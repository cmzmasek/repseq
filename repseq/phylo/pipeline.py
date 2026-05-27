"""Phylogeny orchestrator: short-id rename → MAFFT → IQ-TREE/FastTree → phyloXML.

The pipeline:

1. Skip with a warning if there are fewer than 3 representatives —
   neither IQ-TREE nor FastTree can build a tree on a pair or singleton.
2. Decide protein vs nucleotide from the first representative's
   ``seq_type`` (the upstream pipeline guarantees a single alphabet).
3. Pick the tree builder: IQ-TREE for protein, FastTree for NT, unless
   overridden by ``phylo.tool``.
4. Assign deterministic short ids (``S0001``, ``S0002``, …) so that
   downstream tools — many of which choke on long names, pipes, or
   whitespace — see only clean tokens. The mapping is written to
   ``{prefix}_tree_id_map.tsv`` so the MSA and Newick (which keep the
   short ids) remain readable.
5. Run MAFFT → MSA FASTA (short ids retained).
6. Run IQ-TREE or FastTree → Newick (short ids retained).
7. Render the Newick to phyloXML via the rich writer in
   :mod:`repseq.phylo.phyloxml_writer` — each leaf gets a formatted
   ``<name>``, a ``<taxonomy>`` block, a ``<sequence>`` block with the
   GenBank accession, and ``repseq:``-namespaced ``<property>``
   elements for host, collection_date, country, strain, isolate_id,
   year, and four taxonomy ranks. The tree itself is ladderized and
   the ``<phylogeny>`` element carries a ``<name>`` and ``<description>``
   recording the tools, versions, model, and bootstrap settings used.

Outputs (all under ``{prefix}_*``):
    {prefix}_msa.fasta           aligned MSA, short-id headers
    {prefix}_tree.nwk            tree-builder Newick, short-id leaves
    {prefix}_tree.xml            phyloXML, rich per-leaf annotation
    {prefix}_tree_id_map.tsv     short_id<TAB>accession
    {prefix}_iqtree_summary.txt  IQ-TREE ModelFinder report (IQ-TREE only)

The orchestrator never raises out of the click command — it catches its
own subprocess and conversion errors and reports them to stderr, matching
the existing ``--plot`` behaviour.
"""

from __future__ import annotations

import logging
import sys
import time
from pathlib import Path
from typing import Any, Optional

from Bio import Phylo

from ..models import Sequence, SequenceType
from . import fasttree as fasttree_mod
from . import iqtree as iqtree_mod
from . import mafft as mafft_mod
from .fasttree import FastTreeError, run_fasttree
from .iqtree import IQTreeError, run_iqtree
from .iqtree_parse import (
    format_models_for_description,
    parse_chosen_models,
    write_iqtree_model_file,
)
from .labels import format_leaf_label, labeling_options, pick_format_string
from .lca import (
    annotate_internal_nodes,
    keep_deepest_labels,
    suppress_same_species_pairs,
)
from .coloring import ColorScheme, build_color_scheme
from .mafft import MafftError, run_mafft
from . import trimal as trimal_mod
from .trimal import maybe_trim
from .phyloxml_writer import write_phyloxml
from .rooting import root_tree

logger = logging.getLogger(__name__)


class PhyloError(RuntimeError):
    """Raised when the phylogeny step cannot proceed — surfaced to stderr,
    not propagated, so a tree failure never destroys an otherwise good run."""


_SHORT_ID_FMT = "S{:04d}"


def _pick_tree_tool(cfg: Optional[dict[str, Any]], is_protein: bool) -> str:
    """Choose ``"iqtree"`` or ``"fasttree"`` for this run.

    ``phylo.tool`` is the override (`auto` | `iqtree` | `fasttree`); the
    `auto` default delegates to alphabet — IQ-TREE for protein,
    FastTree for nucleotide. The alphabet check matches the rest of the
    pipeline: protein alignments benefit from IQ-TREE's ModelFinder
    (JTT vs WAG vs LG can change topology), while FastTree's `-nt -gtr`
    is well-understood and orders of magnitude faster on NT.
    """
    if cfg is None:
        return "iqtree" if is_protein else "fasttree"
    tool = (cfg.get("phylo", {}) or {}).get("tool", "auto")
    if tool == "iqtree":
        return "iqtree"
    if tool == "fasttree":
        return "fasttree"
    return "iqtree" if is_protein else "fasttree"


def _use_protein_sequence(reps: list[Sequence], cfg: Optional[dict[str, Any]]) -> bool:
    """True when the phylo input should be each rep's protein_sequence.

    Active when ``clustering.alphabet_for_clustering="protein"`` AND every
    rep actually carries a protein_sequence (a no-resolve fallback or
    missing-marker drop could leave some empty; we never want a
    half-protein, half-NT alignment).
    """
    if cfg is None:
        return False
    if cfg.get("clustering", {}).get("alphabet_for_clustering") != "protein":
        return False
    return all(r.protein_sequence for r in reps)


def _is_protein(reps: list[Sequence]) -> bool:
    """Pick a single alphabet for the whole rep set when reading seq.sequence.

    The upstream pipeline never mixes protein and nucleotide reps in a
    single run, so checking the first one is enough; unknown defaults
    to protein since FastTree's protein model is the default and works
    on anything alphabetic.
    """
    return reps[0].seq_type != SequenceType.NUCLEOTIDE


def _build_id_map(reps: list[Sequence]) -> dict[str, str]:
    """Assign ``S0001``, ``S0002``, … to each representative.

    Returns a short_id → original seq.id mapping. Order is stable
    (the input order), so re-running on the same reps produces
    identical ids — convenient for diff-based comparison.
    """
    return {_SHORT_ID_FMT.format(i + 1): rep.id for i, rep in enumerate(reps)}


def _write_short_id_fasta(
    reps: list[Sequence],
    id_map: dict[str, str],
    path: Path,
    use_protein: bool = False,
    labels: Optional[dict[str, str]] = None,
    bodies: Optional[dict[str, str]] = None,
) -> None:
    """Write a FASTA whose primary header token is the short id.

    The short id is the first (and possibly only) whitespace-separated
    token so phylo binaries (MAFFT, FastTree, IQ-TREE) — which all key
    on the first token — can never confuse identity at a special-char
    or whitespace boundary. When ``labels`` carries a non-empty entry
    for a short id, that label is appended as the FASTA *description*
    (after a single space); MAFFT preserves the description verbatim
    into the output MSA, where it's visible in alignment viewers like
    AliView / Jalview / MEGA but ignored by tree-builders.

    Body selection, in priority order:
      * ``bodies`` (keyed by ``seq.id``) — an explicit per-rep string,
        used by the per-protein-tree path (2F) where the leaf body is a
        specific CDS translation rather than a whole-rep sequence;
      * else ``rep.protein_sequence`` when ``use_protein`` (marker / AA
        concat, 2E protein alphabet);
      * else ``rep.sequence`` (2E nucleotide).
    """
    reverse = {v: k for k, v in id_map.items()}
    path.parent.mkdir(parents=True, exist_ok=True)
    labels = labels or {}
    with open(path, "w") as fh:
        for rep in reps:
            short = reverse[rep.id]
            # Newlines/carriage returns in the label would split the
            # header across lines, breaking the FASTA. Tabs are fine in
            # a FASTA description but normalise to a space for legibility.
            raw_label = labels.get(short, "") or ""
            label = raw_label.replace("\n", " ").replace("\r", " ") \
                             .replace("\t", " ").strip()
            header = f">{short} {label}" if label else f">{short}"
            fh.write(header + "\n")
            if bodies is not None:
                seq = bodies.get(rep.id) or ""
            else:
                seq = (rep.protein_sequence if use_protein else rep.sequence) or ""
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def _write_id_map(id_map: dict[str, str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("short_id\taccession\n")
        for short, original in id_map.items():
            fh.write(f"{short}\t{original}\n")


def _resolved_model(cfg: Optional[dict[str, Any]], tree_tool: str) -> Optional[str]:
    """Return the substitution model recorded in cfg for ``tree_tool``.

    Used only for the phyloXML description block — no runtime effect.
    FastTree's model is inferred from the alphabet (NT → GTR, AA → JTT)
    and we record it as such, since FastTree itself doesn't echo it.
    """
    if cfg is None:
        return None
    phylo_cfg = cfg.get("phylo", {}) or {}
    if tree_tool == "iqtree":
        return (phylo_cfg.get("iqtree", {}) or {}).get("model") or "MFP"
    return None


def _resolved_ufboot(cfg: Optional[dict[str, Any]], tree_tool: str) -> Optional[int]:
    if cfg is None or tree_tool != "iqtree":
        return None
    return (cfg.get("phylo", {}) or {}).get("iqtree", {}).get(
        "ultrafast_bootstrap", 1000
    )


def run_phylogeny(
    representatives: list[Sequence],
    cfg: dict[str, Any],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    """Run the full phylogeny pipeline on a list of representative sequences.

    Returns the list of files written so the caller can append them to its
    output manifest. Raises :class:`PhyloError` when the step cannot
    proceed (under 3 reps, missing binary, or subprocess failure) — the
    caller is expected to catch and report rather than crash the run.
    """
    n = len(representatives)
    if n < 3:
        raise PhyloError(
            f"need >= 3 representatives to build a tree, got {n}"
        )

    use_protein = _use_protein_sequence(representatives, cfg)
    is_protein = use_protein or _is_protein(representatives)
    tree_tool = _pick_tree_tool(cfg, is_protein)
    color_scheme = build_color_scheme(representatives, cfg)

    # Partitioned supermatrix (default on for protein + IQ-TREE runs):
    # align each marker family separately and let IQ-TREE fit a model per
    # partition, instead of gluing the markers into one string and fitting a
    # single model over the lot. Soft-falls back to concat-then-align below
    # (returns None) when the run can't be partitioned — FastTree, no
    # HMM-resolvable families, or fewer than two alignable families.
    part_cfg = ((cfg or {}).get("phylo", {}) or {}).get("partition", {}) or {}
    if part_cfg.get("enabled", True) and use_protein and tree_tool == "iqtree":
        from .partition import build_partitioned_phylogeny
        partitioned = build_partitioned_phylogeny(
            representatives, cfg, out_dir, prefix, color_scheme=color_scheme,
        )
        if partitioned is not None:
            return partitioned

    bodies = {
        rep.id: ((rep.protein_sequence if use_protein else rep.sequence) or "")
        for rep in representatives
    }
    return _build_tree(
        representatives,
        bodies,
        is_protein=is_protein,
        cfg=cfg,
        out_dir=out_dir,
        file_prefix=prefix,
        xml_name_prefix=prefix,
        color_scheme=color_scheme,
        trimal_settings=(cfg.get("phylo", {}) or {}).get("trimal"),
    )


def _build_tree(
    representatives: list[Sequence],
    bodies: dict[str, str],
    *,
    is_protein: bool,
    cfg: dict[str, Any],
    out_dir: Path,
    file_prefix: str,
    xml_name_prefix: str,
    color_scheme: Optional[ColorScheme] = None,
    leaf_protein_ids: Optional[dict[str, set[str]]] = None,
    mafft_extra_args: Optional[list[str]] = None,
    mafft_use_auto: bool = True,
    domain_architecture: bool = False,
    trimal_settings: Optional[dict[str, Any]] = None,
) -> list[Path]:
    """Build one MSA + tree + phyloXML from a list of leaves.

    The shared engine behind both the whole-representative tree (2E,
    :func:`run_phylogeny`) and each per-protein-family tree (2F,
    :func:`repseq.phylo.per_protein.run_per_protein_phylogeny`).

    ``bodies`` maps ``seq.id`` → the AA/NT string to align for that leaf
    (a whole-rep sequence for 2E, a specific CDS translation for 2F).
    Output files are ``{file_prefix}_{msa.fasta,tree.nwk,tree.xml,
    tree_id_map.tsv,iqtree_summary.txt}`` under ``out_dir``;
    ``xml_name_prefix`` becomes the phyloXML ``<phylogeny>`` name.
    ``color_scheme`` (when given) is a shared taxonomy-colour palette
    forwarded to the writer; passing the *same* scheme into every tree of
    a run keeps a given taxon the same colour across 2E and all 2F trees.
    ``leaf_protein_ids`` (per-protein trees) maps ``seq.id`` → the set of
    CDS protein ids that fed this tree, so the writer shows only those as
    ``<sequence type="protein">``; ``None`` (2E) shows the full set.
    ``mafft_extra_args`` / ``mafft_use_auto`` override the MAFFT command
    (the per-protein path passes its L-INS-i args with auto off); ``None``
    falls back to ``phylo.mafft.extra_args`` with ``--auto``.
    ``domain_architecture`` (per-protein trees) makes the writer emit a
    ``<domain_architecture>`` from each shown CDS's HMM hits.
    Raises :class:`PhyloError` on any binary/parse failure.
    """
    tree_tool = _pick_tree_tool(cfg, is_protein)
    id_map = _build_id_map(representatives)

    input_fasta = out_dir / f"{file_prefix}_phylo_input.fasta"
    msa_fasta = out_dir / f"{file_prefix}_msa.fasta"
    newick_path = out_dir / f"{file_prefix}_tree.nwk"
    phyloxml_path = out_dir / f"{file_prefix}_tree.xml"
    id_map_path = out_dir / f"{file_prefix}_tree_id_map.tsv"
    iqtree_summary_path = out_dir / f"{file_prefix}_iqtree_summary.txt"

    # Build a per-rep description label using the same labeling format
    # that drives the tree leaves. This lands as the FASTA description
    # (after the short-id primary token), so AliView/Jalview/MEGA show
    # recognisable names while phylo binaries keep using the safe
    # short id.
    segmented = any(r.id.startswith("CONCAT|") for r in representatives)
    fmt = pick_format_string(cfg, segmented=segmented)
    opts = labeling_options(cfg)
    short_ids = list(id_map.keys())  # insertion-ordered, matches reps
    labels = {
        short: format_leaf_label(rep, fmt, **opts)
        for short, rep in zip(short_ids, representatives)
    }
    _write_short_id_fasta(
        representatives, id_map, input_fasta,
        bodies=bodies, labels=labels,
    )
    _write_id_map(id_map, id_map_path)

    # Resolve the MAFFT command once — used for both the run and the
    # phyloXML provenance description. A caller-supplied override (the
    # per-protein L-INS-i default) drops --auto so its strategy flags
    # take effect; otherwise --auto + phylo.mafft.extra_args (2E).
    if mafft_extra_args is not None:
        mafft_args_used = list(mafft_extra_args)
        mafft_auto_used = mafft_use_auto
    else:
        mafft_cfg = ((cfg or {}).get("phylo", {}).get("mafft", {}) or {})
        mafft_args_used = list(mafft_cfg.get("extra_args", []) or [])
        mafft_auto_used = bool(mafft_cfg.get("use_auto", True))

    try:
        run_mafft(
            input_fasta, msa_fasta, cfg,
            extra_args=mafft_args_used, use_auto=mafft_auto_used,
        )
    except MafftError as exc:
        raise PhyloError(str(exc)) from exc

    written_extras: list[Path] = []

    # Optional trimAl trimming between MAFFT and the tree-builder. On
    # success {file_prefix}_msa.fasta becomes the trimmed (tree-input)
    # alignment and the raw MAFFT output is retained as _msa_untrimmed.fasta;
    # on a soft-fail (disabled / binary missing / degenerate) the MAFFT
    # output stays as _msa.fasta and the tree is built untrimmed.
    trim_note_str: Optional[str] = None
    if trimal_settings and trimal_settings.get("enabled"):
        untrimmed = out_dir / f"{file_prefix}_msa_untrimmed.fasta"
        try:
            if untrimmed.exists():
                untrimmed.unlink()
            msa_fasta.rename(untrimmed)
        except OSError:
            untrimmed = msa_fasta
        if maybe_trim(untrimmed, msa_fasta, cfg, trimal_settings, label=file_prefix):
            written_extras.append(untrimmed)
            trim_note_str = trimal_mod.trim_note(trimal_settings)
        else:
            # Restore the MAFFT output as the canonical (untrimmed) MSA.
            if untrimmed != msa_fasta:
                if msa_fasta.exists():
                    try:
                        msa_fasta.unlink()
                    except OSError:
                        pass
                try:
                    untrimmed.rename(msa_fasta)
                except OSError:
                    pass

    chosen_models: dict[str, str] = {}
    if tree_tool == "iqtree":
        try:
            run_iqtree(
                msa_fasta, newick_path, cfg,
                is_protein=is_protein,
                summary_path=iqtree_summary_path,
            )
        except IQTreeError as exc:
            raise PhyloError(str(exc)) from exc
        if iqtree_summary_path.exists():
            written_extras.append(iqtree_summary_path)
            # ModelFinder records its pick deep in the .iqtree report;
            # extract it for the phyloXML description, _summary.md, and a
            # grep-friendly sidecar. Soft-fail (empty dict) when the
            # format changes — the tree itself is unaffected.
            chosen_models = parse_chosen_models(iqtree_summary_path)
            if chosen_models:
                model_file = out_dir / f"{file_prefix}_iqtree_model.txt"
                write_iqtree_model_file(chosen_models, model_file)
                if model_file.exists():
                    written_extras.append(model_file)
    else:
        try:
            run_fasttree(msa_fasta, newick_path, cfg, is_protein=is_protein)
        except FastTreeError as exc:
            raise PhyloError(str(exc)) from exc

    if tree_tool == "fasttree":
        model_label = "JTT" if is_protein else "GTR"
    else:
        # Prefer the actual ModelFinder pick when we parsed one; fall
        # back to the configured input (typically "MFP") so the
        # description always says *something*.
        model_label = (
            format_models_for_description(chosen_models)
            or _resolved_model(cfg, tree_tool)
        )
    return _finalize_tree(
        representatives=representatives,
        id_map=id_map,
        cfg=cfg,
        out_dir=out_dir,
        file_prefix=file_prefix,
        xml_name_prefix=xml_name_prefix,
        is_protein=is_protein,
        tree_tool=tree_tool,
        newick_path=newick_path,
        msa_fasta=msa_fasta,
        id_map_path=id_map_path,
        model_label=model_label,
        extra_mafft=mafft_args_used,
        written_extras=written_extras,
        color_scheme=color_scheme,
        leaf_protein_ids=leaf_protein_ids,
        domain_architecture=domain_architecture,
        input_fasta=input_fasta,
        trim_note=trim_note_str,
    )


def _finalize_tree(
    *,
    representatives: list[Sequence],
    id_map: dict[str, str],
    cfg: dict[str, Any],
    out_dir: Path,
    file_prefix: str,
    xml_name_prefix: str,
    is_protein: bool,
    tree_tool: str,
    newick_path: Path,
    msa_fasta: Path,
    id_map_path: Path,
    model_label: Optional[str],
    extra_mafft: list[str],
    written_extras: list[Path],
    color_scheme: Optional[ColorScheme] = None,
    leaf_protein_ids: Optional[dict[str, set[str]]] = None,
    domain_architecture: bool = False,
    input_fasta: Optional[Path] = None,
    extra_outputs: Optional[list[Path]] = None,
    trim_note: Optional[str] = None,
) -> list[Path]:
    """Shared post-tree tail: parse Newick → root → LCA → phyloXML.

    Both the concat-then-align builder (:func:`_build_tree`) and the
    partitioned-supermatrix builder
    (:func:`repseq.phylo.partition.build_partitioned_phylogeny`) end here, so
    rooting, LCA labelling, taxonomy colouring, and the rich phyloXML writer
    are identical regardless of how the MSA + Newick were produced. The
    leaves are representatives (``id_map`` insertion order must match
    ``representatives`` 1:1). ``input_fasta`` (when given) is the temp MSA
    input to unlink; ``extra_outputs`` are appended to the returned file
    list (the partitioned path passes its NEXUS + per-family MSAs here).
    Returns ``[msa, nwk, xml, id_map] + written_extras + extra_outputs``.
    """
    phyloxml_path = out_dir / f"{file_prefix}_tree.xml"
    alphabet_label = "protein" if is_protein else "nucleotide"
    extra_tree = list(
        ((cfg or {}).get("phylo", {}).get(tree_tool, {}) or {}).get(
            "extra_args", [],
        )
        or []
    )

    # Load the Newick once. Root and LCA-annotate before handing to
    # the writer so both rooting choice and internal labels make it
    # into the final phyloXML.
    try:
        parsed_tree = Phylo.read(str(newick_path), "newick")
    except Exception as exc:
        raise PhyloError(f"could not parse Newick {newick_path}: {exc}") from exc

    reps_by_short_id = dict(zip(id_map.keys(), representatives))

    rooting_cfg = (cfg or {}).get("phylo", {}).get("rooting", {}) or {}
    rooting_method_req = rooting_cfg.get("method", "auto")
    print(
        f"[phylo] starting rooting (method={rooting_method_req})",
        file=sys.stderr,
    )
    t_root = time.time()
    try:
        parsed_tree, rooting_method_used = root_tree(
            parsed_tree, reps_by_short_id, method=rooting_method_req,
            outgroup=rooting_cfg.get("outgroup"),
            outgroup_rank=rooting_cfg.get("outgroup_rank"),
        )
    except Exception as exc:
        # Rooting is a soft step — fall back to whatever the parser gave us.
        logger.warning("[phylo] rooting failed: %s; leaving tree as parsed", exc)
        rooting_method_used = "none"
    print(
        f"[phylo] rooting finished ({time.time() - t_root:.1f}s, "
        f"used={rooting_method_used})",
        file=sys.stderr,
    )

    lca_cfg = (cfg or {}).get("phylo", {}).get("lca", {}) or {}
    if lca_cfg.get("enabled", True):
        print("[phylo] starting LCA annotation", file=sys.stderr)
        t_lca = time.time()
        try:
            annotate_internal_nodes(
                parsed_tree, reps_by_short_id,
                min_rank=lca_cfg.get("min_rank", "genus"),
                coverage_threshold=lca_cfg.get("coverage_threshold", 0.5),
            )
            keep_deepest_labels(parsed_tree)
            suppress_same_species_pairs(parsed_tree, reps_by_short_id)
        except Exception as exc:
            logger.warning("[phylo] LCA annotation failed: %s", exc)
        print(
            f"[phylo] LCA annotation finished ({time.time() - t_lca:.1f}s)",
            file=sys.stderr,
        )

    try:
        write_phyloxml(
            None,
            phyloxml_path,
            representatives,
            id_map,
            cfg=cfg or {},
            prefix=xml_name_prefix,
            alphabet=alphabet_label,
            msa_tool="MAFFT",
            msa_version=mafft_mod.tool_version(),
            tree_tool="IQ-TREE" if tree_tool == "iqtree" else "FastTree",
            tree_version=(
                iqtree_mod.tool_version(
                    ((cfg or {}).get("phylo", {}).get("iqtree", {}) or {}).get("binary"),
                )
                if tree_tool == "iqtree"
                else fasttree_mod.tool_version()
            ),
            model=model_label,
            ufboot=_resolved_ufboot(cfg, tree_tool),
            extra_msa_args=extra_mafft,
            extra_tree_args=extra_tree,
            tree=parsed_tree,
            rooting_method=rooting_method_used,
            color_scheme=color_scheme,
            leaf_protein_ids=leaf_protein_ids,
            domain_architecture=domain_architecture,
            trim_note=trim_note,
        )
    except Exception as exc:
        raise PhyloError(f"Newick → phyloXML conversion failed: {exc}") from exc

    # The temp input is redundant once the MSA is written.
    if input_fasta is not None:
        try:
            input_fasta.unlink()
        except OSError:
            pass

    return (
        [msa_fasta, newick_path, phyloxml_path, id_map_path]
        + written_extras
        + (extra_outputs or [])
    )
