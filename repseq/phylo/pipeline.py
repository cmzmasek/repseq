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
from pathlib import Path
from typing import Any, Optional

from Bio import Phylo

from ..models import Sequence, SequenceType
from . import fasttree as fasttree_mod
from . import iqtree as iqtree_mod
from . import mafft as mafft_mod
from .fasttree import FastTreeError, run_fasttree
from .iqtree import IQTreeError, run_iqtree
from .lca import (
    annotate_internal_nodes,
    keep_deepest_labels,
    suppress_same_species_pairs,
)
from .mafft import MafftError, run_mafft
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

    Active when ``clustering.alphabet="protein"`` AND every rep actually
    carries a protein_sequence (a no-resolve fallback or missing-marker
    drop could leave some empty; we never want a half-protein, half-NT
    alignment).
    """
    if cfg is None:
        return False
    if cfg.get("clustering", {}).get("alphabet") != "protein":
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
) -> None:
    """Write a FASTA whose header is the short id and nothing else.

    Mirrors ``clustering.mmseqs2._write_id_fasta`` — sole-token headers
    so the alignment/tree tools cannot lose track of identity at any
    whitespace or pipe boundary. With ``use_protein=True`` the body comes
    from ``rep.protein_sequence`` (the marker / per-isolate AA concat).
    """
    reverse = {v: k for k, v in id_map.items()}
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for rep in reps:
            short = reverse[rep.id]
            fh.write(f">{short}\n")
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
    id_map = _build_id_map(representatives)

    input_fasta = out_dir / f"{prefix}_phylo_input.fasta"
    msa_fasta = out_dir / f"{prefix}_msa.fasta"
    newick_path = out_dir / f"{prefix}_tree.nwk"
    phyloxml_path = out_dir / f"{prefix}_tree.xml"
    id_map_path = out_dir / f"{prefix}_tree_id_map.tsv"
    iqtree_summary_path = out_dir / f"{prefix}_iqtree_summary.txt"

    _write_short_id_fasta(representatives, id_map, input_fasta, use_protein=use_protein)
    _write_id_map(id_map, id_map_path)

    try:
        run_mafft(input_fasta, msa_fasta, cfg)
    except MafftError as exc:
        raise PhyloError(str(exc)) from exc

    written_extras: list[Path] = []
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
    else:
        try:
            run_fasttree(msa_fasta, newick_path, cfg, is_protein=is_protein)
        except FastTreeError as exc:
            raise PhyloError(str(exc)) from exc

    alphabet_label = "protein" if is_protein else "nucleotide"
    if tree_tool == "fasttree":
        model_label = "JTT" if is_protein else "GTR"
    else:
        model_label = _resolved_model(cfg, tree_tool)
    extra_mafft = list(
        ((cfg or {}).get("phylo", {}).get("mafft", {}) or {}).get("extra_args", [])
        or []
    )
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

    # _build_id_map walks representatives in order, so id_map's
    # insertion order matches the representatives list 1:1.
    reps_by_short_id = dict(zip(id_map.keys(), representatives))

    rooting_cfg = (cfg or {}).get("phylo", {}).get("rooting", {}) or {}
    rooting_method_req = rooting_cfg.get("method", "auto")
    try:
        parsed_tree, rooting_method_used = root_tree(
            parsed_tree, reps_by_short_id, method=rooting_method_req,
        )
    except Exception as exc:
        # Rooting is a soft step — fall back to whatever the parser gave us.
        logger.warning("[phylo] rooting failed: %s; leaving tree as parsed", exc)
        rooting_method_used = "none"

    lca_cfg = (cfg or {}).get("phylo", {}).get("lca", {}) or {}
    if lca_cfg.get("enabled", True):
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

    try:
        write_phyloxml(
            None,
            phyloxml_path,
            representatives,
            id_map,
            cfg=cfg or {},
            prefix=prefix,
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
        )
    except Exception as exc:
        raise PhyloError(f"Newick → phyloXML conversion failed: {exc}") from exc

    # The temp input is redundant once the MSA is written.
    try:
        input_fasta.unlink()
    except OSError:
        pass

    return [msa_fasta, newick_path, phyloxml_path, id_map_path] + written_extras
