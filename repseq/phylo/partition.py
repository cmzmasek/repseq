"""2E partitioned-supermatrix tree: align each marker family separately,
concatenate column-wise, and let IQ-TREE fit a substitution model per
partition.

This is the default for protein + IQ-TREE runs (``phylo.partition.enabled``,
default true). The concat-then-align path
(:func:`repseq.phylo.pipeline._build_tree`) glues each isolate's marker
proteins into one string and runs **one** MAFFT + **one** model over the
lot — so MAFFT can align L-polymerase residues against M-glycoprotein
residues at the segment seams, and a single substitution model spans
genuinely different protein families. This module instead:

  1. enumerates the marker families (the same HMM tokens used for QC and the
     2F per-protein trees, via
     :func:`repseq.phylo.per_protein.collect_family_specs`);
  2. extracts the satisfying CDS per representative per family and aligns
     each family **separately** (one MAFFT call per family);
  3. concatenates the per-family MSAs **column-wise** into a supermatrix (a
     representative missing a family is gap-filled for that block);
  4. writes a NEXUS partition file declaring one ``charset`` per family;
  5. runs **one IQ-TREE** with ``-p`` / ``-q`` / ``-Q`` (linkage from
     ``phylo.partition.linkage``, default ``-p`` edge-linked proportional)
     so ModelFinder fits a model per partition.

The leaves are representatives — one leaf per rep, exactly as the
whole-genome tree (2E) — so rooting, LCA labelling, colouring, and the rich
phyloXML writer are shared verbatim via
:func:`repseq.phylo.pipeline._finalize_tree`.

Soft-fallback contract: :func:`build_partitioned_phylogeny` returns ``None``
(and :func:`repseq.phylo.pipeline.run_phylogeny` falls back to
concat-then-align) when the run can't be partitioned — the HMM tier didn't
run, or fewer than two families have at least two representatives carrying
them. It raises :class:`repseq.phylo.pipeline.PhyloError` only on a genuine
MAFFT / IQ-TREE failure, matching the rest of the phylo step.

FastTree can't fit partition models, so the dispatcher in
:func:`run_phylogeny` only reaches here when the tool resolves to IQ-TREE;
a FastTree run silently uses concat-then-align.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, Optional

from ..hmm.runner import parse_hmm_token
from ..models import Sequence
from .coloring import ColorScheme
from .iqtree import IQTreeError, run_iqtree
from .iqtree_parse import (
    format_models_for_description,
    parse_chosen_models,
    write_iqtree_model_file,
)
from .labels import format_leaf_label, labeling_options, pick_format_string
from .mafft import MafftError, run_mafft
from . import trimal as trimal_mod
from .trimal import maybe_trim
from .per_protein import (
    _best_satisfying_cds_any,
    _hmm_tier_ran,
    _segment_proteins,
    collect_family_specs,
)
from .pipeline import (
    PhyloError,
    _build_id_map,
    _finalize_tree,
    _write_id_map,
    _write_short_id_fasta,
)

logger = logging.getLogger(__name__)

# IQ-TREE partition linkage flags (Chernomor et al. 2016):
#   -p  edge-linked proportional (one branch-length set + per-partition rate)
#   -q  edge-equal (all partitions share branch lengths)
#   -Q  edge-unlinked (each partition gets independent branch lengths)
_LINKAGE_FLAGS = {"proportional": "-p", "equal": "-q", "unlinked": "-Q"}


def _nexus_safe(name: str) -> str:
    """A NEXUS-identifier-safe charset name (no ``-``, ``.``, or specials).

    The ``--`` multidomain token separator and ``.`` are not legal in a bare
    NEXUS identifier, so collapse any non-word run to a single underscore.
    """
    return re.sub(r"[^A-Za-z0-9_]+", "_", name or "").strip("_") or "part"


def _unique(base: str, used: set[str]) -> str:
    """Disambiguate a NEXUS name that collides after sanitising."""
    name = base
    i = 2
    while name in used:
        name = f"{base}_{i}"
        i += 1
    used.add(name)
    return name


def read_msa(path: Path) -> dict[str, str]:
    """Parse a FASTA MSA into ``{first_header_token: aligned_seq}``.

    Keyed on the short id (the safe first whitespace-separated token the
    MAFFT input carries), which is how the supermatrix rows are matched
    back across the per-family alignments.
    """
    records: dict[str, str] = {}
    cur: Optional[str] = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    records[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        records[cur] = "".join(buf)
    return records


def build_supermatrix(
    fam_msas: list[tuple[str, str, dict[str, str]]],
    short_ids: list[str],
) -> tuple[dict[str, str], list[tuple[str, str, int, int]]]:
    """Concatenate per-family MSAs column-wise into one supermatrix.

    ``fam_msas`` is ``[(family_label, nexus_name, {short_id: aligned_seq}),
    …]`` in partition order. Returns ``(supermatrix, blocks)`` where
    ``supermatrix`` maps each short id to its concatenated row (rows missing
    a family get that block gap-filled with ``-`` × the family's alignment
    width) and ``blocks`` is ``[(family_label, nexus_name, start, end), …]``
    with 1-based inclusive column ranges for the NEXUS charsets.
    """
    blocks: list[tuple[str, str, int, int]] = []
    parts: list[tuple[dict[str, str], int]] = []
    pos = 1
    for label, nexus_name, aligned in fam_msas:
        if not aligned:
            continue
        width = len(next(iter(aligned.values())))
        if width == 0:
            continue
        blocks.append((label, nexus_name, pos, pos + width - 1))
        parts.append((aligned, width))
        pos += width

    supermatrix: dict[str, str] = {}
    for short in short_ids:
        chunks = [aligned.get(short) or ("-" * width) for aligned, width in parts]
        supermatrix[short] = "".join(chunks)
    return supermatrix, blocks


def write_partition_nexus(
    blocks: list[tuple[str, str, int, int]],
    models: dict[str, str],
    path: Path,
) -> bool:
    """Write an IQ-TREE NEXUS partition file. Returns whether it pins models.

    Always emits one ``charset`` per family. The ``charpartition`` is the
    subtle part: IQ-TREE's charpartition syntax requires a concrete
    substitution-model name for **every** partition (``MODEL:charset``) —
    it rejects bare charset names, and ``MFP`` (ModelFinder Plus) is a
    model-*selection strategy*, not a model name, so ``MFP:charset`` makes
    IQ-TREE try to open a file literally named ``MFP`` and abort with
    "File not found MFP". (Empirically confirmed against IQ-TREE 2.3.2.)

    So per-partition ModelFinder cannot be expressed in the charpartition.
    The supported idiom is instead a charsets-only file run with ``-m MFP``
    on the command line, which makes IQ-TREE run ModelFinder independently
    per charset. We therefore:

    * write a ``charpartition`` (with concrete models) **only when every**
      family has a pinned model in ``models`` (keyed by family label) — the
      caller then runs IQ-TREE without ``-m`` so the file's models win;
    * otherwise omit the charpartition entirely (charsets only) — the
      caller passes ``-m MFP`` for per-partition ModelFinder.

    Returns ``True`` in the first case (models pinned in-file), ``False`` in
    the second. When *some but not all* families are pinned, the partial
    pins are dropped with a warning (a charpartition can't mix concrete
    models with per-partition ModelFinder), and the function returns
    ``False``.
    """
    lines = ["#nexus", "begin sets;"]
    for _label, nexus_name, start, end in blocks:
        lines.append(f"    charset {nexus_name} = {start}-{end};")

    pinned = {label: models.get(label) for label, _n, _s, _e in blocks}
    n_pinned = sum(1 for v in pinned.values() if v)
    all_pinned = bool(blocks) and n_pinned == len(blocks)

    if all_pinned:
        assignments = [
            f"{models[label]}:{nexus_name}"
            for label, nexus_name, _start, _end in blocks
        ]
        lines.append(f"    charpartition repseq = {', '.join(assignments)};")
    elif n_pinned:
        logger.warning(
            "[phylo] partition: %d/%d families have a pinned phylo.partition."
            "models entry but the rest do not; IQ-TREE cannot mix pinned "
            "models with per-partition ModelFinder in one charpartition, so "
            "ALL partitions will use ModelFinder (-m MFP). Pin every family "
            "to honour the models map.",
            n_pinned, len(blocks),
        )

    lines.append("end;")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")
    return all_pinned


def _leaf_labels(
    leaf_reps: list[Sequence], id_map: dict[str, str], cfg: dict[str, Any],
) -> dict[str, str]:
    """Short-id → display label, same labelling the 2E/2F leaves use."""
    segmented = any(r.id.startswith("CONCAT|") for r in leaf_reps)
    fmt = pick_format_string(cfg, segmented=segmented)
    opts = labeling_options(cfg)
    return {
        short: format_leaf_label(rep, fmt, **opts)
        for short, rep in zip(id_map.keys(), leaf_reps)
    }


def _write_supermatrix_fasta(
    supermatrix: dict[str, str],
    id_map: dict[str, str],
    labels: dict[str, str],
    path: Path,
) -> None:
    """Write the concatenated supermatrix as ``{prefix}_msa.fasta``.

    Mirrors :func:`repseq.phylo.pipeline._write_short_id_fasta`'s header
    convention: ``>SXXXX <label>`` so viewers show recognisable names while
    the short id stays the safe first token.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for short in id_map:
            raw = labels.get(short, "") or ""
            label = raw.replace("\n", " ").replace("\r", " ").replace("\t", " ").strip()
            fh.write(f">{short} {label}\n" if label else f">{short}\n")
            seq = supermatrix.get(short) or ""
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")


def build_partitioned_phylogeny(
    representatives: list[Sequence],
    cfg: dict[str, Any],
    out_dir: Path,
    prefix: str,
    *,
    color_scheme: Optional[ColorScheme] = None,
) -> Optional[list[Path]]:
    """Build the partitioned-supermatrix tree, or return ``None`` to fall back.

    Returns the output-file list on success; ``None`` when the run can't be
    partitioned (no HMM-resolvable families, fewer than two alignable
    families, or fewer than three placeable representatives), in which case
    :func:`repseq.phylo.pipeline.run_phylogeny` proceeds with
    concat-then-align. Raises :class:`PhyloError` only on a real
    MAFFT / IQ-TREE failure.
    """
    part_cfg = ((cfg or {}).get("phylo", {}) or {}).get("partition", {}) or {}

    specs = collect_family_specs(cfg)
    if not specs:
        logger.info(
            "[phylo] partition: no HMM marker tokens configured (hmms:) — "
            "using concat-then-align"
        )
        return None
    if not _hmm_tier_ran(cfg, representatives):
        logger.info(
            "[phylo] partition: HMM tier did not run this session — "
            "using concat-then-align"
        )
        return None

    # Per family: the satisfying CDS translation per representative. Keep
    # families carried by >= 2 reps (a partition needs >= 2 sequences to be
    # alignable / informative).
    fam_bodies: list[tuple[str, dict[str, str]]] = []
    for family_label, tokens, segment in specs:
        parsed_tokens: list[list[str]] = []
        for token in tokens:
            try:
                parsed_tokens.append(parse_hmm_token(token))
            except ValueError as exc:
                logger.warning("[phylo] partition: skipping token %r: %s", token, exc)
        if not parsed_tokens:
            continue
        bodies: dict[str, str] = {}
        for rep in representatives:
            cds = _best_satisfying_cds_any(
                _segment_proteins(rep, segment), parsed_tokens
            )
            if cds and cds.get("sequence"):
                bodies[rep.id] = cds["sequence"]
        if len(bodies) >= 2:
            fam_bodies.append((family_label, bodies))

    if len(fam_bodies) < 2:
        logger.info(
            "[phylo] partition: %d alignable marker family(ies) (< 2) — "
            "using concat-then-align", len(fam_bodies),
        )
        return None

    # Leaves = reps carrying >= 1 surviving family (a rep with no marker
    # data can't be placed). id_map / leaf_reps stay 1:1 for _finalize_tree.
    carried: set[str] = set()
    for _label, bodies in fam_bodies:
        carried.update(bodies)
    leaf_reps = [r for r in representatives if r.id in carried]
    if len(leaf_reps) < 3:
        logger.info(
            "[phylo] partition: %d placeable representative(s) (< 3) — "
            "using concat-then-align", len(leaf_reps),
        )
        return None

    id_map = _build_id_map(leaf_reps)
    labels = _leaf_labels(leaf_reps, id_map, cfg)
    mafft_cfg = ((cfg or {}).get("phylo", {}).get("mafft", {}) or {})
    mafft_extra = list(mafft_cfg.get("extra_args", []) or [])
    mafft_use_auto = bool(mafft_cfg.get("use_auto", True))
    # The genome-tree trimAl setting governs the per-partition trimming;
    # families are trimmed BEFORE concatenation so the charset ranges
    # reflect the trimmed widths.
    trimal_settings = (cfg.get("phylo", {}) or {}).get("trimal")
    trim_enabled = bool(trimal_settings and trimal_settings.get("enabled"))

    out_dir.mkdir(parents=True, exist_ok=True)

    # One MAFFT (+ optional trimAl) per family. The tree-input alignment is
    # {prefix}_msa_<family>.fasta; when trimming runs, the raw MAFFT output
    # for that family is retained as {prefix}_msa_<family>_untrimmed.fasta.
    fam_msas: list[tuple[str, str, dict[str, str]]] = []        # trimmed (tree input)
    fam_msas_raw: list[tuple[str, str, dict[str, str]]] = []    # untrimmed companion
    extra_files: list[Path] = []
    used_nexus: set[str] = set()
    any_trimmed = False
    for family_label, bodies in fam_bodies:
        fam_reps = [r for r in leaf_reps if r.id in bodies]
        in_fa = out_dir / f"{prefix}_partin_{_nexus_safe(family_label)}.fasta"
        msa_fa = out_dir / f"{prefix}_msa_{family_label}.fasta"
        _write_short_id_fasta(
            fam_reps, id_map, in_fa, bodies=bodies, labels=labels,
        )
        logger.info(
            "[phylo] partition: aligning family %s (%d sequences)…",
            family_label, len(fam_reps),
        )
        try:
            run_mafft(
                in_fa, msa_fa, cfg,
                extra_args=mafft_extra, use_auto=mafft_use_auto,
            )
        except MafftError as exc:
            raise PhyloError(
                f"MAFFT failed on partition {family_label}: {exc}"
            ) from exc
        try:
            in_fa.unlink()
        except OSError:
            pass

        raw_aligned = read_msa(msa_fa)
        if trim_enabled:
            untrimmed_fa = out_dir / f"{prefix}_msa_{family_label}_untrimmed.fasta"
            try:
                if untrimmed_fa.exists():
                    untrimmed_fa.unlink()
                msa_fa.rename(untrimmed_fa)
            except OSError:
                untrimmed_fa = msa_fa
            if maybe_trim(
                untrimmed_fa, msa_fa, cfg, trimal_settings, label=family_label,
            ):
                any_trimmed = True
                extra_files.append(untrimmed_fa)
                aligned = read_msa(msa_fa)
            else:
                # Restore the MAFFT output as the family's tree-input MSA.
                if untrimmed_fa != msa_fa:
                    if msa_fa.exists():
                        try:
                            msa_fa.unlink()
                        except OSError:
                            pass
                    try:
                        untrimmed_fa.rename(msa_fa)
                    except OSError:
                        pass
                aligned = raw_aligned
        else:
            aligned = raw_aligned

        nexus_name = _unique(_nexus_safe(family_label), used_nexus)
        fam_msas.append((family_label, nexus_name, aligned))
        fam_msas_raw.append((family_label, nexus_name, raw_aligned))

    supermatrix, blocks = build_supermatrix(fam_msas, list(id_map.keys()))
    msa_fasta = out_dir / f"{prefix}_msa.fasta"
    _write_supermatrix_fasta(supermatrix, id_map, labels, msa_fasta)

    trim_note = None
    if any_trimmed:
        # Companion untrimmed supermatrix (column ranges differ from the
        # trimmed one, so it carries no charsets — it's an audit artefact).
        raw_super, _ = build_supermatrix(fam_msas_raw, list(id_map.keys()))
        untrimmed_super = out_dir / f"{prefix}_msa_untrimmed.fasta"
        _write_supermatrix_fasta(raw_super, id_map, labels, untrimmed_super)
        extra_files.append(untrimmed_super)
        trim_note = trimal_mod.trim_note(trimal_settings)

    nexus_path = out_dir / f"{prefix}_partition.nex"
    partition_has_models = write_partition_nexus(
        blocks, dict(part_cfg.get("models") or {}), nexus_path
    )

    linkage = part_cfg.get("linkage", "proportional")
    newick_path = out_dir / f"{prefix}_tree.nwk"
    summary_path = out_dir / f"{prefix}_iqtree_summary.txt"
    logger.info(
        "[phylo] partition: %d charsets, linkage=%s (%s) — running IQ-TREE",
        len(blocks), linkage, _LINKAGE_FLAGS.get(linkage, "-p"),
    )
    try:
        run_iqtree(
            msa_fasta, newick_path, cfg,
            is_protein=True,
            summary_path=summary_path,
            partition_file=nexus_path,
            partition_linkage=linkage,
            partition_has_models=partition_has_models,
        )
    except IQTreeError as exc:
        raise PhyloError(str(exc)) from exc

    written_extras = [summary_path] if summary_path.exists() else []
    id_map_path = out_dir / f"{prefix}_tree_id_map.tsv"
    _write_id_map(id_map, id_map_path)

    # Parse the per-partition ModelFinder picks and persist them as a
    # grep-friendly sidecar; the same dict is folded into the
    # phyloXML <description> and the _summary.md Methods section.
    partition_labels = [label for label, _nx, _al in fam_msas]
    chosen_models = (
        parse_chosen_models(summary_path, partition_labels=partition_labels)
        if summary_path.exists()
        else {}
    )
    if chosen_models:
        model_file = out_dir / f"{prefix}_iqtree_model.txt"
        write_iqtree_model_file(chosen_models, model_file)
        if model_file.exists():
            written_extras.append(model_file)

    extra_outputs = (
        [nexus_path]
        + [out_dir / f"{prefix}_msa_{label}.fasta" for label, _nx, _al in fam_msas]
        + extra_files
    )
    linkage_phrase = f"linkage: {_LINKAGE_FLAGS.get(linkage, '-p')} {linkage}"
    if chosen_models:
        picks = format_models_for_description(chosen_models)
        model_label = f"partitioned ({picks}; {linkage_phrase})"
    else:
        model_label = (
            f"partitioned: {len(blocks)} charsets, ModelFinder per partition "
            f"({linkage_phrase})"
        )
    return _finalize_tree(
        representatives=leaf_reps,
        id_map=id_map,
        cfg=cfg,
        out_dir=out_dir,
        file_prefix=prefix,
        xml_name_prefix=prefix,
        is_protein=True,
        tree_tool="iqtree",
        newick_path=newick_path,
        msa_fasta=msa_fasta,
        id_map_path=id_map_path,
        model_label=model_label,
        extra_mafft=mafft_extra,
        written_extras=written_extras,
        color_scheme=color_scheme,
        leaf_protein_ids=None,
        domain_architecture=False,
        input_fasta=None,
        extra_outputs=extra_outputs,
        trim_note=trim_note,
        run_review=True,
        basis_role="genome_partitioned",
        basis_families=[label for label, _nx, _al in fam_msas],
    )
