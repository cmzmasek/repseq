"""2F — one phylogenetic tree per declared HMM domain-architecture token.

Where :func:`repseq.phylo.pipeline.run_phylogeny` (2E) builds a *single*
tree from each representative's whole sequence (or marker-protein concat),
this module builds *one tree per protein family*, where a "family" is an
HMM domain-architecture **token** the user already declared for QC
(``virus.segment_markers[seg].hmms`` / ``virus.cluster_protein[seg].hmms``
for segmented runs, ``clustering.cluster_protein[*].hmms`` otherwise).

For each token we walk the selected representatives, find the CDS that
satisfies the token (via :func:`repseq.hmm.runner.cds_satisfies_token` —
so ``Bunya_G1--Bunya_G2`` is one family and intervening/extra domains are
tolerated exactly as the QC gate tolerates them), and build an MSA + tree
+ phyloXML on those protein translations. The heavy lifting (MAFFT →
IQ-TREE/FastTree → root → LCA → phyloXML) is the same
:func:`repseq.phylo.pipeline._build_tree` engine 2E uses, so rooting, LCA
labelling, ``phylo.tool``, and per-leaf annotation all carry over.

Requirements (soft-fail, mirroring ``--plot`` / ``--phylo``):
  * The HMM tier must have run this session (``hmm.enabled`` + configured
    ``hmms:`` + populated ``hmm_hits`` on the proteins). Without it no CDS
    can satisfy any token; we raise :class:`PhyloError` with that reason.
  * Each family needs at least ``phylo.per_protein.min_taxa`` (default 3,
    the tree-builder floor) representatives carrying the architecture;
    sparser families are skipped with a log note.

Outputs land in a ``{prefix}_per_protein/`` subdirectory, one set per
built family: ``<family>_msa.fasta``, ``<family>_tree.nwk``,
``<family>_tree.xml``, ``<family>_tree_id_map.tsv`` (+
``<family>_iqtree_summary.txt`` for IQ-TREE). The family name is the
sanitised token, prefixed with the segment in segmented mode
(``S_Bunya_nucleocap``) so two segments declaring the same token never
collide.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path
from typing import Any, Optional

from ..config import get_virus_config
from ..hmm.runner import cds_satisfies_token, parse_hmm_token
from ..models import Sequence
from .coloring import build_color_scheme
from .pipeline import PhyloError, _build_tree

logger = logging.getLogger(__name__)


def _sanitize(name: str) -> str:
    """Filesystem-safe family component.

    Keeps alphanumerics and ``-`` / ``_`` / ``.`` (so the multidomain
    ``--`` separator survives intact), collapses any other run — notably
    whitespace and path separators — to a single ``_``.
    """
    return re.sub(r"[^A-Za-z0-9._-]+", "_", (name or "").strip()).strip("_") or "family"


def _resolve_segment_tokens(
    seg_name: str,
    segment_markers: dict,
    cluster_protein_per_seg: dict,
) -> list[str]:
    """Tokens gating a segment — mirrors ``cli._resolve_segment_hmms``.

    Replicated here (rather than imported) because ``cli`` imports the
    phylo package, so importing back would be circular. Lookup order:
    ``segment_markers[seg]`` first (v0.13+ HMM-aware form), then any
    dict-form ``cluster_protein[seg]`` entry's ``hmms``.
    """
    if seg_name in segment_markers:
        spec = segment_markers[seg_name] or {}
        return list(spec.get("hmms") or [])
    if seg_name in cluster_protein_per_seg:
        tokens: list[str] = []
        for entry in cluster_protein_per_seg[seg_name] or []:
            if isinstance(entry, dict):
                tokens.extend(entry.get("hmms") or [])
        return tokens
    return []


def collect_family_specs(cfg: dict[str, Any]) -> list[tuple[str, str, Optional[str]]]:
    """Return ``[(family_label, token, segment_or_None), …]`` for the run.

    ``segment`` scopes the CDS search: in segmented mode a token declared
    on segment ``S`` only ever pulls CDS from the ``S`` segment of each
    isolate; ``None`` (non-segmented) searches the representative's own
    proteins. Duplicate ``(segment, token)`` pairs are dropped, first
    occurrence wins, declaration order is preserved.
    """
    specs: list[tuple[str, str, Optional[str]]] = []
    seen: set[tuple[Optional[str], str]] = set()

    segmented = bool((cfg.get("segmented", {}) or {}).get("enabled"))
    if segmented:
        virus = get_virus_config(cfg) or {}
        segment_markers = virus.get("segment_markers") or {}
        cluster_protein = virus.get("cluster_protein") or {}
        segments = list(virus.get("segments") or [])
        # Include any segment that only appears under the marker maps.
        for extra in (*segment_markers.keys(), *cluster_protein.keys()):
            if extra not in segments:
                segments.append(extra)
        for seg in segments:
            for tok in _resolve_segment_tokens(seg, segment_markers, cluster_protein):
                key = (seg, tok)
                if key in seen:
                    continue
                seen.add(key)
                label = f"{_sanitize(seg)}_{_sanitize(tok)}"
                specs.append((label, tok, seg))
    else:
        for entry in (cfg.get("clustering", {}) or {}).get("cluster_protein", []) or []:
            if not isinstance(entry, dict):
                continue
            for tok in entry.get("hmms") or []:
                key = (None, tok)
                if key in seen:
                    continue
                seen.add(key)
                specs.append((_sanitize(tok), tok, None))
    return specs


def _segment_proteins(rep: Sequence, segment: Optional[str]) -> list[dict]:
    """Proteins to search for ``rep`` under family ``segment`` scope.

    Non-segmented (``segment is None``): the rep's own ``proteins``.
    Segmented: the matching entry in ``concat_segments`` (the per-segment
    Sequences carried on a CONCAT rep); falls back to the rep itself when
    a non-CONCAT rep happens to carry that segment label.
    """
    if segment is None:
        return rep.proteins or []
    for seg in rep.concat_segments or []:
        if seg.segment == segment:
            return seg.proteins or []
    if rep.segment == segment:
        return rep.proteins or []
    return []


def _best_satisfying_cds(proteins: list[dict], token_hmms: list[str]) -> Optional[dict]:
    """The CDS in ``proteins`` that best satisfies the token, or None.

    "Best" = longest translation (mirrors ``select_marker_protein``'s
    default), tie-broken by the lowest worst-domain E-value the token
    matching reports. Proteins without a translation can't seed a tree
    leaf and are skipped.
    """
    candidates: list[tuple[dict, float]] = []
    for prot in proteins:
        seq = prot.get("sequence")
        if not seq:
            continue
        worst_e = cds_satisfies_token(prot.get("hmm_hits") or [], token_hmms)
        if worst_e is not None:
            candidates.append((prot, worst_e))
    if not candidates:
        return None
    candidates.sort(
        key=lambda pe: (-(pe[0].get("length") or len(pe[0]["sequence"])), pe[1])
    )
    return candidates[0][0]


def _min_taxa(cfg: dict[str, Any]) -> int:
    """Per-family minimum, never below the tree-builder floor of 3."""
    raw = ((cfg or {}).get("phylo", {}) or {}).get("per_protein", {}) or {}
    try:
        return max(3, int(raw.get("min_taxa", 3)))
    except (TypeError, ValueError):
        return 3


def _hmm_tier_ran(cfg: dict[str, Any], representatives: list[Sequence]) -> bool:
    """Did the HMM scan actually annotate any protein this session?

    Prefer the runtime flag the pipeline stashes (``_hmm_runtime.active``);
    fall back to scanning for any populated ``hmm_hits`` so a warm cache
    or an unusual entry point is still recognised.
    """
    if ((cfg or {}).get("_hmm_runtime") or {}).get("active"):
        return True
    for rep in representatives:
        pools = [rep.proteins or []]
        pools.extend((seg.proteins or []) for seg in rep.concat_segments or [])
        for pool in pools:
            if any(p.get("hmm_hits") for p in pool):
                return True
    return False


def run_per_protein_phylogeny(
    representatives: list[Sequence],
    cfg: dict[str, Any],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    """Build one tree per declared domain-architecture token.

    Returns the list of files written (across all families). Raises
    :class:`PhyloError` — for the caller to catch and report, never crash —
    when the step can't proceed at all: no tokens configured, the HMM tier
    didn't run, or no family cleared ``min_taxa``.
    """
    specs = collect_family_specs(cfg)
    if not specs:
        raise PhyloError(
            "no HMM marker tokens configured (hmms:) — per-protein trees need "
            "at least one declared domain architecture in segment_markers / "
            "cluster_protein"
        )
    if not _hmm_tier_ran(cfg, representatives):
        raise PhyloError(
            "the HMM tier did not run this session (need hmm.enabled plus "
            "configured hmms:), so no CDS can be assigned to a protein family"
        )

    min_taxa = _min_taxa(cfg)
    sub_dir = out_dir / f"{prefix}_per_protein"
    written: list[Path] = []
    built = 0
    built_labels: list[str] = []
    sparse: list[str] = []

    # Build the leaf-colour palette once over the FULL representative set
    # (not each family's subset) so a taxon keeps the same colour across
    # every per-protein tree and the whole-genome tree (2E) — what makes
    # cross-tree incongruence legible at a glance.
    color_scheme = build_color_scheme(representatives, cfg)

    for family_label, token, segment in specs:
        try:
            parsed = parse_hmm_token(token)
        except ValueError as exc:
            logger.warning("[per-protein] skipping malformed token %r: %s", token, exc)
            continue

        leaf_reps: list[Sequence] = []
        bodies: dict[str, str] = {}
        for rep in representatives:
            cds = _best_satisfying_cds(_segment_proteins(rep, segment), parsed)
            if cds is None:
                continue
            leaf_reps.append(rep)
            bodies[rep.id] = cds["sequence"]

        if len(leaf_reps) < min_taxa:
            sparse.append(f"{family_label} ({len(leaf_reps)})")
            logger.info(
                "[per-protein] family %s: %d representative(s) carry it "
                "(< %d) — skipped", family_label, len(leaf_reps), min_taxa,
            )
            continue

        logger.info(
            "[per-protein] building %s from %d representative(s)…",
            family_label, len(leaf_reps),
        )
        try:
            files = _build_tree(
                leaf_reps,
                bodies,
                is_protein=True,
                cfg=cfg,
                out_dir=sub_dir,
                file_prefix=family_label,
                xml_name_prefix=f"{prefix}_{family_label}",
                color_scheme=color_scheme,
            )
        except PhyloError as exc:
            logger.warning("[per-protein] family %s failed: %s", family_label, exc)
            continue
        written.extend(files)
        built += 1
        built_labels.append(family_label)

    if built == 0:
        detail = f" (sparse families: {', '.join(sparse)})" if sparse else ""
        raise PhyloError(
            f"no protein family had >= {min_taxa} representatives carrying the "
            f"declared architecture — nothing built{detail}"
        )

    inc_path = _write_incongruence_table(
        built_labels, sub_dir, out_dir, prefix, cfg,
    )
    if inc_path is not None:
        written.append(inc_path)

    return written


def _write_incongruence_table(
    built_labels: list[str],
    sub_dir: Path,
    out_dir: Path,
    prefix: str,
    cfg: dict[str, Any],
) -> Optional[Path]:
    """Pairwise unrooted-RF table across the built family trees (+ the 2E
    whole-genome tree when ``--phylo`` wrote one). Soft-fail: returns the
    path written, or None when disabled / too few trees / on any error.

    Reads each tree's on-disk Newick + ``_tree_id_map.tsv`` (the genome
    tree's land in ``out_dir`` from the ``--phylo`` step, which always
    runs before this one), so no tree objects need threading in.
    """
    pp_cfg = ((cfg or {}).get("phylo", {}) or {}).get("per_protein", {}) or {}
    if not pp_cfg.get("incongruence", True):
        return None
    try:
        trees: list[tuple[str, Path, Path]] = [
            (
                label,
                sub_dir / f"{label}_tree.nwk",
                sub_dir / f"{label}_tree_id_map.tsv",
            )
            for label in built_labels
        ]
        genome_nwk = out_dir / f"{prefix}_tree.nwk"
        genome_map = out_dir / f"{prefix}_tree_id_map.tsv"
        if genome_nwk.exists() and genome_map.exists():
            trees.append(("GENOME", genome_nwk, genome_map))
        if len(trees) < 2:
            return None  # need a pair to compare
        from .incongruence import compute_incongruence, write_incongruence_tsv

        rows = compute_incongruence(trees)
        if not rows:
            return None
        inc_path = sub_dir / f"{prefix}_incongruence.tsv"
        write_incongruence_tsv(rows, inc_path)
        logger.info(
            "[per-protein] wrote pairwise incongruence table (%d pair(s)): %s",
            len(rows), inc_path,
        )
        return inc_path
    except Exception as exc:  # never let a metric bug void the trees
        logger.warning("[per-protein] incongruence table failed: %s", exc)
        return None
