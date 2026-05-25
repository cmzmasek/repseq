"""Orchestrator for ``--conservation-heatmap``.

Walks the same marker-spec enumeration that drives the 2F per-protein
trees, locates each family's already-aligned MSA under
``{prefix}_per_protein/``, picks the longest-non-gap row as the
reference, projects that reference rep's HMM hits onto the MSA
columns, and hands off to
:func:`repseq.viz.conservation.write_conservation_heatmap` for the
actual PNG. Outputs land under ``{prefix}_conservation/`` as
``{prefix}_<family>.png``.

**Hard requirement**: ``--per-protein-phylo`` must have run earlier
in the same invocation so each family's ``<family>_msa.fasta`` is on
disk. The flag is intentionally independent — a user who has already
run ``--per-protein-phylo`` once and just wants new heatmaps from
the same outputs can re-invoke with only ``--conservation-heatmap``,
provided the output directory still carries the per-protein MSAs.
Soft-fails (stderr line) per family when its MSA is missing.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Optional

from ..hmm.runner import parse_hmm_token
from ..models import Sequence
from ..phylo.per_protein import (
    _best_satisfying_cds_any,
    _segment_proteins,
    collect_family_specs,
    overlap_tolerance_from_cfg,
)
from .conservation import (
    family_color_hex,
    pick_reference_row,
    project_hits_to_alignment,
    read_msa,
    write_conservation_heatmap,
)

logger = logging.getLogger(__name__)


def _load_id_map(path: Path) -> dict[str, str]:
    """Read ``{prefix}_per_protein/<family>_tree_id_map.tsv`` (short_id → rep_id).

    The file is a two-column TSV with a header — same format
    :func:`repseq.phylo.pipeline._write_id_map` produces for both
    2E and every 2F family. Returns ``{}`` on any read error.
    """
    out: dict[str, str] = {}
    try:
        with open(path) as fh:
            header = fh.readline()  # short_id\taccession
            if not header.startswith("short_id"):
                # Older / unexpected format — try to parse anyway.
                parts = header.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    out[parts[0]] = parts[1]
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0]:
                    out[parts[0]] = parts[1]
    except OSError:
        return {}
    return out


def _reference_hits(
    rep: Sequence,
    tokens: list[str],
    segment: Optional[str],
    tol: int,
) -> tuple[list[dict], Optional[dict]]:
    """Return ``(all_hits_on_satisfying_cds, satisfying_cds_dict)``.

    Re-runs the same CDS-picker used by 2F so the projected hits are
    those of the *same* CDS the heatmap's MSA leaf carries. The full
    ``hmm_hits`` list (passing + non-passing) is returned so the
    domain ribbon can show every architectural box the user might
    recognise — passing-only would hide non-required accessories that
    are still informative.
    """
    parsed = [parse_hmm_token(t) for t in tokens]
    proteins = _segment_proteins(rep, segment)
    cds = _best_satisfying_cds_any(
        proteins, parsed, overlap_tolerance=tol,
    )
    if cds is None:
        return [], None
    return list(cds.get("hmm_hits") or []), cds


def run_conservation_heatmaps(
    representatives: list[Sequence],
    cfg: dict[str, Any],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    """Build one heatmap PNG per declared marker spec.

    Returns the list of paths written. Empty list on configuration
    errors (no specs declared, no 2F output directory found) — the
    caller emits a single stderr line then moves on. Per-family
    failures are logged but do not stop sibling families.
    """
    specs = collect_family_specs(cfg)
    if not specs:
        raise FileNotFoundError(
            "no marker specs configured (hmms: under cluster_protein / "
            "segment_markers) — nothing to plot"
        )
    per_protein_dir = out_dir / f"{prefix}_per_protein"
    if not per_protein_dir.exists():
        raise FileNotFoundError(
            f"{per_protein_dir} not found — --conservation-heatmap "
            "requires --per-protein-phylo to have run first"
        )

    cons_dir = out_dir / f"{prefix}_conservation"
    cons_dir.mkdir(parents=True, exist_ok=True)
    tol = overlap_tolerance_from_cfg(cfg)
    rep_by_id = {r.id: r for r in representatives}
    family_labels = [lab for lab, _t, _s in specs]
    written: list[Path] = []

    for family_label, tokens, segment in specs:
        msa_path = per_protein_dir / f"{family_label}_msa.fasta"
        id_map_path = per_protein_dir / f"{family_label}_tree_id_map.tsv"
        if not msa_path.exists():
            logger.warning(
                "[conservation] %s_msa.fasta not found — skipped "
                "(did --per-protein-phylo build this family?)",
                family_label,
            )
            continue
        msa = read_msa(msa_path)
        if not msa:
            logger.warning("[conservation] %s: empty MSA — skipped", family_label)
            continue

        # Look up reps by their short id via the per-family id map. A
        # missing/unreadable id map is non-fatal — we still build the
        # heatmap, just without the domain ribbon.
        id_map = _load_id_map(id_map_path) if id_map_path.exists() else {}

        short_id = pick_reference_row(msa)
        ref_rep: Optional[Sequence] = None
        if short_id and short_id in id_map:
            ref_rep = rep_by_id.get(id_map[short_id])

        projected: list[dict] = []
        if ref_rep is not None:
            hits, _cds = _reference_hits(ref_rep, tokens, segment, tol)
            if hits and short_id in msa:
                projected = project_hits_to_alignment(msa[short_id], hits)
        else:
            logger.info(
                "[conservation] %s: no reference rep resolvable from id "
                "map — drawing heatmap without domain ribbon",
                family_label,
            )

        color = family_color_hex(family_label, family_labels)
        out_png = cons_dir / f"{prefix}_{family_label}.png"
        try:
            result_path = write_conservation_heatmap(
                msa_path,
                out_png=out_png,
                family_label=family_label,
                family_color=color,
                hmm_hits_on_reference=projected,
            )
        except ImportError as exc:
            # Matplotlib unavailable. Surface once and bail on the
            # whole step — every family would fail the same way.
            raise
        except Exception as exc:
            logger.warning(
                "[conservation] %s failed: %s — skipped", family_label, exc,
            )
            continue
        if result_path is not None:
            written.append(result_path)

    return written
