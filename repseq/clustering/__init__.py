"""Clustering backends.

``run_clustering`` is a thin dispatcher that selects the backend from
``cfg['clustering']['backend']`` and forwards to the matching wrapper.

Modes should import ``run_clustering`` from this package — not from
``mmseqs2`` directly — so swapping backends is a config change, not a
code change. Tests that need to mock clustering should patch the symbol
in the *mode* module's namespace (e.g. ``repseq.modes.global_mode.
run_clustering``), since Python re-exports the dispatcher by name there.
"""

from __future__ import annotations

import time
from typing import Any, Optional

import click

from ..models import Cluster, Sequence
from . import cdhit, mmseqs2


def run_clustering(
    sequences: list[Sequence],
    threshold: float,
    cfg: dict[str, Any],
    tmp_dir: Optional[str] = None,
) -> list[Cluster]:
    """Dispatch to the clustering backend named in ``cfg['clustering']['backend']``."""
    backend = cfg.get("clustering", {}).get("backend", "mmseqs2")
    if backend == "mmseqs2":
        return mmseqs2.run_clustering(sequences, threshold, cfg, tmp_dir=tmp_dir)
    if backend == "cdhit":
        return cdhit.run_clustering(sequences, threshold, cfg, tmp_dir=tmp_dir)
    raise ValueError(f"Unknown clustering backend: {backend!r}")


def min_threshold(cfg: dict[str, Any], sequences: list[Sequence]) -> float:
    """Lowest identity threshold the active backend will accept.

    The binary-search wrappers in modes use this as a stopping floor to
    avoid asking cd-hit for thresholds it would refuse (e.g. below 0.80
    for cd-hit-est).
    """
    backend = cfg.get("clustering", {}).get("backend", "mmseqs2")
    if backend == "cdhit":
        return cdhit.min_threshold(sequences, cfg)
    return 0.0


def compute_diversity_curve(
    sequences: list[Sequence],
    cfg: dict[str, Any],
    label: Optional[str] = None,
) -> Optional[dict[float, Optional[int]]]:
    """Cluster ``sequences`` at each configured "standard" identity threshold
    and return ``{cutoff: n_clusters}``. Reporting-only — does not influence
    representative selection.

    Cutoffs below the active backend's identity floor (``min_threshold``)
    map to ``None`` and the backend is not invoked for them. The TSV
    writer renders ``None`` as ``NA``.

    Returns ``None`` (not an empty dict) when the feature is disabled —
    either because ``clustering.diversity_curve_cutoffs`` is missing /
    empty, or because the caller passed no sequences. ``None`` lets the
    TSV writer distinguish "feature off" from "feature on, all cells
    below floor".

    Configured cutoffs are de-duplicated and sorted high → low before use,
    so a cutoff listed twice costs one clustering pass (not two) and the
    console ordering matches the ``_group_counts.tsv`` column ordering,
    which ``output/report.py`` sorts the same way.

    Progress is echoed to the console in the same shape as the modes'
    binary search (``[label] `` tag, per-pass timing). This step runs one
    full clustering pass per cutoff *after* selection has already settled,
    so on a large group it can take minutes — without the echo it reads as
    a hang. ``label`` is the group name the caller is working on (None for
    a single ungrouped pass). Thresholds print with ``:g`` — the form the
    user wrote in YAML and the form the TSV column names use
    (``n_clusters_0.99``) — so the two surfaces cross-reference directly.
    (The binary search prints ``:.4f`` instead because its thresholds are
    continuous search results, not config constants.)
    """
    raw_cutoffs = cfg.get("clustering", {}).get("diversity_curve_cutoffs", []) or []
    if not raw_cutoffs or not sequences:
        return None
    # De-dup + sort descending. Without the de-dup a cutoff repeated in the
    # config would run the backend twice yet still collapse to ONE key in
    # `out` (and one TSV column), so the announced cutoff count would
    # overstate what the run actually delivers.
    cutoffs = sorted({float(c) for c in raw_cutoffs}, reverse=True)
    if not cutoffs:
        return None
    floor = min_threshold(cfg, sequences)
    tag = f"[{label}] " if label else ""
    click.echo(
        f"    {tag}diversity curve (report only): {len(sequences):,} sequence(s) "
        f"at {len(cutoffs)} standard cutoff(s) ..."
    )
    out: dict[float, Optional[int]] = {}
    for i, c in enumerate(cutoffs):
        if c < floor:
            out[c] = None
            click.echo(
                f"      {tag}cutoff {i + 1}/{len(cutoffs)}: threshold={c:g} "
                f"→ skipped (below the backend's {floor:g} floor)"
            )
            continue
        t0 = time.perf_counter()
        clusters = run_clustering(sequences, c, cfg)
        out[c] = len(clusters)
        click.echo(
            f"      {tag}cutoff {i + 1}/{len(cutoffs)}: threshold={c:g} "
            f"→ {len(clusters):,} cluster(s) [{time.perf_counter() - t0:.1f}s]"
        )
    return out


__all__ = [
    "run_clustering", "min_threshold", "compute_diversity_curve",
    "mmseqs2", "cdhit",
]
