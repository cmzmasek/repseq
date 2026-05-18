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

from typing import Any, Optional

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
    """
    cutoffs = cfg.get("clustering", {}).get("diversity_curve_cutoffs", []) or []
    if not cutoffs or not sequences:
        return None
    floor = min_threshold(cfg, sequences)
    out: dict[float, Optional[int]] = {}
    for raw in cutoffs:
        c = float(raw)
        if c < floor:
            out[c] = None
            continue
        clusters = run_clustering(sequences, c, cfg)
        out[c] = len(clusters)
    return out


__all__ = [
    "run_clustering", "min_threshold", "compute_diversity_curve",
    "mmseqs2", "cdhit",
]
