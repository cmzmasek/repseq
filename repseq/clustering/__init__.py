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
        return cdhit.min_threshold(sequences)
    return 0.0


__all__ = ["run_clustering", "min_threshold", "mmseqs2", "cdhit"]
