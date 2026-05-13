"""Taxonomic Lineage Aware Mode 1: N representatives per taxonomic group.

For each group at the specified rank, binary-search for the MMseqs2 threshold
that yields <= n_per_group representatives, then return those.
"""

from __future__ import annotations

from typing import Any, Literal, Optional

from ..clustering.diversity import select_diverse
from ..clustering.mmseqs2 import run_clustering
from ..models import Cluster, RunResult, Sequence
from ..representative.selector import apply_representative_selection
from .base import BaseMode


OverflowStrategy = Literal["keep", "trim"]


def _group_by_rank(sequences: list[Sequence], rank: str) -> dict[str, list[Sequence]]:
    groups: dict[str, list[Sequence]] = {}
    for seq in sequences:
        if seq.taxonomy:
            label = seq.taxonomy.get_rank(rank)
        else:
            label = None
        key = label or "Unknown"
        groups.setdefault(key, []).append(seq)
    return groups


def _binary_search_threshold(
    sequences: list[Sequence],
    n_target: int,
    cfg: dict[str, Any],
    overflow: OverflowStrategy,
    lo: float = 0.3,
    hi: float = 1.0,
    max_iter: int = 12,
) -> tuple[list[Sequence], float]:
    """Binary-search for threshold that yields <= n_target clusters.

    Returns (representatives, threshold_used).
    """
    best_reps: list[Sequence] = []
    best_threshold = hi

    for _ in range(max_iter):
        mid = (lo + hi) / 2
        clusters = run_clustering(sequences, mid, cfg)
        clusters = apply_representative_selection(clusters, cfg)
        n_reps = len(clusters)

        if n_reps <= n_target:
            best_reps = [c.representative for c in clusters]
            best_threshold = mid
            hi = mid  # try to lower threshold (more clusters) to get closer to n_target
        else:
            lo = mid  # raise threshold (fewer clusters)

        if abs(hi - lo) < 0.005:
            break

    if not best_reps:
        # Could not reach target even at hi=1.0; use all sequences
        best_reps = sequences
        best_threshold = hi

    # If overflow == "trim", diversity-select exactly n_target from best_reps
    if overflow == "trim" and len(best_reps) > n_target:
        seed = cfg.get("seed", 42)
        best_reps = select_diverse(best_reps, n_target, seed=seed)

    return best_reps, best_threshold


class TaxonomicMode1(BaseMode):
    def __init__(
        self,
        cfg: dict[str, Any],
        rank: str,
        n_per_group: int,
        overflow: OverflowStrategy = "keep",
    ) -> None:
        super().__init__(cfg)
        self.rank = rank
        self.n_per_group = n_per_group
        self.overflow = overflow

    def run(self, sequences: list[Sequence]) -> RunResult:
        groups = _group_by_rank(sequences, self.rank)
        all_reps: list[Sequence] = []
        all_clusters: list[Cluster] = []

        for group_label, group_seqs in groups.items():
            if len(group_seqs) <= self.n_per_group:
                # Keep all — no clustering needed
                for seq in group_seqs:
                    all_clusters.append(
                        Cluster(
                            cluster_id=f"{group_label}|{seq.id}",
                            representative=seq,
                        )
                    )
                all_reps.extend(group_seqs)
            else:
                reps, _ = _binary_search_threshold(
                    group_seqs,
                    self.n_per_group,
                    self.cfg,
                    self.overflow,
                )
                all_reps.extend(reps)
                for rep in reps:
                    all_clusters.append(
                        Cluster(
                            cluster_id=f"{group_label}|{rep.id}",
                            representative=rep,
                        )
                    )

        return RunResult(
            mode="taxonomic1",
            representatives=all_reps,
            clusters=all_clusters,
            config_snapshot={
                "rank": self.rank,
                "n_per_group": self.n_per_group,
                "overflow": self.overflow,
            },
        )
