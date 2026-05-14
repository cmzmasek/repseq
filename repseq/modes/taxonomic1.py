"""Taxonomic Lineage Aware Mode 1: N representatives per taxonomic group.

For each group at the specified rank, binary-search for the MMseqs2 threshold
that yields <= n_per_group representatives, then return those.
"""

from __future__ import annotations

import logging
from typing import Any, Literal, Optional

from ..clustering.diversity import select_diverse
from ..clustering.mmseqs2 import run_clustering
from ..models import Cluster, GroupStat, RunResult, Sequence
from ..representative.selector import apply_representative_selection
from .base import BaseMode

logger = logging.getLogger(__name__)

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
    label: str = "",
) -> tuple[list[Sequence], float]:
    """Binary-search for the MMseqs2 threshold that yields ~n_target clusters.

    MMseqs2 ``--min-seq-id`` semantics: a *higher* identity threshold makes
    clustering stricter and produces *more*, smaller clusters; a *lower*
    threshold merges more aggressively into *fewer* clusters. The search
    moves accordingly.

    Cluster count is a step function of the threshold, so an exact landing
    on ``n_target`` is often impossible. The result kept is the one whose
    count is largest while still ``<= n_target`` (closest from below); if
    no threshold gets at or below the target, the closest result from
    above is returned instead. A significant undershoot is logged as a
    warning — silently returning far fewer representatives than requested
    is a real footgun.

    Returns (representatives, threshold_used).
    """
    best_reps: list[Sequence] = []
    best_count = -1
    best_threshold = hi
    fallback_reps: list[Sequence] = []
    fallback_count: Optional[int] = None  # smallest count seen that exceeds n_target

    for _ in range(max_iter):
        mid = (lo + hi) / 2
        clusters = run_clustering(sequences, mid, cfg)
        clusters = apply_representative_selection(clusters, cfg)
        n_reps = len(clusters)

        if n_reps <= n_target:
            # Few enough — keep the count closest to the target (largest
            # that is still <= n_target), then raise the threshold to
            # split further.
            if n_reps > best_count:
                best_reps = [c.representative for c in clusters]
                best_count = n_reps
                best_threshold = mid
            lo = mid
        else:
            # Too many clusters — lower the threshold to merge more.
            if fallback_count is None or n_reps < fallback_count:
                fallback_reps = [c.representative for c in clusters]
                fallback_count = n_reps
            hi = mid

        if n_reps == n_target or abs(hi - lo) < 0.005:
            break

    if not best_reps:
        # No threshold reached at-or-below the target; return the closest
        # result from above (overflow="trim" will pare it down exactly).
        best_reps = fallback_reps or list(sequences)
        best_count = fallback_count if fallback_count is not None else len(sequences)
        best_threshold = hi

    where = f" for '{label}'" if label else ""
    if best_count < n_target:
        logger.warning(
            "Binary search%s: requested %d representatives but cluster count "
            "is a step function of identity; closest achievable is %d at "
            "threshold %.3f.",
            where, n_target, best_count, best_threshold,
        )
    else:
        logger.info(
            "Binary search%s: %d representatives at threshold %.3f.",
            where, best_count, best_threshold,
        )

    # If overflow == "trim", diversity-select exactly n_target from best_reps.
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
        group_stats: list[GroupStat] = []

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
                group_stats.append(GroupStat(
                    grouping=self.rank, group=group_label,
                    n_before=len(group_seqs), n_after=len(group_seqs),
                    clustered=False,
                ))
            else:
                reps, threshold = _binary_search_threshold(
                    group_seqs,
                    self.n_per_group,
                    self.cfg,
                    self.overflow,
                    label=group_label,
                )
                all_reps.extend(reps)
                for rep in reps:
                    all_clusters.append(
                        Cluster(
                            cluster_id=f"{group_label}|{rep.id}",
                            representative=rep,
                        )
                    )
                group_stats.append(GroupStat(
                    grouping=self.rank, group=group_label,
                    n_before=len(group_seqs), n_after=len(reps),
                    clustered=True, cutoff=threshold,
                ))

        return RunResult(
            mode="taxonomic1",
            representatives=all_reps,
            clusters=all_clusters,
            group_stats=group_stats,
            config_snapshot={
                "rank": self.rank,
                "n_per_group": self.n_per_group,
                "overflow": self.overflow,
            },
        )
