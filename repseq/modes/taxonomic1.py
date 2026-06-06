"""Taxonomic Lineage Aware Mode 1: N representatives per taxonomic group.

For each group at the specified rank, binary-search for the MMseqs2 threshold
that yields <= n_per_group representatives, then return those.
"""

from __future__ import annotations

import logging
import time
from typing import Any, Literal, Optional

import click

from ..clustering import compute_diversity_curve, min_threshold, run_clustering
from ..clustering.diversity import select_diverse
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
) -> tuple[list[Cluster], float]:
    """Binary-search for the clustering threshold that yields ~n_target clusters.

    Identity-threshold semantics (same direction for both supported
    backends): a *higher* threshold makes clustering stricter and produces
    *more*, smaller clusters; a *lower* threshold merges more aggressively
    into *fewer* clusters. The search moves accordingly.

    Cluster count is a step function of the threshold, so an exact landing
    on ``n_target`` is often impossible. The result kept is the one whose
    count is largest while still ``<= n_target`` (closest from below); if
    no threshold gets at or below the target, the closest result from
    above is returned instead. A significant undershoot is logged as a
    warning — silently returning far fewer representatives than requested
    is a real footgun.

    The lower bound is clamped to ``min_threshold(cfg, sequences)`` so the
    search never asks the backend for a threshold it would refuse — cd-hit
    rejects ``-c`` below 0.40 (protein) or 0.80 (nucleotide); mmseqs2's
    floor is 0.0, so this is a no-op there.

    Returns ``(clusters, threshold_used)`` — the full ``Cluster`` objects
    with their members intact (representative accessible via
    ``c.representative``), so callers can report accurate cluster sizes
    instead of treating every representative as a singleton.
    """
    floor = min_threshold(cfg, sequences)
    if floor > lo:
        lo = floor
    if lo >= hi:
        raise ValueError(
            f"Binary search cannot run: clustering backend floor ({floor}) "
            f"is at or above the search upper bound ({hi}). Lower the "
            f"backend's threshold floor or switch backends."
        )

    best_clusters: list[Cluster] = []
    best_count = -1
    best_threshold = hi
    fallback_clusters: list[Cluster] = []
    fallback_count: Optional[int] = None  # smallest count seen that exceeds n_target

    tag = f"[{label}] " if label else ""
    click.echo(
        f"  {tag}clustering {len(sequences)} sequences "
        f"(target = {n_target} reps) ..."
    )

    for i in range(max_iter):
        mid = (lo + hi) / 2
        t0 = time.perf_counter()
        clusters = run_clustering(sequences, mid, cfg)
        clusters = apply_representative_selection(clusters, cfg)
        n_reps = len(clusters)
        dt = time.perf_counter() - t0
        click.echo(
            f"    {tag}iter {i + 1}/{max_iter}: threshold={mid:.4f} "
            f"→ {n_reps} cluster(s) [{dt:.1f}s]"
        )

        if n_reps <= n_target:
            # Few enough — keep the count closest to the target (largest
            # that is still <= n_target), then raise the threshold to
            # split further.
            if n_reps > best_count:
                best_clusters = clusters
                best_count = n_reps
                best_threshold = mid
            lo = mid
        else:
            # Too many clusters — lower the threshold to merge more.
            if fallback_count is None or n_reps < fallback_count:
                fallback_clusters = clusters
                fallback_count = n_reps
            hi = mid

        if n_reps == n_target or abs(hi - lo) < 0.005:
            break

    if not best_clusters:
        # No threshold reached at-or-below the target; return the closest
        # result from above (overflow="trim" will pare it down exactly).
        if fallback_clusters:
            best_clusters = fallback_clusters
            best_count = fallback_count if fallback_count is not None else len(fallback_clusters)
        else:
            # Nothing clustered at all — every sequence is its own singleton.
            best_clusters = [
                Cluster(cluster_id=f"cluster_{i + 1:06d}", representative=s)
                for i, s in enumerate(sequences)
            ]
            best_count = len(best_clusters)
        best_threshold = hi

    where = f" for '{label}'" if label else ""
    if best_count < n_target:
        click.echo(
            f"    {tag}settled on {best_count} rep(s) at threshold "
            f"{best_threshold:.4f} (target {n_target} not reachable — "
            f"cluster count is a step function of identity)"
        )
        logger.warning(
            "Binary search%s: requested %d representatives but cluster count "
            "is a step function of identity; closest achievable is %d at "
            "threshold %.3f.",
            where, n_target, best_count, best_threshold,
        )
    else:
        click.echo(
            f"    {tag}settled on {best_count} rep(s) at threshold "
            f"{best_threshold:.4f}"
        )
        logger.info(
            "Binary search%s: %d representatives at threshold %.3f.",
            where, best_count, best_threshold,
        )

    # If overflow == "trim", diversity-select exactly n_target clusters by
    # picking the most diverse representatives, then keeping their clusters
    # (members intact).
    if overflow == "trim" and len(best_clusters) > n_target:
        seed = cfg.get("seed", 42)
        reps = [c.representative for c in best_clusters]
        chosen = select_diverse(reps, n_target, seed=seed, cfg=cfg)
        chosen_ids = {s.id for s in chosen}
        best_clusters = [c for c in best_clusters if c.representative.id in chosen_ids]

    return best_clusters, best_threshold


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
                clusters, threshold = _binary_search_threshold(
                    group_seqs,
                    self.n_per_group,
                    self.cfg,
                    self.overflow,
                    label=group_label,
                )
                for c in clusters:
                    c.cluster_id = f"{group_label}|{c.representative.id}"
                all_clusters.extend(clusters)
                all_reps.extend(c.representative for c in clusters)
                group_stats.append(GroupStat(
                    grouping=self.rank, group=group_label,
                    n_before=len(group_seqs), n_after=len(clusters),
                    clustered=True, cutoff=threshold,
                    cutoff_counts=compute_diversity_curve(group_seqs, self.cfg),
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
