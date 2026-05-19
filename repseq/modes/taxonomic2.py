"""Taxonomic Lineage Aware Mode 2: Hierarchical/nested multi-rank clustering.

Cluster at coarser ranks first (family), then within each group at finer ranks
(genus, species), ensuring representation at every level.
"""

from __future__ import annotations

from typing import Any, Optional

from ..clustering import compute_diversity_curve, run_clustering
from ..models import Cluster, GroupStat, RunResult, Sequence
from ..representative.selector import apply_representative_selection
from .base import BaseMode
from .taxonomic1 import _group_by_rank, _binary_search_threshold


class TaxonomicMode2(BaseMode):
    """
    Hierarchical clustering: apply ranks in order from coarsest to finest.

    cfg example:
        mode: taxonomic2
        ranks:
          - rank: family
            n_per_group: 20
          - rank: genus
            n_per_group: 10
          - rank: species
            n_per_group: 3
        overflow: keep
    """

    def __init__(
        self,
        cfg: dict[str, Any],
        rank_levels: list[dict[str, Any]],
        overflow: str = "keep",
    ) -> None:
        """
        Args:
            rank_levels: List of dicts with keys 'rank' and 'n_per_group',
                         ordered from coarsest to finest.
        """
        super().__init__(cfg)
        self.rank_levels = rank_levels
        self.overflow = overflow

    def run(self, sequences: list[Sequence]) -> RunResult:
        all_reps: list[Sequence] = []
        all_clusters: list[Cluster] = []
        group_stats: list[GroupStat] = []

        # Start with all sequences, progressively refine
        self._recurse(sequences, self.rank_levels, all_reps, all_clusters,
                      group_stats, prefix="")

        return RunResult(
            mode="taxonomic2",
            representatives=all_reps,
            clusters=all_clusters,
            group_stats=group_stats,
            config_snapshot={
                "rank_levels": self.rank_levels,
                "overflow": self.overflow,
            },
        )

    def _recurse(
        self,
        sequences: list[Sequence],
        remaining_levels: list[dict[str, Any]],
        all_reps: list[Sequence],
        all_clusters: list[Cluster],
        group_stats: list[GroupStat],
        prefix: str,
        represented: Optional[dict[str, list[Sequence]]] = None,
    ) -> None:
        # ``represented`` maps a (currently-surviving) sequence id to the
        # list of *other* original sequences it stands for, accumulated
        # across the hierarchical clustering levels above this one. The
        # hierarchical mode thins by re-clustering survivors at each rank,
        # so a leaf representative can stand for many inputs merged over
        # several levels; threading this accumulator lets the leaf cluster
        # report an accurate size instead of always 1.
        if represented is None:
            represented = {}
        if not sequences:
            return

        if not remaining_levels:
            # Leaf: keep all remaining sequences, each carrying the
            # transitive set of inputs it represents as cluster members.
            for seq in sequences:
                members = represented.get(seq.id, [])
                all_reps.append(seq)
                all_clusters.append(Cluster(
                    cluster_id=f"{prefix}|{seq.id}",
                    representative=seq,
                    members=list(members),
                ))
            return

        level = remaining_levels[0]
        rank = level["rank"]
        n_per_group = level["n_per_group"]
        next_levels = remaining_levels[1:]

        groups = _group_by_rank(sequences, rank)

        for group_label, group_seqs in groups.items():
            group_prefix = f"{prefix}/{rank}={group_label}" if prefix else f"{rank}={group_label}"

            if len(group_seqs) <= n_per_group:
                # No clustering at this level; record the pass-through count
                # before recursing into the next rank.
                group_stats.append(GroupStat(
                    grouping=rank, group=group_prefix,
                    n_before=len(group_seqs), n_after=len(group_seqs),
                    clustered=False,
                ))
                self._recurse(group_seqs, next_levels, all_reps, all_clusters,
                              group_stats, group_prefix, represented)
            else:
                clusters, threshold = _binary_search_threshold(
                    group_seqs,
                    n_per_group,
                    self.cfg,
                    self.overflow,  # type: ignore[arg-type]
                    label=group_prefix,
                )
                # Each cluster's representative now also stands for its
                # members (and, transitively, whatever those members already
                # represented from levels above). Carry that forward so the
                # eventual leaf cluster sizes are accurate.
                next_represented = dict(represented)
                reps: list[Sequence] = []
                for c in clusters:
                    rep = c.representative
                    acc = list(represented.get(rep.id, []))
                    for m in c.members:
                        acc.append(m)
                        acc.extend(represented.get(m.id, []))
                    next_represented[rep.id] = acc
                    reps.append(rep)
                group_stats.append(GroupStat(
                    grouping=rank, group=group_prefix,
                    n_before=len(group_seqs), n_after=len(reps),
                    clustered=True, cutoff=threshold,
                    cutoff_counts=compute_diversity_curve(group_seqs, self.cfg),
                ))
                # Recurse into next level with the representatives
                self._recurse(reps, next_levels, all_reps, all_clusters,
                              group_stats, group_prefix, next_represented)
