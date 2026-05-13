"""Taxonomic Lineage Aware Mode 2: Hierarchical/nested multi-rank clustering.

Cluster at coarser ranks first (family), then within each group at finer ranks
(genus, species), ensuring representation at every level.
"""

from __future__ import annotations

from typing import Any

from ..clustering.mmseqs2 import run_clustering
from ..models import Cluster, RunResult, Sequence
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

        # Start with all sequences, progressively refine
        self._recurse(sequences, self.rank_levels, all_reps, all_clusters, prefix="")

        return RunResult(
            mode="taxonomic2",
            representatives=all_reps,
            clusters=all_clusters,
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
        prefix: str,
    ) -> None:
        if not sequences:
            return

        if not remaining_levels:
            # Leaf: keep all remaining sequences
            for seq in sequences:
                all_reps.append(seq)
                all_clusters.append(Cluster(cluster_id=f"{prefix}|{seq.id}", representative=seq))
            return

        level = remaining_levels[0]
        rank = level["rank"]
        n_per_group = level["n_per_group"]
        next_levels = remaining_levels[1:]

        groups = _group_by_rank(sequences, rank)

        for group_label, group_seqs in groups.items():
            group_prefix = f"{prefix}/{rank}={group_label}" if prefix else f"{rank}={group_label}"

            if len(group_seqs) <= n_per_group:
                # Recurse into next level without clustering
                self._recurse(group_seqs, next_levels, all_reps, all_clusters, group_prefix)
            else:
                reps, threshold = _binary_search_threshold(
                    group_seqs,
                    n_per_group,
                    self.cfg,
                    self.overflow,  # type: ignore[arg-type]
                )
                # Recurse into next level with the representatives
                self._recurse(reps, next_levels, all_reps, all_clusters, group_prefix)
