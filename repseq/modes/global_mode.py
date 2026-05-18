"""Global mode: single threshold clustering or select-N across all sequences."""

from __future__ import annotations

import time
from typing import Any, Optional

import click

from ..clustering import compute_diversity_curve, run_clustering
from ..clustering.diversity import select_diverse
from ..models import Cluster, GroupStat, RunResult, Sequence
from ..representative.selector import apply_representative_selection
from .base import BaseMode


class GlobalMode(BaseMode):
    """
    Two sub-modes:
      - threshold: cluster at a fixed identity threshold, return representatives
      - count: select exactly N diverse sequences (MaxMin)
    """

    def __init__(
        self,
        cfg: dict[str, Any],
        threshold: Optional[float] = None,
        n_select: Optional[int] = None,
    ) -> None:
        super().__init__(cfg)
        if threshold is None and n_select is None:
            raise ValueError("GlobalMode requires either threshold or n_select")
        self.threshold = threshold
        self.n_select = n_select

    def run(self, sequences: list[Sequence]) -> RunResult:
        if self.threshold is not None:
            return self._run_threshold(sequences)
        return self._run_count(sequences)

    def _run_threshold(self, sequences: list[Sequence]) -> RunResult:
        click.echo(
            f"  clustering {len(sequences)} sequences at threshold "
            f"{self.threshold:.4f} ..."
        )
        t0 = time.perf_counter()
        clusters = run_clustering(sequences, self.threshold, self.cfg)
        clusters = apply_representative_selection(clusters, self.cfg)
        representatives = [c.representative for c in clusters]
        click.echo(
            f"  → {len(representatives)} rep(s) across {len(clusters)} "
            f"cluster(s) [{time.perf_counter() - t0:.1f}s]"
        )
        return RunResult(
            mode="global:threshold",
            representatives=representatives,
            clusters=clusters,
            group_stats=[GroupStat(
                grouping="global", group="(all)",
                n_before=len(sequences), n_after=len(representatives),
                clustered=True, cutoff=self.threshold,
                cutoff_counts=compute_diversity_curve(sequences, self.cfg),
            )],
            config_snapshot={"threshold": self.threshold},
        )

    def _run_count(self, sequences: list[Sequence]) -> RunResult:
        selected = select_diverse(sequences, self.n_select, seed=self.seed)
        clusters = [
            Cluster(cluster_id=f"div_{i+1:06d}", representative=s)
            for i, s in enumerate(selected)
        ]
        # Diversity selection (MaxMin), not threshold clustering — there is
        # no identity cutoff to report.
        return RunResult(
            mode="global:count",
            representatives=selected,
            clusters=clusters,
            group_stats=[GroupStat(
                grouping="global", group="(all)",
                n_before=len(sequences), n_after=len(selected),
                clustered=False,
            )],
            config_snapshot={"n_select": self.n_select},
        )
