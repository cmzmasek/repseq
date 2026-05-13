"""Global mode: single threshold clustering or select-N across all sequences."""

from __future__ import annotations

from typing import Any, Optional

from ..clustering.diversity import select_diverse
from ..clustering.mmseqs2 import run_clustering
from ..models import Cluster, RunResult, Sequence
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
        clusters = run_clustering(sequences, self.threshold, self.cfg)
        clusters = apply_representative_selection(clusters, self.cfg)
        representatives = [c.representative for c in clusters]
        return RunResult(
            mode="global:threshold",
            representatives=representatives,
            clusters=clusters,
            config_snapshot={"threshold": self.threshold},
        )

    def _run_count(self, sequences: list[Sequence]) -> RunResult:
        selected = select_diverse(sequences, self.n_select, seed=self.seed)
        clusters = [
            Cluster(cluster_id=f"div_{i+1:06d}", representative=s)
            for i, s in enumerate(selected)
        ]
        return RunResult(
            mode="global:count",
            representatives=selected,
            clusters=clusters,
            config_snapshot={"n_select": self.n_select},
        )
