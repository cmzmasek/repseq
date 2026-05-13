"""Host-stratified mode: N representatives per host organism."""

from __future__ import annotations

from typing import Any

from ..models import Cluster, RunResult, Sequence
from .base import BaseMode
from .taxonomic1 import _binary_search_threshold


def _group_by_host(sequences: list[Sequence]) -> dict[str, list[Sequence]]:
    groups: dict[str, list[Sequence]] = {}
    for seq in sequences:
        key = seq.host or "Unknown"
        groups.setdefault(key, []).append(seq)
    return groups


class HostMode(BaseMode):
    def __init__(
        self,
        cfg: dict[str, Any],
        n_per_host: int,
        overflow: str = "keep",
    ) -> None:
        super().__init__(cfg)
        self.n_per_host = n_per_host
        self.overflow = overflow

    def run(self, sequences: list[Sequence]) -> RunResult:
        groups = _group_by_host(sequences)
        all_reps: list[Sequence] = []
        all_clusters: list[Cluster] = []

        for host, group_seqs in groups.items():
            if len(group_seqs) <= self.n_per_host:
                for seq in group_seqs:
                    all_clusters.append(Cluster(cluster_id=f"host={host}|{seq.id}", representative=seq))
                all_reps.extend(group_seqs)
            else:
                reps, _ = _binary_search_threshold(
                    group_seqs, self.n_per_host, self.cfg, self.overflow  # type: ignore[arg-type]
                )
                all_reps.extend(reps)
                for rep in reps:
                    all_clusters.append(Cluster(cluster_id=f"host={host}|{rep.id}", representative=rep))

        return RunResult(
            mode="host",
            representatives=all_reps,
            clusters=all_clusters,
            config_snapshot={"n_per_host": self.n_per_host, "overflow": self.overflow},
        )
