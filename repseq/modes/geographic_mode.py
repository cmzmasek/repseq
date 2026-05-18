"""Geographic mode: N representatives per country/region."""

from __future__ import annotations

from typing import Any

from ..clustering import compute_diversity_curve
from ..models import Cluster, GroupStat, RunResult, Sequence
from .base import BaseMode
from .taxonomic1 import _binary_search_threshold


def _group_by_geography(sequences: list[Sequence]) -> dict[str, list[Sequence]]:
    groups: dict[str, list[Sequence]] = {}
    for seq in sequences:
        # country may be "Country: Region" — use only the country part
        country = seq.country or "Unknown"
        key = country.split(":")[0].strip()
        groups.setdefault(key, []).append(seq)
    return groups


class GeographicMode(BaseMode):
    def __init__(
        self,
        cfg: dict[str, Any],
        n_per_country: int,
        overflow: str = "keep",
    ) -> None:
        super().__init__(cfg)
        self.n_per_country = n_per_country
        self.overflow = overflow

    def run(self, sequences: list[Sequence]) -> RunResult:
        groups = _group_by_geography(sequences)
        all_reps: list[Sequence] = []
        all_clusters: list[Cluster] = []
        group_stats: list[GroupStat] = []

        for country, group_seqs in sorted(groups.items()):
            if len(group_seqs) <= self.n_per_country:
                for seq in group_seqs:
                    all_clusters.append(Cluster(cluster_id=f"geo={country}|{seq.id}", representative=seq))
                all_reps.extend(group_seqs)
                group_stats.append(GroupStat(
                    grouping="country", group=country,
                    n_before=len(group_seqs), n_after=len(group_seqs),
                    clustered=False,
                ))
            else:
                reps, threshold = _binary_search_threshold(
                    group_seqs, self.n_per_country, self.cfg, self.overflow,  # type: ignore[arg-type]
                    label=country,
                )
                all_reps.extend(reps)
                for rep in reps:
                    all_clusters.append(Cluster(cluster_id=f"geo={country}|{rep.id}", representative=rep))
                group_stats.append(GroupStat(
                    grouping="country", group=country,
                    n_before=len(group_seqs), n_after=len(reps),
                    clustered=True, cutoff=threshold,
                    cutoff_counts=compute_diversity_curve(group_seqs, self.cfg),
                ))

        return RunResult(
            mode="geographic",
            representatives=all_reps,
            clusters=all_clusters,
            group_stats=group_stats,
            config_snapshot={"n_per_country": self.n_per_country, "overflow": self.overflow},
        )
