"""Hybrid mode: multi-dimensional stratification combining any grouping keys."""

from __future__ import annotations

from itertools import product
from typing import Any, Optional

from ..clustering import compute_diversity_curve
from ..models import Cluster, GroupStat, RunResult, Sequence
from .base import BaseMode
from .custom_mode import _get_field_value, load_metadata_table
from .taxonomic1 import _binary_search_threshold


# ---------------------------------------------------------------------------
# Multi-key grouping
# ---------------------------------------------------------------------------

def _multi_group(
    sequences: list[Sequence],
    fields: list[str],
    metadata_table: Optional[dict[str, dict[str, str]]],
    field_regexes: Optional[dict[str, str]],
) -> dict[tuple[str, ...], list[Sequence]]:
    """Group sequences by a composite key of multiple fields."""
    groups: dict[tuple[str, ...], list[Sequence]] = {}
    for seq in sequences:
        key_parts: list[str] = []
        for field in fields:
            regex = (field_regexes or {}).get(field)
            key_parts.append(_get_field_value(seq, field, metadata_table, regex))
        key = tuple(key_parts)
        groups.setdefault(key, []).append(seq)
    return groups


# ---------------------------------------------------------------------------
# Mode
# ---------------------------------------------------------------------------

class HybridMode(BaseMode):
    """
    Stratify by any combination of fields (e.g. genus × host × decade) and
    select N representatives per stratum.

    Config example:
        mode: hybrid
        fields:
          - genus
          - host
          - year_window: decade
        n_per_group: 5
        overflow: keep
        metadata_table: /path/to/metadata.tsv
        field_regexes:
          genotype: "genotype[=:]?\\s*(\\w+)"
    """

    def __init__(
        self,
        cfg: dict[str, Any],
        fields: list[str],
        n_per_group: int,
        overflow: str = "keep",
        metadata_table_path: Optional[str] = None,
        field_regexes: Optional[dict[str, str]] = None,
    ) -> None:
        super().__init__(cfg)
        self.fields = fields
        self.n_per_group = n_per_group
        self.overflow = overflow
        self.field_regexes = field_regexes or {}
        self.metadata_table: Optional[dict[str, dict[str, str]]] = None
        if metadata_table_path:
            self.metadata_table = load_metadata_table(metadata_table_path)

    def run(self, sequences: list[Sequence]) -> RunResult:
        groups = _multi_group(
            sequences, self.fields, self.metadata_table, self.field_regexes
        )
        all_reps: list[Sequence] = []
        all_clusters: list[Cluster] = []
        group_stats: list[GroupStat] = []
        grouping = "|".join(self.fields)

        for key_tuple, group_seqs in sorted(groups.items()):
            label = "|".join(f"{f}={v}" for f, v in zip(self.fields, key_tuple))

            if len(group_seqs) <= self.n_per_group:
                for seq in group_seqs:
                    all_clusters.append(
                        Cluster(cluster_id=f"{label}|{seq.id}", representative=seq)
                    )
                all_reps.extend(group_seqs)
                group_stats.append(GroupStat(
                    grouping=grouping, group=label,
                    n_before=len(group_seqs), n_after=len(group_seqs),
                    clustered=False,
                ))
            else:
                clusters, threshold = _binary_search_threshold(
                    group_seqs, self.n_per_group, self.cfg, self.overflow,  # type: ignore[arg-type]
                    label=label,
                )
                for c in clusters:
                    c.cluster_id = f"{label}|{c.representative.id}"
                all_clusters.extend(clusters)
                all_reps.extend(c.representative for c in clusters)
                group_stats.append(GroupStat(
                    grouping=grouping, group=label,
                    n_before=len(group_seqs), n_after=len(clusters),
                    clustered=True, cutoff=threshold,
                    cutoff_counts=compute_diversity_curve(group_seqs, self.cfg),
                ))

        return RunResult(
            mode="hybrid",
            representatives=all_reps,
            clusters=all_clusters,
            group_stats=group_stats,
            config_snapshot={
                "fields": self.fields,
                "n_per_group": self.n_per_group,
                "overflow": self.overflow,
            },
        )
