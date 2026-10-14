"""Hybrid mode: multi-dimensional stratification combining any grouping keys."""

from __future__ import annotations

from typing import Any, Optional

from ..clustering import compute_diversity_curve
from ..errors import InputError
from ..models import Cluster, GroupStat, RunResult, Sequence
from .base import BaseMode
from .custom_mode import (
    FieldValue,
    MetadataTable,
    assignment_row,
    load_metadata_table,
    report_table_coverage,
    resolve_field,
)
from .taxonomic1 import _binary_search_threshold


# ---------------------------------------------------------------------------
# Multi-key grouping
# ---------------------------------------------------------------------------

def _resolve_all(
    sequences: list[Sequence],
    fields: list[str],
    metadata_table: Optional[MetadataTable],
    field_regexes: Optional[dict[str, str]],
) -> list[list[FieldValue]]:
    """Per sequence, the resolved value of each field (in ``fields`` order)."""
    return [
        [
            resolve_field(seq, field, metadata_table, (field_regexes or {}).get(field))
            for field in fields
        ]
        for seq in sequences
    ]


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
        self.metadata_table: Optional[MetadataTable] = None
        if metadata_table_path:
            self.metadata_table = load_metadata_table(metadata_table_path)
            if not any(self.metadata_table.column_for(f) for f in fields):
                raise InputError(
                    f"None of --fields ({', '.join(fields)}) is a column of "
                    f"{self.metadata_table.path.name}.\n"
                    f"       Columns found: "
                    f"{', '.join(self.metadata_table.columns)}\n"
                    f"       Name at least one of those in --fields (upper / "
                    f"lower case does not matter), or drop --metadata-table."
                )

    def run(self, sequences: list[Sequence]) -> RunResult:
        resolved = _resolve_all(
            sequences, self.fields, self.metadata_table, self.field_regexes
        )
        coverage: dict[str, dict[str, int]] = {}
        if self.metadata_table is not None:
            for j, field in enumerate(self.fields):
                if self.metadata_table.column_for(field) is None:
                    continue
                coverage[field] = report_table_coverage(
                    field, sequences, [vals[j] for vals in resolved],
                    self.metadata_table, self.cfg, "hybrid",
                )

        groups: dict[tuple[str, ...], list[Sequence]] = {}
        assignments: list[dict[str, str]] = []
        for seq, vals in zip(sequences, resolved):
            key = tuple(v.value for v in vals)
            groups.setdefault(key, []).append(seq)
            label = "|".join(f"{f}={v}" for f, v in zip(self.fields, key))
            for field, fv in zip(self.fields, vals):
                assignments.append(assignment_row(seq, field, fv, label))

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
                    cutoff_counts=compute_diversity_curve(
                        group_seqs, self.cfg, label=label),
                ))

        snapshot: dict[str, Any] = {
            "fields": self.fields,
            "n_per_group": self.n_per_group,
            "overflow": self.overflow,
        }
        if self.metadata_table is not None:
            snapshot["metadata_table"] = str(self.metadata_table.path)
            snapshot["table_coverage"] = coverage
        return RunResult(
            mode="hybrid",
            representatives=all_reps,
            clusters=all_clusters,
            group_stats=group_stats,
            config_snapshot=snapshot,
            group_assignments=assignments,
        )
