"""Custom metadata mode: group by arbitrary header fields or external metadata table."""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any, Optional

from ..models import Cluster, RunResult, Sequence
from .base import BaseMode
from .taxonomic1 import _binary_search_threshold


# ---------------------------------------------------------------------------
# Metadata table loader
# ---------------------------------------------------------------------------

def load_metadata_table(path: str | Path) -> dict[str, dict[str, str]]:
    """Load a TSV/CSV metadata table keyed by accession.

    Expected columns: accession (required) + any number of metadata columns.
    """
    path = Path(path)
    delimiter = "\t" if path.suffix in (".tsv", ".txt") else ","
    table: dict[str, dict[str, str]] = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        for row in reader:
            key = row.get("accession") or row.get("Accession") or row.get("id")
            if key:
                table[key.strip()] = {k.strip(): v.strip() for k, v in row.items()}
    return table


# ---------------------------------------------------------------------------
# Grouping
# ---------------------------------------------------------------------------

def _get_field_value(
    seq: Sequence,
    field: str,
    metadata_table: Optional[dict[str, dict[str, str]]],
    regex: Optional[str],
) -> str:
    """Extract a grouping value for a sequence."""
    # 1. External metadata table
    if metadata_table and seq.accession and seq.accession in metadata_table:
        value = metadata_table[seq.accession].get(field)
        if value:
            return value

    # 2. Sequence object attribute
    attr = getattr(seq, field, None)
    if attr:
        return str(attr)

    # 3. Taxonomy rank
    if seq.taxonomy:
        rank_val = seq.taxonomy.get_rank(field)
        if rank_val:
            return rank_val

    # 4. Regex against header
    if regex:
        m = re.search(regex, seq.header, re.IGNORECASE)
        if m:
            try:
                return m.group("value")
            except IndexError:
                if m.lastindex:
                    return m.group(1)
                return m.group(0)

    return "Unknown"


def _group_by_field(
    sequences: list[Sequence],
    field: str,
    metadata_table: Optional[dict[str, dict[str, str]]],
    regex: Optional[str],
) -> dict[str, list[Sequence]]:
    groups: dict[str, list[Sequence]] = {}
    for seq in sequences:
        key = _get_field_value(seq, field, metadata_table, regex)
        groups.setdefault(key, []).append(seq)
    return groups


# ---------------------------------------------------------------------------
# Mode
# ---------------------------------------------------------------------------

class CustomMode(BaseMode):
    """Group sequences by a user-defined field and select N representatives per group.

    Args:
        field: Attribute name, taxonomy rank, or column name in metadata table.
        n_per_group: Target representatives per group.
        metadata_table_path: Optional path to TSV/CSV with per-accession metadata.
        field_regex: Optional regex to extract the field value from the header.
        overflow: "keep" or "trim".
    """

    def __init__(
        self,
        cfg: dict[str, Any],
        field: str,
        n_per_group: int,
        metadata_table_path: Optional[str] = None,
        field_regex: Optional[str] = None,
        overflow: str = "keep",
    ) -> None:
        super().__init__(cfg)
        self.field = field
        self.n_per_group = n_per_group
        self.field_regex = field_regex
        self.overflow = overflow
        self.metadata_table: Optional[dict[str, dict[str, str]]] = None
        if metadata_table_path:
            self.metadata_table = load_metadata_table(metadata_table_path)

    def run(self, sequences: list[Sequence]) -> RunResult:
        groups = _group_by_field(
            sequences, self.field, self.metadata_table, self.field_regex
        )
        all_reps: list[Sequence] = []
        all_clusters: list[Cluster] = []

        for group_label, group_seqs in sorted(groups.items()):
            if len(group_seqs) <= self.n_per_group:
                for seq in group_seqs:
                    all_clusters.append(
                        Cluster(cluster_id=f"{self.field}={group_label}|{seq.id}", representative=seq)
                    )
                all_reps.extend(group_seqs)
            else:
                reps, _ = _binary_search_threshold(
                    group_seqs, self.n_per_group, self.cfg, self.overflow,  # type: ignore[arg-type]
                    label=str(group_label),
                )
                all_reps.extend(reps)
                for rep in reps:
                    all_clusters.append(
                        Cluster(cluster_id=f"{self.field}={group_label}|{rep.id}", representative=rep)
                    )

        return RunResult(
            mode="custom",
            representatives=all_reps,
            clusters=all_clusters,
            config_snapshot={
                "field": self.field,
                "n_per_group": self.n_per_group,
                "overflow": self.overflow,
            },
        )
