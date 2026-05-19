"""Time-stratified mode: N representatives per time window."""

from __future__ import annotations

import re
from typing import Any, Optional

from ..clustering import compute_diversity_curve
from ..models import Cluster, GroupStat, RunResult, Sequence
from .base import BaseMode
from .taxonomic1 import _binary_search_threshold

_YEAR_RE = re.compile(r"\b(19\d{2}|20\d{2})\b")


def _extract_year(seq: Sequence) -> Optional[int]:
    # collection_date first, then the strain label (e.g. influenza
    # A/host/place/n/YEAR). The raw header is deliberately NOT scanned: a
    # stray 4-digit number there — a clone ID, an accession, a passage
    # number — is too easily mistaken for a collection year.
    for text in (seq.collection_date or "", seq.strain or ""):
        m = _YEAR_RE.search(text)
        if m:
            return int(m.group(1))
    return None


def _window_label(year: Optional[int], window: str) -> str:
    if year is None:
        return "Unknown"
    if window == "year":
        return str(year)
    if window == "decade":
        return f"{(year // 10) * 10}s"
    # Custom: "5" means 5-year bins
    try:
        size = int(window)
        start = (year // size) * size
        return f"{start}-{start + size - 1}"
    except ValueError:
        return str(year)


def _group_by_time(sequences: list[Sequence], window: str) -> dict[str, list[Sequence]]:
    groups: dict[str, list[Sequence]] = {}
    for seq in sequences:
        year = _extract_year(seq)
        key = _window_label(year, window)
        groups.setdefault(key, []).append(seq)
    return groups


class TimeMode(BaseMode):
    def __init__(
        self,
        cfg: dict[str, Any],
        n_per_window: int,
        window: str = "year",
        overflow: str = "keep",
    ) -> None:
        super().__init__(cfg)
        if window not in ("year", "decade"):
            try:
                size = int(window)
            except (ValueError, TypeError):
                raise ValueError(
                    f"time window must be 'year', 'decade', or a positive "
                    f"integer bin size; got {window!r}"
                )
            if size <= 0:
                raise ValueError(
                    f"time window bin size must be positive; got {window!r}"
                )
        self.n_per_window = n_per_window
        self.window = window
        self.overflow = overflow

    def run(self, sequences: list[Sequence]) -> RunResult:
        groups = _group_by_time(sequences, self.window)
        all_reps: list[Sequence] = []
        all_clusters: list[Cluster] = []
        group_stats: list[GroupStat] = []

        for label, group_seqs in sorted(groups.items()):
            if len(group_seqs) <= self.n_per_window:
                for seq in group_seqs:
                    all_clusters.append(Cluster(cluster_id=f"time={label}|{seq.id}", representative=seq))
                all_reps.extend(group_seqs)
                group_stats.append(GroupStat(
                    grouping="time", group=label,
                    n_before=len(group_seqs), n_after=len(group_seqs),
                    clustered=False,
                ))
            else:
                clusters, threshold = _binary_search_threshold(
                    group_seqs, self.n_per_window, self.cfg, self.overflow,  # type: ignore[arg-type]
                    label=label,
                )
                for c in clusters:
                    c.cluster_id = f"time={label}|{c.representative.id}"
                all_clusters.extend(clusters)
                all_reps.extend(c.representative for c in clusters)
                group_stats.append(GroupStat(
                    grouping="time", group=label,
                    n_before=len(group_seqs), n_after=len(clusters),
                    clustered=True, cutoff=threshold,
                    cutoff_counts=compute_diversity_curve(group_seqs, self.cfg),
                ))

        return RunResult(
            mode="time",
            representatives=all_reps,
            clusters=all_clusters,
            group_stats=group_stats,
            config_snapshot={
                "n_per_window": self.n_per_window,
                "window": self.window,
                "overflow": self.overflow,
            },
        )
