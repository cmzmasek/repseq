"""Time-stratified mode: N representatives per time window."""

from __future__ import annotations

import re
from typing import Any, Optional

from ..models import Cluster, RunResult, Sequence
from .base import BaseMode
from .taxonomic1 import _binary_search_threshold

_YEAR_RE = re.compile(r"\b(19\d{2}|20\d{2})\b")


def _extract_year(seq: Sequence) -> Optional[int]:
    date = seq.collection_date or ""
    m = _YEAR_RE.search(date)
    if m:
        return int(m.group(1))
    # Try strain/header
    for text in (seq.strain or "", seq.header):
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
        self.n_per_window = n_per_window
        self.window = window
        self.overflow = overflow

    def run(self, sequences: list[Sequence]) -> RunResult:
        groups = _group_by_time(sequences, self.window)
        all_reps: list[Sequence] = []
        all_clusters: list[Cluster] = []

        for label, group_seqs in sorted(groups.items()):
            if len(group_seqs) <= self.n_per_window:
                for seq in group_seqs:
                    all_clusters.append(Cluster(cluster_id=f"time={label}|{seq.id}", representative=seq))
                all_reps.extend(group_seqs)
            else:
                reps, _ = _binary_search_threshold(
                    group_seqs, self.n_per_window, self.cfg, self.overflow  # type: ignore[arg-type]
                )
                all_reps.extend(reps)
                for rep in reps:
                    all_clusters.append(Cluster(cluster_id=f"time={label}|{rep.id}", representative=rep))

        return RunResult(
            mode="time",
            representatives=all_reps,
            clusters=all_clusters,
            config_snapshot={
                "n_per_window": self.n_per_window,
                "window": self.window,
                "overflow": self.overflow,
            },
        )
