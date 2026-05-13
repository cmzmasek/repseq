"""Abstract base class for all selection modes."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

from ..models import RunResult, Sequence


class BaseMode(ABC):
    def __init__(self, cfg: dict[str, Any]) -> None:
        self.cfg = cfg

    @abstractmethod
    def run(self, sequences: list[Sequence]) -> RunResult:
        """Execute the selection mode and return a RunResult."""
        ...

    @property
    def seed(self) -> int:
        return self.cfg.get("seed", 42)

    @property
    def threads(self) -> int:
        return self.cfg.get("threads", 4)
