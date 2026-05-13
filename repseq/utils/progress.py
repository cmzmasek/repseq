"""Progress bar utilities."""

from __future__ import annotations

from contextlib import contextmanager
from typing import Optional

from tqdm import tqdm


@contextmanager
def progress_bar(total: int, desc: str, unit: str = "seq", disable: bool = False):
    """Context manager yielding a tqdm progress bar."""
    bar = tqdm(total=total, desc=desc, unit=unit, disable=disable, ncols=80)
    try:
        yield bar
    finally:
        bar.close()
