"""Exceptions raised by the HMM machinery."""
from __future__ import annotations


class HMMError(Exception):
    """Base for HMM-tier failures."""


class HMMScanError(HMMError):
    """hmmscan invocation or output-parsing failed."""


class HMMDatabaseError(HMMError):
    """HMM database is missing, unreadable, or cannot be indexed."""
