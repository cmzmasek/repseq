"""HMM-based marker-protein selection (HMMER hmmscan wrapper).

Public API:
    - ``resolve_database_path``, ``ensure_pressed``, ``db_signature``,
      ``parse_ga_cutoffs``, ``profile_count`` — database management.
    - ``is_available``, ``scan`` — hmmscan invocation.
    - ``BUNDLED_DB_PATH`` — path to the bundled viral-core profile set.
    - ``HMMError``, ``HMMScanError``, ``HMMDatabaseError`` — exceptions.

See ``repseq.clustering.marker`` for how HMM hits feed marker
selection.
"""
from .database import (
    BUNDLED_DB_PATH,
    HMMPRESS_INDEX_SUFFIXES,
    db_signature,
    ensure_pressed,
    has_press_index,
    parse_ga_cutoffs,
    profile_count,
    resolve_database_path,
)
from .errors import HMMDatabaseError, HMMError, HMMScanError
from .hmmscan import is_available, scan
from .runner import (
    CACHE_SOURCE,
    coverage_of,
    get_ga_cutoffs,
    passes_cutoffs,
    scan_proteins,
)

__all__ = [
    "BUNDLED_DB_PATH",
    "CACHE_SOURCE",
    "HMMPRESS_INDEX_SUFFIXES",
    "HMMDatabaseError",
    "HMMError",
    "HMMScanError",
    "coverage_of",
    "db_signature",
    "ensure_pressed",
    "get_ga_cutoffs",
    "has_press_index",
    "is_available",
    "parse_ga_cutoffs",
    "passes_cutoffs",
    "profile_count",
    "resolve_database_path",
    "scan",
    "scan_proteins",
]
