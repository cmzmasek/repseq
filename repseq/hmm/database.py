"""HMM database resolution, indexing, and GA-cutoff parsing.

The database file is a concatenated HMMER3 ``.hmm`` text file (one or
more profiles separated by ``//`` records). ``hmmscan`` requires the
binary indexes ``.h3f / .h3i / .h3m / .h3p`` produced by ``hmmpress``;
this module auto-presses on first use when those are missing.
"""
from __future__ import annotations

import hashlib
import shutil
import subprocess
from pathlib import Path
from typing import Optional

from .errors import HMMDatabaseError

# Bundled set ships under repseq/data/hmms/. Resolved relative to this
# file so editable installs (``pip install -e .``) work without a copy.
BUNDLED_DB_PATH = (
    Path(__file__).resolve().parent.parent
    / "data"
    / "hmms"
    / "repseq_viral_core.hmm"
)

# Suffixes ``hmmpress`` writes alongside the .hmm file.
HMMPRESS_INDEX_SUFFIXES = (".h3f", ".h3i", ".h3m", ".h3p")


def resolve_database_path(user_path: Optional[str]) -> Path:
    """Return the resolved .hmm path; ``None`` → bundled.

    Raises HMMDatabaseError if the resolved path does not exist.
    """
    if user_path is None:
        path = BUNDLED_DB_PATH
        if not path.exists():
            raise HMMDatabaseError(
                f"Bundled HMM database not found at {path}. Either "
                "reinstall repseq with the bundled data files or set "
                "hmm.database to a user-supplied .hmm path in your config."
            )
        return path
    path = Path(user_path).expanduser().resolve()
    if not path.exists():
        raise HMMDatabaseError(f"hmm.database path does not exist: {path}")
    if not path.is_file():
        raise HMMDatabaseError(f"hmm.database is not a file: {path}")
    return path


def has_press_index(db_path: Path) -> bool:
    """True iff every .h3* index file exists alongside db_path."""
    return all(
        Path(str(db_path) + suffix).exists()
        for suffix in HMMPRESS_INDEX_SUFFIXES
    )


def ensure_pressed(db_path: Path) -> None:
    """Auto-run ``hmmpress`` if the .h3* indexes are missing."""
    if has_press_index(db_path):
        return
    if shutil.which("hmmpress") is None:
        raise HMMDatabaseError(
            "hmmpress is not on PATH; cannot index HMM database. Install "
            "HMMER (http://hmmer.org/) or pre-index the database manually."
        )
    proc = subprocess.run(
        ["hmmpress", "-f", str(db_path)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        raise HMMDatabaseError(
            f"hmmpress failed for {db_path} (exit {proc.returncode}):\n"
            f"{proc.stderr.strip()}"
        )


def db_signature(db_path: Path) -> str:
    """Stable cache-key suffix: sha256 of (realpath, mtime, size).

    Cheap even on multi-GB databases — no content read. Different DB
    file (or modified DB) invalidates cached hmmscan results
    automatically.
    """
    p = db_path.resolve()
    stat = p.stat()
    payload = f"{p}|{int(stat.st_mtime)}|{stat.st_size}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def parse_ga_cutoffs(db_path: Path) -> dict[str, Optional[float]]:
    """Parse Pfam GA gathering thresholds from the .hmm headers.

    Each profile may carry a ``GA <seq_ga> <dom_ga>`` line — Pfam-A
    curated profiles always do; user-built ones often don't. The
    DOMAIN GA (second number) is what we apply at ``--domtblout``
    filter time; profiles without GA map to ``None``.
    """
    cutoffs: dict[str, Optional[float]] = {}
    current_name: Optional[str] = None
    current_ga: Optional[float] = None
    try:
        with open(db_path) as fh:
            for line in fh:
                if line.startswith("NAME"):
                    if current_name is not None:
                        cutoffs[current_name] = current_ga
                    parts = line.split(None, 1)
                    current_name = parts[1].strip() if len(parts) == 2 else None
                    current_ga = None
                elif line.startswith("GA ") and current_name is not None:
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            current_ga = float(parts[2].rstrip(";"))
                        except ValueError:
                            current_ga = None
                elif line.startswith("//"):
                    if current_name is not None:
                        cutoffs[current_name] = current_ga
                        current_name = None
                        current_ga = None
        if current_name is not None:
            cutoffs[current_name] = current_ga
    except OSError as e:
        raise HMMDatabaseError(
            f"Failed to read HMM database {db_path}: {e}"
        ) from e
    return cutoffs


def profile_count(db_path: Path) -> int:
    """Count ``NAME`` records in the .hmm file (= number of profiles)."""
    count = 0
    try:
        with open(db_path) as fh:
            for line in fh:
                if line.startswith("NAME"):
                    count += 1
    except OSError:
        pass
    return count
