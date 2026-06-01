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
import tempfile
from pathlib import Path
from typing import Optional

from .errors import HMMDatabaseError

# Bundled HMM sets ship under repseq/data/hmms/. Resolved relative to this
# file so editable installs (``pip install -e .``) work without a copy.
# Family-specific sets live in subdirectories (e.g.
# repseq/data/hmms/Filoviridae/), selectable by bare name in
# ``hmm.database`` (see resolve_database_path).
BUNDLED_HMMS_DIR = Path(__file__).resolve().parent.parent / "data" / "hmms"
BUNDLED_DB_PATH = BUNDLED_HMMS_DIR / "repseq_viral_core.hmm"

# Suffixes ``hmmpress`` writes alongside the .hmm file.
HMMPRESS_INDEX_SUFFIXES = (".h3f", ".h3i", ".h3m", ".h3p")


def resolve_database_path(
    user_path: Optional[str], *, cache_dir: Optional[Path] = None
) -> Path:
    """Return a single .hmm file path to scan against; ``None`` → bundled.

    ``hmm.database`` may be:

    * ``None`` — the bundled ``repseq_viral_core.hmm``.
    * a **file** path — used as-is (the classic single-database form).
    * a **directory** path — every ``*.hmm`` in it is concatenated into one
      combined database (cached under ``cache_dir``), so a user can keep
      family-specific profiles as separate files. The combined file is what
      gets pressed / scanned / GA-parsed downstream.
    * a **bare name** that doesn't resolve as a path — looked up under the
      bundled ``repseq/data/hmms/`` directory, so ``hmm.database:
      Filoviridae`` selects the bundled family set
      ``repseq/data/hmms/Filoviridae/`` (a directory, combined as above).

    Raises HMMDatabaseError if nothing resolves. ``cache_dir`` is where the
    combined directory database is written (defaults to a temp dir); only
    consulted for the directory case.
    """
    if user_path is None:
        if not BUNDLED_DB_PATH.exists():
            raise HMMDatabaseError(
                f"Bundled HMM database not found at {BUNDLED_DB_PATH}. Either "
                "reinstall repseq with the bundled data files or set "
                "hmm.database to a user-supplied .hmm path in your config."
            )
        return BUNDLED_DB_PATH

    path = Path(user_path).expanduser()
    if not path.exists():
        # Fall back to a bundled set selected by bare name
        # (e.g. "Filoviridae" → repseq/data/hmms/Filoviridae).
        bundled = BUNDLED_HMMS_DIR / user_path
        if bundled.exists():
            path = bundled
        else:
            raise HMMDatabaseError(
                f"hmm.database does not exist as a path ({path}) and is not a "
                f"bundled set under {BUNDLED_HMMS_DIR} (looked for "
                f"'{user_path}'). Point it at a .hmm file, a directory of "
                f".hmm files, or a bundled set name."
            )
    path = path.resolve()
    if path.is_dir():
        return _combine_hmm_directory(path, cache_dir)
    if not path.is_file():
        raise HMMDatabaseError(f"hmm.database is not a file or directory: {path}")
    return path


def _combine_hmm_directory(
    dir_path: Path, cache_dir: Optional[Path]
) -> Path:
    """Concatenate every ``*.hmm`` in ``dir_path`` into one combined database.

    The combined file is cached under ``cache_dir`` (or a temp dir) keyed by
    a content signature of the member files, so it is rebuilt only when a
    member is added, removed, or modified — and the same combined file (and
    its pressed indexes) is reused across runs. Members are concatenated in
    sorted filename order for a deterministic result.
    """
    members = sorted(p for p in dir_path.glob("*.hmm") if p.is_file())
    if not members:
        raise HMMDatabaseError(
            f"HMM database directory {dir_path} contains no .hmm files."
        )
    # Content signature: name + size + mtime of every member (cheap, no read).
    sig_payload = "|".join(
        f"{p.name}:{p.stat().st_size}:{int(p.stat().st_mtime)}" for p in members
    ).encode("utf-8")
    sig = hashlib.sha256(sig_payload).hexdigest()[:16]

    base = Path(cache_dir).expanduser() if cache_dir else Path(tempfile.gettempdir())
    combined_dir = base / "hmm_combined"
    combined_dir.mkdir(parents=True, exist_ok=True)
    combined = combined_dir / f"{dir_path.name}_{sig}.hmm"

    if not combined.exists():
        tmp = combined.with_suffix(".hmm.partial")
        with open(tmp, "w") as out:
            for member in members:
                text = member.read_text()
                out.write(text)
                if not text.endswith("\n"):
                    out.write("\n")
        tmp.replace(combined)  # atomic publish so a half-written file is never used
    return combined


def has_press_index(db_path: Path) -> bool:
    """True iff every .h3* index file exists alongside db_path."""
    return all(
        Path(str(db_path) + suffix).exists()
        for suffix in HMMPRESS_INDEX_SUFFIXES
    )


def is_press_index_stale(db_path: Path) -> bool:
    """True iff any .h3* index is older than the .hmm source file.

    Catches the upgrade scenario: a user pulls down a newer bundled
    ``.hmm`` (e.g. v0.34.0 added the Coronaviridae replicase profiles)
    while their existing ``.h3*`` files were pressed against the old
    content. The stale binary indexes would silently mis-index — hmmscan
    either fails to find the new HMMs by name or reads corrupted offsets.
    """
    if not has_press_index(db_path):
        return False
    hmm_mtime = db_path.stat().st_mtime
    for suffix in HMMPRESS_INDEX_SUFFIXES:
        if Path(str(db_path) + suffix).stat().st_mtime < hmm_mtime:
            return True
    return False


def ensure_pressed(db_path: Path) -> None:
    """Auto-run ``hmmpress`` if the .h3* indexes are missing OR stale."""
    if has_press_index(db_path) and not is_press_index_stale(db_path):
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
