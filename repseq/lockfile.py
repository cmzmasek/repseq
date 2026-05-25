"""``{prefix}_lockfile.json`` — the structured reproducibility record.

Where ``{prefix}_summary.md`` records *what was done* in prose for a
methods section, the lockfile records *what came out* in machine-
readable form: the elected representative IDs + their parent
accessions, the post-mutation configuration the pipeline ran under,
the version of every external tool that ran, and a SHA-256
fingerprint of the HMM database. Together with the input FASTA's
sha256, this is enough information to re-fetch the same sequences
months later and re-emit byte-identical representative FASTAs via
``repseq replay``.

Design decisions (locked in v0.30.0):

* **JSON only.** Sequences are not embedded — they're re-fetched via
  the existing NCBI/UniProt cache path, so the lockfile stays small
  and human-readable.
* **Representatives only.** Cluster membership is not pinned;
  ``{prefix}_clusters.tsv`` carries that already, and recording it
  twice would let the two go out of sync.
* **Schema-versioned.** ``schema_version: 1`` is at the top so
  future format changes can hard-fail an old replay rather than
  produce silently-wrong outputs.
* **Sorted-key JSON.** Stable byte output makes the file diff-
  friendly between runs (the only field that changes for an
  identical pipeline rerun is ``created_utc``).
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from . import __version__ as REPSEQ_VERSION
from .models import RunResult, Sequence

# Lockfile JSON shape version. Bump on any breaking change to the
# top-level keys or to a representative entry's shape.
SCHEMA_VERSION = 1


def compute_sha256(path: Path) -> str:
    """Hex sha256 of the file at ``path``, streamed in 1 MiB blocks.

    Used for both input FASTAs and the HMM database. Returns an empty
    string if the file can't be read — the caller decides whether a
    missing fingerprint is fatal (HMM DB → warn loudly at replay time;
    input FASTA → just records absence).
    """
    h = hashlib.sha256()
    try:
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(chunk)
    except OSError:
        return ""
    return h.hexdigest()


def _input_sha256_list(input_paths: list[str]) -> list[dict[str, Any]]:
    """``[{path, sha256, size_bytes}, …]`` for every input FASTA.

    A list rather than a dict so the order — which can affect QC
    counts when duplicates straddle files — is preserved.
    """
    out: list[dict[str, Any]] = []
    for p in input_paths:
        path = Path(p)
        size = path.stat().st_size if path.exists() else None
        out.append({
            "path": str(p),
            "sha256": compute_sha256(path),
            "size_bytes": size,
        })
    return out


def _hmm_db_fingerprint(cfg: dict[str, Any]) -> Optional[dict[str, Any]]:
    """Return the HMM database fingerprint, or ``None`` if the HMM tier
    wasn't configured at all.

    Records:
      * ``path`` — the resolved absolute path
      * ``sha256`` — full content hash (not the cheap mtime/size
        signature ``db_signature`` uses for cache-keying, because
        replay needs to detect *any* content change including a
        retread of the same path)
      * ``bundled`` — true when the bundled db_dir copy was used
      * ``n_profiles`` — count of ``NAME`` records, sanity-check
    """
    hmm_cfg = (cfg or {}).get("hmm", {}) or {}
    # The HMM tier may exist in cfg even when nothing in this run
    # invoked it; we still record the fingerprint because a future
    # replay reproducing the SAME run might.
    try:
        from .hmm.database import (
            BUNDLED_DB_PATH,
            profile_count,
            resolve_database_path,
        )
    except ImportError:
        return None
    user_path = hmm_cfg.get("database")
    try:
        db_path = resolve_database_path(user_path)
    except Exception:
        return None
    return {
        "path": str(db_path),
        "sha256": compute_sha256(db_path),
        "bundled": db_path == BUNDLED_DB_PATH,
        "n_profiles": profile_count(db_path),
    }


def _representative_entry(rep: Sequence) -> dict[str, Any]:
    """Serialise one representative to the lockfile schema.

    Segmented reps are synthetic ``CONCAT|<isolate_id>`` objects
    carrying ``concat_segments`` — we record the per-segment
    accessions so a replay can refetch each segment and rebuild the
    CONCAT. Non-segmented reps just carry an accession.
    """
    if rep.concat_segments:
        seg_accessions: dict[str, Optional[str]] = {}
        for seg in rep.concat_segments:
            if seg.segment:
                seg_accessions[seg.segment] = seg.accession
        return {
            "id": rep.id,
            "kind": "isolate",
            "isolate_id": rep.isolate_id,
            "organism": rep.organism,
            "segment_accessions": seg_accessions,
        }
    return {
        "id": rep.id,
        "kind": "sequence",
        "accession": rep.accession,
        "organism": rep.organism,
    }


def build_lockfile(
    cfg: dict[str, Any],
    result: RunResult,
    input_paths: list[str],
    command: str,
) -> dict[str, Any]:
    """Assemble the lockfile dict.

    The caller is expected to have already finished the run; ``cfg``
    must be the post-mutation cfg (after ``--fast`` overrides,
    ``--alphabet-for-clustering``, etc.) so replay reproduces what
    actually ran rather than what the YAML said.
    """
    # Tool versions: reuse the same probes the summary renderer uses,
    # imported lazily so a render error in the summary can't break
    # lockfile emission.
    try:
        from .output.summary import detect_tool_versions
        tool_versions = detect_tool_versions()
    except Exception:
        tool_versions = {}

    return {
        "schema_version": SCHEMA_VERSION,
        "repseq_version": REPSEQ_VERSION,
        "created_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "mode": result.mode,
        "command": command,
        "python_version": (
            f"{sys.version_info.major}.{sys.version_info.minor}."
            f"{sys.version_info.micro}"
        ),
        "config": cfg,
        "inputs": _input_sha256_list(input_paths),
        "tools": dict(tool_versions),
        "hmm_db": _hmm_db_fingerprint(cfg),
        "representatives": [
            _representative_entry(rep) for rep in result.representatives
        ],
    }


def write_lockfile(lockfile: dict[str, Any], path: Path) -> Path:
    """Serialise the lockfile dict to ``path``.

    Sorted-key JSON with 2-space indentation — stable byte output so
    two lockfiles from runs that differ only in time-of-day diff
    cleanly on ``created_utc`` and nothing else.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(lockfile, sort_keys=True, indent=2, default=str) + "\n",
    )
    return path


class LockfileVersionError(RuntimeError):
    """Raised when a lockfile's ``schema_version`` is incompatible
    with the running repseq's expected version. Caller surfaces as
    actionable error rather than silently misreading old/new fields.
    """


def read_lockfile(path: Path) -> dict[str, Any]:
    """Load and validate a lockfile.

    Hard-fails on missing/unparseable JSON or a ``schema_version``
    higher than this repseq supports (older repseq + newer lockfile —
    we'd misread fields we don't know about). A lower schema version
    is loaded as-is; downstream callers must handle missing optional
    keys defensively.
    """
    with open(path) as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise LockfileVersionError(
            f"{path}: top-level JSON must be an object, got {type(data).__name__}"
        )
    schema = data.get("schema_version")
    if not isinstance(schema, int):
        raise LockfileVersionError(
            f"{path}: missing or non-integer schema_version"
        )
    if schema > SCHEMA_VERSION:
        raise LockfileVersionError(
            f"{path}: schema_version {schema} is newer than this "
            f"repseq supports (max {SCHEMA_VERSION}). Upgrade repseq "
            "to read this lockfile."
        )
    return data
