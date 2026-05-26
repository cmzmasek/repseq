"""Tests for repseq.hmm.database (path resolution, GA parsing, signatures)."""
from __future__ import annotations

import shutil
import time
from pathlib import Path

import pytest

from repseq.hmm.database import (
    BUNDLED_DB_PATH,
    HMMPRESS_INDEX_SUFFIXES,
    db_signature,
    ensure_pressed,
    has_press_index,
    is_press_index_stale,
    parse_ga_cutoffs,
    profile_count,
    resolve_database_path,
)
from repseq.hmm.errors import HMMDatabaseError


_MINI_HMM_SOURCE = """\
HMMER3/f [3.3 | Nov 2019]
NAME  TestProfileA
ACC   PFTEST.1
DESC  Hand-crafted test profile
LENG  10
GA    25.00 25.00
TC    27.00 27.00
NC    20.00 20.00
//
HMMER3/f [3.3 | Nov 2019]
NAME  TestProfileB
ACC   PFTEST.2
DESC  Profile without GA
LENG  10
//
"""


def _write_mini_hmm(path: Path) -> Path:
    path.write_text(_MINI_HMM_SOURCE)
    return path


def test_resolve_database_path_returns_bundled_when_user_path_is_none():
    """None → bundled. Raises if the bundled file is missing."""
    if BUNDLED_DB_PATH.exists():
        path = resolve_database_path(None)
        assert path == BUNDLED_DB_PATH
    else:
        with pytest.raises(HMMDatabaseError, match="Bundled HMM database not found"):
            resolve_database_path(None)


def test_resolve_database_path_accepts_user_path(tmp_path):
    p = _write_mini_hmm(tmp_path / "mini.hmm")
    resolved = resolve_database_path(str(p))
    assert resolved.resolve() == p.resolve()


def test_resolve_database_path_raises_for_missing(tmp_path):
    with pytest.raises(HMMDatabaseError, match="does not exist"):
        resolve_database_path(str(tmp_path / "nope.hmm"))


def test_resolve_database_path_raises_for_directory(tmp_path):
    with pytest.raises(HMMDatabaseError, match="is not a file"):
        resolve_database_path(str(tmp_path))


def test_db_signature_stable_for_same_file(tmp_path):
    p = _write_mini_hmm(tmp_path / "mini.hmm")
    sig1 = db_signature(p)
    sig2 = db_signature(p)
    assert sig1 == sig2
    # 64 hex chars from sha256
    assert len(sig1) == 64


def test_db_signature_changes_when_file_changes(tmp_path):
    p = _write_mini_hmm(tmp_path / "mini.hmm")
    sig1 = db_signature(p)
    # Mutate content + bump mtime (force, since sub-second mtime may not
    # change within the same test run on some filesystems).
    time.sleep(1.1)
    p.write_text(_MINI_HMM_SOURCE + "\n# trailing change\n")
    sig2 = db_signature(p)
    assert sig1 != sig2


def test_parse_ga_cutoffs_extracts_dom_ga(tmp_path):
    """The DOMAIN GA (second number on the GA line) is what we want."""
    p = _write_mini_hmm(tmp_path / "mini.hmm")
    ga = parse_ga_cutoffs(p)
    assert ga["TestProfileA"] == 25.00
    # Profile B has no GA line — must map to None, not raise.
    assert ga["TestProfileB"] is None


def test_profile_count_counts_NAME_records(tmp_path):
    p = _write_mini_hmm(tmp_path / "mini.hmm")
    assert profile_count(p) == 2


def test_has_press_index_detects_present_and_missing(tmp_path):
    p = _write_mini_hmm(tmp_path / "mini.hmm")
    assert has_press_index(p) is False
    # Create empty stub files for every required suffix.
    for sfx in HMMPRESS_INDEX_SUFFIXES:
        Path(str(p) + sfx).touch()
    assert has_press_index(p) is True


def test_ensure_pressed_no_op_when_indexed(tmp_path):
    p = _write_mini_hmm(tmp_path / "mini.hmm")
    for sfx in HMMPRESS_INDEX_SUFFIXES:
        Path(str(p) + sfx).touch()
    # No exception, no work; the function returns None.
    assert ensure_pressed(p) is None


def test_is_press_index_stale_returns_false_when_indexes_are_fresh(tmp_path):
    p = _write_mini_hmm(tmp_path / "mini.hmm")
    # No indexes at all → not stale (there's nothing to be stale).
    assert is_press_index_stale(p) is False
    for sfx in HMMPRESS_INDEX_SUFFIXES:
        Path(str(p) + sfx).touch()
    # Indexes created after the .hmm (or simultaneously) → not stale.
    assert is_press_index_stale(p) is False


def test_is_press_index_stale_detects_outdated_indexes(tmp_path):
    """Upgrade scenario: existing .h3* indexes from the old bundled DB,
    then the user pulls down a newer .hmm. Until repressed the indexes
    point at stale offsets — the staleness check has to catch this so
    ``ensure_pressed`` reruns hmmpress automatically.
    """
    p = _write_mini_hmm(tmp_path / "mini.hmm")
    for sfx in HMMPRESS_INDEX_SUFFIXES:
        Path(str(p) + sfx).touch()
    # Simulate the user pulling a newer .hmm onto a checkout: bump its
    # mtime past the existing indexes.
    idx_mtime = Path(str(p) + HMMPRESS_INDEX_SUFFIXES[0]).stat().st_mtime
    import os
    os.utime(p, (idx_mtime + 10, idx_mtime + 10))
    assert is_press_index_stale(p) is True


def test_ensure_pressed_represses_when_index_is_stale(tmp_path, monkeypatch):
    """``ensure_pressed`` must invoke hmmpress when indexes are stale,
    not just when they're missing — otherwise a v0.34.0 user with the
    pre-v0.34.0 indexes would silently mis-index the new bundled DB.
    Subprocess is mocked so the test stays portable.
    """
    p = _write_mini_hmm(tmp_path / "mini.hmm")
    for sfx in HMMPRESS_INDEX_SUFFIXES:
        Path(str(p) + sfx).touch()
    # Make the .hmm newer than the indexes (the upgrade scenario).
    idx_mtime = Path(str(p) + HMMPRESS_INDEX_SUFFIXES[0]).stat().st_mtime
    import os
    os.utime(p, (idx_mtime + 10, idx_mtime + 10))

    invocations: list[list[str]] = []

    class FakeProc:
        returncode = 0
        stdout = ""
        stderr = ""

    def fake_run(argv, capture_output=False, text=False):
        invocations.append(list(argv))
        return FakeProc()

    monkeypatch.setattr("repseq.hmm.database.shutil.which", lambda _: "/usr/bin/hmmpress")
    monkeypatch.setattr("repseq.hmm.database.subprocess.run", fake_run)
    ensure_pressed(p)
    assert invocations, "ensure_pressed should have invoked hmmpress on stale indexes"
    assert invocations[0][0] == "hmmpress"


@pytest.mark.skipif(shutil.which("hmmpress") is None, reason="hmmpress not on PATH")
def test_ensure_pressed_invokes_hmmpress_when_missing(tmp_path):
    """End-to-end: actually run hmmpress on a hand-built tiny .hmm.

    Skipped on environments without HMMER. Our hand-built profile won't
    pass hmmpress's checks (it's missing the transition tables), so the
    call should RAISE HMMDatabaseError — the important assertion is that
    we attempted the invocation and reported failure cleanly.
    """
    p = _write_mini_hmm(tmp_path / "mini.hmm")
    with pytest.raises(HMMDatabaseError, match="hmmpress failed"):
        ensure_pressed(p)
