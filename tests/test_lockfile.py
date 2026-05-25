"""Lockfile assembly + schema-versioning round-trip.

These tests don't touch the network. They exercise the schema, the
top-level keys, segmented vs non-segmented representative
serialisation, and the read-back validation contract — all parts the
replay subcommand depends on.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from repseq import __version__ as REPSEQ_VERSION
from repseq.lockfile import (
    LockfileVersionError,
    SCHEMA_VERSION,
    build_lockfile,
    compute_sha256,
    read_lockfile,
    write_lockfile,
)
from repseq.models import Cluster, RunResult, Sequence, SequenceSource, SequenceType


def _seq(sid: str, body: str, **kwargs) -> Sequence:
    return Sequence(
        id=sid, header=sid, sequence=body,
        seq_type=SequenceType.NUCLEOTIDE, source=SequenceSource.NCBI,
        accession=sid, **kwargs,
    )


def _result_non_seg() -> RunResult:
    reps = [_seq("AB1234", "ACGT" * 10, organism="Virus A")]
    return RunResult(
        mode="global",
        representatives=reps,
        clusters=[Cluster(cluster_id="c0", representative=reps[0], members=[])],
    )


def _result_segmented() -> RunResult:
    s1 = _seq("L01", "ACGT" * 10, segment="L", isolate_id="iso1")
    s2 = _seq("M01", "ACGT" * 8, segment="M", isolate_id="iso1")
    s3 = _seq("S01", "ACGT" * 6, segment="S", isolate_id="iso1")
    concat = Sequence(
        id="CONCAT|iso1", header="iso1",
        sequence=(s1.sequence + s2.sequence + s3.sequence),
        seq_type=SequenceType.NUCLEOTIDE,
        source=SequenceSource.NCBI,
        organism="Virus B", isolate_id="iso1",
    )
    concat.concat_segments = [s1, s2, s3]
    return RunResult(
        mode="taxonomic1",
        representatives=[concat],
        clusters=[Cluster(cluster_id="c0", representative=concat, members=[])],
    )


# ---------------------------------------------------------------------------
# build_lockfile — top-level shape
# ---------------------------------------------------------------------------

def test_build_lockfile_top_level_keys(tmp_path):
    cfg = {"output": {"dir": str(tmp_path), "prefix": "test"}, "segmented": {"enabled": False}}
    lf = build_lockfile(cfg, _result_non_seg(), ["/path/to/in.fasta"], command="repseq global ...")
    assert lf["schema_version"] == SCHEMA_VERSION
    assert lf["repseq_version"] == REPSEQ_VERSION
    assert lf["mode"] == "global"
    assert lf["command"] == "repseq global ..."
    assert "created_utc" in lf
    assert "python_version" in lf
    assert "config" in lf
    assert "inputs" in lf
    assert "tools" in lf
    assert "representatives" in lf


def test_build_lockfile_records_inputs(tmp_path):
    """Input list carries one entry per input file, with sha256
    computed against the on-disk contents."""
    fa = tmp_path / "in.fasta"
    fa.write_text(">seq1\nACGT\n")
    cfg = {"output": {"dir": str(tmp_path), "prefix": "test"}, "segmented": {"enabled": False}}
    lf = build_lockfile(cfg, _result_non_seg(), [str(fa)], command="x")
    assert len(lf["inputs"]) == 1
    assert lf["inputs"][0]["path"] == str(fa)
    assert lf["inputs"][0]["sha256"] == compute_sha256(fa)
    assert lf["inputs"][0]["size_bytes"] == fa.stat().st_size


def test_build_lockfile_representatives_non_segmented(tmp_path):
    """Non-segmented reps carry kind=sequence + a single accession."""
    cfg = {"output": {"dir": str(tmp_path)}, "segmented": {"enabled": False}}
    lf = build_lockfile(cfg, _result_non_seg(), [], command="x")
    reps = lf["representatives"]
    assert len(reps) == 1
    assert reps[0]["kind"] == "sequence"
    assert reps[0]["accession"] == "AB1234"
    assert reps[0]["organism"] == "Virus A"


def test_build_lockfile_representatives_segmented(tmp_path):
    """Segmented reps carry kind=isolate + a segment→accession map
    (so a replay can re-fetch all three segments and rebuild CONCAT)."""
    cfg = {"output": {"dir": str(tmp_path)}, "segmented": {"enabled": True}}
    lf = build_lockfile(cfg, _result_segmented(), [], command="x")
    reps = lf["representatives"]
    assert len(reps) == 1
    rep = reps[0]
    assert rep["kind"] == "isolate"
    assert rep["isolate_id"] == "iso1"
    assert rep["segment_accessions"] == {"L": "L01", "M": "M01", "S": "S01"}


# ---------------------------------------------------------------------------
# write_lockfile / read_lockfile round-trip
# ---------------------------------------------------------------------------

def test_write_lockfile_produces_stable_sorted_json(tmp_path):
    cfg = {"output": {"dir": str(tmp_path)}, "segmented": {"enabled": False}}
    lf = build_lockfile(cfg, _result_non_seg(), [], command="x")
    path = write_lockfile(lf, tmp_path / "out.json")
    text = path.read_text()
    # Sorted-key invariant: 'config' comes before 'created_utc' before
    # 'inputs' before 'mode' before 'repseq_version' …
    config_pos = text.index('"config"')
    created_pos = text.index('"created_utc"')
    inputs_pos = text.index('"inputs"')
    assert config_pos < created_pos < inputs_pos


def test_read_lockfile_roundtrip(tmp_path):
    cfg = {"output": {"dir": str(tmp_path)}, "segmented": {"enabled": False}}
    lf = build_lockfile(cfg, _result_non_seg(), [], command="x")
    path = write_lockfile(lf, tmp_path / "out.json")
    loaded = read_lockfile(path)
    assert loaded["schema_version"] == SCHEMA_VERSION
    assert loaded["representatives"] == lf["representatives"]


def test_read_lockfile_rejects_future_schema(tmp_path):
    """A lockfile from a future major schema is rejected — we'd
    misread fields we don't know about."""
    path = tmp_path / "bad.json"
    payload = {"schema_version": SCHEMA_VERSION + 1, "repseq_version": "1.0.0"}
    path.write_text(json.dumps(payload))
    with pytest.raises(LockfileVersionError, match="newer than this"):
        read_lockfile(path)


def test_read_lockfile_rejects_non_object(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text("[]")
    with pytest.raises(LockfileVersionError, match="must be an object"):
        read_lockfile(path)


def test_read_lockfile_rejects_missing_schema_version(tmp_path):
    path = tmp_path / "bad.json"
    path.write_text('{"repseq_version": "1.0.0"}')
    with pytest.raises(LockfileVersionError, match="schema_version"):
        read_lockfile(path)


# ---------------------------------------------------------------------------
# compute_sha256
# ---------------------------------------------------------------------------

def test_compute_sha256_known_value(tmp_path):
    """Verify the SHA256 against a known reference for 'hello\\n'."""
    p = tmp_path / "hello.txt"
    p.write_bytes(b"hello\n")
    # Known: sha256('hello\n') = 5891b5b522d5df086d0ff0b110fbd9d21bb4fc7163af34d08286a2e846f6be03
    assert compute_sha256(p) == (
        "5891b5b522d5df086d0ff0b110fbd9d21bb4fc7163af34d08286a2e846f6be03"
    )


def test_compute_sha256_missing_file_returns_empty(tmp_path):
    """Soft-fail for missing files so the lockfile writer doesn't
    crash when an input path was deleted between read and lockfile time."""
    assert compute_sha256(tmp_path / "nope.txt") == ""
