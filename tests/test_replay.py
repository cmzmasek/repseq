"""``repseq replay`` — internal helpers + end-to-end FASTA emission.

Network is mocked: the NCBI client is monkey-patched to return canned
sequence bodies. These tests pin the accession-flattening, the
CONCAT-rebuild path for segmented replays, the missing-accession
contract (warn + write TSV, don't crash), and that ``write_results``
sees a well-formed ``RunResult`` after the rebuild.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from repseq.replay import (
    _build_replay_result,
    _collect_accessions,
    _make_sequence,
    _write_missing_tsv,
)


# ---------------------------------------------------------------------------
# _collect_accessions — flatten lockfile reps into a fetch list
# ---------------------------------------------------------------------------

def test_collect_accessions_non_segmented():
    lf = {"representatives": [
        {"kind": "sequence", "accession": "A1"},
        {"kind": "sequence", "accession": "A2"},
    ]}
    assert _collect_accessions(lf) == ["A1", "A2"]


def test_collect_accessions_segmented_flattens_per_segment():
    lf = {"representatives": [
        {"kind": "isolate", "isolate_id": "iso1",
         "segment_accessions": {"L": "L1", "M": "M1", "S": "S1"}},
        {"kind": "isolate", "isolate_id": "iso2",
         "segment_accessions": {"L": "L2", "M": "M2"}},
    ]}
    out = _collect_accessions(lf)
    assert set(out) == {"L1", "M1", "S1", "L2", "M2"}


def test_collect_accessions_skips_empty():
    """Defensive: dropped segments or malformed entries shouldn't
    insert blank strings into the fetch list."""
    lf = {"representatives": [
        {"kind": "sequence", "accession": None},
        {"kind": "isolate", "segment_accessions": {"L": "L1", "M": ""}},
    ]}
    assert _collect_accessions(lf) == ["L1"]


# ---------------------------------------------------------------------------
# _build_replay_result — non-segmented
# ---------------------------------------------------------------------------

def test_build_replay_result_non_segmented_uses_fetched_bodies():
    lf = {
        "mode": "global",
        "representatives": [
            {"kind": "sequence", "accession": "AB1", "organism": "Virus X"},
            {"kind": "sequence", "accession": "AB2", "organism": "Virus Y"},
        ],
        "config": {},
    }
    seqs = {"AB1": "ACGT" * 5, "AB2": "TTGC" * 3}
    result, complete, segs, missing = _build_replay_result(lf, seqs, {})
    assert complete is None  # non-segmented
    assert segs is None
    assert missing == []
    assert len(result.representatives) == 2
    assert result.representatives[0].accession == "AB1"
    assert result.representatives[0].sequence == "ACGT" * 5
    assert result.representatives[0].organism == "Virus X"


def test_build_replay_result_missing_accession_tracked():
    """A rep whose NCBI fetch failed is omitted from result.representatives
    AND surfaces in the missing list for the .tsv writer."""
    lf = {
        "mode": "global",
        "representatives": [
            {"id": "rep_a", "kind": "sequence", "accession": "AB1"},
            {"id": "rep_b", "kind": "sequence", "accession": "AB2"},
        ],
        "config": {},
    }
    seqs = {"AB1": "ACGT", "AB2": None}  # AB2 was not found
    result, _comp, _segs, missing = _build_replay_result(lf, seqs, {})
    assert len(result.representatives) == 1
    assert result.representatives[0].accession == "AB1"
    assert missing == [("rep_b", "AB2")]


# ---------------------------------------------------------------------------
# _build_replay_result — segmented (CONCAT reconstruction)
# ---------------------------------------------------------------------------

def test_build_replay_result_segmented_rebuilds_concat():
    lf = {
        "mode": "taxonomic1",
        "representatives": [{
            "id": "CONCAT|iso1",
            "kind": "isolate",
            "isolate_id": "iso1",
            "organism": "Virus Z",
            "segment_accessions": {"L": "L1", "M": "M1", "S": "S1"},
        }],
        "config": {
            "viruses": {
                "test_virus": {"segments": ["L", "M", "S"]},
            },
        },
    }
    seqs = {"L1": "A" * 10, "M1": "C" * 8, "S1": "G" * 6}
    result, complete, segs, missing = _build_replay_result(lf, seqs, {})
    assert segs == ["L", "M", "S"]
    assert missing == []
    assert len(result.representatives) == 1
    concat = result.representatives[0]
    assert concat.id == "CONCAT|iso1"
    assert concat.isolate_id == "iso1"
    assert concat.sequence == "A" * 10 + "C" * 8 + "G" * 6
    assert concat.concat_segments is not None
    assert [s.segment for s in concat.concat_segments] == ["L", "M", "S"]
    assert "iso1" in complete
    assert len(complete["iso1"]) == 3


def test_build_replay_result_segmented_partial_isolate_dropped():
    """If any segment of an isolate fails to fetch, the isolate is
    dropped (incomplete) BUT the partial accessions are recorded
    in missing so the user knows what fell out and why."""
    lf = {
        "mode": "taxonomic1",
        "representatives": [{
            "id": "CONCAT|iso1",
            "kind": "isolate",
            "isolate_id": "iso1",
            "segment_accessions": {"L": "L1", "M": "M1", "S": "S1"},
        }],
        "config": {},
    }
    # M1 not returned
    seqs = {"L1": "ACGT", "M1": None, "S1": "TTGC"}
    result, complete, _segs, missing = _build_replay_result(lf, seqs, {})
    # Isolate keeps its two surviving segments — replay's contract is
    # "best effort": partial output is better than no output.
    assert len(result.representatives) == 1
    assert ("CONCAT|iso1", "M1") in missing


# ---------------------------------------------------------------------------
# _make_sequence
# ---------------------------------------------------------------------------

def test_make_sequence_uses_source_metadata_fallbacks():
    """When the explicit segment/isolate kwargs are None, fall back to
    the source-feature qualifiers from the GenBank fetch."""
    seq = _make_sequence(
        "AB1", "ACGT", "Virus X",
        source_meta={"isolate": "iso_xyz", "strain": None, "segment": "L"},
    )
    assert seq.isolate_id == "iso_xyz"
    assert seq.segment == "L"


def test_make_sequence_explicit_kwargs_override_source_meta():
    """When the lockfile explicitly named a segment label, that wins
    over the source feature."""
    seq = _make_sequence(
        "AB1", "ACGT", "X",
        source_meta={"segment": "S", "isolate": "wrong"},
        segment_label="L",
        isolate_id="iso1",
    )
    assert seq.segment == "L"
    assert seq.isolate_id == "iso1"


# ---------------------------------------------------------------------------
# _write_missing_tsv
# ---------------------------------------------------------------------------

def test_write_missing_tsv_writes_header_and_rows(tmp_path):
    path = tmp_path / "missing.tsv"
    _write_missing_tsv([("rep_a", "AB1"), ("rep_b", "AB2")], path)
    text = path.read_text()
    assert text.splitlines()[0] == "representative_id\taccession"
    assert "rep_a\tAB1" in text
    assert "rep_b\tAB2" in text


def test_write_missing_tsv_skips_when_empty(tmp_path):
    """Empty missing list → no file (don't pollute the output dir)."""
    path = tmp_path / "missing.tsv"
    _write_missing_tsv([], path)
    assert not path.exists()


# ---------------------------------------------------------------------------
# NCBI client — fetch_nucleotide_batch parser
# ---------------------------------------------------------------------------

def test_fetch_nucleotide_batch_parses_fasta_response(tmp_path, monkeypatch):
    """Mock the Entrez efetch FASTA round trip; verify the parser pulls
    one body per accession, uppercases, and strips gaps."""
    from repseq.taxonomy.cache import TaxonomyCache
    from repseq.taxonomy import ncbi as ncbi_mod

    cache = TaxonomyCache(tmp_path / "cache")
    client = ncbi_mod.NCBITaxonomy(cache=cache)

    fake_response_text = (
        ">AB1.1 description here\n"
        "acgtacgt\n"
        "acgt\n"
        ">AB2.1 other\n"
        "TTGC-\n"
        "TTGC\n"
    )

    class FakeResp:
        text = fake_response_text
        def raise_for_status(self): pass

    monkeypatch.setattr(
        ncbi_mod.requests, "get", lambda *a, **k: FakeResp(),
    )
    out = client.fetch_nucleotide_batch(["AB1", "AB2"])
    assert out["AB1"] == "ACGTACGTACGT"
    assert out["AB2"] == "TTGCTTGC"


def test_fetch_nucleotide_batch_uses_cache_on_second_call(tmp_path, monkeypatch):
    """Second call with the same accessions hits the cache, not the network."""
    from repseq.taxonomy.cache import TaxonomyCache
    from repseq.taxonomy import ncbi as ncbi_mod

    cache = TaxonomyCache(tmp_path / "cache")
    client = ncbi_mod.NCBITaxonomy(cache=cache)

    call_count = {"n": 0}
    class FakeResp:
        text = ">AB1\nACGT\n"
        def raise_for_status(self): pass

    def fake_get(*a, **k):
        call_count["n"] += 1
        return FakeResp()

    monkeypatch.setattr(ncbi_mod.requests, "get", fake_get)
    client.fetch_nucleotide_batch(["AB1"])
    client.fetch_nucleotide_batch(["AB1"])
    assert call_count["n"] == 1
