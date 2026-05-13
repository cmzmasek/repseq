"""QC pipeline: duplicates, length, ambiguous chars, annotation keywords."""
from __future__ import annotations

from repseq.models import QCReport, SequenceType
from repseq.qc.pipeline import (
    ambiguous_filter,
    annotation_filter,
    length_filter,
    remove_duplicates,
    run_qc,
)


# ---------------------------------------------------------------------------
# remove_duplicates
# ---------------------------------------------------------------------------

def test_remove_duplicates_keeps_first(make_seq):
    a = make_seq("a", "ACGTACGT")
    b = make_seq("b", "ACGTACGT")  # exact dup of a
    c = make_seq("c", "TTTTGGGG")
    report = QCReport()
    kept = remove_duplicates([a, b, c], report)
    ids = [s.id for s in kept]
    assert ids == ["a", "c"]
    assert report.removed_duplicates == 1
    assert b.qc_passed is False
    assert b.qc_fail_reason and "duplicate" in b.qc_fail_reason


# ---------------------------------------------------------------------------
# length_filter
# ---------------------------------------------------------------------------

def test_length_filter_median_percent_drops_short(make_seq):
    # Median of [100,100,100,10] = 100; with 50% cutoff => min_len=50
    seqs = [
        make_seq("s1", "A" * 100),
        make_seq("s2", "A" * 100),
        make_seq("s3", "A" * 100),
        make_seq("short", "A" * 10),
    ]
    report = QCReport()
    kept = length_filter(seqs, {"mode": "median_percent", "min_percent": 50}, report)
    assert {s.id for s in kept} == {"s1", "s2", "s3"}
    assert report.removed_length == 1


def test_length_filter_min_max(make_seq):
    seqs = [
        make_seq("tiny", "A" * 5),
        make_seq("ok", "A" * 50),
        make_seq("huge", "A" * 1000),
    ]
    report = QCReport()
    kept = length_filter(
        seqs, {"mode": "min_max", "min_length": 10, "max_length": 500}, report
    )
    assert [s.id for s in kept] == ["ok"]
    assert report.removed_length == 2


def test_length_filter_empty_input_no_crash():
    report = QCReport()
    assert length_filter([], {"mode": "median_percent", "min_percent": 50}, report) == []


# ---------------------------------------------------------------------------
# ambiguous_filter
# ---------------------------------------------------------------------------

def test_ambiguous_filter_drops_high_ambiguous(make_seq):
    clean = make_seq("clean", "ACGTACGTACGT", seq_type=SequenceType.NUCLEOTIDE)
    dirty = make_seq("dirty", "ACGTNNNNNNNN", seq_type=SequenceType.NUCLEOTIDE)
    report = QCReport()
    kept = ambiguous_filter([clean, dirty], threshold=0.2, report=report)
    assert [s.id for s in kept] == ["clean"]
    assert report.removed_ambiguous == 1


def test_ambiguous_filter_protein_X_is_ambiguous(make_seq):
    s = make_seq("p", "MEEPXXXXX" + "MEEPQSDPSVE", seq_type=SequenceType.PROTEIN)
    report = QCReport()
    kept = ambiguous_filter([s], threshold=0.05, report=report)
    assert kept == []


# ---------------------------------------------------------------------------
# annotation_filter
# ---------------------------------------------------------------------------

def test_annotation_filter_blocks_keyword(make_seq):
    bad = make_seq("p1", "ACGT", header="P1 hypothetical protein [foo]")
    ok = make_seq("p2", "ACGT", header="P2 RNA polymerase [foo]")
    report = QCReport()
    cfg = {"enabled": True, "keywords": ["hypothetical", "MAG:"]}
    kept = annotation_filter([bad, ok], cfg, report)
    assert [s.id for s in kept] == ["p2"]
    assert report.removed_annotation == 1


def test_annotation_filter_disabled_passthrough(make_seq):
    bad = make_seq("p1", "ACGT", header="P1 hypothetical")
    report = QCReport()
    cfg = {"enabled": False, "keywords": ["hypothetical"]}
    kept = annotation_filter([bad], cfg, report)
    assert kept == [bad]
    assert report.removed_annotation == 0


# ---------------------------------------------------------------------------
# end-to-end pipeline
# ---------------------------------------------------------------------------

def test_run_qc_pipeline(make_seq):
    # NB: each sequence has unique content so the duplicate filter only catches
    # `dup` (identical to `ok`), letting `hyp` proceed to the annotation step.
    seqs = [
        make_seq("ok", "ACGT" * 50, seq_type=SequenceType.NUCLEOTIDE,
                 header="ok normal sequence"),
        make_seq("dup", "ACGT" * 50, seq_type=SequenceType.NUCLEOTIDE,
                 header="dup duplicate of ok"),
        make_seq("short", "ACGT", seq_type=SequenceType.NUCLEOTIDE,
                 header="short tiny"),
        make_seq("hyp", "TGCA" * 50, seq_type=SequenceType.NUCLEOTIDE,
                 header="hyp hypothetical protein"),
    ]
    cfg = {
        "qc": {
            "remove_duplicates": True,
            "length_filter": {"mode": "min_max", "min_length": 50},
            "ambiguous_threshold": 0.05,
            "annotation_filter": {"enabled": True, "keywords": ["hypothetical"]},
        }
    }
    kept, report = run_qc(seqs, cfg)
    assert [s.id for s in kept] == ["ok"]
    assert report.total_input == 4
    assert report.passed == 1
    assert report.removed_duplicates == 1
    assert report.removed_length == 1
    assert report.removed_annotation == 1
