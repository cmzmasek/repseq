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


def test_remove_duplicates_keeps_refseq_over_earlier_plain(make_seq):
    # Byte-identical records: the curated RefSeq copy must survive even
    # though a plain duplicate appeared first in the file.
    plain = make_seq("plain", "ACGTACGT")
    refseq = make_seq("refseq", "ACGTACGT", is_refseq=True)
    report = QCReport()
    kept = remove_duplicates([plain, refseq], report)
    assert [s.id for s in kept] == ["refseq"]
    assert plain.qc_passed is False
    assert plain.qc_fail_reason == "exact_duplicate_of:refseq"


def test_remove_duplicates_prefers_reviewed_then_first_seen(make_seq):
    a = make_seq("a", "TTTTGGGG")
    b = make_seq("b", "TTTTGGGG", is_reviewed=True)
    c = make_seq("c", "TTTTGGGG")
    report = QCReport()
    kept = remove_duplicates([a, b, c], report)
    assert [s.id for s in kept] == ["b"]
    assert report.removed_duplicates == 2


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


def test_length_filter_median_percent_caps_long(make_seq):
    # max_percent should also drop oversized records (median = 100).
    seqs = [
        make_seq("s1", "A" * 100),
        make_seq("s2", "A" * 100),
        make_seq("s3", "A" * 100),
        make_seq("big", "A" * 500),
    ]
    report = QCReport()
    cfg = {"mode": "median_percent", "min_percent": 50, "max_percent": 150}
    kept = length_filter(seqs, cfg, report)
    assert {s.id for s in kept} == {"s1", "s2", "s3"}
    assert report.removed_length == 1


def test_length_filter_min_percent_zero_disables_lower_bound(make_seq):
    seqs = [make_seq("s1", "A" * 100), make_seq("tiny", "A" * 2)]
    report = QCReport()
    kept = length_filter(seqs, {"mode": "median_percent", "min_percent": 0}, report)
    assert {s.id for s in kept} == {"s1", "tiny"}
    assert report.removed_length == 0


def test_run_qc_skips_length_filter_in_segmented_mode(make_seq):
    # Regression: in segmented mode the input is a mix of long and short
    # segments. A whole-pool median would drop the short segments, leaving
    # every isolate incomplete. run_qc must skip the global length filter so
    # per-segment segment_lengths bounds can do the job downstream.
    seqs = [
        make_seq("L1", "A" * 6000),
        make_seq("M1", "A" * 4000),
        make_seq("S1", "A" * 900),   # < 50% of median — would be dropped
    ]
    cfg = {"qc": {"length_filter": {"mode": "median_percent", "min_percent": 50}},
           "segmented": {"enabled": True}}
    kept, report = run_qc(seqs, cfg)
    assert {s.id for s in kept} == {"L1", "M1", "S1"}
    assert report.removed_length == 0
    assert report.length_filter_skipped is True


def test_run_qc_skips_dedup_in_segmented_mode(make_seq):
    # Regression: a segment can be byte-identical between two otherwise
    # distinct isolates (a conserved segment). If run_qc deduplicated the
    # segment pool it would drop one copy, and the completeness step would
    # then discard a whole legitimate isolate as "missing" that segment.
    # Dedup must be skipped here and applied to concatenated isolates instead.
    seqs = [
        make_seq("a_S", "ACGTACGTACGT", segment="S"),
        make_seq("b_S", "ACGTACGTACGT", segment="S"),  # identical S, different isolate
        make_seq("a_L", "TTTTGGGGTTTT", segment="L"),
        make_seq("b_L", "CCCCAAAACCCC", segment="L"),
    ]
    cfg = {"qc": {}, "segmented": {"enabled": True}}
    kept, report = run_qc(seqs, cfg)
    assert {s.id for s in kept} == {"a_S", "b_S", "a_L", "b_L"}
    assert report.removed_duplicates == 0
    assert report.dedup_skipped is True


def test_run_qc_applies_dedup_when_not_segmented(make_seq):
    seqs = [
        make_seq("a", "ACGTACGTACGT"),
        make_seq("b", "ACGTACGTACGT"),  # exact dup
        make_seq("c", "TTTTGGGGTTTT"),
    ]
    cfg = {"qc": {}, "segmented": {"enabled": False}}
    kept, report = run_qc(seqs, cfg)
    assert {s.id for s in kept} == {"a", "c"}
    assert report.removed_duplicates == 1
    assert report.dedup_skipped is False


def test_run_qc_applies_length_filter_when_not_segmented(make_seq):
    seqs = [
        make_seq("L1", "A" * 6000),
        make_seq("M1", "A" * 4000),
        make_seq("S1", "A" * 900),
    ]
    cfg = {"qc": {"length_filter": {"mode": "median_percent", "min_percent": 50}},
           "segmented": {"enabled": False}}
    kept, report = run_qc(seqs, cfg)
    assert "S1" not in {s.id for s in kept}
    assert report.removed_length == 1
    assert report.length_filter_skipped is False


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


def test_ambiguous_fraction_protein_excludes_real_residues(make_seq):
    # U (selenocysteine) and O (pyrrolysine) are definite residues, not
    # ambiguity codes — only X/B/Z/J should count toward the fraction.
    s = make_seq("p", "MKUXOJ", seq_type=SequenceType.PROTEIN)
    # ambiguous chars in "MKUXOJ": X and J  => 2 of 6
    assert abs(s.ambiguous_fraction - 2 / 6) < 1e-9


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


def test_annotation_filter_matches_keyword_with_nonword_edge(make_seq):
    # Regression: a keyword ending in a non-word char (e.g. "MAG:") must
    # still match the real NCBI title format "MAG: Genus species ...".
    # The old \b...\b wrapping made the trailing \b after ':' unsatisfiable.
    bad = make_seq("p1", "ACGT", header="CP000001.1 MAG: Escherichia coli genome")
    ok = make_seq("p2", "ACGT", header="P2 RNA polymerase [foo]")
    report = QCReport()
    cfg = {"enabled": True, "keywords": ["MAG:"]}
    kept = annotation_filter([bad, ok], cfg, report)
    assert [s.id for s in kept] == ["p2"]
    assert report.removed_annotation == 1
    assert bad.qc_fail_reason == "annotation_keyword:MAG:"


def test_annotation_filter_word_keyword_still_bounded(make_seq):
    # A plain-word keyword must still require word boundaries — "partial"
    # should not match inside "impartiality".
    inside = make_seq("p1", "ACGT", header="P1 impartiality study")
    real = make_seq("p2", "ACGT", header="P2 partial cds")
    report = QCReport()
    cfg = {"enabled": True, "keywords": ["partial"]}
    kept = annotation_filter([inside, real], cfg, report)
    assert [s.id for s in kept] == ["p1"]


def test_annotation_filter_targets_description_not_structured_header(make_seq):
    # A keyword that appears only in a structured header field (e.g. a
    # UniProt OS= organism name) must NOT trigger removal — matching
    # targets the parsed description.
    seq = make_seq("p1", "ACGT", header="sp|P1|X pol OS=Partial-named organism sp.")
    seq.description = "RNA polymerase"          # clean description
    report = QCReport()
    cfg = {"enabled": True, "keywords": ["partial"]}
    kept = annotation_filter([seq], cfg, report)
    assert [s.id for s in kept] == ["p1"]
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
