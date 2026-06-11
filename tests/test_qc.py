"""QC pipeline: duplicates, length, ambiguous chars, annotation keywords."""
from __future__ import annotations

from repseq.models import QCReport, SequenceType
from repseq.qc.pipeline import (
    ambiguous_filter,
    annotation_filter,
    genome_length_filter,
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

def test_genome_length_filter_drops_short(make_seq):
    seqs = [
        make_seq("ok1", "A" * 100),
        make_seq("ok2", "A" * 100),
        make_seq("short", "A" * 10),
    ]
    report = QCReport()
    kept = genome_length_filter(seqs, {"min": 50}, report)
    assert {s.id for s in kept} == {"ok1", "ok2"}
    assert report.removed_length == 1
    assert "length_too_short:10<50" in {s.qc_fail_reason for s in seqs if not s.qc_passed}


def test_genome_length_filter_min_and_max(make_seq):
    seqs = [
        make_seq("tiny", "A" * 5),
        make_seq("ok", "A" * 50),
        make_seq("huge", "A" * 1000),
    ]
    report = QCReport()
    kept = genome_length_filter(seqs, {"min": 10, "max": 500}, report)
    assert [s.id for s in kept] == ["ok"]
    assert report.removed_length == 2


def test_genome_length_filter_empty_input_no_crash():
    report = QCReport()
    assert genome_length_filter([], {"min": 50}, report) == []


def test_genome_length_filter_no_bounds_is_noop(make_seq):
    # The function-level no-op (config validation separately forbids
    # enabling with no bounds): with neither min nor max, nothing drops.
    seqs = [make_seq("s1", "A" * 100), make_seq("tiny", "A" * 2)]
    report = QCReport()
    kept = genome_length_filter(seqs, {}, report)
    assert {s.id for s in kept} == {"s1", "tiny"}
    assert report.removed_length == 0


def test_run_qc_skips_genome_length_filter_in_segmented_mode(make_seq):
    # Regression: in segmented mode the input is a mix of long and short
    # segments, so the whole-genome length filter must never run there
    # (per-segment segment_lengths bounds do the job downstream). Even with
    # the filter "enabled" the segmented branch wins and skips it.
    seqs = [
        make_seq("L1", "A" * 6000),
        make_seq("M1", "A" * 4000),
        make_seq("S1", "A" * 900),
    ]
    cfg = {"qc": {"genome_length_filter": {"enabled": True, "min": 5000}},
           "segmented": {"enabled": True}}
    kept, report = run_qc(seqs, cfg)
    assert {s.id for s in kept} == {"L1", "M1", "S1"}
    assert report.removed_length == 0
    assert report.length_filter_skipped is True
    # Reason is recorded; the segmented summary footer notes the
    # whole-genome length filter was not applied (per-segment bounds run
    # downstream instead).
    assert report.length_filter_skip_reason == "segmented"
    out = report.summary()
    assert "not applied in segmented mode" in out
    assert "length filters" in out


def test_run_qc_skips_genome_length_filter_when_disabled(make_seq):
    # Non-segmented but the filter is disabled (the default): no drops,
    # length_filter_skipped recorded.
    seqs = [make_seq("L1", "A" * 6000), make_seq("S1", "A" * 50)]
    cfg = {"qc": {"genome_length_filter": {"enabled": False, "min": 5000}},
           "segmented": {"enabled": False}}
    kept, report = run_qc(seqs, cfg)
    assert {s.id for s in kept} == {"L1", "S1"}
    assert report.removed_length == 0
    assert report.length_filter_skipped is True
    # Non-segmented + disabled must NOT claim "segmented mode": the
    # non-segmented summary renders a single per-sequence tally and simply
    # omits the (un-run) whole-genome length-filter line.
    assert report.length_filter_skip_reason == "disabled"
    summary = report.summary()
    assert "segmented mode" not in summary
    assert "Per-record screening  (unit: sequences)" in summary


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


def test_run_qc_applies_genome_length_filter_when_enabled(make_seq):
    seqs = [
        make_seq("L1", "A" * 6000),
        make_seq("M1", "A" * 4000),
        make_seq("S1", "A" * 900),
    ]
    cfg = {"qc": {"genome_length_filter": {"enabled": True, "min": 2000}},
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
            "genome_length_filter": {"enabled": True, "min": 50},
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


def test_qc_summary_segmented_running_tally_and_unit_change():
    """The segmented summary makes the segment->isolate unit change
    explicit: 'input records' (segments) at the top, 'complete isolates
    assembled' as the Phase 2 endpoint, and a reconciliation footer tying
    the two units together — this is the line that dissolves the
    '1.6M segments -> 12K isolates' surprise."""
    from repseq.models import QCReport
    report = QCReport(total_input=100, passed=80)
    report.final_survivors = 12
    report.final_survivors_unit = "isolates"
    report.display["segmented"] = True
    out = report.summary()
    assert "input records" in out
    assert "complete isolates assembled" in out
    assert "100 input segments → 12 complete isolates" in out


def test_qc_summary_omits_final_survivors_when_unset():
    """A QCReport built without setting final_survivors (e.g. legacy
    code paths, unit tests) should not show the new line at all."""
    from repseq.models import QCReport
    report = QCReport(total_input=100, passed=80)
    out = report.summary()
    assert "Final survivors" not in out


def _segmented_report():
    """An Influenza-like segmented QCReport with a fully populated display
    context — the shape the CLI driver builds for the console summary."""
    from repseq.models import QCReport
    r = QCReport(total_input=1_598_432)
    r.removed_annotation = 12_300
    r.removed_ambiguous = 4_210
    r.removed_proteins = 3_120
    r.removed_taxonomy_mismatch = 940
    r.removed_protein_quality = 3
    r.removed_hmm_failed = 16          # 12 pre-assembly + 4 at concat
    r.removed_incomplete_isolates = 1_478_902 + 1
    r.final_survivors = 12_031
    r.final_survivors_unit = "isolates"
    r.display.update({
        "segmented": True,
        "ambiguous_pct": 0.02,
        "protein_quality_pct": 0.05,
        "annotation_keywords": ["MAG:", "synthetic", "fragment", "partial"],
        "segments_entering_assembly": 1_574_950,
        "n_segments": 8,
        "hmm_isolates_pre_assembly": 12,
        "incomplete_segments": 1_478_902,
        "concat_marker_isolates": 5,
        "concat_hmm_isolates": 4,
        "duplicate_isolates": 210,
    })
    return r


def test_qc_summary_segmented_phases_thresholds_and_units():
    out = _segmented_report().summary()
    # Two explicit phases with unit headers.
    assert "Phase 1 · per-record screening  (unit: segments)" in out
    assert "Phase 2 · isolate assembly  (8 segments → 1 isolate; unit: isolates)" in out
    # Live thresholds rendered into the labels.
    assert "ambiguous nucleotides (>2% N)" in out
    assert "ambiguous protein residues (>5% X)" in out
    assert "description keywords (MAG:, synthetic, fragment, …)" in out
    # Isolate-level gate sub-tally (counted in isolates, not segments).
    assert "isolate-level identity gates" in out
    assert "marker fails HMM identity" in out
    # The unit-change reconciliation line + footer.
    assert "segments entering assembly" in out
    assert "complete isolates assembled" in out
    assert "1,598,432 input segments → 12,031 complete isolates" in out


def test_qc_summary_segmented_phase2_separates_missing_from_concat_hmm():
    """Phase 2 attributes the big drop (missing segments) separately from
    the small concatenation/HMM marker loss, with the HMM share noted."""
    out = _segmented_report().summary()
    assert "missing ≥1 of 8 segments / unidentifiable" in out
    assert "≈1,478,902 segments" in out          # the cliff, in segments
    assert "marker lost at concatenation" in out
    assert "(4 via HMM)" in out                   # HMM share of the concat drop
    assert "identical duplicate genomes" in out


def test_qc_summary_running_tally_is_arithmetic():
    """Each Phase 1 survivor number equals input minus the cumulative
    per-segment removals shown above it."""
    out = _segmented_report().summary()
    # input 1,598,432 − 12,300 = 1,586,132 ; − 4,210 = 1,581,922 ; etc.
    assert "1,586,132" in out
    assert "1,581,922" in out
    assert "1,578,802" in out
    assert "1,577,862" in out          # after the cross-segment mismatch drop
    # ground-truth pool entering assembly (after the isolate gates)
    assert "1,574,950" in out


def test_qc_summary_nonsegmented_single_tally_with_bounds():
    from repseq.models import QCReport
    r = QCReport(total_input=84_201)
    r.removed_duplicates = 9_140
    r.removed_length = 1_002
    r.removed_ambiguous = 410
    r.removed_hmm_failed = 330
    r.final_survivors = 73_319
    r.final_survivors_unit = "sequences"
    r.display.update({
        "segmented": False,
        "ambiguous_pct": 0.05,
        "length_bounds": {"min": 1_000, "max": 15_000},
        "annotation_keywords": [],
    })
    out = r.summary()
    assert "Per-record screening  (unit: sequences)" in out
    assert "Phase 1" not in out                    # single phase, no segments
    assert "genome length out of bounds (<1,000 or >15,000 nt)" in out
    assert "exact duplicates" in out
    assert "84,201 input sequences → 73,319 passing QC" in out
