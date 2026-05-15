"""CLI helpers: the closing run summary / no-output warning."""
from __future__ import annotations

from repseq.cli import _final_summary, _handle_segmented
from repseq.models import Cluster, QCReport, RunResult


def _result(make_seq, n_reps=0, n_clusters=0):
    reps = [make_seq(f"r{i}", "ACGT" * 10) for i in range(n_reps)]
    clusters = []
    for i in range(n_clusters):
        rep = make_seq(f"c{i}", "ACGT" * 10)
        clusters.append(Cluster(cluster_id=f"c{i}", representative=rep, members=[]))
    return RunResult(mode="global", representatives=reps, clusters=clusters)


def test_final_summary_reports_success(make_seq, capsys):
    result = _result(make_seq, n_reps=12, n_clusters=4)
    report = QCReport(total_input=100, passed=40)
    _final_summary(result, report, {"segmented": {"enabled": False}})
    out = capsys.readouterr().out
    assert "selected 12 representative sequence(s) across 4 cluster(s)" in out
    assert "40 of 100 input sequences passed QC" in out


def test_final_summary_success_segmented_uses_isolate_wording(make_seq, capsys):
    result = _result(make_seq, n_reps=3)
    report = QCReport(total_input=30, passed=24)
    _final_summary(result, report, {"segmented": {"enabled": True}})
    out = capsys.readouterr().out
    assert "selected 3 isolate(s)" in out


def test_final_summary_warns_when_no_sequences_loaded(make_seq, capsys):
    result = _result(make_seq, n_reps=0)
    report = QCReport(total_input=0, passed=0)
    _final_summary(result, report, {"segmented": {"enabled": False}})
    err = capsys.readouterr().err
    assert "WARNING: no representative sequences" in err
    assert "No sequences were loaded" in err


def test_final_summary_warns_when_qc_removed_everything(make_seq, capsys):
    result = _result(make_seq, n_reps=0)
    report = QCReport(
        total_input=50, passed=0, removed_length=30, removed_annotation=20
    )
    _final_summary(result, report, {"segmented": {"enabled": False}})
    err = capsys.readouterr().err
    assert "QC removed all 50 input sequences" in err
    assert "30 on length" in err
    assert "20 on annotation keywords" in err


def test_final_summary_warns_when_segmented_filter_dropped_everything(make_seq, capsys):
    result = _result(make_seq, n_reps=0)
    report = QCReport(
        total_input=40, passed=32, removed_incomplete_isolates=32
    )
    _final_summary(result, report, {"segmented": {"enabled": True}})
    err = capsys.readouterr().err
    assert "completeness/length filter" in err
    assert "isolate_regex" in err


def test_final_summary_warns_when_selection_produced_nothing(make_seq, capsys):
    result = _result(make_seq, n_reps=0)
    report = QCReport(total_input=40, passed=40)
    _final_summary(result, report, {"segmented": {"enabled": False}})
    err = capsys.readouterr().err
    assert "selection produced" in err
    assert "MMseqs2" in err


# ---------------------------------------------------------------------------
# _handle_segmented — dedup applied to concatenated isolates, not segments
# ---------------------------------------------------------------------------

def _segmented_cfg():
    return {
        "segmented": {
            "enabled": True,
            "virus": "test",
            "viruses": {
                "test": {
                    "expected_segments": 2,
                    "segments": ["S", "L"],
                    "isolate_regex": r"(?P<isolate>iso\d+)",
                }
            },
        }
    }


def test_handle_segmented_keeps_isolates_sharing_a_conserved_segment(make_seq):
    # Regression: iso1 and iso2 have a byte-identical S segment but differ on
    # L. Neither isolate may be dropped — they are not duplicates.
    seqs = [
        make_seq("a_S", "ACGTACGTACGT", header="iso1", segment="S"),
        make_seq("a_L", "TTTTGGGGTTTT", header="iso1", segment="L"),
        make_seq("b_S", "ACGTACGTACGT", header="iso2", segment="S"),  # identical S
        make_seq("b_L", "CCCCAAAACCCC", header="iso2", segment="L"),
    ]
    report = QCReport()
    report.dedup_skipped = True  # run_qc would have set this in segmented mode
    concat, complete, segments = _handle_segmented(seqs, _segmented_cfg(), report)
    assert sorted(s.isolate_id for s in concat) == ["iso1", "iso2"]
    assert set(complete.keys()) == {"iso1", "iso2"}
    assert report.removed_duplicates == 0


def test_handle_segmented_collapses_fully_identical_isolates(make_seq):
    # iso3 is byte-identical to iso1 across *both* segments — a true duplicate
    # isolate, which must collapse to one concatenated representative.
    seqs = [
        make_seq("a_S", "ACGTACGTACGT", header="iso1", segment="S"),
        make_seq("a_L", "TTTTGGGGTTTT", header="iso1", segment="L"),
        make_seq("c_S", "ACGTACGTACGT", header="iso3", segment="S"),
        make_seq("c_L", "TTTTGGGGTTTT", header="iso3", segment="L"),
    ]
    report = QCReport()
    report.dedup_skipped = True
    concat, complete, segments = _handle_segmented(seqs, _segmented_cfg(), report)
    assert [s.isolate_id for s in concat] == ["iso1"]
    assert set(complete.keys()) == {"iso1"}
    assert report.removed_duplicates == 1
