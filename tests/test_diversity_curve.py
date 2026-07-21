"""Tests for repseq.clustering.compute_diversity_curve and the
clustering.diversity_curve_cutoffs config validation."""
from __future__ import annotations

from unittest.mock import patch

from repseq.clustering import compute_diversity_curve
from repseq.config import load_config, validate_config
from repseq.models import Cluster


def _fake_run_clustering(sequences, threshold, cfg, tmp_dir=None):
    """Stand-in for run_clustering: returns a deterministic, threshold-
    dependent number of dummy clusters so the test can verify which
    cutoffs were actually invoked and which were skipped. Uses a simple
    lookup so the test can assert exact counts."""
    counts = {0.99: 9, 0.95: 5, 0.9: 4, 0.8: 3, 0.7: 2}
    n = counts.get(round(threshold, 4), 1)
    return [
        Cluster(cluster_id=f"c{i}", representative=sequences[0])
        for i in range(n)
    ]


def test_compute_diversity_curve_invokes_each_configured_cutoff(make_seq):
    seqs = [make_seq(f"s{i}", "ACGT" * 10) for i in range(5)]
    cfg = {
        "clustering": {
            "backend": "mmseqs2",  # floor = 0.0 so all cutoffs run
            "diversity_curve_cutoffs": [0.99, 0.95, 0.9, 0.8],
        }
    }
    with patch("repseq.clustering.run_clustering", side_effect=_fake_run_clustering):
        out = compute_diversity_curve(seqs, cfg)
    # _fake_run_clustering: 0.99→9, 0.95→5, 0.9→4, 0.8→3.
    assert out == {0.99: 9, 0.95: 5, 0.9: 4, 0.8: 3}


def test_compute_diversity_curve_returns_none_below_backend_floor(make_seq):
    """Cutoffs below the active backend's identity floor (cd-hit-est: 0.80)
    must be reported as None without invoking the backend."""
    seqs = [make_seq("s1", "ACGT" * 10)]
    cfg = {
        "clustering": {
            "backend": "cdhit",
            "alphabet_for_clustering": "nucleotide",  # cd-hit-est, floor 0.80
            "diversity_curve_cutoffs": [0.99, 0.85, 0.70, 0.50],
        }
    }
    called_with: list[float] = []
    def _spy(sequences, threshold, cfg_, tmp_dir=None):
        called_with.append(threshold)
        return [Cluster(cluster_id="c", representative=sequences[0])]
    with patch("repseq.clustering.run_clustering", side_effect=_spy):
        out = compute_diversity_curve(seqs, cfg)
    # 0.70 and 0.50 are below cd-hit-est's 0.80 floor → None, no backend call.
    assert out[0.99] == 1
    assert out[0.85] == 1
    assert out[0.70] is None
    assert out[0.50] is None
    assert sorted(called_with) == [0.85, 0.99]


def test_compute_diversity_curve_disabled_when_no_cutoffs(make_seq):
    """Empty / missing cutoffs list returns None (feature off, not
    'feature on with empty result' — TSV writer uses this to skip
    the curve columns entirely)."""
    seqs = [make_seq("s1", "ACGT" * 10)]
    assert compute_diversity_curve(seqs, {"clustering": {"diversity_curve_cutoffs": []}}) is None
    assert compute_diversity_curve(seqs, {"clustering": {}}) is None
    assert compute_diversity_curve(seqs, {}) is None


def test_compute_diversity_curve_disabled_when_no_sequences():
    """No sequences → None. Avoids a meaningless 0-cluster cell."""
    cfg = {"clustering": {"diversity_curve_cutoffs": [0.95]}}
    assert compute_diversity_curve([], cfg) is None


# ---------------------------------------------------------------------------
# Console progress reporting
#
# This step runs one full clustering pass per cutoff AFTER selection has
# already settled, so on a large group it can take minutes. Without the
# per-cutoff echo it reads as a hang right after the binary search finishes.
# ---------------------------------------------------------------------------


def test_compute_diversity_curve_echoes_progress_per_cutoff(make_seq, capsys):
    """Each cutoff gets a numbered progress line naming the threshold, the
    resulting cluster count, and the elapsed time."""
    seqs = [make_seq(f"s{i}", "ACGT" * 10) for i in range(5)]
    cfg = {
        "clustering": {
            "backend": "mmseqs2",  # floor = 0.0 so all cutoffs run
            "diversity_curve_cutoffs": [0.99, 0.95, 0.9],
        }
    }
    with patch("repseq.clustering.run_clustering", side_effect=_fake_run_clustering):
        compute_diversity_curve(seqs, cfg)
    out = capsys.readouterr().out
    # Header names the step and the amount of work queued up.
    assert "diversity curve (report only): 5 sequence(s) at 3 standard cutoff(s)" in out
    # One line per cutoff, numbered i/N, with count and timing. Thresholds
    # print as `:g` — the YAML form and the TSV column form.
    assert "cutoff 1/3: threshold=0.99 → 9 cluster(s) [" in out
    assert "cutoff 2/3: threshold=0.95 → 5 cluster(s) [" in out
    assert "cutoff 3/3: threshold=0.9 → 4 cluster(s) [" in out


def test_compute_diversity_curve_echo_carries_group_label(make_seq, capsys):
    """When the caller is working through groups, the label is tagged onto
    every line — matching the binary search's `[label] ` convention — so
    interleaved group output stays attributable."""
    seqs = [make_seq("s1", "ACGT" * 10)]
    cfg = {"clustering": {"backend": "mmseqs2", "diversity_curve_cutoffs": [0.95]}}
    with patch("repseq.clustering.run_clustering", side_effect=_fake_run_clustering):
        compute_diversity_curve(seqs, cfg, label="Dengue virus")
    out = capsys.readouterr().out
    assert "[Dengue virus] diversity curve (report only)" in out
    assert "[Dengue virus] cutoff 1/1: threshold=0.95" in out


def test_compute_diversity_curve_echo_untagged_without_label(make_seq, capsys):
    """The ungrouped path (global mode) must emit NO `[...] ` tag.

    Guards the negative case: the label-carrying test above only proves a
    tag appears when asked for, and the per-cutoff assertions elsewhere
    match with or without a prefix. Without this, a leaked/stale label on
    the ungrouped path would keep the whole suite green.
    """
    seqs = [make_seq("s1", "ACGT" * 10)]
    cfg = {"clustering": {"backend": "mmseqs2", "diversity_curve_cutoffs": [0.95]}}
    with patch("repseq.clustering.run_clustering", side_effect=_fake_run_clustering):
        compute_diversity_curve(seqs, cfg)
    out = capsys.readouterr().out
    # `[` does legitimately appear in the trailing `[0.0s]` timing field, so
    # assert on the tag's own shape rather than on the bare character.
    assert "] diversity curve" not in out
    assert "] cutoff" not in out


def test_compute_diversity_curve_echo_reports_skipped_cutoffs(make_seq, capsys):
    """A cutoff below the backend floor is silent work that never happens —
    say so explicitly rather than leaving a gap in the numbering."""
    seqs = [make_seq("s1", "ACGT" * 10)]
    cfg = {
        "clustering": {
            "backend": "cdhit",
            "alphabet_for_clustering": "nucleotide",  # cd-hit-est, floor 0.80
            "diversity_curve_cutoffs": [0.99, 0.50],
        }
    }
    with patch("repseq.clustering.run_clustering", side_effect=_fake_run_clustering):
        compute_diversity_curve(seqs, cfg)
    out = capsys.readouterr().out
    assert "cutoff 2/2: threshold=0.5 → skipped (below the backend's 0.8 floor)" in out


def test_compute_diversity_curve_dedups_repeated_cutoffs(make_seq, capsys):
    """A cutoff repeated in the config must cost ONE clustering pass and be
    announced once.

    `out` is keyed by float, so a duplicate always collapsed to a single
    TSV column — but the backend was still invoked twice and the console
    would announce two cutoffs' worth of work. validate_config accepts
    duplicates, so this is reachable from a hand-edited config.
    """
    seqs = [make_seq("s1", "ACGT" * 10)]
    cfg = {
        "clustering": {
            "backend": "mmseqs2",
            "diversity_curve_cutoffs": [0.95, 0.95, 0.99],
        }
    }
    called: list[float] = []

    def _spy(sequences, threshold, cfg_, tmp_dir=None):
        called.append(threshold)
        return _fake_run_clustering(sequences, threshold, cfg_, tmp_dir)

    with patch("repseq.clustering.run_clustering", side_effect=_spy):
        out_dict = compute_diversity_curve(seqs, cfg)
    assert sorted(called) == [0.95, 0.99]        # 0.95 clustered once, not twice
    assert out_dict == {0.99: 9, 0.95: 5}
    # Announced count matches the columns actually delivered.
    assert "at 2 standard cutoff(s)" in capsys.readouterr().out


def test_compute_diversity_curve_echo_sorts_cutoffs_descending(make_seq, capsys):
    """Console order must match the TSV column order.

    output/report.py builds curve columns with `sorted(..., reverse=True)`,
    so config order must not leak into the console — otherwise a user
    cross-referencing the two maps counts to the wrong columns.
    """
    seqs = [make_seq("s1", "ACGT" * 10)]
    cfg = {
        "clustering": {
            "backend": "mmseqs2",
            "diversity_curve_cutoffs": [0.7, 0.99, 0.9],  # deliberately unsorted
        }
    }
    with patch("repseq.clustering.run_clustering", side_effect=_fake_run_clustering):
        compute_diversity_curve(seqs, cfg)
    out = capsys.readouterr().out
    assert "cutoff 1/3: threshold=0.99" in out
    assert "cutoff 2/3: threshold=0.9 " in out
    assert "cutoff 3/3: threshold=0.7 " in out


def test_compute_diversity_curve_echo_uses_thousands_separators(make_seq, capsys):
    """Large counts get `{:,}`, matching models.py / output/summary.py.

    Real pools reach six figures (a 124,548-isolate influenza run), where
    an unseparated count is the readability problem the v0.61.0 QC-summary
    rewrite existed to fix.
    """
    seqs = [make_seq(f"s{i}", "ACGT" * 10) for i in range(1234)]
    cfg = {"clustering": {"backend": "mmseqs2", "diversity_curve_cutoffs": [0.99]}}

    def _many(sequences, threshold, cfg_, tmp_dir=None):
        return [
            Cluster(cluster_id=f"c{i}", representative=sequences[0])
            for i in range(5678)
        ]

    with patch("repseq.clustering.run_clustering", side_effect=_many):
        compute_diversity_curve(seqs, cfg)
    out = capsys.readouterr().out
    assert "1,234 sequence(s)" in out
    assert "5,678 cluster(s)" in out


# ---------------------------------------------------------------------------
# Config validation
# ---------------------------------------------------------------------------

def test_validate_config_default_diversity_curve_cutoffs():
    cfg = load_config(None)
    assert cfg["clustering"]["diversity_curve_cutoffs"] == [0.99, 0.95, 0.9, 0.8, 0.7]
    assert validate_config(cfg) == []


def test_validate_config_accepts_empty_diversity_curve_cutoffs():
    cfg = load_config(None)
    cfg["clustering"]["diversity_curve_cutoffs"] = []
    assert validate_config(cfg) == []


def test_validate_config_rejects_out_of_range_diversity_cutoff():
    cfg = load_config(None)
    cfg["clustering"]["diversity_curve_cutoffs"] = [0.9, 1.5]  # >1
    errors = validate_config(cfg)
    assert any("diversity_curve_cutoffs" in e for e in errors)


def test_validate_config_rejects_non_numeric_diversity_cutoff():
    cfg = load_config(None)
    cfg["clustering"]["diversity_curve_cutoffs"] = [0.9, "0.7"]
    errors = validate_config(cfg)
    assert any("diversity_curve_cutoffs" in e for e in errors)


def test_validate_config_rejects_zero_diversity_cutoff():
    cfg = load_config(None)
    cfg["clustering"]["diversity_curve_cutoffs"] = [0.0]
    errors = validate_config(cfg)
    assert any("diversity_curve_cutoffs" in e for e in errors)
