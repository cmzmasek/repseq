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
