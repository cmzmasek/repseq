"""Tests for the Methods-section summary writer (repseq.output.summary)."""
from __future__ import annotations

from unittest.mock import patch

import pytest

from repseq.models import Cluster, QCReport, RunResult
from repseq.output.summary import (
    detect_tool_versions,
    render_summary,
    write_summary,
)


def _base_cfg(tmp_path):
    return {
        "output": {"dir": str(tmp_path), "prefix": "test"},
        "clustering": {"backend": "mmseqs2", "alphabet_for_clustering": "protein"},
        "qc": {"length_filter": {"mode": "median_percent", "min_percent": 50},
               "ambiguous_threshold": 0.05},
        "representative": {"priority": ["refseq", "reviewed_uniprot", "longest"]},
        "segmented": {"enabled": False},
    }


def _result(make_seq, n_reps=3, n_clusters=3, mode="global"):
    reps = [make_seq(f"r{i}", "ACGT" * 10) for i in range(n_reps)]
    clusters = [
        Cluster(cluster_id=f"c{i}", representative=reps[i], members=[])
        for i in range(min(n_clusters, n_reps))
    ]
    return RunResult(mode=mode, representatives=reps, clusters=clusters)


def _qc(total_input=100, passed=92, dedup=4, length=2, ambig=2):
    r = QCReport()
    r.total_input = total_input
    r.passed = passed
    r.removed_duplicates = dedup
    r.removed_length = length
    r.removed_ambiguous = ambig
    return r


# ---------------------------------------------------------------------------
# detect_tool_versions
# ---------------------------------------------------------------------------

def test_detect_tool_versions_returns_dict_with_known_keys():
    versions = detect_tool_versions()
    # We don't assert which binaries are installed (CI may have none) —
    # only that every probed tool appears as a key.
    for binary in ("cd-hit", "cd-hit-est", "mmseqs", "mafft", "FastTree", "iqtree2"):
        assert binary in versions
        # Either a string (detected) or None (not on PATH or probe failed).
        assert versions[binary] is None or isinstance(versions[binary], str)


# ---------------------------------------------------------------------------
# render_summary — top-level shape and conditional sections
# ---------------------------------------------------------------------------

def test_render_summary_includes_core_sections(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "# Methods — repseq global selection" in md
    assert "## Input" in md
    assert "## Quality control" in md
    assert "## Representative selection" in md
    assert "## Software and references" in md
    # Phylogeny and segmented sections must NOT appear in this non-segmented,
    # phylo-off run.
    assert "## Phylogenetic inference" not in md
    assert "## Segmented-virus handling" not in md


def test_render_summary_qc_numbers_pulled_from_report(make_seq, tmp_path):
    qc = _qc(total_input=1234, passed=1200, dedup=10, length=20, ambig=4)
    md = render_summary(_base_cfg(tmp_path), qc, _result(make_seq), ["a.fasta"])
    assert "**1,234**" in md           # input count (thousands separator)
    assert "**10** exact duplicates" in md
    assert "**20** outside" in md
    assert "**4** with > 5% ambiguous" in md
    assert "*1,200* passed basic QC" in md


def test_render_summary_segmented_length_filter_uses_per_segment_wording(make_seq, tmp_path):
    """When the segmented per-segment length filter fired, the QC section
    must describe the per-segment drops in isolate units — NOT the
    global ±median-percent window (which never ran in segmented mode).
    Regression guard for the v0.13.x bug where summary.md mislabelled the
    per-segment counter as the global filter."""
    cfg = _base_cfg(tmp_path)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "peribunyaviridae",
        "viruses": {
            "peribunyaviridae": {"expected_segments": 3, "segments": ["S", "M", "L"]},
        },
    }
    qc = _qc(total_input=2000, passed=1900, dedup=0, length=0, ambig=0)
    # Mirror what run_qc + segment_length_filter would set in segmented mode:
    qc.length_filter_skipped = True
    qc.removed_length = 942  # in segments — must NOT appear as the global count
    qc.removed_length_by_segment = {
        "S": {"too_short": 200, "too_long": 0},
        "M": {"too_short": 0,   "too_long": 100},
        "L": {"too_short": 14,  "too_long": 0},
    }
    md = render_summary(cfg, qc, _result(make_seq), ["a.fasta"])
    # The misleading global-filter wording must NOT appear.
    assert "outside the configured length window" not in md
    assert "outside ±50% of the per-rank median length" not in md
    # The actionable per-segment sentence MUST appear, in isolate units
    # (sum of the per-segment counters), with the breakdown spelled out.
    assert "**314** isolate(s) were dropped by the per-segment length filter" in md
    assert "S too short: 200" in md
    assert "M too long: 100" in md
    assert "L too short: 14" in md


def test_render_summary_non_segmented_keeps_global_length_wording(make_seq, tmp_path):
    """Regression guard for the inverse: a non-segmented run with a
    global-filter drop still gets the ±median-percent wording."""
    cfg = _base_cfg(tmp_path)
    qc = _qc(total_input=1000, passed=970, dedup=0, length=30, ambig=0)
    # Non-segmented: length_filter_skipped stays False (default).
    md = render_summary(cfg, qc, _result(make_seq), ["a.fasta"])
    assert "**30** outside the configured length window" in md
    assert "outside ±50% of the per-rank median length" in md
    # And the per-segment sentence must NOT appear when removed_length_by_segment is empty.
    assert "per-segment length filter" not in md


def test_render_summary_segmented_block_appears_only_when_enabled(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "hantaviridae",
        "viruses": {
            "hantaviridae": {"expected_segments": 3, "segments": ["S", "M", "L"]},
        },
    }
    qc = _qc(total_input=500, passed=480)
    qc.final_survivors = 42
    qc.final_survivors_unit = "isolates"
    complete = {f"iso{i}": [None, None, None] for i in range(50)}
    md = render_summary(cfg, qc, _result(make_seq), ["a.fasta"],
                        complete_isolates=complete, segment_names=["S", "M", "L"])
    assert "## Segmented-virus handling" in md
    assert "**3** expected segments" in md
    assert "S, M, L" in md
    assert "**50**" in md  # complete isolates
    assert "**42**" in md  # final survivors


def test_render_summary_phylo_section_appears_when_phylo_ran(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "fasttree", "rooting": {"method": "taxonomy_guided"}}
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "## Phylogenetic inference" in md
    assert "MAFFT" in md
    assert "FastTree" in md
    # IQ-TREE row should still appear in the software table (chosen-or-not).
    assert "IQ-TREE" in md


def test_render_summary_protein_alphabet_described(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    # protein + non-segmented → "marker-protein sequences"
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "marker-protein sequences" in md


def test_render_summary_nucleotide_alphabet_described(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["clustering"]["alphabet_for_clustering"] = "nucleotide"
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "input nucleotide sequences" in md
    # No protein-only language in the selection paragraph.
    assert "marker-protein" not in md.split("## Software")[0]


def test_render_summary_cdhit_backend_named_in_selection(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["clustering"]["backend"] = "cdhit"
    cfg["clustering"]["alphabet_for_clustering"] = "nucleotide"
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    # When NT + cdhit, the dispatcher would pick cd-hit-est.
    assert "cd-hit-est" in md


def test_render_summary_mentions_diversity_curve_when_configured(make_seq, tmp_path):
    """When clustering.diversity_curve_cutoffs is non-empty, the selection
    section names the cutoffs and points the reader at group_counts.tsv."""
    cfg = _base_cfg(tmp_path)
    cfg["clustering"]["diversity_curve_cutoffs"] = [0.99, 0.95, 0.9]
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "fixed identity thresholds (0.99, 0.95, 0.9)" in md
    assert "group_counts.tsv" in md
    assert "did not influence representative selection" in md


def test_render_summary_omits_diversity_curve_sentence_when_disabled(make_seq, tmp_path):
    """Empty cutoff list → no sentence at all (don't pollute Methods with
    a feature that wasn't on)."""
    cfg = _base_cfg(tmp_path)
    cfg["clustering"]["diversity_curve_cutoffs"] = []
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "fixed identity thresholds" not in md
    assert "diagnostic of within-stratum" not in md


def test_render_summary_mentions_hmm_tier_when_active(make_seq, tmp_path):
    """When _hmm_runtime.active is True, the selection section names the
    database, cites HMMER, and lists the cutoffs in effect. Per the
    summary-renderer-drift memory."""
    cfg = _base_cfg(tmp_path)
    cfg["_hmm_runtime"] = {"active": True, "ga_cutoffs": {}, "hmm_cfg": cfg.get("hmm", {})}
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "HMMER hmmscan" in md
    assert "bundled viral-core profile set" in md
    # Coverage cutoff rendered as a percent (default 0.5 → 50%).
    assert "≥ 50%" in md or "50 %" in md


def test_render_summary_omits_hmm_sentence_when_not_active(make_seq, tmp_path):
    """No _hmm_runtime → no HMM sentence (don't add false claims to the
    Methods section if the HMM tier never ran)."""
    cfg = _base_cfg(tmp_path)
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "HMMER hmmscan" not in md
    assert "bundled viral-core" not in md


def test_render_summary_hmm_software_row_only_when_active(make_seq, tmp_path):
    """Software table gains an HMMER row only when the HMM tier fired —
    otherwise it'd clutter every run's Methods section."""
    cfg = _base_cfg(tmp_path)
    md_off = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "HMMER hmmscan" not in md_off
    cfg["_hmm_runtime"] = {"active": True}
    md_on = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "HMMER hmmscan" in md_on


def test_render_summary_software_table_marks_used_and_unused(make_seq, tmp_path):
    """When backend=mmseqs2, cd-hit row must read '(not used)' and vice versa."""
    cfg = _base_cfg(tmp_path)
    with patch("repseq.output.summary.detect_tool_versions",
               return_value={k: "9.9.9" for k in
                             ("cd-hit", "cd-hit-est", "mmseqs", "mafft", "FastTree", "iqtree2")}):
        md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "MMseqs2 | 9.9.9 | Sequence clustering (used)" in md
    assert "cd-hit | 9.9.9 | Sequence clustering (not used)" in md


# ---------------------------------------------------------------------------
# write_summary — file IO + path return
# ---------------------------------------------------------------------------

def test_write_summary_creates_file_at_prefix_summary_md(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    path = write_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert path.exists()
    assert path.name == "test_summary.md"
    content = path.read_text()
    assert content.startswith("# Methods — ")
    assert "Auto-generated" in content
