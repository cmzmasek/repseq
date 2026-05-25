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
        "qc": {"genome_length_filter": {"enabled": True, "min": 9000, "max": 13000},
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


def test_render_summary_non_segmented_describes_absolute_genome_bounds(make_seq, tmp_path):
    """A non-segmented run's length-filter prose must describe the absolute
    nucleotide bounds of qc.genome_length_filter — no median/relative
    wording (that filter was removed)."""
    cfg = _base_cfg(tmp_path)  # genome_length_filter min=9000, max=13000
    qc = _qc(total_input=1000, passed=970, dedup=0, length=30, ambig=0)
    # Non-segmented: length_filter_skipped stays False (default).
    md = render_summary(cfg, qc, _result(make_seq), ["a.fasta"])
    assert "**30** outside the configured whole-genome length bounds" in md
    assert "shorter than 9,000 nt" in md
    assert "longer than 13,000 nt" in md
    # Regression guards: the removed median wording must not reappear.
    assert "median" not in md
    # "per-rank" used to be a regression guard for an older relative-bounds
    # phrasing in this same filter section; v0.23 introduced legitimate
    # "per-rank" prose for the taxonomic / nucleotide-length reports, so
    # we drop the substring guard and rely on the "median" check instead.
    # And the per-segment sentence must NOT appear when removed_length_by_segment is empty.
    assert "per-segment length filter" not in md


def test_render_summary_length_filter_lower_bound_only(make_seq, tmp_path):
    """With only a min bound set, the prose names just the lower bound."""
    cfg = _base_cfg(tmp_path)
    cfg["qc"]["genome_length_filter"] = {"enabled": True, "min": 9000, "max": None}
    qc = _qc(total_input=1000, passed=970, dedup=0, length=30, ambig=0)
    md = render_summary(cfg, qc, _result(make_seq), ["a.fasta"])
    assert "shorter than 9,000 nt" in md
    assert "longer than" not in md


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
    # Colouring is on by default (genus mode).
    assert "coloured by **genus**" in md


def test_render_summary_phylo_describes_partitioned_supermatrix(make_seq, tmp_path):
    """protein + IQ-TREE + partition enabled (the default) → the phylo
    section describes the partitioned-supermatrix analysis, not concat."""
    cfg = _base_cfg(tmp_path)  # alphabet_for_clustering=protein
    cfg["phylo"] = {
        "tool": "iqtree",
        "partition": {"enabled": True, "linkage": "unlinked"},
    }
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "partitioned supermatrix" in md
    assert "per partition" in md
    assert "aligned **separately**" in md
    # Linkage surfaced as words + IQ-TREE flag.
    assert "unlinked" in md and "`-Q`" in md
    assert "{prefix}_partition.nex" in md


def test_render_summary_phylo_describes_trimal_when_enabled(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    # FastTree path → concat branch; trimming enabled on the genome tree.
    cfg["phylo"] = {"tool": "fasttree", "trimal": {"enabled": True, "mode": "automated1"}}
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "trimAl" in md
    assert "`-automated1`" in md
    assert "{prefix}_msa_untrimmed.fasta" in md
    # Software table gains a trimAl row.
    assert "Alignment trimming" in md


def test_render_summary_phylo_no_trimal_prose_when_off(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "fasttree"}  # trimal defaults off
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "trimAl" not in md


def test_render_summary_per_protein_describes_trimal_when_enabled(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {
        "tool": "fasttree",
        "per_protein": {"min_taxa": 3, "trimal": {"enabled": True, "mode": "gappyout"}},
    }
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], per_protein_ran=True)
    assert "trimAl" in md and "`-gappyout`" in md


def test_render_summary_phylo_concat_when_partition_disabled(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "iqtree", "partition": {"enabled": False}}
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "partitioned supermatrix" not in md
    assert "multiple sequence alignment of the" in md


def test_render_summary_phylo_describes_fast_mode_mafft(make_seq, tmp_path):
    """--fast forces MAFFT '--retree 1' with use_auto=False — the summary
    must say so (and flag it as a preliminary-run setting, not pretend
    --auto ran)."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {
        "tool": "fasttree",
        "partition": {"enabled": False},
        "trimal": {"enabled": False},
        "mafft": {"extra_args": ["--retree", "1"], "use_auto": False},
    }
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "`--retree 1`" in md
    assert "single-pass FFT-NS-1" in md
    assert "`--fast` preliminary-run setting" in md
    # And must NOT lie about --auto.
    assert "`--auto`" not in md


def test_render_summary_per_protein_fast_mode_does_not_claim_l_ins_i(make_seq, tmp_path):
    """Per-protein branch must not call --retree 1 'high-accuracy L-INS-i'."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {
        "tool": "fasttree",
        "per_protein": {
            "min_taxa": 3,
            "mafft": {"extra_args": ["--retree", "1"]},
        },
    }
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], per_protein_ran=True,
    )
    assert "high-accuracy L-INS-i" not in md
    assert "`--retree 1`" in md
    assert "`--fast` preliminary-run setting" in md


def test_render_summary_coloring_two_rank_and_disabled(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {
        "tool": "fasttree",
        "coloring": {"enabled": True, "ranks": ["genus", "subgenus"]},
    }
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "coloured by **genus**" in md
    assert "**subgenus** shaded within its parent genus" in md

    cfg["phylo"]["coloring"] = {"enabled": False}
    md_off = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "Tree leaves were coloured" not in md_off


def test_render_summary_per_protein_section_appears(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "fasttree", "per_protein": {"min_taxa": 3}}
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], per_protein_ran=True,
    )
    assert "## Phylogenetic inference" in md
    assert "separate tree was built for each HMM" in md
    # OR-of-architectures prose (one tree per spec, any token satisfies).
    assert "alternative domain architectures" in md
    assert "one tree per `hmms:` spec" in md
    assert "{prefix}_per_protein/" in md
    assert "reassortment" in md
    # MAFFT/tree-builder software rows show even though only 2F ran.
    assert "MAFFT" in md
    # Incongruence table described by default.
    assert "Robinson-Foulds" in md
    assert "{prefix}_per_protein/{prefix}_incongruence.tsv" in md
    # Domain-architecture annotation described by default.
    assert "domain architecture" in md
    assert "<domain_architecture>" in md


def test_render_summary_per_protein_describes_linsi_mafft(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {
        "tool": "fasttree",
        "per_protein": {
            "min_taxa": 3,
            "mafft": {"extra_args": ["--maxiterate", "1000", "--localpair"]},
        },
    }
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], per_protein_ran=True,
    )
    assert "--maxiterate 1000 --localpair" in md
    assert "L-INS-i" in md


def test_render_summary_per_protein_non_linsi_args_not_called_linsi(
    make_seq, tmp_path,
):
    """User passing arbitrary MAFFT args (not the L-INS-i combo) must not
    have the prose claim they ran L-INS-i — the previous wording lied about
    any non-default flag set."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {
        "tool": "fasttree",
        "per_protein": {
            "min_taxa": 3,
            # --maxiterate alone is iterative refinement on FFT-NS, NOT
            # L-INS-i (which strictly requires --localpair).
            "mafft": {"extra_args": ["--maxiterate", "2"]},
        },
    }
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], per_protein_ran=True,
    )
    assert "--maxiterate 2" in md
    assert "L-INS-i" not in md
    assert "user-supplied" in md


def test_render_summary_incongruence_omitted_when_disabled(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {
        "tool": "fasttree",
        "per_protein": {"min_taxa": 3, "incongruence": False},
    }
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], per_protein_ran=True,
    )
    assert "separate tree was built for each HMM" in md
    assert "Robinson-Foulds" not in md


def test_render_summary_no_phylo_section_when_neither_ran(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "## Phylogenetic inference" not in md


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
    """When _hmm_runtime.active is True, the selection section references
    the upstream HMM-QC step and names the database / citation. When
    HMM drops happened, the QC section has a dedicated bullet listing
    the cutoffs in effect. Per the summary-renderer-drift memory."""
    cfg = _base_cfg(tmp_path)
    cfg["_hmm_runtime"] = {"active": True, "ga_cutoffs": {}, "hmm_cfg": cfg.get("hmm", {})}
    qc = _qc()
    qc.removed_hmm_failed = 5
    qc.removed_hmm_by_marker = {"L:RdRP_4": 3, "M:Bunya_G2--Bunya_G1": 2}
    md = render_summary(cfg, qc, _result(make_seq), ["a.fasta"])
    # Selection section references the upstream QC step.
    assert "HMMER hmmscan" in md
    assert "bundled viral-core profile set" in md
    assert "pre-filtered by an HMM-based identity QC step" in md
    # QC section's HMM bullet lists drop count + coverage cutoff.
    assert "**5** isolates were dropped by the HMM-based identity QC" in md or \
           "**5** sequences were dropped by the HMM-based identity QC" in md
    # Coverage cutoff rendered as a percent (default 0.5 → 50%).
    assert "≥ 50%" in md or "50 %" in md


def test_render_summary_omits_hmm_sentence_when_not_active(make_seq, tmp_path):
    """No _hmm_runtime → no HMM sentence (don't add false claims to the
    Methods section if the HMM tier never ran)."""
    cfg = _base_cfg(tmp_path)
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "HMMER hmmscan" not in md
    assert "bundled viral-core" not in md


def test_render_summary_mentions_extra_segments_when_dropped(make_seq, tmp_path):
    """When the extra-segments check drops any isolates, the QC section
    must surface the count (in isolates) and name the config knob the
    user can flip back to 'warn'."""
    cfg = _base_cfg(tmp_path)
    qc = _qc()
    qc.removed_extra_segments = 7
    md = render_summary(cfg, qc, _result(make_seq), ["a.fasta"])
    assert "**7** isolate(s) were dropped by the extra-segments check" in md
    assert "segmented.extra_segments_action: drop" in md
    assert "extra_segments:<extras>" in md


def test_render_summary_omits_extra_segments_when_counter_zero(make_seq, tmp_path):
    """No drops → no sentence. Otherwise every non-segmented run's QC
    section gains an irrelevant bullet."""
    cfg = _base_cfg(tmp_path)
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "extra-segments check" not in md


def test_render_summary_mentions_protein_quality_when_dropped(make_seq, tmp_path):
    """When the protein-quality check drops records, the QC section names
    the count, the threshold, the X/B/Z/J residue set, and the reason
    prefix."""
    cfg = _base_cfg(tmp_path)
    cfg["qc"]["protein_quality"] = {"enabled": True, "max_bad_fraction": 0.05}
    qc = _qc()
    qc.removed_protein_quality = 11
    md = render_summary(cfg, qc, _result(make_seq), ["a.fasta"])
    assert "**11** sequences were dropped by the protein-quality check" in md
    assert "5%" in md
    assert "X/B/Z/J" in md
    assert "protein_quality:" in md


def test_render_summary_omits_protein_quality_when_zero(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "protein-quality check" not in md


def test_render_summary_points_to_taxonomic_report(make_seq, tmp_path):
    """The selection section must point the reader at the new
    taxonomic-diversity report file."""
    cfg = _base_cfg(tmp_path)
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "{prefix}_taxonomic_report.txt" in md
    assert "Taxonomic diversity at each rank before and after clustering" in md


def test_render_summary_points_to_nucleotide_taxonomic_report(make_seq, tmp_path):
    """Summary must also point the reader at the NT length-statistics report."""
    cfg = _base_cfg(tmp_path)
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "{prefix}_nucleotide_taxonomic_report.txt" in md


def test_render_summary_protein_annotation_segmented_uses_segment_wording(make_seq, tmp_path):
    """In segmented mode the protein-annotation drop is per NCBI segment
    record, compared against expected_proteins_per_segment. It must appear
    as its own sentence (not lumped into the 'Initial QC' list) and name
    the per-segment config key."""
    cfg = _base_cfg(tmp_path)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "hantaviridae",
        "viruses": {
            "hantaviridae": {
                "expected_segments": 3,
                "segments": ["S", "M", "L"],
                "expected_proteins_per_segment": {"S": [1, 2], "M": 1, "L": 1},
            },
        },
    }
    qc = _qc()
    qc.removed_proteins = 174
    md = render_summary(cfg, qc, _result(make_seq), ["a.fasta"])
    assert "**174** segment record(s) were dropped by the protein-annotation check" in md
    assert "segmented.viruses.hantaviridae.expected_proteins_per_segment" in md
    assert "protein_count_mismatch:segment=<seg>" in md
    # Must NOT be in the old "Initial QC" list wording.
    assert "failing the protein-annotation check" not in md


def test_render_summary_protein_annotation_non_segmented_uses_min_proteins(make_seq, tmp_path):
    """Non-segmented run with a global min_proteins floor names the floor
    and the matching _qc_removed.tsv reason."""
    cfg = _base_cfg(tmp_path)
    cfg["qc"]["protein_annotation"] = {"enabled": True, "min_proteins": 2}
    qc = _qc()
    qc.removed_proteins = 9
    md = render_summary(cfg, qc, _result(make_seq), ["a.fasta"])
    assert "**9** sequence(s) were dropped by the protein-annotation check" in md
    assert "qc.protein_annotation.min_proteins = 2" in md
    assert "protein_count_below_min:<n><2" in md


def test_render_summary_omits_protein_annotation_when_counter_zero(make_seq, tmp_path):
    """No protein drops → no sentence anywhere in the QC section."""
    cfg = _base_cfg(tmp_path)
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "protein-annotation check" not in md


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


def test_render_summary_phylo_describes_outgroup_rooting(make_seq, tmp_path):
    """method=outgroup with an accession spec → prose names the
    accession and the MRCA fallback, NOT 'MAD/midpoint fallbacks'."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {
        "tool": "fasttree",
        "rooting": {"method": "outgroup", "outgroup": "AB123456"},
    }
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "AB123456" in md
    assert "MRCA" in md
    # The generic MAD/midpoint fallback phrasing must not appear for outgroup.
    assert "minimum ancestor deviation" not in md


def test_render_summary_phylo_describes_outgroup_rank(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {
        "tool": "fasttree",
        "rooting": {
            "method": "outgroup",
            "outgroup_rank": {"family": "Hantaviridae"},
        },
    }
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "Hantaviridae" in md
    assert "family" in md


def test_render_summary_phylo_describes_midpoint_rooting_honestly(
    make_seq, tmp_path,
):
    """method=midpoint has no fallback chain; prose must not claim one."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "fasttree", "rooting": {"method": "midpoint"}}
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "midpoint rooting" in md
    assert "fallback" not in md.lower()


def test_render_summary_describes_per_segment_section_when_2h_ran(
    make_seq, tmp_path,
):
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "fasttree"}
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], per_segment_ran=True,
    )
    assert "nucleotide tree was also built per segment" in md
    assert "{prefix}_per_segment/" in md
    # Reassortment is the headline motivation — verify it's mentioned.
    assert "reassortment" in md
