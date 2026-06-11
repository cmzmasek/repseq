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

def test_at_a_glance_block_present_and_run_specific(make_seq, tmp_path):
    """The lead 'Analysis at a glance' table states the actual decisions
    for this run (dataset type, clustering substrate, tool, counts)."""
    cfg = _base_cfg(tmp_path)
    qc = _qc(total_input=5231)
    md = render_summary(cfg, qc, _result(make_seq, n_reps=4, mode="taxonomic1"),
                        ["a.fasta"])
    assert "## Analysis at a glance" in md
    assert "| Dataset type | Non-segmented" in md
    assert "5,231 input sequences → 4 representative sequences" in md
    # protein alphabet, non-segmented, no concat → single marker substrate.
    assert "single marker protein (amino-acid)" in md
    assert "| Clustering tool | MMseqs2 |" in md
    assert "| Selection mode | `taxonomic1` |" in md


def test_at_a_glance_segmented_and_tree_rows(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["segmented"] = {
        "enabled": True, "virus": "bunya",
        "viruses": {"bunya": {"expected_segments": 3, "segments": ["L", "M", "S"]}},
    }
    cfg["phylo"] = {"tool": "auto", "partition": {"enabled": True},
                    "taxonomy_review": {"enabled": True}}
    md = render_summary(
        cfg, _qc(), _result(make_seq, mode="global"), ["a.fasta"],
        complete_isolates={"iso1": {}}, phylo_ran=True, per_protein_ran=True,
        per_segment_ran=True,
    )
    assert "Segmented virus — 3 segments (L, M, S)" in md
    assert "per-isolate marker-protein concatenation (amino-acid)" in md
    # Whole-genome tree row reflects the partitioned supermatrix substrate.
    assert "partitioned supermatrix of per-marker alignments" in md
    assert "| Per-marker trees |" in md
    assert "| Per-segment NT trees |" in md
    assert "| Taxonomy review | enabled" in md


def test_at_a_glance_names_marker_and_hmm_architecture(make_seq, tmp_path):
    """Clustering substrate + whole-genome tree cells name the actual
    marker(s) and (when the HMM tier ran) the domain architecture, not
    just the substrate kind. Single-marker mode lists declared markers in
    priority order ('first satisfied of …')."""
    cfg = _base_cfg(tmp_path)
    cfg["_hmm_runtime"] = {"active": True}
    cfg["clustering"]["cluster_protein"] = [
        {"name": "Spike", "hmms": ["CoV_S1--CoV_S2", "bCoV_S1_N--CoV_S2"],
         "aliases": ["spike", "surface glycoprotein"]},
        {"name": "Nucleocapsid", "hmms": ["CoV_nucleocap"]},
    ]
    cfg["phylo"] = {"tool": "auto", "partition": {"enabled": True}}
    md = render_summary(cfg, _qc(), _result(make_seq, mode="global"),
                        ["a.fasta"], phylo_ran=True)
    # Marker names present; aliases/synonyms NOT dumped into the cell.
    assert "Spike (HMM: CoV_S1--CoV_S2 OR bCoV_S1_N--CoV_S2)" in md
    assert "Nucleocapsid (HMM: CoV_nucleocap)" in md
    assert "surface glycoprotein" not in md  # synonyms excluded
    assert "first satisfied of Spike" in md  # single-marker priority order


def test_at_a_glance_concat_joins_markers_with_plus(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["_hmm_runtime"] = {"active": True}
    cfg["clustering"]["concatenate_markers"] = True
    cfg["clustering"]["cluster_protein"] = [
        {"name": "Spike", "hmms": ["CoV_S1--CoV_S2"]},
        {"name": "Nucleocapsid", "hmms": ["CoV_nucleocap"]},
    ]
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "Spike (HMM: CoV_S1--CoV_S2) + Nucleocapsid (HMM: CoV_nucleocap)" in md


def test_at_a_glance_omits_architecture_when_hmm_not_used(make_seq, tmp_path):
    """No _hmm_runtime.active → the marker name shows but NOT the HMM
    architecture ('if hmm used')."""
    cfg = _base_cfg(tmp_path)
    cfg["clustering"]["cluster_protein"] = [
        {"name": "Spike", "hmms": ["CoV_S1--CoV_S2"]},
    ]
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "single marker protein (amino-acid) — Spike" in md
    assert "HMM:" not in md.split("## Input")[0]


def test_at_a_glance_segmented_names_per_segment_markers(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["_hmm_runtime"] = {"active": True}
    cfg["segmented"] = {
        "enabled": True, "virus": "bunya",
        "viruses": {"bunya": {
            "expected_segments": 3, "segments": ["L", "M", "S"],
            "segment_markers": {
                "L": {"name": "RdRp", "hmms": ["Bunya_RdRp"]},
                "M": {"name": "GPC", "hmms": ["Bunya_G1--Bunya_G2"]},
                "S": {"name": "N", "hmms": ["Bunya_nucleocap"]},
            },
        }},
    }
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"],
                        complete_isolates={"iso1": {}})
    assert "L: RdRp (HMM: Bunya_RdRp)" in md
    assert "M: GPC (HMM: Bunya_G1--Bunya_G2)" in md
    assert "S: N (HMM: Bunya_nucleocap)" in md


def test_build_provenance_header_one_line(make_seq, tmp_path):
    from repseq.output.summary import build_provenance_header
    cfg = _base_cfg(tmp_path)
    line = build_provenance_header(
        cfg, _result(make_seq, n_reps=7, mode="host"), segmented=False,
    )
    assert line.startswith("# repseq ")
    assert "host selection" in line
    assert "non-segmented" in line
    assert "single marker protein (amino-acid)" in line
    assert "via MMseqs2" in line
    assert "7 representative sequences" in line
    assert "\n" not in line


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
    assert "*1,200* passed this initial screen" in md


def test_render_summary_describes_input_blocklist(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["_excluded_runtime"] = {"audit": [
        {"id": "NC_1.1", "action": "excluded", "detail": "matched on accession"},
        {"id": "GHOST", "action": "unavailable", "detail": "no input sequence matched"},
    ]}
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    # Only the 'excluded' action is counted (not 'unavailable').
    assert "**1** sequence(s) were removed by the user-supplied input blocklist" in md
    assert "`overrides.exclude`" in md
    assert "{prefix}_excluded.tsv" in md


def test_render_summary_omits_blocklist_when_nothing_excluded(make_seq, tmp_path):
    md = render_summary(_base_cfg(tmp_path), _qc(), _result(make_seq), ["a.fasta"])
    assert "input blocklist" not in md


def test_render_summary_diversity_mode_does_not_name_a_clustering_binary(make_seq, tmp_path):
    """global -n (global:count) runs alignment-free MaxMin selection — NO
    clustering binary executes. The summary must not claim/cite that a
    clustering program ran (regression for the cd-hit/MMseqs2 mislabel)."""
    qc = _qc()
    cfg = _base_cfg(tmp_path)
    cfg["clustering"]["backend"] = "cdhit"
    cfg["clustering"]["alphabet_for_clustering"] = "nucleotide"
    result = _result(make_seq, mode="global:count")
    md = render_summary(cfg, qc, result, ["a.fasta"])
    # No false "Clustering was performed ... using **<binary>**" claim/citation.
    assert "Clustering was performed" not in md
    assert "using **cd-hit" not in md
    assert "using **MMseqs2**" not in md
    # Honest MaxMin description instead.
    assert "maximum-diversity sampling" in md
    assert "MaxMin" in md
    # Glance row says no clustering tool.
    assert "| Clustering tool | none — MaxMin diversity selection (alignment-free) |" in md
    # Software table marks BOTH binaries not-used.
    assert md.count("Sequence clustering (not used — MaxMin diversity selection)") == 2
    assert "Sequence clustering (used)" not in md


def test_render_summary_threshold_mode_still_names_the_clustering_tool(make_seq, tmp_path):
    """Regression guard: an actual clustering run (global:threshold) must
    still name and cite the clustering binary it used."""
    qc = _qc()
    cfg = _base_cfg(tmp_path)
    cfg["clustering"]["backend"] = "mmseqs2"
    result = _result(make_seq, mode="global:threshold")
    md = render_summary(cfg, qc, result, ["a.fasta"])
    assert "Clustering was performed" in md
    assert "using **MMseqs2**" in md
    assert "| Clustering tool | MMseqs2 |" in md
    assert "MMseqs2 |" in md and "Sequence clustering (used)" in md


def test_provenance_header_diversity_mode_describes_maxmin(make_seq, tmp_path):
    from repseq.output.summary import build_provenance_header
    result = _result(make_seq, mode="global:count")
    hdr = build_provenance_header(_base_cfg(tmp_path), result, segmented=False)
    assert "via MaxMin diversity (alignment-free)" in hdr
    assert "via MMseqs2" not in hdr and "via cd-hit" not in hdr


def test_render_summary_reports_force_selected(make_seq, tmp_path):
    """When sequences were force-selected, the selection section must say so,
    break down the actions, flag unavailable pins, and point at the TSV."""
    qc = _qc()
    result = _result(make_seq)
    result.force_selected = [
        {"id": "P1", "action": "elected_representative", "detail": "cluster=c1"},
        {"id": "P2", "action": "split_singleton", "detail": "from_cluster=c1"},
        {"id": "D1", "action": "added_representative", "detail": ""},
        {"id": "GHOST", "action": "unavailable", "detail": ""},
    ]
    md = render_summary(_base_cfg(tmp_path), qc, result, ["a.fasta"])
    assert "force-selected" in md
    assert "`overrides.force_select`" in md
    assert "won their cluster's representative slot" in md
    assert "split into singleton clusters" in md
    assert "could not be selected" in md          # the unavailable pin
    assert "{prefix}_force_selected.tsv" in md


def test_render_summary_no_force_select_sentence(make_seq, tmp_path):
    md = render_summary(_base_cfg(tmp_path), _qc(), _result(make_seq), ["a.fasta"])
    assert "force-selected" not in md


def test_render_summary_reports_force_kept_overrides(make_seq, tmp_path):
    """When sequences were force-kept via overrides.protect_qc, the QC
    section must say so, name the stage(s), and point at _overrides.tsv."""
    qc = _qc(total_input=100, passed=98, dedup=0, length=0, ambig=0)
    qc.add_protected("NC_1.1", "ambiguous", "ambiguous_fraction:0.200>0.05")
    qc.add_protected("NC_1.1", "hmm", "hmm_failed:S:Foo")
    qc.add_protected("NC_2.1", "ambiguous", "ambiguous_fraction:0.300>0.05")
    md = render_summary(_base_cfg(tmp_path), qc, _result(make_seq), ["a.fasta"])
    assert "**2** sequence(s) of special importance" in md  # distinct ids
    assert "kept despite failing QC" in md
    assert "`overrides.protect_qc`" in md
    assert "`ambiguous`" in md and "`hmm`" in md
    assert "{prefix}_overrides.tsv" in md


def test_render_summary_no_overrides_no_force_keep_sentence(make_seq, tmp_path):
    qc = _qc()
    md = render_summary(_base_cfg(tmp_path), qc, _result(make_seq), ["a.fasta"])
    assert "kept despite failing QC" not in md


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


def test_render_summary_phylo_surfaces_iqtree_chosen_model(make_seq, tmp_path):
    """Non-partitioned IQ-TREE run: when {prefix}_iqtree_model.txt is on
    disk, the renderer must quote ModelFinder's actual pick rather than
    leaving the generic 'ModelFinder for substitution-model selection'
    line standing alone."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "iqtree", "partition": {"enabled": False}}
    (tmp_path / "test_iqtree_model.txt").write_text("GENOME: LG+I+G4\n")
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "ModelFinder selected **LG+I+G4**" in md


def test_render_summary_phylo_omits_model_clause_when_file_missing(make_seq, tmp_path):
    """No sidecar (e.g. FastTree run, IQ-TREE soft-failed) → no
    pretend pick line."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "iqtree", "partition": {"enabled": False}}
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "ModelFinder selected" not in md


def test_render_summary_phylo_partitioned_lists_per_partition_picks(make_seq, tmp_path):
    """Partitioned IQ-TREE run with a multi-line sidecar: every
    partition's ModelFinder pick should land in the partitioned paragraph."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {
        "tool": "iqtree",
        "partition": {"enabled": True, "linkage": "proportional"},
    }
    (tmp_path / "test_iqtree_model.txt").write_text(
        "CoV_S1: LG+I+G4\nCoV_M: JTT+G4\n"
    )
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "ModelFinder selected:" in md
    assert "CoV_S1=`LG+I+G4`" in md
    assert "CoV_M=`JTT+G4`" in md
    assert "{prefix}_iqtree_model.txt" in md


def test_render_summary_pre_cluster_paragraph_when_run(make_seq, tmp_path):
    """When pre_cluster_ran=True, the phylo section gains a paragraph
    describing the {prefix}_pre_cluster_tree.* outputs with the [repr]
    leaf-prefix convention. The section header appears even if
    pre_cluster_ran is the only flag set."""
    cfg = _base_cfg(tmp_path)
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], pre_cluster_ran=True,
    )
    assert "## Phylogenetic inference" in md
    assert "pre-cluster overview tree" in md
    assert "MAFFT" in md and "`--retree 1`" in md
    assert "FastTree" in md
    assert "midpoint" in md.lower()
    assert "[repr]" in md
    # phylo.newick defaults to false → the .nwk is NOT listed; the phyloXML
    # and id_map are.
    assert "{prefix}_pre_cluster_tree.nwk" not in md
    assert "{prefix}_pre_cluster_tree.xml" in md
    assert "{prefix}_pre_cluster_tree_id_map.tsv" in md
    assert "not retained" in md  # the retention-policy sentence


def test_render_summary_pre_cluster_lists_newick_when_enabled(make_seq, tmp_path):
    """With phylo.newick: true the pre-cluster .nwk is named again."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"newick": True}
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], pre_cluster_ran=True,
    )
    assert "{prefix}_pre_cluster_tree.nwk" in md
    assert "not retained" not in md


def test_render_summary_no_pre_cluster_paragraph_when_off(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True,
    )
    assert "pre-cluster overview tree" not in md


def test_render_summary_msa_conservation_paragraph(make_seq, tmp_path):
    """Any MSA-producing phylo step yields the JSD MSA-conservation
    paragraph by default; disabling the config drops it."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "fasttree"}
    md_on = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True,
    )
    assert "MSA conservation scoring" in md_on
    assert "{prefix}_msa_conservation.tsv" in md_on
    assert "Jensen-Shannon divergence" in md_on
    assert "mean_conservation_core" in md_on
    # Full bibliographic citations for the methods used + the alternative.
    assert "Capra & Singh 2007" in md_on
    assert "Henikoff & Henikoff 1994" in md_on
    assert "Valdar 2002" in md_on

    cfg_off = _base_cfg(tmp_path)
    cfg_off["phylo"] = {"tool": "fasttree", "conservation": {"enabled": False}}
    md_off = render_summary(
        cfg_off, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True,
    )
    assert "MSA conservation scoring" not in md_off


def test_render_summary_monophyly_paragraph(make_seq, tmp_path):
    """A phylo step yields the per-taxon monophyly paragraph; no tree → none."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "fasttree"}
    md_on = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True,
    )
    assert "Per-taxon monophyly" in md_on
    assert "{prefix}_monophyly.tsv" in md_on
    assert "polyphyletic" in md_on
    assert "{prefix}_flags.txt" in md_on
    assert "{prefix}_report.html" in md_on
    assert "support-aware" in md_on
    assert "phylo.monophyly.min_support" in md_on
    # species rank is opt-in: default prose says so and names the knob
    assert "phylo.monophyly.include_species" in md_on
    assert "Species rank is excluded by default" in md_on
    # cross-tree segment-status matrix is named in the monophyly prose
    assert "{prefix}_segment_status_matrix.tsv" in md_on
    assert "single_marker_break" in md_on

    md_off = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "Per-taxon monophyly" not in md_off


def test_render_summary_monophyly_species_on(make_seq, tmp_path):
    """With include_species set, the monophyly prose flips to the included
    wording so a Methods paste reflects what actually ran."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "fasttree", "monophyly": {"include_species": True}}
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True,
    )
    assert "Species-rank rows are included" in md
    assert "pinpoints a reassortant" in md


def test_render_summary_tree_figures_paragraph_on_by_default(make_seq, tmp_path):
    """phylo.pdf defaults to true → the Tree figures paragraph names the
    PDF + PNG outputs and the knob."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "fasttree"}
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "Tree figures" in md
    assert "*_tree.pdf" in md and "*_tree.png" in md
    assert "`phylo.pdf: true`" in md
    assert "--no-pdf" in md


def test_render_summary_tree_figures_paragraph_when_disabled(make_seq, tmp_path):
    """phylo.pdf: false → the paragraph states no figures were written."""
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "fasttree", "pdf": False}
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True)
    assert "Tree figures" in md
    assert "`phylo.pdf: false`" in md
    assert "no `*_tree.pdf`" in md


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


def _poly_cfg(tmp_path, *, whole_tree=False, newick=False):
    cfg = _base_cfg(tmp_path)
    cfg["_hmm_runtime"] = {"active": True}
    cfg["clustering"]["polyprotein"] = [{
        "name": "ORF1ab",
        "peptides": [
            {"name": "nsp3", "hmm": "Macro"},
            {"name": "nsp12", "hmm": "RdRP_1"},
        ],
    }]
    cfg["phylo"] = {
        "tool": "fasttree",
        "newick": newick,
        "per_protein": {"min_taxa": 3, "whole_polyprotein_tree": whole_tree},
    }
    return cfg


def test_render_summary_polyprotein_section_renders(make_seq, tmp_path):
    """Regression: the polyprotein section referenced an out-of-scope
    keep_newick and NameError'd whenever per_protein_ran=True with a spec —
    soft-failed, so the whole _summary.md silently went missing. It must
    render, and name the always-on whole-polyprotein FASTA."""
    cfg = _poly_cfg(tmp_path)
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], per_protein_ran=True,
    )
    assert "## Polyprotein cutting" in md
    assert "tree was also built per peptide" in md
    assert "whole-polyprotein FASTA" in md
    assert "{prefix}_<spec>_polyprotein.fasta" in md


def test_render_summary_whole_polyprotein_tree_prose_when_on(make_seq, tmp_path):
    cfg = _poly_cfg(tmp_path, whole_tree=True)
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], per_protein_ran=True,
    )
    assert "whole-polyprotein tree" in md
    assert "_polyprotein_tree.xml" in md
    assert "across the whole polyprotein end-to-end" in md


def test_render_summary_whole_polyprotein_tree_prose_absent_when_off(make_seq, tmp_path):
    cfg = _poly_cfg(tmp_path, whole_tree=False)
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], per_protein_ran=True,
    )
    assert "## Polyprotein cutting" in md
    assert "whole-polyprotein tree" not in md


def test_render_summary_taxonomy_review_described_when_enabled(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {
        "tool": "fasttree",
        "taxonomy_review": {"enabled": True, "ranks": ["family", "genus", "subgenus"]},
    }
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True,
    )
    assert "Phylogeny-based taxonomy review" in md
    assert "_taxonomy_review.tsv" in md
    assert "never auto-changed" in md


def test_render_summary_taxonomy_review_absent_when_disabled(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["phylo"] = {"tool": "fasttree", "taxonomy_review": {"enabled": False}}
    md = render_summary(
        cfg, _qc(), _result(make_seq), ["a.fasta"], phylo_ran=True,
    )
    assert "Phylogeny-based taxonomy review" not in md


def test_render_summary_no_phylo_section_when_neither_ran(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "## Phylogenetic inference" not in md


def test_render_summary_protein_alphabet_described(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    # protein + non-segmented → "marker-protein sequences"
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "marker-protein sequences" in md


def test_render_summary_concatenate_markers_described(make_seq, tmp_path):
    cfg = _base_cfg(tmp_path)
    cfg["clustering"]["concatenate_markers"] = True
    md = render_summary(cfg, _qc(), _result(make_seq), ["a.fasta"])
    assert "concatenated marker-protein sequences" in md
    assert "concatenate_markers: true" in md


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
