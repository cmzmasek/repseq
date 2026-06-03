"""CLI helpers: the closing run summary / no-output warning."""
from __future__ import annotations

from unittest.mock import MagicMock

import pytest

from repseq.cli import (
    _apply_fast_overrides,
    _check_output_dir,
    _final_summary,
    _handle_segmented,
    _load_and_validate,
    _populate_genbank_isolate_segment,
    _setup_protein_alphabet,
)
from repseq.models import Cluster, QCReport, RunResult, SequenceSource


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
    assert "40 of 100 input sequences passed basic QC" in out


def test_final_summary_appends_final_survivors_when_different(make_seq, capsys):
    """When later QC stages trimmed past basic-QC pass count, the line
    must include the final survivor number so the user isn't misled."""
    result = _result(make_seq, n_reps=2)
    report = QCReport(total_input=100, passed=40)
    report.final_survivors = 12
    report.final_survivors_unit = "isolates"
    _final_summary(result, report, {"segmented": {"enabled": True}})
    out = capsys.readouterr().out
    assert "40 of 100 input sequences passed basic QC" in out
    assert "12 isolates reached selection" in out


def test_final_summary_skips_final_survivors_when_same_as_passed(
    make_seq, capsys
):
    """When no later stage trimmed, don't bother appending — the basic
    QC count is the final count, so the extra clause would be noise."""
    result = _result(make_seq, n_reps=12, n_clusters=4)
    report = QCReport(total_input=100, passed=40)
    report.final_survivors = 40
    report.final_survivors_unit = "sequences"
    _final_summary(result, report, {"segmented": {"enabled": False}})
    out = capsys.readouterr().out
    assert "40 of 100 input sequences passed basic QC." in out
    assert "reached selection" not in out


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
        # These tests exercise nucleotide-level dedup of the concat sequences;
        # protein-alphabet clustering would require marker proteins which the
        # fixture sequences don't carry.
        "clustering": {"alphabet_for_clustering": "nucleotide"},
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


# ---------------------------------------------------------------------------
# _populate_genbank_isolate_segment — GenBank-first, header-regex fallback
# ---------------------------------------------------------------------------

def _cfg_with_toggle(use_genbank_metadata: bool = True) -> dict:
    return {
        "clustering": {"alphabet_for_clustering": "nucleotide"},
        "segmented": {
            "enabled": True,
            "use_genbank_metadata": use_genbank_metadata,
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


def test_populate_genbank_isolate_segment_sets_fields(make_seq):
    a = make_seq("a", "ACGT", source=SequenceSource.NCBI, accession="NC_001.1")
    b = make_seq("b", "ACGT", source=SequenceSource.NCBI, accession="NC_002.1")
    ncbi = MagicMock()
    ncbi.fetch_source_metadata_batch.return_value = {
        "NC_001.1": {"isolate": "SNV-1", "strain": "S1", "segment": "L"},
        "NC_002.1": {"isolate": "SNV-2", "strain": None, "segment": "S"},
    }
    _populate_genbank_isolate_segment([a, b], _cfg_with_toggle(), ncbi)
    assert a.isolate_id == "SNV-1"
    assert a.segment == "L"
    assert a.strain == "S1"
    assert b.isolate_id == "SNV-2"
    assert b.segment == "S"
    assert b.strain is None
    ncbi.fetch_source_metadata_batch.assert_called_once()
    args, _ = ncbi.fetch_source_metadata_batch.call_args
    assert args[0] == ["NC_001.1", "NC_002.1"]


def test_populate_genbank_isolate_segment_falls_back_to_strain(make_seq):
    """When /isolate is absent, /strain takes its place as the isolate id."""
    a = make_seq("a", "ACGT", source=SequenceSource.NCBI, accession="NC_001.1")
    ncbi = MagicMock()
    ncbi.fetch_source_metadata_batch.return_value = {
        "NC_001.1": {"isolate": None, "strain": "Convict Creek 107", "segment": "L"},
    }
    _populate_genbank_isolate_segment([a], _cfg_with_toggle(), ncbi)
    assert a.isolate_id == "Convict Creek 107"
    assert a.segment == "L"


def test_populate_genbank_isolate_segment_preserves_existing_fields(make_seq):
    """If a field is already set, the fetch does not overwrite it."""
    a = make_seq(
        "a", "ACGT", source=SequenceSource.NCBI, accession="NC_001.1",
        isolate_id="prior", segment="prior_seg",
    )
    a.strain = "prior_strain"
    ncbi = MagicMock()
    ncbi.fetch_source_metadata_batch.return_value = {
        "NC_001.1": {"isolate": "NEW", "strain": "NEW_S", "segment": "NEW_SEG"},
    }
    _populate_genbank_isolate_segment([a], _cfg_with_toggle(), ncbi)
    assert a.isolate_id == "prior"
    assert a.segment == "prior_seg"
    assert a.strain == "prior_strain"


def test_populate_genbank_isolate_segment_skips_uniprot_and_no_accession(make_seq):
    """UniProt sequences and sequences without an accession are not fetched."""
    up = make_seq("u", "MEEP", source=SequenceSource.UNIPROT, accession="P12345")
    no_acc = make_seq("n", "ACGT", source=SequenceSource.NCBI)
    no_acc.accession = None
    ncbi = MagicMock()
    _populate_genbank_isolate_segment([up, no_acc], _cfg_with_toggle(), ncbi)
    ncbi.fetch_source_metadata_batch.assert_not_called()


def test_populate_genbank_isolate_segment_no_op_when_toggle_off(make_seq):
    """Toggle set to false bypasses the GenBank fetch entirely."""
    a = make_seq("a", "ACGT", source=SequenceSource.NCBI, accession="NC_001.1")
    ncbi = MagicMock()
    _populate_genbank_isolate_segment([a], _cfg_with_toggle(False), ncbi)
    ncbi.fetch_source_metadata_batch.assert_not_called()
    assert a.isolate_id is None  # untouched — falls through to regex later


def test_populate_genbank_isolate_segment_no_op_without_ncbi(make_seq):
    """--no-resolve leaves ncbi=None — the helper must be a silent no-op."""
    a = make_seq("a", "ACGT", source=SequenceSource.NCBI, accession="NC_001.1")
    _populate_genbank_isolate_segment([a], _cfg_with_toggle(), ncbi=None)
    assert a.isolate_id is None


def test_populate_genbank_isolate_segment_no_op_when_not_segmented(make_seq):
    """Segmented mode disabled → no GenBank fetch, even if ncbi is present."""
    a = make_seq("a", "ACGT", source=SequenceSource.NCBI, accession="NC_001.1")
    ncbi = MagicMock()
    cfg = {"segmented": {"enabled": False}}
    _populate_genbank_isolate_segment([a], cfg, ncbi)
    ncbi.fetch_source_metadata_batch.assert_not_called()


# ---------------------------------------------------------------------------
# _setup_protein_alphabet — auto-trigger fetch + populate protein_sequence
# ---------------------------------------------------------------------------

def _cds(product: str, aa: str):
    return {"protein_id": "P", "product": product, "length": len(aa), "sequence": aa}


def test_setup_protein_alphabet_no_op_when_nucleotide(make_seq):
    a = make_seq("a", "ACGT")
    cfg = {"clustering": {"alphabet_for_clustering": "nucleotide"}}
    out = _setup_protein_alphabet([a], cfg, QCReport(), MagicMock())
    assert out == [a]
    assert a.protein_sequence is None


def test_setup_protein_alphabet_fetches_when_proteins_missing(make_seq):
    a = make_seq("a", "ACGT", source=SequenceSource.NCBI, accession="NC_001.1")
    # Proteins not yet fetched.
    ncbi = MagicMock()
    # attach_proteins will be called; populate proteins via side effect.
    def _attach(seqs, _ncbi):
        for s in seqs:
            s.proteins = [_cds("polymerase", "MMMMM")]
    cfg = {"clustering": {"alphabet_for_clustering": "protein"}, "segmented": {"enabled": False}}
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr("repseq.cli.attach_proteins", _attach)
        out = _setup_protein_alphabet([a], cfg, QCReport(), ncbi)
    assert out == [a]
    assert a.protein_sequence == "MMMMM"


def test_setup_protein_alphabet_aborts_when_no_resolve_and_protein(make_seq):
    a = make_seq("a", "ACGT", source=SequenceSource.NCBI, accession="NC_001.1")
    cfg = {"clustering": {"alphabet_for_clustering": "protein"}, "segmented": {"enabled": False}}
    with pytest.raises(SystemExit) as exc:
        _setup_protein_alphabet([a], cfg, QCReport(), ncbi=None)
    assert exc.value.code == 1


def test_setup_protein_alphabet_concatenates_markers_when_enabled(make_seq):
    """concatenate_markers: true joins one marker CDS per cluster_protein
    spec, in declared order, into seq.protein_sequence."""
    a = make_seq("a", "ACGT", source=SequenceSource.NCBI, accession="NC_001.1")
    a.proteins = [_cds("spike glycoprotein", "S" * 5), _cds("nucleoprotein", "N" * 3)]
    cfg = {
        "clustering": {
            "alphabet_for_clustering": "protein",
            "concatenate_markers": True,
            "cluster_protein": ["spike", "nucleoprotein"],
        },
        "segmented": {"enabled": False},
    }
    out = _setup_protein_alphabet([a], cfg, QCReport(), MagicMock())
    assert out == [a]
    assert a.protein_sequence == "SSSSS" + "NNN"


def test_setup_protein_alphabet_single_marker_when_concat_off(make_seq):
    """Default (flag off): only the first matching spec's marker is used."""
    a = make_seq("a", "ACGT", source=SequenceSource.NCBI, accession="NC_001.1")
    a.proteins = [_cds("spike glycoprotein", "S" * 5), _cds("nucleoprotein", "N" * 3)]
    cfg = {
        "clustering": {
            "alphabet_for_clustering": "protein",
            "cluster_protein": ["spike", "nucleoprotein"],
        },
        "segmented": {"enabled": False},
    }
    out = _setup_protein_alphabet([a], cfg, QCReport(), MagicMock())
    assert a.protein_sequence == "SSSSS"


def test_setup_protein_alphabet_segmented_does_not_set_protein_sequence(make_seq):
    """In segmented mode, _handle_segmented builds the per-isolate concat —
    the per-sequence helper must not touch seq.protein_sequence on segments."""
    a = make_seq("a_S", "ACGT", source=SequenceSource.NCBI, accession="NC_001.1")
    a.proteins = [_cds("nucleoprotein", "NNNN")]
    cfg = {"clustering": {"alphabet_for_clustering": "protein"}, "segmented": {"enabled": True}}
    ncbi = MagicMock()
    out = _setup_protein_alphabet([a], cfg, QCReport(), ncbi)
    assert out == [a]
    assert a.protein_sequence is None


# ---------------------------------------------------------------------------
# _load_and_validate — --alphabet-for-clustering CLI override
# ---------------------------------------------------------------------------

def test_load_and_validate_alphabet_override_applies(tmp_path):
    """--alphabet-for-clustering nucleotide overrides the YAML default."""
    out_dir = tmp_path / "out"
    cfg = _load_and_validate(
        config_path=None,
        output_dir=str(out_dir),
        prefix=None,
        threads=None,
        seed=None,
        alphabet_for_clustering="nucleotide",
    )
    assert cfg["clustering"]["alphabet_for_clustering"] == "nucleotide"


def test_load_and_validate_alphabet_default_unchanged_when_no_override(tmp_path):
    """Without --alphabet-for-clustering the YAML default ('protein') wins."""
    out_dir = tmp_path / "out"
    cfg = _load_and_validate(
        config_path=None,
        output_dir=str(out_dir),
        prefix=None,
        threads=None,
        seed=None,
        alphabet_for_clustering=None,
    )
    assert cfg["clustering"]["alphabet_for_clustering"] == "protein"


def test_load_and_validate_concatenate_markers_override(tmp_path):
    """--concatenate-markers / --no-concatenate-markers overrides the YAML;
    None leaves the default (off) untouched."""
    cfg = _load_and_validate(
        config_path=None, output_dir=str(tmp_path / "o1"), prefix=None,
        threads=None, seed=None, concatenate_markers=True,
    )
    assert cfg["clustering"]["concatenate_markers"] is True

    cfg = _load_and_validate(
        config_path=None, output_dir=str(tmp_path / "o2"), prefix=None,
        threads=None, seed=None, concatenate_markers=False,
    )
    assert cfg["clustering"]["concatenate_markers"] is False

    cfg = _load_and_validate(
        config_path=None, output_dir=str(tmp_path / "o3"), prefix=None,
        threads=None, seed=None, concatenate_markers=None,
    )
    assert cfg["clustering"]["concatenate_markers"] is False  # YAML default


# ---------------------------------------------------------------------------
# --fast: preliminary-run tree pipeline overrides
# ---------------------------------------------------------------------------

def test_apply_fast_overrides_forces_fasttree_and_skips_trimal_partition():
    """--fast forces FastTree, kills trimAl + partitioned supermatrix."""
    cfg = {
        "phylo": {
            "tool": "iqtree",
            "partition": {"enabled": True},
            "trimal": {"enabled": True, "mode": "automated1"},
            "mafft": {"extra_args": ["--maxiterate", "1000"]},
            "per_protein": {
                "trimal": {"enabled": True},
                "mafft": {"extra_args": ["--localpair"]},
            },
        }
    }
    _apply_fast_overrides(cfg)
    assert cfg["phylo"]["tool"] == "fasttree"
    assert cfg["phylo"]["partition"]["enabled"] is False
    assert cfg["phylo"]["trimal"]["enabled"] is False
    assert cfg["phylo"]["per_protein"]["trimal"]["enabled"] is False


def test_apply_fast_overrides_forces_midpoint_rooting():
    """--fast pins midpoint rooting for ALL trees, overriding any config."""
    cfg = {"phylo": {"rooting": {"method": "taxonomy", "outgroup": "AB123456"}}}
    _apply_fast_overrides(cfg)
    assert cfg["phylo"]["rooting"]["method"] == "midpoint"


def test_apply_fast_overrides_sets_midpoint_rooting_from_empty():
    """--fast pins midpoint rooting even when no rooting block is present."""
    cfg = {"phylo": {}}
    _apply_fast_overrides(cfg)
    assert cfg["phylo"]["rooting"]["method"] == "midpoint"


def test_apply_fast_overrides_pins_mafft_retree_one_and_disables_auto():
    """--fast pins MAFFT to '--retree 1' (no --auto) for 2E and 2F."""
    cfg = {"phylo": {}}
    _apply_fast_overrides(cfg)
    assert cfg["phylo"]["mafft"]["extra_args"] == ["--retree", "1"]
    assert cfg["phylo"]["mafft"]["use_auto"] is False
    assert cfg["phylo"]["per_protein"]["mafft"]["extra_args"] == [
        "--retree", "1",
    ]


def test_apply_fast_overrides_replaces_any_user_mafft_extras():
    """An L-INS-i config is wholesale replaced by --retree 1 under --fast."""
    cfg = {
        "phylo": {
            "mafft": {"extra_args": ["--maxiterate", "1000", "--localpair"]},
            "per_protein": {
                "mafft": {"extra_args": ["--maxiterate", "1000", "--localpair"]},
            },
        }
    }
    _apply_fast_overrides(cfg)
    assert cfg["phylo"]["mafft"]["extra_args"] == ["--retree", "1"]
    assert cfg["phylo"]["per_protein"]["mafft"]["extra_args"] == [
        "--retree", "1",
    ]


def test_load_and_validate_fast_flag_applies_overrides(tmp_path):
    """fast=True at the load/validate boundary mutates cfg before run."""
    out_dir = tmp_path / "out"
    cfg = _load_and_validate(
        config_path=None,
        output_dir=str(out_dir),
        prefix=None,
        threads=None,
        seed=None,
        alphabet_for_clustering=None,
        fast=True,
    )
    assert cfg["phylo"]["tool"] == "fasttree"
    assert cfg["phylo"]["partition"]["enabled"] is False
    assert cfg["phylo"]["trimal"]["enabled"] is False
    assert cfg["phylo"]["mafft"]["extra_args"] == ["--retree", "1"]
    assert cfg["phylo"]["mafft"]["use_auto"] is False


def test_load_and_validate_fast_default_off_leaves_yaml_intact(tmp_path):
    """Without --fast the YAML phylo defaults survive."""
    out_dir = tmp_path / "out"
    cfg = _load_and_validate(
        config_path=None,
        output_dir=str(out_dir),
        prefix=None,
        threads=None,
        seed=None,
        alphabet_for_clustering=None,
    )
    # default config ships phylo.tool="auto", partition.enabled=True.
    assert cfg["phylo"]["tool"] == "auto"
    assert cfg["phylo"]["partition"]["enabled"] is True


# ---------------------------------------------------------------------------
# _check_output_dir — refuse to run into a non-empty output directory
# ---------------------------------------------------------------------------

def test_check_output_dir_passes_when_missing(tmp_path):
    cfg = {"output": {"dir": str(tmp_path / "does_not_exist")}}
    _check_output_dir(cfg)  # must not raise


def test_check_output_dir_passes_when_empty(tmp_path):
    cfg = {"output": {"dir": str(tmp_path)}}  # tmp_path exists but is empty
    _check_output_dir(cfg)  # must not raise


def test_check_output_dir_aborts_when_non_empty(tmp_path, capsys):
    (tmp_path / "repseq_representatives.fasta").write_text(">x\nACGT\n")
    cfg = {"output": {"dir": str(tmp_path)}}
    with pytest.raises(SystemExit) as exc:
        _check_output_dir(cfg)
    assert exc.value.code == 1
    err = capsys.readouterr().err
    assert "already exists and is not" in err


def test_check_output_dir_aborts_when_path_is_a_file(tmp_path, capsys):
    file_path = tmp_path / "not_a_dir"
    file_path.write_text("oops")
    cfg = {"output": {"dir": str(file_path)}}
    with pytest.raises(SystemExit) as exc:
        _check_output_dir(cfg)
    assert exc.value.code == 1
    assert "not a directory" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# `repseq stats` — pre-flight input inspection (#5)
# ---------------------------------------------------------------------------

def test_print_input_stats_reports_basics(make_seq, capsys):
    """Stats panel should show count, length distribution, and per-rank
    taxonomy breakdown for taxa actually present in the input."""
    from repseq.cli import _print_input_stats
    from repseq.models import TaxonomyInfo

    seqs = []
    for i in range(5):
        s = make_seq(f"acc{i}", "ACGT" * (250 + 10 * i))
        s.organism = "Foo virus"
        s.host = "human" if i < 3 else None
        s.taxonomy = TaxonomyInfo(genus="Foovirus", family="Fooviridae")
        seqs.append(s)
    _print_input_stats(seqs, top_n=10)
    out = capsys.readouterr().out
    assert "5 sequence(s)" in out
    assert "median" in out
    # Metadata coverage table is present.
    assert "Metadata coverage" in out
    assert "host" in out
    # Taxonomy panel shows the genus we set.
    assert "Foovirus" in out
    assert "Fooviridae" in out


def test_print_input_stats_handles_no_taxonomy(make_seq, capsys):
    from repseq.cli import _print_input_stats

    seqs = [make_seq(f"acc{i}", "ACGT" * 100) for i in range(3)]
    _print_input_stats(seqs, top_n=10)
    out = capsys.readouterr().out
    assert "no taxonomy resolved" in out
