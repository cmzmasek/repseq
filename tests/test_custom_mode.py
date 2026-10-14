"""Tests for the custom / hybrid modes' metadata-table grouping.

Covers table loading, how a row binds to a sequence or a segmented isolate
(accession, isolate name, strain, Nextclade seqName; version- and
case-insensitive), the Unknown / Conflict groups, the coverage report and
its zero-match hard error, and {prefix}_group_assignments.tsv.
"""

from __future__ import annotations

import copy

import pytest

from repseq.config import DEFAULTS
from repseq.errors import InputError
from repseq.models import RunResult, Sequence, SequenceSource, SequenceType
from repseq.modes.custom_mode import (
    CONFLICT_GROUP,
    UNKNOWN_GROUP,
    CustomMode,
    load_metadata_table,
    resolve_field,
)
from repseq.modes.hybrid_mode import HybridMode
from repseq.output.report import write_group_assignments_tsv
from repseq.output.summary import render_summary
from repseq.segmented.completeness import (
    _normalise_isolate_id,
    concatenate_isolate,
)


def _seq(sid, *, accession=None, strain=None, isolate_id=None, segment=None,
         host=None, header=None):
    return Sequence(
        id=sid,
        header=header or sid,
        sequence="ACGT" * 10,
        seq_type=SequenceType.NUCLEOTIDE,
        source=SequenceSource.NCBI,
        accession=accession,
        strain=strain,
        isolate_id=isolate_id,
        segment=segment,
        host=host,
    )


def _isolate(name, accessions, *, host=None, headers=None):
    """A segmented isolate built by the production concatenation code."""
    segments = [
        _seq(acc, accession=acc, strain=name, isolate_id=name,
             segment=str(i + 1), host=host,
             header=(headers[i] if headers else acc))
        for i, acc in enumerate(accessions)
    ]
    return concatenate_isolate(segments, _normalise_isolate_id(name))


def _write(tmp_path, text, name="clades.tsv"):
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _cfg(tmp_path, segmented=False):
    cfg = copy.deepcopy(DEFAULTS)
    cfg["output"]["dir"] = str(tmp_path / "out")
    cfg["output"]["prefix"] = "t"
    cfg["segmented"]["enabled"] = segmented
    return cfg


# ---------------------------------------------------------------------------
# Table loading
# ---------------------------------------------------------------------------

def test_load_tab_table_and_case_insensitive_column(tmp_path):
    t = load_metadata_table(_write(tmp_path, "Accession\tClade\nPQ1.1\t2.3.4.4b\n"))
    assert t.key_columns == ["Accession"]
    assert t.column_for("clade") == "Clade"
    assert t.column_for("genotype") is None


def test_load_sniffs_semicolons_nextclade_csv(tmp_path):
    path = _write(tmp_path, "index;seqName;clade\n0;PQ1.1 Influenza A virus;2.3.4.4b\n",
                  name="nextclade.csv")
    t = load_metadata_table(path)
    assert t.key_columns == ["seqName"]
    assert t.rows_for("PQ1") == [0]


def test_load_strips_excel_byte_order_mark(tmp_path):
    path = tmp_path / "x.csv"
    path.write_text("﻿accession,clade\nPQ1,2.3.4.4b\n", encoding="utf-8")
    t = load_metadata_table(path)
    assert t.key_columns == ["accession"]


@pytest.mark.parametrize("text, needle", [
    ("clade\tsubclade\n2.3.4.4b\tx\n", "no column saying which sequence"),
    ("accession\tclade\n", "no data rows"),
    ("", "empty or has no header"),
])
def test_load_errors_are_plain_english(tmp_path, text, needle):
    with pytest.raises(InputError, match=needle):
        load_metadata_table(_write(tmp_path, text))


def test_load_missing_file(tmp_path):
    with pytest.raises(InputError, match="not found"):
        load_metadata_table(tmp_path / "nope.tsv")


# ---------------------------------------------------------------------------
# Binding a row to a sequence / isolate
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("table_acc, seq_acc", [
    ("PQ123456", "PQ123456.1"),
    ("PQ123456.1", "PQ123456"),
    ("pq123456.2", "PQ123456.1"),
])
def test_accession_binding_ignores_version_and_case(tmp_path, table_acc, seq_acc):
    t = load_metadata_table(_write(tmp_path, f"accession\tclade\n{table_acc}\tC1\n"))
    fv = resolve_field(_seq(seq_acc, accession=seq_acc), "clade", t, None)
    assert (fv.value, fv.matched_by) == ("C1", "table:accession")


@pytest.mark.parametrize("seqname", [
    "PQ1.1 Influenza A virus (A/duck/X/1/2024(H5N1)) segment 4",
    "PQ1.1|A/duck/X/1/2024|H5N1",
])
def test_nextclade_seqname_binds_on_leading_identifier(tmp_path, seqname):
    t = load_metadata_table(_write(tmp_path, f"seqName\tclade\n{seqname}\tC1\n"))
    fv = resolve_field(_seq("PQ1.1", accession="PQ1.1"), "clade", t, None)
    assert fv.value == "C1"


def test_segmented_isolate_binds_through_ha_accession(tmp_path):
    """The regression: a concatenated isolate has accession=None, so the old
    lookup put every isolate in 'Unknown'."""
    iso = _isolate("A/duck/X/1/2024", ["PQ1.1", "PQ2.1", "PQ3.1", "PQ4.1"])
    assert iso.accession is None
    t = load_metadata_table(_write(tmp_path, "accession\tclade\nPQ4\t2.3.4.4b\n"))
    fv = resolve_field(iso, "clade", t, None)
    assert (fv.value, fv.matched_by) == ("2.3.4.4b", "table:accession")
    assert fv.detail == "PQ4.1 (segment 4)"


def test_segmented_isolate_binds_by_strain_name_with_whitespace(tmp_path):
    """GenoFLU-style table keyed by isolate name; flu names contain spaces."""
    iso = _isolate("A/Louisiana/12/2024 OR_IR", ["PQ1.1", "PQ2.1"])
    t = load_metadata_table(_write(
        tmp_path, "strain\tgenotype\nA/Louisiana/12/2024 OR_IR\tD1.1\n"))
    fv = resolve_field(iso, "genotype", t, None)
    assert fv.value == "D1.1"
    assert fv.matched_by in ("table:isolate", "table:strain")


def test_segments_with_different_values_go_to_conflict(tmp_path):
    iso = _isolate("A/duck/X/1/2024", ["PQ1.1", "PQ2.1"])
    t = load_metadata_table(_write(tmp_path, "accession\tclade\nPQ1\tA\nPQ2\tB\n"))
    fv = resolve_field(iso, "clade", t, None)
    assert (fv.value, fv.matched_by) == (CONFLICT_GROUP, "conflict")
    assert "PQ1.1 (segment 1)=A" in fv.detail and "PQ2.1 (segment 2)=B" in fv.detail


def test_rows_agreeing_on_the_value_are_not_a_conflict(tmp_path):
    iso = _isolate("A/duck/X/1/2024", ["PQ1.1", "PQ2.1"])
    t = load_metadata_table(_write(
        tmp_path, "accession\tisolate\tclade\nPQ1\tA/duck/X/1/2024\tA\nPQ2\t\tA\n"))
    assert resolve_field(iso, "clade", t, None).value == "A"


def test_blank_value_falls_through_to_unknown(tmp_path):
    t = load_metadata_table(_write(tmp_path, "accession\tclade\nPQ1\t\nPQ9\tC\n"))
    fv = resolve_field(_seq("PQ1", accession="PQ1"), "clade", t, None)
    assert (fv.value, fv.matched_by) == (UNKNOWN_GROUP, "none")


def test_header_regex_on_isolate_reads_segment_headers():
    iso = _isolate("iso1", ["PQ1", "PQ2"],
                   headers=["PQ1 clade=2.3.4.4b", "PQ2 clade=2.3.4.4b"])
    fv = resolve_field(iso, "clade", None, r"clade=(\S+)")
    assert (fv.value, fv.matched_by) == ("2.3.4.4b", "header_regex")


# ---------------------------------------------------------------------------
# Mode behaviour
# ---------------------------------------------------------------------------

def test_custom_field_must_be_a_table_column(tmp_path):
    path = _write(tmp_path, "accession\tClade\nPQ1\tA\n")
    with pytest.raises(InputError, match="not a column"):
        CustomMode(_cfg(tmp_path), field="genotype", n_per_group=5,
                   metadata_table_path=str(path))


def test_custom_run_groups_isolates_and_reports_coverage(tmp_path, capsys):
    isolates = [
        _isolate("iso1", ["A1", "A2"]),
        _isolate("iso2", ["B1", "B2"]),
        _isolate("iso3", ["C1", "C2"]),
        _isolate("iso4", ["D1", "D2"]),
    ]
    path = _write(tmp_path, "accession\tclade\nA2\tX\nB2\tX\nC1\tY\nD1\tY\nD2\tZ\n")
    mode = CustomMode(_cfg(tmp_path, segmented=True), field="clade",
                      n_per_group=10, metadata_table_path=str(path))
    result = mode.run(isolates)

    groups = {g.group: g.n_before for g in result.group_stats}
    assert groups == {"X": 2, "Y": 1, CONFLICT_GROUP: 1}
    out = capsys.readouterr()
    assert "3 of 4 isolates took their value from clades.tsv (75.0%)" in out.out
    assert "1 with conflicting values -> 'Conflict'" in out.out
    assert "WARNING" in out.err
    assert result.config_snapshot["table_coverage"]["clade"] == {
        "n": 4, "table": 3, "conflict": 1, "unknown": 0, "other": 0}
    assert len(result.group_assignments) == 4


def test_custom_run_stops_when_no_row_matches(tmp_path):
    path = _write(tmp_path, "accession\tclade\nZZ1\tX\n")
    mode = CustomMode(_cfg(tmp_path, segmented=True), field="clade",
                      n_per_group=10, metadata_table_path=str(path))
    with pytest.raises(InputError, match="No row of clades.tsv matched") as exc:
        mode.run([_isolate("iso1", ["A1", "A2"])])
    assert "ZZ1" in str(exc.value) and "A1" in str(exc.value)


def test_custom_run_stops_when_matched_column_is_empty(tmp_path):
    path = _write(tmp_path, "accession\tclade\nA1\t\n")
    mode = CustomMode(_cfg(tmp_path), field="clade", n_per_group=10,
                      metadata_table_path=str(path))
    with pytest.raises(InputError, match="column is empty"):
        mode.run([_seq("A1", accession="A1")])


def test_custom_without_table_records_attribute_assignments(tmp_path):
    seqs = [_seq("A1", accession="A1", host="duck"), _seq("A2", accession="A2")]
    result = CustomMode(_cfg(tmp_path), field="host", n_per_group=10).run(seqs)
    by_acc = {a["accessions"]: a for a in result.group_assignments}
    assert (by_acc["A1"]["value"], by_acc["A1"]["matched_by"]) == ("duck", "attribute")
    assert (by_acc["A2"]["value"], by_acc["A2"]["matched_by"]) == (UNKNOWN_GROUP, "none")
    assert "table_coverage" not in result.config_snapshot


def test_hybrid_mixes_table_column_and_attribute(tmp_path, capsys):
    seqs = [_seq("A1", accession="A1", host="duck"),
            _seq("A2", accession="A2", host="chicken")]
    path = _write(tmp_path, "accession\tclade\nA1\tX\nA2\tX\n")
    result = HybridMode(_cfg(tmp_path), fields=["clade", "host"], n_per_group=10,
                        metadata_table_path=str(path)).run(seqs)
    assert {g.group for g in result.group_stats} == {
        "clade=X|host=chicken", "clade=X|host=duck"}
    assert len(result.group_assignments) == 4        # 2 sequences x 2 fields
    assert set(result.config_snapshot["table_coverage"]) == {"clade"}
    assert "[hybrid] clade: 2 of 2 sequences" in capsys.readouterr().out


def test_hybrid_needs_at_least_one_table_column(tmp_path):
    path = _write(tmp_path, "accession\tclade\nA1\tX\n")
    with pytest.raises(InputError, match="None of --fields"):
        HybridMode(_cfg(tmp_path), fields=["host", "country"], n_per_group=5,
                   metadata_table_path=str(path))


# ---------------------------------------------------------------------------
# {prefix}_group_assignments.tsv
# ---------------------------------------------------------------------------

def test_assignments_tsv_columns_and_representative_flag(tmp_path):
    isolates = [_isolate("iso1", ["A1", "A2"]), _isolate("iso2", ["B1", "B2"])]
    path = _write(tmp_path, "accession\tclade\nA1\tX\nB1\tX\n")
    result = CustomMode(_cfg(tmp_path, segmented=True), field="clade",
                        n_per_group=10, metadata_table_path=str(path)).run(isolates)
    # Simulate the final rep set after force-select changed it.
    result.representatives = [isolates[1]]
    out = tmp_path / "t_group_assignments.tsv"
    assert write_group_assignments_tsv(result, out)
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == [
        "isolate_id", "strain", "accessions", "stratified_by", "value",
        "matched_by", "detail", "stratum", "representative"]
    rows = [dict(zip(lines[0].split("\t"), ln.split("\t"))) for ln in lines[1:]]
    assert rows[0]["isolate_id"] == "iso2"           # representatives first
    assert rows[0]["representative"] == "TRUE"
    assert rows[0]["accessions"] == "B1,B2"
    assert rows[1]["representative"] == "FALSE"
    assert rows[1]["detail"] == "A1 (segment 1)"


def test_assignments_tsv_skipped_for_other_modes(tmp_path):
    assert not write_group_assignments_tsv(RunResult(mode="global"), tmp_path / "x.tsv")
    assert not (tmp_path / "x.tsv").exists()


def test_summary_describes_table_coverage(tmp_path):
    isolates = [_isolate("iso1", ["A1", "A2"]), _isolate("iso2", ["B1", "B2"])]
    path = _write(tmp_path, "accession\tclade\nA1\tX\n")
    cfg = _cfg(tmp_path, segmented=True)
    result = CustomMode(cfg, field="clade", n_per_group=10,
                        metadata_table_path=str(path)).run(isolates)
    from repseq.models import QCReport
    md = render_summary(cfg, QCReport(), result, ["a.fasta"])
    assert "metadata table (`clades.tsv`): **1 of 2** isolates matched" in md
    assert "through any of its segments' accessions" in md
    assert "1 without a label formed an `Unknown` group" in md
    assert "The `Unknown` group was sampled like any other group." in md
    assert "`Conflict`" not in md                     # no conflicts this run
    assert "{prefix}_group_assignments.tsv" in md


def test_summary_names_both_fallback_groups_when_present(tmp_path):
    isolates = [_isolate("iso1", ["A1", "A2"]), _isolate("iso2", ["B1", "B2"]),
                _isolate("iso3", ["C1", "C2"])]
    path = _write(tmp_path, "accession\tclade\nA1\tX\nB1\tX\nB2\tY\n")
    cfg = _cfg(tmp_path, segmented=True)
    result = CustomMode(cfg, field="clade", n_per_group=10,
                        metadata_table_path=str(path)).run(isolates)
    from repseq.models import QCReport
    md = render_summary(cfg, QCReport(), result, ["a.fasta"])
    assert "The `Unknown` and `Conflict` groups were sampled" in md


def test_summary_without_table_only_points_at_assignments(tmp_path):
    cfg = _cfg(tmp_path)
    result = CustomMode(cfg, field="host", n_per_group=10).run(
        [_seq("A1", accession="A1", host="duck")])
    from repseq.models import QCReport
    md = render_summary(cfg, QCReport(), result, ["a.fasta"])
    assert "metadata table" not in md
    assert "The label of every sequence, and how it was assigned" in md
