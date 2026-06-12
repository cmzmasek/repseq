"""Tests for the subtype (serotype) distribution report.

`{prefix}_subtype_report.txt` / `.tsv` mirror one rank of the taxonomic
report but over `seq.subtype` (the viral serotype, e.g. influenza H5N1),
and are emitted only when the representatives carry more than one distinct
subtype.
"""

from pathlib import Path

from repseq.models import Sequence
from repseq.output.report import (
    write_subtype_report,
    write_subtype_report_tsv,
)


def _seq(sid: str, subtype) -> Sequence:
    return Sequence(
        id=sid,
        header=">" + sid,
        sequence="ACGT",
        accession=sid,
        organism="Influenza A virus",
        subtype=subtype,
    )


def _pool():
    return (
        [_seq(f"a{i}", "H3N2") for i in range(40)]
        + [_seq(f"b{i}", "H1N1") for i in range(30)]
        + [_seq(f"c{i}", "H5N1") for i in range(8)]
        + [_seq("d0", "H7N9")]
        + [_seq(f"u{i}", None) for i in range(5)]
    )


def _reps():
    return [
        _seq("r1", "H3N2"),
        _seq("r2", "H1N1"),
        _seq("r3", "H5N1"),
        _seq("r4", "H7N9"),
        _seq("r5", None),
    ]


# ---------------------------------------------------------------------------
# Gating: only when reps carry >1 distinct subtype
# ---------------------------------------------------------------------------

def test_subtype_report_skipped_when_single_subtype_in_reps(tmp_path):
    pool = _pool()
    reps = [_seq("r1", "H3N2"), _seq("r2", "H3N2")]
    txt = tmp_path / "x_subtype_report.txt"
    tsv = tmp_path / "x_subtype_report.tsv"
    assert write_subtype_report(pool, reps, segmented=True, path=txt) is False
    assert write_subtype_report_tsv(pool, reps, path=tsv) is False
    assert not txt.exists()
    assert not tsv.exists()


def test_subtype_report_skipped_when_no_subtype_at_all(tmp_path):
    pool = [_seq(f"u{i}", None) for i in range(10)]
    reps = [_seq("r1", None), _seq("r2", None)]
    txt = tmp_path / "x_subtype_report.txt"
    assert write_subtype_report(pool, reps, segmented=False, path=txt) is False
    assert not txt.exists()


def test_subtype_report_emitted_when_multiple_subtypes_in_reps(tmp_path):
    txt = tmp_path / "x_subtype_report.txt"
    assert write_subtype_report(_pool(), _reps(), segmented=True, path=txt) is True
    assert txt.exists()


# ---------------------------------------------------------------------------
# .txt content
# ---------------------------------------------------------------------------

def test_subtype_report_txt_distribution_and_coverage(tmp_path):
    txt = tmp_path / "x_subtype_report.txt"
    write_subtype_report(_pool(), _reps(), segmented=True, path=txt)
    out = txt.read_text()
    assert "Subtype report" in out
    assert "distinct subtypes    before: 4    after: 4" in out
    # before/after columns for the top subtype
    assert "H3N2" in out and "40" in out
    # coverage line counts only populated subtypes (79 of 84 pool, 4 of 5 reps)
    assert "79 / 84 pool isolates carry a serotype" in out
    assert "4 / 5 representatives" in out
    # the unassigned share is NOT a subtype row
    assert "(unassigned)" not in out
    assert "None" not in out


def test_subtype_report_unit_label_tracks_segmented(tmp_path):
    seg = tmp_path / "seg.txt"
    nonseg = tmp_path / "nonseg.txt"
    write_subtype_report(_pool(), _reps(), segmented=True, path=seg)
    write_subtype_report(_pool(), _reps(), segmented=False, path=nonseg)
    assert "Counting unit: isolates" in seg.read_text()
    assert "Counting unit: sequences" in nonseg.read_text()


def test_subtype_report_excludes_unknown_tokens(tmp_path):
    # "unknown" / "N/A" are not real subtypes — only H3N2 and H1N1 count, so
    # reps still trip the >1 gate but the unknowns never appear.
    pool = [_seq("a", "H3N2"), _seq("b", "H1N1"), _seq("c", "unknown"), _seq("d", "N/A")]
    reps = [_seq("a", "H3N2"), _seq("b", "H1N1"), _seq("c", "unknown")]
    txt = tmp_path / "x.txt"
    assert write_subtype_report(pool, reps, segmented=False, path=txt) is True
    out = txt.read_text()
    assert "distinct subtypes    before: 2    after: 2" in out
    assert "unknown" not in out
    assert "N/A" not in out


def test_subtype_report_truncates_to_max_breakdown(tmp_path):
    pool = [_seq(f"s{i}", f"H{i}N1") for i in range(30)]
    reps = [_seq(f"r{i}", f"H{i}N1") for i in range(30)]
    txt = tmp_path / "x.txt"
    write_subtype_report(pool, reps, segmented=False, path=txt, max_breakdown=20)
    out = txt.read_text()
    assert "top 20 of 30 shown" in out
    assert "+10 more subtypes not shown" in out


# ---------------------------------------------------------------------------
# .tsv schema
# ---------------------------------------------------------------------------

def test_subtype_report_tsv_schema_and_rows(tmp_path):
    tsv = tmp_path / "x_subtype_report.tsv"
    assert write_subtype_report_tsv(_pool(), _reps(), path=tsv) is True
    rows = [ln.split("\t") for ln in tsv.read_text().splitlines()]
    header = rows[0]
    assert header == [
        "report", "rank", "pool", "taxon", "taxon_count", "spec", "metric", "value",
    ]
    body = rows[1:]
    # report / spec columns are constant
    assert all(r[0] == "subtype" for r in body)
    assert all(r[5] == "_subtype" for r in body)
    # distinct_taxa rows for both pools
    assert ["subtype", "subtype", "post_qc", "*ALL*", "4", "_subtype",
            "distinct_taxa", "4"] in body
    assert ["subtype", "subtype", "reps", "*ALL*", "4", "_subtype",
            "distinct_taxa", "4"] in body
    # member_count rows carry the before/after counts
    assert ["subtype", "subtype", "post_qc", "H3N2", "40", "_subtype",
            "member_count", "40"] in body
    assert ["subtype", "subtype", "reps", "H3N2", "1", "_subtype",
            "member_count", "1"] in body


def test_subtype_report_tsv_skipped_when_single_subtype(tmp_path):
    reps = [_seq("r1", "H3N2"), _seq("r2", "H3N2")]
    tsv = tmp_path / "x.tsv"
    assert write_subtype_report_tsv(_pool(), reps, path=tsv) is False
    assert not tsv.exists()
