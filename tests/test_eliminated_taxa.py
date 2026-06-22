"""Pre-QC taxonomic column + the "taxa eliminated by QC" silent-drop alarm.

Covers the v0.67.0 addition: a leading pre-QC column on the taxonomic report
(both .txt and tidy .tsv) and the flags-file / console alarm for an entire
genus-or-higher clade present in the input but wiped out by a QC gate (the
Amarillovirales RdRp-marker episode).
"""
from __future__ import annotations

import csv

from repseq.models import Sequence, SequenceType, TaxonomyInfo
from repseq.output.report import (
    _ALARM_RANKS,
    find_eliminated_taxa,
    tally_taxa_by_rank,
    write_taxonomic_report,
    write_taxonomic_report_tsv,
)
from repseq.output.flags import collect_flags, write_flags_report


def _seq(sid: str, genus: str, family: str = "Testviridae") -> Sequence:
    return Sequence(
        id=sid,
        header=sid,
        sequence="A" * 900,
        seq_type=SequenceType.NUCLEOTIDE,
        accession=sid,
        organism="Test virus",
        taxonomy=TaxonomyInfo(genus=genus, family=family),
    )


# ---------------------------------------------------------------------------
# tally_taxa_by_rank / find_eliminated_taxa
# ---------------------------------------------------------------------------

def test_tally_counts_per_rank_and_excludes_blanks():
    pool = [
        _seq("a", "Orthoflavivirus", "Flaviviridae"),
        _seq("b", "Orthoflavivirus", "Flaviviridae"),
        _seq("c", "Pegivirus", "Hepaciviridae"),
        Sequence(  # blank genus → not a taxon
            id="d", header="d", sequence="A" * 10,
            seq_type=SequenceType.NUCLEOTIDE, accession="d",
            taxonomy=TaxonomyInfo(genus="", family="Hepaciviridae"),
        ),
    ]
    tally = tally_taxa_by_rank(pool)
    assert dict(tally["genus"]) == {"Orthoflavivirus": 2, "Pegivirus": 1}
    assert dict(tally["family"]) == {"Flaviviridae": 2, "Hepaciviridae": 2}


def test_find_eliminated_taxa_flags_only_zero_survivors():
    pre = [
        _seq("a", "Orthoflavivirus", "Flaviviridae"),
        _seq("p1", "Pegivirus", "Hepaciviridae"),
        _seq("p2", "Pegivirus", "Hepaciviridae"),
        _seq("h", "Orthohepacivirus", "Hepaciviridae"),
    ]
    post = [s for s in pre if s.taxonomy.genus != "Pegivirus"]  # Pegivirus gone
    elim = find_eliminated_taxa(tally_taxa_by_rank(pre), post)
    # Pegivirus (genus) eliminated; family Hepaciviridae survives via
    # Orthohepacivirus, so it must NOT be flagged.
    assert elim == [{"rank": "genus", "taxon": "Pegivirus", "pre_qc_count": 2}]


def test_find_eliminated_taxa_ignores_species_and_subgenus():
    # species/subgenus are below the alarm threshold and must never flag.
    assert "species" not in _ALARM_RANKS
    assert "subgenus" not in _ALARM_RANKS
    pre = [Sequence(
        id="x", header="x", sequence="A" * 10,
        seq_type=SequenceType.NUCLEOTIDE, accession="x",
        taxonomy=TaxonomyInfo(species="Gone virus", genus="Survivor"),
    ), _seq("y", "Survivor")]
    post = [_seq("y", "Survivor")]
    elim = find_eliminated_taxa(tally_taxa_by_rank(pre), post)
    # "Gone virus" disappeared at species rank but that rank is not alarmed.
    assert all(e["rank"] != "species" for e in elim)
    assert elim == []  # genus Survivor survives


def test_find_eliminated_taxa_empty_when_no_pre_qc():
    assert find_eliminated_taxa(None, [_seq("a", "X")]) == []
    assert find_eliminated_taxa({}, [_seq("a", "X")]) == []


# ---------------------------------------------------------------------------
# write_taxonomic_report — three-column .txt
# ---------------------------------------------------------------------------

def test_txt_three_columns_show_eliminated_taxon(tmp_path):
    pre_taxa = tally_taxa_by_rank([
        _seq("a", "Orthoflavivirus"), _seq("b", "Orthoflavivirus"),
        _seq("p", "Pegivirus"),
    ])
    post = [_seq("a", "Orthoflavivirus"), _seq("b", "Orthoflavivirus")]
    out = tmp_path / "t_taxonomic_report.txt"
    write_taxonomic_report(
        post, post[:1], segmented=False, path=out,
        pre_qc_by_rank=pre_taxa, pre_qc_total=3,
    )
    txt = out.read_text()
    assert "pre-QC" in txt and "pre-clust" in txt and "reps" in txt
    # The eliminated genus appears with pre-QC 1 and zeros afterwards.
    line = next(l for l in txt.splitlines() if l.strip().startswith("Pegivirus"))
    assert line.split() == ["Pegivirus", "1", "0", "0"]


def test_txt_eliminated_taxon_not_truncated_away(tmp_path):
    # Many surviving genera + one eliminated one with a high pre-QC count.
    pre = [_seq(f"s{i}", f"Genus{i}") for i in range(30)]
    pre += [_seq(f"g{j}", "Doomed") for j in range(99)]  # big pre-QC count
    post = [_seq(f"s{i}", f"Genus{i}") for i in range(30)]  # Doomed gone
    out = tmp_path / "t_taxonomic_report.txt"
    write_taxonomic_report(
        post, post, segmented=False, path=out, max_breakdown=5,
        pre_qc_by_rank=tally_taxa_by_rank(pre), pre_qc_total=len(pre),
    )
    txt = out.read_text()
    # Sorted by pre-QC count: Doomed (99) is the largest, so it survives the
    # top-5 truncation despite having zero post-QC members.
    assert "Doomed" in txt


def test_txt_low_count_eliminated_retained_past_cap(tmp_path):
    # 25 surviving genera each with 10 records (all sort above the eliminated
    # one) + one eliminated genus with a SINGLE pre-QC record — it would fall
    # outside the top-20, but the alarm-rank retention keeps it.
    pre = [_seq(f"s{i:02d}-{k}", f"Genus{i:02d}") for i in range(25) for k in range(10)]
    pre += [_seq("d", "Doomed")]  # pre-QC count 1, far below the cap
    post = [_seq(f"s{i:02d}", f"Genus{i:02d}") for i in range(25)]  # Doomed gone
    out = tmp_path / "t_taxonomic_report.txt"
    write_taxonomic_report(
        post, post, segmented=False, path=out, max_breakdown=20,
        pre_qc_by_rank=tally_taxa_by_rank(pre), pre_qc_total=len(pre),
    )
    txt = out.read_text()
    assert "Doomed" in txt          # retained despite being past the cap
    assert "eliminated" in txt      # the label's "+ N eliminated" note
    # A surviving genus that fell past the cap is NOT force-retained.
    line = next(l for l in txt.splitlines() if l.strip().startswith("Doomed"))
    assert line.split() == ["Doomed", "1", "0", "0"]


def test_txt_segmented_note_warns_about_segment_units(tmp_path):
    pre_taxa = tally_taxa_by_rank([_seq("a", "Hantavirus"), _seq("b", "Hantavirus")])
    out = tmp_path / "t_taxonomic_report.txt"
    write_taxonomic_report(
        [_seq("a", "Hantavirus")], [_seq("a", "Hantavirus")],
        segmented=True, path=out, pre_qc_by_rank=pre_taxa, pre_qc_total=2,
    )
    txt = out.read_text()
    assert "input segments before isolate grouping" in txt
    assert "raw segments" in txt


def test_txt_two_column_backcompat_unchanged(tmp_path):
    out = tmp_path / "t_taxonomic_report.txt"
    write_taxonomic_report(
        [_seq("a", "X")], [_seq("a", "X")], segmented=False, path=out,
    )
    txt = out.read_text()
    assert "before vs after clustering" in txt
    assert "pre-QC" not in txt


# ---------------------------------------------------------------------------
# write_taxonomic_report_tsv — pre_qc pool rows
# ---------------------------------------------------------------------------

def test_tsv_emits_pre_qc_pool_rows(tmp_path):
    pre_taxa = tally_taxa_by_rank([
        _seq("a", "Orthoflavivirus"), _seq("p", "Pegivirus"),
    ])
    post = [_seq("a", "Orthoflavivirus")]
    out = tmp_path / "t_taxonomic_report.tsv"
    write_taxonomic_report_tsv(post, post, path=out, pre_qc_by_rank=pre_taxa)
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    pools = {r["pool"] for r in rows}
    assert pools == {"pre_qc", "post_qc", "reps"}
    # Pegivirus member_count: pre_qc=1, post_qc=0.
    peg = {
        r["pool"]: r["value"]
        for r in rows
        if r["taxon"] == "Pegivirus" and r["metric"] == "member_count"
    }
    assert peg["pre_qc"] == "1" and peg["post_qc"] == "0"


def test_tsv_without_pre_qc_has_no_pre_qc_rows(tmp_path):
    out = tmp_path / "t_taxonomic_report.tsv"
    write_taxonomic_report_tsv([_seq("a", "X")], [_seq("a", "X")], path=out)
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert "pre_qc" not in {r["pool"] for r in rows}


# ---------------------------------------------------------------------------
# flags — the alarm, fired off the TSV
# ---------------------------------------------------------------------------

_TIDY = ("report", "rank", "pool", "taxon", "taxon_count", "spec", "metric", "value")


def _write_tax_tsv(path, rows):
    lines = ["\t".join(_TIDY)] + ["\t".join(map(str, r)) for r in rows]
    path.write_text("\n".join(lines) + "\n")


def test_flags_alarm_fires_without_any_conflict_table(tmp_path):
    # Only the taxonomic report TSV exists (a plain clustering run, no --phylo).
    _write_tax_tsv(tmp_path / "x_taxonomic_report.tsv", [
        ("diversity", "genus", "pre_qc", "Pegivirus", 390, "_diversity", "member_count", 390),
        ("diversity", "genus", "post_qc", "Pegivirus", 0, "_diversity", "member_count", 0),
        ("diversity", "genus", "reps", "Pegivirus", 0, "_diversity", "member_count", 0),
        ("diversity", "genus", "pre_qc", "Orthoflavivirus", 100, "_diversity", "member_count", 100),
        ("diversity", "genus", "post_qc", "Orthoflavivirus", 80, "_diversity", "member_count", 80),
        ("diversity", "genus", "reps", "Orthoflavivirus", 5, "_diversity", "member_count", 5),
    ])
    flags = collect_flags(tmp_path, "x")
    qc = [f for f in flags if f.category == "qc_drop"]
    assert len(qc) == 1
    assert "Pegivirus" in qc[0].message and "390" in qc[0].message
    # The message points at _qc_removed.tsv and does NOT presume the cause.
    assert "_qc_removed.tsv" in qc[0].message
    assert "often an HMM marker that does not cover" not in qc[0].message
    # And the file is written even though no conflict table exists.
    path = write_flags_report(tmp_path, "x")
    assert path is not None
    body = path.read_text()
    assert "Taxa eliminated entirely by QC" in body
    # The QC-drop section comes first (most prominent).
    assert body.index("Taxa eliminated entirely by QC") < len(body)


def test_flags_no_alarm_when_taxon_survives(tmp_path):
    _write_tax_tsv(tmp_path / "x_taxonomic_report.tsv", [
        ("diversity", "genus", "pre_qc", "Survivor", 10, "_diversity", "member_count", 10),
        ("diversity", "genus", "post_qc", "Survivor", 4, "_diversity", "member_count", 4),
        ("diversity", "genus", "reps", "Survivor", 1, "_diversity", "member_count", 1),
    ])
    assert not [f for f in collect_flags(tmp_path, "x") if f.category == "qc_drop"]
    # No conflict tables and no elimination → no file at all.
    assert write_flags_report(tmp_path, "x") is None


def test_flags_alarm_ignores_species_rank(tmp_path):
    _write_tax_tsv(tmp_path / "x_taxonomic_report.tsv", [
        ("diversity", "species", "pre_qc", "Some virus", 5, "_diversity", "member_count", 5),
        ("diversity", "species", "post_qc", "Some virus", 0, "_diversity", "member_count", 0),
    ])
    assert not [f for f in collect_flags(tmp_path, "x") if f.category == "qc_drop"]


def test_console_alarm_fires_even_with_zero_reps(capsys):
    from repseq.cli import _final_summary
    from repseq.models import QCReport, RunResult
    qc = QCReport(total_input=5, passed=5)
    qc.eliminated_taxa = [{"rank": "genus", "taxon": "Pegivirus", "pre_qc_count": 3}]
    result = RunResult(mode="global", representatives=[], clusters=[])
    _final_summary(result, qc, {"segmented": {"enabled": False}})
    err = capsys.readouterr().err
    # No reps selected, yet the silent-drop alarm must still appear.
    assert "ELIMINATED by QC" in err and "Pegivirus" in err


def test_console_alarm_fires_with_reps(capsys):
    from repseq.cli import _final_summary
    from repseq.models import QCReport, RunResult, Sequence, SequenceType
    qc = QCReport(total_input=5, passed=5)
    qc.eliminated_taxa = [{"rank": "family", "taxon": "Pegiviridae", "pre_qc_count": 7}]
    rep = Sequence(id="r", header="r", sequence="A" * 9, seq_type=SequenceType.NUCLEOTIDE)
    result = RunResult(mode="global", representatives=[rep], clusters=[])
    _final_summary(result, qc, {"segmented": {"enabled": False}})
    err = capsys.readouterr().err
    assert "ELIMINATED by QC" in err and "Pegiviridae" in err


def test_html_report_fires_on_qc_drop_without_figures(tmp_path):
    from repseq.output.html_report import write_html_report
    _write_tax_tsv(tmp_path / "x_taxonomic_report.tsv", [
        ("diversity", "genus", "pre_qc", "Pegivirus", 390, "_diversity", "member_count", 390),
        ("diversity", "genus", "post_qc", "Pegivirus", 0, "_diversity", "member_count", 0),
    ])
    # No conflict tables, no tree figures — only the QC-elimination.
    path = write_html_report(tmp_path, "x")
    assert path is not None
    html = path.read_text()
    assert "Taxa eliminated entirely by QC" in html and "Pegivirus" in html
