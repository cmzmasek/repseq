"""Per-label distribution reports (``{prefix}_<label>_report.txt`` / ``.tsv``).

The label analogue of one rank of the taxonomic report: each value of an
external label (a Nextclade HA clade, a GenoFLU genotype) counted pre-QC,
in the post-QC pool and among the representatives. Design points pinned here
(agreed with the designer 2026-10-29):

* a pre-QC column, counted as rows of the label's table matched by the RAW
  input (labels are attached only after QC, so removed records never carry
  one) — WITHHELD when any value has more table-matched post-QC members than
  pre-QC rows, since that proves the pre-QC match was incomplete (an
  isolate-name-keyed table: names are read from GenBank after QC);
* every value listed, no top-N cap (a 47-clade run would otherwise hide 27);
* sorted by pool count; ``Unknown`` / ``Conflict`` below the real values;
* "removed by QC" is a DIRECT check — every input record matched to the value
  is in the QC removal log — never an inference from absence, because a value
  can vanish while its isolates survive (rows disagree → ``Conflict``, or
  assembly kept a different duplicate record than the one the table lists).
  Removed values carry their most common removal reasons, and are flagged in
  ``_flags.txt`` (and the HTML report) and on the console.
"""
from __future__ import annotations

from collections import Counter

from repseq.annotate import (
    label_table_sources,
    resolve_annotate_config,
    tally_pre_qc_labels,
)
from repseq.models import QCReport, Sequence, SequenceSource, SequenceType
from repseq.modes.custom_mode import CustomMode, load_metadata_table
from repseq.modes.hybrid_mode import HybridMode
from repseq.output.flags import collect_flags, write_flags_report
from repseq.output.label_report import (
    build_label_tables,
    label_report_stem,
    render_label_report,
    short_reason,
    table_sourced_counts,
    write_label_reports,
)
from repseq.output.summary import _label_report_sentence


# ---------------------------------------------------------------------------
# Builders (production schema)
# ---------------------------------------------------------------------------

def _seq(sid, *, accession=None, isolate_id=None, labels=None):
    s = Sequence(
        id=sid, header=sid, sequence="ACGT" * 10,
        seq_type=SequenceType.NUCLEOTIDE, source=SequenceSource.NCBI,
        accession=accession or sid, isolate_id=isolate_id,
        organism="Influenza A virus", subtype="H5N1",
    )
    if labels is not None:
        s.group_labels = dict(labels)
    return s


def _isolate(name, accs, clade=None):
    """A segmented CONCAT isolate whose segments carry the accessions."""
    segs = [
        Sequence(
            id=f"{name}_{i}", header=f"{name}_{i}", sequence="ACGT" * 10,
            seq_type=SequenceType.NUCLEOTIDE, source=SequenceSource.NCBI,
            accession=acc, segment=str(i + 1), isolate_id=name,
        )
        for i, acc in enumerate(accs)
    ]
    concat = Sequence(
        id=f"CONCAT|{name}", header=name, sequence="ACGT" * 40,
        seq_type=SequenceType.NUCLEOTIDE, source=SequenceSource.NCBI,
        isolate_id=name, organism="Influenza A virus",
    )
    concat.concat_segments = segs
    if clade is not None:
        concat.group_labels = {"clade": clade}
    return concat


def _table(tmp_path, rows, name="clades.tsv", header="accession\tclade"):
    p = tmp_path / name
    p.write_text(header + "\n" + "\n".join(rows) + "\n")
    return p


def _recs(prefix, n):
    """``n`` pre-QC input records, as the tally stores them."""
    return [(f"{prefix}{i}", f"{prefix}{i}") for i in range(n)]


#: Pre-QC tally of the standard scenario: 2.2 is in the input (3 records) and
#: gone from the pool; the other two clades survive.
PRE = {"clade": {
    "source": "clades.tsv",
    "counts": {"2.3.4.4b": 9, "2.3.2.1a": 4, "2.2": 3},
    "records": {
        "2.3.4.4b": _recs("B", 9), "2.3.2.1a": _recs("A", 4),
        "2.2": _recs("G", 3),
    },
}}

#: QC log removing all three 2.2 records.
REMOVED = [
    {"id": "G0", "reason": "incomplete_isolate:missing_segments:1,2"},
    {"id": "G1", "reason": "incomplete_isolate:missing_segments:3"},
    {"id": "G2", "reason": "marker_length_sibling:4:481aa_below_485"},
]


def _pre(counts, records=None):
    return {"clade": {
        "source": "clades.tsv", "counts": dict(counts),
        "records": records or {},
    }}


def _pool():
    return (
        [_seq(f"b{i}", labels={"clade": "2.3.4.4b"}) for i in range(6)]
        + [_seq(f"a{i}", labels={"clade": "2.3.2.1a"}) for i in range(3)]
        + [_seq("u0", labels={"clade": "Unknown"})]
    )


# ---------------------------------------------------------------------------
# Pre-QC tally
# ---------------------------------------------------------------------------

def test_the_pre_qc_tally_counts_matched_table_rows(tmp_path):
    """Unit = distinct table rows reached by the raw input: a row matched by
    two records counts once, a row no input record reaches does not count,
    and a row with no value in the column is skipped. Each value also keeps
    the input records behind it, for the removal-log check."""
    table = load_metadata_table(_table(tmp_path, [
        "PQ1\t2.3.4.4b", "PQ2\t2.3.4.4b", "PQ3\t2.2", "ZZ9\t1", "PQ4\t",
    ]))
    raw = [
        _seq("PQ1.1", accession="PQ1.1"),
        _seq("dup", accession="PQ1"),
        _seq("PQ2", accession="PQ2"),
        _seq("PQ3", accession="PQ3"),
        _seq("PQ4", accession="PQ4"),
        _seq("other", accession="XX1"),
    ]
    out = tally_pre_qc_labels(raw, {"clade": (table, "clade")})
    assert out["clade"]["source"] == "clades.tsv"
    assert out["clade"]["counts"] == {"2.3.4.4b": 2, "2.2": 1}
    assert sorted(out["clade"]["records"]["2.3.4.4b"]) == [
        ("PQ1.1", "PQ1.1"), ("PQ2", "PQ2"), ("dup", "PQ1"),
    ]


def test_one_file_named_twice_is_walked_once(tmp_path, monkeypatch):
    """The mode's --metadata-table and annotate.table can be the same file
    loaded twice (two objects); the raw input is matched against it once."""
    path = _table(tmp_path, ["PQ1\t2.2\tG1"], header="accession\tclade\tgenotype")
    t1, t2 = load_metadata_table(path), load_metadata_table(path)
    calls = []
    import repseq.modes.custom_mode as cm
    real = cm.matched_rows
    monkeypatch.setattr(cm, "matched_rows", lambda s, t: calls.append(1) or real(s, t))
    out = tally_pre_qc_labels(
        [_seq("PQ1")], {"clade": (t1, "clade"), "genotype": (t2, "genotype")},
    )
    assert out["clade"]["counts"] == {"2.2": 1}
    assert out["genotype"]["counts"] == {"G1": 1}
    assert len(calls) == 1


def test_the_sources_come_from_annotate_and_the_selection_field(tmp_path):
    """annotate.fields and a custom --metadata-table field both feed the
    pre-QC tally; the selection field wins a shared key, because its value
    is the one emitted (set_labels with overwrite=True)."""
    ann = _table(tmp_path, ["PQ1\tA"], name="ann.tsv")
    sel = _table(tmp_path, ["PQ1\tB"], name="sel.tsv")
    cfg = {"annotate": {"table": str(ann), "fields": ["clade"]},
           "segmented": {"enabled": False}}
    resolve_annotate_config(cfg)
    mode = CustomMode(cfg, field="clade", n_per_group=1,
                      metadata_table_path=str(sel))
    sources = label_table_sources(cfg, mode.table_label_sources())
    table, column = sources["clade"]
    assert table.path.name == "sel.tsv" and column == "clade"
    assert label_table_sources(cfg)["clade"][0].path.name == "ann.tsv"


def test_a_reserved_selection_field_has_no_label_source(tmp_path):
    """`--field host` is reported where host already is; no label report."""
    sel = _table(tmp_path, ["PQ1\tAvian"], header="accession\thost")
    mode = CustomMode({}, field="host", n_per_group=1,
                      metadata_table_path=str(sel))
    assert label_table_sources({}, mode.table_label_sources()) == {}


def test_hybrid_offers_only_its_table_columns(tmp_path):
    sel = _table(tmp_path, ["PQ1\t2.2\tG1"], header="accession\tclade\tgenotype")
    mode = HybridMode({}, fields=["clade", "genus", "genotype"], n_per_group=1,
                      metadata_table_path=str(sel))
    assert [f for f, _, _ in mode.table_label_sources()] == ["clade", "genotype"]


# ---------------------------------------------------------------------------
# The tables
# ---------------------------------------------------------------------------

def test_counts_at_the_three_stages():
    pool = _pool()
    reps = [pool[0], pool[1], pool[6], pool[9]]
    (t,) = build_label_tables(pool, reps, PRE, qc_details=REMOVED)
    assert t.pool == Counter({"2.3.4.4b": 6, "2.3.2.1a": 3, "Unknown": 1})
    assert t.reps == Counter({"2.3.4.4b": 2, "2.3.2.1a": 1, "Unknown": 1})
    assert t.pre_qc["2.2"] == 3
    assert t.values() == ["2.3.4.4b", "2.3.2.1a", "2.2"]
    assert t.sentinels() == ["Unknown"]
    assert t.eliminated() == [("2.2", 3)]
    assert t.removal_reasons["2.2"] == Counter({
        "incomplete_isolate:missing_segments": 2, "marker_length": 1,
    })


def test_absence_without_qc_removal_is_not_an_alarm():
    """A value can vanish from the pool while its isolates survive: its
    table rows disagree (Conflict), or assembly kept a different duplicate
    record than the one the table lists (duplicates are discarded silently,
    not logged as QC removals). Only two of the three 2.2 records are in the
    log here, so 2.2 is 'not in pool', never 'removed by QC'."""
    (t,) = build_label_tables(_pool(), [], PRE, qc_details=REMOVED[:2])
    assert t.absent() == ["2.2"]
    assert t.eliminated() == []
    assert t.record_check["2.2"] == (3, 2)
    text = render_label_report(t, segmented=True)
    assert "not in pool, not removed by QC" in text
    assert "WARNING" not in text and "Conflict row" in text


def test_without_a_removal_log_nothing_is_called_removed():
    (t,) = build_label_tables(_pool(), [], PRE)
    assert t.absent() == ["2.2"] and t.eliminated() == []


def test_an_incomplete_pre_qc_match_withholds_the_column():
    """More table-matched post-QC members than pre-QC rows is impossible
    when the pre-QC match was complete, so it means the table could not be
    matched before QC. Showing the column would undercount."""
    (t,) = build_label_tables(_pool(), [], _pre({"2.3.4.4b": 2, "2.2": 5}))
    assert t.pre_qc is None
    assert "withheld" in t.pre_qc_note and "2.3.2.1a" in t.pre_qc_note
    assert t.eliminated() == []


def test_a_fallback_value_does_not_withhold_the_column():
    """custom/hybrid: a value from the attribute / taxonomy / --field-regex
    fallback has no table row by definition. Only TABLE-matched members are
    compared with the pre-QC rows, so the column (and the alarm) stay."""
    pool = _pool() + [_seq("r0", labels={"clade": "from_regex"})]
    assignments = (
        [{"stratified_by": "clade", "value": "2.3.4.4b", "matched_by": "table:accession"}] * 6
        + [{"stratified_by": "clade", "value": "2.3.2.1a", "matched_by": "table:accession"}] * 3
        + [{"stratified_by": "clade", "value": "Unknown", "matched_by": "none"}]
        + [{"stratified_by": "clade", "value": "from_regex", "matched_by": "header_regex"}]
    )
    (t,) = build_label_tables(
        pool, [], PRE, table_counts=table_sourced_counts(assignments),
        qc_details=REMOVED,
    )
    assert t.pre_qc is not None
    assert t.eliminated() == [("2.2", 3)]
    # Without the table-sourced counts the regex value would withhold it.
    (t2,) = build_label_tables(pool, [], PRE, qc_details=REMOVED)
    assert t2.pre_qc is None


def test_table_sourced_counts_cover_every_selection_field():
    rows = [{"stratified_by": "HA clade", "value": "2.2", "matched_by": "attribute"}]
    assert table_sourced_counts(rows) == {"ha_clade": Counter()}


def test_a_label_without_a_table_has_no_pre_qc_column():
    (t,) = build_label_tables(_pool(), [], {})
    assert t.pre_qc is None and "do not come from a metadata table" in t.pre_qc_note


def test_a_total_wipeout_still_reports_and_alarms():
    """No pool, so no labelled record: the key comes from the pre-QC tally,
    and a value whose records were all removed is reported — the run where
    the alarm matters most."""
    pre = _pre({"2.2": 3}, {"2.2": _recs("G", 3)})
    (t,) = build_label_tables([], [], pre, qc_details=REMOVED)
    assert t.key == "clade" and t.eliminated() == [("2.2", 3)]


def test_reserved_keys_get_no_report():
    """The same key set the FASTA tags / tree properties use."""
    pool = [_seq("a", labels={"host": "Avian", "clade": "2.2"})]
    assert [t.key for t in build_label_tables(pool, [])] == ["clade"]


def test_natural_order_breaks_ties():
    pool = [_seq(f"x{v}", labels={"clade": v})
            for v in ("2.3.4.10", "2.3.4.4b", "2.3.4.4", "2.3.4.9")]
    (t,) = build_label_tables(pool, [])
    assert t.values() == ["2.3.4.4", "2.3.4.4b", "2.3.4.9", "2.3.4.10"]


def test_short_reason_drops_per_record_detail():
    assert short_reason("marker_length:4:481aa_below_485") == "marker_length"
    assert short_reason("marker_completeness_sibling:1:no_met_start") == "marker_completeness"
    assert short_reason("incomplete_isolate:missing_segments:1,2") == (
        "incomplete_isolate:missing_segments")
    assert short_reason("segmented_filter:could_not_identify_isolate_or_segment") == (
        "segmented_filter:could_not_identify_isolate_or_segment")
    assert short_reason("hmm_failed:4:Hemagglutinin") == "hmm_failed"
    assert short_reason("ambiguous_fraction:0.02>0.01") == "ambiguous_fraction"
    assert short_reason("protein_count_mismatch:segment=4:got=0") == "protein_count_mismatch"
    assert short_reason("exact_duplicate_of:CONCAT|x") == "exact_duplicate_of"


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def test_the_text_report_lists_every_value_and_marks_removals():
    pool = _pool()
    (t,) = build_label_tables(pool, pool[:2], PRE, qc_details=REMOVED)
    text = render_label_report(t, segmented=True, provenance="# prov")
    lines = text.splitlines()
    assert lines[0] == "# prov" and lines[1] == "Label report: clade"
    assert "Source: clades.tsv" in text
    assert "Counting unit: isolates" in text
    assert "Coverage: 9 / 10 pool isolates carry a clade (90.0%)" in text
    assert "pre-QC: 3    pool: 2    reps: 1" in text
    row = next(ln for ln in lines if ln.strip().startswith("2.2 "))
    assert row.split()[:4] == ["2.2", "3", "0", "0"]
    assert "removed by QC (incomplete_isolate:missing_segments 2, marker_length 1)" in row
    unknown = next(ln for ln in lines if ln.strip().startswith("Unknown"))
    assert unknown.split()[:4] == ["Unknown", "-", "1", "0"]
    assert "WARNING: 1 clade value(s)" in text
    # Unknown sits below the real values.
    assert lines.index(unknown) > lines.index(row)


def test_no_value_is_hidden_behind_a_cap():
    pool = [_seq(f"x{i}", labels={"clade": f"c{i}"}) for i in range(47)]
    (t,) = build_label_tables(pool, [])
    text = render_label_report(t, segmented=False)
    assert all(f"  c{i} " in text for i in range(47))
    assert "Counting unit: sequences" in text


def test_a_withheld_column_says_why_and_has_no_pre_qc_cells():
    (t,) = build_label_tables(_pool(), [], _pre({"2.2": 5}))
    text = render_label_report(t, segmented=True)
    assert "NOTE: pre-QC column withheld" in text
    assert "pre-QC" not in text.split("== clade distribution")[1].splitlines()[2]


def test_the_report_file_name_never_collides_with_another_report():
    assert label_report_stem("clade") == "clade"
    assert label_report_stem("taxonomic") == "label_taxonomic"
    assert label_report_stem("protein_taxonomic") == "label_protein_taxonomic"
    # ...nor with another label: 'label_taxonomic' is not 'taxonomic'.
    assert label_report_stem("label_taxonomic") == "label_label_taxonomic"


# ---------------------------------------------------------------------------
# Writing, TSV, flags
# ---------------------------------------------------------------------------

def _write(tmp_path, pre=PRE, details=REMOVED, prefix="flu"):
    pool = _pool()
    return write_label_reports(
        pool, pool[:2], pre, out_dir=tmp_path, prefix=prefix, segmented=True,
        qc_details=details,
    )


def test_both_files_are_written_and_the_alarm_returned(tmp_path):
    paths, elim, status = _write(tmp_path)
    assert [p.name for p in paths] == ["flu_clade_report.txt", "flu_clade_report.tsv"]
    assert elim == [{
        "label": "clade", "value": "2.2", "pre_qc_count": 3,
        "reasons": ["incomplete_isolate:missing_segments", "marker_length"],
    }]
    assert status == {"clade": {"file": "flu_clade_report.txt", "pre_qc": True}}


def test_an_unlabelled_run_writes_nothing(tmp_path):
    pool = [_seq("a")]
    assert write_label_reports(pool, pool, {}, out_dir=tmp_path, prefix="flu",
                               segmented=False) == ([], [], {})
    assert not list(tmp_path.iterdir())


def test_the_tsv_shares_the_tidy_schema(tmp_path):
    _write(tmp_path)
    lines = (tmp_path / "flu_clade_report.tsv").read_text().splitlines()
    assert lines[0].split("\t") == [
        "report", "rank", "pool", "taxon", "taxon_count", "spec", "metric", "value",
    ]
    rows = [dict(zip(lines[0].split("\t"), ln.split("\t"))) for ln in lines[1:]]
    assert {r["report"] for r in rows} == {"label"}
    assert {r["rank"] for r in rows} == {"clade"}
    assert {r["spec"] for r in rows} == {"_label"}
    get = {(r["pool"], r["taxon"], r["metric"]): r["value"] for r in rows}
    assert get[("pre_qc", "2.2", "member_count")] == "3"
    assert get[("pre_qc", "2.2", "records_matched")] == "3"
    assert get[("pre_qc", "2.2", "records_removed_by_qc")] == "3"
    assert get[("post_qc", "2.2", "member_count")] == "0"
    assert get[("post_qc", "*ALL*", "distinct_taxa")] == "2"
    assert get[("post_qc", "Unknown", "member_count")] == "1"
    # Unknown has no pre-QC meaning (most raw segments bind no row).
    assert ("pre_qc", "Unknown", "member_count") not in get
    # Record counts only where the check ran (absent values).
    assert ("pre_qc", "2.3.4.4b", "records_matched") not in get


def test_flags_name_a_clade_removed_by_qc(tmp_path):
    _write(tmp_path)
    flags = [f for f in collect_flags(tmp_path, "flu") if f.category == "label_drop"]
    assert len(flags) == 1
    assert flags[0].message.startswith("2.2 (clade) was present in the input")
    assert "every input record matched to it was removed by QC" in flags[0].message
    assert "_clade_report.txt" in flags[0].message
    # The alarm alone is enough for a flags file (no --phylo needed).
    text = write_flags_report(tmp_path, "flu").read_text()
    assert "## Labelled groups" in text and "2.2 (clade)" in text


def test_absent_but_not_removed_never_raises_the_flag(tmp_path):
    _write(tmp_path, details=REMOVED[:2])
    assert not [f for f in collect_flags(tmp_path, "flu") if f.category == "label_drop"]


def test_a_withheld_column_never_raises_the_flag(tmp_path):
    _write(tmp_path, pre=_pre({"2.2": 5}, {"2.2": _recs("G", 3)}))
    assert not [f for f in collect_flags(tmp_path, "flu") if f.category == "label_drop"]


def test_a_prefix_with_glob_characters_still_finds_the_report(tmp_path):
    _write(tmp_path, prefix="run[1]")
    flags = [f for f in collect_flags(tmp_path, "run[1]") if f.category == "label_drop"]
    assert len(flags) == 1


def test_other_report_tsvs_are_not_read_as_label_reports(tmp_path):
    """A taxonomic report TSV beside the label report must not be mistaken
    for one, even though both match `{prefix}_*_report.tsv`."""
    (tmp_path / "flu_taxonomic_report.tsv").write_text(
        "report\trank\tpool\ttaxon\ttaxon_count\tspec\tmetric\tvalue\n"
        "label\tclade\tpre_qc\tX\t1\t_label\tmember_count\t1\n"
        "label\tclade\tpre_qc\tX\t1\t_label\trecords_matched\t1\n"
        "label\tclade\tpre_qc\tX\t1\t_label\trecords_removed_by_qc\t1\n"
    )
    assert not [f for f in collect_flags(tmp_path, "flu") if f.category == "label_drop"]


# ---------------------------------------------------------------------------
# Methods prose
# ---------------------------------------------------------------------------

def test_the_summary_names_the_report_and_the_removals():
    qc = QCReport()
    qc.display["label_reports"] = {
        "clade": {"file": "flu_clade_report.txt", "pre_qc": True},
    }
    qc.eliminated_labels = [{"label": "clade", "value": "2.2", "pre_qc_count": 3}]
    text = _label_report_sentence(qc, segmented=True)
    assert "`flu_clade_report.txt`" in text
    assert "a leading pre-QC column" in text
    assert "**1** `clade` value(s) present in the input had no post-QC member" in text


def test_the_summary_says_when_there_is_no_pre_qc_column():
    qc = QCReport()
    qc.display["label_reports"] = {
        "genotype": {"file": "flu_genotype_report.txt", "pre_qc": False},
    }
    text = _label_report_sentence(qc, segmented=False)
    assert "No pre-QC column is given for `genotype`" in text


def test_no_labels_no_summary_sentence():
    assert _label_report_sentence(QCReport(), segmented=True) == ""


# ---------------------------------------------------------------------------
# End-to-end shape of the real case
# ---------------------------------------------------------------------------

def test_a_segmented_isolate_counts_once_pre_qc_and_post_qc(tmp_path):
    """An HA-keyed table, raw segment records pre-QC, one CONCAT per isolate
    post-QC; HA3's isolate was removed by QC, so clade 1 is flagged."""
    table = load_metadata_table(_table(tmp_path, ["HA1\t2.2", "HA2\t2.2", "HA3\t1"]))
    raw = [_seq(a, accession=a) for a in ("PB1", "HA1", "PB2", "HA2", "HA3")]
    pre = tally_pre_qc_labels(raw, {"clade": (table, "clade")})
    pool = [_isolate("A/x/1", ["PB1", "HA1"], clade="2.2")]
    details = [{"id": "HA3", "reason": "incomplete_isolate:missing_segments:1"}]
    (t,) = build_label_tables(pool, pool, pre, qc_details=details)
    assert t.pre_qc == Counter({"2.2": 2, "1": 1})
    assert t.pool["2.2"] == 1
    assert t.eliminated() == [("1", 1)]


def test_a_discarded_duplicate_is_not_a_qc_removal(tmp_path):
    """The reviewer's case: an isolate has two HA records, the table lists
    only the shorter one, and assembly keeps the longer one (silently — no
    QC log entry). The isolate survives as Unknown; its clade is absent from
    the pool but must NOT be reported as removed by QC."""
    table = load_metadata_table(_table(tmp_path, ["HAshort\t2.4"]))
    raw = [_seq("HAshort"), _seq("HAlong"), _seq("PB1")]
    pre = tally_pre_qc_labels(raw, {"clade": (table, "clade")})
    pool = [_isolate("A/x/1", ["PB1", "HAlong"], clade="Unknown")]
    (t,) = build_label_tables(pool, pool, pre, qc_details=[])
    assert t.absent() == ["2.4"]
    assert t.eliminated() == []
