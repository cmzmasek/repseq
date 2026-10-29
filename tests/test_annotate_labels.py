"""External per-sequence labels (``annotate.table`` / ``--annotate-table``).

A clade or genotype call lives in a table some other tool wrote (Nextclade,
GenoFLU), and until now it died inside the ``custom`` / ``hybrid`` mode that
stratified on it: the label never reached a FASTA header or a tree. These
tests pin the three surfaces that now carry it — protein FASTA bracket tags,
phyloXML ``repseq:<field>`` properties, and the monophyly report's non-rank
groupings — plus the binding rules.

Design points pinned here (agreed with the designer 2026-10-29):

* labels are available to EVERY mode through ``annotate.table``, because
  labelling records is independent of choosing them (a global panel run wants
  clade tags too), while ``custom --field clade`` still needs no second table;
* the NUCLEOTIDE FASTA headers are deliberately untouched — they pass the
  GenBank header through verbatim and anything parsing them keeps working;
* tree LABELS are untouched too: the clade reaches the trees as a
  ``<property>``, so no figure or MSA header changes shape;
* ``phylo.monophyly.extra_keys`` assesses such a label (or ``subtype``) as a
  grouping, which is the only way the monophyly report and the segment-status
  matrix say anything about a flat-taxonomy virus like influenza A.

Binding is the shared ``custom_mode`` implementation on purpose: a segmented
CONCAT isolate has no accession of its own, so a table keyed by HA accession
must reach it through a segment — the bug class the ``_norm_id`` invariant
records. Table lookup is the ONLY source here (no attribute / taxonomy /
regex fallback), so a tag always means "this is what the table said".
"""
from __future__ import annotations

import json

import pytest

from repseq.annotate import (
    ANNOTATE_RUNTIME,
    annotate_sequences,
    annotate_stats,
    resolve_annotate_config,
    set_labels,
)
from repseq.errors import InputError
from repseq.lockfile import _representative_entry
from repseq.models import (
    RESERVED_LABEL_KEYS,
    QCReport,
    RunResult,
    Sequence,
    SequenceSource,
    SequenceType,
    TaxonomyInfo,
    sanitize_label_key,
)
from repseq.output.report import (
    _TAX_RANKS,
    _write_protein_fasta_record,
    _write_sliced_peptide_record,
)
from repseq.output.segment_status import write_segment_status_matrix
from repseq.output.summary import _render_selection
from repseq.phylo.monophyly import group_sort_index, write_monophyly_report
from repseq.phylo.phyloxml_writer import write_phyloxml


# ---------------------------------------------------------------------------
# Builders (production schema)
# ---------------------------------------------------------------------------

def _seq(sid, *, accession=None, isolate_id=None, strain=None, clade=None):
    s = Sequence(
        id=sid, header=sid, sequence="ACGT" * 10,
        seq_type=SequenceType.NUCLEOTIDE, source=SequenceSource.NCBI,
        accession=accession or sid, isolate_id=isolate_id, strain=strain,
        organism="Influenza A virus", host="Avian", subtype="H5N1",
        taxonomy=TaxonomyInfo(genus="Alphainfluenzavirus", species="Influenza A virus"),
    )
    if clade is not None:
        s.group_labels = {"clade": clade}
    return s


def _prot(pid, product="hemagglutinin", n=566):
    return {
        "protein_id": pid, "product": product, "length": n,
        "sequence": "M" + "A" * (n - 1), "hmm_hits": [],
    }


def _isolate(name, accs):
    """A segmented CONCAT isolate whose segments carry the accessions."""
    segs = [
        Sequence(
            id=f"{name}_{i}", header=f"{name}_{i}", sequence="ACGT" * 10,
            seq_type=SequenceType.NUCLEOTIDE, source=SequenceSource.NCBI,
            accession=acc, segment=str(i + 1), isolate_id=name,
            proteins=[_prot(f"{acc}_p")],
        )
        for i, acc in enumerate(accs)
    ]
    concat = Sequence(
        id=f"CONCAT|{name}", header=name, sequence="ACGT" * 40,
        seq_type=SequenceType.NUCLEOTIDE, source=SequenceSource.NCBI,
        isolate_id=name, organism="Influenza A virus", subtype="H5N1",
        taxonomy=TaxonomyInfo(genus="Alphainfluenzavirus"),
    )
    concat.concat_segments = segs
    return concat


def _table(tmp_path, rows, name="clades.tsv", header="accession\tclade"):
    p = tmp_path / name
    p.write_text(header + "\n" + "\n".join(rows) + "\n")
    return p


def _cfg(tmp_path, table, fields=("clade",), segmented=False):
    cfg = {
        "annotate": {"table": str(table), "fields": list(fields)},
        "segmented": {"enabled": segmented},
    }
    resolve_annotate_config(cfg)
    return cfg


# ---------------------------------------------------------------------------
# The tag/property key
# ---------------------------------------------------------------------------

def test_a_column_name_becomes_one_stable_token():
    """A table column is whatever the user's tool wrote, but a FASTA tag key
    and an XML property ref need one safe spelling."""
    assert sanitize_label_key("clade") == "clade"
    assert sanitize_label_key("HA clade") == "ha_clade"
    assert sanitize_label_key("Nextclade.clade") == "nextclade_clade"
    # A leading digit would be an illegal QName start.
    assert sanitize_label_key("2024 genotype") == "genotype"
    assert sanitize_label_key("   ") == ""


def test_label_items_keeps_order_and_always_refuses_reserved_keys():
    """A label must never shadow metadata repseq itself emits.

    The skip list used to be derived per surface from the tags that surface
    happened to build, which meant `--annotate-field species` landed a table
    value in the `[species=...]` slot of the NCBI 9-rank ladder on exactly
    the records whose taxonomy had failed to resolve, and
    `--annotate-field strain` reached the FASTA header but never the tree
    (2026-10-29 review). One reserved set, enforced here, for every surface.
    """
    s = _seq("a")
    s.group_labels = {
        "clade": "2.3.4.4b", "Host": "avian", "species": "x", "empty": "",
    }
    assert s.label_items() == [("clade", "2.3.4.4b")]
    # `skip` only ADDS to the reserved set.
    assert s.label_items(skip=["clade"]) == []
    assert _seq("b").label_items() == []


def test_the_reserved_set_covers_what_the_emitters_actually_write():
    """Drift guard: adding a bracket tag or a leaf property without adding
    it to RESERVED_LABEL_KEYS would silently reopen the shadowing hole."""
    from repseq.phylo.phyloxml_writer import _LEAF_PROPERTIES

    for key in _TAX_RANKS:
        assert key in RESERVED_LABEL_KEYS, key
    for _ref, key in _LEAF_PROPERTIES:
        assert key in RESERVED_LABEL_KEYS, key
    # The fixed tag keys of both FASTA writers.
    for key in (
        "organism", "ncbi_taxon_id", "isolate", "segment", "host", "subtype",
        "country", "collection_date", "length", "hmmscan", "parent",
        "polyprotein", "peptide", "peptide_range_aa", "cut_method",
        "matched_arch", "is_representative",
    ):
        assert key in RESERVED_LABEL_KEYS, key


# ---------------------------------------------------------------------------
# Startup validation (fails in seconds, and under --dry-run)
# ---------------------------------------------------------------------------

def test_no_annotate_config_is_a_no_op():
    cfg = {}
    resolve_annotate_config(cfg)
    assert ANNOTATE_RUNTIME not in cfg


def test_fields_without_a_table_is_an_error():
    with pytest.raises(InputError, match="annotate.table is not set"):
        resolve_annotate_config({"annotate": {"fields": ["clade"]}})


def test_a_table_without_fields_is_an_error(tmp_path):
    t = _table(tmp_path, ["PQ1\t2.3.4.4b"])
    with pytest.raises(InputError, match="annotate.fields is"):
        resolve_annotate_config({"annotate": {"table": str(t)}})


def test_a_missing_table_names_the_path():
    with pytest.raises(InputError, match="annotate.table not found"):
        resolve_annotate_config(
            {"annotate": {"table": "/nope/x.tsv", "fields": ["clade"]}}
        )


def test_an_unknown_column_lists_the_real_ones(tmp_path):
    t = _table(tmp_path, ["PQ1\t2.3.4.4b"])
    with pytest.raises(InputError, match="not a column of"):
        resolve_annotate_config(
            {"annotate": {"table": str(t), "fields": ["lineage"]}}
        )


def test_the_runtime_state_stays_json_serialisable(tmp_path):
    """The lockfile keeps `_`-prefixed runtime keys and dumps them with
    default=str, so a parsed-table object here would serialise as
    '<MetadataTable object at 0x...>' and break the lockfile's byte-stability
    between identical runs."""
    cfg = _cfg(tmp_path, _table(tmp_path, ["PQ1\t2.3.4.4b"]))
    json.dumps(cfg[ANNOTATE_RUNTIME])


# ---------------------------------------------------------------------------
# Binding
# ---------------------------------------------------------------------------

def test_labels_bind_by_accession(tmp_path):
    cfg = _cfg(tmp_path, _table(tmp_path, ["PQ1\t2.3.4.4b", "PQ2\t2.3.2.1a"]))
    seqs = [_seq("a", accession="PQ1"), _seq("b", accession="PQ2")]
    annotate_sequences(seqs, cfg)
    assert [s.group_labels["clade"] for s in seqs] == ["2.3.4.4b", "2.3.2.1a"]


def test_an_unbound_sequence_is_labelled_unknown(tmp_path):
    """Not left blank: a missing tag downstream then means "this run had no
    labels", never "this record was not in the table"."""
    cfg = _cfg(tmp_path, _table(tmp_path, ["PQ1\t2.3.4.4b"]))
    seqs = [_seq("a", accession="PQ1"), _seq("b", accession="ZZ9")]
    annotate_sequences(seqs, cfg)
    assert seqs[1].group_labels == {"clade": "Unknown"}
    assert annotate_stats(cfg)["clade"] == {
        "n": 2, "bound": 1, "conflict": 0, "unknown": 1,
    }


def test_disagreeing_rows_are_a_conflict(tmp_path):
    cfg = _cfg(
        tmp_path,
        _table(tmp_path, ["PQ1\t2.3.4.4b", "A/x/1\t2.2"]),
        segmented=True,
    )
    iso = _isolate("A/x/1", ["PQ1", "PQ2"])
    annotate_sequences([iso], cfg)
    assert iso.group_labels == {"clade": "Conflict"}
    assert annotate_stats(cfg)["clade"]["conflict"] == 1


def test_a_segmented_isolate_binds_through_a_segment_accession(tmp_path):
    """The whole point of reusing custom_mode's binding: a CONCAT isolate has
    no accession of its own, so an HA-keyed table must reach it through its
    segments — and both the CONCAT and the segment records must end up
    labelled, because the tree leaf is the CONCAT while the protein-FASTA
    writer sees the segment as the parent."""
    cfg = _cfg(tmp_path, _table(tmp_path, ["PQ4\t2.3.4.4b"]), segmented=True)
    iso = _isolate("A/x/1", ["PQ1", "PQ4"])
    annotate_sequences([iso], cfg)
    assert iso.group_labels == {"clade": "2.3.4.4b"}
    assert all(
        s.group_labels == {"clade": "2.3.4.4b"} for s in iso.concat_segments
    )


def test_an_isolate_name_keyed_table_also_binds(tmp_path):
    cfg = _cfg(
        tmp_path,
        _table(tmp_path, ["A/x/1\t2.3.4.4b"], header="isolate\tclade"),
        segmented=True,
    )
    iso = _isolate("A/x/1", ["PQ1", "PQ2"])
    annotate_sequences([iso], cfg)
    assert iso.group_labels == {"clade": "2.3.4.4b"}


def test_a_table_that_binds_nothing_warns_but_does_not_raise(tmp_path, capsys):
    """Advisory by design: a mis-keyed stratification table is wrong science
    and must stop the run, but a mis-keyed annotation table costs a metadata
    tag — and the bind rate is only knowable hours into a long run."""
    cfg = _cfg(tmp_path, _table(tmp_path, ["ZZ9\t2.3.4.4b"]))
    seqs = [_seq("a", accession="PQ1")]
    annotate_sequences(seqs, cfg)
    err = capsys.readouterr().err
    assert "WARNING" in err and "key column" in err
    assert seqs[0].group_labels == {"clade": "Unknown"}


def test_set_labels_keeps_an_existing_field(tmp_path):
    """The stratification path and the annotate step can both run; whichever
    labelled a field first owns it, so the two never fight."""
    s = _seq("a", clade="from_mode")
    set_labels(s, {"clade": "from_table", "genotype": "B3.13"})
    assert s.group_labels == {"clade": "from_mode", "genotype": "B3.13"}


def test_two_fields_are_both_attached(tmp_path):
    t = _table(
        tmp_path, ["PQ1\t2.3.4.4b\tB3.13"],
        header="accession\tclade\tgenotype",
    )
    cfg = _cfg(tmp_path, t, fields=("clade", "genotype"))
    seqs = [_seq("a", accession="PQ1")]
    annotate_sequences(seqs, cfg)
    assert seqs[0].group_labels == {"clade": "2.3.4.4b", "genotype": "B3.13"}


# ---------------------------------------------------------------------------
# Surface 1: protein FASTA bracket tags
# ---------------------------------------------------------------------------

def _header(seq, prot):
    import io as _io

    fh = _io.StringIO()
    _write_protein_fasta_record(fh, prot, seq, seq.isolate_id)
    return fh.getvalue().splitlines()[0]


def test_the_protein_fasta_header_carries_the_label():
    seq = _seq("a", accession="PQ1", clade="2.3.4.4b")
    header = _header(seq, _prot("PQ1_p"))
    assert "[clade=2.3.4.4b]" in header
    # Beside the biological context it refines, not appended at the end:
    # ... [host=] [subtype=] [clade=] [country=] ... [length=] [parent=]
    assert (
        header.index("[host=")
        < header.index("[subtype=")
        < header.index("[clade=")
        < header.index("[length=")
    )


def test_a_label_named_after_real_metadata_is_not_duplicated():
    """`--annotate-field host` must not put a second, possibly different
    [host=...] in the header."""
    seq = _seq("a", accession="PQ1")
    seq.group_labels = {"host": "something else", "clade": "2.3.4.4b"}
    header = _header(seq, _prot("PQ1_p"))
    assert header.count("[host=") == 1
    assert "[host=Avian]" in header
    assert "[clade=2.3.4.4b]" in header


def test_an_unlabelled_run_leaves_the_header_exactly_as_before():
    seq = _seq("a", accession="PQ1")
    assert "[clade=" not in _header(seq, _prot("PQ1_p"))


# ---------------------------------------------------------------------------
# Surface 2: phyloXML leaf properties (labels/figures deliberately unchanged)
# ---------------------------------------------------------------------------

def _write_xml(tmp_path, name, newick, reps):
    nwk = tmp_path / f"{name}.nwk"
    nwk.write_text(newick)
    id_map = {f"S{i + 1:04d}": r.id for i, r in enumerate(reps)}
    out = tmp_path / f"{name}_tree.xml"
    write_phyloxml(
        nwk, out, reps, id_map, cfg={}, prefix="test",
        alphabet="nucleotide", msa_tool="MAFFT", msa_version="v7",
        tree_tool="FastTree", tree_version="2.1", model="GTR", ufboot=None,
    )
    return out


def test_every_leaf_gets_a_repseq_property(tmp_path):
    reps = [_seq("A", clade="2.3.4.4b"), _seq("B", clade="2.2"), _seq("C")]
    xml = _write_xml(
        tmp_path, "test", "((S0001:1,S0002:1):1,S0003:1);", reps
    ).read_text()
    assert xml.count('ref="repseq:clade"') == 2
    assert "2.3.4.4b" in xml and ">2.2<" in xml


def test_the_leaf_labels_are_untouched(tmp_path):
    """Designer's call: the clade reaches the trees as a property only, so no
    tree figure, MSA header or leaf label changes shape."""
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir()
    b.mkdir()
    plain = _write_xml(a, "t", "(S0001:1,S0002:1);",
                       [_seq("A"), _seq("B")]).read_text()
    labelled = _write_xml(b, "t", "(S0001:1,S0002:1);",
                          [_seq("A", clade="2.2"), _seq("B", clade="2.2")]).read_text()
    names = [ln for ln in plain.splitlines() if "<name>" in ln]
    assert names == [ln for ln in labelled.splitlines() if "<name>" in ln]


def test_a_label_colliding_with_a_metadata_property_is_skipped(tmp_path):
    rep = _seq("A")
    rep.group_labels = {"host": "something else"}
    xml = _write_xml(tmp_path, "test", "(S0001:1,S0002:1);",
                     [rep, _seq("B")]).read_text()
    assert xml.count('ref="repseq:host"') == 2      # the real metadata only
    assert "something else" not in xml


# ---------------------------------------------------------------------------
# Surface 3: monophyly / segment-status on a NON-rank grouping
# ---------------------------------------------------------------------------

def test_extra_keys_are_assessed_like_a_rank(tmp_path):
    """The flat-taxonomy payoff: these leaves are all one species and one
    genus, so every rank row is trivially monophyletic — the clade rows are
    the only informative ones."""
    reps = [
        _seq("A", clade="2.3.4.4b"), _seq("B", clade="2.3.4.4b"),
        _seq("C", clade="2.2"), _seq("D", clade="2.2"),
    ]
    # Clade 2.2 is split by a 2.3.4.4b leaf → not monophyletic.
    _write_xml(tmp_path, "test",
               "(((S0001:1,S0003:1):1,S0002:1):1,S0004:1);", reps)
    out = write_monophyly_report(tmp_path, "test", extra_keys=("clade",))
    rows = [
        dict(zip(out.read_text().splitlines()[0].split("\t"), ln.split("\t")))
        for ln in out.read_text().splitlines()[1:]
    ]
    clade_rows = {r["taxon"]: r for r in rows if r["rank"] == "clade"}
    assert set(clade_rows) == {"2.3.4.4b", "2.2"}
    assert clade_rows["2.2"]["status"] in ("paraphyletic", "polyphyletic")
    # Without the config key, nothing is assessed on it.
    out2 = write_monophyly_report(tmp_path, "nokeys")
    assert not any(
        "clade" in ln.split("\t")[1] for ln in out2.read_text().splitlines()[1:]
    )


def test_extra_key_rows_sort_above_the_rank_rows():
    assert group_sort_index("clade") < group_sort_index("species")
    assert group_sort_index("species") < group_sort_index("genus")


def test_the_segment_status_matrix_picks_up_the_extra_key(tmp_path):
    reps = [
        _seq("A", clade="2.3.4.4b"), _seq("B", clade="2.3.4.4b"),
        _seq("C", clade="2.2"), _seq("D", clade="2.2"),
    ]
    _write_xml(tmp_path, "test", "((S0001:1,S0002:1):1,(S0003:1,S0004:1):1);",
               reps)
    _write_xml(tmp_path, "test_seg4",
               "(((S0001:1,S0003:1):1,S0002:1):1,S0004:1);", reps)
    write_monophyly_report(tmp_path, "test", extra_keys=("clade",))
    out = write_segment_status_matrix(tmp_path, "test")
    assert out is not None
    rows = [ln.split("\t") for ln in out.read_text().splitlines()[1:]]
    assert any(r[0] == "clade" and r[1] == "2.2" for r in rows)


# ---------------------------------------------------------------------------
# Replay must not lose the labels
# ---------------------------------------------------------------------------

def test_the_lockfile_records_the_labels():
    """They came from a user table, so they are not recoverable from the
    record — a replay without them would quietly drop the tags."""
    entry = _representative_entry(_seq("a", accession="PQ1", clade="2.3.4.4b"))
    assert entry["group_labels"] == {"clade": "2.3.4.4b"}
    # Absent (not null) for an unlabelled run, so an old lockfile replays
    # unchanged.
    assert "group_labels" not in _representative_entry(_seq("b"))


def test_a_segmented_lockfile_entry_records_them_too():
    iso = _isolate("A/x/1", ["PQ1", "PQ2"])
    set_labels(iso, {"clade": "2.2"})
    assert _representative_entry(iso)["group_labels"] == {"clade": "2.2"}


# ---------------------------------------------------------------------------
# The Methods prose must not drift
# ---------------------------------------------------------------------------

def test_the_summary_reports_the_annotation(tmp_path):
    cfg = _cfg(tmp_path, _table(tmp_path, ["PQ1\t2.3.4.4b"]))
    annotate_sequences([_seq("a", accession="PQ1"), _seq("b")], cfg)
    text = _render_selection(
        cfg, RunResult(mode="global", representatives=[], clusters=[]),
        QCReport(total_input=2),
    )
    assert "External labels for `clade`" in text
    assert "**1 of 2**" in text
    assert "labelled `Unknown`" in text
    assert "did not affect which sequences were selected" in text


def test_the_summary_names_the_non_rank_groupings(tmp_path):
    from repseq.output.summary import _extra_keys_note

    assert _extra_keys_note({}) == ""
    note = _extra_keys_note(
        {"phylo": {"monophyly": {"extra_keys": ["clade", "subtype"]}}}
    )
    assert "`clade`, `subtype`" in note
    assert "taxonomy is flat" in note


def test_the_shipped_reference_config_documents_both_keys():
    import importlib.resources as resources

    import yaml

    from repseq.config import validate_config

    text = (resources.files("repseq.data") / "default_config.yaml").read_text()
    cfg = yaml.safe_load(text)
    assert cfg["annotate"] == {"table": None, "fields": []}
    assert cfg["phylo"]["monophyly"]["extra_keys"] == []
    assert validate_config(cfg) == []
    # The comments have to explain the flat-taxonomy payoff, or nobody finds
    # the key that makes the monophyly report usable on influenza.
    assert "extra_keys: [subtype]" in text
    assert "annotate.fields" in text


# ---------------------------------------------------------------------------
# Follow-ups from the 2026-10-29 review of this feature
# ---------------------------------------------------------------------------

def test_extra_keys_are_sanitised_like_the_property_they_must_match(tmp_path):
    """The spelling the `--annotate-field` help tells the user to reuse has
    to work.

    `--annotate-field 'HA clade'` emits `repseq:ha_clade`, so matching the
    raw config string found nothing — and produced the silently empty report
    that this key exists to prevent.
    """
    reps = [_seq(c, clade=None) for c in "ABCD"]
    for r, v in zip(reps, ["x", "x", "y", "y"]):
        r.group_labels = {"HA clade": v}
    _write_xml(tmp_path, "t", "(((S0001:1,S0003:1):1,S0002:1):1,S0004:1);", reps)
    for keys in (("HA clade",), ("ha_clade",), ("  HA.clade ",)):
        out = write_monophyly_report(
            tmp_path, f"p{abs(hash(keys)) % 9999}", extra_keys=keys,
        )
        ranks = [ln.split("\t")[1] for ln in out.read_text().splitlines()[1:]]
        assert ranks.count("ha_clade") == 2, keys


def test_an_extra_key_matching_nothing_is_reported(tmp_path, caplog):
    """Zero rows is indistinguishable from "nothing to report", which is the
    misreading the key exists to prevent — so say it, with a hint."""
    reps = [_seq("A", clade="2.2"), _seq("B", clade="2.2")]
    _write_xml(tmp_path, "t", "(S0001:1,S0002:1);", reps)
    with caplog.at_level("WARNING"):
        write_monophyly_report(tmp_path, "t", extra_keys=("clades",))
    assert "matched" not in caplog.text or "no rows" in caplog.text
    assert "'clades'" in caplog.text
    assert "Did you mean 'clade'?" in caplog.text


def test_an_extra_key_naming_a_real_rank_does_not_duplicate_rows(tmp_path):
    reps = [_seq("A"), _seq("B"), _seq("C")]
    _write_xml(tmp_path, "t", "((S0001:1,S0002:1):1,S0003:1);", reps)
    out = write_monophyly_report(tmp_path, "t", extra_keys=("genus", "GENUS"))
    rows = [ln.split("\t") for ln in out.read_text().splitlines()[1:]]
    genus = [r for r in rows if r[1] == "genus"]
    assert len(genus) == len({(r[1], r[2]) for r in genus})


def test_two_fields_that_sanitise_alike_are_refused(tmp_path):
    """Both were accepted, only one could ever be emitted — and the Methods
    prose claimed both had been attached."""
    t = _table(
        tmp_path, ["PQ1\tv1\tv2"], header="accession\tHA clade\tHA-clade",
    )
    with pytest.raises(InputError, match="same label key"):
        resolve_annotate_config(
            {"annotate": {"table": str(t), "fields": ["HA clade", "HA-clade"]}}
        )


def test_a_field_named_after_reserved_metadata_is_refused(tmp_path):
    """Silently dropping it everywhere is worse than refusing it: the user
    asked for a column to be emitted and it never would be."""
    t = _table(tmp_path, ["PQ1\tx"], header="accession\tspecies")
    with pytest.raises(InputError, match="collides with"):
        resolve_annotate_config(
            {"annotate": {"table": str(t), "fields": ["species"]}}
        )


def test_a_label_cannot_shadow_a_rank_tag_without_taxonomy():
    """The FASTA skip list used to come from the tags actually built, and the
    9 rank tags are only built when taxonomy resolved — so the same run put a
    table value in `[species=...]` for unresolved records only."""
    seq = _seq("a", accession="PQ1")
    seq.taxonomy = None
    seq.group_labels = {"species": "MY-TABLE-VALUE", "clade": "2.2"}
    header = _header(seq, _prot("PQ1_p"))
    assert "MY-TABLE-VALUE" not in header
    assert "[clade=2.2]" in header


def test_labels_named_after_leaf_properties_reach_neither_surface(tmp_path):
    """`strain` / `year` / `isolate_id` are leaf properties but were not
    FASTA tag keys, so they appeared in headers and never in trees — and
    `[isolate_id=...]` sat beside a differing `[isolate=...]`."""
    seq = _seq("a", accession="PQ1", isolate_id="realiso", strain="realstrain")
    seq.group_labels = {
        "strain": "LBL", "year": "1999", "isolate_id": "LBL-ISO",
    }
    header = _header(seq, _prot("PQ1_p"))
    assert "LBL" not in header and "1999" not in header
    rep = _seq("A")
    rep.group_labels = {"strain": "LBL"}
    xml = _write_xml(tmp_path, "t", "(S0001:1,S0002:1);",
                     [rep, _seq("B")]).read_text()
    assert "LBL" not in xml


def test_the_sliced_peptide_header_carries_the_label():
    """The peptide writer is a near-copy of the protein writer; the label
    block went into one only, so a polyprotein's whole-ORF FASTA had
    [clade=...] while every per-peptide file beside it did not."""
    from repseq.polyprotein.slicer import SlicedPeptide

    sliced = SlicedPeptide(
        peptide_name="NS5B", parent_protein_id="PQ1_p",
        parent_accession="PQ1", range_aa_from=1, range_aa_to=5,
        length_aa=5, sequence="MAAAA", cut_method_actual="boundary",
        status="ok",
    )
    seq = _seq("a", accession="PQ1", clade="2.3.4.4b")
    fh = __import__("io").StringIO()
    _write_sliced_peptide_record(fh, sliced, "HCV_pp", seq, seq.isolate_id)
    header = fh.getvalue().splitlines()[0]
    assert "[clade=2.3.4.4b]" in header
    assert header.index("[subtype=") < header.index("[clade=")
    assert header.index("[clade=") < header.index("[peptide=")


def test_a_rewritten_table_is_not_served_from_cache(tmp_path):
    """Keyed on the path alone, a second in-process run with a rewritten
    table silently got the first file's rows."""
    t = _table(tmp_path, ["PQ1\tv1"])
    cfg = _cfg(tmp_path, t)
    seqs = [_seq("a", accession="PQ1")]
    annotate_sequences(seqs, cfg)
    assert seqs[0].group_labels == {"clade": "v1"}

    import os
    import time

    time.sleep(0.01)
    t.write_text("accession\tclade\nPQ1\tv2\n")
    os.utime(t, None)
    cfg2 = _cfg(tmp_path, t)
    seqs2 = [_seq("a", accession="PQ1")]
    annotate_sequences(seqs2, cfg2)
    assert seqs2[0].group_labels == {"clade": "v2"}


def test_the_selection_value_outranks_an_annotate_table(tmp_path, capsys):
    """Both can label the same field, and the annotate step runs first — but
    the value that CHOSE the representatives is the one the audit file
    records, so it must be the one the headers and trees show."""
    s = _seq("a", accession="PQ1")
    set_labels(s, {"clade": "from_table"})
    replaced = set_labels(s, {"clade": "from_selection"}, overwrite=True)
    assert s.group_labels == {"clade": "from_selection"}
    assert replaced == ["clade"]
    # Same value from both sources is not a replacement worth reporting.
    assert set_labels(s, {"clade": "from_selection"}, overwrite=True) == []


def test_a_null_segmented_block_does_not_crash(tmp_path):
    """`segmented:` alone is valid YAML for None, and this ran after QC and
    assembly — a traceback there would throw away hours of work."""
    cfg = _cfg(tmp_path, _table(tmp_path, ["PQ1\t2.2"]))
    cfg["segmented"] = None
    annotate_sequences([_seq("a", accession="PQ1")], cfg)


def test_several_columns_resolve_the_matched_rows_once(tmp_path):
    """A segmented CONCAT offers ~34 identifiers, each normalised, so a
    field-major lookup paid for all of that once per column."""
    from repseq.modes import custom_mode

    t = _table(
        tmp_path, ["PQ4\t2.3.4.4b\tB3.13"],
        header="accession\tclade\tgenotype",
    )
    cfg = _cfg(tmp_path, t, fields=("clade", "genotype"), segmented=True)
    iso = _isolate("A/x/1", ["PQ1", "PQ4"])

    calls = {"n": 0}
    real = custom_mode.matched_rows

    def counting(seq, table):
        calls["n"] += 1
        return real(seq, table)

    custom_mode.matched_rows = counting
    try:
        annotate_sequences([iso], cfg)
    finally:
        custom_mode.matched_rows = real
    assert calls["n"] == 1                      # not once per field
    assert iso.group_labels == {"clade": "2.3.4.4b", "genotype": "B3.13"}
