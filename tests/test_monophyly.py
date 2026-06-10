"""Per-taxon monophyly report (2J).

The end-to-end tests drive the REAL phyloXML writer, so the parser is
validated against the exact ``<property>`` format the pipeline emits rather
than a hand-rolled fixture that could drift, then run
``write_monophyly_report`` over the output directory.
"""
from __future__ import annotations

from repseq.models import Sequence, SequenceSource, SequenceType, TaxonomyInfo
from repseq.phylo.monophyly import write_monophyly_report
from repseq.phylo.phyloxml_writer import write_phyloxml


def _seq(sid: str, genus: str) -> Sequence:
    return Sequence(
        id=sid, header=sid, sequence="ACGT" * 10,
        seq_type=SequenceType.NUCLEOTIDE, source=SequenceSource.NCBI,
        accession=sid, organism=genus,
        taxonomy=TaxonomyInfo(genus=genus),
    )


def _seq_sp(sid: str, genus: str, species: str) -> Sequence:
    return Sequence(
        id=sid, header=sid, sequence="ACGT" * 10,
        seq_type=SequenceType.NUCLEOTIDE, source=SequenceSource.NCBI,
        accession=sid, organism=species,
        taxonomy=TaxonomyInfo(genus=genus, species=species),
    )


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


def _rows(tsv_path):
    lines = tsv_path.read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"))) for ln in lines[1:]]


def test_monophyletic_genus(tmp_path):
    reps = [
        _seq("A", "Alphacoronavirus"),
        _seq("B", "Alphacoronavirus"),
        _seq("C", "Alphacoronavirus"),
    ]
    _write_xml(tmp_path, "test", "((S0001:1,S0002:1):1,S0003:1);", reps)
    out = write_monophyly_report(tmp_path, "test")
    assert out is not None
    rows = [r for r in _rows(out) if r["rank"] == "genus"]
    assert len(rows) == 1
    assert rows[0]["taxon"] == "Alphacoronavirus"
    assert rows[0]["status"] == "monophyletic"
    assert rows[0]["n_clusters"] == "1"
    assert rows[0]["n_intruders"] == "0"


def test_polyphyletic_intruders_in_two_clades(tmp_path):
    # (((A,X1),(B,X2)),Y): Alphacoronavirus = {A,B} has Betacoronavirus
    # leaves X1, X2 nested between its members as TWO separate intruder
    # blocks → polyphyletic.
    reps = [
        _seq("A", "Alphacoronavirus"),
        _seq("X1", "Betacoronavirus"),
        _seq("B", "Alphacoronavirus"),
        _seq("X2", "Betacoronavirus"),
        _seq("Y", "Betacoronavirus"),
    ]
    _write_xml(
        tmp_path, "test",
        "(((S0001:1,S0002:1):1,(S0003:1,S0004:1):1):1,S0005:1);", reps,
    )
    rows = {
        r["taxon"]: r
        for r in _rows(write_monophyly_report(tmp_path, "test"))
        if r["rank"] == "genus"
    }
    assert rows["Alphacoronavirus"]["status"] == "polyphyletic"
    assert rows["Alphacoronavirus"]["n_intruders"] == "2"
    assert rows["Alphacoronavirus"]["intruder_clusters"] == "2"
    assert rows["Alphacoronavirus"]["intruder_taxa"] == "Betacoronavirus"
    assert rows["Alphacoronavirus"]["n_clusters"] == "2"


def test_paraphyletic_single_excluded_clade(tmp_path):
    # (((T1,T2),(F1,F2)),T3): GenusT = {T1,T2,T3} is an ancestral grade with a
    # single derived GenusF clade (F1,F2) nested inside → paraphyletic; the
    # carved-out GenusF is itself monophyletic.
    reps = [
        _seq("T1", "Alphacoronavirus"),
        _seq("T2", "Alphacoronavirus"),
        _seq("F1", "Betacoronavirus"),
        _seq("F2", "Betacoronavirus"),
        _seq("T3", "Alphacoronavirus"),
    ]
    _write_xml(
        tmp_path, "test",
        "(((S0001:1,S0002:1):1,(S0003:1,S0004:1):1):1,S0005:1);", reps,
    )
    rows = {
        r["taxon"]: r
        for r in _rows(write_monophyly_report(tmp_path, "test"))
        if r["rank"] == "genus"
    }
    assert rows["Alphacoronavirus"]["status"] == "paraphyletic"
    assert rows["Alphacoronavirus"]["intruder_clusters"] == "1"
    assert rows["Alphacoronavirus"]["intruder_taxa"] == "Betacoronavirus"
    assert rows["Betacoronavirus"]["status"] == "monophyletic"


def test_singleton_taxon_not_reported(tmp_path):
    # Only one Alphacoronavirus leaf → no meaningful monophyly status.
    reps = [
        _seq("A", "Alphacoronavirus"),
        _seq("B", "Betacoronavirus"),
        _seq("C", "Betacoronavirus"),
    ]
    _write_xml(tmp_path, "test", "((S0001:1,S0002:1):1,S0003:1);", reps)
    out = write_monophyly_report(tmp_path, "test")
    taxa = {r["taxon"] for r in _rows(out) if r["rank"] == "genus"}
    assert "Alphacoronavirus" not in taxa     # singleton skipped
    assert "Betacoronavirus" in taxa


def test_multiple_trees_each_assessed(tmp_path):
    reps = [_seq(x, "Alphacoronavirus") for x in ("A", "B", "C")]
    _write_xml(tmp_path, "test", "((S0001:1,S0002:1):1,S0003:1);", reps)
    sub = tmp_path / "test_per_protein"
    sub.mkdir()
    _write_xml(sub, "test_Spike", "((S0001:1,S0002:1):1,S0003:1);", reps)
    out = write_monophyly_report(tmp_path, "test")
    trees = {r["tree"] for r in _rows(out)}
    assert len(trees) == 2   # genome xml + per-protein xml both swept


def test_no_trees_returns_none(tmp_path):
    assert write_monophyly_report(tmp_path, "test") is None


def _genus_rows(tmp_path, min_support):
    return {
        r["taxon"]: r
        for r in _rows(write_monophyly_report(tmp_path, "test",
                                              min_support=min_support))
        if r["rank"] == "genus"
    }


def test_support_aware_collapses_weak_intrusion(tmp_path):
    """(A,(F,B)) with the (F,B) node WEAKLY supported: topology-only calls
    Alphacoronavirus paraphyletic, but at min_support=70 the weak branch
    collapses and it reads as monophyletic (compatibility test)."""
    reps = [
        _seq("A", "Alphacoronavirus"),
        _seq("F", "Betacoronavirus"),
        _seq("B", "Alphacoronavirus"),
    ]
    # FastTree sh_like support 0.30 → 30 in the phyloXML (< 70).
    _write_xml(tmp_path, "test",
               "(S0001:0.1,(S0002:0.1,S0003:0.1)0.30:0.1);", reps)
    assert _genus_rows(tmp_path, 0)["Alphacoronavirus"]["status"] == "paraphyletic"
    aware = _genus_rows(tmp_path, 70)["Alphacoronavirus"]
    assert aware["status"] == "monophyletic"
    assert aware["min_support"] == "70"
    assert aware["n_intruders"] == "0"


def test_support_aware_keeps_strong_intrusion(tmp_path):
    """Same topology but the (F,B) node is STRONGLY supported (0.95 → 95):
    even at min_support=70 it survives the collapse, so the intrusion stands
    and Alphacoronavirus stays non-monophyletic."""
    reps = [
        _seq("A", "Alphacoronavirus"),
        _seq("F", "Betacoronavirus"),
        _seq("B", "Alphacoronavirus"),
    ]
    _write_xml(tmp_path, "test",
               "(S0001:0.1,(S0002:0.1,S0003:0.1)0.95:0.1);", reps)
    assert _genus_rows(tmp_path, 70)["Alphacoronavirus"]["status"] == "paraphyletic"


# --- species rank (opt-in, phylo.monophyly.include_species) ---------------

def test_species_rank_off_by_default(tmp_path):
    """Species rows are emitted only when include_species is set; the default
    sweep stays at subgenus→class (annotation-noise guard)."""
    reps = [
        _seq_sp("A", "Orthobunyavirus", "Orthobunyavirus bunyamweraense"),
        _seq_sp("B", "Orthobunyavirus", "Orthobunyavirus bunyamweraense"),
        _seq_sp("C", "Orthobunyavirus", "Orthobunyavirus cacheense"),
        _seq_sp("D", "Orthobunyavirus", "Orthobunyavirus cacheense"),
    ]
    _write_xml(tmp_path, "test",
               "((S0001:1,S0002:1):1,(S0003:1,S0004:1):1);", reps)
    rows = _rows(write_monophyly_report(tmp_path, "test"))
    assert not any(r["rank"] == "species" for r in rows)
    # genus is still assessed and monophyletic
    assert any(r["rank"] == "genus" and r["status"] == "monophyletic"
               for r in rows)


def test_species_rank_included_when_opted_in(tmp_path):
    reps = [
        _seq_sp("A", "Orthobunyavirus", "Orthobunyavirus bunyamweraense"),
        _seq_sp("B", "Orthobunyavirus", "Orthobunyavirus bunyamweraense"),
        _seq_sp("C", "Orthobunyavirus", "Orthobunyavirus cacheense"),
        _seq_sp("D", "Orthobunyavirus", "Orthobunyavirus cacheense"),
    ]
    _write_xml(tmp_path, "test",
               "((S0001:1,S0002:1):1,(S0003:1,S0004:1):1);", reps)
    sp = {
        r["taxon"]: r
        for r in _rows(write_monophyly_report(
            tmp_path, "test", include_species=True))
        if r["rank"] == "species"
    }
    assert sp["Orthobunyavirus bunyamweraense"]["status"] == "monophyletic"
    assert sp["Orthobunyavirus cacheense"]["status"] == "monophyletic"


def test_species_polyphyly_is_the_reassortment_signal(tmp_path):
    """A species split across two clades (its segments grouping with a
    different species) is the species-level reassortment signal that the
    coarser ranks cannot see — this is the whole point of include_species."""
    reps = [
        _seq_sp("A", "G", "Species X"),
        _seq_sp("B", "G", "Species Y"),
        _seq_sp("C", "G", "Species X"),
        _seq_sp("D", "G", "Species Y"),
    ]
    # ((X,Y),(X,Y)) — each species' two leaves land in different clades.
    _write_xml(tmp_path, "test",
               "((S0001:1,S0002:1):1,(S0003:1,S0004:1):1);", reps)
    sp = {
        r["taxon"]: r
        for r in _rows(write_monophyly_report(
            tmp_path, "test", include_species=True))
        if r["rank"] == "species"
    }
    assert sp["Species X"]["status"] == "polyphyletic"
    assert sp["Species X"]["intruder_clusters"] == "2"
    assert sp["Species Y"]["status"] == "polyphyletic"
