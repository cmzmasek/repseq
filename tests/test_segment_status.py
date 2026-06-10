"""Cross-tree segment-status matrix (Ext-2).

Drives the REAL phyloXML writer + ``write_monophyly_report`` (so the matrix is
fed the exact on-disk ``_monophyly.tsv`` schema it parses in production), then
pivots it with ``write_segment_status_matrix``. The headline assertion is the
``single_marker_break`` localisation: a taxon clean on the whole-genome tree
but broken on exactly one marker tree.
"""
from __future__ import annotations

from repseq.models import Sequence, SequenceSource, SequenceType, TaxonomyInfo
from repseq.phylo.monophyly import write_monophyly_report
from repseq.phylo.phyloxml_writer import write_phyloxml
from repseq.output.segment_status import write_segment_status_matrix


def _seq_sp(sid: str, genus: str, species: str) -> Sequence:
    return Sequence(
        id=sid, header=sid, sequence="ACGT" * 10,
        seq_type=SequenceType.NUCLEOTIDE, source=SequenceSource.NCBI,
        accession=sid, organism=species,
        taxonomy=TaxonomyInfo(genus=genus, species=species),
    )


def _write_xml(out_dir, name, newick, reps):
    nwk = out_dir / f"{name}.nwk"
    nwk.write_text(newick)
    id_map = {f"S{i + 1:04d}": r.id for i, r in enumerate(reps)}
    out = out_dir / f"{name}_tree.xml"
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


# Species X = {A, B}, Species Y = {C, D}. Genome groups each species; the
# marker tree splits X across both clades (its M segment grouped with Y).
_REPS = [
    _seq_sp("A", "G", "G xspecies"),
    _seq_sp("B", "G", "G xspecies"),
    _seq_sp("C", "G", "G yspecies"),
    _seq_sp("D", "G", "G yspecies"),
]
_GENOME_NWK = "((S0001:1,S0002:1):1,(S0003:1,S0004:1):1);"   # X & Y clean
_MARKER_NWK = "((S0001:1,S0003:1):1,(S0002:1,S0004:1):1);"   # X & Y split


def _build(tmp_path, genome_nwk, marker_nwk):
    _write_xml(tmp_path, "test", genome_nwk, _REPS)            # genome tree
    sub = tmp_path / "test_per_protein"
    sub.mkdir()
    _write_xml(sub, "test_M", marker_nwk, _REPS)               # marker tree
    write_monophyly_report(tmp_path, "test", include_species=True)
    return write_segment_status_matrix(tmp_path, "test")


def test_single_marker_break_localises_reassortment(tmp_path):
    matrix = _build(tmp_path, _GENOME_NWK, _MARKER_NWK)
    assert matrix is not None
    sp = {r["taxon"]: r for r in _rows(matrix) if r["rank"] == "species"}
    x = sp["G xspecies"]
    assert x["genome_status"] == "monophyletic"
    assert x["n_nonmono"] == "1"
    assert x["single_marker_break"] == "M"          # the localised call
    assert x["nonmono_trees"] == "M"


def test_clean_taxon_has_no_break(tmp_path):
    # Marker tree identical to genome → every species monophyletic everywhere.
    matrix = _build(tmp_path, _GENOME_NWK, _GENOME_NWK)
    sp = {r["taxon"]: r for r in _rows(matrix) if r["rank"] == "species"}
    x = sp["G xspecies"]
    assert x["n_nonmono"] == "0"
    assert x["single_marker_break"] == ""


def test_break_on_genome_is_not_a_single_marker_call(tmp_path):
    # Genome tree itself splits X → the reference is dirty, so even with the
    # marker also broken we do NOT emit a single_marker_break (no clean ref).
    matrix = _build(tmp_path, _MARKER_NWK, _MARKER_NWK)
    sp = {r["taxon"]: r for r in _rows(matrix) if r["rank"] == "species"}
    x = sp["G xspecies"]
    assert x["genome_status"] != "monophyletic"
    assert x["single_marker_break"] == ""


def test_candidates_sort_to_top_within_rank(tmp_path):
    matrix = _build(tmp_path, _GENOME_NWK, _MARKER_NWK)
    species_rows = [r for r in _rows(matrix) if r["rank"] == "species"]
    # the flagged (single_marker_break non-empty) rows come first
    flagged = [r for r in species_rows if r["single_marker_break"]]
    assert species_rows[0]["single_marker_break"] != ""
    assert len(flagged) >= 1


def test_no_monophyly_report_returns_none(tmp_path):
    assert write_segment_status_matrix(tmp_path, "test") is None
