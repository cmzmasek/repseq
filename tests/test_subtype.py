"""Tests for the viral subtype/serotype field (e.g. influenza H5N1).

Covers the output surfaces that carry it: representative TSVs, protein-FASTA
bracket tags, phyloXML leaf properties, the tree-leaf-label token, and the
segmented CONCAT inheritance. Sourcing from /serotype is covered in
test_protein_qc.py (efetch path) and test_resolver.py (esummary path).
"""

from __future__ import annotations

from pathlib import Path

from repseq.models import Cluster, RunResult, Sequence, SequenceSource, SequenceType
from repseq.output.report import (
    write_proteins_fasta,
    write_representative_isolates_tsv,
    write_representative_sequences_tsv,
)
from repseq.phylo.labels import format_leaf_label
from repseq.phylo.phyloxml_writer import _LEAF_PROPERTIES, _leaf_property_value
from repseq.segmented.completeness import concatenate_isolate


def _seq(sid, **kw):
    return Sequence(
        id=sid, sequence="ACGT" * 10, seq_type=SequenceType.NUCLEOTIDE,
        source=SequenceSource.NCBI, accession=kw.pop("accession", sid),
        description="", header=sid, **kw,
    )


# ---------------------------------------------------------------------------
# Representative TSVs
# ---------------------------------------------------------------------------

def test_subtype_column_in_representative_sequences_tsv(tmp_path: Path):
    rep = _seq("A1", subtype="H5N1", host="duck")
    path = tmp_path / "reps.tsv"
    write_representative_sequences_tsv([rep], path)
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert "subtype" in header
    # Ordered right after host.
    assert header[header.index("host") + 1] == "subtype"
    assert row["subtype"] == "H5N1"


def test_subtype_column_in_representative_isolates_tsv(tmp_path: Path):
    seg = _seq("s4", segment="4", subtype="H5N1")
    concat = Sequence(
        id="CONCAT|iso1", sequence="ACGT" * 10, seq_type=SequenceType.NUCLEOTIDE,
        source=SequenceSource.NCBI, accession=None, description="",
        header="CONCAT|iso1", isolate_id="iso1", subtype="H5N1",
        concat_segments=[seg],
    )
    path = tmp_path / "iso.tsv"
    write_representative_isolates_tsv([concat], path)
    lines = path.read_text().splitlines()
    header = lines[0].split("\t")
    row = dict(zip(header, lines[1].split("\t")))
    assert header[header.index("host") + 1] == "subtype"
    assert row["subtype"] == "H5N1"


def test_subtype_blank_when_absent(tmp_path: Path):
    rep = _seq("A1")  # no subtype
    path = tmp_path / "reps.tsv"
    write_representative_sequences_tsv([rep], path)
    lines = path.read_text().splitlines()
    row = dict(zip(lines[0].split("\t"), lines[1].split("\t")))
    assert row["subtype"] == ""


# ---------------------------------------------------------------------------
# Protein-FASTA bracket tags
# ---------------------------------------------------------------------------

def _rep_with_protein(subtype):
    s = _seq("s1", segment="4", accession="ACC.1", organism="Influenza A virus",
             subtype=subtype)
    s.proteins = [{
        "protein_id": "P_ha", "product": "hemagglutinin",
        "sequence": "M" * 50, "length": 50,
    }]
    return RunResult(mode="x", representatives=[s],
                     clusters=[Cluster(cluster_id="c1", representative=s)])


def test_subtype_in_protein_fasta_header(tmp_path: Path):
    out = tmp_path / "proteins.fasta"
    write_proteins_fasta(_rep_with_protein("H5N1"), complete_isolates=None, path=out)
    header = out.read_text().splitlines()[0]
    assert "[subtype=H5N1]" in header


def test_subtype_omitted_from_header_when_absent(tmp_path: Path):
    out = tmp_path / "proteins.fasta"
    write_proteins_fasta(_rep_with_protein(None), complete_isolates=None, path=out)
    header = out.read_text().splitlines()[0]
    assert "[subtype=" not in header


# ---------------------------------------------------------------------------
# Tree leaf label + phyloXML property
# ---------------------------------------------------------------------------

def test_subtype_leaf_label_token():
    s = _seq("A1", subtype="H5N1", host="duck")
    assert format_leaf_label(s, "{id}|{subtype}|{host}") == "A1|H5N1|duck"


def test_subtype_leaf_label_token_drops_when_empty():
    s = _seq("A1", host="duck")  # no subtype → separator-drop
    assert format_leaf_label(s, "{id}|{subtype}|{host}") == "A1|duck"


def test_subtype_phyloxml_property():
    assert ("repseq:subtype", "subtype") in _LEAF_PROPERTIES
    s = _seq("A1", subtype="H5N1")
    assert _leaf_property_value(s, "subtype") == "H5N1"
    assert _leaf_property_value(_seq("A2"), "subtype") is None


# ---------------------------------------------------------------------------
# Segmented CONCAT inheritance
# ---------------------------------------------------------------------------

def test_concat_inherits_subtype_first_non_empty():
    segs = [_seq("seg4", segment="4", subtype="H5N1"),
            _seq("seg6", segment="6")]
    concat = concatenate_isolate(segs, "iso1")
    assert concat.subtype == "H5N1"


def test_concat_subtype_none_when_no_segment_has_it():
    segs = [_seq("seg4", segment="4"), _seq("seg6", segment="6")]
    concat = concatenate_isolate(segs, "iso1")
    assert concat.subtype is None
