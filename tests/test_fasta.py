"""FASTA parsing and header detection."""
from __future__ import annotations

from pathlib import Path

from repseq.io.fasta import (
    detect_seq_type,
    parse_header,
    read_fasta,
    write_fasta,
)
from repseq.models import SequenceSource, SequenceType


# ---------------------------------------------------------------------------
# detect_seq_type
# ---------------------------------------------------------------------------

def test_detect_seq_type_protein():
    # contains E/F/L/Q — characters that never appear in DNA
    assert detect_seq_type("MEEPQSDPSVEPPLSQETFSDLW") == SequenceType.PROTEIN


def test_detect_seq_type_nucleotide():
    assert detect_seq_type("ATCGATCGATCGATCGATCGATCG") == SequenceType.NUCLEOTIDE


def test_detect_seq_type_empty():
    assert detect_seq_type("") == SequenceType.UNKNOWN


# ---------------------------------------------------------------------------
# parse_header
# ---------------------------------------------------------------------------

def test_parse_header_uniprot_swissprot():
    h = "sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4"
    fields = parse_header(h)
    assert fields["source"] == SequenceSource.UNIPROT
    assert fields["accession"] == "P04637"
    assert fields["organism"] == "Homo sapiens"
    assert fields["taxid"] == 9606
    assert fields["is_reviewed"] is True


def test_parse_header_uniprot_trembl_not_reviewed():
    h = "tr|A0A123|HYP_FAKE Hypothetical protein OS=Foo bar OX=999"
    fields = parse_header(h)
    assert fields["source"] == SequenceSource.UNIPROT
    assert fields["is_reviewed"] is False


def test_parse_header_ncbi_with_organism():
    h = "NP_000537.3 cellular tumor antigen p53 [Homo sapiens]"
    fields = parse_header(h)
    assert fields["source"] == SequenceSource.NCBI_VIRUS or fields["source"] == SequenceSource.NCBI
    # RefSeq is detected regardless of which sub-parser ran
    assert fields["is_refseq"] is True
    assert fields["accession"] == "NP_000537.3"


def test_parse_header_ncbi_refseq_flag():
    h = "NC_026438.1 Influenza A virus segment 1"
    fields = parse_header(h)
    assert fields["is_refseq"] is True


def test_parse_header_ncbi_non_refseq():
    h = "MW626064.1 Influenza A virus polymerase"
    fields = parse_header(h)
    assert fields.get("is_refseq", False) is False


def test_parse_header_unknown_fallback():
    h = "weirdformat_no_match"
    fields = parse_header(h)
    assert fields["source"] == SequenceSource.UNKNOWN
    assert fields["accession"] == "weirdformat_no_match"


# ---------------------------------------------------------------------------
# read_fasta / write_fasta
# ---------------------------------------------------------------------------

FASTA_SAMPLE = """\
>sp|P12345|TEST_HUMAN Test protein OS=Homo sapiens OX=9606
MEEPQSDPSVEPPLSQETFSDLWKLLPEN
NPSL
>NC_026438.1 Influenza A virus segment 1
ATGAAGACTATCATTGCTTTGAGCTACATTTTCTGTCTGGCTCT
"""


def test_read_fasta_roundtrip(tmp_path: Path):
    p = tmp_path / "in.fasta"
    p.write_text(FASTA_SAMPLE)
    seqs = list(read_fasta(p))
    assert len(seqs) == 2

    # First: UniProt SwissProt
    s1 = seqs[0]
    assert s1.source == SequenceSource.UNIPROT
    assert s1.accession == "P12345"
    assert s1.is_reviewed is True
    assert s1.seq_type == SequenceType.PROTEIN
    # multi-line sequence joined
    assert s1.sequence == "MEEPQSDPSVEPPLSQETFSDLWKLLPENNPSL"

    # Second: RefSeq NCBI
    s2 = seqs[1]
    assert s2.is_refseq is True
    assert s2.accession == "NC_026438.1"
    assert s2.seq_type == SequenceType.NUCLEOTIDE


def test_read_fasta_source_override(tmp_path: Path):
    p = tmp_path / "in.fasta"
    p.write_text(FASTA_SAMPLE)
    seqs = list(read_fasta(p, source_override=SequenceSource.NCBI_VIRUS))
    for s in seqs:
        assert s.source == SequenceSource.NCBI_VIRUS


def test_write_fasta_then_read_back(tmp_path: Path, make_seq):
    s = make_seq("X1", "ACGTACGTACGT" * 10, header="X1 some description")
    out = tmp_path / "out.fa"
    write_fasta([s], out, line_width=20)

    # Lines should wrap to width=20
    body_lines = [ln for ln in out.read_text().splitlines() if not ln.startswith(">")]
    assert all(len(ln) <= 20 for ln in body_lines)

    # Re-reading should produce a single sequence with the original content
    [back] = list(read_fasta(out))
    assert back.sequence == s.sequence
