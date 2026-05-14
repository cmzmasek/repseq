"""MMseqs2 wrapper: FASTA-writing and cluster-TSV parsing.

These cover the ID round-trip without needing the mmseqs binary: the
clustering FASTA must use ``seq.id`` as the sole header token so that the
identifiers MMseqs2 emits in its cluster TSV match back to the input
sequences.
"""
from __future__ import annotations

from repseq.clustering.mmseqs2 import _parse_cluster_tsv, _write_id_fasta
from repseq.models import Sequence, SequenceType
from repseq.segmented.completeness import concatenate_isolate


def _seq(sid, seq, header=None, accession=None):
    return Sequence(
        id=sid,
        header=header if header is not None else sid,
        sequence=seq,
        seq_type=SequenceType.NUCLEOTIDE,
        accession=accession,
    )


def test_write_id_fasta_header_token_is_seq_id(tmp_path):
    # Regression: the header written for clustering must be exactly seq.id.
    # MMseqs2 keys its output on the first whitespace token of the header;
    # writing the full descriptive header broke UniProt and CONCAT inputs.
    uniprot = _seq(
        "P12345",
        "MKLPQEFIL",
        header="sp|P12345|VMAIN_HUMAN some protein OS=Homo sapiens",
        accession="P12345",
    )
    concat = concatenate_isolate(
        [_seq("acc1", "ACGT", accession="acc1"),
         _seq("acc2", "TTTT", accession="acc2")],
        "ISO/2009",
    )
    fasta = tmp_path / "input.fasta"
    _write_id_fasta([uniprot, concat], fasta)

    header_tokens = [
        line[1:].strip().split()[0]
        for line in fasta.read_text().splitlines()
        if line.startswith(">")
    ]
    assert header_tokens == ["P12345", "CONCAT|ISO/2009"]


def test_parse_cluster_tsv_matches_id_headers(tmp_path):
    # cluster.tsv produced from an _write_id_fasta input carries seq.id
    # values; the parser must resolve them back to the Sequence objects.
    seqs = [
        _seq("P12345", "MKLPQEFIL", header="sp|P12345|NAME desc", accession="P12345"),
        _seq("P67890", "MKLPQEFAA", header="sp|P67890|NAME2 desc", accession="P67890"),
        _seq("CONCAT|ISO/2009", "ACGTTTTT", header="CONCAT|ISO/2009|acc1|acc2"),
    ]
    tsv = tmp_path / "result_cluster.tsv"
    tsv.write_text(
        "P12345\tP12345\n"
        "P12345\tP67890\n"
        "CONCAT|ISO/2009\tCONCAT|ISO/2009\n"
    )
    clusters = _parse_cluster_tsv(str(tsv), seqs)

    assert len(clusters) == 2
    by_rep = {c.representative.id: c for c in clusters}
    assert set(by_rep) == {"P12345", "CONCAT|ISO/2009"}
    assert [m.id for m in by_rep["P12345"].members] == ["P67890"]
    assert by_rep["CONCAT|ISO/2009"].members == []
