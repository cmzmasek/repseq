"""MMseqs2 wrapper: FASTA-writing and cluster-TSV parsing.

These cover the ID round-trip without needing the mmseqs binary: the
clustering FASTA must use ``seq.id`` as the sole header token so that the
identifiers MMseqs2 emits in its cluster TSV match back to the input
sequences.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from repseq.clustering.mmseqs2 import (
    MMseqs2Error,
    _parse_cluster_tsv,
    _write_id_fasta,
    run_clustering,
)
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


def test_write_id_fasta_alphabet_protein_uses_protein_sequence(tmp_path):
    seq = _seq("CONCAT|iso1", "ACGTACGT")
    seq.protein_sequence = "MEEPMEEP"
    fasta = tmp_path / "input.fasta"
    _write_id_fasta([seq], fasta, alphabet="protein")
    body = "".join(
        line for line in fasta.read_text().splitlines()
        if not line.startswith(">")
    )
    assert body == "MEEPMEEP"


def test_write_id_fasta_alphabet_protein_raises_when_missing(tmp_path):
    seq = _seq("a", "ACGT")  # no protein_sequence
    fasta = tmp_path / "input.fasta"
    with pytest.raises(ValueError, match="protein_sequence"):
        _write_id_fasta([seq], fasta, alphabet="protein")


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


def test_run_clustering_raises_when_round_trip_drops_sequences(tmp_path):
    # Regression: if the cluster TSV references IDs that don't match
    # any input seq.id (e.g. because the input id contained whitespace
    # and MMseqs2 truncated it), the parser silently drops those rows.
    # The wrapper must catch this and raise rather than returning a
    # short/empty cluster list, because callers like the taxonomic1
    # binary search misread an empty result as a successful undershoot
    # and return every input sequence at threshold = 1.0.
    seqs = [
        _seq("seq with space 1", "ACGTACGT"),
        _seq("seq with space 2", "ACGTACGT"),
    ]

    def fake_run(cmd, **kwargs):
        # run_clustering passes positional args
        # [mmseqs, mode, input_fasta, result_prefix, mmseqs_tmp, ...flags]
        result_prefix = cmd[3]
        tsv = Path(result_prefix + "_cluster.tsv")
        # MMseqs2 would emit the truncated tokens (first whitespace word).
        tsv.write_text(
            "seq\tseq\n"
        )
        class _R:
            stderr = ""
            stdout = ""
        return _R()

    with patch("repseq.clustering.mmseqs2._check_mmseqs2", return_value="mmseqs"), \
         patch("repseq.clustering.mmseqs2.subprocess.run", side_effect=fake_run):
        with pytest.raises(MMseqs2Error, match="round-trip"):
            run_clustering(seqs, 0.9, {"temp_dir": str(tmp_path)})
