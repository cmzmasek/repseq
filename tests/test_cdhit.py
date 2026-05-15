"""cd-hit clustering backend: .clstr parsing, dispatch, threshold floor.

Tests run without the cd-hit binary by mocking ``subprocess.run`` to write
synthetic ``.clstr`` output that the parser then reads back. The aim is to
lock the round-trip (input seq.id ↔ .clstr id) and the auto-binary /
auto-word-size choices that depend on input alphabet.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from repseq.clustering import min_threshold, run_clustering as dispatch_run_clustering
from repseq.clustering.cdhit import (
    CDHitError,
    _is_protein,
    _parse_clstr_file,
    _pick_word_size,
    run_clustering,
)
from repseq.models import Sequence, SequenceType


def _seq(sid: str, seq: str, seq_type: SequenceType = SequenceType.PROTEIN) -> Sequence:
    return Sequence(id=sid, header=sid, sequence=seq, seq_type=seq_type, accession=sid)


# ---------------------------------------------------------------------------
# Word-size auto-pick
# ---------------------------------------------------------------------------

@pytest.mark.parametrize(
    "threshold,protein,expected",
    [
        (0.95, True, 5),
        (0.70, True, 5),
        (0.65, True, 4),
        (0.55, True, 3),
        (0.45, True, 2),
        (0.99, False, 10),
        (0.92, False, 8),
        (0.89, False, 7),
        (0.86, False, 6),
        (0.82, False, 5),
    ],
)
def test_pick_word_size_table(threshold, protein, expected):
    assert _pick_word_size(threshold, protein) == expected


def test_pick_word_size_below_floor_raises():
    with pytest.raises(CDHitError, match="floor"):
        _pick_word_size(0.30, protein=True)
    with pytest.raises(CDHitError, match="floor"):
        _pick_word_size(0.70, protein=False)


# ---------------------------------------------------------------------------
# Alphabet detection & min_threshold
# ---------------------------------------------------------------------------

def test_is_protein_all_nucleotide_is_false():
    seqs = [_seq("a", "ACGT", SequenceType.NUCLEOTIDE),
            _seq("b", "ACGT", SequenceType.NUCLEOTIDE)]
    assert _is_protein(seqs) is False


def test_is_protein_any_protein_is_true():
    seqs = [_seq("a", "ACGT", SequenceType.NUCLEOTIDE),
            _seq("b", "MKLV", SequenceType.PROTEIN)]
    assert _is_protein(seqs) is True


def test_min_threshold_dispatch_returns_zero_for_mmseqs2():
    cfg = {"clustering": {"backend": "mmseqs2"}}
    assert min_threshold(cfg, [_seq("a", "MK", SequenceType.PROTEIN)]) == 0.0


def test_min_threshold_dispatch_returns_floor_for_cdhit_protein():
    cfg = {"clustering": {"backend": "cdhit"}}
    assert min_threshold(cfg, [_seq("a", "MK", SequenceType.PROTEIN)]) == 0.40


def test_min_threshold_dispatch_returns_floor_for_cdhit_nucleotide():
    cfg = {"clustering": {"backend": "cdhit"}}
    assert min_threshold(cfg, [_seq("a", "ACGT", SequenceType.NUCLEOTIDE)]) == 0.80


# ---------------------------------------------------------------------------
# .clstr parser
# ---------------------------------------------------------------------------

def test_parse_clstr_file_returns_clusters_with_rep_marked_by_star(tmp_path):
    seqs = [
        _seq("seq_A", "MKLPQE"),
        _seq("seq_B", "MKLPQE"),
        _seq("seq_C", "MMMMMM"),
    ]
    clstr = tmp_path / "result.clstr"
    clstr.write_text(
        ">Cluster 0\n"
        "0\t6aa, >seq_A... *\n"
        "1\t6aa, >seq_B... at 99.93%\n"
        ">Cluster 1\n"
        "0\t6aa, >seq_C... *\n"
    )
    clusters = _parse_clstr_file(str(clstr), seqs)
    assert len(clusters) == 2
    by_rep = {c.representative.id: c for c in clusters}
    assert set(by_rep) == {"seq_A", "seq_C"}
    assert [m.id for m in by_rep["seq_A"].members] == ["seq_B"]
    assert by_rep["seq_C"].members == []


def test_parse_clstr_file_handles_pipe_and_slash_ids(tmp_path):
    # Regression: cd-hit prints '>id...' even with weirdly-punctuated ids;
    # the parser must reliably strip the trailing '...' but keep any
    # internal punctuation (CONCAT|iso pattern, isolate names with slashes).
    seqs = [
        _seq("CONCAT|iso_1", "ACGT"),
        _seq("strain/2009/H1N1", "ACGT"),
    ]
    clstr = tmp_path / "result.clstr"
    clstr.write_text(
        ">Cluster 0\n"
        "0\t4nt, >CONCAT|iso_1... *\n"
        ">Cluster 1\n"
        "0\t4nt, >strain/2009/H1N1... *\n"
    )
    clusters = _parse_clstr_file(str(clstr), seqs)
    rep_ids = sorted(c.representative.id for c in clusters)
    assert rep_ids == ["CONCAT|iso_1", "strain/2009/H1N1"]


# ---------------------------------------------------------------------------
# run_clustering: end-to-end with subprocess mocked
# ---------------------------------------------------------------------------

def _fake_subprocess_writer(clstr_body: str):
    """Build a fake subprocess.run that writes a .clstr file at <output>.clstr."""

    def _run(cmd, **kwargs):
        # cd-hit command: [bin, "-i", in, "-o", out_prefix, "-c", thresh, ...]
        out_prefix = cmd[cmd.index("-o") + 1]
        Path(out_prefix + ".clstr").write_text(clstr_body)
        # Also write the empty output FASTA cd-hit would normally produce;
        # _parse_clstr_file doesn't read it but be tidy.
        Path(out_prefix).write_text("")

        class _R:
            stderr = ""
            stdout = ""

        return _R()

    return _run


def test_run_clustering_round_trip_protein(tmp_path):
    seqs = [
        _seq("P1", "MKLPQEFIL"),
        _seq("P2", "MKLPQEFIL"),
    ]
    clstr_body = (
        ">Cluster 0\n"
        "0\t9aa, >P1... *\n"
        "1\t9aa, >P2... at 100.00%\n"
    )
    with patch("repseq.clustering.cdhit._check_binary", return_value="cd-hit"), \
         patch("repseq.clustering.cdhit.subprocess.run",
               side_effect=_fake_subprocess_writer(clstr_body)):
        clusters = run_clustering(seqs, 0.95, {"temp_dir": str(tmp_path)})

    assert len(clusters) == 1
    assert clusters[0].representative.id == "P1"
    assert [m.id for m in clusters[0].members] == ["P2"]


def test_run_clustering_below_floor_raises(tmp_path):
    # cd-hit-est refuses identity < 0.80; the wrapper should catch this
    # before launching the binary, with a clear message.
    seqs = [_seq("N1", "ACGTACGT", SequenceType.NUCLEOTIDE)]
    with pytest.raises(CDHitError, match=r">= 0\.8"):
        run_clustering(seqs, 0.75, {"temp_dir": str(tmp_path)})


def test_run_clustering_raises_when_round_trip_drops_sequences(tmp_path):
    # If the .clstr references ids that don't match any input, the wrapper
    # must raise — silently returning a short cluster list would let the
    # binary-search caller misread it as a successful undershoot.
    seqs = [
        _seq("seq_X", "MKL"),
        _seq("seq_Y", "MKL"),
    ]
    # Truncated/garbled cluster output that won't round-trip.
    clstr_body = (
        ">Cluster 0\n"
        "0\t3aa, >ghost... *\n"
    )
    with patch("repseq.clustering.cdhit._check_binary", return_value="cd-hit"), \
         patch("repseq.clustering.cdhit.subprocess.run",
               side_effect=_fake_subprocess_writer(clstr_body)):
        with pytest.raises(CDHitError, match="round-trip"):
            run_clustering(seqs, 0.9, {"temp_dir": str(tmp_path)})


def test_run_clustering_picks_cdhit_est_for_nucleotide(tmp_path):
    seqs = [
        _seq("N1", "ACGTACGT", SequenceType.NUCLEOTIDE),
        _seq("N2", "ACGTACGT", SequenceType.NUCLEOTIDE),
    ]
    clstr_body = (
        ">Cluster 0\n"
        "0\t8nt, >N1... *\n"
        "1\t8nt, >N2... at 100.00%\n"
    )
    seen_binary: list[str] = []

    def _writer(cmd, **kwargs):
        seen_binary.append(Path(cmd[0]).name)
        # delegate to the standard fake to write the .clstr
        return _fake_subprocess_writer(clstr_body)(cmd, **kwargs)

    # _check_binary normally checks PATH; stub it so it just echoes back
    # the binary name unchanged so the wrapper's "which binary" decision
    # is visible at cmd[0].
    with patch(
        "repseq.clustering.cdhit._check_binary",
        side_effect=lambda name: name,
    ), patch("repseq.clustering.cdhit.subprocess.run", side_effect=_writer):
        run_clustering(seqs, 0.95, {"temp_dir": str(tmp_path)})

    assert seen_binary == ["cd-hit-est"]


# ---------------------------------------------------------------------------
# Dispatcher
# ---------------------------------------------------------------------------

def test_dispatch_routes_cdhit_backend(tmp_path):
    seqs = [_seq("P1", "MKLPQE")]
    clstr_body = ">Cluster 0\n0\t6aa, >P1... *\n"
    cfg = {"clustering": {"backend": "cdhit"}, "temp_dir": str(tmp_path)}
    with patch("repseq.clustering.cdhit._check_binary", return_value="cd-hit"), \
         patch("repseq.clustering.cdhit.subprocess.run",
               side_effect=_fake_subprocess_writer(clstr_body)):
        clusters = dispatch_run_clustering(seqs, 0.95, cfg)
    assert len(clusters) == 1 and clusters[0].representative.id == "P1"


def test_dispatch_unknown_backend_raises(tmp_path):
    with pytest.raises(ValueError, match="Unknown clustering backend"):
        dispatch_run_clustering([], 0.95, {"clustering": {"backend": "psi-cd-hit"}})
