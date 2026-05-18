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


def test_is_protein_alphabet_override_wins_over_seq_type():
    """clustering.alphabet_for_clustering=protein on NT-typed concat must pick cd-hit."""
    seqs = [_seq("CONCAT|iso1", "ACGT", SequenceType.NUCLEOTIDE)]
    cfg = {"clustering": {"alphabet_for_clustering": "protein"}}
    assert _is_protein(seqs, cfg) is True


def test_is_protein_alphabet_nucleotide_forces_cdhit_est():
    """clustering.alphabet_for_clustering=nucleotide picks cd-hit-est even with protein input."""
    seqs = [_seq("a", "MKLV", SequenceType.PROTEIN)]
    cfg = {"clustering": {"alphabet_for_clustering": "nucleotide"}}
    assert _is_protein(seqs, cfg) is False


def test_min_threshold_floor_follows_alphabet_override():
    """alphabet_for_clustering=protein with NT-typed sequences uses the 0.40 protein floor."""
    seqs = [_seq("a", "ACGT", SequenceType.NUCLEOTIDE)]
    cfg = {"clustering": {"backend": "cdhit", "alphabet_for_clustering": "protein"}}
    assert min_threshold(cfg, seqs) == 0.40


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
    clusters, unmatched = _parse_clstr_file(str(clstr), seqs)
    assert len(clusters) == 2
    assert unmatched == []
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
    clusters, unmatched = _parse_clstr_file(str(clstr), seqs)
    assert unmatched == []
    rep_ids = sorted(c.representative.id for c in clusters)
    assert rep_ids == ["CONCAT|iso_1", "strain/2009/H1N1"]


def test_parse_clstr_file_handles_ellipsis_inside_id(tmp_path):
    # Regression for the round-trip mismatch we hit on a Hantaviridae run:
    # a seq.id containing '...' would be silently dropped because the
    # old non-greedy regex stopped at the first internal '...'. The
    # tightened, end-anchored regex must keep the full id intact.
    seqs = [
        _seq("CONCAT|weird...iso", "MKLPQE"),
        _seq("ordinary_id", "AAAAAA"),
    ]
    clstr = tmp_path / "result.clstr"
    clstr.write_text(
        ">Cluster 0\n"
        "0\t6aa, >CONCAT|weird...iso... *\n"
        ">Cluster 1\n"
        "0\t6aa, >ordinary_id... *\n"
    )
    clusters, unmatched = _parse_clstr_file(str(clstr), seqs)
    assert unmatched == []
    rep_ids = sorted(c.representative.id for c in clusters)
    assert rep_ids == ["CONCAT|weird...iso", "ordinary_id"]


def test_parse_clstr_file_handles_cdhit_est_strand_prefix(tmp_path):
    """Regression for v0.10.2: cd-hit-est member lines carry strand info
    (``at +/99.93%`` or ``at -/99.93%``) that cd-hit protein output does
    not. The old regex matched only ``at 99.93%`` and silently dropped
    every member of every cluster when alphabet=nucleotide — symptom on
    a hantaviridae run was ``91 sequences in, 34 accounted for`` (the
    34 reps had `*`; the 57 members all parsed as no-match)."""
    seqs = [
        _seq("CONCAT|iso1", "ACGTACGTACGT"),
        _seq("CONCAT|iso2", "ACGTACGTACGT"),
        _seq("CONCAT|iso3", "ACGTACGTACGT"),
        _seq("CONCAT|iso4", "ACGTACGTACGT"),
    ]
    clstr = tmp_path / "result.clstr"
    # cd-hit-est nucleotide .clstr — strand info on every member line.
    clstr.write_text(
        ">Cluster 0\n"
        "0\t12nt, >CONCAT|iso1... *\n"
        "1\t12nt, >CONCAT|iso2... at +/99.93%\n"
        "2\t12nt, >CONCAT|iso3... at -/95.00%\n"
        ">Cluster 1\n"
        "0\t12nt, >CONCAT|iso4... *\n"
    )
    clusters, unmatched = _parse_clstr_file(str(clstr), seqs)
    assert unmatched == []
    by_rep = {c.representative.id: c for c in clusters}
    iso1_members = sorted(m.id for m in by_rep["CONCAT|iso1"].members)
    assert iso1_members == ["CONCAT|iso2", "CONCAT|iso3"]
    assert by_rep["CONCAT|iso4"].members == []


def test_parse_clstr_file_still_handles_cdhit_protein_no_strand(tmp_path):
    """Regression guard: the v0.10.2 regex change must not break the
    plain ``at 99.93%`` format that cd-hit protein emits."""
    seqs = [
        _seq("rep", "MKLPQE"),
        _seq("mem", "MKLPQF"),
    ]
    clstr = tmp_path / "result.clstr"
    clstr.write_text(
        ">Cluster 0\n"
        "0\t6aa, >rep... *\n"
        "1\t6aa, >mem... at 99.93%\n"
    )
    clusters, unmatched = _parse_clstr_file(str(clstr), seqs)
    assert unmatched == []
    assert len(clusters) == 1
    assert clusters[0].representative.id == "rep"
    assert [m.id for m in clusters[0].members] == ["mem"]


def test_parse_clstr_file_reports_unmatched_clstr_ids(tmp_path):
    # If cd-hit ever emits an id we can't match back (length cap, truncation,
    # weird transform), _parse_clstr_file must surface it so the caller can
    # build an informative round-trip error instead of silently undercounting.
    seqs = [_seq("seq_A", "MKLPQE"), _seq("seq_B", "MMMMMM")]
    clstr = tmp_path / "result.clstr"
    clstr.write_text(
        ">Cluster 0\n"
        "0\t6aa, >seq_A... *\n"
        "1\t6aa, >totally_unknown... at 99.93%\n"
        ">Cluster 1\n"
        "0\t6aa, >seq_B... *\n"
    )
    clusters, unmatched = _parse_clstr_file(str(clstr), seqs)
    assert unmatched == ["totally_unknown"]
    # The cluster whose member couldn't match should still be kept,
    # just without that member — so the count gap is visible.
    by_rep = {c.representative.id: c for c in clusters}
    assert by_rep["seq_A"].members == []
    assert by_rep["seq_B"].members == []


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
        with pytest.raises(CDHitError) as exc:
            run_clustering(seqs, 0.9, {"temp_dir": str(tmp_path)})
    msg = str(exc.value)
    assert "round-trip" in msg
    # The error must name BOTH sides of the gap so the user can debug:
    # input seqs that never made it into the .clstr, AND .clstr ids that
    # we couldn't match back to any input.
    assert "seq_X" in msg and "seq_Y" in msg
    assert "ghost" in msg


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
