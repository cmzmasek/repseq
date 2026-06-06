"""MaxMin diversity selection."""
from __future__ import annotations

from repseq.clustering.diversity import (
    _containment_distance,
    _jaccard_distance,
    _kmer_set,
    basis_sequence,
    kmer_basis,
    select_diverse,
)


def _cfg(alphabet):
    return {"clustering": {"alphabet_for_clustering": alphabet}}


# ---------------------------------------------------------------------------
# Length-robust containment distance (used by select_diverse)
# ---------------------------------------------------------------------------

def test_containment_distance_substring_is_zero():
    # A short sequence whose k-mers are all contained in a much longer one
    # must score ~0 ("not diverse from it"), even though the length gap is
    # large. Jaccard, by contrast, is badly inflated by that gap.
    prefix = "ACGTTGCAACGTTGCAACGTACGT"
    short = _kmer_set(prefix, k=5)
    long = _kmer_set(prefix + "GGGGCCCCAAAATTTTGGGGCCCCAAAATTTT", k=5)
    assert short.issubset(long)
    assert _containment_distance(short, long) == 0.0
    assert _jaccard_distance(short, long) > 0.4  # length gap inflates Jaccard


def test_containment_distance_disjoint_is_one():
    a = _kmer_set("AAAAAAAAAAAA", k=5)
    b = _kmer_set("CCCCCCCCCCCC", k=5)
    assert _containment_distance(a, b) == 1.0


def test_containment_distance_is_symmetric():
    a = _kmer_set("ACGTTGCAACGTTGCA", k=5)
    b = _kmer_set("ACGTTGCAGGGGCCCC", k=5)
    assert _containment_distance(a, b) == _containment_distance(b, a)


def test_select_diverse_returns_all_when_n_geq_pool(make_seq):
    seqs = [make_seq(f"s{i}", "ACGT" * 5) for i in range(3)]
    out = select_diverse(seqs, n=10, seed=42)
    assert len(out) == 3
    assert {s.id for s in out} == {"s0", "s1", "s2"}


def test_select_diverse_starts_with_longest(make_seq):
    a = make_seq("a", "A" * 20)
    b = make_seq("b", "A" * 50)
    c = make_seq("c", "A" * 30)
    out = select_diverse([a, b, c], n=1, seed=42)
    assert out[0].id == "b"


def test_select_diverse_is_deterministic_with_seed(make_seq):
    seqs = [make_seq(f"s{i}", "ACGT" * 5 + "X" * i) for i in range(10)]
    a = select_diverse(seqs, n=4, seed=123)
    b = select_diverse(seqs, n=4, seed=123)
    assert [s.id for s in a] == [s.id for s in b]


def test_select_diverse_prefers_dissimilar_sequences(make_seq):
    """Given two near-identical sequences and one very different, MaxMin
    should pick one from each cluster rather than both near-identicals."""
    sim1 = make_seq("sim1", "AAAAAAAAAACCCCCCCCCC")
    sim2 = make_seq("sim2", "AAAAAAAAAACCCCCCCCCG")  # 1 char off
    diff = make_seq("diff", "TTTTTTTTTTGGGGGGGGGG")
    out = select_diverse([sim1, sim2, diff], n=2, seed=42)
    ids = {s.id for s in out}
    assert "diff" in ids
    # Must NOT include both near-identical sequences
    assert ids != {"sim1", "sim2"}


def test_select_diverse_empty_input():
    assert select_diverse([], n=5, seed=42) == []


# ---------------------------------------------------------------------------
# Alphabet-aware k-mer basis (rigour audit 2026-06-06): k=5 saturates on
# whole-genome NT, so protein-alphabet runs select on the marker protein at
# k=5 and NT runs use a larger k.
# ---------------------------------------------------------------------------

def test_kmer_basis_protein_when_alphabet_protein_and_markers_present(make_seq):
    seqs = [make_seq("a", "ACGT" * 10), make_seq("b", "ACGT" * 10)]
    for s in seqs:
        s.protein_sequence = "MKVL" * 5
    assert kmer_basis(seqs, _cfg("protein")) == (True, 5)


def test_kmer_basis_nucleotide_alphabet_uses_larger_k(make_seq):
    s = make_seq("a", "ACGT" * 10)
    s.protein_sequence = "MKVL"  # present, but the alphabet is nucleotide
    assert kmer_basis([s], _cfg("nucleotide")) == (False, 11)


def test_kmer_basis_falls_back_to_nt_when_a_marker_is_missing(make_seq):
    a = make_seq("a", "ACGT" * 10)
    a.protein_sequence = "MKVL"
    b = make_seq("b", "ACGT" * 10)  # no protein_sequence populated
    assert kmer_basis([a, b], _cfg("protein")) == (False, 11)


def test_kmer_basis_default_alphabet_is_protein(make_seq):
    a = make_seq("a", "ACGT" * 10)
    a.protein_sequence = "MKVL"
    assert kmer_basis([a], {}) == (True, 5)  # missing key defaults to protein


def test_basis_sequence_switches_representation(make_seq):
    s = make_seq("a", "ACGTACGT")
    s.protein_sequence = "MKVL"
    assert basis_sequence(s, True) == "MKVL"
    assert basis_sequence(s, False) == "ACGTACGT"


def test_select_diverse_protein_basis_uses_marker(make_seq):
    """With cfg=protein, divergence is judged on the marker protein, not NT.

    NT says A==C and B differs; PROTEIN says A==B and C differs. A
    protein-basis MaxMin must therefore pick the protein-distinct C plus
    exactly one of the identical-protein pair {A, B}.
    """
    A = make_seq("A", "AAAAAAAAAAAAAAAA")
    A.protein_sequence = "MKKKKKKKKKKK"
    B = make_seq("B", "CCCCCCCCCCCCCCCC")
    B.protein_sequence = "MKKKKKKKKKKK"        # same protein as A, different NT
    C = make_seq("C", "AAAAAAAAAAAAAAAA")
    C.protein_sequence = "WWWWWWWWWWWW"        # same NT as A, different protein
    out = select_diverse([A, B, C], n=2, seed=42, cfg=_cfg("protein"))
    ids = {s.id for s in out}
    assert "C" in ids
    assert len(ids & {"A", "B"}) == 1


def test_select_diverse_cfgless_keeps_legacy_nt_behaviour(make_seq):
    """No cfg → historical NT/k=5 selection (back-compat)."""
    sim1 = make_seq("sim1", "AAAAAAAAAACCCCCCCCCC")
    sim2 = make_seq("sim2", "AAAAAAAAAACCCCCCCCCG")
    diff = make_seq("diff", "TTTTTTTTTTGGGGGGGGGG")
    out = select_diverse([sim1, sim2, diff], n=2, seed=42)
    assert "diff" in {s.id for s in out}
