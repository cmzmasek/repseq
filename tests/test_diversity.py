"""MaxMin diversity selection."""
from __future__ import annotations

from repseq.clustering.diversity import select_diverse


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
