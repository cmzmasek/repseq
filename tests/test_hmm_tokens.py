"""Tests for v0.14.0 HMM token notation and multidomain matching.

The token notation lives in repseq/hmm/runner.py:
    - ``parse_hmm_token(s)`` splits a token string into HMM names.
    - ``cds_satisfies_token(hits, names)`` returns the worst E-value if
      the CDS's hit list satisfies the (possibly multidomain) token, or
      ``None`` if not.

The C-to-N convention is unconventional (most molecular biology reads
N-to-C left-to-right); these tests pin it explicitly so accidental
flips during future refactors are caught immediately.
"""
from __future__ import annotations

import pytest

from repseq.hmm.runner import cds_satisfies_token, parse_hmm_token


def _hit(target, ali_from, ali_to, dom_evalue=1e-30, passing=True):
    return {
        "target": target,
        "passing": passing,
        "dom_evalue": dom_evalue,
        "ali_from": ali_from,
        "ali_to": ali_to,
        # ali_span included so the test hit shape matches the real parser.
        "ali_span": ali_to - ali_from + 1,
    }


# ---------------------------------------------------------------------------
# parse_hmm_token
# ---------------------------------------------------------------------------

def test_parse_hmm_token_single():
    assert parse_hmm_token("Bunya_G1") == ["Bunya_G1"]


def test_parse_hmm_token_multidomain_two():
    """C-to-N convention: 'A--B' has A C-terminal, B N-terminal."""
    assert parse_hmm_token("Bunya_G1--Bunya_G2") == ["Bunya_G1", "Bunya_G2"]


def test_parse_hmm_token_multidomain_three():
    assert parse_hmm_token("A--B--C") == ["A", "B", "C"]


def test_parse_hmm_token_strips_outer_whitespace():
    """Tolerate ``"A "`` or ``" A"`` so YAML quoting quirks don't bite."""
    assert parse_hmm_token("  A--B  ") == ["A", "B"]


def test_parse_hmm_token_preserves_internal_underscores():
    """Real Pfam names contain underscores; they must not be split."""
    assert parse_hmm_token("Bunya_nucleocap") == ["Bunya_nucleocap"]


def test_parse_hmm_token_rejects_empty_string():
    with pytest.raises(ValueError):
        parse_hmm_token("")
    with pytest.raises(ValueError):
        parse_hmm_token("   ")


def test_parse_hmm_token_rejects_empty_component():
    """A--B with a trailing or leading or internal empty component is a
    user typo, not a single-HMM token; reject it loudly."""
    with pytest.raises(ValueError):
        parse_hmm_token("A--")
    with pytest.raises(ValueError):
        parse_hmm_token("--B")
    with pytest.raises(ValueError):
        parse_hmm_token("A----B")


def test_parse_hmm_token_rejects_non_string():
    with pytest.raises(ValueError):
        parse_hmm_token(None)
    with pytest.raises(ValueError):
        parse_hmm_token(["A", "B"])


# ---------------------------------------------------------------------------
# cds_satisfies_token — single HMM
# ---------------------------------------------------------------------------

def test_satisfies_single_token_with_passing_hit():
    hits = [_hit("RdRP_4", 100, 300)]
    assert cds_satisfies_token(hits, ["RdRP_4"]) == 1e-30


def test_does_not_satisfy_single_token_when_no_passing_hit():
    """Failing hits (passing=False) don't count."""
    hits = [_hit("RdRP_4", 100, 300, passing=False)]
    assert cds_satisfies_token(hits, ["RdRP_4"]) is None


def test_does_not_satisfy_single_token_when_target_absent():
    hits = [_hit("Bunya_G1", 100, 300)]
    assert cds_satisfies_token(hits, ["RdRP_4"]) is None


def test_single_token_picks_best_evalue_when_multiple_hits():
    hits = [
        _hit("RdRP_4", 1, 100, dom_evalue=1e-20),
        _hit("RdRP_4", 200, 300, dom_evalue=1e-40),  # better
    ]
    assert cds_satisfies_token(hits, ["RdRP_4"]) == 1e-40


# ---------------------------------------------------------------------------
# cds_satisfies_token — multidomain with C-to-N ordering
# ---------------------------------------------------------------------------

def test_satisfies_multidomain_in_correct_order():
    """Token 'A--B' = A C-terminal to B. Protein has B at 1-300, A at
    500-800 → A.ali_from(500) > B.ali_to(300) → ordered correctly."""
    hits = [
        _hit("Bunya_G2", 1, 300, dom_evalue=1e-40),  # N-terminal
        _hit("Bunya_G1", 500, 800, dom_evalue=1e-30),  # C-terminal
    ]
    # Worst-of-set returned.
    assert cds_satisfies_token(hits, ["Bunya_G1", "Bunya_G2"]) == 1e-30


def test_multidomain_fails_when_order_is_reversed():
    """Same hits but the token names the WRONG direction: 'B--A' would
    mean B C-terminal to A. With B at 1-300 and A at 500-800, B is
    N-terminal to A, NOT C-terminal — token unsatisfied."""
    hits = [
        _hit("Bunya_G2", 1, 300),
        _hit("Bunya_G1", 500, 800),
    ]
    assert cds_satisfies_token(hits, ["Bunya_G2", "Bunya_G1"]) is None


def test_multidomain_fails_on_overlap():
    """Strict non-overlap required between named domains."""
    hits = [
        _hit("Bunya_G2", 1, 400),
        _hit("Bunya_G1", 300, 700),  # overlaps with G2
    ]
    assert cds_satisfies_token(hits, ["Bunya_G1", "Bunya_G2"]) is None


def test_multidomain_fails_when_one_hmm_missing():
    hits = [_hit("Bunya_G1", 500, 800)]
    assert cds_satisfies_token(hits, ["Bunya_G1", "Bunya_G2"]) is None


def test_multidomain_extras_are_allowed():
    """Extra unrelated HMMs hitting elsewhere on the protein are fine —
    the token only constrains the named domains."""
    hits = [
        _hit("Bunya_G2", 1, 300),
        _hit("Bunya_G1", 500, 800),
        _hit("HMMX", 900, 1100),  # extra, not in the token
    ]
    # Token "Bunya_G1--Bunya_G2" is still satisfied.
    assert cds_satisfies_token(hits, ["Bunya_G1", "Bunya_G2"]) == 1e-30


def test_multidomain_three_hmms_in_correct_order():
    """'A--B--C' = A most C-terminal, C most N-terminal. Protein layout
    N→C: C, B, A."""
    hits = [
        _hit("C", 1, 100),    # most N-terminal
        _hit("B", 200, 400),  # middle
        _hit("A", 500, 700),  # most C-terminal
    ]
    assert cds_satisfies_token(hits, ["A", "B", "C"]) is not None


def test_multidomain_three_hmms_with_middle_overlap_fails():
    hits = [
        _hit("C", 1, 100),
        _hit("B", 50, 400),   # overlaps with C
        _hit("A", 500, 700),
    ]
    assert cds_satisfies_token(hits, ["A", "B", "C"]) is None


def test_multidomain_picks_a_valid_assignment_when_one_hmm_has_multiple_hits():
    """If an HMM has multiple passing hits, the function must find ANY
    consistent assignment that satisfies the order. Here Bunya_G1 has
    two hits — one of them works."""
    hits = [
        _hit("Bunya_G2", 200, 400, dom_evalue=1e-40),
        _hit("Bunya_G1", 1, 100, dom_evalue=1e-20),    # too N-terminal
        _hit("Bunya_G1", 500, 800, dom_evalue=1e-25),  # correct slot
    ]
    # 'A--B' = G1 C-terminal to G2; the 500-800 G1 hit works.
    result = cds_satisfies_token(hits, ["Bunya_G1", "Bunya_G2"])
    assert result == 1e-25  # worst-of-set across chosen hits


def test_multidomain_returns_none_when_no_consistent_assignment():
    """G1 has hits only N-terminal to G2 — no valid C-to-N order."""
    hits = [
        _hit("Bunya_G2", 500, 800),
        _hit("Bunya_G1", 1, 100),    # N-terminal to G2, can't be C-term
        _hit("Bunya_G1", 200, 300),  # also N-terminal to G2
    ]
    assert cds_satisfies_token(hits, ["Bunya_G1", "Bunya_G2"]) is None


def test_satisfies_token_skips_failing_hits():
    """Even when a target HMM has a hit, it doesn't count if passing=False."""
    hits = [
        _hit("Bunya_G2", 1, 300, passing=True),
        _hit("Bunya_G1", 500, 800, passing=False),  # below cutoffs
    ]
    assert cds_satisfies_token(hits, ["Bunya_G1", "Bunya_G2"]) is None


def test_empty_token_list_returns_none():
    """A token with no HMMs is a degenerate case; return None rather
    than 'always satisfied' to avoid masking config bugs."""
    assert cds_satisfies_token([_hit("A", 1, 10)], []) is None


def test_empty_hits_with_nonempty_token():
    assert cds_satisfies_token([], ["RdRP_4"]) is None
    assert cds_satisfies_token(None, ["RdRP_4"]) is None
