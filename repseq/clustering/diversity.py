"""Max-diversity (MaxMin) sequence selection."""

from __future__ import annotations

import random
from typing import Optional

from ..models import Sequence


# ---------------------------------------------------------------------------
# k-mer distance (fast, alignment-free)
# ---------------------------------------------------------------------------

def _kmer_set(sequence: str, k: int = 5) -> set[str]:
    seq = sequence.upper()
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def _jaccard_distance(a: set, b: set) -> float:
    """Return 1 - Jaccard similarity (0 = identical, 1 = no overlap).

    Note: Jaccard divides by the *union*, so it is length-sensitive — a
    short sequence and a long one have an inflated distance purely from
    the size gap. Kept for the illustrative clustering plot; diversity
    selection uses :func:`_containment_distance` instead.
    """
    if not a and not b:
        return 0.0
    intersection = len(a & b)
    union = len(a | b)
    return 1.0 - intersection / union


def _containment_distance(a: set, b: set) -> float:
    """Length-robust k-mer distance: 1 - |A∩B| / min(|A|, |B|).

    Dividing by the *smaller* k-mer set (rather than the union, as Jaccard
    does) removes the length bias: a short sequence whose k-mers are all
    contained in a longer one scores ~0 ("not diverse from it") instead of
    ~1. Without this, MaxMin selection preferentially picks length
    extremes — the shortest and longest sequences — as "most diverse",
    which is an artefact, not biology. Symmetric because ``min`` is.
    """
    if not a or not b:
        return 0.0 if (not a and not b) else 1.0
    intersection = len(a & b)
    return 1.0 - intersection / min(len(a), len(b))


# ---------------------------------------------------------------------------
# MaxMin diversity selection
# ---------------------------------------------------------------------------

def select_diverse(
    sequences: list[Sequence],
    n: int,
    seed: Optional[int] = None,
    k: int = 5,
) -> list[Sequence]:
    """Select n sequences that maximally cover sequence space (MaxMin algorithm).

    Algorithm:
        1. Start with the longest sequence as the first representative.
        2. Maintain a per-sequence distance to the nearest already-selected representative.
        3. Iteratively pick the sequence with the greatest min-distance.

    Distance is the length-robust k-mer containment distance
    (:func:`_containment_distance`), so the selection reflects sequence
    divergence rather than length differences.

    Args:
        sequences: Pool of sequences to select from.
        n: Number of sequences to return.
        seed: Random seed for tie-breaking.
        k: k-mer size for distance computation.

    Returns:
        List of n (or fewer if pool is smaller) selected Sequence objects.
    """
    if not sequences:
        return []
    n = min(n, len(sequences))
    if n == len(sequences):
        return list(sequences)

    rng = random.Random(seed)

    # Pre-compute k-mer sets
    kmer_sets = [_kmer_set(s.sequence, k) for s in sequences]

    # Start with longest sequence
    first_idx = max(range(len(sequences)), key=lambda i: sequences[i].length)
    selected_indices = [first_idx]

    # min_dist[i] = distance from sequence i to its nearest selected representative
    min_dist = [_containment_distance(kmer_sets[i], kmer_sets[first_idx]) for i in range(len(sequences))]
    min_dist[first_idx] = -1.0  # mark as already selected

    for _ in range(n - 1):
        # Find candidate with maximum min-distance (tie-break randomly)
        max_d = max(min_dist)
        candidates = [i for i, d in enumerate(min_dist) if d == max_d]
        chosen = rng.choice(candidates)

        selected_indices.append(chosen)
        min_dist[chosen] = -1.0

        # Update min_dist for remaining sequences
        chosen_kmers = kmer_sets[chosen]
        for i in range(len(sequences)):
            if min_dist[i] < 0:
                continue
            d = _containment_distance(kmer_sets[i], chosen_kmers)
            if d < min_dist[i]:
                min_dist[i] = d

    return [sequences[i] for i in selected_indices]
