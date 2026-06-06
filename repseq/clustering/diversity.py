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
# Representation + k-mer size (shared by diversity selection and the 2G plot)
# ---------------------------------------------------------------------------
#
# k=5 is well-resolved for PROTEIN markers (20^5 ≈ 3.2M possible 5-mers vs
# ~1000 per marker) but SATURATES on nucleotide whole genomes: the 4^5 = 1024
# 5-mer space is essentially fully present in any genome more than a few kb
# long, so every pairwise containment ≈ 1 and all distances collapse to ~0 —
# MaxMin then degenerates to near-random / length-ordered selection. So we
# mirror the clustering alphabet: when alphabet_for_clustering=protein (the
# default) and the reps carry the marker protein that drove clustering, k-mer
# work runs on those amino-acid strings at k=5; otherwise on the nucleotide
# sequence at a larger k. (For highly diverged genomes NO single nucleotide k
# is ideal — small k saturates, large k shares no k-mers across genera — so
# protein clustering is the recommended path for diverse viral families.)
_PROTEIN_KMER = 5
_NUCLEOTIDE_KMER = 11


def kmer_basis(
    sequences: list[Sequence], cfg: Optional[dict]
) -> tuple[bool, int]:
    """Return ``(use_protein, k)`` for k-mer work, from the clustering alphabet.

    ``use_protein`` is True only when ``clustering.alphabet_for_clustering``
    is ``protein`` AND every sequence carries a populated ``protein_sequence``
    (the marker / concat that actually drove clustering) — so diversity
    selection and the clustering plot see the same representation clustering
    did. Otherwise falls back to the nucleotide sequence with the larger
    nucleotide k.
    """
    alphabet = (cfg or {}).get("clustering", {}).get(
        "alphabet_for_clustering", "protein"
    )
    use_protein = (
        alphabet == "protein"
        and bool(sequences)
        and all(getattr(s, "protein_sequence", None) for s in sequences)
    )
    return use_protein, (_PROTEIN_KMER if use_protein else _NUCLEOTIDE_KMER)


def basis_sequence(seq: Sequence, use_protein: bool) -> str:
    """The string fed to k-mer extraction for ``seq`` under the chosen basis."""
    if use_protein:
        return seq.protein_sequence or ""
    return seq.sequence or ""


# ---------------------------------------------------------------------------
# MaxMin diversity selection
# ---------------------------------------------------------------------------

def select_diverse(
    sequences: list[Sequence],
    n: int,
    seed: Optional[int] = None,
    k: Optional[int] = None,
    cfg: Optional[dict] = None,
) -> list[Sequence]:
    """Select n sequences that maximally cover sequence space (MaxMin algorithm).

    Algorithm:
        1. Start with the longest sequence as the first representative.
        2. Maintain a per-sequence distance to the nearest already-selected representative.
        3. Iteratively pick the sequence with the greatest min-distance.

    Distance is the length-robust k-mer containment distance
    (:func:`_containment_distance`), so the selection reflects sequence
    divergence rather than length differences.

    Representation + k follow the clustering alphabet (see :func:`kmer_basis`):
    when ``cfg`` is supplied, a protein-alphabet run selects on the marker
    protein at k=5 and a nucleotide run on the genome at the larger nucleotide
    k. Without ``cfg`` the historical nucleotide/k=5 behaviour is kept. An
    explicit ``k`` always overrides the resolved default.

    Args:
        sequences: Pool of sequences to select from.
        n: Number of sequences to return.
        seed: Random seed for tie-breaking.
        k: k-mer size; ``None`` resolves from the alphabet (see above).
        cfg: full repseq config; enables alphabet-aware representation + k.

    Returns:
        List of n (or fewer if pool is smaller) selected Sequence objects.
    """
    if not sequences:
        return []
    n = min(n, len(sequences))
    if n == len(sequences):
        return list(sequences)

    rng = random.Random(seed)

    if cfg is not None:
        use_protein, basis_k = kmer_basis(sequences, cfg)
    else:
        use_protein, basis_k = False, _PROTEIN_KMER
    if k is None:
        k = basis_k

    # Pre-compute k-mer sets on the chosen representation.
    kmer_sets = [_kmer_set(basis_sequence(s, use_protein), k) for s in sequences]

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
