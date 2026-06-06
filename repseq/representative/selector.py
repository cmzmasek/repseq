"""Representative selection within a cluster.

The representative is the ``max`` of a **lexicographic** sort key built from
``representative.priority`` in the order the user listed it: the first
criterion dominates, each later one breaks ties, and sequence length is the
ultimate tiebreaker. So ``[refseq, reviewed_uniprot, longest]`` prefers a
RefSeq, then (among non-RefSeq) a Swiss-Prot-reviewed entry, then the longest;
``[reviewed_uniprot, refseq, longest]`` flips the first two. (Before v0.53.0
the score was additive with hard-wired weights, so the *order* of refseq vs
reviewed_uniprot was silently ignored — only their presence mattered.)
"""

from __future__ import annotations

from typing import Any

from ..models import Cluster, Sequence


# ---------------------------------------------------------------------------
# Priority scoring (lexicographic, order-honouring)
# ---------------------------------------------------------------------------

def _score(seq: Sequence, priority: list[str]) -> tuple:
    """Return a lexicographic sort key (larger tuple = better representative).

    Each listed criterion contributes one element, **in list order**, so the
    first criterion dominates and later ones break ties:
        ``refseq``           -> 1 if ``seq.is_refseq``   else 0
        ``reviewed_uniprot`` -> 1 if ``seq.is_reviewed`` else 0
        ``longest``          -> ``seq.length``
    ``seq.length`` is always appended last as the ultimate tiebreaker (the
    documented invariant), so selection stays deterministic even when
    ``longest`` is omitted. Unknown criteria can't reach here — config
    validation rejects them — and are skipped defensively.
    """
    key: list[int] = []
    for criterion in priority:
        if criterion == "refseq":
            key.append(1 if seq.is_refseq else 0)
        elif criterion == "reviewed_uniprot":
            key.append(1 if seq.is_reviewed else 0)
        elif criterion == "longest":
            key.append(seq.length)
    key.append(seq.length)  # ultimate tiebreaker, always present
    return tuple(key)


def select_representative(cluster_members: list[Sequence], priority: list[str]) -> Sequence:
    """Return the best representative from a list of sequences."""
    if not cluster_members:
        raise ValueError("Cannot select representative from empty list")
    return max(cluster_members, key=lambda s: _score(s, priority))


def apply_representative_selection(
    clusters: list[Cluster],
    cfg: dict[str, Any],
) -> list[Cluster]:
    """Re-select representatives in each cluster according to priority config.

    MMseqs2 picks longest by default; this step overrides with RefSeq/reviewed
    preference where applicable.
    """
    priority: list[str] = cfg.get("representative", {}).get(
        "priority", ["refseq", "reviewed_uniprot", "longest"]
    )

    for cluster in clusters:
        all_seqs = [cluster.representative] + cluster.members
        best = select_representative(all_seqs, priority)
        if best is not cluster.representative:
            # Swap: put old rep back in members list
            new_members = [s for s in all_seqs if s is not best]
            cluster.representative = best
            cluster.members = new_members

    return clusters
