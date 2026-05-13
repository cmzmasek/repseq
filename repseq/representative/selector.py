"""Representative selection within a cluster: RefSeq > Swiss-Prot > longest."""

from __future__ import annotations

from typing import Any

from ..models import Cluster, Sequence


# ---------------------------------------------------------------------------
# Priority scoring
# ---------------------------------------------------------------------------

_PRIORITY_SCORES = {
    "refseq": 100,
    "reviewed_uniprot": 50,
    "longest": 0,  # tiebreaker — always appended internally
}


def _score(seq: Sequence, priority: list[str]) -> tuple:
    """Return a sort key tuple (higher = better representative)."""
    score = 0
    for criterion in priority:
        if criterion == "refseq" and seq.is_refseq:
            score += _PRIORITY_SCORES["refseq"]
        elif criterion == "reviewed_uniprot" and seq.is_reviewed:
            score += _PRIORITY_SCORES["reviewed_uniprot"]
    # Always use length as final tiebreaker
    return (score, seq.length)


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
