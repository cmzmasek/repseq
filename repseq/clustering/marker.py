"""Marker-protein selection for amino-acid clustering.

When ``clustering.alphabet`` is ``protein``, every sequence must contribute
one amino-acid string to the clustering input. This module picks that
string from the sequence's CDS list (``seq.proteins``).

Selection rule:

* If an alias list is given and any /product matches one of the aliases
  (case-insensitive substring), the *first* alias in the list that matches
  any CDS wins — alias order encodes preference. Among multiple CDSes
  matching the same alias, the longest is kept (defensive against split
  /product strings).
* Otherwise the longest CDS with a /translation is returned.
* If no CDS has a /translation, returns ``None`` — the caller treats this
  as "no marker available" (drops the sequence, or the parent isolate).

Empty alias list is equivalent to no aliases (use longest CDS). ``None``
proteins (never fetched) returns ``None``.
"""

from __future__ import annotations

from typing import Optional

from ..models import QCReport, Sequence


def select_marker_protein(
    proteins: Optional[list[dict]],
    aliases: Optional[list[str]] = None,
) -> Optional[dict]:
    """Return the chosen CDS dict, or None if none qualify."""
    if not proteins:
        return None

    with_seq = [p for p in proteins if p.get("sequence")]
    if not with_seq:
        return None

    if aliases:
        for alias in aliases:
            needle = alias.lower().strip()
            if not needle:
                continue
            matches = [
                p for p in with_seq
                if needle in (p.get("product") or "").lower()
            ]
            if matches:
                return max(matches, key=lambda p: len(p["sequence"]))

    return max(with_seq, key=lambda p: len(p["sequence"]))


def populate_protein_sequences(
    sequences: list[Sequence],
    aliases: Optional[list[str]] = None,
    report: Optional[QCReport] = None,
) -> list[Sequence]:
    """Set ``seq.protein_sequence`` on each sequence to its marker protein.

    Used by non-segmented inputs only — segmented isolates get their
    concatenated marker via ``build_concatenated_sequences``.

    Sequences whose proteins list is ``None`` (never fetched) or has no
    CDS with a /translation are dropped and counted under
    ``report.removed_proteins`` with reason
    ``no_marker_protein_for_clustering``. Returns the surviving list.
    """
    kept: list[Sequence] = []
    for seq in sequences:
        marker = select_marker_protein(seq.proteins, aliases)
        if marker is None:
            if report is not None:
                reason = "no_marker_protein_for_clustering"
                seq.qc_passed = False
                seq.qc_fail_reason = reason
                report.removed_proteins += 1
                report.add_removed(seq.id, reason)
            continue
        seq.protein_sequence = marker["sequence"]
        if marker.get("protein_id"):
            seq.marker_protein_ids = [marker["protein_id"]]
        kept.append(seq)
    return kept
