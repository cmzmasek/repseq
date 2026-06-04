"""Protein-annotation QC: fetch CDS counts from NCBI, filter under-annotated sequences.

Runs after metadata resolution and before the segmented-virus completeness
filter, so segments lacking the expected protein count are removed before
they are counted toward isolate completeness.
"""

from __future__ import annotations

import logging
import sys
from typing import Any, Optional

from ..models import QCReport, Sequence, SequenceSource
from ..overrides import ProtectionPolicy, protected_keep
from ..segmented.completeness import identify_segment
from ..taxonomy.ncbi import NCBITaxonomy

logger = logging.getLogger(__name__)


def _stderr_batch_progress(done: int, total: int, _batch_size: int, cached: int) -> None:
    """Emit a single-line stderr heartbeat each time a GenBank batch lands.

    Prints the cache hit count on the first call (done=0), then ticks once
    per fetched batch. Plain stderr (not tqdm) because batches are slow
    and few, so a single status line per batch is the right granularity.
    """
    if done == 0:
        if total == 0:
            print(
                f"  [proteins] all {cached} accession(s) cached — no fetch needed.",
                file=sys.stderr,
            )
        else:
            print(
                f"  [proteins] {cached} cached, fetching {total} batch(es) of "
                f"up to 200 accession(s) from NCBI ...",
                file=sys.stderr,
            )
        return
    print(f"  [proteins] batch {done}/{total} done.", file=sys.stderr)


# ---------------------------------------------------------------------------
# Protein fetching (populate seq.proteins)
# ---------------------------------------------------------------------------

def attach_proteins(sequences: list[Sequence], ncbi: NCBITaxonomy) -> None:
    """Fetch CDS proteins for all NCBI-sourced sequences and attach in-place.

    UniProt sequences are skipped — they're protein records themselves and
    do not have a GenBank CDS feature table to count.
    """
    accessions: list[str] = []
    seq_by_acc: dict[str, list[Sequence]] = {}
    for seq in sequences:
        if seq.source == SequenceSource.UNIPROT:
            continue
        if not seq.accession:
            continue
        accessions.append(seq.accession)
        seq_by_acc.setdefault(seq.accession, []).append(seq)

    if not accessions:
        return

    proteins_by_acc = ncbi.fetch_proteins_batch(
        accessions, progress=_stderr_batch_progress,
    )
    for acc, proteins in proteins_by_acc.items():
        for seq in seq_by_acc.get(acc, []):
            seq.proteins = proteins


# ---------------------------------------------------------------------------
# Filter
# ---------------------------------------------------------------------------

def filter_by_protein_count(
    sequences: list[Sequence],
    qc_cfg: dict[str, Any],
    virus_cfg: Optional[dict[str, Any]],
    report: QCReport,
    *,
    policy: Optional[ProtectionPolicy] = None,
) -> list[Sequence]:
    """Drop sequences with insufficient protein annotations.

    Two checks are applied (skipped individually if not configured):

    1. ``qc.protein_annotation.min_proteins`` — global floor. A sequence with
       fewer than this many CDS features fails.
    2. ``segmented.viruses.<v>.expected_proteins_per_segment`` — per-segment
       count. Each value may be either an ``int`` (exact count required) or
       a ``list[int]`` (any of the listed counts is acceptable, useful when
       strain variation means a segment may legitimately have N or N+1
       proteins, e.g. PB1 with/without PB1-F2).

    Sequences whose ``proteins`` field is ``None`` (never fetched, e.g.
    UniProt-sourced or no accession) are passed through unchanged.
    """
    pa_cfg = qc_cfg.get("protein_annotation", {}) if qc_cfg else {}
    min_proteins = pa_cfg.get("min_proteins") if pa_cfg.get("enabled") else None

    expected_per_segment: dict[str, int] = {}
    segment_names: list[str] = []
    segment_regex: Optional[str] = None
    segment_aliases: Optional[dict[str, list[str]]] = None
    if virus_cfg:
        expected_per_segment = virus_cfg.get("expected_proteins_per_segment") or {}
        segment_names = virus_cfg.get("segments", [])
        segment_regex = virus_cfg.get("segment_regex")
        segment_aliases = virus_cfg.get("segment_aliases")

    if min_proteins is None and not expected_per_segment:
        return sequences  # nothing to do

    kept: list[Sequence] = []
    for seq in sequences:
        if seq.proteins is None:
            kept.append(seq)
            continue

        n = len(seq.proteins)
        fail_reason: Optional[str] = None

        if min_proteins is not None and n < min_proteins:
            fail_reason = f"protein_count_below_min:{n}<{min_proteins}"

        if not fail_reason and expected_per_segment and segment_names:
            seg = identify_segment(seq, segment_names, segment_regex, segment_aliases)
            if seg and seg in expected_per_segment:
                expected = expected_per_segment[seg]
                if isinstance(expected, list):
                    if n not in expected:
                        fail_reason = (
                            f"protein_count_mismatch:segment={seg}:got={n}:"
                            f"expected_one_of={expected}"
                        )
                else:
                    if n != expected:
                        fail_reason = (
                            f"protein_count_mismatch:segment={seg}:got={n}:expected={expected}"
                        )

        if fail_reason:
            if protected_keep(seq, "protein_count", fail_reason, policy, report):
                kept.append(seq)
                continue
            seq.qc_passed = False
            seq.qc_fail_reason = fail_reason
            report.removed_proteins += 1
            report.add_removed(seq.id, fail_reason)
        else:
            kept.append(seq)

    return kept


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def run_protein_qc(
    sequences: list[Sequence],
    ncbi: NCBITaxonomy,
    cfg: dict[str, Any],
    virus_cfg: Optional[dict[str, Any]],
    report: QCReport,
) -> list[Sequence]:
    """Fetch proteins for the batch and apply the protein-count filter.

    Returns the kept sequences. If neither min_proteins nor
    expected_proteins_per_segment is configured, this is a no-op.
    """
    qc_cfg = cfg.get("qc", {})
    pa_enabled = qc_cfg.get("protein_annotation", {}).get("enabled", False)
    has_per_segment = bool((virus_cfg or {}).get("expected_proteins_per_segment"))

    if not pa_enabled and not has_per_segment:
        return sequences

    attach_proteins(sequences, ncbi)
    policy = ProtectionPolicy.from_cfg(cfg)
    return filter_by_protein_count(
        sequences, qc_cfg, virus_cfg, report, policy=policy
    )
