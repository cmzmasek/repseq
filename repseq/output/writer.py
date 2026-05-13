"""FASTA output writer, including segmented virus multi-file output."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

from ..io.fasta import write_fasta
from ..models import RunResult, Sequence
from ..segmented.completeness import build_concatenated_sequences


def write_results(
    result: RunResult,
    cfg: dict[str, Any],
    complete_isolates: Optional[dict[str, list[Sequence]]] = None,
    segment_names: Optional[list[str]] = None,
) -> list[Path]:
    """Write representative sequences to FASTA files.

    For non-segmented runs: writes a single FASTA file.
    For segmented runs: writes one concatenated FASTA + one FASTA per segment.

    Returns:
        List of paths to written files.
    """
    out_dir = Path(cfg.get("output", {}).get("dir", "./repseq_output"))
    prefix = cfg.get("output", {}).get("prefix", "repseq")
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    if complete_isolates and segment_names:
        # Segmented: write concatenated + individual segment files
        written += _write_segmented(
            result, complete_isolates, segment_names, out_dir, prefix
        )
    else:
        # Standard: single output FASTA
        out_path = out_dir / f"{prefix}_representatives.fasta"
        write_fasta(result.representatives, out_path)
        written.append(out_path)

    return written


def _write_segmented(
    result: RunResult,
    complete_isolates: dict[str, list[Sequence]],
    segment_names: list[str],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    written: list[Path] = []

    # Determine which isolates are representatives
    rep_isolate_ids = {
        seq.isolate_id
        for seq in result.representatives
        if seq.isolate_id
    }
    # Also match by id prefix for concatenated sequences
    for seq in result.representatives:
        if seq.id.startswith("CONCAT|"):
            parts = seq.id.split("|")
            if len(parts) > 1:
                rep_isolate_ids.add(parts[1])

    # Filter isolates to only representative ones
    rep_isolates = {
        iso_id: segs
        for iso_id, segs in complete_isolates.items()
        if iso_id in rep_isolate_ids or not rep_isolate_ids
    }

    # Write concatenated sequences
    concat_seqs = build_concatenated_sequences(rep_isolates)
    concat_path = out_dir / f"{prefix}_concatenated.fasta"
    write_fasta(concat_seqs, concat_path)
    written.append(concat_path)

    # Write one file per segment
    for seg_idx, seg_name in enumerate(segment_names):
        seg_seqs: list[Sequence] = []
        for segs in rep_isolates.values():
            if seg_idx < len(segs):
                seg_seqs.append(segs[seg_idx])
        seg_path = out_dir / f"{prefix}_segment_{seg_name}.fasta"
        write_fasta(seg_seqs, seg_path)
        written.append(seg_path)

    return written
