"""MMseqs2 clustering wrapper."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Optional

from ..models import Cluster, Sequence


class MMseqs2Error(RuntimeError):
    pass


def _check_mmseqs2() -> str:
    path = shutil.which("mmseqs")
    if not path:
        raise MMseqs2Error(
            "mmseqs2 not found in PATH. Install it from https://github.com/soedinglab/MMseqs2"
        )
    return path


def _write_id_fasta(
    sequences: list[Sequence],
    path: Path,
    line_width: int = 70,
    alphabet: str = "nucleotide",
) -> None:
    """Write a FASTA whose header is exactly ``seq.id``.

    MMseqs2 uses the first whitespace-delimited token of each header as the
    sequence identifier in its cluster TSV. Writing the full descriptive
    header (as the general-purpose writer does) would make that token differ
    from ``seq.id`` for UniProt (``sp|ACC|NAME ...``) and concatenated
    segmented isolates (``CONCAT|iso|acc1|acc2``), so the parsed cluster
    members could not be matched back. Using ``seq.id`` as the sole header
    token keeps the round-trip exact.

    With ``alphabet="protein"`` the body comes from ``seq.protein_sequence``
    (the marker protein, or concat thereof in segmented mode) — the
    upstream pipeline guarantees it's populated. Cluster objects returned
    by the backend still carry the original NT-bearing ``Sequence``.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for seq in sequences:
            body = seq.protein_sequence if alphabet == "protein" else seq.sequence
            if not body:
                raise ValueError(
                    f"Sequence {seq.id!r} has no "
                    f"{'protein_sequence' if alphabet == 'protein' else 'sequence'} "
                    f"to write to the clustering input."
                )
            # Defensive: a seq.id with a stray newline or carriage
            # return would split one record into two in the FASTA we
            # hand to MMseqs2, and the post-break fragment would then
            # appear as a phantom "sequence" with junk content.
            safe_id = seq.id.replace("\n", " ").replace("\r", " ")
            fh.write(f">{safe_id}\n")
            for i in range(0, len(body), line_width):
                fh.write(body[i : i + line_width] + "\n")


def run_clustering(
    sequences: list[Sequence],
    threshold: float,
    cfg: dict[str, Any],
    tmp_dir: Optional[str] = None,
) -> list[Cluster]:
    """Cluster sequences with MMseqs2 and return Cluster objects.

    Args:
        sequences: Input sequences (post-QC).
        threshold: Identity threshold (0.0–1.0).
        cfg: Full repseq config dict.
        tmp_dir: Override temp directory.

    Returns:
        List of Cluster objects, one per cluster.
    """
    mmseqs = _check_mmseqs2()
    cluster_cfg = cfg.get("clustering", {})
    mode = cluster_cfg.get("mmseqs2_mode", "easy-linclust")
    coverage = cluster_cfg.get("coverage", 0.8)
    coverage_mode = cluster_cfg.get("coverage_mode", 0)
    extra_args: list[str] = cluster_cfg.get("extra_args", [])
    alphabet = cluster_cfg.get("alphabet", "nucleotide")
    threads = cfg.get("threads", 4)

    work_dir = Path(tmp_dir or cfg.get("temp_dir", "/tmp/repseq"))
    work_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=work_dir, prefix="mmseqs_") as td:
        td = Path(td)
        input_fasta = td / "input.fasta"
        result_prefix = str(td / "result")
        mmseqs_tmp = str(td / "tmp")

        _write_id_fasta(sequences, input_fasta, alphabet=alphabet)

        cmd = [
            mmseqs,
            mode,
            str(input_fasta),
            result_prefix,
            mmseqs_tmp,
            "--min-seq-id", str(threshold),
            "-c", str(coverage),
            "--cov-mode", str(coverage_mode),
            "--threads", str(threads),
        ] + extra_args

        try:
            subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            raise MMseqs2Error(f"MMseqs2 failed:\n{e.stderr}") from e

        clusters = _parse_cluster_tsv(
            tsv_path=result_prefix + "_cluster.tsv",
            sequences=sequences,
        )

    # Defensive sanity check: every input sequence must appear in exactly
    # one cluster (rep or member). A mismatch means the cluster-TSV IDs
    # could not be matched back to the input — typically because the
    # input seq.id contained whitespace and MMseqs2 truncated it. Without
    # this check, _parse_cluster_tsv silently drops the unmatched rows
    # and the binary-search caller misreads the empty result as a
    # successful undershoot, returning all sequences at threshold = 1.0.
    accounted = sum(1 + len(c.members) for c in clusters)
    if accounted != len(sequences):
        raise MMseqs2Error(
            f"Cluster round-trip mismatch: {len(sequences)} sequences in, "
            f"{accounted} accounted for across {len(clusters)} clusters. "
            "Likely an ID/whitespace issue between input FASTA and the "
            "cluster TSV."
        )

    return clusters


def _parse_cluster_tsv(
    tsv_path: str,
    sequences: list[Sequence],
) -> list[Cluster]:
    """Parse MMseqs2 cluster TSV and return Cluster objects.

    TSV format: representative_id <TAB> member_id
    """
    seq_map: dict[str, Sequence] = {s.id: s for s in sequences}
    # Also index by accession in case IDs differ
    for s in sequences:
        if s.accession and s.accession not in seq_map:
            seq_map[s.accession] = s

    raw: dict[str, list[str]] = {}  # rep_id -> [member_ids]

    tsv = Path(tsv_path)
    if not tsv.exists():
        # MMseqs2 may produce _all_seqs.fasta instead; handle gracefully
        raise MMseqs2Error(f"Cluster TSV not found: {tsv_path}")

    with open(tsv) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rep_id, member_id = parts[0], parts[1]
            raw.setdefault(rep_id, [])
            if member_id != rep_id:
                raw[rep_id].append(member_id)

    clusters: list[Cluster] = []
    for i, (rep_id, member_ids) in enumerate(raw.items()):
        rep_seq = seq_map.get(rep_id)
        if rep_seq is None:
            continue
        members = [seq_map[m] for m in member_ids if m in seq_map]
        clusters.append(
            Cluster(
                cluster_id=f"cluster_{i+1:06d}",
                representative=rep_seq,
                members=members,
            )
        )

    return clusters
