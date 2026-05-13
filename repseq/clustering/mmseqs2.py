"""MMseqs2 clustering wrapper."""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Optional

from ..io.fasta import write_fasta
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
    threads = cfg.get("threads", 4)

    work_dir = Path(tmp_dir or cfg.get("temp_dir", "/tmp/repseq"))
    work_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=work_dir, prefix="mmseqs_") as td:
        td = Path(td)
        input_fasta = td / "input.fasta"
        result_prefix = str(td / "result")
        mmseqs_tmp = str(td / "tmp")

        write_fasta(sequences, input_fasta)

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
