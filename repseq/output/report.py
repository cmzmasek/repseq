"""TSV report generation and run log."""

from __future__ import annotations

import datetime
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional

from ..models import QCReport, RunResult, Sequence


# ---------------------------------------------------------------------------
# Run log
# ---------------------------------------------------------------------------

def write_run_log(
    result: RunResult,
    qc_report: QCReport,
    cfg: dict[str, Any],
    input_paths: list[str],
    output_files: list[Path],
    log_path: Path,
) -> None:
    """Write a human-readable run log with all parameters and QC stats."""
    log_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "=" * 72,
        "repseq run log",
        f"Date/time : {datetime.datetime.now().isoformat(timespec='seconds')}",
        f"Python    : {sys.version.split()[0]}",
        f"Platform  : {platform.platform()}",
        f"MMseqs2   : {_mmseqs2_version()}",
        "=" * 72,
        "",
        "INPUT",
        f"  Files   : {', '.join(input_paths)}",
        "",
        "CONFIGURATION",
    ]
    for line in json.dumps(cfg, indent=2, default=str).splitlines():
        lines.append(f"  {line}")
    lines += [
        "",
        "QC SUMMARY",
        qc_report.summary(),
        "",
        "RESULTS",
        f"  Mode              : {result.mode}",
        f"  Representatives   : {len(result.representatives)}",
        f"  Clusters          : {len(result.clusters)}",
        "",
        "OUTPUT FILES",
    ]
    for f in output_files:
        lines.append(f"  {f}")
    lines.append("")

    log_path.write_text("\n".join(lines))


def _mmseqs2_version() -> str:
    try:
        out = subprocess.check_output(["mmseqs", "version"], stderr=subprocess.DEVNULL, text=True)
        return out.strip()
    except Exception:
        return "not found"


# ---------------------------------------------------------------------------
# QC details TSV
# ---------------------------------------------------------------------------

def write_qc_tsv(qc_report: QCReport, path: Path) -> None:
    """Write per-sequence QC removal details to TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("sequence_id\treason\n")
        for entry in qc_report.details:
            fh.write(f"{entry['id']}\t{entry['reason']}\n")


# ---------------------------------------------------------------------------
# Representative metadata TSV
# ---------------------------------------------------------------------------

def write_representative_tsv(representatives: list[Sequence], path: Path) -> None:
    """Write metadata for all representative sequences to TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "accession", "organism", "description", "strain",
        "host", "collection_date", "country", "segment", "isolate_id",
        "seq_type", "length", "is_refseq", "is_reviewed",
        "taxid", "species", "genus", "family", "order",
    ]
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for seq in representatives:
            tax = seq.taxonomy
            row = [
                seq.accession or "",
                seq.organism or "",
                seq.description or "",
                seq.strain or "",
                seq.host or "",
                seq.collection_date or "",
                seq.country or "",
                seq.segment or "",
                seq.isolate_id or "",
                seq.seq_type.value,
                str(seq.length),
                str(seq.is_refseq).lower(),
                str(seq.is_reviewed).lower(),
                str(tax.taxid) if tax and tax.taxid else "",
                tax.species or "" if tax else "",
                tax.genus or "" if tax else "",
                tax.family or "" if tax else "",
                tax.order or "" if tax else "",
            ]
            fh.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# Cluster summary TSV
# ---------------------------------------------------------------------------

def write_cluster_tsv(result: RunResult, path: Path) -> None:
    """Write per-cluster summary to TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("cluster_id\trepresentative_id\trepresentative_organism\tcluster_size\tis_refseq\tis_reviewed\n")
        for cluster in result.clusters:
            rep = cluster.representative
            fh.write(
                f"{cluster.cluster_id}\t{rep.accession or rep.id}\t"
                f"{rep.organism or ''}\t{cluster.size}\t"
                f"{str(rep.is_refseq).lower()}\t{str(rep.is_reviewed).lower()}\n"
            )


# ---------------------------------------------------------------------------
# Write all reports
# ---------------------------------------------------------------------------

def write_all_reports(
    result: RunResult,
    qc_report: QCReport,
    cfg: dict[str, Any],
    input_paths: list[str],
    output_files: list[Path],
) -> None:
    out_dir = Path(cfg.get("output", {}).get("dir", "./repseq_output"))
    prefix = cfg.get("output", {}).get("prefix", "repseq")

    write_run_log(result, qc_report, cfg, input_paths, output_files, out_dir / f"{prefix}_run.log")
    write_qc_tsv(qc_report, out_dir / f"{prefix}_qc_removed.tsv")
    write_representative_tsv(result.representatives, out_dir / f"{prefix}_representatives.tsv")
    write_cluster_tsv(result, out_dir / f"{prefix}_clusters.tsv")
