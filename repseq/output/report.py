"""TSV report generation and run log."""

from __future__ import annotations

import copy
import datetime
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional

import yaml

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
    # Redact secrets before serialising the config into a plaintext log.
    safe_cfg = copy.deepcopy(cfg)
    tax_cfg = safe_cfg.get("taxonomy")
    if isinstance(tax_cfg, dict) and tax_cfg.get("ncbi_api_key"):
        tax_cfg["ncbi_api_key"] = "***redacted***"
    for line in yaml.dump(safe_cfg, default_flow_style=False, sort_keys=False).splitlines():
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
# Per-isolate protein TSV (segmented mode)
# ---------------------------------------------------------------------------

def _write_protein_fasta_record(
    fh,
    prot: dict,
    parent_seq: Sequence,
    isolate_id: Optional[str],
    line_width: int = 70,
) -> None:
    """Emit one FASTA record for a single protein."""
    pid = prot.get("protein_id") or "unknown"
    product = prot.get("product") or ""
    parent_acc = parent_seq.accession or parent_seq.id

    tags: list[str] = []
    if isolate_id:
        tags.append(f"[isolate={isolate_id}]")
    if parent_seq.segment:
        tags.append(f"[segment={parent_seq.segment}]")
    tags.append(f"[parent={parent_acc}]")

    header_parts = [f">{pid}"]
    if product:
        header_parts.append(product)
    header_parts.extend(tags)
    fh.write(" ".join(header_parts) + "\n")

    seq = prot["sequence"]
    for i in range(0, len(seq), line_width):
        fh.write(seq[i : i + line_width] + "\n")


def write_proteins_fasta(
    result: RunResult,
    complete_isolates: Optional[dict[str, list[Sequence]]],
    path: Path,
) -> bool:
    """Write all protein sequences associated with the selected representatives.

    Two paths:
      - Segmented mode (complete_isolates given): emits proteins from every
        segment of every isolate whose CONCAT representative was selected.
      - Non-segmented mode: emits proteins attached directly to each
        representative sequence.

    Skips the file (returns False) when nothing has a populated
    ``proteins[i]["sequence"]`` — e.g. protein QC didn't run, or cached
    entries pre-date the translation-capture change.
    """
    has_any_sequence = False

    if complete_isolates:
        # Selected isolates: extract from CONCAT|<isolate_id> IDs.
        rep_isolate_ids = {
            seq.isolate_id for seq in result.representatives if seq.isolate_id
        }
        for seq in result.representatives:
            if seq.id.startswith("CONCAT|"):
                parts = seq.id.split("|")
                if len(parts) > 1:
                    rep_isolate_ids.add(parts[1])

        targets: list[tuple[Optional[str], Sequence]] = []
        for iso_id, segs in complete_isolates.items():
            if rep_isolate_ids and iso_id not in rep_isolate_ids:
                continue
            for seq in segs:
                targets.append((iso_id, seq))
    else:
        targets = [(None, seq) for seq in result.representatives]

    for _, seq in targets:
        if seq.proteins and any(p.get("sequence") for p in seq.proteins):
            has_any_sequence = True
            break

    if not has_any_sequence:
        return False

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for iso_id, seq in targets:
            if not seq.proteins:
                continue
            for prot in seq.proteins:
                if not prot.get("sequence"):
                    continue
                _write_protein_fasta_record(fh, prot, seq, iso_id)
    return True


def write_isolate_proteins_tsv(
    complete_isolates: dict[str, list[Sequence]],
    path: Path,
) -> bool:
    """Write proteins per segment per isolate, one row per protein.

    Only emits a file when at least one segment has populated `proteins`.
    Returns True if the file was written, False if skipped.
    """
    has_any = any(
        seq.proteins for segs in complete_isolates.values() for seq in segs
    )
    if not has_any:
        return False

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(
            "isolate_id\tsegment\taccession\tprotein_id\tproduct\tlength\n"
        )
        for isolate_id, segs in complete_isolates.items():
            for seq in segs:
                if not seq.proteins:
                    continue
                for prot in seq.proteins:
                    fh.write(
                        f"{isolate_id}\t"
                        f"{seq.segment or ''}\t"
                        f"{seq.accession or seq.id}\t"
                        f"{prot.get('protein_id') or ''}\t"
                        f"{prot.get('product') or ''}\t"
                        f"{prot.get('length') if prot.get('length') is not None else ''}\n"
                    )
    return True


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
    complete_isolates: Optional[dict[str, list[Sequence]]] = None,
) -> None:
    out_dir = Path(cfg.get("output", {}).get("dir", "./repseq_output"))
    prefix = cfg.get("output", {}).get("prefix", "repseq")

    write_run_log(result, qc_report, cfg, input_paths, output_files, out_dir / f"{prefix}_run.log")
    write_qc_tsv(qc_report, out_dir / f"{prefix}_qc_removed.tsv")
    write_representative_tsv(result.representatives, out_dir / f"{prefix}_representatives.tsv")
    write_cluster_tsv(result, out_dir / f"{prefix}_clusters.tsv")
    if complete_isolates:
        write_isolate_proteins_tsv(
            complete_isolates, out_dir / f"{prefix}_isolate_proteins.tsv"
        )
    write_proteins_fasta(
        result, complete_isolates, out_dir / f"{prefix}_proteins.fasta"
    )
