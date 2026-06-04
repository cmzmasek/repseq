"""TSV report generation and run log."""

from __future__ import annotations

import copy
import datetime
import platform
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path
from typing import Any, Optional

import yaml

from ..config import get_virus_config
from ..hmm.runner import coverage_of
from ..models import QCReport, RunResult, Sequence


# A field value bound for a TSV row must contain neither tabs (column
# separator) nor line breaks (record separator); otherwise a single
# record silently splits into extras or drops/gains columns. Tabs and
# line breaks aren't expected in biological metadata, but they sneak in
# from copy-paste, malformed source files, and free-text NCBI fields.
_TSV_UNSAFE_RE = re.compile(r"[\t\r\n]+")


def _tsv_safe(value: Any) -> str:
    """Coerce a field value into a single TSV-safe cell.

    ``None`` becomes the empty string; tabs and line breaks collapse to
    a single space so the surrounding row keeps its column count.
    """
    if value is None:
        return ""
    return _TSV_UNSAFE_RE.sub(" ", str(value))


def _tsv_bool(value: Any) -> str:
    """Format a boolean for a TSV cell as ``TRUE`` / ``FALSE``.

    Uppercase chosen for R interop (R reads native ``TRUE``/``FALSE``
    without coercion) and for visual consistency across every TSV the
    program writes.
    """
    return "TRUE" if bool(value) else "FALSE"


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
        fh.write("accession\treason\n")
        for entry in qc_report.details:
            fh.write(f"{_tsv_safe(entry['id'])}\t{_tsv_safe(entry['reason'])}\n")


def write_overrides_tsv(qc_report: QCReport, path: Path) -> bool:
    """Write the force-keep audit TSV (overrides.protect_qc).

    One row per (sequence, stage) the override list rescued from removal:
    ``id``, the ``stage`` it was protected against, and the
    ``would_be_reason`` it would otherwise have been dropped for. Returns
    ``False`` (writing nothing) when no sequence was protected, so a run
    without overrides leaves no empty file behind.
    """
    if not qc_report.protected:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("id\tstage\twould_be_reason\n")
        for entry in qc_report.protected:
            fh.write(
                f"{_tsv_safe(entry['id'])}\t{_tsv_safe(entry['stage'])}\t"
                f"{_tsv_safe(entry['reason'])}\n"
            )
    return True


def write_force_selected_tsv(result: RunResult, path: Path) -> bool:
    """Write the force-select audit TSV (overrides.force_select).

    One row per pinned sequence: ``id``, the ``action`` taken
    (elected_representative / split_singleton / added_representative /
    already_representative / unavailable), and a ``detail`` string. Returns
    ``False`` (writing nothing) when force-select recorded nothing, so a run
    without the override leaves no empty file behind.
    """
    if not result.force_selected:
        return False
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("id\taction\tdetail\n")
        for entry in result.force_selected:
            fh.write(
                f"{_tsv_safe(entry['id'])}\t{_tsv_safe(entry['action'])}\t"
                f"{_tsv_safe(entry.get('detail', ''))}\n"
            )
    return True


# ---------------------------------------------------------------------------
# Representative metadata TSV
# ---------------------------------------------------------------------------

def write_representative_sequences_tsv(
    representatives: list[Sequence], path: Path
) -> None:
    """Write metadata for non-segmented representative sequences to TSV.

    One row per representative sequence. **Column-identical to**
    ``write_representative_isolates_tsv`` (the segmented rep table) so
    the two files are schema-compatible across modes — a downstream
    consumer reads the same columns regardless of which mode ran. The
    isolate-only columns are blanked or remapped to their per-sequence
    meaning:

    - ``isolate_id``, ``isolate_id_source``, ``n_segments``, ``segments``
      → blank (no isolate / multi-segment concept in non-segmented mode);
    - ``accessions`` → the sequence's single accession;
    - ``total_length_nt`` → the sequence's NT length.

    The per-sequence-only columns the old schema carried (``description``,
    ``segment``, ``molecule_type``, ``length_nt`` under that name) are
    intentionally absent — they have no slot in the shared isolates
    schema. Taxonomic ranks use the shared :data:`_TAX_RANKS` ladder;
    sub-ranks (``subgenus``, ``subfamily``, ``suborder``, ``subclass``)
    only populate when the resolver's lineage map carries them — commonly
    blank for viruses.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "isolate_id", "isolate_id_source", "organism", "strain", "host",
        "collection_date", "country", "n_segments", "segments", "accessions",
        "total_length_nt", "is_refseq", "is_reviewed", "ncbi_taxon_id",
        *_TAX_RANKS,
    ]
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for seq in representatives:
            tax = seq.taxonomy
            tax_cells = [
                _tsv_safe(tax.get_rank(r) if tax else None)
                for r in _TAX_RANKS
            ]
            row = [
                "",  # isolate_id — no isolate concept in non-segmented mode
                "",  # isolate_id_source
                _tsv_safe(seq.organism),
                _tsv_safe(seq.strain),
                _tsv_safe(seq.host),
                _tsv_safe(seq.collection_date),
                _tsv_safe(seq.country),
                "",  # n_segments
                "",  # segments
                _tsv_safe(seq.accession or seq.id),  # accessions (the one)
                _tsv_safe(seq.length),  # total_length_nt
                _tsv_bool(seq.is_refseq),
                _tsv_bool(seq.is_reviewed),
                _tsv_safe(tax.taxid) if tax and tax.taxid else "",
                *tax_cells,
            ]
            fh.write("\t".join(row) + "\n")


def write_representative_isolates_tsv(
    representatives: list[Sequence], path: Path
) -> None:
    """Write metadata for segmented representative isolates to TSV.

    One row per representative isolate (a synthetic
    ``CONCAT|<isolate_id>`` Sequence). Columns that have no isolate-level
    meaning (``accession``, ``segment``, ``description``,
    ``molecule_type``) are replaced by isolate-level equivalents derived
    from ``Sequence.concat_segments``: ``n_segments``, ``segments``
    (comma-joined segment names in concat order), ``accessions``
    (comma-joined per-segment GenBank accessions in concat order), and
    ``total_length_nt`` (sum of per-segment NT lengths). Other metadata
    (organism, strain, host, country, collection_date, taxonomy, RefSeq
    / reviewed flags) is the isolate's, populated by
    ``concatenate_isolate`` via first-non-empty / ``all()`` inheritance
    across segments. Same ``_TAX_RANKS`` ladder as the sequence TSV.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "isolate_id", "isolate_id_source", "organism", "strain", "host",
        "collection_date", "country", "n_segments", "segments", "accessions",
        "total_length_nt", "is_refseq", "is_reviewed", "ncbi_taxon_id",
        *_TAX_RANKS,
    ]
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for seq in representatives:
            tax = seq.taxonomy
            tax_cells = [
                _tsv_safe(tax.get_rank(r) if tax else None)
                for r in _TAX_RANKS
            ]
            segs = seq.concat_segments or []
            seg_names = [s.segment for s in segs if s.segment]
            seg_accs = [s.accession or s.id for s in segs]
            total_nt = sum(s.length for s in segs) if segs else seq.length
            row = [
                _tsv_safe(seq.isolate_id),
                _tsv_safe(seq.isolate_id_source),
                _tsv_safe(seq.organism),
                _tsv_safe(seq.strain),
                _tsv_safe(seq.host),
                _tsv_safe(seq.collection_date),
                _tsv_safe(seq.country),
                _tsv_safe(len(segs) if segs else ""),
                _tsv_safe(",".join(seg_names) if seg_names else ""),
                _tsv_safe(",".join(seg_accs) if seg_accs else ""),
                _tsv_safe(total_nt),
                _tsv_bool(seq.is_refseq),
                _tsv_bool(seq.is_reviewed),
                _tsv_safe(tax.taxid) if tax and tax.taxid else "",
                *tax_cells,
            ]
            fh.write("\t".join(row) + "\n")


# ---------------------------------------------------------------------------
# Per-isolate protein TSV (segmented mode)
# ---------------------------------------------------------------------------

_FASTA_HEADER_UNSAFE_RE = re.compile(r"[\[\]\r\n]+")


def _fasta_safe(value: Any) -> str:
    """Strip characters that would break FASTA-header parsing.

    Square brackets are the bracket-tag delimiter; line breaks split a
    header into two records. Both get collapsed to a single space so the
    tag value stays inline.
    """
    return _FASTA_HEADER_UNSAFE_RE.sub(" ", str(value)).strip()


def _write_protein_fasta_record(
    fh,
    prot: dict,
    parent_seq: Sequence,
    isolate_id: Optional[str],
    line_width: int = 70,
) -> None:
    """Emit one FASTA record for a single protein.

    Header format::

        >{protein_id} {product} [organism=...] [ncbi_taxon_id=...]
        [species=...] [genus=...] [family=...] [order=...] [class=...]
        [isolate=...] [segment=...] [host=...] [country=...]
        [collection_date=...] [length=...] [parent={parent_accession}]

    Tags are NCBI-style ``[key=value]`` and only emitted when the value is
    populated; empty fields are skipped to keep headers short on
    sparse-metadata sequences. Tag order goes from identity (organism +
    taxonomic lineage) to biological context (isolate, segment) to
    collection metadata (host, country, date) to technical (length) to
    provenance (parent). Taxonomy uses the same 9-rank ``_TAX_RANKS``
    ladder as the TSV writers — sub-ranks (``subgenus``, ``subfamily``,
    ``suborder``, ``subclass``) are commonly blank for viruses and are
    therefore omitted via the skip-empty rule.
    """
    pid = prot.get("protein_id") or "unknown"
    product = prot.get("product") or ""
    parent_acc = parent_seq.accession or parent_seq.id
    tax = parent_seq.taxonomy

    # Tag order is deliberate; see docstring. Taxonomy ranks are read
    # via TaxonomyInfo.get_rank() so they include the sub-ranks that
    # only exist in the lineage map.
    tag_specs: list[tuple[str, Any]] = [
        ("organism", parent_seq.organism),
        ("ncbi_taxon_id", tax.taxid if tax else None),
    ]
    if tax:
        tag_specs.extend(
            (rank, tax.get_rank(rank)) for rank in _TAX_RANKS
        )
    tag_specs.extend([
        ("isolate", isolate_id),
        ("segment", parent_seq.segment),
        ("host", parent_seq.host),
        ("country", parent_seq.country),
        ("collection_date", parent_seq.collection_date),
        ("length", prot.get("length")),
        # Matches the TSV ``hmmscan`` column exactly — same
        # ``Name(E=val,cov=val);...`` format, passing hits only,
        # best-E first. Empty (and thus skipped) when no HMM hits or
        # none passed.
        ("hmmscan", _format_hmmscan_cell(prot) or None),
        ("parent", parent_acc),
    ])
    tags = [
        f"[{key}={_fasta_safe(val)}]"
        for key, val in tag_specs
        if val not in (None, "")
    ]

    header_parts = [f">{pid}"]
    if product:
        header_parts.append(_fasta_safe(product))
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


def _pick_marker_cds(
    proteins: list[dict],
    hmm_tokens: list[str],
    aliases: list[str],
    hmm_active: bool,
    *,
    overlap_tolerance: int = 0,
) -> Optional[dict]:
    """Back-compat shim; delegates to :func:`phylo.per_protein.pick_marker_cds`.

    Both writers and the tree-builder share the same picker so a FASTA
    record and the corresponding tree leaf always carry the same CDS.
    """
    from ..phylo.per_protein import pick_marker_cds
    return pick_marker_cds(
        proteins, hmm_tokens, aliases, hmm_active,
        overlap_tolerance=overlap_tolerance,
    )


def _write_specs_to_fastas(
    specs: list[tuple[str, list[str], list[str], Optional[str]]],
    result: RunResult,
    cfg: dict[str, Any],
    sub_dir: Path,
    prefix: str,
    label: str,
) -> list[Path]:
    """Emit one unaligned protein FASTA per spec under ``sub_dir``.

    Shared engine for ``write_per_protein_fastas`` (cluster_protein /
    segment_markers) and ``write_extra_protein_fastas`` (extra_protein).
    Spec format: ``(family_label, hmm_tokens, aliases, segment_or_None)``.
    ``label`` is used only in the soft-fail stderr line.

    Filename: ``{sub_dir}/{prefix}_{family_label}.fasta`` — the prefix
    inside the basename matches the existing improvement-#1 pattern and
    survives a flat-copy out of the subdirectory unambiguously.
    """
    from ..phylo.per_protein import (
        _hmm_tier_ran,
        _segment_proteins,
        overlap_tolerance_from_cfg,
    )

    if not specs:
        return []

    hmm_active = _hmm_tier_ran(cfg, result.representatives)
    tol = overlap_tolerance_from_cfg(cfg)
    # Soft-fail: every spec needs HMM but the HMM tier didn't run.
    if not hmm_active and all(t and not a for _, t, a, _ in specs):
        sys.stderr.write(
            f"[{label} FASTA skipped] HMM tier did not run this session "
            f"and every configured spec is HMM-only — no FASTAs to write\n"
        )
        return []

    written: list[Path] = []
    for family_label, hmm_tokens, aliases, segment in specs:
        records: list[tuple[Optional[str], Sequence, dict]] = []
        for rep in result.representatives:
            # Pick the parent Sequence the CDS lives on so its segment/
            # accession/host/etc. drive the FASTA header tags.
            if segment is None:
                parent = rep
                isolate_id = rep.isolate_id
            else:
                parent = None
                for seg in rep.concat_segments or []:
                    if seg.segment == segment:
                        parent = seg
                        break
                if parent is None and rep.segment == segment:
                    parent = rep
                if parent is None:
                    continue
                isolate_id = parent.isolate_id or rep.isolate_id
                if not isolate_id and rep.id.startswith("CONCAT|"):
                    parts = rep.id.split("|")
                    if len(parts) > 1:
                        isolate_id = parts[1]

            proteins = _segment_proteins(rep, segment) if segment else (rep.proteins or [])
            cds = _pick_marker_cds(
                proteins, hmm_tokens, aliases, hmm_active,
                overlap_tolerance=tol,
            )
            if cds is None:
                continue
            records.append((isolate_id, parent, cds))

        if not records:
            continue

        sub_dir.mkdir(parents=True, exist_ok=True)
        path = sub_dir / f"{prefix}_{family_label}.fasta"
        with open(path, "w") as fh:
            for iso_id, parent_seq, cds in records:
                _write_protein_fasta_record(fh, cds, parent_seq, iso_id)
        written.append(path)

    return written


def write_per_protein_fastas(
    result: RunResult,
    cfg: dict[str, Any],
    complete_isolates: Optional[dict[str, list[Sequence]]],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    """Always-on unaligned protein FASTAs, one per declared marker spec.

    For each entry in ``clustering.cluster_protein`` (non-segmented) or
    ``virus.segment_markers`` / ``virus.cluster_protein`` (segmented), pick
    the satisfying CDS from each representative and write all such CDS to
    ``{prefix}_per_protein_fasta/{prefix}_{family}.fasta``. CDS selection
    uses the HMM gate when ``hmms:`` is declared AND the HMM tier ran
    (same logic as ``--per-protein-phylo``); for alias-only specs (or
    HMM-off runs on alias-bearing specs) it falls back to alias matching
    against ``/product``.

    Headers are byte-identical to the ``_representative_isolate_proteins.fasta``
    / ``_representative_sequence_proteins.fasta`` writer
    (:func:`_write_protein_fasta_record`), so a Spike record in the
    per-protein file and in the all-protein file carry the same tags.

    Specs that no representative satisfies (empty file) are skipped.
    When the HMM tier didn't run AND every spec is HMM-only, prints a
    soft-fail stderr note and returns ``[]``.
    """
    from ..phylo.per_protein import collect_marker_specs

    return _write_specs_to_fastas(
        collect_marker_specs(cfg),
        result, cfg,
        out_dir / f"{prefix}_per_protein_fasta",
        prefix, "per-protein",
    )


def write_extra_protein_fastas(
    result: RunResult,
    cfg: dict[str, Any],
    complete_isolates: Optional[dict[str, list[Sequence]]],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    """Always-on unaligned protein FASTAs for ``extra_protein`` specs.

    Parallel to :func:`write_per_protein_fastas` but driven by
    ``clustering.extra_protein`` (non-segmented) / ``virus.extra_protein``
    (segmented). These proteins are NOT used for clustering or the
    whole-genome tree — they may be sparse across representatives, and
    that is fine: a spec that no rep satisfies simply produces no file.

    Output lives under ``{prefix}_extra_protein_fasta/`` so it's visually
    distinct from the marker-driven per-protein FASTAs. Headers and CDS
    selection chain are identical (HMM-first, alias-fallback).
    """
    from ..phylo.per_protein import collect_extra_specs

    return _write_specs_to_fastas(
        collect_extra_specs(cfg),
        result, cfg,
        out_dir / f"{prefix}_extra_protein_fasta",
        prefix, "extra-protein",
    )


# ---------------------------------------------------------------------------
# Polyprotein-cutting outputs (v0.33.0+)
# ---------------------------------------------------------------------------

def _write_sliced_peptide_record(
    fh,
    sliced,  # SlicedPeptide
    spec_name: str,
    parent_seq: Sequence,
    isolate_id: Optional[str],
    line_width: int = 70,
) -> None:
    """Emit one FASTA record for a mature peptide.

    Headers reuse the existing protein-FASTA bracket-tag set
    (organism + 9-rank taxonomy + isolate + segment + host + country +
    collection_date + length + parent) and append polyprotein-specific
    tags so a downstream consumer can grep them: ``[polyprotein=...]``,
    ``[peptide=...]``, ``[peptide_range_aa=from-to]``, ``[cut_method=...]``.
    The ``>{header_id}`` is built as ``<parent_protein_id>:<peptide_name>``
    so two peptides of the same polyprotein on the same isolate have
    distinct first-token ids (no other writer cares; this just keeps
    things grep-friendly).
    """
    parent_pid = sliced.parent_protein_id or "unknown"
    pep_name = sliced.peptide_name or "peptide"
    parent_acc = sliced.parent_accession or parent_seq.accession or parent_seq.id
    tax = parent_seq.taxonomy

    tag_specs: list[tuple[str, Any]] = [
        ("organism", parent_seq.organism),
        ("ncbi_taxon_id", tax.taxid if tax else None),
    ]
    if tax:
        tag_specs.extend(
            (rank, tax.get_rank(rank)) for rank in _TAX_RANKS
        )
    tag_specs.extend([
        ("isolate", isolate_id),
        ("segment", parent_seq.segment),
        ("host", parent_seq.host),
        ("country", parent_seq.country),
        ("collection_date", parent_seq.collection_date),
        ("length", sliced.length_aa),
        ("polyprotein", spec_name),
        ("peptide", pep_name),
        (
            "peptide_range_aa",
            f"{sliced.range_aa_from}-{sliced.range_aa_to}"
            if sliced.range_aa_from and sliced.range_aa_to else None,
        ),
        ("cut_method", sliced.cut_method_actual or None),
        # When the peptide had multiple alternative architectures, this
        # surfaces which one matched (e.g. `aCoV_NSP1` vs `bCoV_NSP1`).
        # Empty/omitted when the peptide has only one architecture.
        ("matched_arch", sliced.matched_token or None),
        ("parent", parent_acc),
    ])
    tags = [
        f"[{key}={_fasta_safe(val)}]"
        for key, val in tag_specs
        if val not in (None, "")
    ]

    header_id = f"{parent_pid}:{pep_name}"
    header_parts = [f">{_fasta_safe(header_id)}", _fasta_safe(pep_name)]
    header_parts.extend(tags)
    fh.write(" ".join(p for p in header_parts if p) + "\n")

    seq = sliced.sequence
    for i in range(0, len(seq), line_width):
        fh.write(seq[i : i + line_width] + "\n")


def _polyprotein_proteins(rep: Sequence, segment: Optional[str]) -> list[dict]:
    """Proteins to search for ``rep`` under a polyprotein spec.

    Non-segmented (``segment is None``): the rep's own ``proteins``.
    Segmented: the matching entry in ``concat_segments`` (per-segment
    Sequence objects on a CONCAT rep); falls back to the rep itself
    when a non-CONCAT rep happens to carry that segment label.
    """
    if segment is None:
        return rep.proteins or []
    for seg in rep.concat_segments or []:
        if seg.segment == segment:
            return seg.proteins or []
    if rep.segment == segment:
        return rep.proteins or []
    return []


def _polyprotein_parent_seq(rep: Sequence, segment: Optional[str]) -> Sequence:
    """The Sequence whose metadata drives the peptide FASTA header tags.

    Mirrors the segment-scoping in :func:`_polyprotein_proteins`. Falls
    back to ``rep`` itself when the segment isn't found, so the writer
    never crashes on a malformed CONCAT rep — the header will just
    lack segment-specific metadata.
    """
    if segment is None:
        return rep
    for seg in rep.concat_segments or []:
        if seg.segment == segment:
            return seg
    return rep


def write_polyprotein_outputs(
    result: RunResult,
    cfg: dict[str, Any],
    out_dir: Path,
    prefix: str,
) -> list[Path]:
    """Slice each polyprotein spec across every representative.

    One FASTA per peptide of each spec at
    ``{prefix}_polyprotein/{prefix}_{spec_name}_{peptide_name}.fasta``,
    plus one audit TSV per spec at
    ``{prefix}_polyprotein/{prefix}_{spec_name}_peptides.tsv``. The
    audit TSV records every (rep × peptide) attempt with status (``ok``
    / ``missing`` / ``out_of_order`` / ``overlap`` / ``no_parent_cds``),
    so the bench scientist can see at a glance which reps were cleanly
    sliced and which carry caveats.

    Soft-fails (returns ``[]`` with a stderr line) when the HMM tier
    didn't run this session — no peptide HMMs can hit, nothing to do.
    Specs whose every rep yields ``no_parent_cds`` still write the
    audit TSV (so the user sees *why* no FASTAs appeared) but no FASTA.
    """
    from ..polyprotein import collect_polyprotein_specs, slice_polyprotein
    from ..phylo.per_protein import _hmm_tier_ran, overlap_tolerance_from_cfg

    specs = collect_polyprotein_specs(cfg)
    if not specs:
        return []

    hmm_active = _hmm_tier_ran(cfg, result.representatives)
    if not hmm_active:
        sys.stderr.write(
            "[polyprotein cutting skipped] HMM tier did not run this "
            "session — peptide HMMs can't be located on any CDS\n"
        )
        return []

    # Same multidomain overlap tolerance used by the cluster_protein
    # HMM gate (`hmm.multidomain_overlap_tolerance`, default 30 aa) so a
    # multidomain peptide like `CoV_RPol_N--RdRP_1` tolerates Pfam-
    # boundary fuzz at the seam just like the marker gate does.
    overlap_tol = overlap_tolerance_from_cfg(cfg)

    sub_dir = out_dir / f"{prefix}_polyprotein"
    written: list[Path] = []

    for spec in specs:
        records: list[tuple[Sequence, Optional[str], list]] = []
        for rep in result.representatives:
            proteins = _polyprotein_proteins(rep, spec.segment)
            parent_seq = _polyprotein_parent_seq(rep, spec.segment)
            isolate_id = parent_seq.isolate_id or rep.isolate_id
            if not isolate_id and rep.id.startswith("CONCAT|"):
                parts = rep.id.split("|")
                if len(parts) > 1:
                    isolate_id = parts[1]
            _parent_cds, sliced = slice_polyprotein(
                proteins, spec, overlap_tolerance=overlap_tol,
            )
            records.append((parent_seq, isolate_id, sliced))

        # The audit TSV is always written (even when zero sequences came
        # out) so the user sees the per-rep status decisions.
        sub_dir.mkdir(parents=True, exist_ok=True)
        audit_path = sub_dir / f"{prefix}_{spec.file_basename}_peptides.tsv"
        with open(audit_path, "w") as fh:
            fh.write(
                "isolate_id\tpeptide_name\tparent_accession\t"
                "parent_protein_id\trange_aa_from\trange_aa_to\t"
                "length_aa\tcut_method_actual\tmatched_architecture\t"
                "status\tnote\n"
            )
            for parent_seq, iso_id, sliced_list in records:
                effective_iso = iso_id or parent_seq.accession or parent_seq.id
                for s in sliced_list:
                    fh.write(
                        "\t".join([
                            _tsv_safe(effective_iso),
                            _tsv_safe(s.peptide_name),
                            _tsv_safe(s.parent_accession),
                            _tsv_safe(s.parent_protein_id),
                            _tsv_safe(s.range_aa_from or ""),
                            _tsv_safe(s.range_aa_to or ""),
                            _tsv_safe(s.length_aa or ""),
                            _tsv_safe(s.cut_method_actual),
                            _tsv_safe(s.matched_token),
                            _tsv_safe(s.status),
                            _tsv_safe(s.note),
                        ]) + "\n"
                    )
        written.append(audit_path)

        # One FASTA per declared peptide. We only emit a record when the
        # status is "ok" or "overlap" (overlap-affected slices are still
        # usable AA strings — the overlap is noted in the audit TSV).
        peptide_to_records: dict[str, list[tuple[Sequence, Optional[str], object]]] = {
            pep.name: [] for pep in spec.peptides
        }
        for parent_seq, iso_id, sliced_list in records:
            for s in sliced_list:
                if s.status in ("ok", "overlap") and s.sequence:
                    peptide_to_records[s.peptide_name].append((parent_seq, iso_id, s))

        for pep in spec.peptides:
            bucket = peptide_to_records.get(pep.name) or []
            if not bucket:
                continue
            safe_pep = _sanitize_for_path(pep.name)
            path = sub_dir / f"{prefix}_{spec.file_basename}_{safe_pep}.fasta"
            with open(path, "w") as fh:
                for parent_seq, iso_id, sliced in bucket:
                    _write_sliced_peptide_record(
                        fh, sliced, spec.name, parent_seq, iso_id,
                    )
            written.append(path)

    return written


_PATH_SAFE_RE = re.compile(r"[^A-Za-z0-9._-]+")


def _sanitize_for_path(name: str) -> str:
    """Filesystem-safe basename component for peptide filenames.

    Same shape as the per-protein writer's family-name sanitisation:
    keep alphanumerics + ``.`` ``_`` ``-`` (so a peptide name like
    ``NSP3.4`` survives), collapse anything else to ``_``.
    """
    cleaned = _PATH_SAFE_RE.sub("_", (name or "").strip()).strip("_")
    return cleaned or "peptide"


_TAX_RANKS = (
    "species",
    "subgenus",
    "genus",
    "subfamily",
    "family",
    "suborder",
    "order",
    "subclass",
    "class",
)


def _format_hmmscan_cell(prot: dict) -> str:
    """Render passing HMM hits as ``Name(E=val,cov=val);Name(E=val,cov=val)``.

    Sorted by domain E-value (best first). Hits not flagged ``passing``
    are excluded (the v0.13 design choice — show only what passed both
    cutoffs in the output, the raw hit list is internal). Returns ``""``
    when there are no passing hits or no ``hmm_hits`` field.
    """
    hits = prot.get("hmm_hits") or []
    passing = [h for h in hits if h.get("passing")]
    if not passing:
        return ""
    passing.sort(key=lambda h: h.get("dom_evalue", float("inf")))
    parts: list[str] = []
    for h in passing:
        name = h.get("target", "?")
        ev = h.get("dom_evalue")
        # HMM-model coverage (hmm_span / hmm_len), matching the QC gate —
        # not the protein alignment span.
        cov = coverage_of(h)
        parts.append(f"{name}(E={ev:.2g},cov={cov:.2f})")
    return ";".join(parts)


def write_isolate_proteins_tsv(
    complete_isolates: dict[str, list[Sequence]],
    path: Path,
    representative_isolate_ids: Optional[set[str]] = None,
) -> bool:
    """Write proteins per segment per isolate, one row per protein.

    Columns: ``protein_id``, ``product``, ``length_aa`` (protein length,
    amino acids), ``isolate_id``, ``segment``, ``segment_length_nt``
    (nucleotide length of the parent segment), ``accession``,
    ``representative`` (``TRUE`` if the isolate was selected as a
    clustering representative, ``FALSE`` otherwise), then the taxonomic
    ranks ``species``, ``subgenus``, ``genus``, ``subfamily``,
    ``family``, ``suborder``, ``order``, ``subclass``, ``class``.
    Sub-ranks are populated from the resolver's lineage map (NCBI
    ``LineageEx``) and are commonly empty for viruses, where ICTV often
    does not assign every intermediate rank.

    ``representative_isolate_ids`` is the set of isolate ids that survived
    clustering — typically computed by the caller from
    ``RunResult.representatives``. When ``None`` (legacy callers, direct
    unit tests), every row's ``representative`` cell is ``FALSE``.

    Only emits a file when at least one segment has populated `proteins`.
    Returns True if the file was written, False if skipped.
    """
    has_any = any(
        seq.proteins for segs in complete_isolates.values() for seq in segs
    )
    if not has_any:
        return False

    rep_ids = representative_isolate_ids or set()

    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "protein_id\tproduct\tlength_aa\tisolate_id\tisolate_id_source\t"
        "segment\tsegment_length_nt\taccession\trepresentative\thmmscan\t"
        + "\t".join(_TAX_RANKS) + "\n"
    )
    with open(path, "w") as fh:
        fh.write(header)
        for isolate_id, segs in complete_isolates.items():
            is_rep = _tsv_bool(isolate_id in rep_ids)
            # All segments share the same isolate_id_source by construction
            # (it's the provenance of the grouping key, which is the same
            # value for every segment of the isolate). Take the first
            # non-empty in case one segment was populated by the regex
            # fallback while the others weren't tagged yet.
            iso_source = next(
                (s.isolate_id_source for s in segs if s.isolate_id_source),
                None,
            )
            iso_source_cell = _tsv_safe(iso_source)
            for seq in segs:
                if not seq.proteins:
                    continue
                tax = seq.taxonomy
                tax_cells = "\t".join(
                    _tsv_safe(tax.get_rank(r) if tax else None)
                    for r in _TAX_RANKS
                )
                for prot in seq.proteins:
                    fh.write(
                        f"{_tsv_safe(prot.get('protein_id'))}\t"
                        f"{_tsv_safe(prot.get('product'))}\t"
                        f"{_tsv_safe(prot.get('length'))}\t"
                        f"{_tsv_safe(isolate_id)}\t"
                        f"{iso_source_cell}\t"
                        f"{_tsv_safe(seq.segment)}\t"
                        f"{_tsv_safe(seq.length)}\t"
                        f"{_tsv_safe(seq.accession or seq.id)}\t"
                        f"{is_rep}\t"
                        f"{_format_hmmscan_cell(prot)}\t"
                        f"{tax_cells}\n"
                    )
    return True


def write_sequence_proteins_tsv(
    sequences: list[Sequence],
    path: Path,
    representative_ids: Optional[set[str]] = None,
) -> bool:
    """Write proteins per non-segmented sequence, one row per protein.

    The non-segmented counterpart of :func:`write_isolate_proteins_tsv`,
    emitting the **identical column schema** so the per-CDS protein tables
    are joinable / processed the same way across modes:
    ``protein_id``, ``product``, ``length_aa``, ``isolate_id``,
    ``isolate_id_source``, ``segment``, ``segment_length_nt``,
    ``accession``, ``representative``, ``hmmscan``, then the nine
    :data:`_TAX_RANKS`.

    Columns with no non-segmented meaning are blanked: ``isolate_id``,
    ``isolate_id_source``, ``segment``. ``segment_length_nt`` is populated
    with the parent sequence's NT length (in non-segmented mode the whole
    sequence is the "segment"), and ``accession`` is the parent
    sequence's accession.

    ``representative_ids`` is the set of representative sequence ids
    (``Sequence.id``); a row's ``representative`` cell is ``TRUE`` when its
    parent sequence's id is in that set. When ``None`` every cell is
    ``FALSE``.

    Only emits a file when at least one sequence has populated
    ``proteins``. Returns True if written, False if skipped.
    """
    has_any = any(seq.proteins for seq in sequences)
    if not has_any:
        return False

    rep_ids = representative_ids or set()

    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "protein_id\tproduct\tlength_aa\tisolate_id\tisolate_id_source\t"
        "segment\tsegment_length_nt\taccession\trepresentative\thmmscan\t"
        + "\t".join(_TAX_RANKS) + "\n"
    )
    with open(path, "w") as fh:
        fh.write(header)
        for seq in sequences:
            if not seq.proteins:
                continue
            is_rep = _tsv_bool(seq.id in rep_ids)
            tax = seq.taxonomy
            tax_cells = "\t".join(
                _tsv_safe(tax.get_rank(r) if tax else None)
                for r in _TAX_RANKS
            )
            acc_cell = _tsv_safe(seq.accession or seq.id)
            seg_len_cell = _tsv_safe(seq.length)
            for prot in seq.proteins:
                fh.write(
                    f"{_tsv_safe(prot.get('protein_id'))}\t"
                    f"{_tsv_safe(prot.get('product'))}\t"
                    f"{_tsv_safe(prot.get('length'))}\t"
                    f"\t"  # isolate_id — no isolate concept here
                    f"\t"  # isolate_id_source
                    f"\t"  # segment
                    f"{seg_len_cell}\t"
                    f"{acc_cell}\t"
                    f"{is_rep}\t"
                    f"{_format_hmmscan_cell(prot)}\t"
                    f"{tax_cells}\n"
                )
    return True


# ---------------------------------------------------------------------------
# Cluster summary TSV
# ---------------------------------------------------------------------------

def _cluster_purity(values: list[Optional[str]]) -> tuple[Optional[float], int]:
    """Return ``(majority_fraction, n_distinct)`` for a cluster's per-member
    values at one rank.

    ``majority_fraction`` is the fraction of *populated* values that match
    the most common one (1.0 = every populated member agrees, < 1.0 =
    cluster spans multiple taxa at this rank). ``n_distinct`` counts how
    many distinct non-empty values appear; 1 = pure, > 1 = mixed.
    Members with empty/None values are excluded from both numerator and
    denominator — a member with no resolved genus shouldn't count against
    the cluster's genus purity. Returns ``(None, 0)`` when *no* member
    carries a value at this rank (the cluster is unclassifiable here).
    """
    populated = [v for v in values if v]
    if not populated:
        return None, 0
    counts = Counter(populated)
    top_count = max(counts.values())
    return top_count / len(populated), len(counts)


def write_cluster_tsv(result: RunResult, path: Path) -> None:
    """Write per-cluster summary to TSV with diagnostic purity columns.

    Columns: ``cluster_id``, ``accession`` (the representative's),
    ``organism``, ``cluster_size``, ``is_refseq``, ``is_reviewed``, then
    diagnostic columns computed across every member (representative +
    cluster.members):
        - ``species_purity`` / ``species_distinct`` — majority-species
          fraction and distinct-species count. ``species_purity < 1.0``
          means the cluster mixes species; ``species_distinct > 1`` is
          the same signal in raw form. Cells are blank when no member
          carries a species label.
        - ``genus_purity`` / ``genus_distinct`` — same at the genus rank.
        - ``host_purity`` / ``host_distinct`` — same at the host field
          (the literal ``/host`` qualifier, not taxonomy).
        - ``member_length_min`` / ``member_length_max`` /
          ``member_length_median`` — NT-length distribution of the
          cluster's members in nucleotides, useful for catching
          length-outlier inclusions.

    The purity columns let a bench scientist immediately spot clusters
    that span suspiciously many genera (likely the clustering threshold
    was too lax for that stratum) or hosts (cross-host contamination /
    mis-annotation). ``> 1`` distinct genera in a single cluster is
    almost always worth investigating; a single distinct species is the
    expected case at high identity thresholds.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(
            "cluster_id\taccession\torganism\tcluster_size\t"
            "is_refseq\tis_reviewed\t"
            "species_purity\tspecies_distinct\t"
            "genus_purity\tgenus_distinct\t"
            "host_purity\thost_distinct\t"
            "member_length_min\tmember_length_max\tmember_length_median\n"
        )
        for cluster in result.clusters:
            rep = cluster.representative
            members = [rep] + list(cluster.members)
            species = [
                (m.taxonomy.get_rank("species") if m.taxonomy else None)
                for m in members
            ]
            genera = [
                (m.taxonomy.get_rank("genus") if m.taxonomy else None)
                for m in members
            ]
            hosts = [m.host for m in members]
            sp_purity, sp_distinct = _cluster_purity(species)
            gn_purity, gn_distinct = _cluster_purity(genera)
            ho_purity, ho_distinct = _cluster_purity(hosts)
            lengths = [m.length for m in members if m.length]
            if lengths:
                import statistics
                len_min = str(min(lengths))
                len_max = str(max(lengths))
                len_med = str(int(round(statistics.median(lengths))))
            else:
                len_min = len_max = len_med = ""

            def _fmt_purity(p: Optional[float]) -> str:
                return f"{p:.3f}" if p is not None else ""

            fh.write(
                f"{_tsv_safe(cluster.cluster_id)}\t"
                f"{_tsv_safe(rep.accession or rep.id)}\t"
                f"{_tsv_safe(rep.organism)}\t{cluster.size}\t"
                f"{_tsv_bool(rep.is_refseq)}\t{_tsv_bool(rep.is_reviewed)}\t"
                f"{_fmt_purity(sp_purity)}\t{sp_distinct}\t"
                f"{_fmt_purity(gn_purity)}\t{gn_distinct}\t"
                f"{_fmt_purity(ho_purity)}\t{ho_distinct}\t"
                f"{len_min}\t{len_max}\t{len_med}\n"
            )


# ---------------------------------------------------------------------------
# Per-group selection counts TSV
# ---------------------------------------------------------------------------

def write_group_counts_tsv(result: RunResult, path: Path) -> bool:
    """Write per-group before/after selection counts to TSV.

    One row per stratum at whatever dimension the mode stratified on
    (taxonomic rank, host, time window, country, custom field, hybrid
    field combination, or the whole dataset for ``global``).
    ``stratum_size_before`` is the number of sequences entering the
    stratum — in segmented runs these are the concatenated per-isolate
    sequences. ``cutoff`` is the identity threshold the binary search
    settled on for that stratum, left blank when the stratum was small
    enough to keep without clustering (``clustered`` is then ``FALSE``).

    When ``clustering.diversity_curve_cutoffs`` is configured, each
    clustered stratum also carries cluster counts at those fixed
    cutoffs in trailing ``n_clusters_<c>`` columns (e.g.
    ``n_clusters_0.99``). Cells are ``NA`` for cutoffs below the
    backend's identity floor (cd-hit-est < 0.80, cd-hit protein <
    0.40); rows with ``clustered=FALSE`` leave every curve cell
    empty. Columns appear only when at least one row has
    ``cutoff_counts`` populated — the schema adapts to the data.

    Returns False without writing when the mode recorded no group stats.
    """
    if not result.group_stats:
        return False

    # Column set is the union of cutoffs seen across all rows, sorted
    # high → low so the most stringent threshold is leftmost. Reading
    # from the data (not from cfg) keeps the schema consistent with
    # what was actually computed for this run.
    curve_cutoffs: set[float] = set()
    for gs in result.group_stats:
        if gs.cutoff_counts:
            curve_cutoffs.update(gs.cutoff_counts.keys())
    sorted_cutoffs = sorted(curve_cutoffs, reverse=True)
    curve_cols = [f"n_clusters_{c:g}" for c in sorted_cutoffs]

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        header = (
            "stratified_by\tstratum\tstratum_size_before\t"
            "stratum_size_after\tclustered\tcutoff"
        )
        if curve_cols:
            header += "\t" + "\t".join(curve_cols)
        fh.write(header + "\n")
        for gs in result.group_stats:
            cutoff = f"{gs.cutoff:.4f}" if gs.cutoff is not None else ""
            row = (
                f"{_tsv_safe(gs.grouping)}\t{_tsv_safe(gs.group)}\t"
                f"{gs.n_before}\t{gs.n_after}\t"
                f"{_tsv_bool(gs.clustered)}\t{cutoff}"
            )
            if curve_cols:
                cells = []
                for c in sorted_cutoffs:
                    if gs.cutoff_counts is None:
                        cells.append("")
                    else:
                        v = gs.cutoff_counts.get(c)
                        cells.append("NA" if v is None else str(v))
                row += "\t" + "\t".join(cells)
            fh.write(row + "\n")
    return True


# ---------------------------------------------------------------------------
# Write all reports
# ---------------------------------------------------------------------------

def _seq_rank_value(seq: Sequence, rank: str) -> Optional[str]:
    """Return the taxon name a sequence carries at ``rank``, or ``None``.

    An empty string is treated as absent (not a taxon) so blank lineage
    cells never inflate a distinct-taxa count or appear in a breakdown.
    """
    tax = seq.taxonomy
    if not tax:
        return None
    return tax.get_rank(rank) or None


def write_taxonomic_report(
    before_seqs: list[Sequence],
    after_seqs: list[Sequence],
    segmented: bool,
    path: Path,
    max_breakdown: int = 20,
    provenance: Optional[str] = None,
) -> None:
    """Write ``{prefix}_taxonomic_report.txt``: taxonomic diversity at each
    of the nine :data:`_TAX_RANKS`, before vs after clustering.

    "Before" is the post-QC pool fed to the clustering step (one synthetic
    CONCAT isolate per record in segmented mode, one sequence per record
    otherwise); "after" is the selected representatives. The counting unit
    is therefore **isolates** in segmented mode and **sequences** otherwise.
    Empty / missing rank values are not a taxon and are excluded from every
    count.

    Two sections:

    1. Rank diversity — distinct non-empty taxa before and after, one row
       per rank.
    2. Per-rank breakdown — for every rank with at least one distinct
       taxon, each taxon name with its before / after unit count, sorted
       by member count (the *before* unit count) descending. Ranks with
       more than ``max_breakdown`` distinct taxa show only the top
       ``max_breakdown`` by member count, with a note in the rank label.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    unit = "isolates" if segmented else "sequences"

    before_by_rank: dict[str, Counter] = {}
    after_by_rank: dict[str, Counter] = {}
    for rank in _TAX_RANKS:
        before_by_rank[rank] = Counter(
            v for s in before_seqs if (v := _seq_rank_value(s, rank))
        )
        after_by_rank[rank] = Counter(
            v for s in after_seqs if (v := _seq_rank_value(s, rank))
        )

    lines: list[str] = []
    if provenance:
        lines.append(provenance)
    lines.append("Taxonomic report")
    lines.append(f"Generated: {datetime.date.today().isoformat()}")
    lines.append(
        f"Counting unit: {unit} "
        f"(before = post-QC pool fed to clustering, {len(before_seqs)} {unit}; "
        f"after = representatives, {len(after_seqs)} {unit})"
    )
    lines.append("")
    lines.append("== Rank diversity (distinct taxa before vs after clustering) ==")
    lines.append("")

    rank_w = max(len(r) for r in _TAX_RANKS)
    header = f"{'rank':<{rank_w}}  {'before':>7}  {'after':>7}"
    lines.append(header)
    lines.append("-" * len(header))
    for rank in _TAX_RANKS:
        nb = len(before_by_rank[rank])
        na = len(after_by_rank[rank])
        lines.append(f"{rank:<{rank_w}}  {nb:>7}  {na:>7}")
    lines.append("")

    lines.append(
        f"== Per-rank breakdown (per-taxon member counts, sorted by member "
        f"count; ranks with > {max_breakdown} distinct taxa show the top "
        f"{max_breakdown}) =="
    )
    lines.append("")
    printed_any = False
    for rank in _TAX_RANKS:
        before_c = before_by_rank[rank]
        n_distinct = len(before_c)
        if n_distinct == 0:
            continue
        printed_any = True
        after_c = after_by_rank[rank]
        names = sorted(before_c, key=lambda n: (-before_c[n], n))
        if n_distinct > max_breakdown:
            names = names[:max_breakdown]
            label = (
                f"{rank} ({n_distinct} distinct, top {max_breakdown} by "
                f"member count shown):"
            )
        else:
            label = f"{rank} ({n_distinct} distinct):"
        lines.append(label)
        name_w = max(max(len(n) for n in names), len("name"))
        sub_header = f"  {'name':<{name_w}}  {'before':>7}  {'after':>7}"
        lines.append(sub_header)
        lines.append("  " + "-" * (len(sub_header) - 2))
        for n in names:
            lines.append(
                f"  {n:<{name_w}}  {before_c[n]:>7}  {after_c.get(n, 0):>7}"
            )
        lines.append("")
    if not printed_any:
        lines.append("(no taxonomy available at any rank)")
        lines.append("")

    path.write_text("\n".join(lines).rstrip() + "\n")


# Ranks for the protein taxonomic report — `_TAX_RANKS` without ``species``.
# Species-level diversity is reported in `_taxonomic_report.txt`; the protein
# report starts at subgenus where coverage gaps between marker and accessory
# proteins actually start to inform biology (within a species you'd expect
# every isolate to carry every protein, modulo CDS-annotation noise).
_PROTEIN_REPORT_RANKS = tuple(r for r in _TAX_RANKS if r != "species")


def _protein_report_specs(
    cfg: dict[str, Any],
) -> list[tuple[str, list[str], list[str], Optional[str], bool]]:
    """Specs in the order the protein report's columns will appear.

    Returns ``[(family_label, hmm_tokens, aliases, segment_or_None,
    is_cluster), …]`` — cluster_protein / segment_markers entries first
    (rendered with a trailing ``*`` in the header to flag "drives the
    clustering and the whole-genome tree"), then extra_protein entries
    in declaration order.
    """
    from ..phylo.per_protein import collect_extra_specs, collect_marker_specs

    specs: list[tuple[str, list[str], list[str], Optional[str], bool]] = []
    for label, tokens, aliases, segment in collect_marker_specs(cfg):
        specs.append((label, tokens, aliases, segment, True))
    for label, tokens, aliases, segment in collect_extra_specs(cfg):
        specs.append((label, tokens, aliases, segment, False))
    return specs


def _proteins_at_segment(seq: Sequence, segment: Optional[str]) -> list[dict]:
    """Proteins to search for ``seq`` under a spec's ``segment`` scope.

    Non-segmented (``segment is None``): the sequence's own ``proteins``.
    Segmented: the matching entry in ``concat_segments``; falls back to
    the sequence itself if it directly carries that segment label
    (resilient against non-CONCAT segmented input).
    """
    if segment is None:
        return seq.proteins or []
    for sub in seq.concat_segments or []:
        if sub.segment == segment:
            return sub.proteins or []
    if seq.segment == segment:
        return seq.proteins or []
    return []


def _quartile_summary(values: list[int]) -> Optional[tuple[int, int, int, int, int]]:
    """``(min, max, median, Q3-Q1, n)`` for a list of ints, or ``None`` if empty.

    Q3-Q1 (IQR) is the third minus first quartile of the distribution.
    ``statistics.quantiles`` requires n>=2; for a single value the IQR is
    zero by definition (one-point distribution). ``n`` is the sample
    size so the user can judge how trustworthy the quartiles are.
    """
    if not values:
        return None
    import statistics

    vmin = min(values)
    vmax = max(values)
    median = int(round(statistics.median(values)))
    if len(values) >= 2:
        q1, _q2, q3 = statistics.quantiles(values, n=4)
        iqr = int(round(q3 - q1))
    else:
        iqr = 0
    return vmin, vmax, median, iqr, len(values)


def _format_coverage_cell(count: int, total: int) -> str:
    """Render one (count, total) coverage cell as ``"<count> <pct>%"``.

    ``total == 0`` is the "no items in this taxon" case (shouldn't happen
    in practice — taxa are derived from the same items being counted —
    but guarded so a division-by-zero can't void the report).
    """
    if total <= 0:
        return f"{count} ---"
    pct = round(100.0 * count / total)
    return f"{count} {pct}%"


def _format_length_cell(stats: Optional[tuple[int, int, int, int, int]]) -> str:
    """Render a ``(min, max, median, Q3-Q1, n)`` cell, or ``"n/a"`` when empty."""
    if stats is None:
        return "n/a"
    vmin, vmax, median, iqr, n = stats
    return f"{vmin}, {vmax}, {median}, {iqr}, {n}"


def _coverage_data_per_taxon(
    seqs: list[Sequence],
    specs: list[tuple[str, list[str], list[str], Optional[str], bool]],
    hmm_active: bool,
    rank: str,
    overlap_tolerance: int = 0,
) -> tuple[dict[str, int], dict[str, list[list[int]]]]:
    """For each taxon at ``rank``, count items and gather per-spec AA lengths.

    Returns ``(taxon_totals, taxon_lengths)``:

    * ``taxon_totals[taxon]`` = number of items in that taxon (the
      denominator of the coverage percentage).
    * ``taxon_lengths[taxon][spec_idx]`` = list of AA lengths of the
      satisfying CDS, one per item that carries the spec. Items without
      a satisfying CDS contribute *nothing* to the list (so the count
      from ``len(taxon_lengths[t][i])`` is the numerator of coverage).

    Items whose taxonomy lacks a value at ``rank`` are excluded (mirrors
    ``_taxonomic_report.txt`` behaviour). The picker chain is the same
    one used for the per-protein FASTAs and trees, so a coverage hit
    here corresponds exactly to a record in the matching artifact.
    """
    from ..phylo.per_protein import pick_marker_cds

    taxon_totals: dict[str, int] = {}
    taxon_lengths: dict[str, list[list[int]]] = {}
    n_specs = len(specs)

    for seq in seqs:
        taxon = _seq_rank_value(seq, rank)
        if not taxon:
            continue
        taxon_totals[taxon] = taxon_totals.get(taxon, 0) + 1
        if taxon not in taxon_lengths:
            taxon_lengths[taxon] = [[] for _ in range(n_specs)]
        for i, (_label, tokens, aliases, segment, _is_cluster) in enumerate(specs):
            proteins = _proteins_at_segment(seq, segment)
            cds = pick_marker_cds(
                proteins, tokens, aliases, hmm_active,
                overlap_tolerance=overlap_tolerance,
            )
            if cds is None:
                continue
            seq_str = cds.get("sequence") or ""
            length = cds.get("length") or len(seq_str)
            if not length:
                continue
            taxon_lengths[taxon][i].append(int(length))

    return taxon_totals, taxon_lengths


def _format_coverage_table(
    indent: str,
    title: str,
    totals: dict[str, int],
    lengths: dict[str, list[list[int]]],
    spec_headers: list[str],
    max_breakdown: int,
) -> list[str]:
    """Render one coverage sub-table.

    Sorted by ``total`` desc (matches ``_taxonomic_report.txt``); ranks
    with more than ``max_breakdown`` distinct taxa get a ``"+N more taxa
    not shown"`` summary line so the column widths stay terminal-friendly.
    """
    lines: list[str] = [f"{indent}{title}"]
    # Truncate by member count (the rank header advertises "top N by
    # member count") so a long-tail rank doesn't drop high-count taxa
    # just for being late in the alphabet; then sort the kept taxa
    # alphabetically for display.
    taxa = sorted(totals, key=lambda t: (-totals[t], t))
    truncated = 0
    if len(taxa) > max_breakdown:
        truncated = len(taxa) - max_breakdown
        taxa = taxa[:max_breakdown]
    taxa = sorted(taxa)
    if not taxa:
        lines.append(f"{indent}  (no taxa)")
        return lines

    cells: dict[tuple[str, int], str] = {}
    for taxon in taxa:
        total = totals[taxon]
        for i, _h in enumerate(spec_headers):
            count = len(lengths[taxon][i])
            cells[(taxon, i)] = _format_coverage_cell(count, total)

    taxon_w = max(max(len(t) for t in taxa), len("taxon"))
    count_w = max(max(len(str(totals[t])) for t in taxa), len("count"))
    spec_ws: list[int] = []
    for i, header in enumerate(spec_headers):
        w = max(
            len(header),
            *(len(cells[(t, i)]) for t in taxa),
        )
        spec_ws.append(w)

    sub_indent = indent + "  "
    header_parts = [f"{'taxon':<{taxon_w}}", f"{'count':>{count_w}}"]
    for i, header in enumerate(spec_headers):
        header_parts.append(f"{header:>{spec_ws[i]}}")
    header_line = "  ".join(header_parts)
    lines.append(sub_indent + header_line)
    lines.append(sub_indent + "-" * len(header_line))
    for taxon in taxa:
        row_parts = [f"{taxon:<{taxon_w}}", f"{totals[taxon]:>{count_w}}"]
        for i in range(len(spec_headers)):
            row_parts.append(f"{cells[(taxon, i)]:>{spec_ws[i]}}")
        lines.append(sub_indent + "  ".join(row_parts))
    if truncated:
        lines.append(sub_indent + f"... +{truncated} more taxa not shown")
    return lines


def _format_length_table(
    indent: str,
    title: str,
    totals: dict[str, int],
    lengths: dict[str, list[list[int]]],
    spec_headers: list[str],
    max_breakdown: int,
) -> list[str]:
    """Render one length-statistics sub-table (``min, max, median, Q3-Q1, n``)."""
    lines: list[str] = [f"{indent}{title}"]
    # Same truncate-then-alphabetise pattern as the coverage table —
    # see _format_coverage_table for rationale.
    taxa = sorted(totals, key=lambda t: (-totals[t], t))
    truncated = 0
    if len(taxa) > max_breakdown:
        truncated = len(taxa) - max_breakdown
        taxa = taxa[:max_breakdown]
    taxa = sorted(taxa)
    if not taxa:
        lines.append(f"{indent}  (no taxa)")
        return lines

    cells: dict[tuple[str, int], str] = {}
    for taxon in taxa:
        for i in range(len(spec_headers)):
            cells[(taxon, i)] = _format_length_cell(
                _quartile_summary(lengths[taxon][i])
            )

    taxon_w = max(max(len(t) for t in taxa), len("taxon"))
    spec_ws: list[int] = []
    for i, header in enumerate(spec_headers):
        w = max(
            len(header),
            *(len(cells[(t, i)]) for t in taxa),
        )
        spec_ws.append(w)

    sub_indent = indent + "  "
    header_parts = [f"{'taxon':<{taxon_w}}"]
    for i, header in enumerate(spec_headers):
        header_parts.append(f"{header:>{spec_ws[i]}}")
    header_line = "  ".join(header_parts)
    lines.append(sub_indent + header_line)
    lines.append(sub_indent + "-" * len(header_line))
    for taxon in taxa:
        row_parts = [f"{taxon:<{taxon_w}}"]
        for i in range(len(spec_headers)):
            row_parts.append(f"{cells[(taxon, i)]:>{spec_ws[i]}}")
        lines.append(sub_indent + "  ".join(row_parts))
    if truncated:
        lines.append(sub_indent + f"... +{truncated} more taxa not shown")
    return lines


def _format_hmm_architecture_section(
    specs: list[tuple[str, list[str], list[str], Optional[str], bool]],
    cfg: dict[str, Any],
) -> list[str]:
    """Architecture + cutoff-policy summary for HMM-bearing specs."""
    hmm_specs = [s for s in specs if s[1]]  # non-empty hmm_tokens
    if not hmm_specs:
        return []

    hmm_cfg = cfg.get("hmm", {}) or {}
    use_ga = bool(hmm_cfg.get("use_ga_when_available", True))
    default_e = hmm_cfg.get("default_evalue", 1.0e-5)
    rel_len = hmm_cfg.get("relative_length_cutoff", 0.5)
    cutoff_summary = (
        f"GA cutoffs when available, else E≤{default_e:g}"
        if use_ga else f"E≤{default_e:g}"
    )

    lines: list[str] = []
    lines.append("== HMM marker architectures ==")
    lines.append("")
    lines.append(
        "Tokens within a spec are alternative architectures (OR): a CDS "
        "satisfying ANY token satisfies the spec. Multidomain tokens "
        "(\"A--B--C\") list HMMs in N-to-C order on the protein."
    )
    lines.append("")
    for label, tokens, _aliases, segment, is_cluster in hmm_specs:
        star = "*" if is_cluster else ""
        scope = f" (segment {segment})" if segment else ""
        joined = " OR ".join(tokens) if len(tokens) > 1 else tokens[0]
        lines.append(f"{label}{star}{scope}: {joined}")
        lines.append(f"  ({cutoff_summary}; coverage ≥ {rel_len:g})")
    lines.append("")
    return lines


def write_protein_taxonomic_report(
    before_seqs: list[Sequence],
    after_seqs: list[Sequence],
    cfg: dict[str, Any],
    segmented: bool,
    path: Path,
    provenance: Optional[str] = None,
) -> bool:
    """Write ``{prefix}_protein_taxonomic_report.txt``.

    For each rank from ``subgenus`` up to ``class`` (skipping
    ``species`` — within-species coverage gaps are dominated by
    annotation noise, not biology), emits two coverage sub-tables
    (post-QC pool + representatives) and two length-statistics
    sub-tables for every declared protein. Cluster-driving markers
    (``cluster_protein`` / ``segment_markers``) get a trailing ``*``
    in the header — they participate in clustering AND the
    whole-genome tree. ``extra_protein`` entries follow.

    Returns False (writes nothing) when no specs are configured —
    the file would carry only structural noise. The picker chain
    matches the per-protein/extra-protein FASTAs and trees so a
    coverage hit here corresponds exactly to a record in those
    artifacts.
    """
    from ..phylo.per_protein import _hmm_tier_ran, overlap_tolerance_from_cfg

    specs = _protein_report_specs(cfg)
    if not specs:
        return False

    pr_cfg = cfg.get("output", {}).get("protein_report", {}) or {}
    max_breakdown = int(pr_cfg.get("max_breakdown", 20))

    hmm_active = _hmm_tier_ran(cfg, after_seqs)
    tol = overlap_tolerance_from_cfg(cfg)
    spec_headers = [
        f"{label}{'*' if is_cluster else ''}"
        for label, _t, _a, _s, is_cluster in specs
    ]
    unit = "isolates" if segmented else "sequences"

    lines: list[str] = []
    if provenance:
        lines.append(provenance)
    lines.append("Protein taxonomic report")
    lines.append(f"Generated: {datetime.date.today().isoformat()}")
    lines.append(
        f"Counting unit: {unit} "
        f"(post-QC pool fed to clustering = {len(before_seqs)} {unit}; "
        f"representatives = {len(after_seqs)} {unit})"
    )
    lines.append(
        "Coverage cell format: <count> <%>; length cell format: "
        "min, max, median, Q3-Q1, n (amino acids; n = number of "
        "items contributing the length)."
    )
    lines.append(
        "An asterisk (*) marks proteins also used for clustering and the "
        "whole-genome tree; the rest are accessory proteins declared via "
        "`extra_protein:` and reported here only."
    )
    lines.append("")

    rank_printed_any = False
    for rank in _PROTEIN_REPORT_RANKS:
        b_totals, b_lengths = _coverage_data_per_taxon(
            before_seqs, specs, hmm_active, rank, overlap_tolerance=tol,
        )
        a_totals, a_lengths = _coverage_data_per_taxon(
            after_seqs, specs, hmm_active, rank, overlap_tolerance=tol,
        )
        n_distinct = len(b_totals)
        if n_distinct == 0:
            continue
        rank_printed_any = True
        # Use the same rank label style as _taxonomic_report.txt so the
        # two reports key on the same headings ("subgenus (18 distinct):").
        if n_distinct > max_breakdown:
            header = (
                f"{rank} ({n_distinct} distinct, top {max_breakdown} by "
                f"member count shown):"
            )
        else:
            header = f"{rank} ({n_distinct} distinct):"
        lines.append(header)
        lines.extend(_format_coverage_table(
            indent="  ",
            title=f"coverage (post-QC pool, N={sum(b_totals.values())})",
            totals=b_totals, lengths=b_lengths,
            spec_headers=spec_headers,
            max_breakdown=max_breakdown,
        ))
        lines.append("")
        if a_totals:
            lines.extend(_format_coverage_table(
                indent="  ",
                title=f"coverage (representatives, N={sum(a_totals.values())})",
                totals=a_totals, lengths=a_lengths,
                spec_headers=spec_headers,
                max_breakdown=max_breakdown,
            ))
            lines.append("")
        lines.extend(_format_length_table(
            indent="  ",
            title="protein length statistics [min, max, median, Q3-Q1, n] (post-QC pool)",
            totals=b_totals, lengths=b_lengths,
            spec_headers=spec_headers,
            max_breakdown=max_breakdown,
        ))
        lines.append("")
        if a_totals:
            lines.extend(_format_length_table(
                indent="  ",
                title="protein length statistics [min, max, median, Q3-Q1, n] (representatives)",
                totals=a_totals, lengths=a_lengths,
                spec_headers=spec_headers,
                max_breakdown=max_breakdown,
            ))
            lines.append("")

    if not rank_printed_any:
        lines.append("(no taxonomy available at any rank from "
                     "subgenus to class)")
        lines.append("")

    lines.extend(_format_hmm_architecture_section(specs, cfg))

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n")
    return True


def _polyprotein_coverage_data_per_taxon(
    seqs: list[Sequence],
    spec,
    rank: str,
    overlap_tolerance: int,
) -> tuple[dict[str, int], dict[str, list[list[int]]]]:
    """Per-taxon coverage + length data for one polyprotein spec.

    Slices every sequence under the spec and aggregates by taxon at
    ``rank``. ``taxon_lengths[taxon][peptide_idx]`` is the list of
    peptide AA lengths for the items that yielded a usable slice
    (``status in {"ok", "overlap"}`` — the same statuses
    ``write_polyprotein_outputs`` actually writes to the per-peptide
    FASTAs, so coverage here corresponds to records on disk). Items
    whose taxonomy lacks a value at ``rank`` are excluded, mirroring
    ``_taxonomic_report.txt`` and the protein report.
    """
    from ..polyprotein import slice_polyprotein

    taxon_totals: dict[str, int] = {}
    taxon_lengths: dict[str, list[list[int]]] = {}
    n_peptides = len(spec.peptides)
    pep_index = {pep.name: i for i, pep in enumerate(spec.peptides)}

    for seq in seqs:
        taxon = _seq_rank_value(seq, rank)
        if not taxon:
            continue
        taxon_totals[taxon] = taxon_totals.get(taxon, 0) + 1
        if taxon not in taxon_lengths:
            taxon_lengths[taxon] = [[] for _ in range(n_peptides)]
        proteins = _polyprotein_proteins(seq, spec.segment)
        if not proteins:
            continue
        _parent, sliced = slice_polyprotein(
            proteins, spec, overlap_tolerance=overlap_tolerance,
        )
        for s in sliced:
            if s.status not in ("ok", "overlap"):
                continue
            length = s.length_aa or 0
            if not length:
                continue
            idx = pep_index.get(s.peptide_name)
            if idx is None:
                continue
            taxon_lengths[taxon][idx].append(int(length))

    return taxon_totals, taxon_lengths


def write_polyprotein_taxonomic_report(
    before_seqs: list[Sequence],
    after_seqs: list[Sequence],
    cfg: dict[str, Any],
    segmented: bool,
    path: Path,
    provenance: Optional[str] = None,
) -> bool:
    """Write ``{prefix}_polyprotein_taxonomic_report.txt``.

    Per-rank coverage + length sub-tables for the mature peptides of
    every declared ``polyprotein:`` spec — the sliced-peptide analogue
    of :func:`write_protein_taxonomic_report`. Each spec gets its own
    section; columns are the declared peptides in N→C order, kept as
    columns even when zero items carry the peptide (a declared-but-
    never-located peptide is itself a useful audit signal).

    Returns False (writes nothing) when:

    * no ``polyprotein:`` spec is declared, OR
    * the HMM tier didn't run this session (peptide HMMs can't be
      located on any CDS — same soft-fail posture as
      :func:`write_polyprotein_outputs`).

    The picker chain (and the "what counts as covered" rule —
    ``status in {"ok", "overlap"}``) matches the per-peptide FASTAs,
    so a coverage hit here corresponds exactly to a record in the
    matching artifact under ``{prefix}_polyprotein/``.
    """
    from ..polyprotein import collect_polyprotein_specs
    from ..phylo.per_protein import _hmm_tier_ran, overlap_tolerance_from_cfg

    specs = collect_polyprotein_specs(cfg)
    if not specs:
        return False
    if not _hmm_tier_ran(cfg, after_seqs):
        return False

    pr_cfg = cfg.get("output", {}).get("protein_report", {}) or {}
    max_breakdown = int(pr_cfg.get("max_breakdown", 20))
    tol = overlap_tolerance_from_cfg(cfg)
    unit = "isolates" if segmented else "sequences"

    lines: list[str] = []
    if provenance:
        lines.append(provenance)
    lines.append("Polyprotein taxonomic report")
    lines.append(f"Generated: {datetime.date.today().isoformat()}")
    lines.append(
        f"Counting unit: {unit} "
        f"(post-QC pool fed to clustering = {len(before_seqs)} {unit}; "
        f"representatives = {len(after_seqs)} {unit})"
    )
    lines.append(
        "Coverage cell format: <count> <%>; length cell format: "
        "min, max, median, Q3-Q1, n (amino acids; n = number of "
        "items contributing the length)."
    )
    lines.append(
        "Coverage counts peptides with status `ok` or `overlap` — the "
        "same statuses written to the per-peptide FASTAs under "
        "`{prefix}_polyprotein/`. Peptides with status `missing`, "
        "`out_of_order`, or `no_parent_cds` are not counted as covered."
    )
    lines.append("")

    for spec in specs:
        scope = f" (segment {spec.segment})" if spec.segment else ""
        lines.append(f"========== {spec.name}{scope} ==========")
        lines.append("")
        peptide_headers = [pep.name for pep in spec.peptides]

        rank_printed_any = False
        for rank in _PROTEIN_REPORT_RANKS:
            b_totals, b_lengths = _polyprotein_coverage_data_per_taxon(
                before_seqs, spec, rank, overlap_tolerance=tol,
            )
            a_totals, a_lengths = _polyprotein_coverage_data_per_taxon(
                after_seqs, spec, rank, overlap_tolerance=tol,
            )
            n_distinct = len(b_totals)
            if n_distinct == 0:
                continue
            rank_printed_any = True
            if n_distinct > max_breakdown:
                header = (
                    f"{rank} ({n_distinct} distinct, top {max_breakdown} by "
                    f"member count shown):"
                )
            else:
                header = f"{rank} ({n_distinct} distinct):"
            lines.append(header)
            lines.extend(_format_coverage_table(
                indent="  ",
                title=f"coverage (post-QC pool, N={sum(b_totals.values())})",
                totals=b_totals, lengths=b_lengths,
                spec_headers=peptide_headers,
                max_breakdown=max_breakdown,
            ))
            lines.append("")
            if a_totals:
                lines.extend(_format_coverage_table(
                    indent="  ",
                    title=f"coverage (representatives, N={sum(a_totals.values())})",
                    totals=a_totals, lengths=a_lengths,
                    spec_headers=peptide_headers,
                    max_breakdown=max_breakdown,
                ))
                lines.append("")
            lines.extend(_format_length_table(
                indent="  ",
                title="peptide length statistics [min, max, median, Q3-Q1, n] (post-QC pool)",
                totals=b_totals, lengths=b_lengths,
                spec_headers=peptide_headers,
                max_breakdown=max_breakdown,
            ))
            lines.append("")
            if a_totals:
                lines.extend(_format_length_table(
                    indent="  ",
                    title="peptide length statistics [min, max, median, Q3-Q1, n] (representatives)",
                    totals=a_totals, lengths=a_lengths,
                    spec_headers=peptide_headers,
                    max_breakdown=max_breakdown,
                ))
                lines.append("")

        if not rank_printed_any:
            lines.append("(no taxonomy available at any rank from "
                         "subgenus to class)")
            lines.append("")

        # Architecture readout per spec: each peptide's OR-alternatives
        # joined with " OR " (so a bench scientist reading the report
        # can tell at a glance what was declared as the locator for
        # each row), plus the cleavage motif when set.
        lines.append("== Peptide architectures ==")
        lines.append("")
        for pep in spec.peptides:
            tokens = pep.hmms or []
            joined = " OR ".join(tokens) if len(tokens) > 1 else (
                tokens[0] if tokens else "(no HMM token declared)"
            )
            motif = f"  [cleavage_motif={pep.cleavage_motif}]" if pep.cleavage_motif else ""
            lines.append(f"{pep.name}: {joined}{motif}")
        lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n")
    return True


# ---------------------------------------------------------------------------
# Tidy long-format TSV companions to the four `_*_report.txt` files.
# Shared 8-column schema across all four so a user can concat them in
# pandas / R / Excel pivot tables and key on the `report` column to
# distinguish sources. Rows are pure observations — one numeric per
# row, named by `metric`. No multi-table layout, no headers, no banner
# rows. The `.txt` files remain the primary humans-first format; these
# `.tsv` files are the machine-readable companion.

# The full canonical column set, in fixed order. Every writer below
# emits exactly this header; downstream consumers can rely on it.
_TIDY_TSV_COLUMNS = (
    "report", "rank", "pool", "taxon", "taxon_count", "spec", "metric", "value",
)


def _emit_tidy_rows(
    fh,
    report: str,
    rank: str,
    pool: str,
    taxon: str,
    taxon_count: int,
    spec: str,
    metric_value_pairs: list[tuple[str, Any]],
) -> None:
    """Write one TSV row per (metric, value) pair sharing the same key.

    All ints are written as decimal integers; floats are formatted to
    four decimal places trimmed of trailing zeros so a typical
    coverage_pct like ``91.6667`` stays exact-ish without `.0` noise
    on round percentages. Strings (taxon names, spec names) are passed
    through ``_tsv_safe`` to scrub embedded tabs/newlines.
    """
    for metric, value in metric_value_pairs:
        if isinstance(value, float):
            v = f"{value:.4f}".rstrip("0").rstrip(".") if value else "0"
        elif isinstance(value, int):
            v = str(value)
        else:
            v = _tsv_safe(value)
        fh.write(
            "\t".join([
                _tsv_safe(report),
                _tsv_safe(rank),
                _tsv_safe(pool),
                _tsv_safe(taxon),
                str(int(taxon_count)),
                _tsv_safe(spec),
                _tsv_safe(metric),
                v,
            ]) + "\n"
        )


def _emit_length_rows(
    fh,
    report: str,
    rank: str,
    pool: str,
    taxon: str,
    taxon_count: int,
    spec: str,
    lengths: list[int],
) -> None:
    """Emit length_n and (when n>0) length_min/max/median/iqr rows.

    Zero-length lists emit ``length_n=0`` and nothing else — we won't
    fake min/max/median values for a taxon that didn't contribute any
    length to this spec, since downstream aggregators would treat them
    as real zeros.
    """
    stats = _quartile_summary(lengths)
    if stats is None:
        _emit_tidy_rows(
            fh, report, rank, pool, taxon, taxon_count, spec,
            [("length_n", 0)],
        )
        return
    vmin, vmax, median, iqr, n = stats
    _emit_tidy_rows(
        fh, report, rank, pool, taxon, taxon_count, spec, [
            ("length_min", vmin),
            ("length_max", vmax),
            ("length_median", median),
            ("length_iqr", iqr),
            ("length_n", n),
        ],
    )


def write_taxonomic_report_tsv(
    before_seqs: list[Sequence],
    after_seqs: list[Sequence],
    path: Path,
) -> bool:
    """Tidy long-format TSV companion to ``_taxonomic_report.txt``.

    One row per (rank, pool, taxon, metric). Two metrics only:

    * ``distinct_taxa`` — number of distinct populated values at that
      rank in that pool. Emitted with ``taxon='*ALL*'`` and
      ``taxon_count`` equal to the value.
    * ``member_count`` — per-taxon item count (the same numbers shown
      in the per-rank breakdown of the `.txt`). One row per
      (rank, pool, taxon) actually populated; emits both pools for
      every taxon seen in *either* pool (so reps-only rare survivors
      and pool-only drops both show up).

    All nine ``_TAX_RANKS`` are covered. Empty / missing rank values
    are excluded from every count (same rule as the `.txt`).
    """
    before_by_rank: dict[str, Counter] = {}
    after_by_rank: dict[str, Counter] = {}
    for rank in _TAX_RANKS:
        before_by_rank[rank] = Counter(
            v for s in before_seqs if (v := _seq_rank_value(s, rank))
        )
        after_by_rank[rank] = Counter(
            v for s in after_seqs if (v := _seq_rank_value(s, rank))
        )

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(_TIDY_TSV_COLUMNS) + "\n")
        for rank in _TAX_RANKS:
            nb = len(before_by_rank[rank])
            na = len(after_by_rank[rank])
            _emit_tidy_rows(
                fh, "diversity", rank, "post_qc", "*ALL*", nb,
                "_diversity", [("distinct_taxa", nb)],
            )
            _emit_tidy_rows(
                fh, "diversity", rank, "reps", "*ALL*", na,
                "_diversity", [("distinct_taxa", na)],
            )
            # Per-taxon member counts — union of pool keys so neither
            # pool's exclusive members are lost.
            names: set[str] = set(before_by_rank[rank]) | set(after_by_rank[rank])
            for name in sorted(names):
                nb_c = before_by_rank[rank].get(name, 0)
                na_c = after_by_rank[rank].get(name, 0)
                _emit_tidy_rows(
                    fh, "diversity", rank, "post_qc", name, nb_c,
                    "_diversity", [("member_count", nb_c)],
                )
                _emit_tidy_rows(
                    fh, "diversity", rank, "reps", name, na_c,
                    "_diversity", [("member_count", na_c)],
                )
    return True


def write_protein_taxonomic_report_tsv(
    before_seqs: list[Sequence],
    after_seqs: list[Sequence],
    cfg: dict[str, Any],
    path: Path,
) -> bool:
    """Tidy long-format TSV companion to ``_protein_taxonomic_report.txt``.

    One row per (rank, pool, taxon, spec, metric). Metrics:
    ``coverage_count``, ``coverage_pct``, ``length_min``,
    ``length_max``, ``length_median``, ``length_iqr``, ``length_n``.
    Length metrics other than ``length_n`` are emitted only when at
    least one item contributed a length (so a 0% coverage row carries
    ``coverage_count=0, coverage_pct=0, length_n=0`` and nothing else
    — no fake zeros that would skew downstream aggregations).
    Returns False (writes nothing) when no protein specs are configured.
    """
    from ..phylo.per_protein import _hmm_tier_ran, overlap_tolerance_from_cfg

    specs = _protein_report_specs(cfg)
    if not specs:
        return False

    hmm_active = _hmm_tier_ran(cfg, after_seqs)
    tol = overlap_tolerance_from_cfg(cfg)

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(_TIDY_TSV_COLUMNS) + "\n")
        for rank in _PROTEIN_REPORT_RANKS:
            b_totals, b_lengths = _coverage_data_per_taxon(
                before_seqs, specs, hmm_active, rank, overlap_tolerance=tol,
            )
            a_totals, a_lengths = _coverage_data_per_taxon(
                after_seqs, specs, hmm_active, rank, overlap_tolerance=tol,
            )
            for pool_label, totals, lengths in (
                ("post_qc", b_totals, b_lengths),
                ("reps", a_totals, a_lengths),
            ):
                for taxon in sorted(totals):
                    total = totals[taxon]
                    for i, (label, _t, _a, _s, is_cluster) in enumerate(specs):
                        spec_name = f"{label}*" if is_cluster else label
                        ls = lengths[taxon][i]
                        cov_count = len(ls)
                        cov_pct = (100.0 * cov_count / total) if total else 0.0
                        _emit_tidy_rows(
                            fh, "protein", rank, pool_label, taxon, total,
                            spec_name, [
                                ("coverage_count", cov_count),
                                ("coverage_pct", cov_pct),
                            ],
                        )
                        _emit_length_rows(
                            fh, "protein", rank, pool_label, taxon, total,
                            spec_name, ls,
                        )
    return True


def write_polyprotein_taxonomic_report_tsv(
    before_seqs: list[Sequence],
    after_seqs: list[Sequence],
    cfg: dict[str, Any],
    path: Path,
) -> bool:
    """Tidy long-format TSV companion to ``_polyprotein_taxonomic_report.txt``.

    Same schema as the protein tidy TSV; ``spec`` carries the
    ``<polyprotein_name>:<peptide_name>`` composite (e.g. ``ORF1ab:NSP1``)
    so the spec column stays a single string while still distinguishing
    peptides across multiple specs. Same metric set + same
    length-skipping rule (``length_n=0`` only on 0-coverage rows).
    Returns False when no ``polyprotein:`` spec is declared OR the HMM
    tier didn't run (parity with the `.txt` writer).
    """
    from ..polyprotein import collect_polyprotein_specs
    from ..phylo.per_protein import _hmm_tier_ran, overlap_tolerance_from_cfg

    specs = collect_polyprotein_specs(cfg)
    if not specs:
        return False
    if not _hmm_tier_ran(cfg, after_seqs):
        return False
    tol = overlap_tolerance_from_cfg(cfg)

    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(_TIDY_TSV_COLUMNS) + "\n")
        for spec in specs:
            for rank in _PROTEIN_REPORT_RANKS:
                b_totals, b_lengths = _polyprotein_coverage_data_per_taxon(
                    before_seqs, spec, rank, overlap_tolerance=tol,
                )
                a_totals, a_lengths = _polyprotein_coverage_data_per_taxon(
                    after_seqs, spec, rank, overlap_tolerance=tol,
                )
                for pool_label, totals, lengths in (
                    ("post_qc", b_totals, b_lengths),
                    ("reps", a_totals, a_lengths),
                ):
                    for taxon in sorted(totals):
                        total = totals[taxon]
                        for i, pep in enumerate(spec.peptides):
                            spec_name = f"{spec.name}:{pep.name}"
                            ls = lengths[taxon][i]
                            cov_count = len(ls)
                            cov_pct = (100.0 * cov_count / total) if total else 0.0
                            _emit_tidy_rows(
                                fh, "polyprotein", rank, pool_label, taxon, total,
                                spec_name, [
                                    ("coverage_count", cov_count),
                                    ("coverage_pct", cov_pct),
                                ],
                            )
                            _emit_length_rows(
                                fh, "polyprotein", rank, pool_label, taxon, total,
                                spec_name, ls,
                            )
    return True


def write_nucleotide_taxonomic_report_tsv(
    before_seqs: list[Sequence],
    after_seqs: list[Sequence],
    cfg: dict[str, Any],
    segmented: bool,
    path: Path,
) -> bool:
    """Tidy long-format TSV companion to ``_nucleotide_taxonomic_report.txt``.

    Length metrics only (no coverage — every passing entity carries
    every required nucleotide unit by construction). ``spec`` is
    ``genome`` non-segmented; ``<segment>`` / ``total`` segmented. One
    row per (rank, pool, taxon, spec, metric).
    """
    segment_order = _resolve_segment_order(cfg, before_seqs)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(_TIDY_TSV_COLUMNS) + "\n")
        for rank in _PROTEIN_REPORT_RANKS:
            b_totals, b_lengths, headers = _nucleotide_lengths_per_taxon(
                before_seqs, rank, segmented, segment_order,
            )
            a_totals, a_lengths, _ = _nucleotide_lengths_per_taxon(
                after_seqs, rank, segmented, segment_order,
            )
            for pool_label, totals, lengths in (
                ("post_qc", b_totals, b_lengths),
                ("reps", a_totals, a_lengths),
            ):
                for taxon in sorted(totals):
                    total = totals[taxon]
                    for i, spec_name in enumerate(headers):
                        ls = lengths[taxon][i]
                        _emit_length_rows(
                            fh, "nucleotide", rank, pool_label, taxon, total,
                            spec_name, ls,
                        )
    return True


def _nucleotide_lengths_per_taxon(
    seqs: list[Sequence],
    rank: str,
    segmented: bool,
    segment_order: list[str],
) -> tuple[dict[str, int], dict[str, list[list[int]]], list[str]]:
    """For each taxon at ``rank``, gather NT lengths column-by-column.

    Non-segmented: one column ``"genome"`` with each sequence's NT length.
    Segmented: one column per segment in ``segment_order`` plus a trailing
    ``"total"`` column (sum of segment lengths per CONCAT isolate). The
    total is read from ``seq.length`` directly — it equals the sum of
    ``concat_segments`` NT lengths regardless of clustering alphabet,
    since ``build_concatenated_sequences`` always sets ``seq.sequence``
    to the joined per-segment NT.

    Items whose taxonomy lacks a value at ``rank`` are excluded (same
    rule as ``_taxonomic_report.txt`` and the protein report). Counting
    unit is the input sequence — one row per ``Sequence`` (isolate in
    segmented mode, sequence otherwise).
    """
    if segmented:
        spec_headers = list(segment_order) + ["total"]
        seg_index = {name: i for i, name in enumerate(segment_order)}
    else:
        spec_headers = ["genome"]
        seg_index = {}

    taxon_totals: dict[str, int] = {}
    taxon_lengths: dict[str, list[list[int]]] = {}
    n_cols = len(spec_headers)

    for seq in seqs:
        taxon = _seq_rank_value(seq, rank)
        if not taxon:
            continue
        taxon_totals[taxon] = taxon_totals.get(taxon, 0) + 1
        if taxon not in taxon_lengths:
            taxon_lengths[taxon] = [[] for _ in range(n_cols)]
        if segmented:
            for sub in seq.concat_segments or []:
                col = seg_index.get(sub.segment) if sub.segment else None
                if col is None:
                    continue
                length = sub.length
                if length:
                    taxon_lengths[taxon][col].append(int(length))
            # Trailing total column: the joined NT length on the CONCAT
            # itself. Equal to sum(sub.length for sub in concat_segments)
            # but cheaper and resilient against missing segment labels.
            if seq.length:
                taxon_lengths[taxon][-1].append(int(seq.length))
        else:
            if seq.length:
                taxon_lengths[taxon][0].append(int(seq.length))

    return taxon_totals, taxon_lengths, spec_headers


def _resolve_segment_order(
    cfg: dict[str, Any],
    seqs: list[Sequence],
) -> list[str]:
    """Pick the segment column order for the nucleotide report.

    Prefers the active virus block's ``segments:`` list (the canonical
    declared order — S/M/L for hanta, HA/NA/NS/… for influenza). Falls
    back to the order segments first appear across the input sequences
    so the report still works on a partially configured run.
    """
    virus_cfg = get_virus_config(cfg)
    if virus_cfg:
        seg = virus_cfg.get("segments")
        if isinstance(seg, list) and seg:
            return [str(s) for s in seg]
    seen: list[str] = []
    seen_set: set[str] = set()
    for seq in seqs:
        for sub in seq.concat_segments or []:
            if sub.segment and sub.segment not in seen_set:
                seen.append(sub.segment)
                seen_set.add(sub.segment)
    return seen


def _best_hit_on_profile(
    proteins: list[dict], profile: str,
) -> Optional[dict]:
    """Return the hit dict with the lowest dom_evalue for ``profile`` across
    every CDS in ``proteins``, or ``None`` if no CDS has a hit on it.

    The pre-computed ``passing`` flag is preserved on the returned hit;
    callers use it to decide whether the cutoff was cleared. Picking the
    single best hit is the right diagnostic level — the user wants
    "what's the BEST evidence for this marker on this isolate", not the
    cross-product of every weak match.
    """
    best: Optional[dict] = None
    for prot in proteins or []:
        for h in prot.get("hmm_hits") or []:
            if h.get("target") != profile:
                continue
            if best is None or h.get("dom_evalue", float("inf")) < best.get(
                "dom_evalue", float("inf"),
            ):
                best = h
    return best


def write_hmm_diagnostic_tsv(
    result: RunResult,
    cfg: dict[str, Any],
    complete_isolates: Optional[dict[str, list[Sequence]]],
    path: Path,
) -> bool:
    """Write a per-(representative, profile) HMM diagnostic table.

    One row for every (representative, marker_spec, declared HMM
    profile) combination — including profiles that no CDS hit at all
    (rendered as empty cells, ``hit=FALSE``). This surfaces the silent-
    failure pattern the basic QC gate hides: "65 isolates had RdRP_4
    hits at E=2e-5 but the default cutoff was 1e-5 — consider relaxing".

    Columns:
        - ``isolate_id`` / ``accession`` — the rep's identity.
        - ``segment`` — empty in non-segmented mode.
        - ``spec_name`` — the marker spec's ``name`` (or its single
          token).
        - ``profile`` — the individual HMM profile name (one of the
          spec's tokens after splitting multidomain ``--``).
        - ``hit`` — TRUE/FALSE: did any CDS in scope have a hit on this
          profile?
        - ``best_dom_evalue`` / ``best_dom_score`` / ``best_coverage`` —
          the best (lowest-E) hit's stats. Empty when ``hit=FALSE``.
        - ``ga_cutoff`` — the profile's Pfam GA threshold from the .hmm
          headers (empty when the profile has no curated GA).
        - ``cutoff_used`` — ``GA`` (when ``use_ga_when_available`` and a
          curated GA exists) or ``default_evalue=<value>``.
        - ``passing`` — TRUE/FALSE: did the best hit clear both the
          similarity gate AND the coverage gate? Empty when no hit.

    Only emits a file when the HMM tier ran this session (the
    `_hmm_runtime.active` flag) AND at least one marker spec declares
    HMMs. Returns True on write, False on skip.
    """
    from ..phylo.per_protein import (
        _hmm_tier_ran,
        _segment_proteins,
        collect_marker_specs,
    )
    from ..hmm.runner import parse_hmm_token

    if not _hmm_tier_ran(cfg, result.representatives):
        return False
    specs = collect_marker_specs(cfg)
    if not specs:
        return False

    hmm_rt = (cfg.get("_hmm_runtime") or {})
    hmm_cfg = hmm_rt.get("hmm_cfg") or (cfg.get("hmm") or {})
    ga_cutoffs = hmm_rt.get("ga_cutoffs") or {}
    use_ga = bool(hmm_cfg.get("use_ga_when_available", True))
    default_e = float(hmm_cfg.get("default_evalue", 1.0e-5))

    # For each spec, collect the SET of distinct profile names across all
    # its tokens (multidomain "A--B--C" contributes A, B, C). Token-level
    # OR / order info isn't relevant at the diagnostic level — the user
    # wants "did this profile hit, and was it below the gate" per profile.
    spec_profiles: list[tuple[str, list[str], list[str], Optional[str]]] = []
    for label, tokens, _aliases, segment in specs:
        profiles: list[str] = []
        seen: set[str] = set()
        for tok in tokens:
            try:
                for name in parse_hmm_token(tok):
                    if name and name not in seen:
                        profiles.append(name)
                        seen.add(name)
            except ValueError:
                continue
        if profiles:
            spec_profiles.append((label, tokens, profiles, segment))
    if not spec_profiles:
        return False

    # Build the rep list. Segmented runs read the per-segment proteins
    # via concat_segments; non-segmented runs use the rep's own proteins.
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(
            "isolate_id\taccession\tsegment\tspec_name\tprofile\t"
            "hit\tbest_dom_evalue\tbest_dom_score\tbest_coverage\t"
            "ga_cutoff\tcutoff_used\tpassing\n"
        )
        for rep in result.representatives:
            iso_id = rep.isolate_id or ""
            # For non-segmented reps, accession is meaningful; for CONCAT
            # reps it's None — emit the isolate_id alone.
            acc = rep.accession or (rep.id if not rep.id.startswith("CONCAT|") else "")
            for spec_label, _tokens, profiles, segment in spec_profiles:
                proteins = _segment_proteins(rep, segment) if segment else (rep.proteins or [])
                # Identify the parent segment's accession when applicable
                # (more useful than the CONCAT id).
                if segment is not None:
                    parent_acc = ""
                    for sub in rep.concat_segments or []:
                        if sub.segment == segment:
                            parent_acc = sub.accession or sub.id or ""
                            break
                    seg_cell = segment
                    acc_cell = parent_acc or acc
                else:
                    seg_cell = ""
                    acc_cell = acc
                for profile in profiles:
                    ga = ga_cutoffs.get(profile)
                    if use_ga and ga is not None:
                        cutoff_used = f"GA={ga:g}"
                    else:
                        cutoff_used = f"default_evalue={default_e:g}"
                    best = _best_hit_on_profile(proteins, profile)
                    if best is None:
                        fh.write(
                            f"{_tsv_safe(iso_id)}\t{_tsv_safe(acc_cell)}\t"
                            f"{_tsv_safe(seg_cell)}\t{_tsv_safe(spec_label)}\t"
                            f"{_tsv_safe(profile)}\t{_tsv_bool(False)}\t"
                            f"\t\t\t"
                            f"{_tsv_safe(ga) if ga is not None else ''}\t"
                            f"{_tsv_safe(cutoff_used)}\t\n"
                        )
                        continue
                    cov = coverage_of(best)
                    fh.write(
                        f"{_tsv_safe(iso_id)}\t{_tsv_safe(acc_cell)}\t"
                        f"{_tsv_safe(seg_cell)}\t{_tsv_safe(spec_label)}\t"
                        f"{_tsv_safe(profile)}\t{_tsv_bool(True)}\t"
                        f"{best.get('dom_evalue', ''):.3g}\t"
                        f"{best.get('dom_score', ''):.2f}\t"
                        f"{cov:.3f}\t"
                        f"{_tsv_safe(ga) if ga is not None else ''}\t"
                        f"{_tsv_safe(cutoff_used)}\t"
                        f"{_tsv_bool(bool(best.get('passing')))}\n"
                    )
    return True


def write_nucleotide_taxonomic_report(
    before_seqs: list[Sequence],
    after_seqs: list[Sequence],
    cfg: dict[str, Any],
    segmented: bool,
    path: Path,
    provenance: Optional[str] = None,
) -> bool:
    """Write ``{prefix}_nucleotide_taxonomic_report.txt``.

    Parallel to :func:`write_protein_taxonomic_report` but reports NT
    length distributions, not protein coverage. For each rank from
    ``subgenus`` up to ``class`` (skipping ``species`` — within-species
    length variation is dominated by sequencing/assembly noise), emits
    two length-statistics sub-tables (post-QC pool + representatives)
    with one column per segment plus a ``total`` column in segmented
    mode, or a single ``genome`` column in non-segmented mode. Cells
    are ``min, max, median, Q3-Q1, n`` in nucleotides; ``n`` is the
    number of items contributing the length.

    No coverage tables — every passing entity carries every required
    nucleotide unit by construction (segmented completeness QC
    guarantees all ``expected_segments`` are present), so coverage
    would always be 100% and add no signal.
    """
    pr_cfg = cfg.get("output", {}).get("protein_report", {}) or {}
    max_breakdown = int(pr_cfg.get("max_breakdown", 20))

    if segmented:
        segment_order = _resolve_segment_order(cfg, before_seqs)
    else:
        segment_order = []

    unit = "isolates" if segmented else "sequences"

    lines: list[str] = []
    if provenance:
        lines.append(provenance)
    lines.append("Nucleotide taxonomic report")
    lines.append(f"Generated: {datetime.date.today().isoformat()}")
    lines.append(
        f"Counting unit: {unit} "
        f"(post-QC pool fed to clustering = {len(before_seqs)} {unit}; "
        f"representatives = {len(after_seqs)} {unit})"
    )
    lines.append(
        "Length cell format: min, max, median, Q3-Q1, n (nucleotides; "
        "n = number of items contributing the length)."
    )
    if segmented:
        lines.append(
            "Per-segment columns plus a trailing `total` column "
            "(sum of segment lengths per isolate, i.e. the concatenated "
            "genome length)."
        )
    lines.append("")

    rank_printed_any = False
    for rank in _PROTEIN_REPORT_RANKS:
        b_totals, b_lengths, spec_headers = _nucleotide_lengths_per_taxon(
            before_seqs, rank, segmented, segment_order,
        )
        a_totals, a_lengths, _ = _nucleotide_lengths_per_taxon(
            after_seqs, rank, segmented, segment_order,
        )
        n_distinct = len(b_totals)
        if n_distinct == 0:
            continue
        rank_printed_any = True
        if n_distinct > max_breakdown:
            header = (
                f"{rank} ({n_distinct} distinct, top {max_breakdown} by "
                f"member count shown):"
            )
        else:
            header = f"{rank} ({n_distinct} distinct):"
        lines.append(header)
        lines.extend(_format_length_table(
            indent="  ",
            title="nucleotide length statistics [min, max, median, Q3-Q1, n] (post-QC pool)",
            totals=b_totals, lengths=b_lengths,
            spec_headers=spec_headers,
            max_breakdown=max_breakdown,
        ))
        lines.append("")
        if a_totals:
            lines.extend(_format_length_table(
                indent="  ",
                title="nucleotide length statistics [min, max, median, Q3-Q1, n] (representatives)",
                totals=a_totals, lengths=a_lengths,
                spec_headers=spec_headers,
                max_breakdown=max_breakdown,
            ))
            lines.append("")

    if not rank_printed_any:
        lines.append("(no taxonomy available at any rank from "
                     "subgenus to class)")
        lines.append("")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines).rstrip() + "\n")
    return True


def write_all_reports(
    result: RunResult,
    qc_report: QCReport,
    cfg: dict[str, Any],
    input_paths: list[str],
    output_files: list[Path],
    complete_isolates: Optional[dict[str, list[Sequence]]] = None,
    pre_clustering_sequences: Optional[list[Sequence]] = None,
) -> None:
    out_dir = Path(cfg.get("output", {}).get("dir", "./repseq_output"))
    prefix = cfg.get("output", {}).get("prefix", "repseq")

    write_run_log(result, qc_report, cfg, input_paths, output_files, out_dir / f"{prefix}_run.log")
    write_qc_tsv(qc_report, out_dir / f"{prefix}_qc_removed.tsv")
    # Force-keep audit (overrides.protect_qc): only written when something
    # was actually protected, so a normal run leaves no empty file.
    overrides_path = out_dir / f"{prefix}_overrides.tsv"
    if write_overrides_tsv(qc_report, overrides_path):
        output_files.append(overrides_path)
    # Force-select audit (overrides.force_select): only written when at
    # least one sequence was pinned.
    force_selected_path = out_dir / f"{prefix}_force_selected.tsv"
    if write_force_selected_tsv(result, force_selected_path):
        output_files.append(force_selected_path)
    # Mode-aware: segmented runs produce one row per representative
    # isolate (CONCAT|<isolate_id> Sequence with concat_segments
    # populated); non-segmented runs produce one row per representative
    # sequence. The two writers emit different column sets — each
    # tailored to the entity its rows actually represent.
    if complete_isolates:
        write_representative_isolates_tsv(
            result.representatives,
            out_dir / f"{prefix}_representative_isolates.tsv",
        )
    else:
        write_representative_sequences_tsv(
            result.representatives,
            out_dir / f"{prefix}_representative_sequences.tsv",
        )
        # Per-CDS protein tables — non-segmented counterparts of the
        # segmented _isolate_proteins.tsv / _representative_isolate_proteins.tsv
        # pair, sharing their exact column schema. The "all post-QC" file
        # needs the pre-clustering pool; the rep-only file is filtered from
        # result.representatives.
        rep_seq_ids: set[str] = {seq.id for seq in result.representatives}
        if pre_clustering_sequences is not None:
            write_sequence_proteins_tsv(
                pre_clustering_sequences,
                out_dir / f"{prefix}_sequence_proteins.tsv",
                representative_ids=rep_seq_ids,
            )
        write_sequence_proteins_tsv(
            result.representatives,
            out_dir / f"{prefix}_representative_sequence_proteins.tsv",
            representative_ids=rep_seq_ids,
        )
    write_cluster_tsv(result, out_dir / f"{prefix}_clusters.tsv")
    write_group_counts_tsv(result, out_dir / f"{prefix}_group_counts.tsv")
    # HMM diagnostic table — emitted only when the HMM tier actually ran
    # AND at least one marker spec declares HMMs. Per-(rep, profile) rows
    # let the user spot near-miss patterns ("65 isolates barely missed
    # the cutoff at E=2e-5") that the binary pass/fail gate hides.
    try:
        if write_hmm_diagnostic_tsv(
            result, cfg, complete_isolates,
            out_dir / f"{prefix}_hmm_diagnostic.tsv",
        ):
            output_files.append(out_dir / f"{prefix}_hmm_diagnostic.tsv")
    except Exception as exc:
        sys.stderr.write(f"[HMM diagnostic skipped] {exc}\n")
    if complete_isolates:
        rep_isolate_ids: set[str] = {
            seq.isolate_id for seq in result.representatives if seq.isolate_id
        }
        # Segmented reps are emitted as synthetic CONCAT|<isolate_id>;
        # mirror the parse used by write_proteins_fasta so the column is
        # correct even if a rep lacks an isolate_id attribute.
        for seq in result.representatives:
            if seq.id.startswith("CONCAT|"):
                parts = seq.id.split("|")
                if len(parts) > 1:
                    rep_isolate_ids.add(parts[1])
        write_isolate_proteins_tsv(
            complete_isolates,
            out_dir / f"{prefix}_isolate_proteins.tsv",
            representative_isolate_ids=rep_isolate_ids,
        )
        # Row-filtered companion: same schema as _isolate_proteins.tsv,
        # but only proteins of representative isolates. `representative`
        # column is retained (all values TRUE) so the file format is
        # byte-identical column-wise.
        rep_only_isolates = {
            iso_id: segs
            for iso_id, segs in complete_isolates.items()
            if iso_id in rep_isolate_ids
        }
        if rep_only_isolates:
            write_isolate_proteins_tsv(
                rep_only_isolates,
                out_dir / f"{prefix}_representative_isolate_proteins.tsv",
                representative_isolate_ids=rep_isolate_ids,
            )
    # Mode-aware filename parallel to the rep TSV split:
    #   segmented      -> _representative_isolate_proteins.fasta
    #   non-segmented  -> _representative_sequence_proteins.fasta
    # Each name mirrors its TSV companion exactly.
    proteins_basename = (
        "representative_isolate_proteins" if complete_isolates
        else "representative_sequence_proteins"
    )
    write_proteins_fasta(
        result, complete_isolates,
        out_dir / f"{prefix}_{proteins_basename}.fasta",
    )
    # Always-on per-marker unaligned FASTAs (one per declared marker
    # spec — Spike, Nucleocapsid, etc.). Headers match the
    # all-proteins file above; the difference is one file per marker
    # with only the satisfying CDS from each rep. Independent of
    # `--per-protein-phylo`: those build trees over the same CDS, this
    # always emits the raw protein set so a user can re-align or
    # post-process without re-extracting.
    try:
        pp_files = write_per_protein_fastas(
            result, cfg, complete_isolates, out_dir, prefix,
        )
        output_files.extend(pp_files)
    except Exception as exc:
        sys.stderr.write(f"[per-protein FASTA skipped] {exc}\n")
    # Always-on extra-protein FASTAs (one per `extra_protein` entry).
    # Same selection chain and headers as the per-protein FASTAs above,
    # but a separate output subdirectory so accessory / sparse proteins
    # don't visually mix with the clustering markers.
    try:
        ep_files = write_extra_protein_fastas(
            result, cfg, complete_isolates, out_dir, prefix,
        )
        output_files.extend(ep_files)
    except Exception as exc:
        sys.stderr.write(f"[extra-protein FASTA skipped] {exc}\n")
    # Always-on polyprotein cutting (v0.33.0+). For each declared
    # clustering.polyprotein / virus.polyprotein spec, slice every
    # representative's polyprotein CDS into mature peptides via one HMM
    # per peptide. Lives in its own {prefix}_polyprotein/ subdirectory
    # (peptide FASTAs + audit TSV); polyprotein still drives clustering
    # and the whole-genome tree, so this is purely additive.
    try:
        poly_files = write_polyprotein_outputs(
            result, cfg, out_dir, prefix,
        )
        output_files.extend(poly_files)
    except Exception as exc:
        sys.stderr.write(f"[polyprotein cutting skipped] {exc}\n")
