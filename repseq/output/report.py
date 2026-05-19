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


# ---------------------------------------------------------------------------
# Representative metadata TSV
# ---------------------------------------------------------------------------

def write_representative_sequences_tsv(
    representatives: list[Sequence], path: Path
) -> None:
    """Write metadata for non-segmented representative sequences to TSV.

    One row per representative sequence. Taxonomic ranks use the shared
    :data:`_TAX_RANKS` ladder so this file and
    ``_isolate_proteins.tsv`` are joinable on the same column names.
    Sub-ranks (``subgenus``, ``subfamily``, ``suborder``, ``subclass``)
    only populate when the resolver's lineage map carries them —
    commonly blank for viruses.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "accession", "organism", "description", "strain",
        "host", "collection_date", "country", "segment", "isolate_id",
        "molecule_type", "length_nt", "is_refseq", "is_reviewed",
        "ncbi_taxon_id", *_TAX_RANKS,
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
                _tsv_safe(seq.accession),
                _tsv_safe(seq.organism),
                _tsv_safe(seq.description),
                _tsv_safe(seq.strain),
                _tsv_safe(seq.host),
                _tsv_safe(seq.collection_date),
                _tsv_safe(seq.country),
                _tsv_safe(seq.segment),
                _tsv_safe(seq.isolate_id),
                _tsv_safe(seq.seq_type.value),
                _tsv_safe(seq.length),
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


# ---------------------------------------------------------------------------
# Cluster summary TSV
# ---------------------------------------------------------------------------

def write_cluster_tsv(result: RunResult, path: Path) -> None:
    """Write per-cluster summary to TSV."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("cluster_id\taccession\torganism\tcluster_size\tis_refseq\tis_reviewed\n")
        for cluster in result.clusters:
            rep = cluster.representative
            fh.write(
                f"{_tsv_safe(cluster.cluster_id)}\t{_tsv_safe(rep.accession or rep.id)}\t"
                f"{_tsv_safe(rep.organism)}\t{cluster.size}\t"
                f"{_tsv_bool(rep.is_refseq)}\t{_tsv_bool(rep.is_reviewed)}\n"
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
    write_cluster_tsv(result, out_dir / f"{prefix}_clusters.tsv")
    write_group_counts_tsv(result, out_dir / f"{prefix}_group_counts.tsv")
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
