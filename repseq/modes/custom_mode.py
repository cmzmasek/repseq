"""Custom metadata mode: group by arbitrary header fields or external metadata table."""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import click

from ..clustering import compute_diversity_curve
from ..errors import InputError
from ..annotate import set_labels
from ..models import Cluster, GroupStat, RunResult, Sequence
from ..overrides import _norm_id, _strip_version
from .base import BaseMode
from .taxonomic1 import _binary_search_threshold


#: Group label for a sequence / isolate with no value for the field.
UNKNOWN_GROUP = "Unknown"
#: Group label for an isolate whose identifiers reach table rows that carry
#: different values (e.g. its HA and NA accessions list different clades).
CONFLICT_GROUP = "Conflict"

# Header columns that say which sequence / isolate a row describes, matched
# case-insensitively. ``seqname`` is Nextclade's column, which holds the whole
# FASTA header line (only its leading identifier is used).
_KEY_COLUMNS = ("accession", "id", "isolate", "isolate_id", "strain", "seqname")

# Sequence fields offered for matching, as (kind shown in the audit, attribute).
_SEQ_IDENTIFIERS = (
    ("accession", "accession"),
    ("id", "id"),
    ("isolate", "isolate_id"),
    ("strain", "strain"),
)


def _match_key(value: Optional[str]) -> Optional[str]:
    """Matching key for an identifier: the same normalisation the overrides
    use (case-, whitespace- and pipe-insensitive), minus any accession
    version, so ``PQ123456`` in a table binds ``PQ123456.1`` in the input."""
    norm = _norm_id(value)
    if not norm:
        return None
    return _strip_version(norm) or None


# ---------------------------------------------------------------------------
# Metadata table loader
# ---------------------------------------------------------------------------

class MetadataTable:
    """A user-supplied TSV/CSV of per-sequence or per-isolate metadata.

    Every row is indexed under each of its key-column values, so a table
    can be keyed by accession, by isolate / strain name, or both.
    """

    def __init__(
        self,
        path: Path,
        columns: list[str],
        rows: list[dict[str, str]],
        key_columns: list[str],
    ) -> None:
        self.path = path
        self.columns = columns
        self.rows = rows
        self.key_columns = key_columns
        self._index: dict[str, list[int]] = {}
        for i, row in enumerate(rows):
            for col in key_columns:
                for key in _row_keys(col, row.get(col, "")):
                    hits = self._index.setdefault(key, [])
                    if not hits or hits[-1] != i:
                        hits.append(i)

    def column_for(self, field: str) -> Optional[str]:
        """The table column named ``field``, ignoring upper / lower case."""
        wanted = field.strip().lower()
        for col in self.columns:
            if col.lower() == wanted:
                return col
        return None

    def rows_for(self, raw_id: str) -> list[int]:
        key = _match_key(raw_id)
        return self._index.get(key, []) if key else []

    def example_keys(self, n: int = 3) -> list[str]:
        out: list[str] = []
        for row in self.rows:
            for col in self.key_columns:
                cell = (row.get(col) or "").strip()
                if cell:
                    out.append(cell.split()[0])
                    break
            if len(out) >= n:
                break
        return out


def _row_keys(column: str, cell: str) -> set[str]:
    cell = (cell or "").strip()
    if not cell:
        return set()
    if column.lower() == "seqname":
        # Nextclade: the full header, e.g. "PQ123456.1 Influenza A virus ..."
        # or "PQ123456.1|A/duck/...". Index the first word (what repseq uses
        # as the sequence id) and its first pipe-delimited field (the
        # accession).
        first_word = cell.split()[0]
        candidates = {first_word, first_word.split("|")[0]}
    else:
        candidates = {cell}
    return {k for k in (_match_key(c) for c in candidates) if k}


def _sniff_delimiter(header_line: str, path: Path) -> str:
    counts = {d: header_line.count(d) for d in ("\t", ";", ",")}
    best = max(counts, key=lambda d: counts[d])
    if counts[best] > 0:
        return best
    return "\t" if path.suffix.lower() in (".tsv", ".txt") else ","


def load_metadata_table(path: str | Path) -> MetadataTable:
    """Load a TSV/CSV metadata table (tab, comma or semicolon separated).

    The table needs at least one column saying which sequence or isolate a
    row describes — ``accession``, ``id``, ``isolate``, ``isolate_id``,
    ``strain`` or Nextclade's ``seqName`` (case-insensitive) — plus any
    number of metadata columns. Raises :class:`InputError` with a plain-English
    message when the file is missing, unreadable, empty, or has no key column.
    """
    path = Path(path)
    if not path.is_file():
        raise InputError(
            f"Metadata table not found: {path}\n"
            f"       Check the --metadata-table path."
        )
    try:
        # utf-8-sig drops the byte-order mark Excel writes at the start of
        # a CSV export, which would otherwise glue itself to the first
        # column name.
        with open(path, newline="", encoding="utf-8-sig") as fh:
            header_line = fh.readline()
            fh.seek(0)
            reader = csv.DictReader(fh, delimiter=_sniff_delimiter(header_line, path))
            columns = [c.strip() for c in (reader.fieldnames or []) if c]
            rows = [
                {k.strip(): (v or "").strip() for k, v in row.items()
                 if isinstance(k, str) and isinstance(v, (str, type(None)))}
                for row in reader
            ]
    except (OSError, UnicodeDecodeError, csv.Error) as exc:
        raise InputError(
            f"Could not read metadata table {path}: {exc}\n"
            f"       Save it as plain-text TSV or CSV (UTF-8) and try again."
        ) from exc

    if not columns:
        raise InputError(
            f"Metadata table {path} is empty or has no header row.\n"
            f"       The first line must name the columns, e.g. "
            f"'accession<TAB>clade'."
        )
    key_columns = [c for c in columns if c.lower() in _KEY_COLUMNS]
    if not key_columns:
        raise InputError(
            f"Metadata table {path} has no column saying which sequence each "
            f"row describes.\n"
            f"       Columns found: {', '.join(columns)}\n"
            f"       Add (or rename one to) 'accession', 'isolate', 'strain', "
            f"'id' or 'seqName' (upper / lower case does not matter)."
        )
    if not rows:
        raise InputError(
            f"Metadata table {path} has a header row but no data rows."
        )
    return MetadataTable(path, columns, rows, key_columns)


# ---------------------------------------------------------------------------
# Field lookup
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FieldValue:
    """The group label for one sequence and where it came from.

    ``matched_by`` is ``table:<kind>`` (``kind`` = accession / id / isolate /
    strain — the identifier that reached the table row), ``attribute``,
    ``taxonomy``, ``header_regex``, ``conflict`` or ``none``. ``detail`` names
    the identifier(s) that matched, for the audit file.
    """

    value: str
    matched_by: str
    detail: str = ""


def _identifiers(seq: Sequence) -> list[tuple[str, str, Optional[str]]]:
    """``(kind, raw value, segment)`` for each identifier a table row could
    be keyed on.

    A concatenated segmented isolate has no accession of its own (its id is
    synthetic), so it offers its isolate name / strain plus every
    per-segment record's identifiers — that is how a table keyed by HA
    accession reaches the isolate.
    """
    out: list[tuple[str, str, Optional[str]]] = []
    if seq.concat_segments:
        for kind, attr in (("isolate", "isolate_id"), ("strain", "strain")):
            value = getattr(seq, attr, None)
            if value:
                out.append((kind, str(value), None))
        for seg in seq.concat_segments:
            for kind, attr in _SEQ_IDENTIFIERS:
                value = getattr(seg, attr, None)
                if value:
                    out.append((kind, str(value), seg.segment or "?"))
    else:
        for kind, attr in _SEQ_IDENTIFIERS:
            value = getattr(seq, attr, None)
            if value:
                out.append((kind, str(value), None))
    return out


def matched_rows(
    seq: Sequence, table: MetadataTable,
) -> list[tuple[int, str, str]]:
    """``(row index, identifier kind, where)`` for every table row this
    sequence's identifiers reach, in identifier order, de-duplicated.

    Split out so a MULTI-column lookup derives the sequence's identifiers and
    walks the index ONCE: a segmented CONCAT offers ~34 identifier strings
    (isolate, strain, and every segment's accession/id/isolate/strain), each
    normalised through ``_norm_id`` / ``_strip_version``, so resolving two
    annotate columns per isolate used to do all of that twice for rows it had
    already found (2026-10-29 review). ``resolve_field`` and
    :func:`_table_value` go through here too, so there is one implementation.
    """
    out: list[tuple[int, str, str]] = []
    seen_rows: set[int] = set()
    for kind, raw, segment in _identifiers(seq):
        for i in table.rows_for(raw):
            if i in seen_rows:
                continue
            seen_rows.add(i)
            where = raw if segment is None else f"{raw} (segment {segment})"
            out.append((i, kind, where))
    return out


def value_from_rows(
    rows: list[tuple[int, str, str]], column: str, table: MetadataTable,
) -> Optional[FieldValue]:
    """One column's :class:`FieldValue` off rows from :func:`matched_rows`."""
    found: list[tuple[str, str, str]] = []  # (value, kind, where)
    for i, kind, where in rows:
        value = (table.rows[i].get(column) or "").strip()
        if value:
            found.append((value, kind, where))
    if not found:
        return None
    if len({v for v, _, _ in found}) == 1:
        value, kind, where = found[0]
        return FieldValue(value, f"table:{kind}", where)
    return FieldValue(
        CONFLICT_GROUP, "conflict",
        "; ".join(f"{where}={value}" for value, _, where in found),
    )


def _table_value(
    seq: Sequence, column: str, table: MetadataTable,
) -> Optional[FieldValue]:
    """One column's value for one sequence (the single-column shorthand)."""
    return value_from_rows(matched_rows(seq, table), column, table)


def _regex_group(m: re.Match) -> str:
    try:
        return m.group("value")
    except IndexError:
        return m.group(1) if m.lastindex else m.group(0)


def _regex_value(seq: Sequence, regex: str) -> Optional[FieldValue]:
    # A concatenated isolate's own header is synthetic ("CONCAT|..."), so
    # search the original per-segment headers instead.
    headers = (
        [s.header for s in seq.concat_segments]
        if seq.concat_segments else [seq.header]
    )
    values = []
    for header in headers:
        m = re.search(regex, header or "", re.IGNORECASE)
        if m and _regex_group(m):
            values.append(_regex_group(m))
    if not values:
        return None
    if len(set(values)) == 1:
        return FieldValue(values[0], "header_regex")
    return FieldValue(CONFLICT_GROUP, "conflict", "; ".join(dict.fromkeys(values)))


def resolve_field(
    seq: Sequence,
    field: str,
    metadata_table: Optional[MetadataTable],
    regex: Optional[str],
) -> FieldValue:
    """Group label for ``seq``: metadata table → sequence attribute →
    taxonomy rank → header regex → ``Unknown``."""
    # 1. External metadata table
    if metadata_table is not None:
        column = metadata_table.column_for(field)
        if column is not None:
            hit = _table_value(seq, column, metadata_table)
            if hit is not None:
                return hit

    # 2. Sequence object attribute
    attr = getattr(seq, field, None)
    if attr:
        return FieldValue(str(attr), "attribute")

    # 3. Taxonomy rank
    if seq.taxonomy:
        rank_val = seq.taxonomy.get_rank(field)
        if rank_val:
            return FieldValue(rank_val, "taxonomy")

    # 4. Regex against header
    if regex:
        hit = _regex_value(seq, regex)
        if hit is not None:
            return hit

    return FieldValue(UNKNOWN_GROUP, "none")


# ---------------------------------------------------------------------------
# Coverage reporting + audit rows (shared with hybrid mode)
# ---------------------------------------------------------------------------

def report_table_coverage(
    field: str,
    sequences: list[Sequence],
    values: list[FieldValue],
    table: MetadataTable,
    cfg: dict[str, Any],
    tag: str,
) -> dict[str, int]:
    """Print how many sequences took ``field`` from the table; raise when
    none did.

    A table that matches nothing almost always has the wrong key column (or
    describes different sequences), and silently putting the whole pool in
    one ``Unknown`` group would hide that. Returns the counts for the run
    summary.
    """
    unit = "isolates" if cfg.get("segmented", {}).get("enabled") else "sequences"
    n = len(values)
    counts = {
        "n": n,
        "table": sum(1 for v in values if v.matched_by.startswith("table:")),
        "conflict": sum(1 for v in values if v.matched_by == "conflict"),
        "unknown": sum(1 for v in values if v.matched_by == "none"),
    }
    counts["other"] = n - counts["table"] - counts["conflict"] - counts["unknown"]

    if n and counts["table"] == 0 and counts["conflict"] == 0:
        column = table.column_for(field) or field
        n_rows_hit = sum(
            1 for seq in sequences
            if any(table.rows_for(raw) for _, raw, _ in _identifiers(seq))
        )
        if n_rows_hit:
            raise InputError(
                f"Rows of {table.path.name} matched {n_rows_hit:,} of your "
                f"{n:,} {unit}, but the '{column}' column is empty for all "
                f"of them.\n"
                f"       Fill in that column, or pick another one for "
                f"--field."
            )
        # One identifier of each kind per sequence, so both an accession and
        # an isolate name are shown for a segmented isolate.
        examples = []
        for seq in sequences[:2]:
            shown_kinds: set[str] = set()
            for kind, raw, _ in _identifiers(seq):
                if kind not in shown_kinds and kind != "id":
                    shown_kinds.add(kind)
                    examples.append(f"{raw} ({kind})")
        raise InputError(
            f"No row of {table.path.name} matched any of your {n:,} {unit}, "
            f"so every one would land in the '{UNKNOWN_GROUP}' group for "
            f"'{field}'.\n"
            f"       Rows are matched on the table's "
            f"{', '.join(table.key_columns)} column(s); in segmented runs "
            f"any segment's accession reaches its isolate.\n"
            f"       Example table keys:      "
            f"{', '.join(table.example_keys()) or '(none)'}\n"
            f"       Example ids in your data: "
            f"{', '.join(examples) or '(none)'}\n"
            f"       Check that the table describes the same sequences as "
            f"your input FASTA."
        )

    pct = 100.0 * counts["table"] / n if n else 0.0
    parts = [
        f"[{tag}] {field}: {counts['table']:,} of {n:,} {unit} took their "
        f"value from {table.path.name} ({pct:.1f}%)"
    ]
    if counts["other"]:
        parts.append(
            f"{counts['other']:,} from a sequence attribute, taxonomy rank "
            f"or header regex"
        )
    if counts["unknown"]:
        parts.append(f"{counts['unknown']:,} with no value -> '{UNKNOWN_GROUP}'")
    if counts["conflict"]:
        parts.append(
            f"{counts['conflict']:,} with conflicting values -> "
            f"'{CONFLICT_GROUP}'"
        )
    click.echo("; ".join(parts))
    if counts["conflict"]:
        prefix = cfg.get("output", {}).get("prefix", "repseq")
        click.echo(
            f"WARNING: {counts['conflict']:,} {unit} reach table rows with "
            f"different '{field}' values and were grouped as "
            f"'{CONFLICT_GROUP}'. The rows with matched_by=conflict in "
            f"{prefix}_group_assignments.tsv list each one.",
            err=True,
        )
    return counts


def assignment_row(seq: Sequence, field: str, fv: FieldValue, stratum: str) -> dict[str, str]:
    """One ``{prefix}_group_assignments.tsv`` row (``representative`` is
    filled in by the writer, after force-select has run)."""
    segs = seq.concat_segments or []
    accessions = (
        ",".join(s.accession or s.id for s in segs) if segs
        else (seq.accession or seq.id)
    )
    return {
        "seq_id": seq.id,
        # Non-segmented rows leave isolate_id blank, as the rep TSVs do.
        "isolate_id": (seq.isolate_id or "") if segs else "",
        "strain": seq.strain or "",
        "accessions": accessions,
        "stratified_by": field,
        "value": fv.value,
        "matched_by": fv.matched_by,
        "detail": fv.detail,
        "stratum": stratum,
    }


# ---------------------------------------------------------------------------
# Mode
# ---------------------------------------------------------------------------

class CustomMode(BaseMode):
    """Group sequences by a user-defined field and select N representatives per group.

    Args:
        field: Attribute name, taxonomy rank, or column name in metadata table.
        n_per_group: Target representatives per group.
        metadata_table_path: Optional path to TSV/CSV with per-sequence or
            per-isolate metadata.
        field_regex: Optional regex to extract the field value from the header.
        overflow: "keep" or "trim".
    """

    def __init__(
        self,
        cfg: dict[str, Any],
        field: str,
        n_per_group: int,
        metadata_table_path: Optional[str] = None,
        field_regex: Optional[str] = None,
        overflow: str = "keep",
    ) -> None:
        super().__init__(cfg)
        self.field = field
        self.n_per_group = n_per_group
        self.field_regex = field_regex
        self.overflow = overflow
        self.metadata_table: Optional[MetadataTable] = None
        if metadata_table_path:
            self.metadata_table = load_metadata_table(metadata_table_path)
            if self.metadata_table.column_for(field) is None:
                raise InputError(
                    f"--field '{field}' is not a column of "
                    f"{self.metadata_table.path.name}.\n"
                    f"       Columns found: "
                    f"{', '.join(self.metadata_table.columns)}\n"
                    f"       Use one of those for --field (upper / lower case "
                    f"does not matter), or drop --metadata-table to group by "
                    f"a sequence attribute or taxonomy rank."
                )

    def table_label_sources(self) -> list[tuple[str, MetadataTable, str]]:
        """``(field, table, column)`` when the selection field comes from
        ``--metadata-table`` — what the label report's pre-QC column needs
        (``annotate.tally_pre_qc_labels``)."""
        if self.metadata_table is None:
            return []
        column = self.metadata_table.column_for(self.field)
        return [(self.field, self.metadata_table, column)] if column else []

    def run(self, sequences: list[Sequence]) -> RunResult:
        values = [
            resolve_field(seq, self.field, self.metadata_table, self.field_regex)
            for seq in sequences
        ]
        coverage = None
        if self.metadata_table is not None:
            coverage = report_table_coverage(
                self.field, sequences, values, self.metadata_table,
                self.cfg, "custom",
            )
        groups: dict[str, list[Sequence]] = {}
        assignments: list[dict[str, str]] = []
        relabelled: list[str] = []
        for seq, fv in zip(sequences, values):
            groups.setdefault(fv.value, []).append(seq)
            assignments.append(assignment_row(seq, self.field, fv, fv.value))
            # The stratification value is also an external label, so a
            # `--field clade` run gets [clade=...] FASTA tags and
            # repseq:clade tree properties without a second table through
            # annotate.table. A field naming metadata the output already
            # carries (`--field host`) is dropped at emission by
            # Sequence.label_items, never duplicated. overwrite=True because
            # THIS is the value that chose the representatives: if an
            # annotate table labelled the same field first, the audit file
            # and the headers would otherwise disagree.
            relabelled += set_labels(
                seq, {self.field: fv.value}, overwrite=True,
            )

        # An annotate table had labelled the same field with a different
        # value; the stratification value replaced it, because that is what
        # chose the representatives. Say so rather than letting the two
        # sources disagree silently.
        for fname in dict.fromkeys(relabelled):
            click.echo(
                f"[annotate] '{fname}' was labelled by both annotate.table "
                f"and the selection field; the selection value is the one "
                f"emitted, so the FASTA tags and tree properties match "
                "{prefix}_group_assignments.tsv.",
                err=True,
            )

        all_reps: list[Sequence] = []
        all_clusters: list[Cluster] = []
        group_stats: list[GroupStat] = []

        for group_label, group_seqs in sorted(groups.items()):
            if len(group_seqs) <= self.n_per_group:
                for seq in group_seqs:
                    all_clusters.append(
                        Cluster(cluster_id=f"{self.field}={group_label}|{seq.id}", representative=seq)
                    )
                all_reps.extend(group_seqs)
                group_stats.append(GroupStat(
                    grouping=self.field, group=str(group_label),
                    n_before=len(group_seqs), n_after=len(group_seqs),
                    clustered=False,
                ))
            else:
                search = _binary_search_threshold(
                    group_seqs, self.n_per_group, self.cfg, self.overflow,  # type: ignore[arg-type]
                    label=str(group_label),
                )
                clusters = search.clusters
                for c in clusters:
                    c.cluster_id = f"{self.field}={group_label}|{c.representative.id}"
                all_clusters.extend(clusters)
                all_reps.extend(c.representative for c in clusters)
                group_stats.append(GroupStat(
                    grouping=self.field, group=str(group_label),
                    n_before=len(group_seqs), n_after=len(clusters),
                    clustered=True, cutoff=search.threshold,
                    trimmed_from=search.trimmed_from, note=search.note,
                    cutoff_counts=compute_diversity_curve(
                        group_seqs, self.cfg, label=str(group_label)),
                ))

        snapshot: dict[str, Any] = {
            "field": self.field,
            "n_per_group": self.n_per_group,
            "overflow": self.overflow,
        }
        if self.metadata_table is not None:
            snapshot["metadata_table"] = str(self.metadata_table.path)
            snapshot["table_coverage"] = {self.field: coverage}
        return RunResult(
            mode="custom",
            representatives=all_reps,
            clusters=all_clusters,
            group_stats=group_stats,
            config_snapshot=snapshot,
            group_assignments=assignments,
        )
