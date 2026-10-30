"""Per-label distribution report: ``{prefix}_<label>_report.txt`` / ``.tsv``.

One report per external label — an ``annotate.fields`` column (``clade``,
``genotype``) or a ``custom`` / ``hybrid`` selection field — giving each
value's count at three stages, the label analogue of one rank of the
taxonomic report:

* **pre-QC** — rows of the label's metadata table matched by at least one
  input record before any QC removal (only when the values come from a
  table; see below);
* **pool** — the post-QC pool fed to selection (isolates when segmented);
* **reps** — the representatives.

**Why a pre-QC column needs a table, and a guard.** Labels are attached
after QC and isolate assembly (step 9b), so the input records QC removed
never carry one. The pre-QC count is therefore taken separately, at the QC
chokepoint, by matching the RAW input records against the table
(``annotate.tally_pre_qc_labels``). Its unit is table rows: for a Nextclade
table with one row per HA record that is "HA records in the input". The
records are matched by accession / id only at that point — isolate names
are read from GenBank after QC — so an isolate-name-keyed table matches
little or nothing pre-QC. Rather than print undercounts, the column is
WITHHELD whenever any value has more table-sourced post-QC members than
pre-QC table rows, which cannot happen when the pre-QC match was complete.
"Table-sourced" matters for ``custom`` / ``hybrid``: a value their
attribute / taxonomy / ``--field-regex`` fallback produced has no table row
by definition and must not be read as a failed match.

**The alarm is a direct check, not an inference.** A value with no post-QC
member is reported as **removed by QC** only when EVERY input record matched
to it is in the QC removal log (``QCReport.details``). Merely being absent
from the pool is not enough, because a value can vanish while its isolates
survive: their table rows disagree (they sit under ``Conflict``), or
assembly kept a different duplicate record of the segment, one the table
does not list (duplicates are discarded silently — they are not QC
removals). Such a value is listed as "not in the pool" instead, with no
alarm. The removal reasons of the matched records are summarised beside
each removed value, so the alarm says WHY. The alarm is returned for the
console final summary, and :mod:`repseq.output.flags` re-derives it from the
TSV (``records_matched`` vs ``records_removed_by_qc``) for ``_flags.txt``
and the HTML report.

Which labels get a report: every key :meth:`Sequence.label_items` would emit,
i.e. the same set that reaches the FASTA tags and tree properties. A
selection field named after metadata repseq already writes (``--field
host``, ``--field genus``) is reported where that metadata already is.
"""

from __future__ import annotations

import datetime
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

from ..models import Sequence, sanitize_label_key
from ..modes.custom_mode import CONFLICT_GROUP, UNKNOWN_GROUP

#: The two sentinel values the binding code assigns. Never a "value" for
#: counting distinct groups or for the alarm; listed below the real values.
_SENTINELS = (UNKNOWN_GROUP, CONFLICT_GROUP)

#: Stems already taken by other ``{prefix}_<stem>_report.*`` outputs. A label
#: whose key would land on one of them gets a ``label_`` prefix instead, so a
#: column named ``taxonomic`` can never overwrite the taxonomic report.
_TAKEN_STEMS = frozenset({
    "taxonomic", "protein_taxonomic", "nucleotide_taxonomic",
    "polyprotein_taxonomic", "subtype",
})

#: The tidy-TSV ``report`` / ``spec`` values for these files.
TIDY_REPORT = "label"
TIDY_SPEC = "_label"


def label_report_stem(key: str) -> str:
    """The ``<stem>`` in ``{prefix}_<stem>_report.txt`` for label ``key``.

    A key that is another report's stem gets a ``label_`` prefix, and so does
    any key already starting with ``label_`` — otherwise ``taxonomic`` and a
    second label literally named ``label_taxonomic`` would share one file.
    With both rules the mapping is one-to-one.
    """
    if key in _TAKEN_STEMS or key.startswith("label_"):
        return f"label_{key}"
    return key


def is_label_report_stem(stem: str) -> bool:
    """True unless ``stem`` belongs to one of the other report files."""
    return stem not in _TAKEN_STEMS


def table_sourced_counts(group_assignments: Optional[Iterable[dict]]) -> dict[str, Counter]:
    """``{label key: Counter(value -> pool members)}`` counting only the
    members whose ``custom`` / ``hybrid`` value came from the metadata table.

    Every selection field gets an entry (possibly empty), so the pre-QC guard
    never falls back to the full pool count for it — that count includes
    values the attribute / taxonomy / regex fallback produced.
    """
    out: dict[str, Counter] = {}
    for row in group_assignments or ():
        key = sanitize_label_key(row.get("stratified_by", ""))
        if not key:
            continue
        counts = out.setdefault(key, Counter())
        if str(row.get("matched_by", "")).startswith("table:"):
            counts[row.get("value", "")] += 1
    return out


_REASON_WORD = re.compile(r"^[a-z][a-z_]{3,}$")


def short_reason(reason: str) -> str:
    """A QC removal reason without its per-record detail.

    ``marker_length:4:481aa_below_485`` → ``marker_length``;
    ``incomplete_isolate:missing_segments:1,2`` →
    ``incomplete_isolate:missing_segments``; a ``_sibling`` reason (the
    record was fine, another segment of its isolate failed) folds into the
    stage it names, since the stage is what removed the isolate.
    """
    parts = str(reason or "").split(":")
    head = re.split(r"[=<>]", parts[0], maxsplit=1)[0]
    if head.endswith("_sibling"):
        head = head[: -len("_sibling")]
    if len(parts) > 1 and _REASON_WORD.match(parts[1]):
        head += ":" + parts[1]
    return head or "unknown"


@dataclass
class LabelTable:
    """Per-value counts for one label at the three stages."""

    key: str
    source: Optional[str]           # table file name, when known
    pre_qc: Optional[Counter]       # None when unavailable or withheld
    pre_qc_note: Optional[str]      # why it is missing, when it is
    pool: Counter
    reps: Counter
    n_pool: int
    n_reps: int
    #: value -> (records matched pre-QC, records in the QC removal log), for
    #: values with no pool member.
    record_check: dict[str, tuple[int, int]] = field(default_factory=dict)
    #: value -> Counter(short removal reason) for values removed by QC.
    removal_reasons: dict[str, Counter] = field(default_factory=dict)

    def values(self) -> list[str]:
        """Real values (no sentinels), largest pool first, then pre-QC,
        then a natural name order so ``2.3.4.4`` sits before ``2.3.4.4b``."""
        names = set(self.pool) | set(self.reps) | set(self.pre_qc or ())
        names -= set(_SENTINELS)
        pre = self.pre_qc or Counter()
        return sorted(
            names,
            key=lambda v: (-self.pool[v], -pre[v], _natural_key(v)),
        )

    def sentinels(self) -> list[str]:
        return [s for s in _SENTINELS if self.pool[s] or self.reps[s]]

    def absent(self) -> list[str]:
        """Values with pre-QC rows but no pool member (shown column only)."""
        if self.pre_qc is None:
            return []
        return [
            v for v in self.values()
            if self.pre_qc[v] > 0 and self.pool[v] == 0
        ]

    def removed_by_qc(self, value: str) -> bool:
        matched, removed = self.record_check.get(value, (0, 0))
        return matched > 0 and removed == matched

    def eliminated(self) -> list[tuple[str, int]]:
        """``(value, pre-QC rows)`` for values QC removed in full: absent
        from the pool AND every matched input record in the removal log."""
        return [
            (v, self.pre_qc[v]) for v in self.absent() if self.removed_by_qc(v)
        ]


def _natural_key(value: str) -> list:
    return [
        (0, int(tok), "") if tok.isdigit() else (1, 0, tok.lower())
        for tok in re.split(r"(\d+)", value) if tok
    ]


def _labels(seq: Sequence) -> dict[str, str]:
    return dict(seq.label_items())


def _removal_log(candidates: set[str], qc_details) -> dict[str, str]:
    """``{record id: reason}`` for the candidate ids found in the QC log.

    One pass over the log, comparing raw strings: the log stores
    ``seq.id`` (or ``seq.accession``) and the pre-QC tally stored the same
    attributes of the same objects, so no normalisation is needed — and only
    candidates are kept, so a 1.4M-row log costs no memory.
    """
    out: dict[str, str] = {}
    if not candidates:
        return out
    for entry in qc_details or ():
        rid = entry.get("id")
        if rid in candidates and rid not in out:
            out[rid] = entry.get("reason", "")
    return out


def build_label_tables(
    pool: list[Sequence],
    reps: list[Sequence],
    pre_qc_labels: Optional[dict] = None,
    *,
    table_counts: Optional[dict[str, Counter]] = None,
    qc_details=None,
) -> list[LabelTable]:
    """One :class:`LabelTable` per label key, in declaration order.

    ``pre_qc_labels`` is ``QCReport.pre_qc_labels``:
    ``{key: {"source", "counts": {value: rows}, "records": {value:
    [(id, accession), ...]}}}``. ``table_counts`` is
    :func:`table_sourced_counts` (custom / hybrid only), ``qc_details`` the
    QC removal log. A key known only pre-QC (a total wipe-out leaves no
    labelled pool) still gets a table, so the alarm can fire on exactly the
    run where it matters most.
    """
    pre_qc_labels = pre_qc_labels or {}
    table_counts = table_counts or {}
    pool_labels = [_labels(s) for s in pool]
    rep_labels = [_labels(s) for s in reps]

    keys: list[str] = []
    for labels in pool_labels:
        for key in labels:
            if key not in keys:
                keys.append(key)
    for key in pre_qc_labels:
        if key not in keys:
            keys.append(key)

    tables: list[LabelTable] = []
    for key in keys:
        # A record without this key counts as Unknown: annotate and the
        # selection modes label every pool record, so this only guards odd
        # inputs, and an Unknown row is more honest than a silent gap.
        pool_c = Counter(lb.get(key, UNKNOWN_GROUP) for lb in pool_labels)
        reps_c = Counter(lb.get(key, UNKNOWN_GROUP) for lb in rep_labels)
        info = pre_qc_labels.get(key) or {}
        pre: Optional[Counter] = None
        note: Optional[str] = None
        if "counts" in info:
            pre = Counter({
                v: int(n) for v, n in (info.get("counts") or {}).items()
                if v not in _SENTINELS
            })
            sourced = table_counts.get(key, pool_c)
            short = sorted(
                v for v, n in sourced.items()
                if v not in _SENTINELS and n > pre[v]
            )
            if short:
                note = (
                    f"pre-QC column withheld: {len(short)} value(s) have more "
                    f"post-QC members matched through the table than table "
                    f"rows matched before QC (e.g. {short[0]}: "
                    f"{pre[short[0]]} pre-QC vs {sourced[short[0]]} "
                    f"post-QC), so the pre-QC match was incomplete. Input "
                    f"records are matched by accession / id before QC; a "
                    f"table keyed by isolate name or strain can only be "
                    f"matched after GenBank metadata is read."
                )
                pre = None
        else:
            note = (
                "no pre-QC column: these values do not come from a metadata "
                "table, so the input records removed by QC cannot be counted."
            )
        tables.append(LabelTable(
            key=key, source=info.get("source"), pre_qc=pre, pre_qc_note=note,
            pool=pool_c, reps=reps_c, n_pool=len(pool), n_reps=len(reps),
        ))

    # The direct check: are ALL the input records matched to an absent value
    # in the QC removal log? One pass over the log for every label at once.
    wanted: dict[tuple[str, str], list] = {}
    for t in tables:
        records = (pre_qc_labels.get(t.key) or {}).get("records") or {}
        for v in t.absent():
            wanted[(t.key, v)] = list(records.get(v) or ())
    candidates = {
        rid for recs in wanted.values() for rec in recs for rid in rec if rid
    }
    log = _removal_log(candidates, qc_details)
    by_key = {t.key: t for t in tables}
    for (key, value), recs in wanted.items():
        reasons: Counter = Counter()
        removed = 0
        for rec in recs:
            hit = next((log[rid] for rid in rec if rid and rid in log), None)
            if hit is not None:
                removed += 1
                reasons[short_reason(hit)] += 1
        t = by_key[key]
        t.record_check[value] = (len(recs), removed)
        if recs and removed == len(recs):
            t.removal_reasons[value] = reasons
    return tables


def _pct(a: int, b: int) -> str:
    return f"{(100.0 * a / b):.1f}%" if b else "n/a"


def _reasons_text(reasons: Counter, n: int = 2) -> str:
    top = reasons.most_common(n)
    text = ", ".join(f"{r} {c}" for r, c in top)
    if len(reasons) > n:
        text += ", ..."
    return text


def render_label_report(
    table: LabelTable, segmented: bool, provenance: Optional[str] = None,
) -> str:
    """The plain-text report for one label."""
    unit = "isolates" if segmented else "sequences"
    key = table.key
    pre = table.pre_qc
    values = table.values()
    sentinels = table.sentinels()

    lines: list[str] = []
    if provenance:
        lines.append(provenance)
    lines.append(f"Label report: {key}")
    lines.append(f"Generated: {datetime.date.today().isoformat()}")
    if table.source:
        lines.append(f"Source: {table.source}")
    lines.append(
        f"Counting unit: {unit} (pool = post-QC pool fed to selection, "
        f"{table.n_pool:,}; reps = representatives, {table.n_reps:,})"
    )
    if pre is not None:
        lines.append(
            "pre-QC = rows of the table matched by at least one input record "
            "before any QC removal. Unit: table rows (for a table with one row "
            "per HA record, that is HA records), so compare presence across "
            "the pre-QC / pool boundary rather than exact magnitudes."
        )
    elif table.pre_qc_note:
        lines.append(f"NOTE: {table.pre_qc_note}")

    labelled_pool = table.n_pool - sum(table.pool[s] for s in _SENTINELS)
    labelled_reps = table.n_reps - sum(table.reps[s] for s in _SENTINELS)
    lines.append(
        f"Coverage: {labelled_pool:,} / {table.n_pool:,} pool {unit} carry a "
        f"{key} ({_pct(labelled_pool, table.n_pool)}); {labelled_reps:,} / "
        f"{table.n_reps:,} representatives "
        f"({_pct(labelled_reps, table.n_reps)})"
    )
    lines.append("")

    n_pre = sum(1 for v in values if pre is not None and pre[v] > 0)
    n_pool = sum(1 for v in values if table.pool[v] > 0)
    n_reps = sum(1 for v in values if table.reps[v] > 0)
    lines.append(f"== Distinct {key} values ==")
    lines.append("")
    distinct = f"  pool: {n_pool}    reps: {n_reps}"
    if pre is not None:
        distinct = f"  pre-QC: {n_pre}  " + distinct
    lines.append(distinct)
    lines.append("")

    lines.append(
        f"== {key} distribution (sorted by pool count; all values shown) =="
    )
    lines.append("")
    name_w = max([len(key)] + [len(v) for v in values + sentinels])
    cols = (["pre-QC"] if pre is not None else []) + ["pool", "reps"]
    header = f"  {key:<{name_w}}" + "".join(f"  {c:>8}" for c in cols)
    rule = "  " + "-" * (len(header) - 2)
    lines.append(header)
    lines.append(rule)

    def _row(name: str, pre_cell: str, note: str = "") -> str:
        cells = ([pre_cell] if pre is not None else []) + [
            str(table.pool[name]), str(table.reps[name]),
        ]
        row = f"  {name:<{name_w}}" + "".join(f"  {c:>8}" for c in cells)
        return row + (f"   {note}" if note else "")

    absent = set(table.absent())
    n_not_in_pool = 0
    for v in values:
        note = ""
        if v in absent:
            if table.removed_by_qc(v):
                note = (
                    f"removed by QC "
                    f"({_reasons_text(table.removal_reasons.get(v, Counter()))})"
                )
            else:
                note = "not in pool, not removed by QC (see note below)"
                n_not_in_pool += 1
        lines.append(_row(v, str(pre[v]) if pre is not None else "", note))
    if sentinels:
        lines.append(rule)
        for s in sentinels:
            note = (
                "no table row matched" if s == UNKNOWN_GROUP
                else "table rows disagreed"
            )
            lines.append(_row(s, "-", note))
    lines.append("")

    elim = table.eliminated()
    if elim:
        lines.append(
            f"WARNING: {len(elim)} {key} value(s) were present in the input "
            f"but every input record matched to them was removed by QC "
            f"(\"removed by QC\" above, with the most common reasons). See "
            f"_qc_removed.tsv for the per-record reason."
        )
        lines.append("")
    if n_not_in_pool:
        lines.append(
            f"NOTE: {n_not_in_pool} {key} value(s) are absent from the pool "
            f"although not all of their input records were removed by QC: "
            f"their {unit} survived but carry another label — the Conflict "
            f"row when their table rows disagree, or Unknown / another value "
            f"when assembly kept a different duplicate record of the segment "
            f"than the one the table lists."
        )
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def write_label_report_tsv(table: LabelTable, path: Path) -> None:
    """Tidy long-format companion, sharing the 8-column schema of the other
    report TSVs (``report`` = ``label``, ``rank`` = the label key, ``spec`` =
    ``_label``). ``pre_qc`` rows appear only when the pre-QC column is shown,
    and each value absent from the pool also gets ``records_matched`` /
    ``records_removed_by_qc`` on its ``pre_qc`` row — equal counts are what
    the ``_flags.txt`` alarm keys on."""
    from .report import _TIDY_TSV_COLUMNS, _emit_tidy_rows

    values = table.values()
    sentinels = table.sentinels()
    pre = table.pre_qc
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(_TIDY_TSV_COLUMNS) + "\n")
        pools: list[tuple[str, Counter]] = []
        if pre is not None:
            pools.append(("pre_qc", pre))
        pools += [("post_qc", table.pool), ("reps", table.reps)]
        for pool_name, counts in pools:
            n = sum(1 for v in values if counts[v] > 0)
            _emit_tidy_rows(
                fh, TIDY_REPORT, table.key, pool_name, "*ALL*", n, TIDY_SPEC,
                [("distinct_taxa", n)],
            )
        for v in values:
            for pool_name, counts in pools:
                metrics: list[tuple[str, int]] = [("member_count", counts[v])]
                if pool_name == "pre_qc" and v in table.record_check:
                    matched, removed = table.record_check[v]
                    metrics += [
                        ("records_matched", matched),
                        ("records_removed_by_qc", removed),
                    ]
                _emit_tidy_rows(
                    fh, TIDY_REPORT, table.key, pool_name, v, counts[v],
                    TIDY_SPEC, metrics,
                )
        for s in sentinels:
            for pool_name, counts in pools[-2:]:
                _emit_tidy_rows(
                    fh, TIDY_REPORT, table.key, pool_name, s, counts[s],
                    TIDY_SPEC, [("member_count", counts[s])],
                )


def write_label_reports(
    pool: list[Sequence],
    reps: list[Sequence],
    pre_qc_labels: Optional[dict],
    *,
    out_dir: Path,
    prefix: str,
    segmented: bool,
    provenance: Optional[str] = None,
    table_counts: Optional[dict[str, Counter]] = None,
    qc_details=None,
) -> tuple[list[Path], list[dict], dict[str, dict]]:
    """Write one ``.txt`` + ``.tsv`` per label.

    Returns ``(paths written, eliminated, status)``: ``eliminated`` holds
    ``{label, value, pre_qc_count, reasons}`` per value QC removed in full,
    for the console final summary (``reasons``: most common short removal
    reasons, most frequent first); ``status`` is ``{key: {"file": <txt
    name>, "pre_qc": shown?}}``, for the Methods prose, which has to say
    whether the pre-QC column could be given. Writes nothing for a run with
    no labels.
    """
    tables = build_label_tables(
        pool, reps, pre_qc_labels,
        table_counts=table_counts, qc_details=qc_details,
    )
    paths: list[Path] = []
    eliminated: list[dict] = []
    status: dict[str, dict] = {}
    out_dir = Path(out_dir)
    for table in tables:
        stem = f"{prefix}_{label_report_stem(table.key)}_report"
        txt = out_dir / f"{stem}.txt"
        tsv = out_dir / f"{stem}.tsv"
        txt.parent.mkdir(parents=True, exist_ok=True)
        txt.write_text(render_label_report(table, segmented, provenance))
        write_label_report_tsv(table, tsv)
        paths += [txt, tsv]
        eliminated += [
            {
                "label": table.key, "value": v, "pre_qc_count": n,
                "reasons": [
                    r for r, _ in
                    table.removal_reasons.get(v, Counter()).most_common(3)
                ],
            }
            for v, n in table.eliminated()
        ]
        status[table.key] = {
            "file": txt.name, "pre_qc": table.pre_qc is not None,
        }
    return paths, eliminated, status
