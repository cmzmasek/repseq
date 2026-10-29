"""External per-sequence labels: ``annotate.table`` / ``--annotate-table``.

An influenza HA clade, a GenoFLU genotype, a vaccine-era bucket — metadata
that lives in a table some other tool produced (Nextclade, GenoFLU, a
curator's spreadsheet) rather than in GenBank. This step attaches those
values to the sequences and the rest of the pipeline emits them: a
``[clade=2.3.4.4b]`` tag in every protein FASTA header, a
``repseq:clade`` ``<property>`` on every tree leaf, and — with
``phylo.monophyly.extra_keys`` — per-clade monophyly rows in
``{prefix}_monophyly.tsv`` and the segment-status matrix.

**Why it is separate from ``custom`` / ``hybrid`` stratification.** Those
modes already resolve such a column, but only the mode that *stratifies* on
it ever learns the value, and a panel run is usually `global` or
`taxonomic1`. Annotation decouples "label these records" from "choose these
records", so a global 0.97 panel can carry clade tags while still spending
its budget on diversity. A ``custom --field clade`` run needs no second
table: the mode writes its own stratification labels through the same
:func:`set_labels`, so both paths produce identical headers.

**Table lookup only, deliberately.** ``custom_mode.resolve_field`` falls back
through sequence attribute → taxonomy rank → header regex before giving up;
this step does not. A ``[clade=...]`` tag has to mean "this is what the table
said", so a sequence the table does not bind is labelled ``Unknown`` and one
with disagreeing rows is ``Conflict`` — the same two sentinels the
stratification path uses, through the same binding code
(``custom_mode._table_value``), so a table keyed by HA accession reaches a
segmented isolate exactly as it does under ``--field``. That shared binding
is the load-bearing part: a CONCAT isolate has no accession of its own, and
duplicating the identifier logic here is precisely the drift trap the
``_norm_id`` invariant records.

**A table that binds nothing WARNS here; the stratification path raises.**
Keying a stratification table wrong silently collapses selection into one
``Unknown`` group, which is wrong science and must stop the run. A mis-keyed
annotation table costs a metadata tag and nothing else — and the bind rate
is only knowable after isolate assembly, hours into a long run, so aborting
there would throw away real work for a cosmetic loss.

Config is validated at STARTUP (``resolve_annotate_config``, called from
``cli._load_and_validate``), so a missing file or an unknown column fails in
seconds and under ``--dry-run``; the labels are applied after
``_handle_segmented``, once the CONCAT isolates exist.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import click

from .errors import InputError
from .models import RESERVED_LABEL_KEYS, Sequence, sanitize_label_key

#: Runtime state key. Holds only JSON-serialisable facts (the table path, the
#: resolved column names, bind counts) — never the parsed table object, which
#: the lockfile would render as ``<MetadataTable object at 0x...>`` and break
#: its byte-stability between identical runs.
ANNOTATE_RUNTIME = "_annotate_runtime"

#: Parsed tables, so startup validation and the later apply step read the
#: file once. Keyed on (path, mtime, size) rather than the path alone: the
#: same path can hold different content within one process — a rewritten
#: table between two in-process runs, or two tests using one filename — and a
#: path-only key silently served the first file's rows (2026-10-29 review).
_TABLE_CACHE: dict[tuple[str, int, int], Any] = {}


def _load_table(path: str):
    from .modes.custom_mode import load_metadata_table

    resolved = Path(path).expanduser()
    try:
        st = resolved.stat()
        key = (str(resolved), st.st_mtime_ns, st.st_size)
    except OSError:
        # Unreadable: let load_metadata_table raise the friendly error.
        return load_metadata_table(str(resolved))
    if key not in _TABLE_CACHE:
        table = load_metadata_table(str(resolved))
        # One table per path: a new stat key replaces the old entry rather
        # than accumulating parsed copies of a table that can hold a row per
        # sequence of a 1.6M-record input.
        for stale in [k for k in _TABLE_CACHE if k[0] == key[0]]:
            del _TABLE_CACHE[stale]
        _TABLE_CACHE[key] = table
    return _TABLE_CACHE[key]


def annotate_fields(cfg: dict) -> list[str]:
    """The configured ``annotate.fields``, in declaration order."""
    block = (cfg or {}).get("annotate") or {}
    raw = block.get("fields") or []
    return [str(f) for f in raw if str(f).strip()]


def resolve_annotate_config(cfg: dict) -> None:
    """Validate ``annotate.*`` and cache what the apply step needs.

    Raises :class:`InputError` for the three mistakes a user can make here —
    fields without a table, a table without fields, and a column the table
    does not have — each naming the next step. No-op when the block is empty.
    """
    block = (cfg or {}).get("annotate") or {}
    table_path = block.get("table")
    fields = annotate_fields(cfg)

    if not table_path and not fields:
        return
    if fields and not table_path:
        raise InputError(
            f"annotate.fields names {', '.join(repr(f) for f in fields)} but "
            f"annotate.table is not set, so there is no table to read them "
            f"from.\n"
            f"       Set annotate.table (or --annotate-table) to your "
            f"metadata file, or drop annotate.fields."
        )
    if table_path and not fields:
        raise InputError(
            f"annotate.table is set to {table_path} but annotate.fields is "
            f"empty, so no column would be attached.\n"
            f"       Name at least one column, e.g. --annotate-field clade."
        )

    path = Path(str(table_path)).expanduser()
    if not path.is_file():
        raise InputError(
            f"annotate.table not found: {table_path}\n"
            f"       Provide a readable TSV/CSV with one row per sequence or "
            f"isolate."
        )

    table = _load_table(str(path))
    columns: dict[str, str] = {}
    keys_seen: dict[str, str] = {}
    for field in fields:
        column = table.column_for(field)
        if column is None:
            available = ", ".join(sorted(table.columns)[:12]) or "(none)"
            raise InputError(
                f"annotate.fields names '{field}', which is not a column of "
                f"{path.name}.\n"
                f"       Columns in that file: {available}\n"
                f"       Column names are matched case-insensitively."
            )
        key = sanitize_label_key(field)
        if not key:
            raise InputError(
                f"annotate.fields names '{field}', which has no letters or "
                f"digits to form a FASTA tag / tree property name.\n"
                f"       Rename the column, or name a different one."
            )
        if key in RESERVED_LABEL_KEYS:
            raise InputError(
                f"annotate.fields names '{field}', which collides with "
                f"metadata repseq already writes (`{key}`): the label would "
                f"be silently dropped from every FASTA header and tree "
                f"rather than overwrite the real value.\n"
                f"       Rename that column in your table (e.g. "
                f"'table_{key}') and name the new column instead."
            )
        if key in keys_seen:
            raise InputError(
                f"annotate.fields names both '{keys_seen[key]}' and "
                f"'{field}', which become the same label key `{key}`, so "
                f"only one of them could ever be emitted.\n"
                f"       Name one of them, or rename the column in your "
                f"table so the two differ in more than punctuation or case."
            )
        keys_seen[key] = field
        columns[field] = column

    cfg[ANNOTATE_RUNTIME] = {
        "table": str(path),
        "fields": fields,
        "columns": columns,
    }


def set_labels(
    seq: Sequence, labels: dict[str, str], *, overwrite: bool = False,
) -> list[str]:
    """Merge ``labels`` into ``seq.group_labels`` and its segment records.

    The protein-FASTA writer receives the SEGMENT record as the parent, and a
    segmented tree leaf is the CONCAT, so both have to carry the label or one
    of the two surfaces comes out blank.

    ``overwrite=False`` (the annotate step) keeps a value already present.
    ``overwrite=True`` is for a ``custom`` / ``hybrid`` STRATIFICATION value,
    which must win: it is the label that actually chose the representatives,
    so `{prefix}_group_assignments.tsv` and the FASTA / tree labels have to
    agree with it. The annotate step runs first, so without this a run
    passing both ``--field clade`` and ``--annotate-field clade`` reported
    one table's clades in its audit file and the other's in every header
    (2026-10-29 review). Returns the fields whose differing value was
    replaced, so the caller can say so.
    """
    if not labels:
        return []
    replaced: list[str] = []
    for key, value in labels.items():
        old = (seq.group_labels or {}).get(key)
        if overwrite and old is not None and old != value:
            replaced.append(key)
    for target in (seq, *(seq.concat_segments or [])):
        current = target.group_labels
        if current is None:
            target.group_labels = dict(labels)
        elif overwrite:
            current.update(labels)
        else:
            for key, value in labels.items():
                current.setdefault(key, value)
    return replaced


def annotate_sequences(sequences: list[Sequence], cfg: dict) -> None:
    """Attach the configured table columns to every sequence, in place.

    Runs after ``_handle_segmented`` so a segmented run labels the assembled
    CONCAT isolates (and, through :func:`set_labels`, their segment records).
    Prints one bind-rate line per field and stashes the counts for the
    Methods summary. No-op unless ``annotate.table`` is configured.
    """
    runtime = (cfg or {}).get(ANNOTATE_RUNTIME) or {}
    columns = runtime.get("columns") or {}
    if not columns or not sequences:
        return

    from .modes.custom_mode import (
        CONFLICT_GROUP,
        UNKNOWN_GROUP,
        matched_rows,
        value_from_rows,
    )

    table = _load_table(runtime["table"])
    unit = (
        "isolates" if (cfg.get("segmented") or {}).get("enabled")
        else "sequences"
    )
    counts = {f: {"bound": 0, "conflict": 0} for f in columns}

    # Sequence-major, and the matched rows are resolved ONCE per sequence:
    # every configured column lives in the same rows, so a field-major loop
    # re-derived each isolate's ~34 identifiers and re-walked the table index
    # per column (2026-10-29 review). One `set_labels` call per sequence too.
    for seq in sequences:
        labels: dict[str, str] = {}
        rows = matched_rows(seq, table)
        for field, column in columns.items():
            hit = value_from_rows(rows, column, table)
            if hit is None:
                labels[field] = UNKNOWN_GROUP
            else:
                labels[field] = hit.value
                if hit.value == CONFLICT_GROUP:
                    counts[field]["conflict"] += 1
                else:
                    counts[field]["bound"] += 1
        set_labels(seq, labels)

    n = len(sequences)
    stats: dict[str, dict[str, int]] = {}
    for field in columns:
        bound = counts[field]["bound"]
        conflict = counts[field]["conflict"]
        stats[field] = {
            "n": n, "bound": bound, "conflict": conflict,
            "unknown": n - bound - conflict,
        }
        pct = 100.0 * bound / n if n else 0.0
        note = f", {conflict} conflicting" if conflict else ""
        click.echo(
            f"[annotate] {field}: {bound:,} of {n:,} {unit} labelled from "
            f"{Path(runtime['table']).name} ({pct:.1f}%{note})."
        )
        if bound == 0 and conflict == 0:
            # Advisory by design — see the module docstring. Loud, because a
            # silent run of Unknown tags looks like working output.
            click.echo(
                f"[annotate] WARNING: no {unit} matched a row of "
                f"{Path(runtime['table']).name} for '{field}', so every "
                f"label is '{UNKNOWN_GROUP}'. The table's key column "
                f"probably does not hold the identifiers in your FASTA "
                f"(accession, isolate name or strain). Selection and every "
                f"other output are unaffected.",
                err=True,
            )

    runtime["stats"] = stats


def annotate_stats(cfg: dict) -> Optional[dict]:
    """``{field: counts}`` from the last apply, for the Methods summary."""
    runtime = (cfg or {}).get(ANNOTATE_RUNTIME) or {}
    stats = runtime.get("stats")
    return dict(stats) if stats else None
