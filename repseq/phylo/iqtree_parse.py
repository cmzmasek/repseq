"""Parse the ModelFinder pick(s) out of IQ-TREE's ``.iqtree`` report.

IQ-TREE writes its model-selection result into the human-readable
``.iqtree`` summary file (the one repseq copies to
``{prefix}_iqtree_summary.txt``). The model is buried in a long report,
so we extract it here for use in:

* the phyloXML ``<description>`` block (currently records the literal
  ``"MFP"`` input rather than the chosen output),
* a small ``{prefix}_iqtree_model.txt`` sidecar that bench scientists
  can grep without learning the IQ-TREE report layout,
* the ``_summary.md`` Methods section.

IQ-TREE 2.x uses two report shapes, both handled here:

**Non-partitioned (single ``-m MFP`` run):**

::

    ModelFinder
    -----------

    Best-fit model according to BIC: LG+I+G4

Some older versions emit ``Best-fit model: LG+I+G4 chosen according to BIC``;
both spellings are matched.

**Partitioned (``-p partition.nex``, MFP per charset):**

::

    List of best-fit models per partition:

      ID  Model           Speed  Parameters
       1  LG+I+G4         1.000  LG+I{0.123}+G4{0.456}
       2  JTT+G4          0.500  ...

We parse the table in ``ID`` order and zip it back to the charset names
the caller passed in (the partition path knows them from the NEXUS file
it wrote, so we don't have to scrape them out of the report).

All functions return ``{}`` on any parse problem so the caller can soft-fail
(a missing model annotation is cosmetic; failing the tree over it would be
worse than not annotating).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional


_BEST_FIT_PATTERNS = (
    # IQ-TREE 2.x canonical
    re.compile(r"^\s*Best-fit model according to BIC:\s*(\S+)\s*$", re.MULTILINE),
    # Older / alternate phrasing
    re.compile(r"^\s*Best-fit model:\s*(\S+)\s+chosen according to", re.MULTILINE),
)

_PARTITION_TABLE_HEADER = re.compile(
    r"^\s*ID\s+Model\b", re.MULTILINE,
)

_PARTITION_ROW = re.compile(
    r"^\s*(\d+)\s+(\S+)\s+",  # ID  Model  rest...
)


def _parse_non_partitioned(text: str) -> Optional[str]:
    """Return the single chosen model or ``None`` if not found."""
    for pat in _BEST_FIT_PATTERNS:
        m = pat.search(text)
        if m:
            return m.group(1).strip()
    return None


def _parse_partition_table(text: str) -> list[tuple[int, str]]:
    """Return ``[(id, model), ...]`` from the per-partition best-fit table.

    Empty list when the table header is absent or no rows parse.
    """
    m = _PARTITION_TABLE_HEADER.search(text)
    if not m:
        return []
    # Rows follow the header. Stop at the first blank line or any line
    # that doesn't start with an integer ID (sentinel lines like
    # "Best-fit partition model selected" terminate the block).
    rows: list[tuple[int, str]] = []
    for line in text[m.end():].splitlines():
        if not line.strip():
            if rows:
                break
            continue
        row_m = _PARTITION_ROW.match(line)
        if not row_m:
            if rows:
                break
            continue
        rows.append((int(row_m.group(1)), row_m.group(2).strip()))
    return rows


def parse_chosen_models(
    iqtree_report: Path,
    *,
    partition_labels: Optional[list[str]] = None,
) -> dict[str, str]:
    """Extract the ModelFinder pick(s) from an IQ-TREE ``.iqtree`` report.

    * Non-partitioned run (``partition_labels`` is ``None``): returns
      ``{"GENOME": <model>}`` or ``{}`` if the model can't be located.
    * Partitioned run (``partition_labels`` provided, in the order the
      partitions were declared to IQ-TREE): returns ``{label: model}``
      for as many rows as parse cleanly. Extra labels (no row) and
      extra rows (no label) are silently dropped — both directions are
      cosmetic and a partial annotation is better than none.
    """
    try:
        text = Path(iqtree_report).read_text(errors="replace")
    except OSError:
        return {}

    if partition_labels is None:
        model = _parse_non_partitioned(text)
        return {"GENOME": model} if model else {}

    rows = _parse_partition_table(text)
    if not rows:
        return {}
    # Zip ID-ordered rows to caller-supplied labels in order.
    out: dict[str, str] = {}
    by_index = {pid: model for pid, model in rows}
    for i, label in enumerate(partition_labels, start=1):
        model = by_index.get(i)
        if model:
            out[label] = model
    return out


def format_models_for_description(chosen: dict[str, str]) -> Optional[str]:
    """Render ``chosen`` for the phyloXML ``<description>`` block.

    * Empty → ``None`` (caller omits the field).
    * Single entry → just the model string (``"LG+I+G4"``).
    * Multiple entries → ``"label1:model1|label2:model2|..."``,
      partition order preserved.
    """
    if not chosen:
        return None
    if len(chosen) == 1:
        return next(iter(chosen.values()))
    return "|".join(f"{label}:{model}" for label, model in chosen.items())


def write_iqtree_model_file(chosen: dict[str, str], path: Path) -> None:
    """Write a flat ``<label>: <model>`` per line to ``path``.

    No-op when ``chosen`` is empty (we don't create a stub file just to
    document a parse failure).
    """
    if not chosen:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for label, model in chosen.items():
            fh.write(f"{label}: {model}\n")
