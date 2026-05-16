"""Per-leaf display label formatting for phylogenetic-tree output.

A bench scientist staring at a tree of accession numbers learns nothing;
labels need to combine taxonomy with isolate context. The formatter is
driven by a Python ``str.format``-style template (e.g.
``"{species}|{id}|{host}"``) so different downstream uses (segmented vs
non-segmented, manuscript figure vs working tree) can pick their own
shape without code changes.

Two design rules worth pinning:

1. **Strain-isolate fallback.** In segmented mode the leaf *is* the
   isolate, so the user's choice of ``{strain}`` should not produce
   ``{species}||{host}`` when the GenBank source feature has no
   ``/strain``. When a placeholder resolves to empty, we look at the
   isolate-id field; if a strain placeholder was requested and the
   isolate_id is present, we substitute the isolate_id. This means a
   single config template works for both segmented and non-segmented
   runs.

2. **Separator-drop on empty.** When a placeholder resolves to empty
   *and* ``keep_separator_on_empty=False``, we drop the preceding
   separator too so we never produce ``Hantaan virus||1976`` — the tree
   viewer would render the double-pipe verbatim. The drop is
   conservative: only the single separator character immediately before
   the empty placeholder is removed.
"""

from __future__ import annotations

import re
from typing import Optional

from ..models import Sequence


# Placeholders that *should* be present on every well-resolved leaf.
# Values for missing placeholders come from None, the empty string, or
# common "I don't know" tokens that NCBI's metadata uses ("unknown",
# "n/a", "not collected"). All such values trip separator-drop and the
# strain→isolate_id fallback.
_EMPTY_TOKENS = {"", "unknown", "n/a", "na", "none", "null", "not collected"}


def _is_empty(value: Optional[str]) -> bool:
    if value is None:
        return True
    return str(value).strip().lower() in _EMPTY_TOKENS


def _parse_year(date_str: Optional[str]) -> Optional[str]:
    """Extract a 4-digit year from a free-text collection_date.

    GenBank dates come in many shapes — ``"1976"``, ``"1976-04-12"``,
    ``"Apr-1976"``, ``"04-Apr-1976"`` — and we just want the year
    component for labels. First 4-digit run wins; returns None if none
    found.
    """
    if not date_str:
        return None
    m = re.search(r"\b(\d{4})\b", str(date_str))
    return m.group(1) if m else None


def _resolve_field(seq: Sequence, field: str) -> Optional[str]:
    """Look up a single placeholder value on a Sequence/TaxonomyInfo.

    Unknown placeholders return None (the formatter then either keeps
    the literal ``{name}`` or drops it, depending on policy — see
    ``format_leaf_label``).
    """
    f = field.lower()
    if f == "id":
        return seq.id
    if f == "accession":
        return seq.accession or seq.id
    if f == "host":
        return seq.host
    if f == "strain":
        return seq.strain
    if f == "isolate_id":
        return seq.isolate_id
    if f == "country":
        return seq.country
    if f == "date":
        return seq.collection_date
    if f == "year":
        return _parse_year(seq.collection_date)
    if f == "organism":
        return seq.organism
    if f in ("species", "genus", "family", "subfamily", "order", "suborder",
            "class", "subclass", "phylum", "kingdom", "superkingdom",
            "subgenus"):
        if seq.taxonomy is None:
            return None
        return seq.taxonomy.get_rank(f)
    return None


_PLACEHOLDER_RE = re.compile(r"\{([a-zA-Z_][a-zA-Z0-9_]*)\}")


def format_leaf_label(
    seq: Sequence,
    format_str: str,
    *,
    replace_whitespace: bool = True,
    keep_separator_on_empty: bool = False,
) -> str:
    """Render a display label for ``seq`` using ``format_str``.

    Placeholders are ``{name}`` tokens. Supported names: ``species``,
    ``genus``, ``subgenus``, ``family``, ``subfamily``, ``order``,
    ``class``, ``phylum``, ``id``, ``accession``, ``host``, ``strain``,
    ``isolate_id``, ``country``, ``date``, ``year``, ``organism``.

    Empty resolutions trip two behaviours:

    * If the placeholder is ``{strain}`` and the strain is unknown but
      ``isolate_id`` is present, the isolate_id is substituted — so a
      single template works for both segmented (isolate-keyed) and
      non-segmented (accession-keyed) input.
    * Otherwise, when ``keep_separator_on_empty`` is False, the single
      separator character immediately before the placeholder is
      removed along with the placeholder itself, so the rendered
      string never contains a ``||`` or trailing ``|``.

    With ``replace_whitespace=True`` (the default), each resolved value
    has its internal whitespace runs replaced with ``_`` — this keeps
    the label round-trippable through Newick / phyloXML viewers that
    treat whitespace as a token boundary. The separator characters of
    the format string itself are left untouched.
    """
    # Two-pass: first compute resolved values per placeholder, then
    # stitch back together so separator-drop has full context.
    matches = list(_PLACEHOLDER_RE.finditer(format_str))
    if not matches:
        return format_str

    # Pre-resolve every placeholder.
    resolved: list[tuple[re.Match[str], Optional[str]]] = []
    for m in matches:
        field = m.group(1)
        value = _resolve_field(seq, field)
        if field.lower() == "strain" and _is_empty(value):
            isolate = seq.isolate_id
            if not _is_empty(isolate):
                value = isolate
        if _is_empty(value):
            value = None
        else:
            value = str(value).strip()
            if replace_whitespace:
                value = re.sub(r"\s+", "_", value)
        resolved.append((m, value))

    # Stitch. We walk through the format string and, when we hit an
    # empty placeholder, drop one separator character on either side
    # so the rendered string never contains ``||`` or a leading/trailing
    # ``|``. The previous-separator drop covers ``{a}|{b}`` with b
    # empty; the next-separator drop (carried via ``drop_next_sep``)
    # covers the symmetric ``{a}|{b}`` with a empty.
    sep_chars = "|;,/_-:"
    out: list[str] = []
    cursor = 0
    drop_next_sep = False
    for m, value in resolved:
        chunk = format_str[cursor:m.start()]
        if drop_next_sep and chunk and chunk[0] in sep_chars:
            chunk = chunk[1:]
            drop_next_sep = False
        out.append(chunk)
        if value is None:
            if not keep_separator_on_empty:
                # Backward drop: take the preceding separator off the
                # output we've already emitted. If there is no
                # backward separator (e.g. the placeholder is at
                # column 0), arm the next-chunk drop instead — never
                # both, or "{a}|{b}|{c}" with empty b would collapse
                # to "ac" instead of "a|c".
                dropped_backward = False
                if out:
                    last = out[-1]
                    if last and last[-1] in sep_chars:
                        out[-1] = last[:-1]
                        dropped_backward = True
                if not dropped_backward:
                    drop_next_sep = True
        else:
            out.append(value)
            drop_next_sep = False
        cursor = m.end()
    tail = format_str[cursor:]
    if drop_next_sep and tail and tail[0] in sep_chars:
        tail = tail[1:]
    out.append(tail)

    return "".join(out)


def pick_format_string(cfg: Optional[dict], segmented: bool) -> str:
    """Return the label format string to use for this run.

    Defaults from ``cfg["phylo"]["labeling"]`` (``format``,
    ``segmented_format``). If ``segmented`` is True and
    ``segmented_format`` is set, that wins; otherwise ``format``;
    otherwise the hardcoded fallback ``"{species}|{id}|{host}"``.
    """
    labeling = ((cfg or {}).get("phylo", {}) or {}).get("labeling", {}) or {}
    if segmented:
        fmt = labeling.get("segmented_format") or labeling.get("format")
    else:
        fmt = labeling.get("format")
    return fmt or "{species}|{id}|{host}"


def labeling_options(cfg: Optional[dict]) -> dict:
    """Return ``replace_whitespace`` and ``keep_separator_on_empty`` for
    a run, falling back to the defaults documented in
    ``config/default_config.yaml``."""
    labeling = ((cfg or {}).get("phylo", {}) or {}).get("labeling", {}) or {}
    return {
        "replace_whitespace": labeling.get("replace_whitespace", True),
        "keep_separator_on_empty": labeling.get("keep_separator_on_empty", False),
    }
