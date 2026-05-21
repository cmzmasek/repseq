"""Taxonomy-driven leaf colouring for the phyloXML output (2E / 2F).

Each external node (leaf) is given a
``<property ref="style:font_color" datatype="xsd:token"
applies_to="node">#RRGGBB</property>`` element so a viewer
(Archaeopteryx in particular) tints the leaf label by taxonomy.

Two modes, selected by the number of ranks in ``phylo.coloring.ranks``:

* **one rank** (default ``[genus]``): each distinct value of the rank
  gets its own hue, spaced around the colour wheel by the golden angle
  (137.5°) so consecutive taxa land far apart regardless of how many
  there are.
* **two ranks** (e.g. ``[genus, subgenus]``): the first (parent) rank
  fixes a base hue as above; the second (child) rank fans its values
  across a narrow hue band *around* the parent's base hue, so every
  subgenus reads as a variation of its genus colour. The band half-width
  auto-shrinks when many parents crowd the wheel, so a child band never
  bleeds into a neighbouring parent's hue.

Missing taxonomy renders as ``missing_color`` (medium grey by default).
"Missing" means empty/None or one of a small sentinel set (``unknown``,
``na``, ``?``, …), matched case-insensitively. In two-rank mode a present
parent with a *missing* child takes the parent's undifferentiated base
colour — it still has a genus, so only a missing **parent** goes grey.

The palette is keyed by taxon *name* and is built once over the full
representative set, then shared across the whole-genome tree (2E) and
every per-protein tree (2F). That keeps a given genus the same colour in
every tree of a run, which is what makes eyeballing reassortment /
topological incongruence across the per-protein trees actually possible.
"""

from __future__ import annotations

import colorsys
import re
from dataclasses import dataclass
from typing import Any, Optional

from ..models import Sequence

# The golden angle in degrees — successive multiples are maximally
# spread on the [0, 360) circle, so taxon i and taxon i+1 are always far
# apart in hue no matter the total count.
GOLDEN_ANGLE = 137.50776405003785

# Taxon strings that mean "no usable value" → rendered as missing_color.
_MISSING_TOKENS = {
    "",
    "unknown",
    "unclassified",
    "unidentified",
    "unassigned",
    "na",
    "n/a",
    "none",
    "null",
    "nan",
    "?",
    "-",
    ".",
}

DEFAULT_MISSING_COLOR = "#808080"  # medium grey
_HEX_RE = re.compile(r"^#[0-9A-Fa-f]{6}$")

# A child hue-band never spreads wider than this half-width, even when a
# parent has just one or two crowded neighbours — keeps subgenera reading
# as one genus rather than several.
_MAX_BAND_HALF_DEG = 22.0
# …and a child band is never allowed to occupy more than this fraction of
# the gap to the nearest neighbouring parent hue.
_BAND_GAP_FRACTION = 0.35


def _normalize(value: Optional[str]) -> Optional[str]:
    """Return a usable taxon string, or None for empty/sentinel values."""
    if value is None:
        return None
    s = str(value).strip()
    if s.lower() in _MISSING_TOKENS:
        return None
    return s


def _hsv_to_hex(h_deg: float, s: float, v: float) -> str:
    """HSV (hue in degrees, s/v in [0,1]) → ``#RRGGBB`` upper-case."""
    r, g, b = colorsys.hsv_to_rgb((h_deg % 360.0) / 360.0, s, v)
    return "#{:02X}{:02X}{:02X}".format(
        round(r * 255), round(g * 255), round(b * 255)
    )


def _band_half_width(parent_hues: list[float]) -> float:
    """Half-width of each parent's child hue-band, in degrees.

    Bounded by both an absolute cap and a fraction of the smallest gap
    between any two parent hues (including the wrap-around gap), so child
    bands can't overlap a neighbouring parent however many parents there
    are.
    """
    if len(parent_hues) < 2:
        return _MAX_BAND_HALF_DEG
    hues = sorted(parent_hues)
    gaps = [hues[i + 1] - hues[i] for i in range(len(hues) - 1)]
    gaps.append(360.0 - hues[-1] + hues[0])  # wrap-around
    return min(_MAX_BAND_HALF_DEG, _BAND_GAP_FRACTION * min(gaps))


def _fan_offsets(k: int, band_half: float) -> list[float]:
    """``k`` hue offsets evenly spread over ``[-band_half, +band_half]``.

    A single child sits exactly on the parent's base hue (offset 0).
    """
    if k <= 1:
        return [0.0] * max(k, 0)
    step = (2.0 * band_half) / (k - 1)
    return [-band_half + i * step for i in range(k)]


@dataclass
class ColorScheme:
    """Resolves a per-leaf ``#RRGGBB`` from one or two taxonomy ranks.

    Built by :func:`build_color_scheme` over the full representative set;
    :meth:`color_for` is then called once per leaf at serialisation time.
    """

    parent_rank: str
    child_rank: Optional[str]
    saturation: float
    value: float
    missing_color: str
    parent_hue: dict[str, float]
    child_color: dict[tuple[str, str], str]

    def _rank_value(self, seq: Sequence, rank: str) -> Optional[str]:
        if seq.taxonomy is None:
            return None
        return _normalize(seq.taxonomy.get_rank(rank))

    def color_for(self, seq: Sequence) -> str:
        """The hex colour for one leaf — always a value (grey if missing)."""
        parent = self._rank_value(seq, self.parent_rank)
        if parent is None or parent not in self.parent_hue:
            return self.missing_color
        base = self.parent_hue[parent]
        base_hex = _hsv_to_hex(base, self.saturation, self.value)
        if self.child_rank is None:
            return base_hex
        child = self._rank_value(seq, self.child_rank)
        if child is None:
            # Known parent, no child → undifferentiated parent colour.
            return base_hex
        return self.child_color.get((parent, child), base_hex)


def build_color_scheme(
    representatives: list[Sequence], cfg: dict[str, Any]
) -> Optional[ColorScheme]:
    """Build the shared leaf-colour palette for a run, or None if disabled.

    Returns None when ``phylo.coloring.enabled`` is false (so the writer
    emits no colour property at all). Otherwise returns a
    :class:`ColorScheme` even if nothing resolves — in that case every
    leaf falls through to ``missing_color``.
    """
    cc = ((cfg or {}).get("phylo", {}) or {}).get("coloring", {}) or {}
    if not cc.get("enabled", True):
        return None

    ranks = list(cc.get("ranks") or ["genus"])
    if not ranks:
        return None
    parent_rank = ranks[0]
    child_rank = ranks[1] if len(ranks) > 1 else None

    try:
        saturation = float(cc.get("saturation", 0.65))
        value = float(cc.get("value", 0.90))
    except (TypeError, ValueError):
        saturation, value = 0.65, 0.90
    missing = cc.get("missing_color") or DEFAULT_MISSING_COLOR

    parents: set[str] = set()
    children_by_parent: dict[str, set[str]] = {}
    for seq in representatives:
        if seq.taxonomy is None:
            continue
        parent = _normalize(seq.taxonomy.get_rank(parent_rank))
        if parent is None:
            continue
        parents.add(parent)
        if child_rank is not None:
            child = _normalize(seq.taxonomy.get_rank(child_rank))
            if child is not None:
                children_by_parent.setdefault(parent, set()).add(child)

    sorted_parents = sorted(parents)
    parent_hue = {
        parent: (i * GOLDEN_ANGLE) % 360.0
        for i, parent in enumerate(sorted_parents)
    }

    child_color: dict[tuple[str, str], str] = {}
    if child_rank is not None and sorted_parents:
        band_half = _band_half_width(list(parent_hue.values()))
        for parent in sorted_parents:
            kids = sorted(children_by_parent.get(parent, set()))
            base = parent_hue[parent]
            for kid, offset in zip(kids, _fan_offsets(len(kids), band_half)):
                child_color[(parent, kid)] = _hsv_to_hex(
                    base + offset, saturation, value
                )

    return ColorScheme(
        parent_rank=parent_rank,
        child_rank=child_rank,
        saturation=saturation,
        value=value,
        missing_color=missing,
        parent_hue=parent_hue,
        child_color=child_color,
    )


def is_valid_hex_color(value: Any) -> bool:
    """True for a ``#RRGGBB`` string (used by config validation)."""
    return isinstance(value, str) and bool(_HEX_RE.match(value))
