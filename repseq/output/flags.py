"""Plain-English analysis flags — a synthesis of the conflict reports.

The interesting anomalies a run can surface are scattered across several
machine-readable tables: ``{prefix}_monophyly.tsv`` (per-taxon mono/para/poly),
``{prefix}_per_protein/{prefix}_incongruence.tsv`` (tree-vs-tree RF distance),
and ``{prefix}_taxonomy_review.tsv`` (per-leaf rank conflicts, when the opt-in
review ran). This module reads whichever of those exist and distils them into
one human-readable ``{prefix}_flags.txt`` — the "what's worth a look" summary a
bench scientist can skim without joining four TSVs by hand.

It is a pure post-hoc synthesis (no new computation), runs whenever any source
table is present, and soft-fails to nothing otherwise. ``collect_flags``
returns the structured flags so the HTML report can reuse them.

Flags are heuristics over thresholds documented inline; absence of a flag is
not proof of cleanliness — the source tables remain authoritative.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

# Normalised Robinson-Foulds at/above which two marker trees are called
# incongruent enough to flag. RF is in [0, 1]; small values are common noise.
_RF_FLAG = 0.2
# Cap the incongruence flags so a many-marker run doesn't flood the report.
_MAX_INCONGRUENCE_FLAGS = 15


@dataclass
class Flag:
    """One analysis flag. ``severity`` ∈ {``warn``, ``info``}; ``category`` ∈
    {``reassortment``, ``monophyly``, ``taxonomy``}."""
    severity: str
    category: str
    message: str


def _read_tsv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    try:
        lines = path.read_text().splitlines()
    except OSError:
        return []
    if not lines:
        return []
    header = lines[0].split("\t")
    rows: list[dict] = []
    for line in lines[1:]:
        if not line.strip():
            continue
        rows.append(dict(zip(header, line.split("\t"))))
    return rows


def _short_tree(tree: str, prefix: str) -> str:
    """A friendly label for a monophyly ``tree`` relpath.

    ``coronavirus_per_protein/Spike_tree.xml`` → ``Spike``;
    ``coronavirus_tree.xml`` → ``genome``.
    """
    name = tree.rsplit("/", 1)[-1]
    if name.endswith("_tree.xml"):
        name = name[: -len("_tree.xml")]
    pfx = f"{prefix}_"
    if name.startswith(pfx):
        name = name[len(pfx):]
    return name or "genome"


def _monophyly_flags(out_dir: Path, prefix: str) -> list[Flag]:
    rows = _read_tsv(out_dir / f"{prefix}_monophyly.tsv")
    if not rows:
        return []
    genome_file = f"{prefix}_tree.xml"
    by_taxon: dict[tuple[str, str], dict[str, dict]] = {}
    for r in rows:
        by_taxon.setdefault((r["rank"], r["taxon"]), {})[r["tree"]] = r

    out: list[Flag] = []
    for (rank, taxon), trees in sorted(by_taxon.items()):
        genome = trees.get(genome_file)
        markers = {t: rr for t, rr in trees.items() if t != genome_file}

        if genome is not None and genome["status"] == "polyphyletic":
            out.append(Flag(
                "warn", "monophyly",
                f"{taxon} ({rank}) is polyphyletic on the whole-genome tree "
                f"({genome['n_leaves']} members in {genome['n_clusters']} "
                f"blocks; intruders: "
                f"{genome['intruder_taxa'] or 'unlabelled leaves'}).",
            ))
        elif genome is not None and genome["status"] == "paraphyletic":
            out.append(Flag(
                "info", "monophyly",
                f"{taxon} ({rank}) is paraphyletic on the whole-genome tree "
                f"(intruders: {genome['intruder_taxa'] or 'unlabelled leaves'}).",
            ))

        # Reassortment / recombination: clean on the genome tree but broken on
        # a marker tree — the per-marker signal.
        if genome is not None and genome["status"] == "monophyletic":
            broken = sorted(
                _short_tree(t, prefix)
                for t, rr in markers.items()
                if rr["status"] != "monophyletic"
            )
            if broken:
                out.append(Flag(
                    "warn", "reassortment",
                    f"{taxon} ({rank}) is monophyletic on the whole-genome "
                    f"tree but not on marker tree(s): {', '.join(broken)} — "
                    f"possible reassortment / recombination.",
                ))
    return out


def _incongruence_flags(out_dir: Path, prefix: str) -> list[Flag]:
    rows = _read_tsv(
        out_dir / f"{prefix}_per_protein" / f"{prefix}_incongruence.tsv"
    )
    scored: list[tuple[float, dict]] = []
    for r in rows:
        raw = (r.get("norm_rf") or "").strip()
        if raw in ("", "NA"):
            continue
        try:
            val = float(raw)
        except ValueError:
            continue
        if val >= _RF_FLAG:
            scored.append((val, r))

    out: list[Flag] = []
    for val, r in sorted(scored, key=lambda x: -x[0])[:_MAX_INCONGRUENCE_FLAGS]:
        strength = "strongly" if val >= 0.5 else "notably"
        out.append(Flag(
            "warn", "reassortment",
            f"Trees {r['tree_a']} and {r['tree_b']} are {strength} incongruent "
            f"(normalised RF = {val:.2f} over {r['n_common_taxa']} shared "
            f"taxa) — possible reassortment / recombination.",
        ))
    return out


def _taxonomy_review_flags(out_dir: Path, prefix: str) -> list[Flag]:
    rows = _read_tsv(out_dir / f"{prefix}_taxonomy_review.tsv")
    out: list[Flag] = []
    for r in rows:
        cur = (r.get("current_value") or "").strip()
        sug = (r.get("suggested_value") or "").strip()
        # A genuine conflict: a populated label the tree neighbourhood
        # disagrees with (blank-fill imputations are not flagged here).
        if cur and sug and cur != sug:
            out.append(Flag(
                "warn", "taxonomy",
                f"{r.get('accession', '?')} ({r.get('organism', '')}): "
                f"{r.get('rank', '?')} labelled '{cur}' but the tree "
                f"neighbourhood suggests '{sug}' "
                f"(support {r.get('clade_support', '?')}, "
                f"purity {r.get('clade_purity', '?')}).",
            ))
    return out


def collect_flags(out_dir: Path, prefix: str) -> list[Flag]:
    """Gather all flags from whichever source tables are present."""
    return (
        _monophyly_flags(out_dir, prefix)
        + _incongruence_flags(out_dir, prefix)
        + _taxonomy_review_flags(out_dir, prefix)
    )


_SECTIONS = [
    ("reassortment", "Reassortment / recombination signals"),
    ("monophyly", "Taxonomy ↔ tree conflicts (non-monophyletic taxa)"),
    ("taxonomy", "Per-isolate taxonomy conflicts"),
]


def write_flags_report(out_dir: Path, prefix: str) -> Optional[Path]:
    """Write ``{prefix}_flags.txt`` summarising the conflict tables.

    Returns the path, or None when no source table exists (nothing to
    synthesise). A clean run still gets a file saying so, so a user can tell
    "no conflicts" from "report not produced".
    """
    sources = (
        (out_dir / f"{prefix}_monophyly.tsv").exists()
        or (out_dir / f"{prefix}_per_protein" / f"{prefix}_incongruence.tsv").exists()
        or (out_dir / f"{prefix}_taxonomy_review.tsv").exists()
    )
    if not sources:
        return None

    flags = collect_flags(out_dir, prefix)
    path = out_dir / f"{prefix}_flags.txt"
    lines = [
        f"# {prefix} — analysis flags",
        "",
        "Taxonomy / tree conflicts worth a look, distilled from "
        "_monophyly.tsv,",
        "_incongruence.tsv, and _taxonomy_review.tsv. These are heuristics — "
        "absence of",
        "a flag is not proof of cleanliness; the source tables are "
        "authoritative.",
        "",
    ]
    if not flags:
        lines.append("No flags raised. ✓")
    else:
        for category, title in _SECTIONS:
            group = [f for f in flags if f.category == category]
            if not group:
                continue
            lines.append(f"## {title}")
            for f in group:
                marker = "⚠ " if f.severity == "warn" else "– "
                lines.append(f"{marker}{f.message}")
            lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n")
    return path
