"""Plain-English analysis flags — a synthesis of the conflict reports.

The interesting anomalies a run can surface are scattered across several
machine-readable tables: ``{prefix}_taxonomic_report.tsv`` (pre-QC vs post-QC
taxon counts — a whole genus/family eliminated by QC), ``{prefix}_monophyly.tsv``
(per-taxon mono/para/poly), ``{prefix}_per_protein/{prefix}_incongruence.tsv``
(tree-vs-tree RF distance), and ``{prefix}_taxonomy_review.tsv`` (per-leaf rank
conflicts, when the opt-in review ran). This module reads whichever of those
exist and distils them into one human-readable ``{prefix}_flags.txt`` — the
"what's worth a look" summary a bench scientist can skim without joining the
TSVs by hand.

It is a pure post-hoc synthesis (no new computation), runs whenever any source
table is present OR a QC elimination is detected, and soft-fails to nothing
otherwise. ``collect_flags`` returns the structured flags so the HTML report
can reuse them.

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
    {``qc_drop``, ``reassortment``, ``monophyly``, ``taxonomy``}."""
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


def _qc_elimination_flags(out_dir: Path, prefix: str) -> list[Flag]:
    """Genus-and-higher taxa present in the input but eliminated by QC.

    Reads the pre-QC vs post-QC ``member_count`` rows of
    ``{prefix}_taxonomic_report.tsv`` (the ``pre_qc`` rows exist only when
    metadata was resolved) and flags any taxon, at a rank in
    :data:`repseq.output.report._ALARM_RANKS` (genus and up), whose pre-QC
    count is > 0 but whose post-QC count is 0 — a whole genus / family wiped
    out by a QC gate (the silent-drop danger). The most prominent flag class,
    so it's collected first and given the leading section.
    """
    rows = _read_tsv(out_dir / f"{prefix}_taxonomic_report.tsv")
    if not rows:
        return []
    from .report import _ALARM_RANKS
    alarm = set(_ALARM_RANKS)
    counts: dict[tuple[str, str], dict[str, int]] = {}
    for r in rows:
        if r.get("report") != "diversity" or r.get("metric") != "member_count":
            continue
        rank, taxon = r.get("rank", ""), r.get("taxon", "")
        if rank not in alarm or not taxon or taxon == "*ALL*":
            continue
        try:
            val = int(float(r.get("value") or 0))
        except ValueError:
            continue
        counts.setdefault((rank, taxon), {})[r.get("pool", "")] = val

    eliminated = [
        (rank, taxon, pools.get("pre_qc", 0))
        for (rank, taxon), pools in counts.items()
        if pools.get("pre_qc", 0) > 0 and pools.get("post_qc", 0) == 0
    ]
    rank_order = {r: i for i, r in enumerate(_ALARM_RANKS)}
    eliminated.sort(key=lambda e: (rank_order.get(e[0], 99), -e[2], e[1]))
    return [
        Flag(
            "warn", "qc_drop",
            f"{taxon} ({rank}) was present in the input ({pre} pre-QC "
            f"record(s)) but has ZERO survivors after QC — the entire {rank} "
            f"was eliminated by QC. See _qc_removed.tsv for the per-record "
            f"reason (e.g. an HMM marker, a missing segment, length, or "
            f"ambiguous residues).",
        )
        for rank, taxon, pre in eliminated
    ]


def _polyprotein_wall_flags(
    out_dir: Path, prefix: str, cfg: Optional[dict] = None
) -> list[Flag]:
    """Polyprotein peptide-coverage "walls of zeros".

    Reads the ``reps`` / ``coverage_count`` rows of
    ``{prefix}_polyprotein_taxonomic_report.tsv`` (``spec`` column carries the
    ``<spec>:<peptide>`` composite) and delegates to the shared
    :func:`repseq.output.report.compute_polyprotein_walls`: a clade whose
    home polyprotein spec leaves most peptides at 0 % coverage (mistuned
    profiles), or that no spec slices at all. Thresholds come from
    ``output.polyprotein_report.wall_warning.*``; the alarm is skippable via
    its ``enabled`` flag.
    """
    ww = (
        (cfg or {}).get("output", {}).get("polyprotein_report", {}) or {}
    ).get("wall_warning", {}) or {}
    if not ww.get("enabled", True):
        return []
    rank = ww.get("rank", "genus")
    wall_fraction = ww.get("wall_fraction", 0.6)
    min_reps = ww.get("min_reps", 3)

    rows = _read_tsv(out_dir / f"{prefix}_polyprotein_taxonomic_report.tsv")
    if not rows:
        return []

    coverage: dict[tuple[str, str], dict[str, int]] = {}
    totals: dict[str, int] = {}
    spec_peptides: dict[str, list[str]] = {}
    for r in rows:
        if r.get("report") != "polyprotein" or r.get("pool") != "reps":
            continue
        if r.get("rank") != rank or r.get("metric") != "coverage_count":
            continue
        spec_full = r.get("spec", "")
        if ":" not in spec_full:
            continue
        spec_name, pep = spec_full.split(":", 1)
        taxon = r.get("taxon", "")
        if not taxon or taxon == "*ALL*":
            continue
        try:
            cnt = int(float(r.get("value") or 0))
            tot = int(float(r.get("taxon_count") or 0))
        except ValueError:
            continue
        coverage.setdefault((spec_name, taxon), {})[pep] = cnt
        totals[taxon] = tot
        peps = spec_peptides.setdefault(spec_name, [])
        if pep not in peps:
            peps.append(pep)

    if not coverage:
        return []
    from .report import compute_polyprotein_walls
    walls = compute_polyprotein_walls(
        coverage, totals, spec_peptides,
        rank=rank, wall_fraction=wall_fraction, min_reps=min_reps,
    )
    out: list[Flag] = []
    for w in walls:
        if w["kind"] == "unsliced_taxon":
            out.append(Flag(
                "warn", "polyprotein_wall",
                f"{w['taxon']} ({w['rank']}): {w['n_reps']} representative(s) "
                f"cluster but NO polyprotein spec slices any peptide for them — "
                f"the peptide profiles don't cover this clade. See "
                f"_polyprotein_taxonomic_report.txt.",
            ))
        else:
            out.append(Flag(
                "warn", "polyprotein_wall",
                f"{w['taxon']} ({w['rank']}): the '{w['spec']}' polyprotein "
                f"spec covers only {w['covered_peptides']}/{w['total_peptides']} "
                f"of its peptides for this {w['rank']} "
                f"({w['zero_peptides']} at 0 % coverage across {w['n_reps']} "
                f"representative(s)) — the peptide profiles may not fit this "
                f"clade. See _polyprotein_taxonomic_report.txt.",
            ))
    return out


def collect_flags(
    out_dir: Path, prefix: str, cfg: Optional[dict] = None
) -> list[Flag]:
    """Gather all flags from whichever source tables are present."""
    return (
        _qc_elimination_flags(out_dir, prefix)
        + _polyprotein_wall_flags(out_dir, prefix, cfg)
        + _monophyly_flags(out_dir, prefix)
        + _incongruence_flags(out_dir, prefix)
        + _taxonomy_review_flags(out_dir, prefix)
    )


_SECTIONS = [
    ("qc_drop", "Taxa eliminated entirely by QC"),
    ("polyprotein_wall", "Polyprotein peptide-coverage walls"),
    ("reassortment", "Reassortment / recombination signals"),
    ("monophyly", "Taxonomy ↔ tree conflicts (non-monophyletic taxa)"),
    ("taxonomy", "Per-isolate taxonomy conflicts"),
]


def write_flags_report(
    out_dir: Path, prefix: str, cfg: Optional[dict] = None
) -> Optional[Path]:
    """Write ``{prefix}_flags.txt`` summarising the conflict tables.

    Returns the path, or None when there is nothing to synthesise (no
    conflict table AND no QC-elimination AND no polyprotein-coverage wall).
    A clean run that *does* have a conflict table still gets a file saying
    "No flags raised", so a user can tell "no conflicts" from "report not
    produced".
    """
    flags = collect_flags(out_dir, prefix, cfg)
    standalone = [
        f for f in flags if f.category in ("qc_drop", "polyprotein_wall")
    ]
    has_conflict_tables = (
        (out_dir / f"{prefix}_monophyly.tsv").exists()
        or (out_dir / f"{prefix}_per_protein" / f"{prefix}_incongruence.tsv").exists()
        or (out_dir / f"{prefix}_taxonomy_review.tsv").exists()
    )
    # The QC-elimination and polyprotein-wall alarms fire even on a plain
    # clustering run (no --phylo, so no conflict tables) — that's exactly when
    # an inexperienced user can't otherwise see a whole genus vanish or a
    # clade's peptides silently go uncovered.
    if not has_conflict_tables and not standalone:
        return None

    path = out_dir / f"{prefix}_flags.txt"
    lines = [
        f"# {prefix} — analysis flags",
        "",
        "Anomalies worth a look, distilled from the taxonomic report "
        "(pre-QC vs",
        "post-QC), _monophyly.tsv, _incongruence.tsv, and "
        "_taxonomy_review.tsv. These",
        "are heuristics — absence of a flag is not proof of cleanliness; the "
        "source",
        "tables are authoritative.",
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
