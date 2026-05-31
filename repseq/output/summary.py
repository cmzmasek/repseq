"""Auto-generated Methods-section starter for each repseq run.

Writes ``{prefix}_summary.md`` describing what the pipeline did,
in scientific prose suitable as a starting point for a methods
section. Numbers come from the live ``QCReport`` and ``RunResult``;
tool versions are detected by invoking each external binary's
``--version`` (or equivalent) once at write time.
"""

from __future__ import annotations

import datetime as _dt
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Optional

from .. import __version__ as REPSEQ_VERSION
from ..models import QCReport, RunResult


# ---------------------------------------------------------------------------
# Tool-version detection
# ---------------------------------------------------------------------------

# Tools we may invoke and how to coax a version string out of them. Each
# entry: (binary_name, argv_after_binary, regex_to_extract). The regex must
# have group(1) = version string, or None to use the first matching line.
_TOOL_PROBES: tuple[tuple[str, list[str], Optional[str]], ...] = (
    # cd-hit prints "CD-HIT version 4.8.1 (built on ...)" to stdout on -h.
    ("cd-hit",     ["-h"],          r"CD-HIT version\s+(\S+)"),
    ("cd-hit-est", ["-h"],          r"CD-HIT version\s+(\S+)"),
    # mmseqs version prints "16.747c6"-like single line to stdout.
    ("mmseqs",     ["version"],     r"^(\S+)$"),
    # mafft --version writes "v7.520 (2023/Mar/16)" to stderr.
    ("mafft",      ["--version"],   r"v?(\d[\d.]*)"),
    # trimal --version prints "trimAl v1.4.rev15 build[...]" to stdout.
    ("trimal",     ["--version"],   r"trimAl\s+v?(\S+)"),
    # FastTree with no args prints usage with "FastTree 2.1.11" on stderr.
    ("FastTree",   [],              r"FastTree\s+(?:Version\s+)?(\d[\d.]*)"),
    # iqtree2 --version: "IQ-TREE multicore version 2.2.5 ..."
    ("iqtree2",    ["--version"],   r"IQ-TREE.+version\s+(\S+)"),
    # hmmscan -h prints "# HMMER 3.3.2 (Nov 2020); http://hmmer.org/" on stdout.
    ("hmmscan",    ["-h"],          r"HMMER\s+(\S+)"),
)


def _probe_version(binary: str, args: list[str], pattern: Optional[str]) -> Optional[str]:
    """Best-effort: invoke the binary, grep its output for a version string.

    Returns None if the binary is not on PATH or invocation fails. Both
    stdout and stderr are scanned because different tools write to
    different streams (mafft → stderr, cd-hit → stdout, FastTree → stderr).
    """
    path = shutil.which(binary)
    if not path:
        return None
    try:
        # FastTree with no args exits non-zero (usage); accept that.
        proc = subprocess.run(
            [path, *args],
            capture_output=True,
            text=True,
            timeout=5,
            check=False,
        )
    except (subprocess.TimeoutExpired, OSError):
        return None
    haystack = (proc.stdout or "") + "\n" + (proc.stderr or "")
    if pattern is None:
        return haystack.splitlines()[0].strip() if haystack.strip() else None
    m = re.search(pattern, haystack, re.MULTILINE)
    return m.group(1) if m else None


_VERSION_CACHE: dict[str, Optional[str]] = {}


def detect_tool_versions() -> dict[str, Optional[str]]:
    """Detect installed versions of every external tool repseq may use.

    Results cached for the lifetime of the process so a single run
    pays for the probes only once. Missing binaries return None.
    """
    if _VERSION_CACHE:
        return dict(_VERSION_CACHE)
    for binary, args, pattern in _TOOL_PROBES:
        _VERSION_CACHE[binary] = _probe_version(binary, args, pattern)
    return dict(_VERSION_CACHE)


def _python_version() -> str:
    v = sys.version_info
    return f"{v.major}.{v.minor}.{v.micro}"


def _biopython_version() -> Optional[str]:
    try:
        import Bio
        return getattr(Bio, "__version__", None)
    except ImportError:
        return None


# ---------------------------------------------------------------------------
# Citations (hardcoded — these don't change)
# ---------------------------------------------------------------------------

_CITATIONS: dict[str, str] = {
    "repseq":      "https://github.com/cmzmasek/repseq",
    "python":      "https://www.python.org",
    "biopython":   "Cock et al. 2009, *Bioinformatics* 25(11):1422-1423",
    "cd-hit":      "Fu, Niu, Zhu, Wu & Li 2012, *Bioinformatics* 28(23):3150-3152",
    "cd-hit-est":  "Fu, Niu, Zhu, Wu & Li 2012, *Bioinformatics* 28(23):3150-3152",
    "mmseqs":      "Steinegger & Söding 2017, *Nat. Biotechnol.* 35(11):1026-1028",
    "mafft":       "Katoh & Standley 2013, *Mol. Biol. Evol.* 30(4):772-780",
    "trimal":      "Capella-Gutiérrez, Silla-Martínez & Gabaldón 2009, *Bioinformatics* 25(15):1972-1973",
    "FastTree":    "Price, Dehal & Arkin 2010, *PLoS ONE* 5(3):e9490",
    "iqtree2":     "Minh et al. 2020, *Mol. Biol. Evol.* 37(5):1530-1534",
    "ncbi":        "Sayers et al. 2022, *Nucleic Acids Res.* 50(D1):D20-D26",
    "uniprot":     "UniProt Consortium 2023, *Nucleic Acids Res.* 51(D1):D523-D531",
    "hmmer":       "Eddy 2011, *PLoS Comput. Biol.* 7(10):e1002195",
    "pfam":        "Mistry et al. 2021, *Nucleic Acids Res.* 49(D1):D412-D419",
}


# ---------------------------------------------------------------------------
# Section renderers
# ---------------------------------------------------------------------------

def _fmt_int(n: Optional[int]) -> str:
    return f"{n:,}" if isinstance(n, int) else "—"


def _render_header(cfg: dict, command: str, mode: str) -> str:
    when = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    virus = cfg.get("segmented", {}).get("virus") or "input dataset"
    title = f"Methods — repseq {mode} selection ({virus})"
    cmd_block = f"\n*Run command:* `{command}`" if command else ""
    return (
        f"# {title}\n\n"
        f"*Generated by repseq v{REPSEQ_VERSION} on {when}.*"
        f"{cmd_block}\n"
    )


def _render_input(qc_report: QCReport, input_paths: list[str]) -> str:
    n = _fmt_int(qc_report.total_input)
    n_files = len(input_paths)
    file_phrase = (
        f"a single FASTA file" if n_files == 1 else f"{n_files} FASTA files"
    )
    return (
        f"## Input\n\n"
        f"A total of **{n}** sequences were read from {file_phrase} and "
        f"used as input. Taxonomy, GenBank source-feature qualifiers "
        f"(`/isolate`, `/strain`, `/segment`, `/host`, `/country`, "
        f"`/collection_date`), and CDS annotations were retrieved from "
        f"NCBI Entrez and UniProt REST where applicable and cached "
        f"locally in a SQLite database with TTL eviction "
        f"({_CITATIONS['ncbi']}; {_CITATIONS['uniprot']}).\n"
    )


def _render_qc(qc_report: QCReport, cfg: dict) -> str:
    qc_cfg = cfg.get("qc", {}) or {}
    parts: list[str] = []
    if qc_report.removed_duplicates:
        parts.append(f"**{_fmt_int(qc_report.removed_duplicates)}** exact duplicates")
    # Describe the whole-genome length filter only when it actually ran
    # (non-segmented + enabled). In segmented mode it's skipped
    # (length_filter_skipped) and removed_length is overloaded by the
    # per-segment filter in segmented/completeness.py; that case is
    # described separately below.
    if qc_report.removed_length and not qc_report.length_filter_skipped:
        glf = qc_cfg.get("genome_length_filter", {}) or {}
        mn = glf.get("min")
        mx = glf.get("max")
        bound_bits: list[str] = []
        if mn is not None:
            bound_bits.append(f"shorter than {_fmt_int(mn)} nt")
        if mx is not None:
            bound_bits.append(f"longer than {_fmt_int(mx)} nt")
        length_note = f" ({' or '.join(bound_bits)})" if bound_bits else ""
        parts.append(
            f"**{_fmt_int(qc_report.removed_length)}** outside the configured "
            f"whole-genome length bounds{length_note}"
        )
    if qc_report.removed_ambiguous:
        thr = qc_cfg.get("ambiguous_threshold", 0.05)
        parts.append(
            f"**{_fmt_int(qc_report.removed_ambiguous)}** with > {thr:.0%} "
            f"ambiguous characters"
        )
    basic = (
        f"Initial quality control removed "
        f"{', '.join(parts) if parts else 'no'} sequence(s) "
        f"(*{_fmt_int(qc_report.passed)}* passed basic QC)."
    )

    extra_lines: list[str] = []
    # The protein-annotation check is sequence-level QC but is segmented-aware:
    # in segmented mode it compares each record's GenBank CDS count against the
    # per-segment expected count (and so depends on seq.segment being populated
    # upstream by _populate_genbank_isolate_segment). It does NOT run during
    # the segmented-handling stage despite needing segment labels — surfacing
    # it as its own sentence here (rather than in the "Initial QC" list) keeps
    # the units honest: per NCBI segment record in segmented mode, per sequence
    # otherwise.
    if qc_report.removed_proteins:
        seg_cfg = cfg.get("segmented", {}) or {}
        virus_name = seg_cfg.get("virus") or ""
        virus_cfg = (seg_cfg.get("viruses") or {}).get(virus_name, {}) if virus_name else {}
        has_per_seg = bool(virus_cfg.get("expected_proteins_per_segment"))
        pa_cfg = (qc_cfg.get("protein_annotation") or {})
        min_proteins = pa_cfg.get("min_proteins") if pa_cfg.get("enabled") else None
        if seg_cfg.get("enabled") and has_per_seg:
            cfg_path = (
                f"segmented.viruses.{virus_name}.expected_proteins_per_segment"
                if virus_name
                else "segmented.viruses.*.expected_proteins_per_segment"
            )
            extra_lines.append(
                f"**{_fmt_int(qc_report.removed_proteins)}** segment record(s) "
                f"were dropped by the protein-annotation check — their GenBank "
                f"CDS count did not match `{cfg_path}` for their segment. "
                f"Failing records appear in `_qc_removed.tsv` with reason "
                f"`protein_count_mismatch:segment=<seg>:got=<n>:"
                f"expected_one_of=[...]`."
            )
        elif min_proteins is not None:
            extra_lines.append(
                f"**{_fmt_int(qc_report.removed_proteins)}** sequence(s) were "
                f"dropped by the protein-annotation check "
                f"(`qc.protein_annotation.min_proteins = {min_proteins}`). "
                f"Failing records appear in `_qc_removed.tsv` with reason "
                f"`protein_count_below_min:<n><{min_proteins}`."
            )
        else:
            extra_lines.append(
                f"**{_fmt_int(qc_report.removed_proteins)}** record(s) were "
                f"dropped by the protein-annotation check (see `_qc_removed.tsv` "
                f"for per-record reasons starting with `protein_count_`)."
            )
    if qc_report.removed_protein_quality:
        pq_cfg = (qc_cfg.get("protein_quality") or {})
        mbf = pq_cfg.get("max_bad_fraction", 0.05)
        unit = "isolates" if cfg.get("segmented", {}).get("enabled") else "sequences"
        extra_lines.append(
            f"**{_fmt_int(qc_report.removed_protein_quality)}** {unit} were "
            f"dropped by the protein-quality check — a CDS protein carried "
            f"more than **{mbf:.0%}** ambiguous residues (X/B/Z/J; the "
            f"amino-acid analogue of the nucleotide ambiguity filter), which "
            f"fails its segment and drops the whole "
            f"{'isolate' if unit == 'isolates' else 'sequence'}. Dropped "
            f"records appear in `_qc_removed.tsv` with reason "
            f"`protein_quality:...`."
        )
    if qc_report.removed_length_by_segment:
        per_seg = qc_report.removed_length_by_segment
        # The per-segment counter is in isolates (one bad segment drops the
        # whole isolate), so summing too_short + too_long gives the total
        # isolates lost — the actionable number for a bench scientist. The
        # raw qc_report.removed_length counter is in segments and would
        # mislead the reader, so we deliberately do not surface it here.
        total_isolates = sum(
            c.get("too_short", 0) + c.get("too_long", 0)
            for c in per_seg.values()
        )
        breakdown_bits: list[str] = []
        for seg_name in sorted(per_seg.keys()):
            counts = per_seg[seg_name]
            if counts.get("too_short"):
                breakdown_bits.append(f"{seg_name} too short: {counts['too_short']}")
            if counts.get("too_long"):
                breakdown_bits.append(f"{seg_name} too long: {counts['too_long']}")
        breakdown = "; ".join(breakdown_bits)
        extra_lines.append(
            f"**{_fmt_int(total_isolates)}** isolate(s) were dropped by the "
            f"per-segment length filter (configured bounds in "
            f"`segmented.viruses.*.segment_lengths`): {breakdown}."
        )
    if qc_report.removed_incomplete_isolates:
        extra_lines.append(
            f"**{_fmt_int(qc_report.removed_incomplete_isolates)}** isolates "
            f"were dropped during segmented-completeness filtering (missing "
            f"one or more expected segments, or missing a marker protein "
            f"required for amino-acid clustering)."
        )
    if qc_report.removed_taxonomy_mismatch:
        rank = (
            (cfg.get("segmented", {}) or {})
            .get("taxonomy_consistency", {})
            .get("rank", "species")
        )
        extra_lines.append(
            f"**{_fmt_int(qc_report.removed_taxonomy_mismatch)}** segments "
            f"(belonging to isolates whose segments resolved to different "
            f"taxa at the **{rank}** rank) were removed by the taxonomy-"
            f"consistency check — typically reassortants or `/isolate` "
            f"identifier collisions across distinct viruses."
        )
    if qc_report.removed_strain_collisions:
        extra_lines.append(
            f"**{_fmt_int(qc_report.removed_strain_collisions)}** records "
            f"were dropped by the strain-collision detector "
            f"(`/strain`-derived `isolate_id`s that collided across "
            f"different submitters)."
        )
    if qc_report.removed_extra_segments:
        extra_lines.append(
            f"**{_fmt_int(qc_report.removed_extra_segments)}** isolate(s) "
            f"were dropped by the extra-segments check "
            f"(`segmented.extra_segments_action: drop`) — their segment "
            f"set carried names outside the configured `segments` list, "
            f"and the entire isolate was removed (every segment recorded "
            f"in `_qc_removed.tsv` with reason `extra_segments:<extras>`)."
        )
    if qc_report.removed_hmm_failed:
        hmm_rt = cfg.get("_hmm_runtime", {}) or {}
        hmm_cfg = hmm_rt.get("hmm_cfg") or (cfg.get("hmm", {}) or {})
        breakdown_bits: list[str] = []
        for key, count in sorted(
            qc_report.removed_hmm_by_marker.items(),
            key=lambda kv: -kv[1],
        )[:5]:
            breakdown_bits.append(f"`{key}` ({_fmt_int(count)})")
        breakdown = ", ".join(breakdown_bits) if breakdown_bits else "—"
        unit = (
            "isolates" if cfg.get("segmented", {}).get("enabled") else "sequences"
        )
        rel = hmm_cfg.get("relative_length_cutoff", 0.5)
        ev = hmm_cfg.get("default_evalue", 1.0e-5)
        use_ga = hmm_cfg.get("use_ga_when_available", True)
        gate_phrase = (
            f"the Pfam gathering threshold (GA) where available "
            f"(otherwise E ≤ {ev:g})"
            if use_ga
            else f"E ≤ {ev:g}"
        )
        extra_lines.append(
            f"**{_fmt_int(qc_report.removed_hmm_failed)}** {unit} were "
            f"dropped by the HMM-based identity QC (HMMER hmmscan; "
            f"{_CITATIONS['hmmer']}) — each marker's configured HMM(s) "
            f"had to hit a CDS with {gate_phrase} AND the alignment span "
            f"had to cover ≥ {rel:.0%} of the HMM model length. Top "
            f"reasons: {breakdown}."
        )

    extra_block = ("\n\n" + " ".join(extra_lines)) if extra_lines else ""
    return f"## Quality control\n\n{basic}{extra_block}\n"


def _render_segmented(cfg: dict, qc_report: QCReport,
                      complete_isolates: Optional[dict],
                      segment_names: Optional[list]) -> str:
    if not cfg.get("segmented", {}).get("enabled"):
        return ""
    virus_cfg = (cfg.get("segmented", {}).get("viruses") or {}).get(
        cfg.get("segmented", {}).get("virus") or "", {}
    )
    expected = virus_cfg.get("expected_segments")
    segs = segment_names or virus_cfg.get("segments") or []
    seg_str = ", ".join(segs) if segs else "the expected segments"
    n_complete = len(complete_isolates) if complete_isolates else None
    n_final = qc_report.final_survivors if qc_report.final_survivors is not None else None
    return (
        f"## Segmented-virus handling\n\n"
        f"Sequences were processed in segmented-virus mode with "
        f"**{expected}** expected segments per isolate ({seg_str}). "
        f"Records were grouped by isolate using the GenBank `/isolate` "
        f"qualifier as the primary key, with `/strain` as a fallback "
        f"(provenance is recorded in the `isolate_id_source` column of "
        f"the per-isolate output TSV). "
        f"**{_fmt_int(n_complete)}** complete isolates carrying all "
        f"**{expected}** expected segments survived completeness "
        f"filtering, of which **{_fmt_int(n_final)}** reached the "
        f"clustering step after all QC stages.\n"
    )


def _render_selection(cfg: dict, result: RunResult, qc_report: QCReport) -> str:
    cluster_cfg = cfg.get("clustering", {}) or {}
    backend = cluster_cfg.get("backend", "mmseqs2")
    alphabet = cluster_cfg.get("alphabet_for_clustering", "protein")
    segmented = bool(cfg.get("segmented", {}).get("enabled"))
    n_input_to_clust = qc_report.final_survivors if qc_report.final_survivors is not None else None
    n_reps = len(result.representatives)
    n_clusters = len(result.clusters)
    priority = ", ".join(cfg.get("representative", {}).get("priority", []))

    if backend == "cdhit":
        tool_name = "cd-hit" if alphabet == "protein" else "cd-hit-est"
        tool_cite = _CITATIONS["cd-hit"]
    else:
        tool_name = "MMseqs2"
        tool_cite = _CITATIONS["mmseqs"]

    if alphabet == "protein" and segmented:
        input_desc = (
            f"**{_fmt_int(n_input_to_clust)} per-isolate concatenated "
            f"marker-protein sequences** (one marker per segment — longest "
            f"CDS by default, overridable via `cluster_protein` aliases — "
            f"joined in canonical segment order)"
        )
    elif alphabet == "protein" and cluster_cfg.get("concatenate_markers"):
        input_desc = (
            f"**{_fmt_int(n_input_to_clust)} concatenated marker-protein "
            f"sequences** (the marker CDS from every `cluster_protein` "
            f"spec joined in declared order — `concatenate_markers: true`)"
        )
    elif alphabet == "protein":
        input_desc = (
            f"**{_fmt_int(n_input_to_clust)} marker-protein sequences** "
            f"(longest CDS per record by default, overridable via "
            f"`cluster_protein` aliases)"
        )
    elif segmented:
        input_desc = (
            f"**{_fmt_int(n_input_to_clust)} per-isolate concatenated "
            f"nucleotide sequences** (per-segment sequences joined "
            f"head-to-tail in canonical segment order)"
        )
    else:
        input_desc = (
            f"**{_fmt_int(n_input_to_clust)} input nucleotide sequences**"
        )

    mode_desc = _describe_mode(result.mode, cfg)
    rep_unit = "representative isolates" if segmented else "representative sequences"

    diversity_cutoffs = cluster_cfg.get("diversity_curve_cutoffs") or []
    if diversity_cutoffs:
        cutoffs_str = ", ".join(f"{c:g}" for c in sorted(diversity_cutoffs, reverse=True))
        diversity_sentence = (
            f" For comparison and as a diagnostic of within-stratum "
            f"sequence diversity, the same clustering backend was also "
            f"run at fixed identity thresholds ({cutoffs_str}) for each "
            f"clustered stratum; the resulting cluster counts are "
            f"reported in `{{prefix}}_group_counts.tsv` and did not "
            f"influence representative selection."
        )
    else:
        diversity_sentence = ""

    hmm_sentence = ""
    hmm_rt = cfg.get("_hmm_runtime", {}) or {}
    if hmm_rt.get("active"):
        hmm_cfg = cfg.get("hmm", {}) or {}
        user_db = hmm_cfg.get("database")
        db_desc = (
            "the bundled viral-core profile set (Pfam-A subset; "
            f"{_CITATIONS['pfam']})"
            if not user_db
            else f"a user-supplied HMM database (`{user_db}`)"
        )
        # The drop count is already detailed in the Quality control section's
        # HMM-QC bullet; here we just reference the upstream gate and note
        # how it shaped the input to clustering.
        hmm_sentence = (
            f" The clustering input was pre-filtered by an HMM-based "
            f"identity QC step "
            f"(HMMER hmmscan v{detect_tool_versions().get('hmmscan') or '?'}; "
            f"{_CITATIONS['hmmer']}) against {db_desc} — see the Quality "
            f"control section above for the per-marker drop breakdown. "
            f"When `alphabet_for_clustering=protein`, the marker CDS for "
            f"each surviving sequence/segment is then chosen as the "
            f"longest CDS satisfying any of the configured HMM tokens."
        )
        if cluster_cfg.get("concatenate_markers") and not segmented:
            hmm_sentence += (
                " With `concatenate_markers: true`, one such marker CDS is "
                "selected per `cluster_protein` spec and their amino-acid "
                "sequences are concatenated in declared spec order to form "
                "the clustering string; a sequence missing any required "
                "marker is dropped."
            )

    return (
        f"## Representative selection\n\n"
        f"Clustering was performed on {input_desc} using **{tool_name}** "
        f"({tool_cite}). {mode_desc}{hmm_sentence} Within each cluster the "
        f"representative was chosen by the configured priority "
        f"(**{priority}**, with sequence length as the final tiebreaker). "
        f"The final set contains **{_fmt_int(n_reps)} {rep_unit}** "
        f"across **{_fmt_int(n_clusters)} cluster(s)**.{diversity_sentence} "
        f"Taxonomic diversity at each rank before and after clustering "
        f"(distinct species, genera, families, etc., with per-taxon "
        f"{rep_unit.split()[-1]} counts — the top 20 by member count for "
        f"high-diversity ranks) is reported in "
        f"`{{prefix}}_taxonomic_report.txt`, a parallel "
        f"per-protein coverage and length-statistics report (one column "
        f"per declared marker / extra_protein, four sub-tables per rank "
        f"from subgenus to class, plus a trailing HMM-architecture "
        f"section) is written to `{{prefix}}_protein_taxonomic_report.txt`, "
        f"and a corresponding per-rank **NT length-statistics report** "
        f"(per-segment columns + a `total` column in segmented mode, a "
        f"single `genome` column otherwise) is written to "
        f"`{{prefix}}_nucleotide_taxonomic_report.txt`. "
        f"Each of these four `.txt` reports is accompanied by a "
        f"machine-readable **tidy long-format TSV** "
        f"(`{{prefix}}_*_taxonomic_report.tsv`) — one row per "
        f"(rank, pool, taxon, spec, metric) observation, with a "
        f"shared 8-column schema across all four files so downstream "
        f"analysis (Excel pivot tables, R/pandas) can read them as-is. "
        f"For each declared marker spec "
        f"(`clustering.cluster_protein` in non-segmented mode, "
        f"`virus.segment_markers` / `virus.cluster_protein` in segmented "
        f"mode), the satisfying CDS from each representative was also "
        f"written, unaligned, to "
        f"`{{prefix}}_per_protein_fasta/{{prefix}}_<marker>.fasta` "
        f"(one file per marker, headers identical to the all-proteins "
        f"FASTA above) so the per-marker protein set can be re-aligned "
        f"or rebuilt into trees without re-extracting. "
        f"Any **accessory proteins** declared in `extra_protein:` "
        f"(proteins that are sparse across taxa and therefore unsuitable "
        f"for clustering — e.g. coronavirus ORF7) were extracted with "
        f"the same selection chain and written to "
        f"`{{prefix}}_extra_protein_fasta/{{prefix}}_<name>.fasta`; "
        f"under `--per-protein-phylo` a separate tree per entry is also "
        f"emitted into `{{prefix}}_extra_protein/`.\n"
    )


def _describe_mode(mode: str, cfg: dict) -> str:
    """One-sentence prose description of the selection mode."""
    descriptions = {
        "global":      "A single global clustering pass was applied.",
        "taxonomic1":  "Sequences were stratified by taxonomic rank and a "
                       "binary search over the identity threshold selected a "
                       "target number of representatives per group.",
        "taxonomic2":  "Sequences were stratified hierarchically across "
                       "multiple taxonomic ranks; representatives were "
                       "selected at each nested level.",
        "host":        "Sequences were stratified by host organism and a "
                       "binary search over the identity threshold selected a "
                       "target number of representatives per host.",
        "temporal":    "Sequences were stratified by collection time window "
                       "and a binary search over the identity threshold "
                       "selected a target number of representatives per "
                       "window.",
        "geographic":  "Sequences were stratified by country of origin and a "
                       "binary search over the identity threshold selected a "
                       "target number of representatives per country.",
        "custom":      "Sequences were stratified by a user-defined metadata "
                       "field and a binary search over the identity threshold "
                       "selected a target number of representatives per "
                       "group.",
        "hybrid":      "Sequences were stratified across multiple metadata "
                       "dimensions simultaneously; representatives were "
                       "selected per multi-dimensional group.",
    }
    return descriptions.get(mode, f"Selection mode: `{mode}`.")


def _read_iqtree_model_file(cfg: dict) -> dict[str, str]:
    """Load ``{prefix}_iqtree_model.txt`` into ``{label: model}``.

    The file is written by the phylo pipeline after parsing the actual
    ModelFinder picks out of IQ-TREE's ``.iqtree`` report. Returns ``{}``
    when the file is absent (FastTree run, IQ-TREE soft-failed, or the
    parser couldn't locate the picks) — the renderer falls back to the
    generic "ModelFinder for substitution-model selection" prose.
    """
    try:
        out_dir = Path(cfg["output"]["dir"])
        prefix = cfg["output"].get("prefix", "repseq")
    except (KeyError, TypeError):
        return {}
    path = out_dir / f"{prefix}_iqtree_model.txt"
    if not path.exists():
        return {}
    chosen: dict[str, str] = {}
    try:
        for line in path.read_text().splitlines():
            line = line.strip()
            if not line or ":" not in line:
                continue
            label, model = line.split(":", 1)
            label = label.strip()
            model = model.strip()
            if label and model:
                chosen[label] = model
    except OSError:
        return {}
    return chosen


def _render_phylo(cfg: dict, result: RunResult, phylo_ran: bool,
                  versions: dict, per_protein_ran: bool = False,
                  per_segment_ran: bool = False,
                  conservation_ran: bool = False,
                  pre_cluster_ran: bool = False) -> str:
    if (not phylo_ran and not per_protein_ran
            and not per_segment_ran and not conservation_ran
            and not pre_cluster_ran):
        return ""
    n_reps = len(result.representatives)
    alphabet = cfg.get("clustering", {}).get("alphabet_for_clustering", "protein")
    phylo_cfg = cfg.get("phylo", {}) or {}
    tool_pref = phylo_cfg.get("tool", "auto")
    # Match the runtime auto-pick: protein → IQ-TREE, NT → FastTree.
    if tool_pref == "auto":
        chosen = "iqtree2" if alphabet == "protein" else "FastTree"
    else:
        chosen = "iqtree2" if tool_pref == "iqtree" else "FastTree"
    chosen_models = _read_iqtree_model_file(cfg) if phylo_ran else {}
    if chosen == "iqtree2":
        # Append the actual ModelFinder pick when we have it — buried in
        # the .iqtree report otherwise. For non-partitioned the file has
        # a single GENOME entry; for partitioned the per-partition picks
        # land in the partitioned paragraph below instead.
        model_pick_clause = ""
        if chosen_models and len(chosen_models) == 1 and "GENOME" in chosen_models:
            model_pick_clause = (
                f" ModelFinder selected **{chosen_models['GENOME']}** "
                f"as the best-fit substitution model."
            )
        tree_sentence = (
            f"A maximum-likelihood tree was inferred with **IQ-TREE** "
            f"v{versions.get('iqtree2') or '?'} ({_CITATIONS['iqtree2']}) "
            f"using ModelFinder for substitution-model selection and "
            f"ultrafast bootstrap approximation (UFBoot).{model_pick_clause}"
        )
    else:
        flag_note = "`-nt -gtr`" if alphabet != "protein" else "default protein model"
        tree_sentence = (
            f"An approximate-maximum-likelihood tree was inferred with "
            f"**FastTree** v{versions.get('FastTree') or '?'} "
            f"({_CITATIONS['FastTree']}, {flag_note})."
        )
    rooting_cfg = (phylo_cfg.get("rooting", {}) or {})
    rooting = rooting_cfg.get("method", "auto")
    # Describe the rooting fallback chain honestly per method — the old
    # phrasing claimed "MAD and midpoint as fallbacks" for every method,
    # which lies for `outgroup` (falls back to midpoint), `none` (no
    # fallback at all), and `midpoint` (no fallback either).
    if rooting == "outgroup":
        og = rooting_cfg.get("outgroup")
        ogr = rooting_cfg.get("outgroup_rank")
        og_bits: list[str] = []
        if og:
            og_str = og if isinstance(og, str) else ", ".join(og)
            og_bits.append(f"accession(s) `{og_str}`")
        if ogr:
            og_bits.append(
                "rank "
                + ", ".join(f"{r}=`{t}`" for r, t in ogr.items())
            )
        target = " and ".join(og_bits) if og_bits else "user-specified outgroup"
        rooting_clause = (
            f"on the {target} (a multi-leaf spec roots at the MRCA), "
            f"with midpoint as a fallback if the spec matches no rep"
        )
    elif rooting == "midpoint":
        rooting_clause = "using midpoint rooting"
    elif rooting == "none":
        rooting_clause = (
            "left as parsed (no rooting applied — the input is assumed "
            "to be rooted already)"
        )
    else:
        rooting_clause = (
            f"using the **{rooting}** strategy (with minimum ancestor "
            f"deviation and midpoint as fallbacks where applicable)"
        )
    # Partitioned supermatrix is the default for protein + IQ-TREE runs; it
    # only applies when the HMM tier resolves >= 2 marker families, else the
    # pipeline falls back to concat-then-align (described in the else branch).
    part_cfg = phylo_cfg.get("partition", {}) or {}
    partition_active = (
        bool(part_cfg.get("enabled", True))
        and alphabet == "protein"
        and chosen == "iqtree2"
    )
    # Describe the actual MAFFT command rather than hardcoding `--auto`,
    # so `--fast` (`--retree 1`, no `--auto`) and any user override get
    # an honest write-up.
    mafft_cfg = phylo_cfg.get("mafft", {}) or {}
    mafft_use_auto = bool(mafft_cfg.get("use_auto", True))
    mafft_extra = list(mafft_cfg.get("extra_args", []) or [])
    if mafft_use_auto:
        genome_mafft_phrase = "`--auto`"
        if mafft_extra:
            genome_mafft_phrase = f"`--auto {' '.join(mafft_extra)}`"
    else:
        if mafft_extra == ["--retree", "1"]:
            genome_mafft_phrase = (
                "`--retree 1` (single-pass FFT-NS-1; the `--fast` "
                "preliminary-run setting — switch this off for the "
                "publication run)"
            )
        elif mafft_extra:
            genome_mafft_phrase = f"`{' '.join(mafft_extra)}`"
        else:
            genome_mafft_phrase = "(no `--auto`)"
    # Optional trimAl trimming, described where it sits (between MAFFT and
    # the tree). Whole-genome (phylo.trimal) and per-protein
    # (phylo.per_protein.trimal) have independent switches.
    _gt = phylo_cfg.get("trimal", {}) or {}
    if _gt.get("enabled"):
        _gm = _gt.get("mode", "automated1")
        genome_trim_clause = (
            f"The alignment was then trimmed with **trimAl** "
            f"({_CITATIONS.get('trimal', 'Capella-Gutiérrez et al. 2009')}, "
            f"`-{_gm}`) to drop poorly-aligned / gap-rich columns before tree "
            f"inference, with the untrimmed MAFFT alignment retained as "
            f"`{{prefix}}_msa_untrimmed.fasta`. "
        )
        genome_trim_partition_clause = (
            f"Each per-family alignment was trimmed with **trimAl** "
            f"({_CITATIONS.get('trimal', 'Capella-Gutiérrez et al. 2009')}, "
            f"`-{_gm}`) before concatenation, so the partition column ranges "
            f"reflect the trimmed widths (untrimmed alignments retained). "
        )
    else:
        genome_trim_clause = ""
        genome_trim_partition_clause = ""
    paragraphs = ["## Phylogenetic inference\n"]
    if phylo_ran and partition_active:
        linkage = part_cfg.get("linkage", "proportional")
        linkage_flag = {
            "proportional": "-p", "equal": "-q", "unlinked": "-Q",
        }.get(linkage, "-p")
        # Surface the per-partition ModelFinder picks if available. We
        # exclude any spurious "GENOME" key (only emitted for the
        # non-partitioned single-MFP path) so a stale sidecar from a
        # prior run can't poison this paragraph.
        partition_picks = {
            k: v for k, v in chosen_models.items() if k != "GENOME"
        }
        if partition_picks:
            picks_str = ", ".join(
                f"{k}=`{v}`" for k, v in partition_picks.items()
            )
            partition_picks_clause = (
                f" ModelFinder selected: {picks_str} (also written to "
                f"`{{prefix}}_iqtree_model.txt`)."
            )
        else:
            partition_picks_clause = ""
        paragraphs.append(
            f"The whole-genome tree of the **{_fmt_int(n_reps)}** "
            f"representatives was built as a **partitioned supermatrix**. "
            f"Each declared marker family (the `hmms:` domain-architecture "
            f"tokens used for QC) was aligned **separately** with **MAFFT** "
            f"v{versions.get('mafft') or '?'} ({_CITATIONS['mafft']}, "
            f"{genome_mafft_phrase}). {genome_trim_partition_clause}The per-family "
            f"alignments were concatenated "
            f"column-wise into one supermatrix. A maximum-likelihood tree was "
            f"then inferred with **IQ-TREE** v{versions.get('iqtree2') or '?'} "
            f"({_CITATIONS['iqtree2']}) fitting a substitution model **per "
            f"partition** (ModelFinder) under **{linkage}**-linkage branch "
            f"lengths (`{linkage_flag}`) with ultrafast bootstrap (UFBoot)."
            f"{partition_picks_clause} "
            f"This avoids aligning unrelated proteins across segment seams and "
            f"forcing a single model across distinct protein families. The "
            f"tree was rooted {rooting_clause} and internal nodes "
            f"annotated with the last common ancestor (LCA) of their "
            f"descendants. Outputs include the supermatrix alignment "
            f"(`{{prefix}}_msa.fasta`), the per-family alignments "
            f"(`{{prefix}}_msa_<family>.fasta`), a NEXUS partition file "
            f"(`{{prefix}}_partition.nex`), the Newick, and the PhyloXML. When "
            f"fewer than two marker families are resolvable, the pipeline "
            f"falls back to a single concatenated-marker alignment under one "
            f"model.\n"
        )
    elif phylo_ran:
        paragraphs.append(
            f"A multiple sequence alignment of the **{_fmt_int(n_reps)}** "
            f"representative sequences was built with **MAFFT** "
            f"v{versions.get('mafft') or '?'} ({_CITATIONS['mafft']}, "
            f"{genome_mafft_phrase}). {genome_trim_clause}{tree_sentence} The tree was rooted "
            f"{rooting_clause} and internal nodes "
            f"were annotated with the last common ancestor (LCA) of their "
            f"descendants. The final tree was emitted as PhyloXML alongside "
            f"the underlying Newick file and FASTA alignment.\n"
        )
    if per_protein_ran:
        pp_cfg = phylo_cfg.get("per_protein", {}) or {}
        min_taxa = max(3, int(pp_cfg.get("min_taxa", 3) or 3))
        pp_mafft = list((pp_cfg.get("mafft", {}) or {}).get("extra_args", []) or [])
        pp_trim = pp_cfg.get("trimal", {}) or {}
        pp_trim_clause = (
            f" (then trimmed with **trimAl** `-{pp_trim.get('mode', 'automated1')}`)"
            if pp_trim.get("enabled") else ""
        )
        # Recognise specific common strategies so the prose stays honest:
        # `--retree 1` is the `--fast` preliminary-run setting, and
        # `--maxiterate ... --localpair` is L-INS-i. Anything else is
        # described as "user-supplied flags" — describing arbitrary args as
        # L-INS-i would be a lie. (L-INS-i strictly requires --localpair;
        # --maxiterate alone is iterative-refinement on FFT-NS, not L-INS-i.)
        if pp_mafft == ["--retree", "1"]:
            align_sentence = (
                f"aligned with MAFFT (`--retree 1`; single-pass FFT-NS-1, "
                f"the `--fast` preliminary-run setting — switch this off for "
                f"the publication run){pp_trim_clause} and inferred with "
                f"{chosen}"
            )
        elif pp_mafft and "--localpair" in pp_mafft and "--maxiterate" in pp_mafft:
            align_sentence = (
                f"aligned with MAFFT (`{' '.join(pp_mafft)}`; high-accuracy "
                f"L-INS-i for these single-gene sets){pp_trim_clause} and "
                f"inferred with {chosen}"
            )
        elif pp_mafft:
            align_sentence = (
                f"aligned with MAFFT (`{' '.join(pp_mafft)}`; user-supplied "
                f"per-protein MAFFT flags){pp_trim_clause} and inferred with "
                f"{chosen}"
            )
        else:
            align_sentence = (
                f"aligned{pp_trim_clause} and inferred with the same "
                f"MAFFT/{chosen} pipeline"
            )
        paragraphs.append(
            f"In addition, a **separate tree was built for each HMM "
            f"marker** declared for quality control (one tree per `hmms:` "
            f"spec). A marker's `hmms:` list holds **alternative domain "
            f"architectures** (e.g. coronavirus Spike as "
            f"`CoV_S1--CoV_S2` *or* `bCoV_S1_N--bCoV_S1_RBD--CoV_S2`); a CDS "
            f"satisfying **any** of them counts, so divergent forms of the "
            f"same protein land in one tree. For every marker, the "
            f"satisfying CDS translation was taken from each representative "
            f"that carries the architecture and {align_sentence}, with the "
            f"same rooting and "
            f"LCA annotation as above. Markers carried by fewer than "
            f"**{min_taxa}** representatives were skipped. Each tree was "
            f"emitted as PhyloXML with its alignment, Newick, and id-map "
            f"into the `{{prefix}}_per_protein/` subdirectory (named after "
            f"the marker's `name:`, or its domain-architecture token when "
            f"unnamed)."
            + (
                " Each leaf's protein carries its HMM **domain "
                "architecture** (a phyloXML `<domain_architecture>` of every "
                "domain hit with its E-value), which viewers such as "
                "Archaeopteryx render as domain boxes."
                if pp_cfg.get("domain_architecture", True) else ""
            )
            + f" Incongruence "
            f"between these single-marker trees (e.g. an L-segment "
            f"polymerase tree disagreeing with an M-segment glycoprotein "
            f"tree) is the expected signature of reassortment.\n"
        )
        if (phylo_cfg.get("per_protein", {}) or {}).get("incongruence", True):
            paragraphs.append(
                f"That incongruence was quantified with pairwise **unrooted "
                f"Robinson-Foulds (RF) distances** between the marker trees "
                f"(and the whole-genome tree, when `--phylo` also ran), "
                f"scored on each pair's shared taxa and written to "
                f"`{{prefix}}_per_protein/{{prefix}}_incongruence.tsv` "
                f"(RF = 0 for identical unrooted topologies, higher for "
                f"more topological disagreement; normalised by the maximum "
                f"RF for the shared-taxon count).\n"
            )
    if per_segment_ran:
        paragraphs.append(
            f"A separate **nucleotide tree was also built per segment** "
            f"(one tree per declared segment, over the raw per-segment NT "
            f"of each representative isolate) into the "
            f"`{{prefix}}_per_segment/` subdirectory. This complements the "
            f"per-marker (per-protein) trees above: reassortment may not "
            f"surface in any single marker (one CDS per segment) but shows "
            f"up as topological incongruence between the per-segment trees "
            f"themselves. Each segment needed at least "
            f"**{max(3, int((phylo_cfg.get('per_protein', {}) or {}).get('min_taxa', 3) or 3))}** "
            f"carriers to be built. Same MAFFT / tree-builder / rooting / "
            f"LCA / colour palette as 2E and 2F.\n"
        )
    coloring_cfg = (phylo_cfg.get("coloring", {}) or {})
    if coloring_cfg.get("enabled", True):
        cranks = list(coloring_cfg.get("ranks") or ["genus"])
        if len(cranks) >= 2:
            rank_sentence = (
                f"by **{cranks[0]}** (a distinct hue per {cranks[0]}), with "
                f"**{cranks[1]}** shaded within its parent {cranks[0]}'s hue"
            )
        else:
            rank_sentence = f"by **{cranks[0]}** (a distinct hue per {cranks[0]})"
        paragraphs.append(
            f"Tree leaves were coloured {rank_sentence}; unresolved taxa "
            f"were left a neutral grey. The same colour palette was applied "
            f"across every tree, and is stored in the PhyloXML as a node "
            f"font-colour property readable by tree viewers such as "
            f"Archaeopteryx.\n"
        )
    if pre_cluster_ran:
        paragraphs.append(
            "A **pre-cluster overview tree** of every post-QC sequence "
            "(one leaf per CONCAT isolate in segmented mode) was "
            "inferred with a single-pass pipeline — **MAFFT** "
            f"v{versions.get('mafft') or '?'} ({_CITATIONS['mafft']}) "
            "with `--retree 1` (FFT-NS-1, no `--auto`), followed by "
            f"**FastTree** v{versions.get('FastTree') or '?'} "
            f"({_CITATIONS['FastTree']}), midpoint-rooted, with no "
            "LCA annotation, trimAl, or bootstrap. The phyloXML "
            "`<name>` of each representative leaf is prefixed with "
            "`[repr] ` so the elected representatives can be picked "
            "out at a glance against the broader diversity of the "
            "input pool. Outputs: `{prefix}_pre_cluster_tree.nwk`, "
            "`{prefix}_pre_cluster_tree.xml`, and "
            "`{prefix}_pre_cluster_tree_id_map.tsv` (with an "
            "`is_rep` column for grep-without-XML use).\n"
        )
    if conservation_ran:
        paragraphs.append(
            "Per-marker **conservation plots** were rendered to "
            "`{prefix}_conservation/{prefix}_<family>.png`, one PNG per "
            "declared HMM marker spec. Each figure stacks two metric "
            "line charts — per-column **Shannon entropy** (bits, gaps "
            "excluded) and **fraction matching consensus** (non-gap "
            "rows whose residue equals the column mode), both smoothed "
            "with a **15-residue centered sliding window** so single-"
            "column spikes (especially from low-coverage columns) "
            "don't drown out real conservation patterns — over a "
            "**domain-architecture ribbon** drawn from the HMM hits on "
            "the longest satisfying CDS across representatives, "
            "projected from ungapped CDS coordinates onto MSA columns "
            "and labelled with the HMM profile name. Domain boxes use "
            "a stable golden-angle family colour so the same marker is "
            "visually identifiable across runs. Computed in-process "
            "from the per-protein MSAs already written by "
            "`--per-protein-phylo` — no fresh alignment is run.\n"
        )
    return "\n".join(paragraphs)


def _collect_polyprotein_specs_for_summary(cfg: dict) -> list[dict]:
    """Lightweight enumeration of polyprotein specs for the summary prose.

    Returns one dict per declared spec with the fields the renderer
    needs: ``name``, ``segment`` (or ``None``), ``cut_strategy``,
    ``peptides`` (list of peptide-name strings). Lives here rather than
    pulling :func:`repseq.polyprotein.collect_polyprotein_specs` so this
    module stays free of the cycle (summary is consumed by cli; the
    polyprotein package eventually imports config which eventually
    imports summary in some entry points). The shape is enough for
    prose; the writer module does the real heavy lifting.
    """
    out: list[dict] = []
    segmented = bool((cfg.get("segmented", {}) or {}).get("enabled"))
    if segmented:
        seg = cfg.get("segmented", {}) or {}
        virus_name = seg.get("virus")
        viruses = seg.get("viruses", {}) or {}
        virus = viruses.get(virus_name) or {}
        per_seg = virus.get("polyprotein") or {}
        seg_order = list(virus.get("segments") or [])
        for extra in per_seg.keys():
            if extra not in seg_order:
                seg_order.append(extra)
        for seg_name in seg_order:
            for entry in per_seg.get(seg_name, []) or []:
                if not isinstance(entry, dict):
                    continue
                name = (entry.get("name") or "").strip()
                if not name:
                    continue
                peptides = [
                    (p.get("name") or "").strip()
                    for p in (entry.get("peptides") or [])
                    if isinstance(p, dict) and (p.get("name") or "").strip()
                ]
                if len(peptides) < 2:
                    continue
                has_motif = any(
                    isinstance(p, dict) and (p.get("cleavage_motif") or "")
                    for p in (entry.get("peptides") or [])
                )
                default_strat = "motif" if has_motif else "bisect"
                out.append({
                    "name": name,
                    "segment": seg_name,
                    "cut_strategy": entry.get("cut_strategy") or default_strat,
                    "peptides": peptides,
                })
    else:
        for entry in (cfg.get("clustering", {}) or {}).get("polyprotein", []) or []:
            if not isinstance(entry, dict):
                continue
            name = (entry.get("name") or "").strip()
            if not name:
                continue
            peptides = [
                (p.get("name") or "").strip()
                for p in (entry.get("peptides") or [])
                if isinstance(p, dict) and (p.get("name") or "").strip()
            ]
            if len(peptides) < 2:
                continue
            has_motif = any(
                isinstance(p, dict) and (p.get("cleavage_motif") or "")
                for p in (entry.get("peptides") or [])
            )
            default_strat = "motif" if has_motif else "bisect"
            out.append({
                "name": name,
                "segment": None,
                "cut_strategy": entry.get("cut_strategy") or default_strat,
                "peptides": peptides,
            })
    return out


def _render_polyprotein(cfg: dict, per_protein_ran: bool = False) -> str:
    """Polyprotein cutting section.

    Fires when at least one ``clustering.polyprotein`` /
    ``virus.polyprotein`` spec is declared AND the HMM tier was active
    this session (otherwise :func:`write_polyprotein_outputs` soft-fails
    with a stderr line and emits no FASTAs, so the summary should match).

    When ``per_protein_ran`` is True the section also describes the
    per-peptide phylogenetic trees (one per spec × peptide) that
    ``--per-protein-phylo`` produced alongside the peptide FASTAs.
    """
    specs = _collect_polyprotein_specs_for_summary(cfg)
    if not specs:
        return ""
    hmm_active = bool((cfg.get("_hmm_runtime", {}) or {}).get("active"))
    if not hmm_active:
        return ""

    strategies_used = sorted({s["cut_strategy"] for s in specs})
    strat_phrases = {
        "boundary": (
            "**boundary** (each peptide spans its HMM hit's "
            "`ali_from`..`ali_to` verbatim; deterministic but lossy "
            "at the seams)"
        ),
        "bisect": (
            "**bisect** (cuts placed at the midpoint between adjacent "
            "located peptides; no residues dropped between two adjacent "
            "hits, but cut sites are geometric rather than biological. "
            "When a peptide is missing, the flanking peptides keep their "
            "HMM-hit boundary on the missing side — the gap is left "
            "unassigned rather than absorbed by either neighbour)"
        ),
        "motif": (
            "**motif** (cuts snap to the user-declared cleavage motif "
            "— e.g. `LQ` for coronavirus 3CL, `Q` for picornavirus 3C "
            "— within ±`motif_window_aa` of the bisect point, falling "
            "back to bisect when no motif lies in the window)"
        ),
    }
    strat_clause = "; ".join(
        strat_phrases.get(s, s) for s in strategies_used
    )
    bullet_lines = []
    for s in specs:
        scope = f" (segment {s['segment']})" if s["segment"] else ""
        pep_list = ", ".join(s["peptides"])
        bullet_lines.append(
            f"* **{s['name']}**{scope} — {len(s['peptides'])} peptides "
            f"in N→C order: {pep_list}; cut strategy `{s['cut_strategy']}`."
        )
    bullets = "\n".join(bullet_lines)

    peptide_trees_clause = ""
    if per_protein_ran:
        peptide_trees_clause = (
            "\n\nWith `--per-protein-phylo` active, **one phylogenetic "
            "tree was also built per peptide** (one tree per declared "
            "spec × peptide combination) using the same MAFFT / "
            "IQ-TREE-or-FastTree / rooting / LCA / colour palette as "
            "the per-marker trees (`phylo.per_protein` settings apply "
            "verbatim). Peptide trees land alongside the FASTAs under "
            "`{prefix}_polyprotein/` with `{prefix}_<spec>_<peptide>_msa.fasta`, "
            "`_tree.nwk`, `_tree.xml`, `_tree_id_map.tsv` basenames. "
            "Sparse peptides (fewer than "
            f"`phylo.per_protein.min_taxa` `ok`-status slices) are "
            "skipped. Peptide trees are intentionally **kept out** of "
            "`{prefix}_per_protein/{prefix}_incongruence.tsv` (which "
            "compares whole-genome markers) — peptide-vs-peptide "
            "comparison within a polyprotein answers a different "
            "scientific question."
        )
    return (
        "## Polyprotein cutting\n\n"
        f"For each declared polyprotein spec, the elected representatives' "
        f"polyprotein CDS was sliced into its mature peptides using one "
        f"HMM per peptide. The parent CDS was identified by counting "
        f"peptide-HMM hits (with a configurable minimum number of distinct "
        f"matches), and the protein was cut according to the configured "
        f"strategy: {strat_clause}.\n\n"
        f"Declared specs:\n\n"
        f"{bullets}\n\n"
        f"Outputs are written under `{{prefix}}_polyprotein/` — one FASTA "
        f"per peptide of each spec, plus a `{{prefix}}_<spec>_peptides.tsv` "
        f"audit table recording every (representative × peptide) attempt "
        f"with status (`ok`, `missing`, `out_of_order`, `overlap`, "
        f"`no_parent_cds`) so the bench scientist can see which cuts are "
        f"clean and which carry caveats. Polyprotein cutting is purely "
        f"additive — the polyprotein itself still drives clustering and "
        f"the whole-genome tree. A per-rank **peptide coverage and "
        f"length-statistics report** (one column per declared peptide of "
        f"each spec, four sub-tables per rank from subgenus to class — "
        f"the sliced-peptide analogue of "
        f"`{{prefix}}_protein_taxonomic_report.txt`) is written to "
        f"`{{prefix}}_polyprotein_taxonomic_report.txt`, accompanied by "
        f"its tidy long-format TSV companion "
        f"`{{prefix}}_polyprotein_taxonomic_report.tsv` for downstream "
        f"analysis (Excel pivot tables, R/pandas)."
        f"{peptide_trees_clause}\n"
    )


def _render_software(cfg: dict, versions: dict, phylo_ran: bool) -> str:
    backend = cfg.get("clustering", {}).get("backend", "mmseqs2")
    alphabet = cfg.get("clustering", {}).get("alphabet_for_clustering", "protein")
    cdhit_binary = "cd-hit" if alphabet == "protein" else "cd-hit-est"
    hmm_active = bool((cfg.get("_hmm_runtime", {}) or {}).get("active"))

    def _row(name: str, version: Optional[str], role: str, cite: str) -> str:
        v = version if version else "—"
        return f"| {name} | {v} | {role} | {cite} |"

    rows = [
        _row("repseq",  REPSEQ_VERSION,            "Pipeline orchestration",  _CITATIONS["repseq"]),
        _row("Python",  _python_version(),         "Runtime",                 _CITATIONS["python"]),
        _row("BioPython", _biopython_version(),    "FASTA / GenBank parsing", _CITATIONS["biopython"]),
    ]
    if backend == "mmseqs2":
        rows.append(_row("MMseqs2", versions.get("mmseqs"), "Sequence clustering (used)",      _CITATIONS["mmseqs"]))
        rows.append(_row("cd-hit",  versions.get("cd-hit"), "Sequence clustering (not used)",  _CITATIONS["cd-hit"]))
    else:
        rows.append(_row(cdhit_binary, versions.get(cdhit_binary), "Sequence clustering (used)", _CITATIONS["cd-hit"]))
        rows.append(_row("MMseqs2",    versions.get("mmseqs"),     "Sequence clustering (not used)", _CITATIONS["mmseqs"]))
    if hmm_active:
        rows.append(_row("HMMER hmmscan", versions.get("hmmscan"), "HMM-based marker selection (used)", _CITATIONS["hmmer"]))
        rows.append(_row("Pfam-A", "—", "HMM profile source (bundled subset)", _CITATIONS["pfam"]))
    if phylo_ran:
        rows.append(_row("MAFFT",    versions.get("mafft"),    "Multiple sequence alignment", _CITATIONS["mafft"]))
        _phylo_cfg = cfg.get("phylo", {}) or {}
        _trim_on = (
            (_phylo_cfg.get("trimal", {}) or {}).get("enabled")
            or ((_phylo_cfg.get("per_protein", {}) or {}).get("trimal", {}) or {}).get("enabled")
        )
        if _trim_on:
            rows.append(_row("trimAl", versions.get("trimal"), "Alignment trimming (when enabled)", _CITATIONS["trimal"]))
        rows.append(_row("IQ-TREE",  versions.get("iqtree2"),  "Maximum-likelihood tree (when chosen)", _CITATIONS["iqtree2"]))
        rows.append(_row("FastTree", versions.get("FastTree"), "Approximate-ML tree (when chosen)",     _CITATIONS["FastTree"]))
    rows.append(_row("NCBI E-utilities", "—", "Taxonomy / CDS retrieval", _CITATIONS["ncbi"]))
    rows.append(_row("UniProt REST API", "—", "Reviewed-entry metadata",  _CITATIONS["uniprot"]))

    body = "\n".join(rows)
    return (
        "## Software and references\n\n"
        "| Tool | Version | Role | Reference |\n"
        "| --- | --- | --- | --- |\n"
        f"{body}\n"
    )


def _render_footer() -> str:
    return (
        "---\n\n"
        "*Auto-generated. The numbers above reflect what the pipeline "
        "actually did on this run. Treat this as a Methods-section "
        "starting point — edit for prose flow before submission, and "
        "verify the citations against the journal's preferred format.*\n"
    )


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------

def render_summary(
    cfg: dict[str, Any],
    qc_report: QCReport,
    result: RunResult,
    input_paths: list[str],
    complete_isolates: Optional[dict] = None,
    segment_names: Optional[list[str]] = None,
    phylo_ran: bool = False,
    per_protein_ran: bool = False,
    per_segment_ran: bool = False,
    conservation_ran: bool = False,
    pre_cluster_ran: bool = False,
    command: str = "",
) -> str:
    """Build the full Markdown summary as a single string."""
    versions = detect_tool_versions()
    any_phylo = (
        phylo_ran or per_protein_ran or per_segment_ran or pre_cluster_ran
    )
    sections = [
        _render_header(cfg, command, result.mode),
        _render_input(qc_report, input_paths),
        _render_qc(qc_report, cfg),
        _render_segmented(cfg, qc_report, complete_isolates, segment_names),
        _render_selection(cfg, result, qc_report),
        _render_phylo(
            cfg, result, phylo_ran, versions,
            per_protein_ran, per_segment_ran,
            conservation_ran=conservation_ran,
            pre_cluster_ran=pre_cluster_ran,
        ),
        _render_polyprotein(cfg, per_protein_ran=per_protein_ran),
        _render_software(cfg, versions, any_phylo),
        _render_footer(),
    ]
    # Drop empty sections (segmented + phylo are conditional).
    return "\n".join(s for s in sections if s.strip())


def write_summary(
    cfg: dict[str, Any],
    qc_report: QCReport,
    result: RunResult,
    input_paths: list[str],
    complete_isolates: Optional[dict] = None,
    segment_names: Optional[list[str]] = None,
    phylo_ran: bool = False,
    per_protein_ran: bool = False,
    per_segment_ran: bool = False,
    conservation_ran: bool = False,
    pre_cluster_ran: bool = False,
    command: str = "",
) -> Path:
    """Render the summary and write it to ``{prefix}_summary.md``.

    Returns the path on success. Caller is responsible for soft-fail
    wrapping if a render error must not kill the parent run.
    """
    out_dir = Path(cfg["output"]["dir"])
    prefix = cfg["output"].get("prefix", "repseq")
    path = out_dir / f"{prefix}_summary.md"
    md = render_summary(
        cfg, qc_report, result, input_paths,
        complete_isolates=complete_isolates,
        segment_names=segment_names,
        phylo_ran=phylo_ran,
        per_protein_ran=per_protein_ran,
        per_segment_ran=per_segment_ran,
        conservation_ran=conservation_ran,
        pre_cluster_ran=pre_cluster_ran,
        command=command,
    )
    path.write_text(md)
    return path
