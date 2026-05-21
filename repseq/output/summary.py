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
        f"`{{prefix}}_taxonomic_report.txt`.\n"
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


def _render_phylo(cfg: dict, result: RunResult, phylo_ran: bool,
                  versions: dict, per_protein_ran: bool = False) -> str:
    if not phylo_ran and not per_protein_ran:
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
    if chosen == "iqtree2":
        tree_sentence = (
            f"A maximum-likelihood tree was inferred with **IQ-TREE** "
            f"v{versions.get('iqtree2') or '?'} ({_CITATIONS['iqtree2']}) "
            f"using ModelFinder for substitution-model selection and "
            f"ultrafast bootstrap approximation (UFBoot)."
        )
    else:
        flag_note = "`-nt -gtr`" if alphabet != "protein" else "default protein model"
        tree_sentence = (
            f"An approximate-maximum-likelihood tree was inferred with "
            f"**FastTree** v{versions.get('FastTree') or '?'} "
            f"({_CITATIONS['FastTree']}, {flag_note})."
        )
    rooting = (phylo_cfg.get("rooting", {}) or {}).get("method", "taxonomy_guided")
    paragraphs = ["## Phylogenetic inference\n"]
    if phylo_ran:
        paragraphs.append(
            f"A multiple sequence alignment of the **{_fmt_int(n_reps)}** "
            f"representative sequences was built with **MAFFT** "
            f"v{versions.get('mafft') or '?'} ({_CITATIONS['mafft']}, "
            f"`--auto`). {tree_sentence} The tree was rooted using the "
            f"**{rooting}** strategy (with minimum ancestor deviation and "
            f"midpoint as fallbacks where applicable) and internal nodes "
            f"were annotated with the last common ancestor (LCA) of their "
            f"descendants. The final tree was emitted as PhyloXML alongside "
            f"the underlying Newick file and FASTA alignment.\n"
        )
    if per_protein_ran:
        pp_cfg = phylo_cfg.get("per_protein", {}) or {}
        min_taxa = max(3, int(pp_cfg.get("min_taxa", 3) or 3))
        pp_mafft = list((pp_cfg.get("mafft", {}) or {}).get("extra_args", []) or [])
        if pp_mafft:
            align_sentence = (
                f"aligned with MAFFT (`{' '.join(pp_mafft)}`; high-accuracy "
                f"L-INS-i for these single-gene sets) and inferred with "
                f"{chosen}"
            )
        else:
            align_sentence = f"aligned and inferred with the same MAFFT/{chosen} pipeline"
        paragraphs.append(
            f"In addition, a **separate tree was built for each HMM "
            f"domain-architecture marker** declared for quality control "
            f"(the `hmms:` tokens; e.g. `Bunya_G1--Bunya_G2`). For every "
            f"marker, the satisfying CDS translation was taken from each "
            f"representative that carries that architecture and {align_sentence}, "
            f"with the same rooting and "
            f"LCA annotation as above. Markers carried by fewer than "
            f"**{min_taxa}** representatives were skipped. Each tree was "
            f"emitted as PhyloXML with its alignment, Newick, and id-map "
            f"into the `{{prefix}}_per_protein/` subdirectory. Incongruence "
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
    return "\n".join(paragraphs)


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
    command: str = "",
) -> str:
    """Build the full Markdown summary as a single string."""
    versions = detect_tool_versions()
    sections = [
        _render_header(cfg, command, result.mode),
        _render_input(qc_report, input_paths),
        _render_qc(qc_report, cfg),
        _render_segmented(cfg, qc_report, complete_isolates, segment_names),
        _render_selection(cfg, result, qc_report),
        _render_phylo(cfg, result, phylo_ran, versions, per_protein_ran),
        _render_software(cfg, versions, phylo_ran or per_protein_ran),
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
        command=command,
    )
    path.write_text(md)
    return path
