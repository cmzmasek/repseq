"""YAML configuration loading, validation, and defaults."""

from __future__ import annotations

import copy
import os
import re
from pathlib import Path
from typing import Any, Optional

import yaml

from .errors import ConfigError


# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

DEFAULTS: dict[str, Any] = {
    "cache_dir": "~/.repseq/cache",
    "temp_dir": "/tmp/repseq",
    "threads": 4,
    "seed": 42,
    "qc": {
        "remove_duplicates": True,
        "genome_length_filter": {
            # Absolute nucleotide-length bounds on whole input sequences.
            # NON-SEGMENTED MODE ONLY — enabling this with segmented.enabled
            # is a config error (segmented runs use per-segment bounds via
            # segmented.viruses.<v>.segment_lengths instead). Default off:
            # the user must opt in and supply at least one of min/max.
            "enabled": False,
            "min": None,   # drop sequences shorter than this many nt
            "max": None,   # drop sequences longer than this many nt
        },
        "ambiguous_threshold": 0.05,
        "protein_annotation": {
            # Drop sequences whose NCBI GenBank record has fewer than
            # min_proteins annotated CDS features. Requires network access
            # (skipped automatically with --no-resolve).
            "enabled": False,
            "min_proteins": 1,
        },
        "protein_quality": {
            # Amino-acid analogue of ambiguous_threshold. Drop a CDS
            # protein when the fraction of ambiguous residues (X/B/Z/J) in
            # its translation exceeds max_bad_fraction; a bad protein fails
            # its segment, which drops the whole isolate (segmented mode) or
            # the sequence (non-segmented). An empty/absent translation
            # counts as fully bad. Network-dependent: when enabled it
            # force-fetches GenBank CDS translations if no earlier step
            # already did. Skipped under --no-resolve.
            "enabled": False,
            "max_bad_fraction": 0.05,
        },
        "annotation_filter": {
            "enabled": True,
            "keywords": [
                "MAG:",
                "metagenome-assembled",
                "synthetic",
                "artificial",
                "fragment",
                "partial",
                "environmental sample",
                "uncultured",
                "unclassified",
                "unidentified",
                "hypothetical",
            ],
        },
    },
    "segmented": {
        "enabled": False,
        "virus": None,
        "viruses": {},
        # When true (the default), repseq fetches the GenBank source feature
        # for each NCBI-sourced sequence and uses its /isolate, /strain, and
        # /segment qualifiers in preference to the header-regex parse. The
        # regex still runs as a fallback for sequences without an accession,
        # for UniProt input, when --no-resolve is set, or when the GenBank
        # record lacks the qualifier. Set to false to bypass the GenBank
        # lookup entirely (header-regex only).
        "use_genbank_metadata": True,
        # Drop any segmented isolate whose segments disagree on the
        # taxonomic rank named in ``rank``. Reassortment between
        # different parent species is real biology for many segmented
        # viruses (peribunyaviruses, orthomyxoviruses, …), so the user
        # who turns this on is asking for monophyletic-at-this-rank
        # isolates only. Off-by-default would let those isolates
        # through unflagged, so the default is on. Missing values are
        # ignored — an isolate is only dropped when *populated* labels
        # disagree.
        "taxonomy_consistency": {
            "enabled": True,
            "rank": "species",
        },
        # What to do when the strain-collision detector finds two or more
        # distinct accessions sharing the same strain-derived isolate_id
        # AND the same segment — the over-merge signature of the
        # /strain → isolate_id fallback in _populate_genbank_isolate_segment.
        # "warn" (default) prints one line per collision to stderr and
        # lets the pipeline continue (the over-merged isolate keeps the
        # longest sequence per segment, the rest get dedup-dropped
        # downstream). "drop" removes every accession involved in any
        # collision before the completeness filter runs, adds them to
        # _qc_removed.tsv with reason "strain_collision:<segment>", and
        # increments QCReport.removed_strain_collisions.
        "strain_collision_action": "warn",
        # What to do when an isolate's segments include names outside
        # the configured ``segments`` list (e.g. an L/M/S virus that
        # has a fourth segment, or a non-canonical segment identifier
        # that ``identify_segment`` couldn't map). "warn" (default)
        # prints one line per affected isolate to stderr and pipes the
        # isolate through with just its expected segments — today's
        # silent prune becomes visible without changing analysis.
        # "drop" removes the entire isolate (every segment lands in
        # _qc_removed.tsv with reason "extra_segments:<extras>") and
        # increments QCReport.removed_extra_segments (units: isolates).
        "extra_segments_action": "warn",
    },
    "clustering": {
        "backend": "mmseqs2",              # "mmseqs2" | "cdhit"
        # Alphabet fed to the clustering backend.
        #
        # IMPORTANT: this setting only chooses what the clustering backend
        # (mmseqs2 / cd-hit) sees. It does NOT disable GenBank CDS download,
        # protein-count QC (qc.protein_annotation), or the
        # virus.expected_proteins_per_segment check — those run on every
        # isolate regardless of this value.
        #
        #   "protein"    — cluster on amino acid sequences (recommended for
        #                  diverged virus families). Non-segmented: the
        #                  marker protein (longest CDS, or the first
        #                  matching alias in `cluster_protein`).
        #                  Segmented: in-order concat of each segment's
        #                  marker protein. Triggers a one-shot GenBank CDS
        #                  fetch when proteins aren't already cached.
        #   "nucleotide" — cluster on the raw nucleotide sequence.
        #                  Non-segmented: the input FASTA sequence as-is.
        #                  Segmented: concat of all segments in
        #                  `segments` order.
        "alphabet_for_clustering": "protein",
        # Non-segmented marker-protein override. Alias list, matched
        # case-insensitively as substrings against /product. First alias
        # that matches a CDS wins; if no aliases match (or the list is
        # empty), the longest CDS on the sequence is used.
        "cluster_protein": [],
        # Additional proteins to extract and (optionally) phylogenize but
        # NOT use for clustering or the whole-genome tree. Same schema as
        # cluster_protein (list of {name, aliases?, hmms?} dicts). For
        # every entry, a single CDS per representative (HMM-gated when
        # `hmms` is declared and the HMM tier ran; alias fallback against
        # /product otherwise) is written to
        # {prefix}_extra_protein_fasta/{prefix}_<name>.fasta. With
        # --per-protein-phylo, one tree per entry is also emitted into
        # {prefix}_extra_protein/. Intended for proteins that are sparse
        # across taxa (e.g. coronavirus ORF7) — clustering them would
        # systematically drop representatives that lack the protein.
        "extra_protein": [],
        # Diagnostic: for each stratum where binary-search clustering ran,
        # also run the clustering backend at each of these identity
        # thresholds and report the cluster count as additional columns
        # (n_clusters_0.99, n_clusters_0.95, ...) in
        # {prefix}_group_counts.tsv. Reporting-only — does NOT influence
        # representative selection. Cutoffs below the backend's identity
        # floor (cd-hit-est: 0.80, cd-hit protein: 0.40, MMseqs2: 0) are
        # reported as NA. Set to [] to disable (e.g. on huge runs where
        # the extra clustering work would dominate runtime).
        "diversity_curve_cutoffs": [0.99, 0.95, 0.9, 0.8, 0.7],
        "mmseqs2_mode": "easy-linclust",   # "easy-linclust" | "easy-cluster"
        "coverage": 0.8,
        "coverage_mode": 0,
        "extra_args": [],
        "cdhit": {
            # Binary auto-selected from input alphabet:
            #   protein  -> cd-hit
            #   nucleic  -> cd-hit-est
            # Override to pin a specific path or variant.
            "binary": None,
            # Word size (-n). None = auto-pick from threshold per the cd-hit
            # user guide; cd-hit refuses out-of-range -n for a given -c.
            "word_size": None,
            # Shorter-sequence coverage (-aS); only honoured when
            # global_alignment is False (cd-hit's -G 0).
            "coverage": 0.8,
            # -G: True = global identity (cd-hit default), False = local.
            "global_alignment": True,
            # -g: False = greedy (fast, default), True = accurate
            # (slower; compares each input against every existing cluster).
            "accurate": False,
            # -M: memory cap in MB. 0 = unlimited (cd-hit default).
            "memory_mb": 0,
            # Raw cd-hit flags appended verbatim (e.g. ["-s", "0.8"]).
            "extra_args": [],
        },
    },
    "hmm": {
        # HMM-based marker-protein selection. When enabled AND a marker
        # has `hmms: [...]` configured AND hmmscan (HMMER) is on PATH,
        # HMM hits become the AUTHORITATIVE gate for that marker — a CDS
        # whose /product matches an alias but FAILS the HMM check is
        # rejected, and the segment / sequence is dropped. Markers
        # without `hmms` fall through to the legacy alias → longest
        # chain. Soft-fails (warns and falls back) when hmmscan is
        # missing or the database is unavailable.
        "enabled": True,
        # null = use the bundled viral-core set (`repseq/data/hmms/
        # repseq_viral_core.hmm`, 19 Pfam-A profiles for RdRp,
        # nucleocapsid, glycoprotein, protease, etc.). Absolute path = user-
        # supplied .hmm file; auto-`hmmpress`-ed on first use if the
        # .h3* index files are missing.
        "database": None,
        # E-value cutoff used when a profile has no curated Pfam GA
        # (gathering threshold). Hits with E ≤ this value pass.
        "default_evalue": 1.0e-5,
        # When true, use each profile's GA cutoff when available; fall
        # back to default_evalue otherwise. When false, always use
        # default_evalue regardless of GA availability.
        "use_ga_when_available": True,
        # Length cutoff: ali_span / hmm_model_length must be ≥ this.
        # Guards against tiny single-domain hits being accepted as full
        # marker matches. Range (0, 1].
        "relative_length_cutoff": 0.5,
        # Multidomain-token overlap tolerance, in amino acids. Consecutive
        # named domains in an N-to-C token (`A--B--C`) may overlap by up
        # to this many residues at their boundary AND still satisfy the
        # token. The walk also requires forward progression — the next
        # hit must start strictly C-terminal to the prior hit's start AND
        # end strictly C-terminal to the prior hit's end — so a single
        # fused region can never masquerade as two domains, and full
        # containment is always rejected. Set 0 for strict non-overlap
        # (pre-v0.22 behaviour). 30 (default) is comfortable for typical
        # Pfam viral-glycoprotein profiles whose boundaries don't align
        # exactly to the molecular cleavage site (e.g. coronavirus
        # S1/S2 across the furin site, bunya G1/G2).
        "multidomain_overlap_tolerance": 30,
        # null = use cfg.threads.
        "threads": None,
    },
    "representative": {
        "priority": ["refseq", "reviewed_uniprot", "longest"],
    },
    "phylo": {
        # Optional MSA + phylogeny step. Triggered with --phylo on any
        # mode subcommand; skipped automatically if fewer than 3
        # representatives survive selection.
        #
        # Tree-builder selection:
        #   "auto"     — IQ-TREE for protein alignments, FastTree for
        #                nucleotide. Default; matches what each tool is
        #                best at.
        #   "iqtree"   — always IQ-TREE (slower; ModelFinder + UFBoot)
        #   "fasttree" — always FastTree (faster; approximate-ML)
        "tool": "auto",
        "mafft": {
            # Raw mafft flags appended to "mafft --auto --thread N <input>".
            # Examples: ["--maxiterate", "1000"] for L-INS-i; or
            # ["--retree", "1"] for a faster pass on a very large input.
            "extra_args": [],
            # Pass --auto to MAFFT. Set false to drop it and rely on
            # extra_args alone (used by the --fast CLI flag to force
            # single-pass FFT-NS-1).
            "use_auto": True,
        },
        # Optional alignment trimming (trimAl) between MAFFT and the
        # tree-builder, for the whole-genome tree (2E). OFF by default.
        # In partitioned mode each per-family alignment is trimmed BEFORE
        # concatenation (so the partition charset ranges stay valid).
        # Soft-fails (loud warning + builds on the UNTRIMMED alignment)
        # when trimal is missing, errors, or strips the alignment to
        # nothing. The per-protein trees have their own knob at
        # phylo.per_protein.trimal.
        "trimal": {
            "enabled": False,
            # trimAl method → the `-<mode>` flag. "automated1" is trimAl's
            # heuristic best-for-ML-trees default; other column-trimming
            # methods: gappyout, strict, strictplus, nogaps, noallgaps.
            # Threshold trimming (-gt / -st / -cons) goes in extra_args.
            "mode": "automated1",
            # Raw trimal flags appended verbatim, e.g. ["-gt", "0.8"].
            "extra_args": [],
        },
        "fasttree": {
            # Raw FastTree flags appended to its argv. The protein /
            # nucleotide model is picked automatically from the rep
            # alphabet (default JTT for protein, -nt -gtr for nucleotide).
            "extra_args": [],
        },
        "iqtree": {
            # Binary auto-detected: tries iqtree2 first, then iqtree.
            # Set to a name or absolute path to pin a specific build.
            "binary": None,
            # Substitution model. "MFP" runs ModelFinder Plus (recommended
            # for protein — JTT/WAG/LG/etc. tested by BIC). Pin a model
            # like "LG+G4" or "JTT+G4" for a faster fixed-model pass.
            "model": "MFP",
            # Ultrafast bootstrap replicates. 0 disables bootstrap (faster).
            # IQ-TREE recommends >= 1000 when reporting branch support.
            "ultrafast_bootstrap": 1000,
            # Raw flags appended verbatim, e.g. ["-alrt", "1000"] for SH-aLRT.
            "extra_args": [],
        },
        # Partitioned-supermatrix tree (the principled multi-marker
        # analysis). When enabled AND the run resolves to protein + IQ-TREE
        # AND the HMM tier resolved >= 2 marker families, the whole-genome
        # tree (2E) is built by aligning each marker family separately, then
        # concatenating the per-family MSAs column-wise into a supermatrix
        # and letting IQ-TREE fit a model PER partition (-p/-q/-Q). This
        # replaces gluing the markers into one string + one model + one
        # MAFFT (which can align unrelated proteins across segment seams).
        # Soft-falls back to that concat-then-align path when the run can't
        # be partitioned (FastTree, no HMM families, < 2 families).
        "partition": {
            # Master switch. Default ON for protein + IQ-TREE runs; set
            # false to force the legacy concat-then-align behaviour.
            "enabled": True,
            # IQ-TREE partition linkage (Chernomor et al. 2016):
            #   "proportional" — -p, edge-linked proportional: one shared
            #                    branch-length set + a per-partition rate
            #                    multiplier. The standard default.
            #   "equal"        — -q, edge-equal: all partitions share branch
            #                    lengths (fewest parameters).
            #   "unlinked"     — -Q, edge-unlinked: each partition gets its
            #                    own branch lengths (most flexible; many more
            #                    parameters — useful when segments may have
            #                    different histories, e.g. reassortment).
            "linkage": "proportional",
            # Optional per-family substitution-model pin, keyed by the family
            # LABEL (segment_token, e.g. "L_RdRP_4" or "M_Bunya_G1--Bunya_G2";
            # the segment prefix is dropped in non-segmented runs). A family
            # absent here gets "MFP" → ModelFinder picks its model. Example:
            #   models: {L_RdRP_4: "LG+G4", S_Bunya_nucleocap: "WAG+G4"}
            "models": {},
        },
        # Per-leaf display labels on the phyloXML tree.
        # Supported placeholders: {species}, {genus}, {subgenus},
        # {subfamily}, {family}, {order}, {class}, {phylum}, {id},
        # {accession}, {host}, {strain}, {isolate_id}, {country},
        # {date}, {year}, {organism}.
        "labeling": {
            # Default for non-segmented runs.
            "format": "{species}|{id}|{host}",
            # Used when segmented.enabled is true. Falls back to {format}
            # if null. When {strain} is requested but the GenBank record
            # has no /strain qualifier, the writer substitutes
            # {isolate_id} (which segmented mode always has) so the
            # label never collapses to ``...||host``.
            "segmented_format": "{species}|{strain}|{host}",
            # Replace internal whitespace runs in each placeholder value
            # with underscores. Keeps the label round-trippable through
            # tree viewers that treat whitespace as a token boundary.
            "replace_whitespace": True,
            # When a placeholder resolves to empty (and isn't a {strain}
            # that can fall back to {isolate_id}), drop the placeholder
            # AND the single separator character immediately before it,
            # so the rendered label never contains ``||`` or trailing
            # ``|``. Set true to keep all separators verbatim.
            "keep_separator_on_empty": False,
        },
        # Tree rooting. Tree-building tools (FastTree, IQ-TREE) produce
        # unrooted trees; the post-processing step picks a root before
        # the writer ladderizes + serialises.
        #
        # method:
        #   auto       — try taxonomy-guided → MAD → midpoint, first
        #                success wins. The most-likely-correct default.
        #   taxonomy   — root at the branch that maximises mean LCA
        #                specificity of internal clades against the
        #                resolved NCBI lineages. Falls through to
        #                midpoint if no leaves carry lineage data.
        #   mad        — Minimal Ancestor Deviation (Tria et al. 2017).
        #                Pure-Python implementation; robust when
        #                taxonomy is sparse.
        #   midpoint   — Bio.Phylo's root_at_midpoint. Last-resort
        #                fallback; always succeeds.
        #   none       — leave the tree as parsed (use when the input
        #                is already rooted, e.g. by an outgroup).
        "rooting": {
            "method": "auto",
        },
        # Internal-node LCA labels. After rooting, every internal
        # clade is labelled with the lowest common ancestor of its
        # terminals (read from the resolved NCBI lineage) and the
        # label's rank attached as a PhyloXML <rank>.
        "lca": {
            # Master switch.
            "enabled": True,
            # Leaves whose lineage doesn't reach this rank are
            # excluded from the LCA *vote* — they stay on the tree
            # but don't pull internal labels toward an over-coarse
            # taxon. "none" disables the gate. Default "genus" suits
            # viral data, where lots of leaves lack species-level
            # classification.
            "min_rank": "genus",
            # An internal node is annotated only if at least this
            # fraction of its terminals carry usable lineage data.
            # Guards against a handful of well-annotated leaves
            # dictating the label of a much larger bare clade.
            "coverage_threshold": 0.5,
        },
        # PhyloXML writer knobs.
        "phyloxml": {
            # Override the <confidence type="..."> attribute. ``auto``
            # picks ``sh_like`` for FastTree and ``ufboot`` for IQ-TREE,
            # which matches what each tool actually produces by default.
            # Set explicitly if you pass non-default tree args (e.g.
            # IQ-TREE ``-b`` for classical bootstrap, in which case use
            # ``bootstrap``; ``-alrt`` only, use ``sh_alrt``).
            "confidence_type": "auto",
        },
        # Taxonomy-driven leaf colouring of the phyloXML output. Each
        # external node gets a <property ref="style:font_color"
        # datatype="xsd:token" applies_to="node">#RRGGBB</property> so a
        # viewer (Archaeopteryx) tints the leaf label by taxonomy. The
        # palette is shared across the whole-genome tree (2E) and every
        # per-protein tree (2F), so a taxon is the same colour everywhere.
        "coloring": {
            # Master switch. Default ON, colouring by genus.
            "enabled": True,
            # One or two taxonomy ranks (from the 9-rank ladder or
            # phylum/kingdom/superkingdom). ONE rank → each value gets a
            # distinct hue. TWO ranks [parent, child] → the child fans
            # out across hues of its parent's base colour (e.g.
            # [genus, subgenus]: subgenera as shades-of-genus). List the
            # coarser rank first.
            "ranks": ["genus"],
            # HSV saturation / value for the generated palette, in [0,1].
            # Defaults give legible mid-tone colours on a white canvas.
            "saturation": 0.65,
            "value": 0.90,
            # Colour for leaves whose rank is empty / unknown / na / ? /
            # etc. — a medium grey by default. Must be #RRGGBB.
            "missing_color": "#808080",
        },
        # Per-protein trees (2F), triggered with --per-protein-phylo.
        # One tree per declared HMM domain-architecture token (the same
        # hmms: used for QC). Requires the HMM tier to have run.
        "per_protein": {
            # A family is built only when at least this many
            # representatives carry the architecture. Never drops below
            # 3 (the tree-builder floor) regardless of what's set here.
            "min_taxa": 3,
            # After building the per-protein trees, write a pairwise
            # unrooted Robinson-Foulds distance table
            # ({prefix}_incongruence.tsv) quantifying topological
            # incongruence between the marker trees (and the
            # whole-genome tree from --phylo, when present). High RF
            # between marker trees is the reassortment/recombination
            # signal. Needs >= 2 trees to compare; soft-fails otherwise.
            "incongruence": True,
            # MAFFT args for the per-protein (single-gene) alignments.
            # Default is empty → MAFFT --auto (fast; size-adaptive), the
            # same as the whole-genome tree. For a high-accuracy
            # publication run set this to ["--maxiterate", "1000",
            # "--localpair"] (L-INS-i) — these single-gene alignments are
            # small enough to afford it, but it's noticeably slower, so
            # it's opt-in. A non-empty list is passed to MAFFT WITHOUT
            # --auto (so the explicit strategy takes effect).
            "mafft": {
                "extra_args": [],
            },
            # Emit a phyloXML <domain_architecture> on each per-protein
            # tree leaf, built from that CDS's HMM hits (every hit, with
            # its E-value as the domain confidence — so Archaeopteryx's
            # interactive E-value slider can filter them). Default on.
            "domain_architecture": True,
            # Optional trimAl trimming for the per-protein (single-gene)
            # alignments — independent of the whole-genome phylo.trimal.
            # OFF by default; same shape (mode → -<mode>, default
            # automated1; extra_args raw passthrough). Soft-fails to the
            # untrimmed alignment exactly like phylo.trimal.
            "trimal": {
                "enabled": False,
                "mode": "automated1",
                "extra_args": [],
            },
        },
    },
    "taxonomy": {
        "ncbi_email": None,
        "ncbi_api_key": None,
        "cache_ttl_days": 30,
    },
    "output": {
        "dir": "./repseq_output",
        "prefix": "repseq",
        # {prefix}_protein_taxonomic_report.txt — protein coverage and length
        # statistics per taxonomic rank (subgenus → class). max_breakdown caps
        # the rows per rank table; extra taxa are summarised on a "+N more"
        # line so column widths stay readable on a terminal. Set higher for
        # publication runs where the full per-rank table is expected.
        "protein_report": {
            "max_breakdown": 20,
        },
    },
}


# ---------------------------------------------------------------------------
# Required fields for specific scenarios
# ---------------------------------------------------------------------------

SEGMENTED_VIRUS_REQUIRED_FIELDS = ["expected_segments", "segments", "isolate_regex"]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _deep_merge(base: dict, override: dict) -> dict:
    """Recursively merge override into base, returning a new dict."""
    result = copy.deepcopy(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def _expand_paths(cfg: dict) -> dict:
    """Expand ~ in path fields."""
    for key in ("cache_dir", "temp_dir"):
        if cfg.get(key):
            cfg[key] = str(Path(cfg[key]).expanduser())
    if cfg.get("output", {}).get("dir"):
        cfg["output"]["dir"] = str(Path(cfg["output"]["dir"]).expanduser())
    return cfg


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def _format_yaml_error(path: Path, err: yaml.YAMLError) -> str:
    """Turn a PyYAML parse error into a one-glance, line-located message."""
    mark = getattr(err, "problem_mark", None)
    loc = f" at line {mark.line + 1}, column {mark.column + 1}" if mark else ""
    problem = getattr(err, "problem", None) or str(err).splitlines()[0]
    return (
        f"Config file {path} is not valid YAML{loc}: {problem}.\n"
        f"       Check the indentation (use spaces, not tabs) and quoting."
    )


def load_config(path: Optional[str | Path] = None) -> dict[str, Any]:
    """Load config from YAML file, merging over defaults.

    Raises :class:`~repseq.errors.ConfigError` (rendered without a traceback
    at the CLI boundary) for the mistakes a user is likely to make: a path
    that doesn't exist, a file that isn't valid YAML, or a top level that
    isn't a mapping.
    """
    cfg = copy.deepcopy(DEFAULTS)
    if path is not None:
        path = Path(path)
        if not path.exists():
            raise ConfigError(
                f"Config file not found: {path}\n"
                f"       Check the path you passed to -c/--config."
            )
        try:
            with open(path) as fh:
                user_cfg = yaml.safe_load(fh)
        except yaml.YAMLError as e:
            raise ConfigError(_format_yaml_error(path, e)) from e
        except OSError as e:
            raise ConfigError(
                f"Could not read config file {path}: {e.strerror or e}."
            ) from e
        if user_cfg is None:
            user_cfg = {}
        if not isinstance(user_cfg, dict):
            raise ConfigError(
                f"Config file {path} must contain a mapping (key: value pairs) "
                f"at the top level, but found a {type(user_cfg).__name__}.\n"
                f"       Check the file's structure — a repseq config is a set "
                f"of named sections, not a bare list or value."
            )
        cfg = _deep_merge(cfg, user_cfg)

    # Environment variable overrides
    if os.environ.get("REPSEQ_NCBI_EMAIL"):
        cfg["taxonomy"]["ncbi_email"] = os.environ["REPSEQ_NCBI_EMAIL"]
    if os.environ.get("REPSEQ_NCBI_API_KEY"):
        cfg["taxonomy"]["ncbi_api_key"] = os.environ["REPSEQ_NCBI_API_KEY"]

    cfg = _expand_paths(cfg)
    return cfg


def _validate_hmm_tokens(hmms: Any, path: str) -> tuple[list[str], list[str]]:
    """Validate an ``hmms:`` list of token strings.

    Each token is either a single HMM name (``"Name"``) or a multidomain
    spec joined with ``--`` (``"A--B--C"``, HMMs listed in N-to-C order).
    Multiple tokens in one list are **alternative architectures (OR)** — a
    CDS satisfying any one of them satisfies the marker. Returns
    ``(errors, validated_tokens)``. Invalid tokens are dropped from
    ``validated_tokens`` so the caller can still check the "at least one of
    aliases / hmms" invariant.
    """
    from .hmm.runner import parse_hmm_token

    errs: list[str] = []
    validated: list[str] = []
    if not isinstance(hmms, list):
        errs.append(
            f"{path} must be a list of HMM token strings "
            "(single 'Name' or multidomain 'A--B--C')"
        )
        return errs, validated
    for i, token in enumerate(hmms):
        if not isinstance(token, str):
            errs.append(f"{path}[{i}] must be a string, got {type(token).__name__}")
            continue
        try:
            parse_hmm_token(token)
        except ValueError as e:
            errs.append(f"{path}[{i}]: {e}")
            continue
        validated.append(token)
    return errs, validated


def _validate_marker_entry(entry: Any, path: str) -> list[str]:
    """Validate one cluster_protein marker entry.

    Each entry is either a non-empty alias string (legacy, alias-only)
    or a dict with required ``name`` and at least one of ``aliases``
    (list of non-empty strings) or ``hmms`` (list of HMM token strings,
    where a token is either a single HMM name like ``"Name"`` or a
    multidomain spec like ``"A--B"`` in N-to-C order; multiple tokens are
    alternative architectures, OR).
    """
    errs: list[str] = []
    if isinstance(entry, str):
        if not entry.strip():
            errs.append(f"{path}: alias string must be non-empty")
        return errs
    if not isinstance(entry, dict):
        errs.append(
            f"{path} must be an alias string or a dict "
            "{name, aliases?, hmms?}"
        )
        return errs
    name = entry.get("name")
    if not isinstance(name, str) or not name.strip():
        errs.append(f"{path}: dict-form entry must include a non-empty 'name'")
    aliases = entry.get("aliases", [])
    if not isinstance(aliases, list) or not all(
        isinstance(a, str) and a.strip() for a in aliases
    ):
        errs.append(f"{path}: 'aliases' must be a list of non-empty strings")
        aliases = []
    token_errs, hmms = _validate_hmm_tokens(entry.get("hmms", []), f"{path}.hmms")
    errs.extend(token_errs)
    if not aliases and not hmms:
        errs.append(
            f"{path}: dict-form entry must define at least one of "
            "'aliases' or 'hmms' (otherwise the marker can't be matched)"
        )
    return errs


def _validate_segment_markers(
    sm: Any, virus_name: str, seg_names: set[str]
) -> list[str]:
    """Validate per-virus segment_markers block.

    Shape: ``{segment_name: {aliases: [...], hmms: [...]}}``. Each
    segment-spec must define at least one of aliases / hmms. Coexists
    with the legacy per-segment ``cluster_protein`` block; when both
    define a marker for the same segment, ``segment_markers`` wins.
    """
    errs: list[str] = []
    if not isinstance(sm, dict):
        errs.append(
            f"segmented.viruses.{virus_name}.segment_markers must be a "
            "mapping of segment-name → {aliases: [...], hmms: [...]}"
        )
        return errs
    for seg_name, spec in sm.items():
        prefix = f"segmented.viruses.{virus_name}.segment_markers.{seg_name}"
        if seg_name not in seg_names:
            errs.append(
                f"segmented.viruses.{virus_name}.segment_markers: "
                f"unknown segment '{seg_name}'"
            )
        if not isinstance(spec, dict):
            errs.append(
                f"{prefix} must be a dict with 'aliases' and/or 'hmms' keys"
            )
            continue
        aliases = spec.get("aliases", [])
        if not isinstance(aliases, list) or not all(
            isinstance(a, str) and a.strip() for a in aliases
        ):
            errs.append(f"{prefix}.aliases must be a list of non-empty strings")
            aliases = []
        token_errs, hmms = _validate_hmm_tokens(spec.get("hmms", []), f"{prefix}.hmms")
        errs.extend(token_errs)
        if not aliases and not hmms:
            errs.append(
                f"{prefix}: must define at least one of 'aliases' or 'hmms'"
            )
    return errs


def validate_config(cfg: dict[str, Any]) -> list[str]:
    """Return a list of validation error messages (empty = valid)."""
    errors: list[str] = []

    # Genome length filter (non-segmented only). The old median-percent
    # filter (qc.length_filter) was removed in favour of explicit absolute
    # bounds — reject the renamed key so configs migrate consciously rather
    # than silently losing their length filtering.
    qc = cfg.get("qc", {}) or {}
    if "length_filter" in qc:
        errors.append(
            "qc.length_filter was renamed to qc.genome_length_filter and the "
            "median-percent mode was removed. Use "
            "qc.genome_length_filter: {enabled: true, min: <nt>, max: <nt>} "
            "with absolute nucleotide bounds (non-segmented mode only)."
        )
    glf = qc.get("genome_length_filter", {}) or {}
    if "enabled" in glf and not isinstance(glf.get("enabled"), bool):
        errors.append("qc.genome_length_filter.enabled must be a boolean")
    if glf.get("enabled"):
        if cfg.get("segmented", {}).get("enabled"):
            errors.append(
                "qc.genome_length_filter.enabled cannot be true when "
                "segmented.enabled is true — the whole-genome length filter "
                "only applies to non-segmented runs. For segmented viruses, "
                "set per-segment bounds via "
                "segmented.viruses.<virus>.segment_lengths instead."
            )
        mn = glf.get("min")
        mx = glf.get("max")
        if mn is None and mx is None:
            errors.append(
                "qc.genome_length_filter.enabled is true but neither min nor "
                "max is set — supply at least one absolute nucleotide bound."
            )
        for label, val in (("min", mn), ("max", mx)):
            if val is not None and (not isinstance(val, int) or isinstance(val, bool) or val <= 0):
                errors.append(
                    f"qc.genome_length_filter.{label} must be a positive integer "
                    "(nucleotides) or null"
                )
        if (
            isinstance(mn, int) and not isinstance(mn, bool)
            and isinstance(mx, int) and not isinstance(mx, bool)
            and mn > mx
        ):
            errors.append(
                "qc.genome_length_filter.min must be <= qc.genome_length_filter.max"
            )

    # Ambiguous threshold
    thresh = cfg.get("qc", {}).get("ambiguous_threshold")
    if not isinstance(thresh, (int, float)) or not (0 <= thresh <= 1):
        errors.append("qc.ambiguous_threshold must be a number between 0 and 1")

    # Protein annotation QC
    pa = cfg.get("qc", {}).get("protein_annotation", {})
    if pa.get("enabled"):
        mp = pa.get("min_proteins")
        if not isinstance(mp, int) or mp < 0:
            errors.append("qc.protein_annotation.min_proteins must be a non-negative integer")

    # Protein quality QC
    pq = cfg.get("qc", {}).get("protein_quality", {})
    if pq.get("enabled"):
        mbf = pq.get("max_bad_fraction")
        if not isinstance(mbf, (int, float)) or not (0 <= mbf <= 1):
            errors.append(
                "qc.protein_quality.max_bad_fraction must be a number between 0 and 1"
            )

    # Segmented virus
    seg = cfg.get("segmented", {})
    if "use_genbank_metadata" in seg and not isinstance(
        seg["use_genbank_metadata"], bool
    ):
        errors.append("segmented.use_genbank_metadata must be a boolean")

    sca = seg.get("strain_collision_action", "warn")
    if sca not in ("warn", "drop"):
        errors.append(
            f"segmented.strain_collision_action '{sca}' is not supported "
            "(use 'warn' or 'drop')"
        )

    esa = seg.get("extra_segments_action", "warn")
    if esa not in ("warn", "drop"):
        errors.append(
            f"segmented.extra_segments_action '{esa}' is not supported "
            "(use 'warn' or 'drop')"
        )

    tc = seg.get("taxonomy_consistency", {}) or {}
    if "enabled" in tc and not isinstance(tc["enabled"], bool):
        errors.append("segmented.taxonomy_consistency.enabled must be a boolean")
    # The rank must be one a resolved TaxonomyInfo can actually answer.
    # Standard fields on TaxonomyInfo plus arbitrary entries from
    # ``lineage`` are both accepted by ``get_rank``, so we only
    # validate against the named-attribute set — anything else still
    # works, just relies on the lineage map being populated.
    valid_consistency_ranks = {
        "species", "subgenus", "genus", "subfamily", "family",
        "suborder", "order", "subclass", "class", "phylum",
        "kingdom", "superkingdom",
    }
    rank = tc.get("rank", "species")
    if not isinstance(rank, str) or rank.lower() not in valid_consistency_ranks:
        errors.append(
            f"segmented.taxonomy_consistency.rank '{rank}' is not supported "
            f"(use one of {sorted(valid_consistency_ranks)})"
        )
    if seg.get("enabled"):
        virus_name = seg.get("virus")
        if not virus_name:
            errors.append("segmented.virus must be set when segmented.enabled is true")
        else:
            viruses = seg.get("viruses", {})
            if virus_name not in viruses:
                errors.append(
                    f"segmented.virus '{virus_name}' not found in segmented.viruses"
                )
            else:
                vdef = viruses[virus_name]
                for field in SEGMENTED_VIRUS_REQUIRED_FIELDS:
                    if field not in vdef:
                        errors.append(
                            f"segmented.viruses.{virus_name} missing required field '{field}'"
                        )
                # Segment aliases: optional dict[canonical → list[str]]
                aliases = vdef.get("segment_aliases")
                if aliases is not None:
                    if not isinstance(aliases, dict):
                        errors.append(
                            f"segmented.viruses.{virus_name}.segment_aliases "
                            f"must be a mapping of canonical-segment-name → list of strings"
                        )
                    else:
                        seg_names = set(vdef.get("segments", []))
                        for canonical, syns in aliases.items():
                            if canonical not in seg_names:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"segment_aliases: unknown segment '{canonical}'"
                                )
                            if not isinstance(syns, list) or not all(
                                isinstance(s, str) and s.strip() for s in syns
                            ):
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"segment_aliases.{canonical} "
                                    f"must be a list of non-empty strings"
                                )

                cp = vdef.get("cluster_protein")
                if cp is not None:
                    if not isinstance(cp, dict):
                        errors.append(
                            f"segmented.viruses.{virus_name}.cluster_protein "
                            f"must be a mapping of segment-name → list of "
                            f"alias strings and/or {{name, aliases?, hmms?}} dicts"
                        )
                    else:
                        seg_names = set(vdef.get("segments", []))
                        for seg_name, entries in cp.items():
                            if seg_name not in seg_names:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"cluster_protein: unknown segment '{seg_name}'"
                                )
                            if not isinstance(entries, list):
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"cluster_protein.{seg_name} "
                                    f"must be a list of alias strings and/or "
                                    f"{{name, aliases?, hmms?}} dicts"
                                )
                            else:
                                for j, entry in enumerate(entries):
                                    errors.extend(_validate_marker_entry(
                                        entry,
                                        f"segmented.viruses.{virus_name}."
                                        f"cluster_protein.{seg_name}[{j}]",
                                    ))

                sm = vdef.get("segment_markers")
                if sm is not None:
                    errors.extend(_validate_segment_markers(
                        sm, virus_name, set(vdef.get("segments", []))
                    ))

                ep = vdef.get("extra_protein")
                if ep is not None:
                    if not isinstance(ep, dict):
                        errors.append(
                            f"segmented.viruses.{virus_name}.extra_protein "
                            f"must be a mapping of segment-name → list of "
                            f"{{name, aliases?, hmms?}} dicts"
                        )
                    else:
                        seg_names = set(vdef.get("segments", []))
                        seen_extra: set[str] = set()
                        for seg_name, entries in ep.items():
                            if seg_name not in seg_names:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"extra_protein: unknown segment "
                                    f"'{seg_name}'"
                                )
                            if not isinstance(entries, list):
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"extra_protein.{seg_name} must be a "
                                    f"list of {{name, aliases?, hmms?}} dicts"
                                )
                                continue
                            for j, entry in enumerate(entries):
                                ipath = (
                                    f"segmented.viruses.{virus_name}."
                                    f"extra_protein.{seg_name}[{j}]"
                                )
                                if isinstance(entry, str):
                                    errors.append(
                                        f"{ipath} must be a dict "
                                        f"{{name, aliases?, hmms?}}, not a "
                                        f"bare alias string"
                                    )
                                    continue
                                errors.extend(
                                    _validate_marker_entry(entry, ipath)
                                )
                                if isinstance(entry, dict):
                                    nm = (entry.get("name") or "").strip()
                                    if nm:
                                        # Names must be unique across the
                                        # whole virus (not just within one
                                        # segment) — the output basename is
                                        # the name alone, no segment prefix.
                                        key = nm.lower()
                                        if key in seen_extra:
                                            errors.append(
                                                f"{ipath}: duplicate name "
                                                f"'{nm}' — extra_protein "
                                                f"names must be unique "
                                                f"across all segments "
                                                f"(they key the output "
                                                f"filenames)"
                                            )
                                        seen_extra.add(key)

                epps = vdef.get("expected_proteins_per_segment")
                if epps is not None:
                    if not isinstance(epps, dict):
                        errors.append(
                            f"segmented.viruses.{virus_name}.expected_proteins_per_segment "
                            f"must be a mapping of segment-name → int or list[int]"
                        )
                    else:
                        seg_names = set(vdef.get("segments", []))
                        for seg_name, count in epps.items():
                            if seg_name not in seg_names:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"expected_proteins_per_segment: unknown segment '{seg_name}'"
                                )
                            # Allow int (exact) or list[int] (any-of).
                            # bool is a subclass of int in Python — reject it explicitly.
                            valid = (
                                (isinstance(count, int) and not isinstance(count, bool) and count >= 0)
                                or (
                                    isinstance(count, list)
                                    and len(count) > 0
                                    and all(
                                        isinstance(x, int)
                                        and not isinstance(x, bool)
                                        and x >= 0
                                        for x in count
                                    )
                                )
                            )
                            if not valid:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"expected_proteins_per_segment.{seg_name} "
                                    f"must be a non-negative integer or a non-empty "
                                    f"list of non-negative integers"
                                )

                sl = vdef.get("segment_lengths")
                if sl is not None:
                    if not isinstance(sl, dict):
                        errors.append(
                            f"segmented.viruses.{virus_name}.segment_lengths "
                            f"must be a mapping of segment-name → {{min: N, max: M}}"
                        )
                    else:
                        seg_names = set(vdef.get("segments", []))
                        for seg_name, bounds in sl.items():
                            if seg_name not in seg_names:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"segment_lengths: unknown segment '{seg_name}'"
                                )
                            if not isinstance(bounds, dict):
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"segment_lengths.{seg_name} must be a dict "
                                    f"with optional 'min' and/or 'max' integer keys"
                                )
                            else:
                                mn = bounds.get("min")
                                mx = bounds.get("max")
                                if mn is not None and (
                                    not isinstance(mn, int) or isinstance(mn, bool) or mn < 0
                                ):
                                    errors.append(
                                        f"segmented.viruses.{virus_name}."
                                        f"segment_lengths.{seg_name}.min must be a "
                                        f"non-negative integer"
                                    )
                                if mx is not None and (
                                    not isinstance(mx, int) or isinstance(mx, bool) or mx < 0
                                ):
                                    errors.append(
                                        f"segmented.viruses.{virus_name}."
                                        f"segment_lengths.{seg_name}.max must be a "
                                        f"non-negative integer"
                                    )
                                if (
                                    mn is not None and mx is not None
                                    and isinstance(mn, int) and isinstance(mx, int)
                                    and mn >= mx
                                ):
                                    errors.append(
                                        f"segmented.viruses.{virus_name}."
                                        f"segment_lengths.{seg_name}: min must be "
                                        f"less than max"
                                    )

    # Clustering backend
    backend = cfg.get("clustering", {}).get("backend")
    if backend not in ("mmseqs2", "cdhit"):
        errors.append(
            f"clustering.backend '{backend}' is not supported "
            f"(use 'mmseqs2' or 'cdhit')"
        )

    alphabet = cfg.get("clustering", {}).get("alphabet_for_clustering", "protein")
    if alphabet not in ("protein", "nucleotide"):
        errors.append(
            f"clustering.alphabet_for_clustering '{alphabet}' is not supported "
            f"(use 'protein' or 'nucleotide')"
        )

    cluster_protein_global = cfg.get("clustering", {}).get("cluster_protein", [])
    if not isinstance(cluster_protein_global, list):
        errors.append(
            "clustering.cluster_protein must be a list of alias strings "
            "and/or {name, aliases?, hmms?} dicts"
        )
    else:
        for i, entry in enumerate(cluster_protein_global):
            errors.extend(
                _validate_marker_entry(entry, f"clustering.cluster_protein[{i}]")
            )

    extra_protein_global = cfg.get("clustering", {}).get("extra_protein", [])
    if not isinstance(extra_protein_global, list):
        errors.append(
            "clustering.extra_protein must be a list of "
            "{name, aliases?, hmms?} dicts"
        )
    else:
        seen_extra_names: set[str] = set()
        for i, entry in enumerate(extra_protein_global):
            # extra_protein entries MUST be dicts with a name — they drive
            # output filenames, so a bare alias string (which has no name)
            # has nothing to key the artifact on.
            if isinstance(entry, str):
                errors.append(
                    f"clustering.extra_protein[{i}] must be a dict "
                    f"{{name, aliases?, hmms?}}, not a bare alias string"
                )
                continue
            errors.extend(
                _validate_marker_entry(entry, f"clustering.extra_protein[{i}]")
            )
            if isinstance(entry, dict):
                name = (entry.get("name") or "").strip()
                if name:
                    key = name.lower()
                    if key in seen_extra_names:
                        errors.append(
                            f"clustering.extra_protein[{i}]: duplicate name "
                            f"'{name}' — extra_protein names must be unique "
                            f"(they key the output filenames)"
                        )
                    seen_extra_names.add(key)

    diversity_cutoffs = cfg.get("clustering", {}).get("diversity_curve_cutoffs", [])
    if not isinstance(diversity_cutoffs, list):
        errors.append(
            "clustering.diversity_curve_cutoffs must be a list of floats in (0, 1]"
        )
    else:
        for v in diversity_cutoffs:
            if not isinstance(v, (int, float)) or not (0.0 < float(v) <= 1.0):
                errors.append(
                    f"clustering.diversity_curve_cutoffs entry {v!r} is not a "
                    f"number in (0, 1]"
                )
                break

    mmseqs2_mode = cfg.get("clustering", {}).get("mmseqs2_mode")
    if mmseqs2_mode not in ("easy-linclust", "easy-cluster"):
        errors.append(
            f"clustering.mmseqs2_mode must be 'easy-linclust' or 'easy-cluster', got '{mmseqs2_mode}'"
        )

    cdhit_cfg = cfg.get("clustering", {}).get("cdhit", {}) or {}
    if cdhit_cfg:
        ws = cdhit_cfg.get("word_size")
        if ws is not None and (
            not isinstance(ws, int) or isinstance(ws, bool) or not (2 <= ws <= 11)
        ):
            errors.append(
                "clustering.cdhit.word_size must be an integer in [2, 11] or null "
                "(auto-pick from threshold)"
            )
        cov = cdhit_cfg.get("coverage", 0.8)
        if not isinstance(cov, (int, float)) or not (0 <= cov <= 1):
            errors.append("clustering.cdhit.coverage must be a number between 0 and 1")
        mem = cdhit_cfg.get("memory_mb", 0)
        if not isinstance(mem, int) or isinstance(mem, bool) or mem < 0:
            errors.append(
                "clustering.cdhit.memory_mb must be a non-negative integer "
                "(0 = unlimited)"
            )
        for flag in ("global_alignment", "accurate"):
            if flag in cdhit_cfg and not isinstance(cdhit_cfg[flag], bool):
                errors.append(f"clustering.cdhit.{flag} must be a boolean")
        extra = cdhit_cfg.get("extra_args", [])
        if not isinstance(extra, list) or not all(isinstance(x, str) for x in extra):
            errors.append("clustering.cdhit.extra_args must be a list of strings")

    # Phylo
    phylo_cfg = cfg.get("phylo", {}) or {}
    tool = phylo_cfg.get("tool", "auto")
    if tool not in ("auto", "iqtree", "fasttree"):
        errors.append(
            f"phylo.tool '{tool}' is not supported "
            f"(use 'auto', 'iqtree', or 'fasttree')"
        )
    for tool_name in ("mafft", "fasttree", "iqtree"):
        tool_cfg = phylo_cfg.get(tool_name, {}) or {}
        extra = tool_cfg.get("extra_args", [])
        if not isinstance(extra, list) or not all(isinstance(x, str) for x in extra):
            errors.append(f"phylo.{tool_name}.extra_args must be a list of strings")

    iq_cfg = phylo_cfg.get("iqtree", {}) or {}
    if "model" in iq_cfg and not isinstance(iq_cfg["model"], str):
        errors.append("phylo.iqtree.model must be a string (e.g. 'MFP', 'LG+G4')")
    if "binary" in iq_cfg and iq_cfg["binary"] is not None and not isinstance(
        iq_cfg["binary"], str
    ):
        errors.append("phylo.iqtree.binary must be a string or null")
    ufb = iq_cfg.get("ultrafast_bootstrap", 0)
    if not isinstance(ufb, int) or isinstance(ufb, bool) or ufb < 0:
        errors.append(
            "phylo.iqtree.ultrafast_bootstrap must be a non-negative integer "
            "(0 disables; IQ-TREE recommends >= 1000 for interpretable support)"
        )

    labeling_cfg = phylo_cfg.get("labeling", {}) or {}
    for key in ("format", "segmented_format"):
        if key in labeling_cfg and labeling_cfg[key] is not None and not isinstance(
            labeling_cfg[key], str
        ):
            errors.append(f"phylo.labeling.{key} must be a string or null")
    for key in ("replace_whitespace", "keep_separator_on_empty"):
        if key in labeling_cfg and not isinstance(labeling_cfg[key], bool):
            errors.append(f"phylo.labeling.{key} must be a boolean")

    phyloxml_cfg = phylo_cfg.get("phyloxml", {}) or {}
    ct = phyloxml_cfg.get("confidence_type", "auto")
    if ct not in ("auto", "sh_like", "sh_alrt", "ufboot", "bootstrap"):
        errors.append(
            f"phylo.phyloxml.confidence_type '{ct}' is not supported "
            "(use 'auto', 'sh_like', 'sh_alrt', 'ufboot', or 'bootstrap')"
        )

    coloring_cfg = phylo_cfg.get("coloring", {}) or {}
    if "enabled" in coloring_cfg and not isinstance(coloring_cfg["enabled"], bool):
        errors.append("phylo.coloring.enabled must be a boolean")
    # Ranks usable for colouring — the 9-rank ladder plus the coarse
    # ranks TaxonomyInfo.get_rank answers from its named attributes.
    valid_color_ranks = {
        "species", "subgenus", "genus", "subfamily", "family",
        "suborder", "order", "subclass", "class", "phylum",
        "kingdom", "superkingdom",
    }
    cranks = coloring_cfg.get("ranks", ["genus"])
    if not isinstance(cranks, list) or not (1 <= len(cranks) <= 2):
        errors.append(
            "phylo.coloring.ranks must be a list of one or two rank names "
            "(one rank = distinct hue per value; two = [parent, child] with "
            "the child shaded within its parent's hue)"
        )
    else:
        for r in cranks:
            if not isinstance(r, str) or r.lower() not in valid_color_ranks:
                errors.append(
                    f"phylo.coloring.ranks contains unsupported rank '{r}' "
                    f"(use one of {sorted(valid_color_ranks)})"
                )
    for key in ("saturation", "value"):
        v = coloring_cfg.get(key, 0.65 if key == "saturation" else 0.90)
        if not isinstance(v, (int, float)) or isinstance(v, bool) or not (0 <= v <= 1):
            errors.append(f"phylo.coloring.{key} must be a number between 0 and 1")
    mc = coloring_cfg.get("missing_color", "#808080")
    if not (isinstance(mc, str) and re.fullmatch(r"#[0-9A-Fa-f]{6}", mc)):
        errors.append(
            "phylo.coloring.missing_color must be a '#RRGGBB' hex colour string"
        )

    pp_cfg = phylo_cfg.get("per_protein", {}) or {}
    mt = pp_cfg.get("min_taxa", 3)
    if not isinstance(mt, int) or isinstance(mt, bool) or mt < 1:
        errors.append(
            "phylo.per_protein.min_taxa must be a positive integer "
            "(values below 3 are clamped to 3 at runtime)"
        )
    if "incongruence" in pp_cfg and not isinstance(pp_cfg["incongruence"], bool):
        errors.append("phylo.per_protein.incongruence must be a boolean")
    pp_mafft = pp_cfg.get("mafft", {}) or {}
    pp_mafft_extra = pp_mafft.get("extra_args", [])
    if not isinstance(pp_mafft_extra, list) or not all(
        isinstance(x, str) for x in pp_mafft_extra
    ):
        errors.append("phylo.per_protein.mafft.extra_args must be a list of strings")
    if "domain_architecture" in pp_cfg and not isinstance(
        pp_cfg["domain_architecture"], bool
    ):
        errors.append("phylo.per_protein.domain_architecture must be a boolean")

    # trimAl blocks (whole-genome phylo.trimal + per-protein
    # phylo.per_protein.trimal share the same shape).
    _TRIMAL_MODES = (
        "automated1", "gappyout", "strict", "strictplus", "nogaps", "noallgaps",
    )

    def _check_trimal_block(block: dict, where: str) -> None:
        if "enabled" in block and not isinstance(block["enabled"], bool):
            errors.append(f"{where}.enabled must be a boolean")
        mode = block.get("mode", "automated1")
        if mode not in _TRIMAL_MODES:
            errors.append(
                f"{where}.mode '{mode}' is not supported "
                f"(use one of {list(_TRIMAL_MODES)}; threshold trimming goes "
                f"in {where}.extra_args)"
            )
        ta_extra = block.get("extra_args", [])
        if not isinstance(ta_extra, list) or not all(
            isinstance(x, str) for x in ta_extra
        ):
            errors.append(f"{where}.extra_args must be a list of strings")

    _check_trimal_block(phylo_cfg.get("trimal", {}) or {}, "phylo.trimal")
    _check_trimal_block(
        pp_cfg.get("trimal", {}) or {}, "phylo.per_protein.trimal",
    )

    part_cfg = phylo_cfg.get("partition", {}) or {}
    if "enabled" in part_cfg and not isinstance(part_cfg["enabled"], bool):
        errors.append("phylo.partition.enabled must be a boolean")
    plink = part_cfg.get("linkage", "proportional")
    if plink not in ("proportional", "equal", "unlinked"):
        errors.append(
            f"phylo.partition.linkage '{plink}' is not supported "
            "(use 'proportional', 'equal', or 'unlinked')"
        )
    pmodels = part_cfg.get("models", {})
    if not isinstance(pmodels, dict) or not all(
        isinstance(k, str) and isinstance(v, str) for k, v in pmodels.items()
    ):
        errors.append(
            "phylo.partition.models must be a mapping of family-label → "
            "model string"
        )

    rooting_cfg = phylo_cfg.get("rooting", {}) or {}
    rmethod = rooting_cfg.get("method", "auto")
    if rmethod not in ("auto", "taxonomy", "mad", "midpoint", "none"):
        errors.append(
            f"phylo.rooting.method '{rmethod}' is not supported "
            "(use 'auto', 'taxonomy', 'mad', 'midpoint', or 'none')"
        )

    lca_cfg = phylo_cfg.get("lca", {}) or {}
    if "enabled" in lca_cfg and not isinstance(lca_cfg["enabled"], bool):
        errors.append("phylo.lca.enabled must be a boolean")
    valid_min_ranks = {
        "none", "superkingdom", "realm", "kingdom", "subkingdom",
        "phylum", "subphylum", "class", "subclass",
        "order", "suborder", "family", "subfamily",
        "genus", "subgenus", "species",
    }
    mr = lca_cfg.get("min_rank", "genus")
    if mr not in valid_min_ranks:
        errors.append(
            f"phylo.lca.min_rank '{mr}' is not supported "
            f"(use one of {sorted(valid_min_ranks)})"
        )
    ctv = lca_cfg.get("coverage_threshold", 0.5)
    if not isinstance(ctv, (int, float)) or isinstance(ctv, bool) or not (0 <= ctv <= 1):
        errors.append(
            "phylo.lca.coverage_threshold must be a number between 0 and 1"
        )

    # HMM block
    hmm = cfg.get("hmm", {}) or {}
    if "enabled" in hmm and not isinstance(hmm["enabled"], bool):
        errors.append("hmm.enabled must be a boolean")
    if "use_ga_when_available" in hmm and not isinstance(
        hmm["use_ga_when_available"], bool
    ):
        errors.append("hmm.use_ga_when_available must be a boolean")
    db = hmm.get("database")
    if db is not None and not isinstance(db, str):
        errors.append(
            "hmm.database must be null (use bundled set) or a path string"
        )
    ev = hmm.get("default_evalue", 1.0e-5)
    if not isinstance(ev, (int, float)) or isinstance(ev, bool) or ev <= 0:
        errors.append(
            "hmm.default_evalue must be a positive number (e.g. 1.0e-5)"
        )
    rc = hmm.get("relative_length_cutoff", 0.5)
    if not isinstance(rc, (int, float)) or isinstance(rc, bool) or not (0 < rc <= 1):
        errors.append(
            "hmm.relative_length_cutoff must be a number in (0, 1] "
            "(fraction of HMM model length the alignment must span)"
        )
    ht = hmm.get("threads")
    if ht is not None and (
        not isinstance(ht, int) or isinstance(ht, bool) or ht < 1
    ):
        errors.append(
            "hmm.threads must be null (use cfg.threads) or a positive integer"
        )
    mot = hmm.get("multidomain_overlap_tolerance", 30)
    if not isinstance(mot, int) or isinstance(mot, bool) or mot < 0:
        errors.append(
            "hmm.multidomain_overlap_tolerance must be a non-negative integer "
            "(amino acids; 0 = strict non-overlap, 30 = default)"
        )

    # Representative priority
    priority = cfg.get("representative", {}).get("priority", [])
    valid_priorities = {"refseq", "reviewed_uniprot", "longest"}
    for p in priority:
        if p not in valid_priorities:
            errors.append(
                f"representative.priority contains unknown value '{p}'. "
                f"Valid values: {sorted(valid_priorities)}"
            )
    if "longest" not in priority:
        errors.append("representative.priority must include 'longest' as a fallback")

    pr_cfg = cfg.get("output", {}).get("protein_report", {}) or {}
    max_breakdown = pr_cfg.get("max_breakdown", 20)
    if (
        not isinstance(max_breakdown, int)
        or isinstance(max_breakdown, bool)
        or max_breakdown < 1
    ):
        errors.append(
            "output.protein_report.max_breakdown must be a positive integer"
        )

    return errors


def get_virus_config(cfg: dict[str, Any]) -> Optional[dict[str, Any]]:
    """Return the active virus config dict, or None if segmented mode is off."""
    seg = cfg.get("segmented", {})
    if not seg.get("enabled"):
        return None
    virus_name = seg.get("virus")
    return seg.get("viruses", {}).get(virus_name)
