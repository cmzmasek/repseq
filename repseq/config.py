"""YAML configuration loading, validation, and defaults."""

from __future__ import annotations

import copy
import os
from pathlib import Path
from typing import Any, Optional

import yaml


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
        "length_filter": {
            "mode": "median_percent",   # "median_percent" | "min_max"
            "min_percent": 50,          # used when mode == median_percent; 0 disables lower bound
            "max_percent": None,        # optional upper cap, as % of median (median_percent mode)
            "min_length": None,         # used when mode == min_max
            "max_length": None,
        },
        "ambiguous_threshold": 0.05,
        "protein_annotation": {
            # Drop sequences whose NCBI GenBank record has fewer than
            # min_proteins annotated CDS features. Requires network access
            # (skipped automatically with --no-resolve).
            "enabled": False,
            "min_proteins": 1,
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
    },
    "clustering": {
        "backend": "mmseqs2",              # "mmseqs2" | "cdhit"
        # Alphabet fed to the clustering backend. "protein" (default) feeds
        # the marker protein (longest CDS, or the first matching alias in
        # `cluster_protein`) — for segmented input, an in-order concat of
        # each segment's marker. "nucleotide" feeds the raw sequence (the
        # pre-0.6 behaviour). "auto" picks protein when any sequence has
        # CDS info available, else nucleotide. Triggers a one-shot
        # GenBank CDS fetch when proteins aren't already cached.
        "alphabet": "protein",
        # Non-segmented marker-protein override. Alias list, matched
        # case-insensitively as substrings against /product. First alias
        # that matches a CDS wins; if no aliases match (or the list is
        # empty), the longest CDS on the sequence is used.
        "cluster_protein": [],
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
    },
    "taxonomy": {
        "ncbi_email": None,
        "ncbi_api_key": None,
        "cache_ttl_days": 30,
    },
    "output": {
        "dir": "./repseq_output",
        "prefix": "repseq",
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

def load_config(path: Optional[str | Path] = None) -> dict[str, Any]:
    """Load config from YAML file, merging over defaults."""
    cfg = copy.deepcopy(DEFAULTS)
    if path is not None:
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"Config file not found: {path}")
        with open(path) as fh:
            user_cfg = yaml.safe_load(fh) or {}
        cfg = _deep_merge(cfg, user_cfg)

    # Environment variable overrides
    if os.environ.get("REPSEQ_NCBI_EMAIL"):
        cfg["taxonomy"]["ncbi_email"] = os.environ["REPSEQ_NCBI_EMAIL"]
    if os.environ.get("REPSEQ_NCBI_API_KEY"):
        cfg["taxonomy"]["ncbi_api_key"] = os.environ["REPSEQ_NCBI_API_KEY"]

    cfg = _expand_paths(cfg)
    return cfg


def validate_config(cfg: dict[str, Any]) -> list[str]:
    """Return a list of validation error messages (empty = valid)."""
    errors: list[str] = []

    # Length filter
    lf = cfg.get("qc", {}).get("length_filter", {})
    mode = lf.get("mode")
    if mode not in ("median_percent", "min_max"):
        errors.append(
            f"qc.length_filter.mode must be 'median_percent' or 'min_max', got '{mode}'"
        )
    if mode == "median_percent":
        pct = lf.get("min_percent")
        if not isinstance(pct, (int, float)) or not (0 <= pct <= 100):
            errors.append(
                "qc.length_filter.min_percent must be a number between 0 and 100 "
                "(0 disables the lower bound)"
            )
        max_pct = lf.get("max_percent")
        if max_pct is not None:
            if not isinstance(max_pct, (int, float)) or max_pct <= 0:
                errors.append(
                    "qc.length_filter.max_percent must be a positive number "
                    "(percent of median length)"
                )
            elif isinstance(pct, (int, float)) and max_pct <= pct:
                errors.append(
                    "qc.length_filter.max_percent must be greater than min_percent"
                )
    if mode == "min_max":
        mn = lf.get("min_length")
        mx = lf.get("max_length")
        if mn is not None and mx is not None and mn >= mx:
            errors.append("qc.length_filter.min_length must be less than max_length")

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
                            f"must be a mapping of segment-name → list of strings"
                        )
                    else:
                        seg_names = set(vdef.get("segments", []))
                        for seg_name, aliases in cp.items():
                            if seg_name not in seg_names:
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"cluster_protein: unknown segment '{seg_name}'"
                                )
                            if not isinstance(aliases, list) or not all(
                                isinstance(s, str) and s.strip() for s in aliases
                            ):
                                errors.append(
                                    f"segmented.viruses.{virus_name}."
                                    f"cluster_protein.{seg_name} "
                                    f"must be a list of non-empty strings"
                                )

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

    alphabet = cfg.get("clustering", {}).get("alphabet", "protein")
    if alphabet not in ("protein", "nucleotide", "auto"):
        errors.append(
            f"clustering.alphabet '{alphabet}' is not supported "
            f"(use 'protein', 'nucleotide', or 'auto')"
        )

    cluster_protein_global = cfg.get("clustering", {}).get("cluster_protein", [])
    if not isinstance(cluster_protein_global, list) or not all(
        isinstance(s, str) and s.strip() for s in cluster_protein_global
    ):
        errors.append(
            "clustering.cluster_protein must be a list of non-empty strings "
            "(case-insensitive substring aliases for /product)"
        )

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

    return errors


def get_virus_config(cfg: dict[str, Any]) -> Optional[dict[str, Any]]:
    """Return the active virus config dict, or None if segmented mode is off."""
    seg = cfg.get("segmented", {})
    if not seg.get("enabled"):
        return None
    virus_name = seg.get("virus")
    return seg.get("viruses", {}).get(virus_name)
