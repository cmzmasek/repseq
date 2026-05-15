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
    },
    "clustering": {
        "backend": "mmseqs2",              # "mmseqs2" | "cdhit"
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
    for tool in ("mafft", "fasttree"):
        tool_cfg = phylo_cfg.get(tool, {}) or {}
        extra = tool_cfg.get("extra_args", [])
        if not isinstance(extra, list) or not all(isinstance(x, str) for x in extra):
            errors.append(f"phylo.{tool}.extra_args must be a list of strings")

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
