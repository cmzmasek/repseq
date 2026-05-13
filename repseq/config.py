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
            "min_percent": 50,          # used when mode == median_percent
            "min_length": None,         # used when mode == min_max
            "max_length": None,
        },
        "ambiguous_threshold": 0.05,
        "annotation_filter": {
            "enabled": True,
            "keywords": [
                "MAG",
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
        "backend": "mmseqs2",
        "mmseqs2_mode": "easy-linclust",   # "easy-linclust" | "easy-cluster"
        "coverage": 0.8,
        "coverage_mode": 0,
        "extra_args": [],
    },
    "representative": {
        "priority": ["refseq", "reviewed_uniprot", "longest"],
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
        if not isinstance(pct, (int, float)) or not (0 < pct < 100):
            errors.append(
                "qc.length_filter.min_percent must be a number between 0 and 100"
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

    # Clustering backend
    backend = cfg.get("clustering", {}).get("backend")
    if backend not in ("mmseqs2",):
        errors.append(f"clustering.backend '{backend}' is not supported (use 'mmseqs2')")

    mmseqs2_mode = cfg.get("clustering", {}).get("mmseqs2_mode")
    if mmseqs2_mode not in ("easy-linclust", "easy-cluster"):
        errors.append(
            f"clustering.mmseqs2_mode must be 'easy-linclust' or 'easy-cluster', got '{mmseqs2_mode}'"
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
