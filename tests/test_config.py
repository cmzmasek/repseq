"""Config loading, defaults merging, validation."""
from __future__ import annotations

from pathlib import Path

import yaml

from repseq.config import DEFAULTS, get_virus_config, load_config, validate_config


def test_load_config_returns_defaults_when_no_file():
    cfg = load_config(None)
    assert cfg["threads"] == DEFAULTS["threads"]
    assert cfg["qc"]["length_filter"]["mode"] == "median_percent"


def test_load_config_merges_user_overrides(tmp_path: Path):
    user = {"threads": 16, "qc": {"ambiguous_threshold": 0.1}}
    cfg_path = tmp_path / "c.yaml"
    cfg_path.write_text(yaml.dump(user))

    cfg = load_config(cfg_path)
    assert cfg["threads"] == 16
    # Overridden inside a nested dict
    assert cfg["qc"]["ambiguous_threshold"] == 0.1
    # Untouched defaults survive the merge
    assert cfg["qc"]["length_filter"]["mode"] == "median_percent"


def test_load_config_env_var_override(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("REPSEQ_NCBI_EMAIL", "test@example.com")
    cfg = load_config(None)
    assert cfg["taxonomy"]["ncbi_email"] == "test@example.com"


def test_validate_config_accepts_defaults():
    errors = validate_config(load_config(None))
    assert errors == []


def test_validate_config_rejects_bad_length_mode():
    cfg = load_config(None)
    cfg["qc"]["length_filter"]["mode"] = "bogus"
    errors = validate_config(cfg)
    assert any("length_filter.mode" in e for e in errors)


def test_validate_config_rejects_out_of_range_ambiguous():
    cfg = load_config(None)
    cfg["qc"]["ambiguous_threshold"] = 1.5
    errors = validate_config(cfg)
    assert any("ambiguous_threshold" in e for e in errors)


def test_validate_config_requires_longest_in_priority():
    cfg = load_config(None)
    cfg["representative"]["priority"] = ["refseq", "reviewed_uniprot"]
    errors = validate_config(cfg)
    assert any("longest" in e for e in errors)


def test_validate_config_segmented_requires_virus_name():
    cfg = load_config(None)
    cfg["segmented"]["enabled"] = True
    cfg["segmented"]["virus"] = None
    errors = validate_config(cfg)
    assert any("segmented.virus" in e for e in errors)


def test_validate_config_segmented_requires_fields_for_named_virus():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "fluA",
        "viruses": {"fluA": {"segments": ["HA", "NA"]}},  # missing isolate_regex/expected_segments
    }
    errors = validate_config(cfg)
    # Two missing required fields => two errors
    assert sum("missing required field" in e for e in errors) == 2


def test_get_virus_config_returns_none_when_disabled():
    cfg = load_config(None)
    assert get_virus_config(cfg) is None


def test_get_virus_config_returns_named_block():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "fluA",
        "viruses": {
            "fluA": {
                "expected_segments": 2,
                "segments": ["HA", "NA"],
                "isolate_regex": r"(?P<isolate>X)",
            }
        },
    }
    v = get_virus_config(cfg)
    assert v["segments"] == ["HA", "NA"]
