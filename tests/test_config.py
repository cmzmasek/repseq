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


def test_validate_config_accepts_min_percent_zero_and_max_percent():
    cfg = load_config(None)
    cfg["qc"]["length_filter"]["min_percent"] = 0    # disables lower bound
    cfg["qc"]["length_filter"]["max_percent"] = 200  # optional upper cap
    assert validate_config(cfg) == []


def test_validate_config_rejects_max_percent_not_above_min():
    cfg = load_config(None)
    cfg["qc"]["length_filter"]["min_percent"] = 80
    cfg["qc"]["length_filter"]["max_percent"] = 50
    errors = validate_config(cfg)
    assert any("max_percent" in e for e in errors)


def test_validate_config_rejects_out_of_range_ambiguous():
    cfg = load_config(None)
    cfg["qc"]["ambiguous_threshold"] = 1.5
    errors = validate_config(cfg)
    assert any("ambiguous_threshold" in e for e in errors)


def test_validate_config_accepts_cdhit_backend():
    cfg = load_config(None)
    cfg["clustering"]["backend"] = "cdhit"
    assert validate_config(cfg) == []


def test_validate_config_rejects_unknown_backend():
    cfg = load_config(None)
    cfg["clustering"]["backend"] = "psi-cd-hit"
    errors = validate_config(cfg)
    assert any("clustering.backend" in e for e in errors)


def test_validate_config_rejects_cdhit_word_size_out_of_range():
    cfg = load_config(None)
    cfg["clustering"]["backend"] = "cdhit"
    cfg["clustering"]["cdhit"]["word_size"] = 12
    errors = validate_config(cfg)
    assert any("word_size" in e for e in errors)


def test_validate_config_rejects_cdhit_negative_memory():
    cfg = load_config(None)
    cfg["clustering"]["backend"] = "cdhit"
    cfg["clustering"]["cdhit"]["memory_mb"] = -1
    errors = validate_config(cfg)
    assert any("memory_mb" in e for e in errors)


def test_validate_config_rejects_cdhit_non_bool_global():
    cfg = load_config(None)
    cfg["clustering"]["backend"] = "cdhit"
    cfg["clustering"]["cdhit"]["global_alignment"] = "yes"
    errors = validate_config(cfg)
    assert any("global_alignment" in e for e in errors)


def test_validate_config_accepts_phylo_extra_args():
    cfg = load_config(None)
    cfg["phylo"]["mafft"]["extra_args"] = ["--maxiterate", "1000"]
    cfg["phylo"]["fasttree"]["extra_args"] = ["-fastest"]
    assert validate_config(cfg) == []


def test_validate_config_rejects_phylo_non_string_extra_args():
    cfg = load_config(None)
    cfg["phylo"]["mafft"]["extra_args"] = [123]
    errors = validate_config(cfg)
    assert any("phylo.mafft.extra_args" in e for e in errors)


def test_validate_config_default_use_genbank_metadata_is_true():
    cfg = load_config(None)
    assert cfg["segmented"]["use_genbank_metadata"] is True
    assert validate_config(cfg) == []


def test_validate_config_rejects_non_bool_use_genbank_metadata():
    cfg = load_config(None)
    cfg["segmented"]["use_genbank_metadata"] = "yes"
    errors = validate_config(cfg)
    assert any("use_genbank_metadata" in e for e in errors)


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


def test_validate_config_segment_aliases_accepted():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "fluA",
        "viruses": {
            "fluA": {
                "expected_segments": 2,
                "segments": ["HA", "NA"],
                "isolate_regex": r"(?P<isolate>X)",
                "segment_aliases": {
                    "HA": ["hemagglutinin", "haemagglutinin"],
                    "NA": ["neuraminidase"],
                },
            }
        },
    }
    assert validate_config(cfg) == []


def test_validate_config_segment_aliases_rejects_unknown_canonical():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "v",
        "viruses": {
            "v": {
                "expected_segments": 1,
                "segments": ["HA"],
                "isolate_regex": r"(?P<isolate>X)",
                "segment_aliases": {"NOPE": ["x"]},
            }
        },
    }
    errors = validate_config(cfg)
    assert any("segment_aliases" in e and "NOPE" in e for e in errors)


def test_validate_config_segment_aliases_rejects_empty_alias():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "v",
        "viruses": {
            "v": {
                "expected_segments": 1,
                "segments": ["HA"],
                "isolate_regex": r"(?P<isolate>X)",
                "segment_aliases": {"HA": ["valid", ""]},
            }
        },
    }
    errors = validate_config(cfg)
    assert any("segment_aliases.HA" in e for e in errors)


def test_validate_config_expected_proteins_per_segment_accepts_int_and_list():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "fluA",
        "viruses": {
            "fluA": {
                "expected_segments": 3,
                "segments": ["HA", "PB1", "NS"],
                "isolate_regex": r"(?P<isolate>X)",
                "expected_proteins_per_segment": {
                    "HA": 1,        # int — exact
                    "PB1": [1, 2],  # list — any of
                    "NS": [1, 2, 3],
                },
            }
        },
    }
    assert validate_config(cfg) == []


def test_validate_config_rejects_bad_expected_proteins_value():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "fluA",
        "viruses": {
            "fluA": {
                "expected_segments": 1,
                "segments": ["HA"],
                "isolate_regex": r"(?P<isolate>X)",
                "expected_proteins_per_segment": {"HA": "two"},
            }
        },
    }
    errors = validate_config(cfg)
    assert any("expected_proteins_per_segment.HA" in e for e in errors)


def test_validate_config_rejects_negative_in_list():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "fluA",
        "viruses": {
            "fluA": {
                "expected_segments": 1,
                "segments": ["HA"],
                "isolate_regex": r"(?P<isolate>X)",
                "expected_proteins_per_segment": {"HA": [1, -1]},
            }
        },
    }
    errors = validate_config(cfg)
    assert any("expected_proteins_per_segment.HA" in e for e in errors)


def test_validate_config_rejects_empty_list():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "fluA",
        "viruses": {
            "fluA": {
                "expected_segments": 1,
                "segments": ["HA"],
                "isolate_regex": r"(?P<isolate>X)",
                "expected_proteins_per_segment": {"HA": []},
            }
        },
    }
    errors = validate_config(cfg)
    assert any("expected_proteins_per_segment.HA" in e for e in errors)


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
