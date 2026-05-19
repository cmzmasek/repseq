"""Config loading, defaults merging, validation."""
from __future__ import annotations

from pathlib import Path

import yaml

from repseq.config import DEFAULTS, get_virus_config, load_config, validate_config


def test_load_config_returns_defaults_when_no_file():
    cfg = load_config(None)
    assert cfg["threads"] == DEFAULTS["threads"]
    assert cfg["qc"]["genome_length_filter"]["enabled"] is False


def test_load_config_merges_user_overrides(tmp_path: Path):
    user = {"threads": 16, "qc": {"ambiguous_threshold": 0.1}}
    cfg_path = tmp_path / "c.yaml"
    cfg_path.write_text(yaml.dump(user))

    cfg = load_config(cfg_path)
    assert cfg["threads"] == 16
    # Overridden inside a nested dict
    assert cfg["qc"]["ambiguous_threshold"] == 0.1
    # Untouched defaults survive the merge
    assert cfg["qc"]["genome_length_filter"]["enabled"] is False


def test_load_config_env_var_override(tmp_path: Path, monkeypatch):
    monkeypatch.setenv("REPSEQ_NCBI_EMAIL", "test@example.com")
    cfg = load_config(None)
    assert cfg["taxonomy"]["ncbi_email"] == "test@example.com"


def test_validate_config_accepts_defaults():
    errors = validate_config(load_config(None))
    assert errors == []


def test_validate_config_rejects_legacy_length_filter_key():
    """The renamed qc.length_filter must be rejected so configs migrate
    consciously instead of silently losing length filtering."""
    cfg = load_config(None)
    cfg["qc"]["length_filter"] = {"mode": "median_percent", "min_percent": 50}
    errors = validate_config(cfg)
    assert any("qc.length_filter was renamed" in e for e in errors)


def test_validate_config_accepts_genome_length_filter_bounds():
    cfg = load_config(None)
    cfg["qc"]["genome_length_filter"] = {"enabled": True, "min": 9000, "max": 13000}
    assert validate_config(cfg) == []


def test_validate_config_genome_length_filter_disabled_ignores_bounds():
    cfg = load_config(None)
    # Disabled → invalid bounds are not checked (the filter won't run).
    cfg["qc"]["genome_length_filter"] = {"enabled": False, "min": 13000, "max": 9000}
    assert validate_config(cfg) == []


def test_validate_config_rejects_genome_length_filter_with_segmented():
    cfg = load_config(None)
    cfg["qc"]["genome_length_filter"] = {"enabled": True, "min": 9000}
    cfg["segmented"]["enabled"] = True
    errors = validate_config(cfg)
    assert any("cannot be true when" in e for e in errors)


def test_validate_config_rejects_genome_length_filter_enabled_without_bounds():
    cfg = load_config(None)
    cfg["qc"]["genome_length_filter"] = {"enabled": True, "min": None, "max": None}
    errors = validate_config(cfg)
    assert any("neither min nor max" in e for e in errors)


def test_validate_config_rejects_genome_length_filter_min_above_max():
    cfg = load_config(None)
    cfg["qc"]["genome_length_filter"] = {"enabled": True, "min": 13000, "max": 9000}
    errors = validate_config(cfg)
    assert any("min must be <=" in e for e in errors)


def test_validate_config_rejects_out_of_range_ambiguous():
    cfg = load_config(None)
    cfg["qc"]["ambiguous_threshold"] = 1.5
    errors = validate_config(cfg)
    assert any("ambiguous_threshold" in e for e in errors)


def test_validate_config_rejects_out_of_range_protein_quality():
    cfg = load_config(None)
    cfg["qc"]["protein_quality"]["enabled"] = True
    cfg["qc"]["protein_quality"]["max_bad_fraction"] = 1.5
    errors = validate_config(cfg)
    assert any("protein_quality.max_bad_fraction" in e for e in errors)


def test_validate_config_accepts_valid_protein_quality():
    cfg = load_config(None)
    cfg["qc"]["protein_quality"]["enabled"] = True
    cfg["qc"]["protein_quality"]["max_bad_fraction"] = 0.1
    assert validate_config(cfg) == []


def test_validate_config_ignores_protein_quality_when_disabled():
    cfg = load_config(None)
    cfg["qc"]["protein_quality"]["enabled"] = False
    cfg["qc"]["protein_quality"]["max_bad_fraction"] = 9.0  # invalid but unused
    assert validate_config(cfg) == []


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
    cfg["phylo"]["iqtree"]["extra_args"] = ["-alrt", "1000"]
    assert validate_config(cfg) == []


def test_validate_config_rejects_phylo_non_string_extra_args():
    cfg = load_config(None)
    cfg["phylo"]["mafft"]["extra_args"] = [123]
    errors = validate_config(cfg)
    assert any("phylo.mafft.extra_args" in e for e in errors)


def test_validate_config_default_phylo_tool_is_auto():
    cfg = load_config(None)
    assert cfg["phylo"]["tool"] == "auto"
    assert cfg["phylo"]["iqtree"]["model"] == "MFP"
    assert cfg["phylo"]["iqtree"]["ultrafast_bootstrap"] == 1000
    assert validate_config(cfg) == []


def test_validate_config_accepts_phylo_tool_iqtree_or_fasttree():
    cfg = load_config(None)
    for tool in ("iqtree", "fasttree", "auto"):
        cfg["phylo"]["tool"] = tool
        assert validate_config(cfg) == []


def test_validate_config_rejects_unknown_phylo_tool():
    cfg = load_config(None)
    cfg["phylo"]["tool"] = "raxml"
    errors = validate_config(cfg)
    assert any("phylo.tool" in e for e in errors)


def test_validate_config_rejects_negative_ufboot():
    cfg = load_config(None)
    cfg["phylo"]["iqtree"]["ultrafast_bootstrap"] = -1
    errors = validate_config(cfg)
    assert any("ultrafast_bootstrap" in e for e in errors)


def test_validate_config_rejects_non_int_ufboot():
    cfg = load_config(None)
    cfg["phylo"]["iqtree"]["ultrafast_bootstrap"] = "1000"
    errors = validate_config(cfg)
    assert any("ultrafast_bootstrap" in e for e in errors)


def test_validate_config_rejects_non_string_iqtree_model():
    cfg = load_config(None)
    cfg["phylo"]["iqtree"]["model"] = 42
    errors = validate_config(cfg)
    assert any("phylo.iqtree.model" in e for e in errors)


def test_validate_config_rejects_non_string_iqtree_binary():
    cfg = load_config(None)
    cfg["phylo"]["iqtree"]["binary"] = 7
    errors = validate_config(cfg)
    assert any("phylo.iqtree.binary" in e for e in errors)


# ---------------------------------------------------------------------------
# phylo.labeling and phylo.phyloxml (Pass A — rich phyloXML annotation)
# ---------------------------------------------------------------------------

def test_validate_config_default_phylo_labeling():
    cfg = load_config(None)
    assert cfg["phylo"]["labeling"]["format"] == "{species}|{id}|{host}"
    assert cfg["phylo"]["labeling"]["segmented_format"] == "{species}|{strain}|{host}"
    assert cfg["phylo"]["labeling"]["replace_whitespace"] is True
    assert cfg["phylo"]["labeling"]["keep_separator_on_empty"] is False
    assert validate_config(cfg) == []


def test_validate_config_rejects_non_string_label_format():
    cfg = load_config(None)
    cfg["phylo"]["labeling"]["format"] = 42
    errors = validate_config(cfg)
    assert any("phylo.labeling.format" in e for e in errors)


def test_validate_config_accepts_null_segmented_format():
    """Null segmented_format means: fall back to the regular format."""
    cfg = load_config(None)
    cfg["phylo"]["labeling"]["segmented_format"] = None
    assert validate_config(cfg) == []


def test_validate_config_rejects_non_bool_replace_whitespace():
    cfg = load_config(None)
    cfg["phylo"]["labeling"]["replace_whitespace"] = "yes"
    errors = validate_config(cfg)
    assert any("phylo.labeling.replace_whitespace" in e for e in errors)


def test_validate_config_default_phyloxml():
    cfg = load_config(None)
    assert cfg["phylo"]["phyloxml"]["confidence_type"] == "auto"
    assert validate_config(cfg) == []


def test_validate_config_accepts_all_confidence_types():
    cfg = load_config(None)
    for ct in ("auto", "sh_like", "sh_alrt", "ufboot", "bootstrap"):
        cfg["phylo"]["phyloxml"]["confidence_type"] = ct
        assert validate_config(cfg) == [], f"rejected {ct}"


def test_validate_config_rejects_unknown_confidence_type():
    cfg = load_config(None)
    cfg["phylo"]["phyloxml"]["confidence_type"] = "raxml-bootstrap"
    errors = validate_config(cfg)
    assert any("confidence_type" in e for e in errors)


# ---------------------------------------------------------------------------
# phylo.rooting and phylo.lca (Pass B — rooting + internal LCA)
# ---------------------------------------------------------------------------

def test_validate_config_default_rooting():
    cfg = load_config(None)
    assert cfg["phylo"]["rooting"]["method"] == "auto"
    assert validate_config(cfg) == []


def test_validate_config_accepts_all_rooting_methods():
    cfg = load_config(None)
    for method in ("auto", "taxonomy", "mad", "midpoint", "none"):
        cfg["phylo"]["rooting"]["method"] = method
        assert validate_config(cfg) == [], f"rejected {method}"


def test_validate_config_rejects_unknown_rooting_method():
    cfg = load_config(None)
    cfg["phylo"]["rooting"]["method"] = "raxml-root"
    errors = validate_config(cfg)
    assert any("phylo.rooting.method" in e for e in errors)


def test_validate_config_default_lca():
    cfg = load_config(None)
    assert cfg["phylo"]["lca"]["enabled"] is True
    assert cfg["phylo"]["lca"]["min_rank"] == "genus"
    assert cfg["phylo"]["lca"]["coverage_threshold"] == 0.5
    assert validate_config(cfg) == []


def test_validate_config_rejects_non_bool_lca_enabled():
    cfg = load_config(None)
    cfg["phylo"]["lca"]["enabled"] = "yes"
    errors = validate_config(cfg)
    assert any("phylo.lca.enabled" in e for e in errors)


def test_validate_config_accepts_lca_min_rank_options():
    cfg = load_config(None)
    for r in ("none", "family", "genus", "species", "order"):
        cfg["phylo"]["lca"]["min_rank"] = r
        assert validate_config(cfg) == [], f"rejected {r}"


def test_validate_config_rejects_unknown_min_rank():
    cfg = load_config(None)
    cfg["phylo"]["lca"]["min_rank"] = "isolate"
    errors = validate_config(cfg)
    assert any("phylo.lca.min_rank" in e for e in errors)


def test_validate_config_rejects_out_of_range_coverage_threshold():
    cfg = load_config(None)
    cfg["phylo"]["lca"]["coverage_threshold"] = 1.5
    errors = validate_config(cfg)
    assert any("coverage_threshold" in e for e in errors)
    cfg["phylo"]["lca"]["coverage_threshold"] = -0.1
    errors = validate_config(cfg)
    assert any("coverage_threshold" in e for e in errors)


def test_validate_config_default_alphabet_is_protein():
    cfg = load_config(None)
    assert cfg["clustering"]["alphabet_for_clustering"] == "protein"
    assert validate_config(cfg) == []


def test_validate_config_accepts_alphabet_nucleotide():
    cfg = load_config(None)
    cfg["clustering"]["alphabet_for_clustering"] = "nucleotide"
    assert validate_config(cfg) == []


def test_validate_config_rejects_unknown_alphabet():
    cfg = load_config(None)
    cfg["clustering"]["alphabet_for_clustering"] = "dna"
    errors = validate_config(cfg)
    assert any("clustering.alphabet_for_clustering" in e for e in errors)


def test_validate_config_rejects_auto_alphabet():
    """v0.10.0 dropped the 'auto' alphabet — it's no longer accepted."""
    cfg = load_config(None)
    cfg["clustering"]["alphabet_for_clustering"] = "auto"
    errors = validate_config(cfg)
    assert any("clustering.alphabet_for_clustering" in e for e in errors)


def test_validate_config_rejects_non_list_cluster_protein():
    cfg = load_config(None)
    cfg["clustering"]["cluster_protein"] = "polymerase"  # string, not list
    errors = validate_config(cfg)
    assert any("clustering.cluster_protein" in e for e in errors)


def test_validate_config_segmented_cluster_protein_accepted():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "v",
        "viruses": {
            "v": {
                "expected_segments": 2,
                "segments": ["L", "S"],
                "isolate_regex": r"(?P<isolate>X)",
                "cluster_protein": {
                    "L": ["polymerase", "L protein"],
                    "S": ["nucleoprotein"],
                },
            }
        },
    }
    assert validate_config(cfg) == []


def test_validate_config_segmented_cluster_protein_rejects_unknown_segment():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "v",
        "viruses": {
            "v": {
                "expected_segments": 1,
                "segments": ["L"],
                "isolate_regex": r"(?P<isolate>X)",
                "cluster_protein": {"NOPE": ["x"]},
            }
        },
    }
    errors = validate_config(cfg)
    assert any("cluster_protein" in e and "NOPE" in e for e in errors)


def test_validate_config_segmented_cluster_protein_rejects_empty_alias():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True,
        "virus": "v",
        "viruses": {
            "v": {
                "expected_segments": 1,
                "segments": ["L"],
                "isolate_regex": r"(?P<isolate>X)",
                "cluster_protein": {"L": ["valid", ""]},
            }
        },
    }
    errors = validate_config(cfg)
    assert any("cluster_protein.L" in e for e in errors)


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


# ---------------------------------------------------------------------------
# v0.13: dict-form cluster_protein + segment_markers + hmm block
# ---------------------------------------------------------------------------

def test_validate_config_accepts_dict_form_cluster_protein():
    cfg = load_config(None)
    cfg["clustering"]["cluster_protein"] = [
        "polymerase",                                          # legacy string
        {"name": "Spike", "aliases": ["spike"], "hmms": ["Corona_S1", "Corona_S2"]},
        {"name": "N", "hmms": ["CoV_nucleocap"]},              # hmms only
        {"name": "M", "aliases": ["membrane"]},                # aliases only
    ]
    assert validate_config(cfg) == []


def test_validate_config_rejects_dict_form_with_no_aliases_or_hmms():
    cfg = load_config(None)
    cfg["clustering"]["cluster_protein"] = [{"name": "empty"}]
    errs = validate_config(cfg)
    assert any("at least one of 'aliases' or 'hmms'" in e for e in errs)


def test_validate_config_rejects_dict_form_missing_name():
    cfg = load_config(None)
    cfg["clustering"]["cluster_protein"] = [{"aliases": ["x"]}]
    errs = validate_config(cfg)
    assert any("non-empty 'name'" in e for e in errs)


def test_validate_config_accepts_segment_markers():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True, "virus": "fluA",
        "viruses": {"fluA": {
            "expected_segments": 1, "segments": ["HA"],
            "isolate_regex": r"(?P<isolate>X)",
            "segment_markers": {
                "HA": {"aliases": ["hemagglutinin"], "hmms": ["Hemagglutinin"]},
            },
        }},
    }
    assert validate_config(cfg) == []


def test_validate_config_rejects_segment_markers_unknown_segment():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True, "virus": "fluA",
        "viruses": {"fluA": {
            "expected_segments": 1, "segments": ["HA"],
            "isolate_regex": r"(?P<isolate>X)",
            "segment_markers": {"NOPE": {"hmms": ["X"]}},
        }},
    }
    errs = validate_config(cfg)
    assert any("unknown segment 'NOPE'" in e for e in errs)


def test_validate_config_hmm_block_defaults_validate():
    cfg = load_config(None)
    assert validate_config(cfg) == []
    assert cfg["hmm"]["enabled"] is True


def test_validate_config_rejects_hmm_default_evalue_zero():
    cfg = load_config(None)
    cfg["hmm"]["default_evalue"] = 0
    assert any(
        "hmm.default_evalue" in e for e in validate_config(cfg)
    )


def test_validate_config_rejects_hmm_relative_length_out_of_range():
    cfg = load_config(None)
    cfg["hmm"]["relative_length_cutoff"] = 1.5
    assert any(
        "hmm.relative_length_cutoff" in e for e in validate_config(cfg)
    )


def test_validate_config_accepts_user_hmm_database_path():
    """validate_config does NOT check the path exists — that's runtime."""
    cfg = load_config(None)
    cfg["hmm"]["database"] = "/path/to/my.hmm"
    assert validate_config(cfg) == []


def test_validate_config_accepts_extra_segments_action_default():
    cfg = load_config(None)
    assert validate_config(cfg) == []
    assert cfg["segmented"]["extra_segments_action"] == "warn"


def test_validate_config_accepts_extra_segments_action_drop():
    cfg = load_config(None)
    cfg["segmented"]["extra_segments_action"] = "drop"
    assert validate_config(cfg) == []


def test_validate_config_rejects_unknown_extra_segments_action():
    cfg = load_config(None)
    cfg["segmented"]["extra_segments_action"] = "ignore"
    errs = validate_config(cfg)
    assert any("extra_segments_action" in e for e in errs)
