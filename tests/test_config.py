"""Config loading, defaults merging, validation."""
from __future__ import annotations

from pathlib import Path

import yaml

from repseq import __version__ as REPSEQ_VERSION
from repseq.config import (
    DEFAULTS,
    effective_config_filename,
    get_virus_config,
    load_config,
    sanitize_config,
    validate_config,
    write_effective_config,
)


def test_sanitize_config_blanks_secrets_keeps_key():
    cfg = {"taxonomy": {"ncbi_email": "x@y.com", "ncbi_api_key": "SECRET",
                        "rate_limit": 3}}
    out = sanitize_config(cfg)
    assert out["taxonomy"]["ncbi_email"] is None
    assert out["taxonomy"]["ncbi_api_key"] is None
    assert out["taxonomy"]["rate_limit"] == 3  # non-secret untouched


def test_sanitize_config_drops_private_keys_by_default():
    cfg = {"a": 1, "_hmm_runtime": {"active": True},
           "_taxonomy_review": {"verdicts": [1]}}
    assert sanitize_config(cfg) == {"a": 1}


def test_sanitize_config_keeps_private_when_requested():
    cfg = {"a": 1, "_hmm_runtime": {"active": True},
           "taxonomy": {"ncbi_api_key": "SECRET"}}
    out = sanitize_config(cfg, drop_private=False)
    assert "_hmm_runtime" in out                       # kept
    assert out["taxonomy"]["ncbi_api_key"] is None      # secret still blanked


def test_sanitize_config_is_non_mutating():
    cfg = {"taxonomy": {"ncbi_api_key": "SECRET"}, "_x": 1}
    sanitize_config(cfg)
    assert cfg["taxonomy"]["ncbi_api_key"] == "SECRET"  # original intact
    assert cfg["_x"] == 1


def test_effective_config_filename_uses_underscored_version():
    name = effective_config_filename("cov")
    assert name == f"cov_config_repseq{REPSEQ_VERSION.replace('.', '_')}.yaml"


def test_write_effective_config_no_comments_no_secrets_reloadable(tmp_path: Path):
    cfg = load_config(None)
    cfg["taxonomy"]["ncbi_email"] = "x@y.com"
    cfg["taxonomy"]["ncbi_api_key"] = "SECRETKEY123"
    cfg["_hmm_runtime"] = {"active": True, "ga_cutoffs": {"RdRp": 25.0}}
    path = tmp_path / effective_config_filename("test")
    write_effective_config(cfg, path)
    text = path.read_text()
    assert "SECRETKEY123" not in text and "x@y.com" not in text
    assert "_hmm_runtime" not in text
    assert "ncbi_api_key: null" in text   # key kept, value gone
    # No YAML comment lines (a '#' inside a quoted hex colour value is fine).
    assert not any(ln.lstrip().startswith("#") for ln in text.splitlines())
    # Fully resolved (all settings present) AND re-loadable as a config.
    reloaded = load_config(path)
    assert set(reloaded) >= set(DEFAULTS)
    assert reloaded["taxonomy"]["ncbi_api_key"] is None


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


def test_validate_config_per_protein_mafft_default_is_auto():
    cfg = load_config(None)
    # Default is empty → MAFFT --auto (fast); L-INS-i is opt-in.
    assert cfg["phylo"]["per_protein"]["mafft"]["extra_args"] == []
    assert validate_config(cfg) == []


def test_validate_config_rejects_bad_per_protein_mafft_args():
    cfg = load_config(None)
    cfg["phylo"]["per_protein"]["mafft"]["extra_args"] = "--localpair"
    errors = validate_config(cfg)
    assert any("per_protein.mafft.extra_args" in e for e in errors)


def test_validate_config_rejects_bad_per_protein_incongruence():
    cfg = load_config(None)
    cfg["phylo"]["per_protein"]["incongruence"] = "yes"
    errors = validate_config(cfg)
    assert any("per_protein.incongruence" in e for e in errors)


def test_validate_config_per_protein_domain_architecture_default_on():
    cfg = load_config(None)
    assert cfg["phylo"]["per_protein"]["domain_architecture"] is True
    assert validate_config(cfg) == []


def test_validate_config_rejects_bad_per_protein_domain_architecture():
    cfg = load_config(None)
    cfg["phylo"]["per_protein"]["domain_architecture"] = "yes"
    errors = validate_config(cfg)
    assert any("per_protein.domain_architecture" in e for e in errors)


def test_validate_config_whole_polyprotein_tree_default_off():
    cfg = load_config(None)
    assert cfg["phylo"]["per_protein"]["whole_polyprotein_tree"] is False
    assert validate_config(cfg) == []


def test_validate_config_rejects_bad_whole_polyprotein_tree():
    cfg = load_config(None)
    cfg["phylo"]["per_protein"]["whole_polyprotein_tree"] = "yes"
    errors = validate_config(cfg)
    assert any("per_protein.whole_polyprotein_tree" in e for e in errors)


def test_validate_config_trimal_default_off_and_valid():
    cfg = load_config(None)
    assert cfg["phylo"]["trimal"]["enabled"] is False
    assert cfg["phylo"]["trimal"]["mode"] == "automated1"
    assert cfg["phylo"]["per_protein"]["trimal"]["enabled"] is False
    assert validate_config(cfg) == []


def test_validate_config_rejects_bad_trimal_mode():
    cfg = load_config(None)
    cfg["phylo"]["trimal"]["mode"] = "supertrim"
    errors = validate_config(cfg)
    assert any("phylo.trimal.mode" in e for e in errors)


def test_validate_config_rejects_bad_trimal_extra_args():
    cfg = load_config(None)
    cfg["phylo"]["trimal"]["extra_args"] = "-gt 0.8"
    errors = validate_config(cfg)
    assert any("phylo.trimal.extra_args" in e for e in errors)


def test_validate_config_validates_per_protein_trimal():
    cfg = load_config(None)
    cfg["phylo"]["per_protein"]["trimal"]["mode"] = "nope"
    errors = validate_config(cfg)
    assert any("phylo.per_protein.trimal.mode" in e for e in errors)


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


def test_validate_config_rejects_phylo_mafft_use_auto_non_bool():
    cfg = load_config(None)
    cfg["phylo"]["mafft"]["use_auto"] = "yes"
    errors = validate_config(cfg)
    assert any("phylo.mafft.use_auto" in e for e in errors)


def test_validate_config_accepts_phylo_mafft_use_auto_bool():
    cfg = load_config(None)
    for v in (True, False):
        cfg["phylo"]["mafft"]["use_auto"] = v
        assert validate_config(cfg) == []


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


def test_validate_config_accepts_outgroup_rooting_with_accession():
    cfg = load_config(None)
    cfg["phylo"]["rooting"]["method"] = "outgroup"
    cfg["phylo"]["rooting"]["outgroup"] = "AB123456"
    assert validate_config(cfg) == []


def test_validate_config_accepts_outgroup_rooting_with_clade_list():
    cfg = load_config(None)
    cfg["phylo"]["rooting"]["method"] = "outgroup"
    cfg["phylo"]["rooting"]["outgroup"] = ["AB1", "AB2"]
    assert validate_config(cfg) == []


def test_validate_config_accepts_outgroup_rooting_with_rank():
    cfg = load_config(None)
    cfg["phylo"]["rooting"]["method"] = "outgroup"
    cfg["phylo"]["rooting"]["outgroup_rank"] = {"family": "Hantaviridae"}
    assert validate_config(cfg) == []


def test_validate_config_rejects_outgroup_method_without_target():
    cfg = load_config(None)
    cfg["phylo"]["rooting"]["method"] = "outgroup"
    errors = validate_config(cfg)
    assert any(
        "phylo.rooting.method='outgroup'" in e for e in errors
    )


def test_validate_config_rejects_outgroup_non_string():
    cfg = load_config(None)
    cfg["phylo"]["rooting"]["outgroup"] = 42
    errors = validate_config(cfg)
    assert any("phylo.rooting.outgroup" in e for e in errors)


def test_validate_config_rejects_outgroup_rank_non_mapping():
    cfg = load_config(None)
    cfg["phylo"]["rooting"]["outgroup_rank"] = ["family"]
    errors = validate_config(cfg)
    assert any("phylo.rooting.outgroup_rank" in e for e in errors)


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


# ---------------------------------------------------------------------------
# Unknown-key audit (a typo'd / misplaced key is a hard error)
# ---------------------------------------------------------------------------

def test_unknown_top_level_key_rejected():
    cfg = load_config(None)
    cfg["foo"] = [{"bar": True}]
    errs = validate_config(cfg)
    assert any("unknown config key 'foo'" in e for e in errs)


def test_unknown_section_key_rejected_with_suggestion():
    cfg = load_config(None)
    cfg["phylo"]["per_protei"] = {}            # typo for per_protein
    errs = validate_config(cfg)
    assert any(
        "unknown config key 'phylo.per_protei'" in e and "per_protein" in e
        for e in errs
    )


def test_unknown_nested_closed_schema_key_rejected():
    cfg = load_config(None)
    cfg["hmm"]["default_eval"] = 1.0e-5        # typo for default_evalue
    errs = validate_config(cfg)
    assert any("unknown config key 'hmm.default_eval'" in e for e in errs)
    assert any("default_evalue" in e for e in errs)


def test_unknown_key_on_polyprotein_spec_rejected():
    """The motivating bug: whole_polyprotein_tree belongs under
    phylo.per_protein, not on a clustering.polyprotein spec."""
    cfg = load_config(None)
    cfg["clustering"]["polyprotein"] = [{
        "name": "ORF1ab",
        "whole_polyprotein_tree": True,        # silently ignored before
        "peptides": [
            {"name": "NSP1", "hmm": "A"},
            {"name": "NSP2", "hmm": "B"},
        ],
    }]
    errs = validate_config(cfg)
    assert any(
        "clustering.polyprotein[0]" in e
        and "unknown key 'whole_polyprotein_tree'" in e
        for e in errs
    )


def test_unknown_key_on_marker_spec_rejected():
    cfg = load_config(None)
    cfg["clustering"]["cluster_protein"] = [
        {"name": "Spike", "hmms": ["CoV_S1"], "hmm": ["typo"]},  # 'hmm' singular
    ]
    errs = validate_config(cfg)
    assert any(
        "clustering.cluster_protein[0]" in e and "unknown key 'hmm'" in e
        for e in errs
    )


def test_unknown_key_on_peptide_rejected():
    cfg = load_config(None)
    cfg["clustering"]["polyprotein"] = [{
        "name": "ORF1ab",
        "peptides": [
            {"name": "NSP1", "hmm": "A", "motif": "GG"},  # typo for cleavage_motif
            {"name": "NSP2", "hmm": "B"},
        ],
    }]
    errs = validate_config(cfg)
    assert any(
        "peptides[0]" in e and "unknown key 'motif'" in e
        and "cleavage_motif" in e
        for e in errs
    )


def test_unknown_key_on_virus_block_rejected():
    cfg = load_config(None)
    cfg["segmented"] = {
        "enabled": True, "virus": "fluA",
        "viruses": {"fluA": {
            "expected_segments": 1, "segments": ["HA"],
            "isolate_regex": r"(?P<isolate>X)",
            "segment_legnths": {},             # typo for segment_lengths
        }},
    }
    errs = validate_config(cfg)
    assert any(
        "segmented.viruses.fluA" in e and "unknown key 'segment_legnths'" in e
        for e in errs
    )


def test_underscore_prefixed_keys_skipped():
    """The sanctioned annotation / runtime-state escape hatch."""
    cfg = load_config(None)
    cfg["_my_note"] = "anything"
    cfg["qc"]["_reminder"] = "tighten later"
    assert validate_config(cfg) == []


def test_injected_verbose_key_not_flagged():
    """cli.py sets cfg['verbose'] before validation; it is not in DEFAULTS."""
    cfg = load_config(None)
    cfg["verbose"] = True
    assert validate_config(cfg) == []


def test_user_keyed_maps_not_audited():
    """Family labels / virus names / rank keys are DATA, not schema keys."""
    cfg = load_config(None)
    cfg["phylo"]["partition"]["models"] = {"L_RdRP": "LG+G4"}
    cfg["phylo"]["rooting"]["method"] = "outgroup"
    cfg["phylo"]["rooting"]["outgroup_rank"] = {"family": "Hantaviridae"}
    cfg["segmented"] = {
        "enabled": False,
        "viruses": {"whatever_name_i_like": {
            "expected_segments": 1, "segments": ["S"],
            "isolate_regex": r"(?P<isolate>X)",
        }},
    }
    assert validate_config(cfg) == []


def test_renamed_length_filter_single_message_not_doubled():
    """qc.length_filter gets its tailored migration message, NOT also a
    generic 'unknown config key' line."""
    cfg = load_config(None)
    cfg["qc"]["length_filter"] = {"mode": "median_percent"}
    errs = validate_config(cfg)
    assert any("renamed to qc.genome_length_filter" in e for e in errs)
    assert not any("unknown config key 'qc.length_filter'" in e for e in errs)


def test_shipped_example_configs_validate_clean():
    """Drift guard: the documented schema and the segmented HMM example must
    pass the unknown-key audit. Add a key to DEFAULTS or a spec without
    updating the matching _ALLOWED_*_KEYS set and this fails."""
    config_dir = Path(__file__).resolve().parent.parent / "config"
    for rel in ("default_config.yaml", "examples/alphainfluenzavirus_hmm.yaml"):
        cfg = load_config(config_dir / rel)
        assert validate_config(cfg) == [], f"{rel} did not validate clean"


def test_phylo_monophyly_min_support_default_is_70():
    cfg = load_config(None)
    assert cfg["phylo"]["monophyly"]["min_support"] == 70
    assert validate_config(cfg) == []


def test_phylo_monophyly_min_support_rejects_out_of_range():
    cfg = load_config(None)
    cfg["phylo"]["monophyly"]["min_support"] = 150
    assert any("phylo.monophyly.min_support" in e for e in validate_config(cfg))


def test_phylo_monophyly_min_support_zero_ok():
    cfg = load_config(None)
    cfg["phylo"]["monophyly"]["min_support"] = 0
    assert validate_config(cfg) == []
