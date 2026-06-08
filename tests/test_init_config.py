"""`repseq init-config` — the config wizard now emits the COMPLETE,
fully-commented reference config (every current setting + inline docs) with
the user's short-Q&A answers overlaid. These tests pin that contract: the
output validates clean, carries every advanced section, preserves the
comments, and round-trips the wizard's answers (incl. a segmented virus)."""
from __future__ import annotations

from pathlib import Path

from click.testing import CliRunner

from repseq.cli import (
    _overlay_value,
    _reference_config_text,
    _render_segmented_virus,
    main,
)
from repseq.config import load_config, validate_config


# --- the bundled reference -------------------------------------------------

def test_reference_config_is_bundled_and_matches_repo():
    """The wizard reads the packaged copy; it must be the single source of
    truth (byte-identical to the repo file it was moved to)."""
    text = _reference_config_text()
    repo = (
        Path(__file__).resolve().parent.parent
        / "repseq" / "data" / "default_config.yaml"
    )
    assert text == repo.read_text()
    assert len(text.splitlines()) > 500  # the full commented schema, not a stub


# --- overlay helpers -------------------------------------------------------

def test_overlay_value_replaces_unique_anchor():
    assert _overlay_value("threads: 4\n", "threads: 4", "threads: 8") == "threads: 8\n"


def test_overlay_value_leaves_text_when_anchor_absent():
    """Template drift must never corrupt the file — a missing anchor is a
    no-op, so the user still gets a complete file with the default in place."""
    txt = "threads: 4\n"
    assert _overlay_value(txt, "nonexistent: x", "nonexistent: y") == txt


def test_overlay_value_leaves_text_when_anchor_not_unique():
    txt = "enabled: false\nenabled: false\n"
    assert _overlay_value(txt, "enabled: false", "enabled: true") == txt


def test_render_segmented_virus_indents_under_viruses():
    block = _render_segmented_virus(
        "flu", {"expected_segments": 2, "segments": ["L", "S"], "isolate_regex": "x"}
    )
    assert block.startswith("  viruses:\n")
    assert "    flu:\n" in block
    assert "      expected_segments: 2" in block


# --- end-to-end wizard runs ------------------------------------------------

def _run(tmp_path: Path, answers: str) -> Path:
    out = tmp_path / "cfg.yaml"
    runner = CliRunner()
    result = runner.invoke(main, ["init-config", "-o", str(out)], input=answers)
    assert result.exit_code == 0, result.output
    return out


def test_init_config_nonsegmented_validates_and_has_every_section(tmp_path):
    out = _run(tmp_path, "./myout\nmyrep\n16\nme@x.org\n\n\nn\n")
    cfg = load_config(out)
    assert validate_config(cfg) == []

    # overlaid answers round-tripped
    assert cfg["threads"] == 16
    assert cfg["output"]["prefix"] == "myrep"
    assert cfg["taxonomy"]["ncbi_email"] == "me@x.org"

    # every advanced section the old sparse wizard dropped is present
    assert "hmm" in cfg and cfg["hmm"]["enabled"] is True
    assert "phylo" in cfg and "partition" in cfg["phylo"]
    assert cfg["clustering"]["alphabet_for_clustering"] == "protein"
    assert "cluster_protein" in cfg["clustering"]
    assert "priority" in cfg["representative"]
    assert cfg["overrides"]["exclude"] is not None
    assert cfg["qc"]["protein_quality"]["enabled"] is False
    assert cfg["segmented"]["taxonomy_consistency"]["rank"] == "species"


def test_init_config_output_is_richly_commented(tmp_path):
    """'explains everything well' — the generated file keeps the inline docs
    (not a stripped key dump). Spot-check a comment from a few sections."""
    out = _run(tmp_path, "\n\n\n\n\n\nn\n")
    text = out.read_text()
    assert text.count("#") > 300
    assert "# repseq default configuration" in text
    assert "Unknown keys are a HARD ERROR" in text


def test_init_config_blank_ncbi_email_stays_null(tmp_path):
    out = _run(tmp_path, "\n\n\n\n\n\nn\n")
    cfg = load_config(out)
    assert cfg["taxonomy"]["ncbi_email"] is None


def test_init_config_cdhit_backend_overlaid(tmp_path):
    out = _run(tmp_path, "\n\n\n\n\ncdhit\nn\n")
    cfg = load_config(out)
    assert cfg["clustering"]["backend"] == "cdhit"
    assert validate_config(cfg) == []


def test_init_config_segmented_virus_validates(tmp_path):
    answers = "\n\n\n\n\n\ny\ntestvirus\n3\nL,M,S\n\nn\nn\ny\n"
    out = _run(tmp_path, answers)
    cfg = load_config(out)
    assert validate_config(cfg) == []
    seg = cfg["segmented"]
    assert seg["enabled"] is True
    assert seg["virus"] == "testvirus"
    v = seg["viruses"]["testvirus"]
    assert v["expected_segments"] == 3
    assert v["segments"] == ["L", "M", "S"]


def test_init_config_segmented_defined_but_not_enabled(tmp_path):
    """User can store a virus definition without flipping enabled on (activate
    per-run with --segmented)."""
    answers = "\n\n\n\n\n\ny\ntestvirus\n3\nL,M,S\n\nn\nn\nn\n"
    out = _run(tmp_path, answers)
    cfg = load_config(out)
    assert validate_config(cfg) == []
    assert cfg["segmented"]["enabled"] is False
    assert "testvirus" in cfg["segmented"]["viruses"]


def test_init_config_segmented_with_protein_counts(tmp_path):
    answers = (
        "\n\n\n\n\n\ny\nflu\n2\nPB1,PB2\n\n"  # dir..backend, segmented y, name, n, segs, regex
        "y\n1,2\n1\n"                            # protein counts: yes, PB1=1,2  PB2=1
        "n\ny\n"                                 # length bounds no, enable yes
    )
    out = _run(tmp_path, answers)
    cfg = load_config(out)
    assert validate_config(cfg) == []
    epps = cfg["segmented"]["viruses"]["flu"]["expected_proteins_per_segment"]
    assert epps == {"PB1": [1, 2], "PB2": 1}
