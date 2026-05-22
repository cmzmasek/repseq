"""Friendly, traceback-free errors for common user mistakes.

Covers the three categories: config-file problems (load_config), input-FASTA
problems (_preflight_input), and the CLI boundary that renders them — plus
the external-tool classifier the boundary uses.
"""
from __future__ import annotations

from pathlib import Path

import pytest
from click.testing import CliRunner

from repseq.cli import _as_external_tool_error, _preflight_input, main
from repseq.config import load_config
from repseq.errors import ConfigError, InputError, RepseqError


# ---------------------------------------------------------------------------
# Config file
# ---------------------------------------------------------------------------

def test_load_config_missing_file_raises_config_error():
    with pytest.raises(ConfigError) as exc:
        load_config("/no/such/repseq_config.yaml")
    assert "Config file not found" in str(exc.value)
    assert "-c/--config" in str(exc.value)


def test_load_config_malformed_yaml_raises_with_location(tmp_path: Path):
    p = tmp_path / "bad.yaml"
    # Unbalanced bracket → PyYAML parse error with a problem_mark.
    p.write_text("threads: 4\nclustering: [unterminated\n")
    with pytest.raises(ConfigError) as exc:
        load_config(p)
    msg = str(exc.value)
    assert "is not valid YAML" in msg
    assert "line" in msg  # the problem_mark location made it in


def test_load_config_non_mapping_root_raises(tmp_path: Path):
    p = tmp_path / "list.yaml"
    p.write_text("- just\n- a\n- list\n")
    with pytest.raises(ConfigError) as exc:
        load_config(p)
    assert "must contain a mapping" in str(exc.value)
    assert "list" in str(exc.value)


def test_load_config_empty_file_is_ok(tmp_path: Path):
    """An empty config means 'all defaults', not an error."""
    p = tmp_path / "empty.yaml"
    p.write_text("")
    cfg = load_config(p)
    assert cfg["threads"] == 4  # defaults applied


def test_config_error_is_repseq_error():
    assert issubclass(ConfigError, RepseqError)
    assert issubclass(InputError, RepseqError)


# ---------------------------------------------------------------------------
# Input FASTA preflight
# ---------------------------------------------------------------------------

def test_preflight_missing_file():
    with pytest.raises(InputError) as exc:
        _preflight_input("/no/such/seqs.fasta")
    assert "Input file not found" in str(exc.value)


def test_preflight_directory(tmp_path: Path):
    with pytest.raises(InputError) as exc:
        _preflight_input(str(tmp_path))
    assert "is a directory" in str(exc.value)


def test_preflight_empty_file(tmp_path: Path):
    p = tmp_path / "empty.fasta"
    p.write_text("")
    with pytest.raises(InputError) as exc:
        _preflight_input(str(p))
    assert "empty" in str(exc.value)


def test_preflight_not_fasta(tmp_path: Path):
    p = tmp_path / "notes.txt"
    p.write_text("this is not fasta\njust some text\n")
    with pytest.raises(InputError) as exc:
        _preflight_input(str(p))
    assert "does not look like FASTA" in str(exc.value)


def test_preflight_valid_fasta_passes(tmp_path: Path):
    p = tmp_path / "ok.fasta"
    p.write_text(">seq1\nACGTACGT\n")
    _preflight_input(str(p))  # no exception


def test_preflight_tolerates_leading_blank_lines(tmp_path: Path):
    p = tmp_path / "lead.fasta"
    p.write_text("\n\n>seq1\nACGT\n")
    _preflight_input(str(p))  # whitespace before the first '>' is fine


# ---------------------------------------------------------------------------
# External-tool classifier
# ---------------------------------------------------------------------------

def test_as_external_tool_error_recognises_mmseqs():
    from repseq.clustering.mmseqs2 import MMseqs2Error

    e = MMseqs2Error("mmseqs2 not found in PATH")
    assert _as_external_tool_error(e) is e


def test_as_external_tool_error_recognises_cdhit():
    from repseq.clustering.cdhit import CDHitError

    e = CDHitError("cd-hit not found in PATH")
    assert _as_external_tool_error(e) is e


def test_as_external_tool_error_ignores_other():
    assert _as_external_tool_error(ValueError("boom")) is None
    assert _as_external_tool_error(RuntimeError("generic")) is None


# ---------------------------------------------------------------------------
# CLI boundary: friendly rendering, no traceback
# ---------------------------------------------------------------------------

def test_cli_missing_config_is_friendly():
    runner = CliRunner()
    result = runner.invoke(
        main,
        ["global", "-c", "/no/such/cfg.yaml", "-i", "x.fasta", "-T", "0.9"],
    )
    assert result.exit_code == 1
    assert "Error: Config file not found" in result.output
    assert "Traceback" not in result.output


def test_cli_missing_input_is_friendly(tmp_path: Path):
    runner = CliRunner()
    result = runner.invoke(
        main,
        ["global", "-i", str(tmp_path / "nope.fasta"), "-T", "0.9", "--no-resolve"],
    )
    assert result.exit_code == 1
    assert "Error: Input file not found" in result.output
    assert "Traceback" not in result.output


def test_cli_unexpected_error_keeps_traceback(monkeypatch, tmp_path: Path):
    """A genuine bug must still surface as a traceback (by design)."""
    p = tmp_path / "ok.fasta"
    p.write_text(">s1\nACGT\n")

    def _boom(*a, **k):
        raise ValueError("synthetic internal bug")

    monkeypatch.setattr("repseq.cli._load_sequences", _boom)
    runner = CliRunner()
    result = runner.invoke(
        main,
        ["global", "-i", str(p), "-T", "0.9", "--no-resolve"],
    )
    assert result.exit_code != 0
    # CliRunner stores the propagated exception rather than swallowing it.
    assert isinstance(result.exception, ValueError)
