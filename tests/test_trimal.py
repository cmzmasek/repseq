"""trimAl wrapper: the maybe_trim soft-fail contract + provenance helpers.

subprocess / PATH are mocked so the disabled / missing / success / failure /
degenerate branches are locked without a real trimal binary.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from repseq.phylo.trimal import (
    _alignment_width,
    maybe_trim,
    tool_version,
    trim_note,
)


def _write_aln(path: Path, cols: int = 6, n: int = 3) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for i in range(n):
            fh.write(f">S{i + 1:04d}\n{'M' * cols}\n")


def _fake_trimal(out_cols: int = 3):
    """A subprocess.run stub that writes a trimmed alignment to -out."""
    captured: dict = {}

    def _run(cmd, **kwargs):
        captured["cmd"] = list(cmd)
        out = Path(cmd[cmd.index("-out") + 1])
        n_in = cmd[cmd.index("-in") + 1]
        n = sum(1 for line in Path(n_in).read_text().splitlines() if line.startswith(">"))
        with open(out, "w") as fh:
            for i in range(n):
                fh.write(f">S{i + 1:04d}\n{'M' * out_cols}\n")

        class _R:
            stdout = ""
            stderr = ""
        return _R()

    return _run, captured


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def test_alignment_width_reads_first_record(tmp_path):
    aln = tmp_path / "a.fasta"
    _write_aln(aln, cols=42)
    assert _alignment_width(aln) == 42


def test_trim_note_formats_version_and_mode():
    with patch("repseq.phylo.trimal.tool_version", return_value="trimAl v1.4.rev15"):
        assert trim_note({"mode": "automated1"}) == "trimAl v1.4.rev15 -automated1"
        assert trim_note({"mode": "gappyout", "extra_args": ["-keepheader"]}) == \
            "trimAl v1.4.rev15 -gappyout -keepheader"


def test_tool_version_unknown_when_missing(monkeypatch):
    monkeypatch.setattr("repseq.phylo.trimal.shutil.which", lambda _n: None)
    assert tool_version() == "unknown"


# ---------------------------------------------------------------------------
# maybe_trim
# ---------------------------------------------------------------------------

def test_maybe_trim_disabled_is_noop(tmp_path):
    src = tmp_path / "in.fasta"
    _write_aln(src)
    out = tmp_path / "out.fasta"
    assert maybe_trim(src, out, {}, {"enabled": False}) is False
    assert not out.exists()


def test_maybe_trim_missing_binary_warns_and_returns_false(tmp_path, capsys):
    src = tmp_path / "in.fasta"
    _write_aln(src)
    out = tmp_path / "out.fasta"
    with patch("repseq.phylo.trimal.shutil.which", return_value=None):
        ok = maybe_trim(src, out, {}, {"enabled": True, "mode": "automated1"},
                        label="genome")
    assert ok is False
    err = capsys.readouterr().err
    assert "trimal requested" in err and "UNTRIMMED" in err


def test_maybe_trim_success(tmp_path):
    src = tmp_path / "in.fasta"
    _write_aln(src, cols=10)
    out = tmp_path / "out.fasta"
    fake_run, captured = _fake_trimal(out_cols=4)
    with patch("repseq.phylo.trimal.shutil.which", return_value="/fake/trimal"), \
         patch("repseq.phylo.trimal.subprocess.run", side_effect=fake_run):
        ok = maybe_trim(src, out, {}, {"enabled": True, "mode": "gappyout"})
    assert ok is True
    assert _alignment_width(out) == 4
    # Mode → -gappyout flag; -in / -out present.
    assert "-gappyout" in captured["cmd"]
    assert "-in" in captured["cmd"] and "-out" in captured["cmd"]


def test_maybe_trim_extra_args_appended(tmp_path):
    src = tmp_path / "in.fasta"
    _write_aln(src)
    out = tmp_path / "out.fasta"
    fake_run, captured = _fake_trimal()
    with patch("repseq.phylo.trimal.shutil.which", return_value="/fake/trimal"), \
         patch("repseq.phylo.trimal.subprocess.run", side_effect=fake_run):
        maybe_trim(src, out, {}, {"enabled": True, "mode": "automated1",
                                  "extra_args": ["-gt", "0.8"]})
    cmd = captured["cmd"]
    assert "-gt" in cmd and cmd[cmd.index("-gt") + 1] == "0.8"


def test_maybe_trim_degenerate_output_falls_back(tmp_path, capsys):
    src = tmp_path / "in.fasta"
    _write_aln(src)
    out = tmp_path / "out.fasta"
    fake_run, _ = _fake_trimal(out_cols=0)  # trimal strips everything
    with patch("repseq.phylo.trimal.shutil.which", return_value="/fake/trimal"), \
         patch("repseq.phylo.trimal.subprocess.run", side_effect=fake_run):
        ok = maybe_trim(src, out, {}, {"enabled": True, "mode": "strict"})
    assert ok is False
    assert "degenerate" in capsys.readouterr().err


def test_maybe_trim_subprocess_failure_falls_back(tmp_path, capsys):
    import subprocess as _sp
    src = tmp_path / "in.fasta"
    _write_aln(src)
    out = tmp_path / "out.fasta"

    def _boom(cmd, **kwargs):
        raise _sp.CalledProcessError(1, cmd, output="", stderr="trimal: bad flag")

    with patch("repseq.phylo.trimal.shutil.which", return_value="/fake/trimal"), \
         patch("repseq.phylo.trimal.subprocess.run", side_effect=_boom):
        ok = maybe_trim(src, out, {}, {"enabled": True, "mode": "automated1"})
    assert ok is False
    err = capsys.readouterr().err
    assert "trimal failed" in err and "bad flag" in err
