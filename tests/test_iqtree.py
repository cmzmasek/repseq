"""IQ-TREE wrapper: argv construction, output capture, UFBoot auto-skip.

Subprocess is mocked: each test stubs ``subprocess.run`` to inspect the
argv IQ-TREE would have received and to drop the expected output files
(``.treefile``, ``.iqtree``) into the work directory so the wrapper's
file-copy logic runs end-to-end.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from repseq.phylo.iqtree import IQTreeError, _check_iqtree, _count_msa_records, run_iqtree


def _write_msa(path: Path, n: int = 3) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for i in range(n):
            fh.write(f">S{i + 1:04d}\nMKLPQE\n")


def _fake_run_writes_outputs(treefile_body: str = "(A:0.1,B:0.1,C:0.1);\n",
                             summary_body: str = "Best model: LG+G4\n"):
    """Build a fake run_streaming that drops the IQ-TREE output files.

    Matches the helper's signature: ``run_streaming(argv, *, stdout_file=None,
    cwd=None, check=True, stream_prefix="", stderr_dest=None) -> str``.
    """
    captured: dict = {}

    def _run(argv, **kwargs):
        captured["cmd"] = list(argv)
        # Find --prefix's value; the wrapper writes outputs to <prefix>.*.
        prefix = None
        for i, token in enumerate(argv):
            if token == "--prefix" and i + 1 < len(argv):
                prefix = Path(argv[i + 1])
                break
        assert prefix is not None, "wrapper must pass --prefix"
        Path(str(prefix) + ".treefile").write_text(treefile_body)
        Path(str(prefix) + ".iqtree").write_text(summary_body)
        return ""  # run_streaming returns the buffered stderr string

    return _run, captured


def test_check_iqtree_prefers_iqtree2(monkeypatch):
    """When both binaries are on PATH, iqtree2 wins."""
    def _which(name):
        return f"/fake/bin/{name}" if name in ("iqtree2", "iqtree") else None
    monkeypatch.setattr("repseq.phylo.iqtree.shutil.which", _which)
    assert _check_iqtree() == "/fake/bin/iqtree2"


def test_check_iqtree_falls_back_to_iqtree(monkeypatch):
    """When only iqtree exists, that's returned."""
    def _which(name):
        return "/fake/bin/iqtree" if name == "iqtree" else None
    monkeypatch.setattr("repseq.phylo.iqtree.shutil.which", _which)
    assert _check_iqtree() == "/fake/bin/iqtree"


def test_check_iqtree_honours_override(monkeypatch):
    """An explicit override is tried first."""
    def _which(name):
        return "/custom/path/myiqtree" if name == "myiqtree" else None
    monkeypatch.setattr("repseq.phylo.iqtree.shutil.which", _which)
    assert _check_iqtree("myiqtree") == "/custom/path/myiqtree"


def test_check_iqtree_raises_when_none_present(monkeypatch):
    monkeypatch.setattr("repseq.phylo.iqtree.shutil.which", lambda _name: None)
    with pytest.raises(IQTreeError, match="IQ-TREE not found"):
        _check_iqtree()


def test_count_msa_records(tmp_path):
    msa = tmp_path / "msa.fasta"
    _write_msa(msa, n=5)
    assert _count_msa_records(msa) == 5


def test_run_iqtree_prints_start_and_finish_to_stderr(tmp_path, capsys):
    """Bench-scientist progress: '[phylo] starting IQ-TREE (...)' before
    the subprocess call, '[phylo] IQ-TREE finished (Xs)' after success.
    Path-bearing args (-s INPUT, --prefix PFX) are stripped from the
    display string."""
    msa = tmp_path / "msa.fasta"
    _write_msa(msa, n=4)
    out = tmp_path / "out.nwk"
    cfg = {"temp_dir": str(tmp_path), "seed": 7, "threads": 2}

    fake_run, _ = _fake_run_writes_outputs()
    with patch("repseq.phylo.iqtree.run_streaming", side_effect=fake_run), \
         patch("repseq.phylo.iqtree._check_iqtree", return_value="/fake/iqtree2"):
        run_iqtree(msa, out, cfg, is_protein=True)

    err = capsys.readouterr().err
    assert "[phylo] starting IQ-TREE (" in err
    assert "[phylo] IQ-TREE finished (" in err
    assert "/fake/iqtree2" not in err
    assert "-s " not in err
    assert "--prefix" not in err
    assert "-m MFP" in err
    assert "-T 2" in err
    assert "-B 1000" in err


def test_run_iqtree_writes_treefile_and_summary(tmp_path):
    msa = tmp_path / "msa.fasta"
    _write_msa(msa, n=4)
    out = tmp_path / "out_tree.nwk"
    summary = tmp_path / "out_iqtree_summary.txt"
    cfg = {"temp_dir": str(tmp_path), "seed": 7, "threads": 2}

    fake_run, captured = _fake_run_writes_outputs(
        treefile_body="(A,B,C,D);\n", summary_body="Best model: WAG+G4\n",
    )
    with patch("repseq.phylo.iqtree.run_streaming", side_effect=fake_run), \
         patch("repseq.phylo.iqtree._check_iqtree", return_value="/fake/iqtree2"):
        run_iqtree(msa, out, cfg, is_protein=True, summary_path=summary)

    assert out.read_text() == "(A,B,C,D);\n"
    assert summary.read_text() == "Best model: WAG+G4\n"
    cmd = captured["cmd"]
    # Argv sanity: required flags present.
    assert "/fake/iqtree2" in cmd
    assert "-s" in cmd
    assert "-m" in cmd and "MFP" in cmd
    assert "-T" in cmd and "2" in cmd
    assert "-seed" in cmd and "7" in cmd
    assert "-B" in cmd and "1000" in cmd          # default UFBoot=1000
    assert "--quiet" in cmd
    assert "--redo" in cmd


def test_run_iqtree_disables_ufboot_when_fewer_than_four_seqs(tmp_path, capsys):
    """3 sequences → UFBoot skipped (IQ-TREE rejects bootstrap < 4)."""
    msa = tmp_path / "msa.fasta"
    _write_msa(msa, n=3)
    out = tmp_path / "out.nwk"
    cfg = {"temp_dir": str(tmp_path), "seed": 7, "threads": 1}

    fake_run, captured = _fake_run_writes_outputs()
    with patch("repseq.phylo.iqtree.run_streaming", side_effect=fake_run), \
         patch("repseq.phylo.iqtree._check_iqtree", return_value="/fake/iqtree2"):
        run_iqtree(msa, out, cfg, is_protein=True)

    cmd = captured["cmd"]
    assert "-B" not in cmd  # bootstrap dropped
    err = capsys.readouterr().err
    assert "skipping ultrafast bootstrap" in err


def test_run_iqtree_honours_explicit_model_and_extra_args(tmp_path):
    msa = tmp_path / "msa.fasta"
    _write_msa(msa, n=4)
    out = tmp_path / "out.nwk"
    cfg = {
        "temp_dir": str(tmp_path),
        "seed": 1, "threads": 1,
        "phylo": {"iqtree": {
            "model": "LG+G4",
            "ultrafast_bootstrap": 500,
            "extra_args": ["-alrt", "1000"],
        }},
    }

    fake_run, captured = _fake_run_writes_outputs()
    with patch("repseq.phylo.iqtree.run_streaming", side_effect=fake_run), \
         patch("repseq.phylo.iqtree._check_iqtree", return_value="/fake/iqtree2"):
        run_iqtree(msa, out, cfg, is_protein=True)

    cmd = captured["cmd"]
    assert "LG+G4" in cmd
    assert "500" in cmd      # UFBoot replicates
    # Extra args appended verbatim.
    i = cmd.index("-alrt")
    assert cmd[i + 1] == "1000"


def test_run_iqtree_can_disable_ufboot_via_config(tmp_path):
    msa = tmp_path / "msa.fasta"
    _write_msa(msa, n=10)  # enough seqs that the auto-skip would not fire
    out = tmp_path / "out.nwk"
    cfg = {
        "temp_dir": str(tmp_path), "seed": 1, "threads": 1,
        "phylo": {"iqtree": {"ultrafast_bootstrap": 0}},
    }

    fake_run, captured = _fake_run_writes_outputs()
    with patch("repseq.phylo.iqtree.run_streaming", side_effect=fake_run), \
         patch("repseq.phylo.iqtree._check_iqtree", return_value="/fake/iqtree2"):
        run_iqtree(msa, out, cfg, is_protein=True)

    cmd = captured["cmd"]
    assert "-B" not in cmd


def test_run_iqtree_summary_path_optional(tmp_path):
    """No summary_path → wrapper doesn't try to copy the .iqtree file."""
    msa = tmp_path / "msa.fasta"
    _write_msa(msa, n=4)
    out = tmp_path / "out.nwk"
    cfg = {"temp_dir": str(tmp_path), "seed": 1, "threads": 1}

    fake_run, _ = _fake_run_writes_outputs()
    with patch("repseq.phylo.iqtree.run_streaming", side_effect=fake_run), \
         patch("repseq.phylo.iqtree._check_iqtree", return_value="/fake/iqtree2"):
        run_iqtree(msa, out, cfg, is_protein=True)  # no summary_path

    assert out.exists()


def test_run_iqtree_raises_when_subprocess_fails(tmp_path):
    msa = tmp_path / "msa.fasta"
    _write_msa(msa, n=4)
    out = tmp_path / "out.nwk"
    cfg = {"temp_dir": str(tmp_path), "seed": 1, "threads": 1}

    from repseq.utils.subprocess_stream import StreamedProcessError

    def _fail(argv, **kwargs):
        raise StreamedProcessError(2, argv, output=None, stderr="boom")

    with patch("repseq.phylo.iqtree.run_streaming", side_effect=_fail), \
         patch("repseq.phylo.iqtree._check_iqtree", return_value="/fake/iqtree2"):
        with pytest.raises(IQTreeError, match="boom"):
            run_iqtree(msa, out, cfg, is_protein=True)


def test_run_iqtree_partition_uses_linkage_flag_and_drops_model(tmp_path):
    """A partition_file switches IQ-TREE into partitioned mode: the linkage
    flag (-p/-q/-Q) carries the NEXUS, -m is dropped (models live in the
    charpartition), and the partition-file path is copied into the workdir."""
    msa = tmp_path / "msa.fasta"
    _write_msa(msa, n=4)
    part = tmp_path / "partition.nex"
    part.write_text("#nexus\nbegin sets;\n  charset a = 1-3;\nend;\n")
    out = tmp_path / "out.nwk"
    cfg = {"temp_dir": str(tmp_path), "seed": 7, "threads": 2}

    fake_run, captured = _fake_run_writes_outputs()
    with patch("repseq.phylo.iqtree.run_streaming", side_effect=fake_run), \
         patch("repseq.phylo.iqtree._check_iqtree", return_value="/fake/iqtree2"):
        run_iqtree(
            msa, out, cfg, is_protein=True,
            partition_file=part, partition_linkage="unlinked",
        )

    cmd = captured["cmd"]
    assert "-Q" in cmd                       # unlinked → -Q
    assert "-m" not in cmd                   # model dropped under partitions
    # The flag is immediately followed by a NEXUS path inside the workdir.
    nexus_arg = cmd[cmd.index("-Q") + 1]
    assert nexus_arg.endswith("partition.nex")
    assert out.read_text()                   # tree copied out


def test_run_iqtree_partition_linkage_defaults_to_proportional(tmp_path):
    msa = tmp_path / "msa.fasta"
    _write_msa(msa, n=4)
    part = tmp_path / "partition.nex"
    part.write_text("#nexus\nbegin sets;\n  charset a = 1-3;\nend;\n")
    out = tmp_path / "out.nwk"
    cfg = {"temp_dir": str(tmp_path), "seed": 1, "threads": 1}

    fake_run, captured = _fake_run_writes_outputs()
    with patch("repseq.phylo.iqtree.run_streaming", side_effect=fake_run), \
         patch("repseq.phylo.iqtree._check_iqtree", return_value="/fake/iqtree2"):
        run_iqtree(msa, out, cfg, is_protein=True, partition_file=part)

    assert "-p" in captured["cmd"]           # proportional default


def test_run_iqtree_raises_when_treefile_missing(tmp_path):
    """IQ-TREE exits 0 but writes no .treefile → clear error, not silent
    success."""
    msa = tmp_path / "msa.fasta"
    _write_msa(msa, n=4)
    out = tmp_path / "out.nwk"
    cfg = {"temp_dir": str(tmp_path), "seed": 1, "threads": 1}

    def _silent_success(argv, **kwargs):
        # run_streaming returns the buffered stderr text; an empty string
        # means a quiet successful run.
        return ""

    with patch("repseq.phylo.iqtree.run_streaming", side_effect=_silent_success), \
         patch("repseq.phylo.iqtree._check_iqtree", return_value="/fake/iqtree2"):
        with pytest.raises(IQTreeError, match="did not produce a .treefile"):
            run_iqtree(msa, out, cfg, is_protein=True)
