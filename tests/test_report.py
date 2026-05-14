"""Report writers: run-log secret redaction, per-group counts TSV."""
from __future__ import annotations

from repseq.models import GroupStat, QCReport, RunResult
from repseq.output.report import write_group_counts_tsv, write_run_log


def test_write_run_log_redacts_ncbi_api_key(tmp_path):
    # The full config is dumped into the plaintext run log; a configured
    # NCBI API key must be redacted so it isn't leaked into the log file.
    cfg = {
        "output": {"dir": str(tmp_path), "prefix": "x"},
        "taxonomy": {"ncbi_email": "me@example.org", "ncbi_api_key": "SECRET123"},
    }
    log_path = tmp_path / "x_run.log"
    write_run_log(
        RunResult(mode="test"),
        QCReport(),
        cfg,
        input_paths=["in.fasta"],
        output_files=[],
        log_path=log_path,
    )
    text = log_path.read_text()
    assert "SECRET123" not in text
    assert "***redacted***" in text
    # Non-secret config is still recorded, and the caller's dict is untouched.
    assert "me@example.org" in text
    assert cfg["taxonomy"]["ncbi_api_key"] == "SECRET123"


# ---------------------------------------------------------------------------
# Per-group selection counts TSV
# ---------------------------------------------------------------------------

def test_write_group_counts_tsv(tmp_path):
    result = RunResult(
        mode="taxonomic1",
        group_stats=[
            GroupStat(grouping="genus", group="Alphainfluenzavirus",
                      n_before=1487, n_after=5, clustered=True, cutoff=0.8342),
            GroupStat(grouping="genus", group="Betainfluenzavirus",
                      n_before=3, n_after=3, clustered=False),
        ],
    )
    path = tmp_path / "x_group_counts.tsv"
    assert write_group_counts_tsv(result, path) is True
    lines = path.read_text().splitlines()
    assert lines[0] == "grouping\tgroup\tn_before\tn_after\tclustered\tcutoff"
    assert lines[1] == "genus\tAlphainfluenzavirus\t1487\t5\ttrue\t0.8342"
    # Group kept without clustering: cutoff column left blank.
    assert lines[2] == "genus\tBetainfluenzavirus\t3\t3\tfalse\t"


def test_write_group_counts_tsv_skips_when_no_stats(tmp_path):
    # A mode that recorded no group stats writes no file.
    path = tmp_path / "x_group_counts.tsv"
    assert write_group_counts_tsv(RunResult(mode="test"), path) is False
    assert not path.exists()
