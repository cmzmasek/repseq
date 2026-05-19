"""Tests for repseq.hmm.hmmscan (subprocess wrapper + domtblout parser)."""
from __future__ import annotations

import shutil
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from repseq.hmm.errors import HMMScanError
from repseq.hmm.hmmscan import _parse_domtblout, is_available, scan


# A realistic hmmscan --domtblout output. Columns (whitespace-separated):
#   target  t-acc  tlen  query  q-acc  qlen  E  score  bias  #  of  cEv  iEv
#   dom-score  dom-bias  hmm-from  hmm-to  ali-from  ali-to  env-from  env-to
#   acc  description
_SAMPLE_DOMTBL = """\
#                                                                            --- full sequence --- -------------- this domain -------------   hmm coord   ali coord   env coord
# target name        accession   tlen query name           accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ----------  ---- -------------------- ---------- -----  --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
Flu_NP               PF00506.24   498 flu_np               -            498   1.9e-79  257.6   0.1   1   1   2.1e-83   1.9e-79  257.6   0.1     1   498     1   498     1   498 1.00 Influenza virus nucleoprotein
RdRP_4               PF02123.22   465 sars_rdrp            -            900   3.5e-12   42.1   0.0   1   1   4.0e-16   3.5e-12   42.1   0.0    10   200   100   300   100   300 0.95 Viral RNA-directed RNA-polymerase
"""


def test_parse_domtblout_reads_two_hits(tmp_path):
    p = tmp_path / "domtbl.tsv"
    p.write_text(_SAMPLE_DOMTBL)
    hits = _parse_domtblout(p)
    assert set(hits.keys()) == {"flu_np", "sars_rdrp"}
    h = hits["flu_np"][0]
    assert h["target"] == "Flu_NP"
    assert h["target_acc"] == "PF00506.24"
    assert h["hmm_len"] == 498
    assert h["query_len"] == 498
    assert h["dom_evalue"] == pytest.approx(1.9e-79)
    assert h["dom_score"] == pytest.approx(257.6)
    assert h["hmm_from"] == 1
    assert h["hmm_to"] == 498
    assert h["ali_from"] == 1
    assert h["ali_to"] == 498
    assert h["ali_span"] == 498


def test_parse_domtblout_handles_partial_alignment(tmp_path):
    p = tmp_path / "domtbl.tsv"
    p.write_text(_SAMPLE_DOMTBL)
    hits = _parse_domtblout(p)
    rdrp = hits["sars_rdrp"][0]
    assert rdrp["ali_span"] == 300 - 100 + 1
    # Coverage check the marker selector cares about — short alignment
    # span over a long HMM model should give low coverage.
    assert rdrp["ali_span"] / rdrp["hmm_len"] < 0.5


def test_parse_domtblout_skips_comments_and_blank_lines(tmp_path):
    p = tmp_path / "domtbl.tsv"
    p.write_text("# this is a comment\n\n# another\n")
    assert _parse_domtblout(p) == {}


def test_parse_domtblout_skips_malformed_rows(tmp_path):
    p = tmp_path / "domtbl.tsv"
    # First row missing many columns — should be silently skipped, not crash.
    p.write_text(
        "broken row only two fields\n"
        + _SAMPLE_DOMTBL.splitlines()[3] + "\n"  # one good row
    )
    hits = _parse_domtblout(p)
    assert set(hits.keys()) == {"flu_np"}


def test_scan_empty_query_dict_returns_empty():
    assert scan(Path("/dev/null"), {}, threads=1) == {}


def test_scan_raises_when_hmmscan_missing(tmp_path):
    """When hmmscan is not on PATH the wrapper raises HMMScanError
    rather than running anything; the caller (cli._run_hmm_scan) is
    responsible for translating that into a soft-fail."""
    with patch("repseq.hmm.hmmscan.is_available", return_value=False):
        with pytest.raises(HMMScanError, match="not on PATH"):
            scan(Path("/dev/null"), {"q1": "MMMM"}, threads=1)


def test_scan_writes_query_fasta_and_parses_output(tmp_path):
    """Mock subprocess.run: confirm the wrapper builds the right argv
    and parses the resulting --domtblout. The Path(tmp_path) is a real
    temp dir so the tempfile creation inside scan() works."""
    sample = _SAMPLE_DOMTBL

    def fake_run(argv, capture_output=True, text=True):
        # Last two args are the database and the queries fasta.
        out_path = None
        for i, a in enumerate(argv):
            if a == "--domtblout":
                out_path = Path(argv[i + 1])
        assert out_path is not None
        out_path.write_text(sample)
        mock = MagicMock()
        mock.returncode = 0
        mock.stdout = ""
        mock.stderr = ""
        return mock

    db = tmp_path / "fake.hmm"
    db.write_text("HMMER3/f\nNAME Foo\n//\n")
    with patch("repseq.hmm.hmmscan.is_available", return_value=True), \
         patch("repseq.hmm.hmmscan.subprocess.run", side_effect=fake_run):
        hits = scan(db, {"flu_np": "MAS...", "sars_rdrp": "SAD..."}, threads=2)
    assert "flu_np" in hits
    assert "sars_rdrp" in hits


def test_scan_propagates_nonzero_exit_as_hmmscanerror(tmp_path):
    def fake_run(argv, capture_output=True, text=True):
        mock = MagicMock()
        mock.returncode = 1
        mock.stdout = ""
        mock.stderr = "fatal: bad database"
        return mock

    db = tmp_path / "fake.hmm"
    db.write_text("HMMER3/f\nNAME Foo\n//\n")
    with patch("repseq.hmm.hmmscan.is_available", return_value=True), \
         patch("repseq.hmm.hmmscan.subprocess.run", side_effect=fake_run):
        with pytest.raises(HMMScanError, match="bad database"):
            scan(db, {"q": "M"}, threads=1)
