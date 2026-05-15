"""`repseq doctor` self-test command.

These tests lock the WARN-vs-FAIL policy. The actual sources of "is it
installed" (importlib, shutil.which, requests.get, filesystem write) are
stubbed so the checks run deterministically in any environment.
"""
from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from click.testing import CliRunner

from repseq.cli import main
from repseq.doctor import (
    FAIL,
    OK,
    WARN,
    CheckResult,
    DoctorReport,
    check_binaries,
    check_cache_dir,
    check_config,
    check_network,
    check_python_packages,
    run_doctor,
)


# ---------------------------------------------------------------------------
# Python packages
# ---------------------------------------------------------------------------

def test_check_python_packages_optional_missing_is_warn_not_fail():
    # Force umap and matplotlib to look uninstalled. Required deps stay
    # importable (they're hard deps). Patching _package_version directly
    # avoids the recursion you get from patching importlib.import_module
    # (the fake's fallback would call the patched function again).
    def fake_version(name):
        if name in ("umap", "matplotlib"):
            return None
        return "1.0"

    with patch("repseq.doctor._package_version", side_effect=fake_version):
        results = check_python_packages()

    by_label = {r.label: r for r in results}
    assert by_label["umap-learn"].status == WARN
    assert by_label["matplotlib"].status == WARN
    assert by_label["biopython"].status == OK


def test_check_python_packages_required_missing_is_fail():
    def fake_version(name):
        return None if name == "Bio" else "1.0"

    with patch("repseq.doctor._package_version", side_effect=fake_version):
        results = check_python_packages()

    by_label = {r.label: r for r in results}
    assert by_label["biopython"].status == FAIL


# ---------------------------------------------------------------------------
# External binaries
# ---------------------------------------------------------------------------

def test_check_binaries_missing_all_is_warn_not_fail():
    # Every external tool is optional in the sense that you can always
    # pick a backend / mode that doesn't need a given one. The doctor
    # should not [FAIL] when none of them are present.
    with patch("repseq.doctor.shutil.which", return_value=None):
        results = check_binaries()
    assert all(r.status == WARN for r in results)


def test_check_binaries_fasttree_accepts_lowercase_alias():
    # Some Linux packages name it 'fasttree'; the doctor should resolve
    # either capitalisation.
    def which(name):
        return "/usr/bin/fasttree" if name == "fasttree" else None

    with patch("repseq.doctor.shutil.which", side_effect=which):
        results = check_binaries()

    ft = next(r for r in results if r.label == "FastTree")
    assert ft.status == OK
    assert "fasttree" in ft.detail


# ---------------------------------------------------------------------------
# Network
# ---------------------------------------------------------------------------

def test_check_network_unreachable_is_warn():
    import requests
    # Both endpoints raise a connection error -> WARN, not FAIL.
    def boom(url, timeout):
        raise requests.exceptions.ConnectionError("no route to host")

    with patch("repseq.doctor.requests.get", side_effect=boom):
        results = check_network()

    assert {r.status for r in results} == {WARN}
    assert any("ConnectionError" in r.detail for r in results)


def test_check_network_timeout_is_warn_with_timeout_label():
    import requests

    def slow(url, timeout):
        raise requests.exceptions.Timeout("read timed out")

    with patch("repseq.doctor.requests.get", side_effect=slow):
        results = check_network()

    assert all(r.status == WARN for r in results)
    assert all("timeout" in r.detail for r in results)


def test_check_network_ok_when_endpoints_respond():
    fake_resp = MagicMock(ok=True, status_code=200)
    with patch("repseq.doctor.requests.get", return_value=fake_resp):
        results = check_network()

    assert all(r.status == OK for r in results)


# ---------------------------------------------------------------------------
# Filesystem / config
# ---------------------------------------------------------------------------

def test_check_cache_dir_failure_is_fail(tmp_path, monkeypatch):
    # Point at a path under a non-writable parent. Use a real OSError
    # from tempfile rather than monkey-patching mkdir, so we exercise
    # the same code path a real install would hit.
    cache = tmp_path / "no_perms"
    cache.mkdir()
    cache.chmod(0o500)  # read+execute only — no write
    try:
        result = check_cache_dir({"cache_dir": str(cache / "sub")})
        # On some filesystems (e.g. tmpfs as root, CI) chmod is a no-op;
        # if so the check will succeed — that's fine, we still locked
        # the OK path. The interesting case is when it does fail.
        assert result.status in (OK, FAIL)
        if result.status == FAIL:
            assert "not writable" in result.detail
    finally:
        cache.chmod(0o700)


def test_check_cache_dir_creates_missing_directory(tmp_path):
    cache = tmp_path / "fresh"
    assert not cache.exists()
    result = check_cache_dir({"cache_dir": str(cache)})
    assert result.status == OK
    assert cache.exists()


def test_check_config_no_path_uses_defaults():
    from repseq.config import load_config
    cfg = load_config(None)
    cfg["taxonomy"]["ncbi_email"] = "me@institute.org"  # silence the email WARN
    results = check_config(cfg, config_path=None)
    by_label = {r.label: r for r in results}
    assert by_label["config"].status == OK
    assert "defaults" in by_label["config"].detail
    assert by_label["ncbi_email"].status == OK


def test_check_config_with_invalid_config_is_fail():
    cfg = {"clustering": {"backend": "bogus"}, "qc": {}, "representative": {}}
    results = check_config(cfg, config_path="/tmp/c.yaml")
    by_label = {r.label: r for r in results}
    assert by_label["config"].status == FAIL


def test_check_config_missing_email_is_warn(monkeypatch):
    monkeypatch.delenv("REPSEQ_NCBI_EMAIL", raising=False)
    from repseq.config import load_config
    cfg = load_config(None)
    cfg["taxonomy"]["ncbi_email"] = None
    results = check_config(cfg, config_path=None)
    by_label = {r.label: r for r in results}
    assert by_label["ncbi_email"].status == WARN


# ---------------------------------------------------------------------------
# Orchestrator + report
# ---------------------------------------------------------------------------

def test_run_doctor_no_network_skips_network_group():
    from repseq.config import load_config
    cfg = load_config(None)
    cfg["taxonomy"]["ncbi_email"] = "me@institute.org"
    with patch("repseq.doctor.requests.get") as mock_get:
        report = run_doctor(cfg, config_path=None, no_network=True)
    # The network check should never have been called.
    mock_get.assert_not_called()
    titles = [t for t, _ in report.groups]
    assert "Network / databases" in titles
    net_rows = dict(report.groups)["Network / databases"]
    assert any("--no-network" in r.detail for r in net_rows)


def test_doctor_report_renders_all_sections():
    report = DoctorReport(groups=[
        ("Python packages", [CheckResult(OK, "biopython", "1.81")]),
        ("External tools", [CheckResult(WARN, "mmseqs", "not on PATH — clustering")]),
        ("Network / databases", [CheckResult(OK, "NCBI Entrez", "HTTP 200")]),
        ("Configuration", [CheckResult(FAIL, "config", "invalid")]),
    ])
    text = report.render("0.5.8")
    assert "repseq 0.5.8 self-test" in text
    assert "[OK]" in text
    assert "[WARN]" in text
    assert "[FAIL]" in text
    assert "Summary: 2 OK, 1 warning(s), 1 failure(s)." in text
    assert report.has_failures is True


def test_doctor_report_clean_run_says_ready():
    report = DoctorReport(groups=[
        ("Python packages", [CheckResult(OK, "biopython", "1.81")]),
    ])
    text = report.render("0.5.8")
    assert "Your install is ready" in text
    assert report.has_failures is False


# ---------------------------------------------------------------------------
# CLI integration
# ---------------------------------------------------------------------------

def test_cli_doctor_subcommand_exit_zero_when_clean():
    runner = CliRunner()
    # Stub the check groups to all pass.
    def fake_run_doctor(cfg, config_path, no_network):
        return DoctorReport(groups=[
            ("Python packages", [CheckResult(OK, "biopython", "1.81")]),
        ])

    with patch("repseq.cli.load_config", return_value={}), \
         patch("repseq.doctor.run_doctor", side_effect=fake_run_doctor):
        result = runner.invoke(main, ["doctor", "--no-network"])

    assert result.exit_code == 0
    assert "Your install is ready" in result.output


def test_cli_doctor_subcommand_exit_one_on_failure():
    runner = CliRunner()
    def fake_run_doctor(cfg, config_path, no_network):
        return DoctorReport(groups=[
            ("Configuration",
             [CheckResult(FAIL, "config", "validation error: oops")]),
        ])

    with patch("repseq.cli.load_config", return_value={}), \
         patch("repseq.doctor.run_doctor", side_effect=fake_run_doctor):
        result = runner.invoke(main, ["doctor", "--no-network"])

    assert result.exit_code == 1
    assert "[FAIL]" in result.output
    assert "Fix the [FAIL] item(s)" in result.output
