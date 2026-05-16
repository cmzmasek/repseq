"""Self-test for `repseq doctor`.

The bench scientists who run this aren't going to chase a stack trace
when something's missing — they want a plain English answer to "is this
installation ready to use, and if not, what should I install?". The
``doctor`` command emits a categorised pass/warn/fail report and exits
non-zero only when something *required* is missing or broken.

WARN vs FAIL policy:
  * required Python packages (biopython, click, pyyaml, requests) → FAIL
    if missing — but they're hard dependencies so this almost never
    fires; included anyway because the failure mode would otherwise be
    a bare ImportError at start-up.
  * optional Python packages (umap-learn, matplotlib) → WARN — only used
    by ``--plot``.
  * every external clustering / phylogeny binary (mmseqs, cd-hit,
    cd-hit-est, mafft, FastTree, iqtree2) → WARN. None of them is
    strictly required: you can pick the backend you have installed,
    or use diversity-only modes that don't shell out.
  * network checks (NCBI Entrez + UniProt REST) → WARN when unreachable;
    you can still run with ``--no-resolve``.
  * cache directory writable → FAIL if not, since every resolved run
    needs to write there.
  * config (if a config file was passed) → FAIL on validation errors;
    WARN if no NCBI email is configured (works but you'll be rate-limited).
"""

from __future__ import annotations

import importlib
import importlib.metadata
import importlib.util
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Optional

import requests


OK, WARN, FAIL = "OK", "WARN", "FAIL"


@dataclass
class CheckResult:
    """One row of the report.

    ``label`` is the short name shown in the left column; ``detail`` is
    the longer human-readable note shown after the tag. ``status`` is
    one of ``OK``/``WARN``/``FAIL``.
    """
    status: str
    label: str
    detail: str = ""


# ---------------------------------------------------------------------------
# Python packages
# ---------------------------------------------------------------------------

# Each entry: (import name, display name shown in the report — typically the
# pip-install name so bench scientists can grep for what to install).
_REQUIRED_PACKAGES: tuple[tuple[str, str], ...] = (
    ("Bio", "biopython"),
    ("click", "click"),
    ("yaml", "PyYAML"),
    ("requests", "requests"),
)
_OPTIONAL_PACKAGES: tuple[tuple[str, str, str], ...] = (
    ("umap", "umap-learn", "only needed for --plot (pip install 'repseq[viz]')"),
    ("matplotlib", "matplotlib", "only needed for --plot (pip install 'repseq[viz]')"),
)


# Import-name → pip-install-name mapping for the version lookup.
_IMPORT_TO_DIST: dict[str, str] = {
    "Bio": "biopython",
    "yaml": "PyYAML",
    "umap": "umap-learn",
}


def _package_version(name: str) -> Optional[str]:
    """Return the installed version of an actually-importable package, or None.

    We attempt the real import (rather than just ``find_spec``) so that a
    broken install — e.g. a numpy/scipy ABI mismatch that lets the module
    be found but raises at load time — is reported as "not installed"
    instead of silently passing. ``importlib.metadata`` then supplies the
    version without going through ``module.__version__`` (which click 9.x
    and others have deprecated).
    """
    import warnings as _warnings
    try:
        with _warnings.catch_warnings():
            _warnings.simplefilter("ignore", DeprecationWarning)
            importlib.import_module(name)
    except ImportError:
        return None
    dist_name = _IMPORT_TO_DIST.get(name, name)
    try:
        return importlib.metadata.version(dist_name)
    except importlib.metadata.PackageNotFoundError:
        return "?"


def check_python_packages() -> list[CheckResult]:
    results: list[CheckResult] = []
    for import_name, display_name in _REQUIRED_PACKAGES:
        ver = _package_version(import_name)
        if ver is None:
            results.append(CheckResult(FAIL, display_name, "not installed (required)"))
        else:
            results.append(CheckResult(OK, display_name, ver))
    for import_name, display_name, note in _OPTIONAL_PACKAGES:
        ver = _package_version(import_name)
        if ver is None:
            results.append(CheckResult(WARN, display_name, f"not installed — {note}"))
        else:
            results.append(CheckResult(OK, display_name, ver))
    return results


# ---------------------------------------------------------------------------
# External binaries
# ---------------------------------------------------------------------------

# Each entry: (binary name, reason it's needed)
_EXTERNAL_BINARIES: tuple[tuple[str, str], ...] = (
    ("mmseqs",     "clustering (default backend)"),
    ("cd-hit",     "clustering (alternative backend, protein)"),
    ("cd-hit-est", "clustering (alternative backend, nucleotide)"),
    ("mafft",      "MSA, --phylo only"),
    ("FastTree",   "phylogeny, --phylo only on nucleotide (also accepts 'fasttree')"),
    ("iqtree2",    "phylogeny, --phylo only on protein (also accepts 'iqtree')"),
)


def _which_either(*names: str) -> Optional[str]:
    """Return the first name on PATH (handles FastTree / fasttree case-flip)."""
    for n in names:
        path = shutil.which(n)
        if path:
            return path
    return None


def check_binaries() -> list[CheckResult]:
    results: list[CheckResult] = []
    for name, note in _EXTERNAL_BINARIES:
        if name == "FastTree":
            path = _which_either("FastTree", "fasttree")
        elif name == "iqtree2":
            path = _which_either("iqtree2", "iqtree")
        else:
            path = shutil.which(name)
        if path:
            results.append(CheckResult(OK, name, path))
        else:
            results.append(CheckResult(WARN, name, f"not on PATH — {note}"))
    return results


# ---------------------------------------------------------------------------
# Network / external databases
# ---------------------------------------------------------------------------

_NCBI_PING_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/einfo.fcgi"
_UNIPROT_PING_URL = "https://rest.uniprot.org/uniprotkb/P12345.json"
_NETWORK_TIMEOUT = 5.0  # seconds — long enough for a slow DNS, short enough not to hang


def _ping(url: str, timeout: float = _NETWORK_TIMEOUT) -> tuple[bool, str]:
    """GET ``url`` with a short timeout; return (ok, status-or-error-string)."""
    try:
        resp = requests.get(url, timeout=timeout)
        if resp.ok:
            return True, f"HTTP {resp.status_code}"
        return False, f"HTTP {resp.status_code}"
    except requests.exceptions.Timeout:
        return False, f"timeout after {timeout:.0f}s"
    except requests.exceptions.RequestException as exc:
        return False, type(exc).__name__


def check_network() -> list[CheckResult]:
    results: list[CheckResult] = []
    for label, url in (("NCBI Entrez", _NCBI_PING_URL),
                       ("UniProt REST", _UNIPROT_PING_URL)):
        ok, msg = _ping(url)
        if ok:
            results.append(CheckResult(OK, label, f"{url} ({msg})"))
        else:
            results.append(CheckResult(
                WARN, label,
                f"{url} unreachable ({msg}); you can still run with --no-resolve",
            ))
    return results


# ---------------------------------------------------------------------------
# Filesystem / configuration
# ---------------------------------------------------------------------------

def check_cache_dir(cfg: dict[str, Any]) -> CheckResult:
    """Verify the cache directory exists (or can be created) and is writable."""
    raw = cfg.get("cache_dir") or "~/.repseq/cache"
    cache = Path(raw).expanduser()
    try:
        cache.mkdir(parents=True, exist_ok=True)
        # touch a probe file to confirm write access; some shared filesystems
        # let you mkdir but not write inside the directory.
        with tempfile.NamedTemporaryFile(dir=cache, delete=True):
            pass
    except OSError as exc:
        return CheckResult(
            FAIL, "cache directory",
            f"{cache} not writable ({exc.strerror or exc}); "
            "set cache_dir: in your config to a writable path",
        )
    return CheckResult(OK, "cache directory", str(cache))


def check_config(cfg: dict[str, Any], config_path: Optional[str]) -> list[CheckResult]:
    """Validate the config and surface whether an NCBI email is set."""
    from .config import validate_config  # local to avoid circular import

    results: list[CheckResult] = []
    errors = validate_config(cfg)
    if config_path is None:
        results.append(CheckResult(
            OK, "config", "using built-in defaults (no -c passed)",
        ))
    elif errors:
        results.append(CheckResult(
            FAIL, "config",
            f"{config_path} has {len(errors)} error(s): " + "; ".join(errors),
        ))
    else:
        results.append(CheckResult(OK, "config", f"{config_path} validates"))

    email = (cfg.get("taxonomy", {}) or {}).get("ncbi_email")
    if not email and not os.environ.get("REPSEQ_NCBI_EMAIL"):
        results.append(CheckResult(
            WARN, "ncbi_email",
            "not set — NCBI rate-limits anonymous traffic; "
            "set taxonomy.ncbi_email in your config or export REPSEQ_NCBI_EMAIL",
        ))
    else:
        results.append(CheckResult(OK, "ncbi_email", email or "REPSEQ_NCBI_EMAIL env"))
    return results


# ---------------------------------------------------------------------------
# Orchestrator + report formatting
# ---------------------------------------------------------------------------

@dataclass
class DoctorReport:
    groups: list[tuple[str, list[CheckResult]]]

    def render(self, version: str) -> str:
        lines: list[str] = []
        header = f"repseq {version} self-test"
        lines.append(header)
        lines.append("=" * len(header))
        lines.append("")
        label_width = max(
            (len(r.label) for _, rows in self.groups for r in rows),
            default=12,
        )
        for title, rows in self.groups:
            lines.append(title)
            for r in rows:
                tag = f"[{r.status}]".ljust(7)
                lines.append(f"  {tag} {r.label.ljust(label_width)}  {r.detail}")
            lines.append("")
        # Tally + closing line.
        flat = [r for _, rows in self.groups for r in rows]
        n_ok = sum(1 for r in flat if r.status == OK)
        n_warn = sum(1 for r in flat if r.status == WARN)
        n_fail = sum(1 for r in flat if r.status == FAIL)
        lines.append(
            f"Summary: {n_ok} OK, {n_warn} warning(s), {n_fail} failure(s)."
        )
        if n_fail == 0:
            lines.append(
                "Your install is ready. Warnings are optional pieces — "
                "install them only if you need the corresponding feature."
            )
        else:
            lines.append(
                "Fix the [FAIL] item(s) above before running a real job."
            )
        return "\n".join(lines)

    @property
    def has_failures(self) -> bool:
        return any(r.status == FAIL for _, rows in self.groups for r in rows)


def run_doctor(
    cfg: dict[str, Any],
    config_path: Optional[str],
    no_network: bool,
) -> DoctorReport:
    """Run every check group and return the assembled report."""
    groups: list[tuple[str, list[CheckResult]]] = []
    groups.append(("Python packages", check_python_packages()))
    groups.append(("External tools", check_binaries()))
    if no_network:
        groups.append((
            "Network / databases",
            [CheckResult(WARN, "skipped", "--no-network passed")],
        ))
    else:
        groups.append(("Network / databases", check_network()))
    groups.append(("Configuration", [check_cache_dir(cfg), *check_config(cfg, config_path)]))
    return DoctorReport(groups=groups)
