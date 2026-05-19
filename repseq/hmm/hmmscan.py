"""hmmscan subprocess wrapper and ``--domtblout`` parser.

Batches arbitrarily many protein queries into a single hmmscan
invocation so the start-up cost is paid once per run, not per CDS.
"""
from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from .errors import HMMScanError


def is_available() -> bool:
    """True iff ``hmmscan`` is on PATH."""
    return shutil.which("hmmscan") is not None


def scan(
    db_path: Path,
    queries: dict[str, str],
    threads: int = 1,
    extra_args: Optional[list[str]] = None,
) -> dict[str, list[dict]]:
    """Run hmmscan on a batch of protein queries.

    Args:
        db_path: pressed HMM database (.h3* indexes must exist).
        queries: ``{query_id: protein_sequence}``. IDs must be unique and
            free of whitespace; the caller is responsible for that.
        threads: ``--cpu N``. 0 lets HMMER auto-pick.
        extra_args: appended verbatim to argv.

    Returns:
        ``{query_id: [hit_dict, ...]}`` where each hit_dict has keys
        ``target, target_acc, hmm_len, query_len, evalue, full_score,
        dom_evalue, dom_score, hmm_from, hmm_to, ali_from, ali_to,
        ali_span``. Queries with no hits are absent from the dict.
    """
    if not queries:
        return {}
    if not is_available():
        raise HMMScanError("hmmscan is not on PATH")
    extra = list(extra_args) if extra_args else []

    with tempfile.TemporaryDirectory(prefix="repseq_hmm_") as td:
        qfa = Path(td) / "queries.faa"
        out = Path(td) / "domtbl.tsv"
        with open(qfa, "w") as fh:
            for qid, seq in queries.items():
                fh.write(f">{qid}\n{seq}\n")
        argv = [
            "hmmscan",
            "--cpu", str(max(0, threads)),
            "--domtblout", str(out),
            "--noali",
            "--notextw",
            *extra,
            str(db_path),
            str(qfa),
        ]
        proc = subprocess.run(argv, capture_output=True, text=True)
        if proc.returncode != 0:
            raise HMMScanError(
                f"hmmscan exited {proc.returncode}\n"
                f"stderr: {proc.stderr.strip()}\n"
                f"stdout: {proc.stdout.strip()}"
            )
        return _parse_domtblout(out)


def _parse_domtblout(path: Path) -> dict[str, list[dict]]:
    """Parse hmmscan ``--domtblout`` (22 whitespace-separated columns).

    Columns: target, t-acc, tlen, query, q-acc, qlen, E-value, score,
    bias, #, of, c-Evalue, i-Evalue, dom-score, dom-bias, hmm-from,
    hmm-to, ali-from, ali-to, env-from, env-to, acc, description.
    """
    hits: dict[str, list[dict]] = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            # split(None, 22) keeps the trailing description (col 23) intact;
            # we ignore it but it keeps positional indexes stable.
            parts = line.split(None, 22)
            if len(parts) < 22:
                continue
            try:
                hit = {
                    "target": parts[0],
                    "target_acc": parts[1] if parts[1] != "-" else "",
                    "hmm_len": int(parts[2]),
                    "query_len": int(parts[5]),
                    "evalue": float(parts[6]),
                    "full_score": float(parts[7]),
                    "dom_evalue": float(parts[12]),
                    "dom_score": float(parts[13]),
                    "hmm_from": int(parts[15]),
                    "hmm_to": int(parts[16]),
                    "ali_from": int(parts[17]),
                    "ali_to": int(parts[18]),
                }
                hit["ali_span"] = hit["ali_to"] - hit["ali_from"] + 1
            except (ValueError, IndexError):
                continue
            hits.setdefault(parts[3], []).append(hit)
    return hits
