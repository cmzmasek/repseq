"""cd-hit clustering wrapper.

Alternative to MMseqs2 with a near-identical interface (``run_clustering``)
so it can be dispatched behind ``clustering.backend``.

cd-hit and cd-hit-est diverge from MMseqs2 in three places worth knowing:

* Two binaries — ``cd-hit`` for protein, ``cd-hit-est`` for nucleotide.
  The wrapper picks based on the input alphabet.
* Word size (``-n``) is not a free parameter: cd-hit refuses thresholds
  that don't match the table in its user guide. We auto-pick a valid ``-n``
  from the requested identity unless one is set in config.
* Identity floors: cd-hit refuses ``-c`` below 0.40 (protein) or 0.80
  (nucleotide). ``min_threshold()`` exposes the floor so callers (the
  binary-search wrapper in modes) can stop before hitting it.

Output format is the ``.clstr`` text file:

    >Cluster 0
    0\t1500aa, >seq_A... *
    1\t1499aa, >seq_B... at 99.93%

The ``*`` marks the representative. With ``-d 0`` the full sequence id
appears between ``>`` and ``...``, so the round-trip is exact.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Optional

from ..models import Cluster, Sequence, SequenceType
from .mmseqs2 import _write_id_fasta


class CDHitError(RuntimeError):
    pass


# Word-size selection tables. Each entry: (low_inclusive, high_inclusive, default_n).
# Ranges are taken from the cd-hit user guide (v4.8.1).
_PROTEIN_WORD_SIZES = (
    (0.70, 1.001, 5),
    (0.60, 0.70, 4),
    (0.50, 0.60, 3),
    (0.40, 0.50, 2),
)
_NUCLEOTIDE_WORD_SIZES = (
    (0.95, 1.001, 10),
    (0.90, 0.95, 8),
    (0.88, 0.90, 7),
    (0.85, 0.88, 6),
    (0.80, 0.85, 5),
)

# Identity floors enforced by the binary itself.
_PROTEIN_MIN = 0.40
_NUCLEOTIDE_MIN = 0.80


def _is_protein(sequences: list[Sequence], cfg: Optional[dict[str, Any]] = None) -> bool:
    """Choose the binary by clustering alphabet, else by majority sequence type.

    ``cfg['clustering']['alphabet']`` is consulted first: ``protein`` forces
    ``cd-hit`` regardless of ``seq.seq_type`` (the upstream pipeline may
    populate ``seq.protein_sequence`` on a nucleotide-typed concat record).
    ``nucleotide`` forces ``cd-hit-est`` for the same reason.

    Without a config override we fall back to the alphabet of the records
    themselves: bias to ``cd-hit-est`` only when every sequence is explicitly
    nucleotide, so any ambiguity falls back to the protein binary (which
    has the wider acceptable threshold range).
    """
    if cfg is not None:
        alphabet = cfg.get("clustering", {}).get("alphabet")
        if alphabet == "protein":
            return True
        if alphabet == "nucleotide":
            return False
    return not all(s.seq_type == SequenceType.NUCLEOTIDE for s in sequences)


def _pick_word_size(threshold: float, protein: bool) -> int:
    table = _PROTEIN_WORD_SIZES if protein else _NUCLEOTIDE_WORD_SIZES
    for lo, hi, n in table:
        if lo <= threshold < hi:
            return n
    floor = _PROTEIN_MIN if protein else _NUCLEOTIDE_MIN
    raise CDHitError(
        f"cd-hit identity threshold {threshold} is below the supported "
        f"floor of {floor} for "
        f"{'protein' if protein else 'nucleotide'} clustering."
    )


def min_threshold(
    sequences: list[Sequence],
    cfg: Optional[dict[str, Any]] = None,
) -> float:
    """Return the lowest identity threshold cd-hit will accept for this input."""
    return _PROTEIN_MIN if _is_protein(sequences, cfg) else _NUCLEOTIDE_MIN


def _check_binary(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise CDHitError(
            f"{name} not found in PATH. Install cd-hit from "
            "https://github.com/weizhongli/cdhit"
        )
    return path


def run_clustering(
    sequences: list[Sequence],
    threshold: float,
    cfg: dict[str, Any],
    tmp_dir: Optional[str] = None,
) -> list[Cluster]:
    """Cluster sequences with cd-hit / cd-hit-est and return Cluster objects.

    Mirrors ``mmseqs2.run_clustering``: same call signature, same return
    type, same round-trip sanity check. Sequence IDs are written as the
    sole FASTA header token (via ``_write_id_fasta``) so cd-hit's truncation
    at whitespace cannot corrupt the .clstr round-trip.
    """
    cluster_cfg = cfg.get("clustering", {})
    cdhit_cfg = cluster_cfg.get("cdhit", {}) or {}
    alphabet = cluster_cfg.get("alphabet", "nucleotide")

    protein = _is_protein(sequences, cfg)
    floor = _PROTEIN_MIN if protein else _NUCLEOTIDE_MIN
    if threshold < floor:
        raise CDHitError(
            f"cd-hit ({'cd-hit' if protein else 'cd-hit-est'}) requires "
            f"identity threshold >= {floor}; got {threshold}."
        )

    binary_name = cdhit_cfg.get("binary") or ("cd-hit" if protein else "cd-hit-est")
    binary = _check_binary(binary_name)

    word_size = cdhit_cfg.get("word_size")
    if word_size is None:
        word_size = _pick_word_size(threshold, protein=protein)

    coverage = cdhit_cfg.get("coverage", 0.8)
    # cd-hit: -G 1 (default) = global identity; -G 0 = local. Coverage
    # bounds only apply in local mode, so flip -G accordingly.
    global_alignment = bool(cdhit_cfg.get("global_alignment", True))
    accurate = bool(cdhit_cfg.get("accurate", False))
    memory_mb = int(cdhit_cfg.get("memory_mb", 0))
    extra_args: list[str] = list(cdhit_cfg.get("extra_args", []))

    threads = cfg.get("threads", 4)

    work_dir = Path(tmp_dir or cfg.get("temp_dir", "/tmp/repseq"))
    work_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(dir=work_dir, prefix="cdhit_") as td:
        td = Path(td)
        input_fasta = td / "input.fasta"
        output_prefix = td / "result"

        _write_id_fasta(sequences, input_fasta, alphabet=alphabet)

        cmd = [
            binary,
            "-i", str(input_fasta),
            "-o", str(output_prefix),
            "-c", str(threshold),
            "-n", str(word_size),
            "-T", str(threads),
            "-M", str(memory_mb),
            "-d", "0",  # keep full sequence id in .clstr (no truncation)
            "-G", "1" if global_alignment else "0",
            "-g", "1" if accurate else "0",
        ]
        if not global_alignment:
            cmd.extend(["-aS", str(coverage)])
        cmd.extend(extra_args)

        try:
            subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            raise CDHitError(f"{binary_name} failed:\n{e.stderr}") from e

        clusters = _parse_clstr_file(
            clstr_path=str(output_prefix) + ".clstr",
            sequences=sequences,
        )

    accounted = sum(1 + len(c.members) for c in clusters)
    if accounted != len(sequences):
        raise CDHitError(
            f"Cluster round-trip mismatch: {len(sequences)} sequences in, "
            f"{accounted} accounted for across {len(clusters)} clusters. "
            "Likely an ID/whitespace issue between input FASTA and the "
            ".clstr file."
        )

    return clusters


# Matches a member line: "<index>\t<len>(aa|nt), >SEQID... *" or "...at 99.93%"
# The id is everything between '>' and the final '...'. cd-hit always emits
# the '...' even when -d 0 keeps the full id.
_CLSTR_MEMBER_RE = re.compile(r">(?P<id>.+?)\.\.\.\s*(?P<rep>\*)?")


def _parse_clstr_file(
    clstr_path: str,
    sequences: list[Sequence],
) -> list[Cluster]:
    """Parse a cd-hit ``.clstr`` file into Cluster objects.

    .clstr format:
        >Cluster 0
        0\t1500aa, >seq_A... *
        1\t1499aa, >seq_B... at 99.93%

    The line starting with ``>Cluster`` opens a new cluster; subsequent
    indented lines list members, one of which is marked ``*`` as the
    representative.
    """
    clstr = Path(clstr_path)
    if not clstr.exists():
        raise CDHitError(f"cd-hit cluster file not found: {clstr_path}")

    seq_map: dict[str, Sequence] = {s.id: s for s in sequences}
    for s in sequences:
        if s.accession and s.accession not in seq_map:
            seq_map[s.accession] = s

    raw: list[tuple[Optional[str], list[str]]] = []  # (rep_id, member_ids)
    cur_rep: Optional[str] = None
    cur_members: list[str] = []

    def _flush() -> None:
        if cur_rep is not None or cur_members:
            raw.append((cur_rep, list(cur_members)))

    with open(clstr) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">Cluster"):
                _flush()
                cur_rep = None
                cur_members = []
                continue
            m = _CLSTR_MEMBER_RE.search(line)
            if not m:
                continue
            mid = m.group("id")
            if m.group("rep") == "*":
                cur_rep = mid
            else:
                cur_members.append(mid)
    _flush()

    clusters: list[Cluster] = []
    for i, (rep_id, member_ids) in enumerate(raw):
        if rep_id is None:
            # Defensive: cd-hit always writes a representative, but skip
            # gracefully if a malformed cluster slips through.
            continue
        rep_seq = seq_map.get(rep_id)
        if rep_seq is None:
            continue
        members = [seq_map[m] for m in member_ids if m in seq_map]
        clusters.append(
            Cluster(
                cluster_id=f"cluster_{i+1:06d}",
                representative=rep_seq,
                members=members,
            )
        )

    return clusters
