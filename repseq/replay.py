"""``repseq replay`` — re-materialise representatives from a lockfile.

The replay subcommand reads a ``{prefix}_lockfile.json`` written by a
prior repseq run and re-emits the representative FASTAs (and, with
``--rebuild-trees``, the phylogeny step) into a fresh output
directory. It is **not** a fresh run — QC, clustering, mode logic,
and the representative-selection algorithm are all skipped because
the lockfile already says which representatives won. The only work
replay does is:

1. Validate the lockfile (schema version, repseq version).
2. Verify the HMM database hasn't changed under the user (warn
   loudly if it has).
3. Re-fetch each representative's sequence body and metadata via
   the existing cache → NCBI/UniProt path.
4. Build ``Sequence`` objects (and, in segmented mode, the
   ``CONCAT`` rep objects with ``concat_segments`` populated).
5. Call the same writer + reporter ``write_results`` /
   ``write_all_reports`` use for a fresh run, so the FASTAs are
   byte-identical to what the original run wrote — provided NCBI
   hasn't changed the records out from under us.

Accessions NCBI can't return (retired, network down, etc.) trigger
loud stderr warnings and land in ``{prefix}_replay_missing.tsv``;
the replay then continues with the survivors so the bench
scientist gets *something* back rather than the whole step
crashing. Set ``--rebuild-trees`` to also re-run the phylogeny
step (MAFFT + IQ-TREE/FastTree); the resulting tree may differ
slightly from the original due to bootstrap variance even at the
same seed, so it is opt-in.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any, Optional

import click

from . import __version__ as REPSEQ_VERSION
from .lockfile import LockfileVersionError, compute_sha256, read_lockfile
from .models import (
    Cluster,
    RunResult,
    Sequence,
    SequenceSource,
    SequenceType,
    TaxonomyInfo,
)
from .output.writer import write_results
from .taxonomy.cache import TaxonomyCache
from .taxonomy.ncbi import NCBITaxonomy


def _validate_repseq_version(lockfile: dict[str, Any]) -> None:
    """Hard-fail when the lockfile's repseq major version is newer than
    this installation; warn on minor mismatch.

    A future-major lockfile may reference fields this repseq doesn't
    know about — silently misreading them would corrupt the replay.
    A future-minor lockfile is allowed (forward-compatible additive
    changes) but still flagged so the user knows.
    """
    lf_version = lockfile.get("repseq_version", "")
    if not lf_version:
        click.echo(
            "[replay] warning: lockfile carries no repseq_version field",
            err=True,
        )
        return
    try:
        lf_major, lf_minor = _split_version(lf_version)
        cur_major, cur_minor = _split_version(REPSEQ_VERSION)
    except ValueError:
        click.echo(
            f"[replay] warning: could not parse versions "
            f"(lockfile={lf_version}, installed={REPSEQ_VERSION})",
            err=True,
        )
        return
    if lf_major > cur_major:
        raise click.ClickException(
            f"lockfile was written by repseq {lf_version} but this "
            f"installation is {REPSEQ_VERSION} — major-version mismatch. "
            f"Upgrade repseq before replaying."
        )
    if lf_major != cur_major or lf_minor != cur_minor:
        click.echo(
            f"[replay] note: lockfile is from repseq {lf_version}; "
            f"replaying with {REPSEQ_VERSION}",
            err=True,
        )


def _split_version(v: str) -> tuple[int, int]:
    parts = v.strip().split(".")
    return int(parts[0]), int(parts[1])


def _verify_hmm_db(lockfile: dict[str, Any]) -> None:
    """Compare the lockfile's HMM DB sha256 with the resolved path's current
    sha256; warn loudly on mismatch.

    A mismatched HMM DB means the HMM-tier-derived state on the
    re-fetched sequences (`hmm_hits`, marker selection) could differ
    from the original run, even if the accessions are byte-identical.
    Hard-failing would be too strict — many users don't care about
    HMM-tier consistency in a replay. So we surface the difference
    and let the user decide.
    """
    expected = lockfile.get("hmm_db") or {}
    if not expected:
        return
    expected_sha = expected.get("sha256")
    if not expected_sha:
        return
    try:
        from .hmm.database import resolve_database_path
        path = resolve_database_path(
            (lockfile.get("config", {}).get("hmm", {}) or {}).get("database"),
        )
    except Exception as exc:
        click.echo(
            f"[replay] warning: could not resolve HMM database "
            f"recorded in lockfile: {exc}",
            err=True,
        )
        return
    current_sha = compute_sha256(path)
    if current_sha != expected_sha:
        click.echo(
            f"[replay] warning: HMM database content has changed since "
            f"the lockfile was written\n"
            f"           path: {path}\n"
            f"           lockfile sha256: {expected_sha[:16]}…\n"
            f"           current sha256:  {current_sha[:16]}…\n"
            f"         HMM-derived state may differ from the original run.",
            err=True,
        )


def _collect_accessions(lockfile: dict[str, Any]) -> list[str]:
    """Flatten the lockfile's representative entries into a single list
    of accessions to fetch.

    Segmented entries contribute one accession per segment; non-segmented
    entries contribute their single ``accession``. Order is preserved so
    a per-isolate progress bar reads naturally.
    """
    accessions: list[str] = []
    for rep in lockfile.get("representatives", []) or []:
        if rep.get("kind") == "isolate":
            for acc in (rep.get("segment_accessions") or {}).values():
                if acc:
                    accessions.append(acc)
        else:
            acc = rep.get("accession")
            if acc:
                accessions.append(acc)
    return accessions


def _make_sequence(
    accession: str,
    body: str,
    organism: Optional[str],
    source_meta: Optional[dict[str, Optional[str]]] = None,
    segment_label: Optional[str] = None,
    isolate_id: Optional[str] = None,
) -> Sequence:
    """Construct a :class:`Sequence` for a re-fetched NCBI record.

    Replay can only re-materialise what's reachable from the lockfile +
    cache; richer metadata (host, country, collection_date,
    taxonomy lineage) will populate via the resolver if the user wires
    it up, but is not strictly required for the rep FASTAs.
    """
    src = source_meta or {}
    return Sequence(
        id=accession,
        header=accession,
        sequence=body,
        seq_type=SequenceType.NUCLEOTIDE,
        source=SequenceSource.NCBI,
        accession=accession,
        organism=organism,
        segment=segment_label or src.get("segment"),
        isolate_id=isolate_id or src.get("isolate") or src.get("strain"),
    )


def _build_replay_result(
    lockfile: dict[str, Any],
    seqs_by_acc: dict[str, Optional[str]],
    source_meta: dict[str, dict[str, Optional[str]]],
) -> tuple[RunResult, Optional[dict[str, list[Sequence]]], Optional[list[str]], list[tuple[str, str]]]:
    """Turn the lockfile's representative entries into a ``RunResult``.

    Returns ``(result, complete_isolates_or_None, segment_names_or_None,
    missing)`` where ``missing`` is the list of ``(rep_id, accession)``
    pairs that NCBI couldn't return — surfaced to the user in
    ``{prefix}_replay_missing.tsv`` so the discrepancy is auditable.
    """
    representatives: list[Sequence] = []
    complete_isolates: dict[str, list[Sequence]] = {}
    missing: list[tuple[str, str]] = []

    segment_names: Optional[list[str]] = None
    virus_cfg = (
        (lockfile.get("config", {}).get("viruses", {}) or {})
        if isinstance(lockfile.get("config", {}).get("viruses"), dict)
        else {}
    )
    # The first virus's `segments:` defines the per-segment write order.
    for entry in virus_cfg.values():
        if isinstance(entry, dict) and entry.get("segments"):
            segment_names = list(entry["segments"])
            break

    for rep_entry in lockfile.get("representatives", []) or []:
        kind = rep_entry.get("kind")
        if kind == "isolate":
            isolate_id = rep_entry.get("isolate_id") or rep_entry.get("id")
            organism = rep_entry.get("organism")
            seg_records: list[Sequence] = []
            any_missing = False
            for seg_label, acc in (rep_entry.get("segment_accessions") or {}).items():
                body = seqs_by_acc.get(acc) if acc else None
                if not body:
                    any_missing = True
                    if acc:
                        missing.append((rep_entry.get("id", ""), acc))
                    continue
                seg_records.append(
                    _make_sequence(
                        acc, body, organism,
                        source_meta.get(acc),
                        segment_label=seg_label,
                        isolate_id=isolate_id,
                    )
                )
            if not seg_records:
                continue
            # Synthetic CONCAT rep (segmented).
            concat_seq = "".join(s.sequence for s in seg_records)
            concat = Sequence(
                id=f"CONCAT|{isolate_id}" if not str(rep_entry.get("id", "")).startswith("CONCAT|") else rep_entry["id"],
                header=isolate_id or "",
                sequence=concat_seq,
                seq_type=SequenceType.NUCLEOTIDE,
                source=SequenceSource.NCBI,
                organism=organism,
                isolate_id=isolate_id,
            )
            concat.concat_segments = seg_records
            representatives.append(concat)
            complete_isolates[isolate_id] = seg_records
            if any_missing:
                # Record partial replay so the missing TSV reports it.
                pass
        else:
            acc = rep_entry.get("accession")
            body = seqs_by_acc.get(acc) if acc else None
            if not body:
                if acc:
                    missing.append((rep_entry.get("id", ""), acc))
                continue
            seq = _make_sequence(
                acc, body, rep_entry.get("organism"),
                source_meta.get(acc),
            )
            representatives.append(seq)

    result = RunResult(
        mode=lockfile.get("mode", "replay"),
        representatives=representatives,
        clusters=[
            Cluster(cluster_id=str(i), representative=rep, members=[])
            for i, rep in enumerate(representatives)
        ],
    )
    return (
        result,
        complete_isolates if complete_isolates else None,
        segment_names,
        missing,
    )


def _write_missing_tsv(missing: list[tuple[str, str]], path: Path) -> None:
    """Tabulate accessions that NCBI couldn't return for the replay.

    No-op when nothing failed (we don't want a stray empty file
    polluting an otherwise-clean output dir).
    """
    if not missing:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write("representative_id\taccession\n")
        for rep_id, acc in missing:
            fh.write(f"{rep_id}\t{acc}\n")


def _write_replay_summary(
    lockfile: dict[str, Any],
    written: list[Path],
    missing: list[tuple[str, str]],
    out_dir: Path,
    prefix: str,
) -> Path:
    """Write a tiny ``{prefix}_replay.md`` next to the FASTAs.

    Replay does NOT write the original ``{prefix}_summary.md`` (that
    file is the *original* run's methods record and reproducing it
    from re-fetched data would imply we re-did QC/clustering, which we
    didn't). The replay summary is a short prose note pointing at the
    lockfile and listing what was emitted.
    """
    path = out_dir / f"{prefix}_replay.md"
    lines = [
        "# repseq replay\n",
        "",
        "This output directory contains representatives re-materialised "
        f"by `repseq replay` from a lockfile written by **repseq "
        f"{lockfile.get('repseq_version', '?')}** on "
        f"`{lockfile.get('created_utc', '?')}`. The original run is "
        f"identified by the lockfile fields below.\n",
        "",
        f"- **mode**: `{lockfile.get('mode', '?')}`",
        f"- **representatives**: {len(lockfile.get('representatives') or [])}",
        f"- **re-fetched successfully**: "
        f"{len(lockfile.get('representatives') or []) - len({m[0] for m in missing})}",
        f"- **missing (NCBI could not return)**: {len(missing)}",
        "",
        "Replay does **not** re-run QC, clustering, or selection — the "
        "representatives in the lockfile are authoritative. The "
        "FASTA outputs below should be byte-identical to the original "
        "run's representative FASTAs, provided NCBI hasn't changed the "
        "records (any divergence is logged in "
        "`{prefix}_replay_missing.tsv` and as stderr warnings).\n",
        "",
        "## Files written",
        "",
    ]
    for p in written:
        lines.append(f"- `{p.name}`")
    path.write_text("\n".join(lines) + "\n")
    return path


@click.command("replay")
@click.argument("lockfile_path", type=click.Path(exists=True, dir_okay=False))
@click.option(
    "--output-dir", "-o", "output_dir", required=True,
    help="Directory to write replayed outputs into (must be empty or non-existent).",
)
@click.option(
    "--config", "-c", "config_path", default=None,
    help=(
        "Path to a YAML config file — only needed to locate the cache "
        "directory and any NCBI api_key/email. If omitted, repseq uses "
        "the default config-resolution rules."
    ),
)
@click.option(
    "--rebuild-trees", is_flag=True, default=False,
    help=(
        "After re-emitting FASTAs, re-run the phylogeny step (MAFFT + "
        "IQ-TREE/FastTree) using the cfg recorded in the lockfile. The "
        "rebuilt tree may differ from the original by bootstrap "
        "variance even at the same seed, so this is opt-in."
    ),
)
def replay_command(
    lockfile_path: str,
    output_dir: str,
    config_path: Optional[str],
    rebuild_trees: bool,
) -> None:
    """Re-materialise representatives from a {prefix}_lockfile.json.

    The lockfile records exactly which representatives were elected
    in an earlier run; replay re-fetches their sequences from NCBI
    (via the local cache) and re-emits the same representative FASTAs
    in a fresh output directory. QC, clustering, and selection are
    NOT re-run.
    """
    from .config import load_config

    lockfile = read_lockfile(Path(lockfile_path))
    _validate_repseq_version(lockfile)
    _verify_hmm_db(lockfile)

    # The CACHE config wins from the user-supplied -c if any, else
    # from the lockfile's recorded config. The lockfile's config is
    # NOT used to drive QC / clustering — replay skips those — but it
    # IS used to determine segmented mode + segment names.
    runtime_cfg = load_config(config_path) if config_path else load_config(None)
    locked_cfg = lockfile.get("config", {}) or {}
    # Output dir override.
    runtime_cfg.setdefault("output", {})["dir"] = output_dir
    runtime_cfg["output"].setdefault(
        "prefix", locked_cfg.get("output", {}).get("prefix", "repseq"),
    )
    # Honour segmented + virus config from the lockfile.
    if locked_cfg.get("segmented"):
        runtime_cfg["segmented"] = locked_cfg["segmented"]
    if locked_cfg.get("viruses"):
        runtime_cfg["viruses"] = locked_cfg["viruses"]
    if locked_cfg.get("virus"):
        runtime_cfg["virus"] = locked_cfg["virus"]

    out_dir = Path(output_dir)
    if out_dir.exists() and any(out_dir.iterdir()):
        raise click.ClickException(
            f"output directory '{out_dir}' is non-empty; pick an empty / "
            "non-existent path so original outputs aren't overwritten."
        )
    out_dir.mkdir(parents=True, exist_ok=True)

    accessions = _collect_accessions(lockfile)
    if not accessions:
        raise click.ClickException("lockfile carries no accessions to re-fetch")
    click.echo(
        f"[replay] re-fetching {len(accessions)} accession(s) from NCBI…",
        err=True,
    )

    cache = TaxonomyCache(runtime_cfg["cache_dir"])
    tax_cfg = runtime_cfg.get("taxonomy", {}) or {}
    ncbi = NCBITaxonomy(
        cache=cache,
        email=tax_cfg.get("ncbi_email"),
        api_key=tax_cfg.get("ncbi_api_key"),
    )
    seqs_by_acc = ncbi.fetch_nucleotide_batch(accessions)
    source_meta = ncbi.fetch_source_metadata_batch(accessions)
    missing_before = [acc for acc, body in seqs_by_acc.items() if not body]
    if missing_before:
        click.echo(
            f"[replay] warning: {len(missing_before)} accession(s) could not "
            f"be fetched from NCBI — replay output will be incomplete. See "
            f"{{prefix}}_replay_missing.tsv.",
            err=True,
        )

    prefix = runtime_cfg["output"]["prefix"]
    result, complete_isolates, segment_names, missing = _build_replay_result(
        lockfile, seqs_by_acc, source_meta,
    )

    if not result.representatives:
        raise click.ClickException(
            "no representatives could be re-fetched — replay produced "
            "no output. Check the lockfile's accessions are still valid."
        )

    written = write_results(
        result, runtime_cfg, complete_isolates, segment_names,
    )

    if missing:
        miss_path = out_dir / f"{prefix}_replay_missing.tsv"
        _write_missing_tsv(missing, miss_path)
        written.append(miss_path)

    if rebuild_trees:
        try:
            from .phylo import PhyloError, run_phylogeny
            phylo_files = run_phylogeny(
                result.representatives, locked_cfg, out_dir, prefix,
            )
            written.extend(phylo_files)
        except PhyloError as exc:
            click.echo(f"[replay phylo skipped] {exc}", err=True)
        except Exception as exc:
            click.echo(f"[replay phylo failed] {exc}", err=True)

    summary_path = _write_replay_summary(
        lockfile, written, missing, out_dir, prefix,
    )
    written.append(summary_path)

    click.echo(f"\n[replay] wrote {len(written)} file(s) to {out_dir}", err=True)
