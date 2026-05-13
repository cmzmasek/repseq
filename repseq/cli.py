"""repseq command-line interface."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import click
import yaml

from .config import load_config, validate_config, get_virus_config
from .io.fasta import read_fasta
from .models import SequenceSource
from .models import RunResult
from .output.report import write_all_reports
from .output.writer import write_results
from .qc.pipeline import run_qc
from .segmented.completeness import (
    build_concatenated_sequences,
    filter_complete_isolates,
)
from .taxonomy.cache import TaxonomyCache
from .taxonomy.ncbi import NCBITaxonomy
from .taxonomy.resolver import MetadataResolver
from .taxonomy.uniprot import UniProtAPI
from .utils.progress import progress_bar


# ---------------------------------------------------------------------------
# Shared options
# ---------------------------------------------------------------------------

def _shared_options(fn):
    fn = click.option("--config", "-c", "config_path", default=None,
                      help="Path to YAML config file.")(fn)
    fn = click.option("--input", "-i", "input_paths", multiple=True, required=True,
                      help="Input FASTA file(s). Repeat for multiple files.")(fn)
    fn = click.option("--output-dir", "-o", default=None,
                      help="Output directory (overrides config).")(fn)
    fn = click.option("--prefix", default=None,
                      help="Output file prefix (overrides config).")(fn)
    fn = click.option("--threads", "-t", default=None, type=int,
                      help="Number of threads (overrides config).")(fn)
    fn = click.option("--seed", default=None, type=int,
                      help="Random seed (overrides config).")(fn)
    fn = click.option("--segmented", is_flag=True, default=False,
                      help="Enable segmented virus mode.")(fn)
    fn = click.option("--dry-run", is_flag=True, default=False,
                      help="Validate inputs and config without running.")(fn)
    fn = click.option("--no-resolve", is_flag=True, default=False,
                      help="Skip database metadata resolution.")(fn)
    fn = click.option("--overflow", type=click.Choice(["keep", "trim"]), default="keep",
                      help="Behaviour when a group exceeds n-per-group.")(fn)
    fn = click.option(
        "--source", "source_override",
        type=click.Choice(["auto", "uniprot", "ncbi", "ncbi_virus"]),
        default="auto",
        help=(
            "Force input source instead of auto-detecting from headers. "
            "Use 'ncbi_virus' for NCBI Virus FASTA downloads, 'uniprot' for "
            "UniProt FASTA, 'ncbi' for standard NCBI nucleotide/protein FASTA. "
            "Default: auto (detect from header format)."
        ),
    )(fn)
    return fn


def _load_and_validate(config_path, output_dir, prefix, threads, seed) -> dict:
    cfg = load_config(config_path)
    if output_dir:
        cfg["output"]["dir"] = output_dir
    if prefix:
        cfg["output"]["prefix"] = prefix
    if threads:
        cfg["threads"] = threads
    if seed is not None:
        cfg["seed"] = seed
    errors = validate_config(cfg)
    if errors:
        for e in errors:
            click.echo(f"[config error] {e}", err=True)
        sys.exit(1)
    return cfg


_SOURCE_MAP = {
    "uniprot": SequenceSource.UNIPROT,
    "ncbi": SequenceSource.NCBI,
    "ncbi_virus": SequenceSource.NCBI_VIRUS,
}


def _load_sequences(input_paths: tuple[str, ...], source_override: str = "auto") -> list:
    override = _SOURCE_MAP.get(source_override)  # None when "auto"
    sequences = []
    for path in input_paths:
        click.echo(f"Reading {path} ...")
        sequences.extend(read_fasta(path, source_override=override))
    click.echo(f"Loaded {len(sequences)} sequences.")
    if override:
        click.echo(f"Source override: {source_override}")
    return sequences


def _resolve_metadata(sequences, cfg, no_resolve):
    if no_resolve:
        return
    tax_cfg = cfg.get("taxonomy", {})
    cache = TaxonomyCache(cfg["cache_dir"], ttl_days=tax_cfg.get("cache_ttl_days", 30))
    ncbi = NCBITaxonomy(cache, email=tax_cfg.get("ncbi_email"), api_key=tax_cfg.get("ncbi_api_key"))
    uniprot = UniProtAPI(cache)
    resolver = MetadataResolver(cache, ncbi, uniprot, threads=cfg.get("threads", 4))
    with progress_bar(len(sequences), "Resolving metadata") as bar:
        resolver.resolve_batch(sequences, progress=bar)
    if resolver.failures:
        click.echo(
            f"  Metadata resolution: {len(resolver.failures)} of {len(sequences)} "
            f"failed (continuing with header-derived metadata).",
            err=True,
        )


def _run_qc(sequences, cfg):
    click.echo("Running QC ...")
    passed, qc_report = run_qc(sequences, cfg)
    click.echo(qc_report.summary())
    return passed, qc_report


def _handle_segmented(sequences, cfg, qc_report):
    virus_cfg = get_virus_config(cfg)
    if not virus_cfg:
        return sequences, None, None
    click.echo("Applying segmented virus completeness filter ...")
    kept, complete_isolates = filter_complete_isolates(sequences, virus_cfg, qc_report)
    click.echo(f"  Complete isolates : {len(complete_isolates)}")
    click.echo(f"  Individual seqs   : {len(kept)}")
    concat_seqs = build_concatenated_sequences(complete_isolates)
    return concat_seqs, complete_isolates, virus_cfg.get("segments")


def _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names):
    out_files = write_results(result, cfg, complete_isolates, segment_names)
    write_all_reports(result, qc_report, cfg, list(input_paths), out_files)
    click.echo(f"\nOutput written to: {cfg['output']['dir']}")
    for f in out_files:
        click.echo(f"  {f.name}")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

@click.group()
@click.version_option(package_name="repseq")
def main():
    """repseq — representative sequence selection for large bioinformatics datasets."""
    pass


# ---------------------------------------------------------------------------
# run global
# ---------------------------------------------------------------------------

@main.command("global")
@_shared_options
@click.option("--threshold", "-T", default=None, type=float,
              help="Identity threshold for clustering (e.g. 0.90).")
@click.option("--n-select", "-n", default=None, type=int,
              help="Number of representative sequences to select.")
def run_global(config_path, input_paths, output_dir, prefix, threads, seed,
               segmented, dry_run, no_resolve, overflow, source_override, threshold, n_select):
    """Global mode: cluster at a threshold or select N diverse sequences."""
    if threshold is None and n_select is None:
        raise click.UsageError("Provide --threshold or --n-select.")

    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed)
    if segmented:
        cfg["segmented"]["enabled"] = True

    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)

    from .modes.global_mode import GlobalMode
    mode = GlobalMode(cfg, threshold=threshold, n_select=n_select)
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names)


# ---------------------------------------------------------------------------
# run taxonomic1
# ---------------------------------------------------------------------------

@main.command("taxonomic1")
@_shared_options
@click.option("--rank", "-r", required=True,
              help="Taxonomic rank to group by (e.g. genus, family).")
@click.option("--n-per-group", "-n", required=True, type=int,
              help="Target representatives per taxonomic group.")
def run_taxonomic1(config_path, input_paths, output_dir, prefix, threads, seed,
                   segmented, dry_run, no_resolve, overflow, source_override, rank, n_per_group):
    """Taxonomic mode 1: N representatives per taxonomic rank group."""
    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)

    from .modes.taxonomic1 import TaxonomicMode1
    mode = TaxonomicMode1(cfg, rank=rank, n_per_group=n_per_group, overflow=overflow)
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names)


# ---------------------------------------------------------------------------
# run taxonomic2
# ---------------------------------------------------------------------------

@main.command("taxonomic2")
@_shared_options
@click.option("--rank-levels", "-r", required=True,
              help='JSON list of {rank, n_per_group} dicts. E.g. \'[{"rank":"family","n_per_group":20},{"rank":"genus","n_per_group":5}]\'')
def run_taxonomic2(config_path, input_paths, output_dir, prefix, threads, seed,
                   segmented, dry_run, no_resolve, overflow, source_override, rank_levels):
    """Taxonomic mode 2: hierarchical multi-rank nested clustering."""
    import json as _json
    try:
        rank_levels_parsed = _json.loads(rank_levels)
    except Exception:
        raise click.UsageError("--rank-levels must be valid JSON.")

    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)

    from .modes.taxonomic2 import TaxonomicMode2
    mode = TaxonomicMode2(cfg, rank_levels=rank_levels_parsed, overflow=overflow)
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names)


# ---------------------------------------------------------------------------
# run host
# ---------------------------------------------------------------------------

@main.command("host")
@_shared_options
@click.option("--n-per-host", "-n", required=True, type=int,
              help="Target representatives per host organism.")
def run_host(config_path, input_paths, output_dir, prefix, threads, seed,
             segmented, dry_run, no_resolve, overflow, source_override, n_per_host):
    """Host-stratified mode: N representatives per host organism."""
    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)

    from .modes.host_mode import HostMode
    mode = HostMode(cfg, n_per_host=n_per_host, overflow=overflow)
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names)


# ---------------------------------------------------------------------------
# run time
# ---------------------------------------------------------------------------

@main.command("time")
@_shared_options
@click.option("--n-per-window", "-n", required=True, type=int,
              help="Target representatives per time window.")
@click.option("--window", default="year",
              help='Time window: "year", "decade", or a number (e.g. "5" for 5-year bins).')
def run_time(config_path, input_paths, output_dir, prefix, threads, seed,
             segmented, dry_run, no_resolve, overflow, source_override, n_per_window, window):
    """Time-stratified mode: N representatives per time window."""
    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)

    from .modes.time_mode import TimeMode
    mode = TimeMode(cfg, n_per_window=n_per_window, window=window, overflow=overflow)
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names)


# ---------------------------------------------------------------------------
# run geographic
# ---------------------------------------------------------------------------

@main.command("geographic")
@_shared_options
@click.option("--n-per-country", "-n", required=True, type=int,
              help="Target representatives per country.")
def run_geographic(config_path, input_paths, output_dir, prefix, threads, seed,
                   segmented, dry_run, no_resolve, overflow, source_override, n_per_country):
    """Geographic mode: N representatives per country."""
    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)

    from .modes.geographic_mode import GeographicMode
    mode = GeographicMode(cfg, n_per_country=n_per_country, overflow=overflow)
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names)


# ---------------------------------------------------------------------------
# run custom
# ---------------------------------------------------------------------------

@main.command("custom")
@_shared_options
@click.option("--field", "-f", required=True,
              help="Field name to group by (attribute, taxonomy rank, or metadata column).")
@click.option("--n-per-group", "-n", required=True, type=int,
              help="Target representatives per group.")
@click.option("--metadata-table", default=None,
              help="Path to TSV/CSV metadata table with accession column.")
@click.option("--field-regex", default=None,
              help="Regex to extract the field value from FASTA headers.")
def run_custom(config_path, input_paths, output_dir, prefix, threads, seed,
               segmented, dry_run, no_resolve, overflow, source_override, field, n_per_group,
               metadata_table, field_regex):
    """Custom metadata mode: group by any field or metadata table column."""
    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)

    from .modes.custom_mode import CustomMode
    mode = CustomMode(
        cfg, field=field, n_per_group=n_per_group,
        metadata_table_path=metadata_table,
        field_regex=field_regex,
        overflow=overflow,
    )
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names)


# ---------------------------------------------------------------------------
# run hybrid
# ---------------------------------------------------------------------------

@main.command("hybrid")
@_shared_options
@click.option("--fields", "-f", required=True,
              help="Comma-separated list of fields to combine (e.g. genus,host,decade).")
@click.option("--n-per-group", "-n", required=True, type=int,
              help="Target representatives per stratum.")
@click.option("--metadata-table", default=None,
              help="Path to TSV/CSV metadata table with accession column.")
def run_hybrid(config_path, input_paths, output_dir, prefix, threads, seed,
               segmented, dry_run, no_resolve, overflow, source_override, fields, n_per_group,
               metadata_table):
    """Hybrid mode: multi-dimensional stratification (e.g. genus × host × year)."""
    field_list = [f.strip() for f in fields.split(",")]
    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)

    from .modes.hybrid_mode import HybridMode
    mode = HybridMode(
        cfg, fields=field_list, n_per_group=n_per_group,
        overflow=overflow, metadata_table_path=metadata_table,
    )
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names)


# ---------------------------------------------------------------------------
# cache management
# ---------------------------------------------------------------------------

@main.group("cache")
def cache_group():
    """Manage the taxonomy/metadata cache."""
    pass


@cache_group.command("stats")
@click.option("--config", "-c", "config_path", default=None)
def cache_stats(config_path):
    """Show cache statistics."""
    cfg = load_config(config_path)
    cache = TaxonomyCache(cfg["cache_dir"])
    stats = cache.stats()
    click.echo(f"Cache: {stats['db_path']}")
    click.echo(f"Total entries : {stats['total_entries']}")
    click.echo(f"DB size       : {stats['db_size_mb']} MB")
    for src, count in stats.get("by_source", {}).items():
        click.echo(f"  {src}: {count}")


@cache_group.command("clear")
@click.option("--config", "-c", "config_path", default=None)
@click.option("--source", default=None, help="Clear only a specific source.")
@click.confirmation_option(prompt="This will delete cached data. Continue?")
def cache_clear(config_path, source):
    """Clear the cache (all sources or a specific one)."""
    cfg = load_config(config_path)
    cache = TaxonomyCache(cfg["cache_dir"])
    n = cache.clear(source)
    click.echo(f"Deleted {n} cache entries.")


@cache_group.command("purge-expired")
@click.option("--config", "-c", "config_path", default=None)
def cache_purge(config_path):
    """Remove expired cache entries."""
    cfg = load_config(config_path)
    cache = TaxonomyCache(cfg["cache_dir"])
    n = cache.purge_expired()
    click.echo(f"Purged {n} expired cache entries.")


# ---------------------------------------------------------------------------
# Config wizard
# ---------------------------------------------------------------------------

@main.command("init-config")
@click.option("--output", "-o", default="repseq_config.yaml",
              help="Path to write the generated config file.")
def init_config(output):
    """Interactive wizard to generate a repseq YAML config file."""
    click.echo("repseq config wizard\n")

    cfg = {}

    cfg["cache_dir"] = click.prompt("Cache directory", default="~/.repseq/cache")
    cfg["temp_dir"] = click.prompt("Temp directory", default="/tmp/repseq")
    cfg["threads"] = click.prompt("Number of threads", default=4, type=int)
    cfg["seed"] = click.prompt("Random seed", default=42, type=int)

    # QC
    cfg["qc"] = {}
    cfg["qc"]["remove_duplicates"] = click.confirm("Remove exact duplicates?", default=True)

    lf_mode = click.prompt(
        "Length filter mode", type=click.Choice(["median_percent", "min_max"]),
        default="median_percent"
    )
    cfg["qc"]["length_filter"] = {"mode": lf_mode}
    if lf_mode == "median_percent":
        cfg["qc"]["length_filter"]["min_percent"] = click.prompt(
            "Minimum percent of median length", default=50, type=int
        )
    else:
        cfg["qc"]["length_filter"]["min_length"] = click.prompt("Minimum length", default=None, type=int)
        cfg["qc"]["length_filter"]["max_length"] = click.prompt("Maximum length", default=None, type=int)

    cfg["qc"]["ambiguous_threshold"] = click.prompt(
        "Max ambiguous character fraction (0-1)", default=0.05, type=float
    )

    cfg["qc"]["annotation_filter"] = {
        "enabled": click.confirm("Enable annotation keyword filter?", default=True),
        "keywords": [
            "MAG:", "metagenome-assembled", "synthetic", "artificial",
            "fragment", "partial", "environmental sample",
            "uncultured", "unclassified", "unidentified", "hypothetical",
        ],
    }

    # Segmented
    if click.confirm("\nConfigure segmented virus mode?", default=False):
        virus_name = click.prompt("Virus name (e.g. influenza_a)")
        n_seg = click.prompt("Number of expected segments", type=int)
        segment_list = click.prompt(
            f"Segment names in order (comma-separated, {n_seg} expected)"
        )
        segments = [s.strip() for s in segment_list.split(",")]
        isolate_regex = click.prompt(
            "Regex to extract isolate ID from header",
            default=r"(?P<isolate>[A-Za-z]+/[^/]+/[^/]+/[^/]+/\d{4})"
        )
        cfg["segmented"] = {
            "enabled": False,
            "virus": virus_name,
            "viruses": {
                virus_name: {
                    "expected_segments": n_seg,
                    "segments": segments,
                    "isolate_regex": isolate_regex,
                }
            },
        }
    else:
        cfg["segmented"] = {"enabled": False}

    # Taxonomy
    cfg["taxonomy"] = {}
    ncbi_email = click.prompt("NCBI email (for Entrez API, leave blank to skip)", default="")
    if ncbi_email:
        cfg["taxonomy"]["ncbi_email"] = ncbi_email
    ncbi_key = click.prompt("NCBI API key (optional, leave blank to skip)", default="")
    if ncbi_key:
        cfg["taxonomy"]["ncbi_api_key"] = ncbi_key
    cfg["taxonomy"]["cache_ttl_days"] = click.prompt("Cache TTL (days)", default=30, type=int)

    # Clustering
    cfg["clustering"] = {
        "backend": "mmseqs2",
        "mmseqs2_mode": click.prompt(
            "MMseqs2 mode", type=click.Choice(["easy-linclust", "easy-cluster"]),
            default="easy-linclust"
        ),
        "coverage": click.prompt("Coverage threshold", default=0.8, type=float),
        "coverage_mode": click.prompt("Coverage mode (0-4)", default=0, type=int),
    }

    # Output
    cfg["output"] = {
        "dir": click.prompt("Output directory", default="./repseq_output"),
        "prefix": click.prompt("Output file prefix", default="repseq"),
    }

    out_path = Path(output)
    with open(out_path, "w") as fh:
        yaml.dump(cfg, fh, default_flow_style=False, sort_keys=False)

    click.echo(f"\nConfig written to: {out_path}")


if __name__ == "__main__":
    main()
