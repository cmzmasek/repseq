"""repseq command-line interface."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional

import click
import yaml

from . import __version__
from .config import load_config, validate_config, get_virus_config
from .errors import InputError, RepseqError
from .io.fasta import read_fasta
from .models import SequenceSource
from .models import RunResult
from .overrides import ProtectionPolicy, protected_keep, resolve_ids, resolve_stages
from .output.report import (
    write_all_reports,
    write_nucleotide_taxonomic_report,
    write_nucleotide_taxonomic_report_tsv,
    write_polyprotein_taxonomic_report,
    write_polyprotein_taxonomic_report_tsv,
    write_protein_taxonomic_report,
    write_protein_taxonomic_report_tsv,
    write_taxonomic_report,
    write_taxonomic_report_tsv,
)
from .output.writer import write_results
from .clustering.marker import populate_protein_sequences
from .qc.pipeline import remove_duplicates, run_qc
from .qc.protein_qc import _stderr_batch_progress, attach_proteins, run_protein_qc
from .segmented.completeness import (
    build_concatenated_sequences,
    detect_strain_collisions,
    filter_complete_isolates,
    segment_length_filter,
)
from .segmented.taxonomy_consistency import filter_taxonomy_consistent_isolates
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
        "--plot", is_flag=True, default=False,
        help=(
            "Render a 2-D scatter of the clustering result to "
            "{prefix}_clustering.png. Uses UMAP when 'umap-learn' is "
            "installed (the [viz-umap] extra) and falls back to a "
            "numpy-only classical MDS (PCoA) on a k-mer Jaccard distance "
            "otherwise; the chosen method appears in the figure title. "
            "Requires the [viz] extras for matplotlib: pip install "
            "'repseq[viz]' (or 'repseq[viz-umap]' for the UMAP upgrade)."
        ),
    )(fn)
    fn = click.option(
        "--phylo", is_flag=True, default=False,
        help=(
            "After selection, build an MSA (MAFFT) and an ML tree over the "
            "representatives, then write {prefix}_msa.fasta, "
            "{prefix}_tree.xml (phyloXML with taxonomy + "
            "metadata properties and taxonomy-driven leaf colours), and "
            "{prefix}_tree_id_map.tsv. The plain-text {prefix}_tree.nwk is "
            "off by default (phylo.newick / --newick to keep it). Tree "
            "builder is chosen by phylo.tool "
            "(default 'auto': IQ-TREE with ModelFinder + UFBoot for protein "
            "alignments, FastTree -nt -gtr for nucleotide). For protein + "
            "IQ-TREE with >=2 HMM-resolvable marker families, builds a "
            "partitioned supermatrix by default (phylo.partition). Optional "
            "trimAl trimming via phylo.trimal (off by default). Requires "
            "'mafft' on PATH, plus 'iqtree2' and/or 'FastTree' as needed; "
            "skipped with a warning if fewer than 3 representatives."
        ),
    )(fn)
    fn = click.option(
        "--per-protein-phylo", is_flag=True, default=False,
        help=(
            "Build one tree per declared protein family AND polyprotein "
            "peptide over the representatives. Three groups, three "
            "distinct subdirectories: cluster_protein / segment_markers "
            "trees → {prefix}_per_protein/; extra_protein trees → "
            "{prefix}_extra_protein/; polyprotein peptide trees → "
            "{prefix}_polyprotein/ (alongside the peptide FASTAs the "
            "polyprotein cutter wrote). Requires the HMM tier "
            "(hmm.enabled + configured hmms:) and 'mafft' plus a tree "
            "builder; each family/peptide needs >= "
            "phylo.per_protein.min_taxa representatives. All three groups "
            "share the same MAFFT / IQ-TREE-or-FastTree / rooting / LCA "
            "/ colour-palette settings under phylo.per_protein. Runs "
            "alone or alongside --phylo."
        ),
    )(fn)
    fn = click.option(
        "--per-segment-phylo", is_flag=True, default=False,
        help=(
            "Segmented mode only. Build one nucleotide tree per segment "
            "over the representative isolates, into the "
            "{prefix}_per_segment/ subdirectory. Each isolate contributes "
            "its raw per-segment NT (from concat_segments). Complements "
            "--per-protein-phylo: reassortment can show up as topological "
            "incongruence between the per-segment trees themselves, even "
            "when no single marker tree captures it. Needs 'mafft' plus "
            "a tree builder; each segment needs >= phylo.per_protein.min_taxa "
            "representatives. Runs alone or alongside --phylo / --per-protein-phylo."
        ),
    )(fn)
    fn = click.option(
        "--newick/--no-newick", "newick", default=None,
        help=(
            "Keep the plain-text Newick (*_tree.nwk) for every tree built "
            "this run (whole-genome, per-protein / extra / segment, "
            "pre-cluster, partition). OFF by default to reduce output "
            "clutter — the annotated phyloXML (*_tree.xml) is a topological "
            "superset and the *_tree_id_map.tsv that decodes the retained "
            "*_msa.fasta is kept regardless. The Newick is still built "
            "internally (the phyloXML is parsed from it; the incongruence "
            "table reads it), then dropped at the end unless kept. "
            "Overrides phylo.newick in the YAML; --no-newick forces it off."
        ),
    )(fn)
    fn = click.option(
        "--fast", "fast", is_flag=True, default=False,
        help=(
            "Preliminary-run mode for fast tree building. Overrides the "
            "YAML phylo: settings and forces, for every tree built in "
            "this run (whole-genome 2E AND per-protein 2F): "
            "FastTree (skips IQ-TREE / ModelFinder / UFBoot), MAFFT "
            "'--retree 1' single-pass progressive alignment (drops "
            "'--auto' and any L-INS-i flags), no trimAl trimming, no "
            "partitioned supermatrix, and midpoint rooting (skips the "
            "taxonomy-guided / MAD rooting chain). Use this to iterate "
            "quickly during config tuning; switch it off for the final "
            "publication run."
        ),
    )(fn)
    fn = click.option(
        "--verbose", "verbose", is_flag=True, default=False,
        help=(
            "Stream the live stderr (heartbeat) of MAFFT, IQ-TREE and "
            "FastTree to the terminal as it arrives. Without --verbose "
            "only one start line and one finish line per subprocess "
            "appears, so the log stays compact. The on-failure error "
            "message still carries the full buffered stderr either way."
        ),
    )(fn)
    fn = click.option(
        "--pre-cluster-tree", "pre_cluster_tree", is_flag=True,
        default=False,
        help=(
            "Build a rough overview tree of EVERY post-QC sequence "
            "(one leaf per CONCAT isolate in segmented mode) BEFORE "
            "clustering, with representative leaves prefixed '[repr] ' "
            "in the phyloXML <name> so you can see at a glance where "
            "the elected reps land in the broader diversity. Pipeline "
            "is hard-coded for speed regardless of the rest of phylo: "
            "MAFFT --retree 1, FastTree, midpoint root only, no LCA, "
            "no trimAl, no bootstrap. Outputs: "
            "{prefix}_pre_cluster_tree.xml + _tree_id_map.tsv "
            "(short_id\\taccession\\tis_rep); the .nwk is kept only with "
            "phylo.newick / --newick. Can also be enabled via "
            "phylo.pre_cluster_tree.enabled: true in the YAML."
        ),
    )(fn)
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
    fn = click.option(
        "--alphabet-for-clustering",
        type=click.Choice(["protein", "nucleotide"]),
        default=None,
        help=(
            "Alphabet fed to the clustering backend (overrides config). "
            "'protein' clusters on amino acid sequences (default; recommended "
            "for diverged virus families); 'nucleotide' clusters on raw NT. "
            "ONLY affects clustering — GenBank CDS download, protein-count QC, "
            "and the per-segment expected-proteins check still run on every "
            "isolate regardless of this value."
        ),
    )(fn)
    fn = click.option(
        "--concatenate-markers/--no-concatenate-markers",
        "concatenate_markers",
        default=None,
        help=(
            "Non-segmented protein clustering: cluster on the concatenation "
            "of the marker CDS from EVERY cluster_protein spec (in declared "
            "order, e.g. Spike+Nucleocapsid) instead of the single first "
            "matching marker. Overrides clustering.concatenate_markers. A "
            "record missing any required marker is dropped. No effect in "
            "segmented or nucleotide mode. Does NOT change the whole-genome "
            "tree — use phylo.tool: auto/iqtree for a partitioned multi-gene "
            "tree. Default: from config (off)."
        ),
    )(fn)
    fn = click.option(
        "--protect-ids", "protect_ids", default=None,
        help=(
            "Path to a file of accessions / isolate-ids of special "
            "importance (one per line; '#' comments and blank lines "
            "ignored). These sequences are FORCE-KEPT through the QC "
            "removal stages they would otherwise fail; the file is unioned "
            "with any overrides.ids / overrides.ids_file in the config and "
            "turns overrides.protect_qc ON for the run. Matching is by "
            "accession (non-segmented) / isolate_id (segmented), case- and "
            "version-insensitive (NC_045512 matches NC_045512.2). Rescued "
            "records — and the reason each would have been dropped — are "
            "written to {prefix}_overrides.tsv. To protect only specific QC "
            "stages (not all), set overrides.protect_stages in the config. "
            "Note: force-keep guarantees a sequence survives QC, NOT that it "
            "becomes a representative (use --pin-ids for that)."
        ),
    )(fn)
    fn = click.option(
        "--pin-ids", "pin_ids", default=None,
        help=(
            "Path to a file of accessions / isolate-ids to FORCE-SELECT as "
            "representatives (one per line; '#' comments and blank lines "
            "ignored). Guarantees each named sequence appears in the output "
            "regardless of clustering: a pinned sequence wins its cluster's "
            "representative slot, pins colliding in one cluster are split "
            "into singletons, and diversity-deselected pins are added as "
            "singletons. The file is unioned with any overrides.ids / "
            "overrides.ids_file in the config and turns overrides.force_select "
            "ON for the run. Audit in {prefix}_force_selected.tsv. Force-"
            "select needs the sequence to survive QC first — combine with "
            "--protect-ids (or list the same ids under both) for a 'present "
            "no matter what' guarantee. Matching is by accession / "
            "isolate_id, case- and version-insensitive."
        ),
    )(fn)
    return fn


def _load_and_validate(config_path, output_dir, prefix, threads, seed,
                       alphabet_for_clustering=None, concatenate_markers=None,
                       fast=False, verbose=False, protect_ids=None,
                       pin_ids=None, newick=None) -> dict:
    cfg = load_config(config_path)
    if output_dir:
        cfg["output"]["dir"] = output_dir
    if prefix:
        cfg["output"]["prefix"] = prefix
    if threads:
        cfg["threads"] = threads
    if seed is not None:
        cfg["seed"] = seed
    if alphabet_for_clustering is not None:
        cfg.setdefault("clustering", {})["alphabet_for_clustering"] = alphabet_for_clustering
    if concatenate_markers is not None:
        cfg.setdefault("clustering", {})["concatenate_markers"] = bool(concatenate_markers)
    if newick is not None:
        cfg.setdefault("phylo", {})["newick"] = bool(newick)
    # Stash --verbose on the cfg so the MAFFT / IQ-TREE / FastTree
    # wrappers can read it via cfg.get("verbose", False) without having
    # to thread an extra arg through every call site.
    cfg["verbose"] = bool(verbose)
    # --protect-ids / --pin-ids: merge each CLI file into overrides.ids and
    # turn the matching flag on for the run (an explicit request, so it
    # overrides a config that left the flag off). Merged before validation
    # so the injected ids are validated too. Both feed the one id list, per
    # the "one list, two independent flags" model.
    _merge_override_id_file(cfg, protect_ids, "--protect-ids", "protect_qc")
    _merge_override_id_file(cfg, pin_ids, "--pin-ids", "force_select")
    errors = validate_config(cfg)
    if errors:
        for e in errors:
            click.echo(f"[config error] {e}", err=True)
        sys.exit(1)
    if fast:
        _apply_fast_overrides(cfg)
    _resolve_overrides(cfg)
    _check_output_dir(cfg)
    return cfg


def _merge_override_id_file(cfg: dict, path, flag_name: str, enable_key: str) -> None:
    """Union an id file from a CLI flag into overrides.ids and flip a flag.

    Shared by ``--protect-ids`` (enable_key=``protect_qc``) and ``--pin-ids``
    (enable_key=``force_select``). No-op when ``path`` is falsy; friendly
    exit on a missing file.
    """
    if not path:
        return
    from pathlib import Path as _Path

    if not _Path(path).expanduser().is_file():
        click.echo(
            f"[error] {flag_name} file not found: {path}. Provide a "
            f"readable file with one accession/isolate-id per line.",
            err=True,
        )
        sys.exit(1)
    from .overrides import load_ids_file

    ov = cfg.setdefault("overrides", {})
    ov["ids"] = list(ov.get("ids") or []) + load_ids_file(
        str(_Path(path).expanduser())
    )
    ov[enable_key] = True


def _resolve_overrides(cfg: dict) -> None:
    """Resolve the overrides id list once and cache it on the cfg.

    Reads ``overrides.ids_file`` from disk (friendly error if missing),
    unions it with the inline ``overrides.ids``, normalises + version-
    augments the set, and stashes ``cfg["_overrides_runtime"]`` so every QC
    stage reads a pre-resolved policy (one file read per run, not per
    stage). Warns when force-keep is requested but no ids are listed —
    that protects nothing, and a silent no-op would be surprising.
    """
    from .overrides import resolve_raw_ids

    ov = cfg.get("overrides") or {}
    protect_qc = bool(ov.get("protect_qc", False))
    force_select = bool(ov.get("force_select", False))
    path = ov.get("ids_file")
    if path:
        from pathlib import Path as _Path

        if not _Path(path).expanduser().is_file():
            click.echo(
                f"[config error] overrides.ids_file not found: {path}. "
                f"Provide a readable file (one accession/isolate-id per line) "
                f"or remove the key.",
                err=True,
            )
            sys.exit(1)
        ov = {**ov, "ids_file": str(_Path(path).expanduser())}
    ids = resolve_ids(ov)
    stages = resolve_stages(ov.get("protect_stages", "all"))
    cfg["_overrides_runtime"] = {
        "ids": ids,
        "stages": stages,
        "protect_qc": protect_qc,
        "force_select": force_select,
        "raw_ids": resolve_raw_ids(ov),
    }
    if (protect_qc or force_select) and not ids:
        which = " / ".join(
            w for w, on in (("protect_qc", protect_qc),
                            ("force_select", force_select)) if on
        )
        click.echo(
            f"[overrides] {which} set but no overrides.ids / "
            f"overrides.ids_file were provided — nothing will be "
            f"{'protected/selected' if protect_qc and force_select else 'applied'}.",
            err=True,
        )
    elif ids:
        if protect_qc:
            n_stages = "all" if stages == resolve_stages("all") else len(stages)
            click.echo(
                f"[overrides] force-keep active: {len(ids)} id-form(s) "
                f"protected across {n_stages} QC stage(s).",
                err=True,
            )
        if force_select:
            click.echo(
                f"[overrides] force-select active: {len(ids)} id-form(s) "
                f"guaranteed in the representatives.",
                err=True,
            )


def _apply_fast_overrides(cfg: dict) -> None:
    """Apply --fast overrides to the loaded config, in place.

    Forces a preliminary-run tree pipeline regardless of YAML: FastTree,
    MAFFT ``--retree 1`` only (no ``--auto``, no L-INS-i), no trimAl, no
    partitioned supermatrix, midpoint rooting (skips the taxonomy-guided
    / MAD chain — all trees in the run are midpoint-rooted). Touches both
    the whole-genome (2E) and per-protein (2F) paths so a run can iterate
    fast end-to-end. The summary renderer reads the post-mutation cfg, so
    {prefix}_summary.md will describe what actually ran.
    """
    phylo = cfg.setdefault("phylo", {})
    phylo["tool"] = "fasttree"
    phylo.setdefault("partition", {})["enabled"] = False
    phylo.setdefault("trimal", {})["enabled"] = False
    phylo.setdefault("rooting", {})["method"] = "midpoint"
    mafft = phylo.setdefault("mafft", {})
    mafft["extra_args"] = ["--retree", "1"]
    mafft["use_auto"] = False
    pp = phylo.setdefault("per_protein", {})
    pp.setdefault("trimal", {})["enabled"] = False
    pp_mafft = pp.setdefault("mafft", {})
    pp_mafft["extra_args"] = ["--retree", "1"]
    click.echo(
        "[--fast] forcing FastTree, MAFFT '--retree 1' (no --auto), "
        "no trimAl, no partitioned supermatrix, midpoint rooting — "
        "overriding phylo: settings from YAML.",
        err=True,
    )


def _check_output_dir(cfg: dict) -> None:
    """Abort before doing any work if the output directory already exists and
    is not empty, so a previous run's files are never silently overwritten or
    mixed with a new run's results."""
    out_dir = Path(cfg["output"]["dir"])
    if not out_dir.exists():
        return
    if not out_dir.is_dir():
        click.echo(
            f"[error] output path '{out_dir}' already exists and is not a "
            f"directory. Choose another location with --output-dir.",
            err=True,
        )
        sys.exit(1)
    if any(out_dir.iterdir()):
        click.echo(
            f"[error] output directory '{out_dir}' already exists and is not "
            f"empty.\n        Remove it, empty it, or pick a different one "
            f"with --output-dir so results from separate runs are not mixed.",
            err=True,
        )
        sys.exit(1)


_SOURCE_MAP = {
    "uniprot": SequenceSource.UNIPROT,
    "ncbi": SequenceSource.NCBI,
    "ncbi_virus": SequenceSource.NCBI_VIRUS,
}


def _preflight_input(path: str) -> None:
    """Fail fast with a friendly message on a bad input FASTA path.

    Catches the mistakes that would otherwise surface as a raw traceback
    from ``open()`` deep inside ``read_fasta`` — a missing file, a
    directory, a permission problem, an empty file, or a file that plainly
    isn't FASTA. Raises :class:`~repseq.errors.InputError` (rendered without
    a traceback at the CLI boundary).
    """
    p = Path(path)
    if not p.exists():
        raise InputError(
            f"Input file not found: {path}\n"
            f"       Check the path(s) you passed to -i/--input."
        )
    if p.is_dir():
        raise InputError(
            f"Input path is a directory, not a FASTA file: {path}\n"
            f"       Pass the FASTA file itself (repeat -i for multiple files)."
        )
    try:
        if p.stat().st_size == 0:
            raise InputError(
                f"Input file is empty: {path}\n"
                f"       It contains no sequences — check that the download or "
                f"export actually produced data."
            )
        with open(p) as fh:
            head = fh.read(4096)
    except PermissionError as e:
        raise InputError(
            f"Input file is not readable (permission denied): {path}\n"
            f"       Check the file's permissions."
        ) from e
    except UnicodeDecodeError as e:
        raise InputError(
            f"Input file is not plain-text FASTA: {path}\n"
            f"       If it is gzip-compressed (.gz), decompress it first "
            f"(repseq reads uncompressed FASTA)."
        ) from e
    except OSError as e:
        raise InputError(
            f"Could not read input file {path}: {e.strerror or e}."
        ) from e
    stripped = head.lstrip()
    if not stripped.startswith(">"):
        raise InputError(
            f"Input file does not look like FASTA (no '>' header found): {path}\n"
            f"       repseq reads FASTA. If this is a different format (or a "
            f"gzip-compressed .gz file), convert or decompress it first."
        )


def _load_sequences(input_paths: tuple[str, ...], source_override: str = "auto") -> list:
    override = _SOURCE_MAP.get(source_override)  # None when "auto"
    # Validate every path up front so we don't print "Reading ..." for the
    # first file and then crash on the second.
    for path in input_paths:
        _preflight_input(path)
    sequences = []
    for path in input_paths:
        click.echo(f"Reading {path} ...")
        sequences.extend(read_fasta(path, source_override=override))
    click.echo(f"Loaded {len(sequences)} sequences.")
    if override:
        click.echo(f"Source override: {source_override}")
    return sequences


def _resolve_metadata(sequences, cfg, no_resolve):
    """Resolve metadata; return the NCBI client (or None if --no-resolve)
    so downstream steps like protein QC can reuse it."""
    if no_resolve:
        return None
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
    return ncbi


def _run_qc(sequences, cfg):
    click.echo("Running QC ...")
    passed, qc_report = run_qc(sequences, cfg)
    # NOTE: the QC summary is printed once at the end (after the
    # extended QC and segmented-completeness steps mutate the same
    # counters), not here. Printing it now would show 0 for every
    # field a later step is responsible for filling.
    return passed, qc_report


def _run_protein_qc(sequences, cfg, qc_report, ncbi):
    """Optional NCBI-backed protein-count filter. No-op if disabled or
    if no NCBI client is available (e.g. --no-resolve)."""
    if ncbi is None:
        return sequences
    virus_cfg = get_virus_config(cfg)
    qc_cfg = cfg.get("qc", {})
    pa_enabled = qc_cfg.get("protein_annotation", {}).get("enabled", False)
    has_per_segment = bool((virus_cfg or {}).get("expected_proteins_per_segment"))
    if not pa_enabled and not has_per_segment:
        return sequences
    click.echo("Running protein-annotation QC ...")
    kept = run_protein_qc(sequences, ncbi, cfg, virus_cfg, qc_report)
    click.echo(f"  Removed (proteins): {qc_report.removed_proteins}")
    return kept


def _any_marker_has_hmms(cfg: dict) -> bool:
    """True iff any configured marker spec carries an ``hmms`` list.

    Used to skip the hmmscan step entirely when no marker would consume
    the hits (the scan is purely diagnostic without configured HMMs).
    Covers cluster_protein, segment_markers, AND extra_protein — the last
    can declare HMMs without any of the others (e.g. extracting a sparse
    accessory protein on a run that clusters by NT only).
    """
    for entry in (cfg.get("clustering", {}).get("cluster_protein", []) or []):
        if isinstance(entry, dict) and (entry.get("hmms") or []):
            return True
    for entry in (cfg.get("clustering", {}).get("extra_protein", []) or []):
        if isinstance(entry, dict) and (entry.get("hmms") or []):
            return True
    # Polyprotein specs (v0.33.0+): each peptide carries either `hmm:`
    # (single token, legacy form) or `hmms:` (list of alternative
    # architectures, v0.34.0 OR form). Either one means the HMM tier
    # needs to fire so the slicer has hits to work with.
    for entry in (cfg.get("clustering", {}).get("polyprotein", []) or []):
        if isinstance(entry, dict):
            for pep in (entry.get("peptides") or []):
                if not isinstance(pep, dict):
                    continue
                if (pep.get("hmm") or "").strip():
                    return True
                hmms = pep.get("hmms") or []
                if isinstance(hmms, list) and any(
                    isinstance(t, str) and t.strip() for t in hmms
                ):
                    return True
    virus_cfg = get_virus_config(cfg)
    if virus_cfg:
        for entries in (virus_cfg.get("cluster_protein") or {}).values():
            for entry in (entries or []):
                if isinstance(entry, dict) and (entry.get("hmms") or []):
                    return True
        for spec in (virus_cfg.get("segment_markers") or {}).values():
            if isinstance(spec, dict) and (spec.get("hmms") or []):
                return True
        for entries in (virus_cfg.get("extra_protein") or {}).values():
            for entry in (entries or []):
                if isinstance(entry, dict) and (entry.get("hmms") or []):
                    return True
        for entries in (virus_cfg.get("polyprotein") or {}).values():
            for entry in (entries or []):
                if isinstance(entry, dict):
                    for pep in (entry.get("peptides") or []):
                        if not isinstance(pep, dict):
                            continue
                        if (pep.get("hmm") or "").strip():
                            return True
                        hmms = pep.get("hmms") or []
                        if isinstance(hmms, list) and any(
                            isinstance(t, str) and t.strip() for t in hmms
                        ):
                            return True
    return False


def _collect_config_hmm_names(cfg: dict) -> set[str]:
    """Return the set of individual HMM profile names referenced in config tokens.

    Splits every token string (single ``"Name"`` or multidomain ``"A--B--C"``)
    into its component names so callers can cross-check against the DB.
    Malformed tokens are silently skipped (config validation already flagged them).
    """
    from .hmm.runner import parse_hmm_token

    names: set[str] = set()

    def _add_tokens(token_list):
        for token in token_list or []:
            if not isinstance(token, str):
                continue
            try:
                names.update(parse_hmm_token(token))
            except ValueError:
                pass

    for entry in (cfg.get("clustering", {}).get("cluster_protein", []) or []):
        if isinstance(entry, dict):
            _add_tokens(entry.get("hmms", []))
    for entry in (cfg.get("clustering", {}).get("extra_protein", []) or []):
        if isinstance(entry, dict):
            _add_tokens(entry.get("hmms", []))
    # Polyprotein peptide HMMs (v0.34.0+: each peptide's locator is
    # either `hmm:` (single token) or `hmms:` (list of alternative
    # tokens). Parse each token via parse_hmm_token so we collect every
    # component HMM name across all alternatives.
    def _collect_peptide_tokens(pep: dict) -> list[str]:
        tokens: list[str] = []
        single = (pep.get("hmm") or "").strip()
        if single:
            tokens.append(single)
        for t in (pep.get("hmms") or []):
            if isinstance(t, str) and t.strip():
                tokens.append(t.strip())
        return tokens

    for entry in (cfg.get("clustering", {}).get("polyprotein", []) or []):
        if isinstance(entry, dict):
            for pep in (entry.get("peptides") or []):
                if isinstance(pep, dict):
                    for token in _collect_peptide_tokens(pep):
                        try:
                            names.update(parse_hmm_token(token))
                        except ValueError:
                            pass
    virus_cfg = get_virus_config(cfg)
    if virus_cfg:
        for entries in (virus_cfg.get("cluster_protein") or {}).values():
            for entry in (entries or []):
                if isinstance(entry, dict):
                    _add_tokens(entry.get("hmms", []))
        for spec in (virus_cfg.get("segment_markers") or {}).values():
            if isinstance(spec, dict):
                _add_tokens(spec.get("hmms", []))
        for entries in (virus_cfg.get("extra_protein") or {}).values():
            for entry in (entries or []):
                if isinstance(entry, dict):
                    _add_tokens(entry.get("hmms", []))
        for entries in (virus_cfg.get("polyprotein") or {}).values():
            for entry in (entries or []):
                if isinstance(entry, dict):
                    for pep in (entry.get("peptides") or []):
                        if isinstance(pep, dict):
                            for token in _collect_peptide_tokens(pep):
                                try:
                                    names.update(parse_hmm_token(token))
                                except ValueError:
                                    pass
    return names


def _run_hmm_scan(sequences, cfg, ncbi) -> None:
    """Run hmmscan over every CDS protein and stash results on cfg.

    Soft-fails: any error (binary missing, DB unreadable, hmmscan
    nonzero exit) emits one stderr line and leaves
    ``cfg["_hmm_runtime"]["active"]`` False so downstream marker
    selection falls back to alias/longest.

    Annotates ``seq.proteins[*]["hmm_hits"]`` with the raw hit list per
    CDS; the marker selector applies the configured E-value / coverage
    cutoffs at selection time.
    """
    cfg.setdefault("_hmm_runtime", {})["active"] = False
    hmm_cfg = cfg.get("hmm", {}) or {}
    if not hmm_cfg.get("enabled", True):
        return
    if not _any_marker_has_hmms(cfg):
        return

    from . import hmm as hmm_pkg

    if not hmm_pkg.is_available():
        click.echo(
            "[hmm] hmmscan not on PATH — HMM tier skipped; marker "
            "selection will fall back to alias/longest.",
            err=True,
        )
        return

    proteins_to_scan: dict[str, str] = {}
    sid_to_loc: dict[str, tuple[int, int]] = {}
    for si, seq in enumerate(sequences):
        if not seq.proteins:
            continue
        for pi, p in enumerate(seq.proteins):
            aa = p.get("sequence")
            if not aa:
                continue
            sid = f"S{si:06d}P{pi:03d}"
            proteins_to_scan[sid] = aa
            sid_to_loc[sid] = (si, pi)

    if not proteins_to_scan:
        return

    import time
    cache = ncbi.cache if ncbi is not None else None
    t0 = time.time()
    try:
        results = hmm_pkg.scan_proteins(proteins_to_scan, cfg, cache=cache)
        ga_cutoffs = hmm_pkg.get_ga_cutoffs(cfg)
    except hmm_pkg.HMMDatabaseError as e:
        click.echo(
            f"[hmm] database error: {e}; HMM tier skipped (falling back "
            "to alias/longest marker selection).",
            err=True,
        )
        return
    except hmm_pkg.HMMScanError as e:
        click.echo(
            f"[hmm] hmmscan failed: {e}; HMM tier skipped (falling back "
            "to alias/longest marker selection).",
            err=True,
        )
        return

    unknown_hmms = sorted(_collect_config_hmm_names(cfg) - set(ga_cutoffs.keys()))
    if unknown_hmms:
        click.echo(
            f"[hmm] WARNING: {len(unknown_hmms)} HMM name(s) referenced in config "
            "are not present in the database — they will never match any CDS "
            "(check for typos):\n"
            + "\n".join(f"  {n}" for n in unknown_hmms),
            err=True,
        )

    # Annotate each hit with a `passing: bool` once here so downstream
    # consumers (marker selector, isolate_proteins TSV/FASTA writers,
    # summary renderer) don't each re-apply the cutoffs in slightly
    # different ways.
    from .hmm.runner import passes_cutoffs as _pc
    use_ga = hmm_cfg.get("use_ga_when_available", True)
    default_ev = hmm_cfg.get("default_evalue", 1.0e-5)
    rel_len = hmm_cfg.get("relative_length_cutoff", 0.5)
    n_with_hits = 0
    n_with_passing = 0
    for sid, (si, pi) in sid_to_loc.items():
        hits = results.get(sid, [])
        for hit in hits:
            hit["passing"] = _pc(hit, ga_cutoffs, default_ev, rel_len, use_ga)
        sequences[si].proteins[pi]["hmm_hits"] = hits
        if hits:
            n_with_hits += 1
        if any(h.get("passing") for h in hits):
            n_with_passing += 1

    click.echo(
        f"[hmm] scanned {len(proteins_to_scan)} CDS across {len(sequences)} "
        f"sequences against {len(ga_cutoffs)} profiles "
        f"({n_with_passing} CDS with ≥1 passing hit, "
        f"{n_with_hits} with ≥1 raw hit) in {time.time() - t0:.1f}s"
    )
    cfg["_hmm_runtime"] = {
        "active": True,
        "ga_cutoffs": ga_cutoffs,
        "hmm_cfg": hmm_cfg,
    }


def _run_hmm_qc(sequences, cfg, qc_report, ncbi):
    """HMM-based QC step: verify each segment/sequence carries the expected
    marker proteins by HMM identity, regardless of clustering alphabet.

    This is the v0.14.0 promotion of the HMM tier from a marker-selection
    helper (v0.13.0, only ran when ``alphabet=protein``) to a first-class
    QC stage that fires whenever ``hmm.enabled=true`` AND any configured
    marker spec carries ``hmms:[…]``. The purpose is identity verification
    (protein X is actually protein X, structurally) instead of trusting
    GenBank ``/product`` strings, which are inconsistent across viral
    families and submitters.

    Pipeline placement: after ``_populate_genbank_isolate_segment`` and
    ``_filter_taxonomy_consistent``, before ``_setup_protein_alphabet``.
    Auto-fetches proteins from GenBank if not already attached.

    Semantic: each ``hmms:`` list element is a TOKEN string ("Name" or
    "A--B--C" multidomain in N-to-C order). A CDS satisfies a token when
    every named HMM has a passing hit AND those hits appear in N-to-C
    order along the protein. Tokens in one spec are ALTERNATIVE
    architectures (OR): a segment passes when at least one CDS satisfies
    ANY one of the spec's tokens, and fails only when none of them is.

    Drops:
        - Segmented: an isolate is dropped when ANY of its segments
          fails its spec. Counter: ``removed_hmm_failed`` (one bump per
          dropped isolate); ``removed_hmm_by_marker`` breaks down by
          "{segment}:{unmatched-token(s)}" (the alternatives joined with
          "|" when the spec declared more than one).
        - Non-segmented: a sequence is dropped when its spec has no
          satisfying CDS. Same counter, key = spec name.

    Soft-fails: missing ``hmmscan`` binary, unreadable DB, or
    ``hmm.enabled=false`` all bypass the step with a single stderr line.
    Marker-selection downstream falls back to alias/longest.
    """
    if not _any_marker_has_hmms(cfg):
        return sequences
    if not (cfg.get("hmm", {}) or {}).get("enabled", True):
        return sequences

    # Auto-fetch CDS proteins if QC didn't already do it.
    needs_proteins = any(seq.proteins is None for seq in sequences)
    if needs_proteins:
        if ncbi is None:
            click.echo(
                "[hmm] HMM QC requires GenBank CDS fetch, but --no-resolve "
                "was set. HMM gate skipped — set hmm.enabled=false to "
                "silence this warning, or drop --no-resolve.",
                err=True,
            )
            return sequences
        click.echo("Fetching CDS proteins for HMM-based QC ...")
        attach_proteins(sequences, ncbi)

    # Run scan + stash runtime context on cfg.
    _run_hmm_scan(sequences, cfg, ncbi)

    hmm_rt = cfg.get("_hmm_runtime", {}) or {}
    if not hmm_rt.get("active"):
        # Scan soft-failed (already printed a warning in _run_hmm_scan).
        return sequences

    segmented = bool(cfg.get("segmented", {}).get("enabled"))
    if segmented:
        return _run_hmm_qc_segmented(sequences, cfg, qc_report)
    return _run_hmm_qc_non_segmented(sequences, cfg, qc_report)


def _resolve_segment_hmms(
    seg_name: str,
    segment_markers: dict,
    cluster_protein_per_seg: dict,
) -> tuple[Optional[str], list[str]]:
    """Return ``(spec_name, token_list)`` for the segment's HMM gate.

    Lookup order: ``segment_markers[seg]`` first (the v0.13+ HMM-aware
    form). If it declares ``hmms:``, those tokens win. If it exists
    but its ``hmms:`` is empty (alias-only segment_markers entry),
    fall through to per-segment ``cluster_protein[seg]`` (legacy: any
    dict-form entry's ``hmms`` pulled into a flat token list) so a user
    who split aliases into ``segment_markers`` and HMMs into
    ``cluster_protein`` still gets the HMM gate fired. Returns
    ``(name, [])`` when no HMM gate is configured for the segment —
    that segment isn't QC'd by HMMs (alias-only / longest fallback in
    marker selection).
    """
    if seg_name in segment_markers:
        spec = segment_markers[seg_name] or {}
        tokens = list(spec.get("hmms") or [])
        if tokens:
            return seg_name, tokens
        # Fall through: segment_markers is alias-only for this segment —
        # consult cluster_protein for HMM tokens.
    if seg_name in cluster_protein_per_seg:
        entries = cluster_protein_per_seg[seg_name] or []
        tokens = []
        spec_name: Optional[str] = None
        for entry in entries:
            if isinstance(entry, dict):
                spec_name = spec_name or entry.get("name") or seg_name
                tokens.extend(entry.get("hmms") or [])
        return spec_name, tokens
    return None, []


def _segment_fails_hmm_gate(seq, tokens: list[str], overlap_tolerance: int = 0) -> Optional[str]:
    """Gate a segment against its ``hmms:`` token list. Return ``None`` when
    the segment passes, or a token string naming what was expected when it
    fails.

    The tokens in one ``hmms:`` list are **alternative architectures (OR)**:
    the segment passes as soon as *any* token is satisfied by some CDS, and
    fails only when *none* of them is. On failure the returned string names
    the unsatisfied alternative(s) — a single token, or the alternatives
    joined with ``|`` when more than one was declared — for the
    ``_qc_removed.tsv`` reason and the per-marker counter key.

    ``overlap_tolerance`` is forwarded to :func:`cds_satisfies_token` so
    Pfam-boundary fuzz at multidomain seams doesn't drop biologically
    valid isolates from the QC pool.
    """
    from .hmm.runner import cds_satisfies_token, parse_hmm_token

    proteins = seq.proteins or []
    declared: list[str] = []
    for token in tokens:
        try:
            parsed = parse_hmm_token(token)
        except ValueError:
            continue  # a malformed token can't be satisfied; skip as a candidate
        declared.append(token)
        if any(
            cds_satisfies_token(
                p.get("hmm_hits") or [], parsed,
                overlap_tolerance=overlap_tolerance,
            ) is not None
            for p in proteins
        ):
            return None  # at least one alternative architecture present → pass
    if not declared:
        # Nothing parseable to satisfy → fail, naming whatever was declared.
        return tokens[0] if tokens else "?"
    return "|".join(declared) if len(declared) > 1 else declared[0]


def _run_hmm_qc_segmented(sequences, cfg, qc_report):
    """Per-isolate HMM QC for segmented runs.

    Group by ``isolate_id``; for each isolate, check each segment's spec.
    If any segment fails its spec, drop the whole isolate (every segment
    seq is recorded in _qc_removed.tsv — the one that actually failed
    gets the structured reason, siblings get a sibling-dropped reason
    referencing the primary failure).

    Sequences without an ``isolate_id`` (UniProt input, missing GenBank
    qualifier, --no-resolve survivors) are skipped here — the
    completeness step's regex fallback may still group them, and
    ``build_concatenated_sequences`` applies the gate again as a
    backstop for that path.
    """
    virus_cfg = get_virus_config(cfg) or {}
    segment_markers = virus_cfg.get("segment_markers") or {}
    cluster_protein_per_seg = virus_cfg.get("cluster_protein") or {}
    overlap_tol = int(
        cfg.get("hmm", {}).get("multidomain_overlap_tolerance", 30)
    )

    by_isolate: dict[str, list] = {}
    for seq in sequences:
        if seq.isolate_id:
            by_isolate.setdefault(seq.isolate_id, []).append(seq)

    # Map isolate_id -> list of (segment, failed_token) tuples in input order.
    failures_by_isolate: dict[str, list[tuple[str, str]]] = {}
    for isolate_id, segs in by_isolate.items():
        for seq in segs:
            seg_name = seq.segment
            if not seg_name:
                continue
            _, tokens = _resolve_segment_hmms(
                seg_name, segment_markers, cluster_protein_per_seg
            )
            if not tokens:
                continue  # no HMM gate for this segment
            failed_token = _segment_fails_hmm_gate(seq, tokens, overlap_tol)
            if failed_token is not None:
                failures_by_isolate.setdefault(isolate_id, []).append(
                    (seg_name, failed_token)
                )

    # Force-keep whitelist: drop protected isolates from the failing set
    # (naming any one segment protects the whole isolate). Record each
    # protected isolate's segments with the would-be HMM-failure reason.
    policy = ProtectionPolicy.from_cfg(cfg)
    if policy.enabled:
        for isolate_id in list(failures_by_isolate):
            segs = by_isolate.get(isolate_id, [])
            if policy.protects_any(segs, "hmm"):
                primary = failures_by_isolate.pop(isolate_id)[0]
                reason = f"hmm_failed:{primary[0]}:{primary[1]}"
                for seq in segs:
                    qc_report.add_protected(seq.id, "hmm", reason)

    if not failures_by_isolate:
        return sequences

    failing_ids = set(failures_by_isolate.keys())
    kept: list = []
    n_dropped_segs = 0
    for seq in sequences:
        if seq.isolate_id and seq.isolate_id in failing_ids:
            failures = failures_by_isolate[seq.isolate_id]
            own = next(((s, t) for (s, t) in failures if s == seq.segment), None)
            if own is not None:
                reason = f"hmm_failed:{own[0]}:{own[1]}"
            else:
                primary = failures[0]
                reason = (
                    f"hmm_failed_sibling:{primary[0]}:{primary[1]}"
                )
            seq.qc_passed = False
            seq.qc_fail_reason = reason
            qc_report.add_removed(seq.id, reason)
            n_dropped_segs += 1
        else:
            kept.append(seq)

    # Counter bumps: one per dropped isolate; per-marker key is
    # "{segment}:{first_failing_token}" so the breakdown reflects which
    # combination is responsible (e.g. "L:RdRP_4--Mononeg_RNA_pol").
    for isolate_id, failures in failures_by_isolate.items():
        qc_report.removed_hmm_failed += 1
        first_seg, first_token = failures[0]
        key = f"{first_seg}:{first_token}"
        qc_report.removed_hmm_by_marker[key] = (
            qc_report.removed_hmm_by_marker.get(key, 0) + 1
        )

    top = ", ".join(
        f"{k}={v}"
        for k, v in sorted(
            qc_report.removed_hmm_by_marker.items(),
            key=lambda kv: -kv[1],
        )[:5]
    )
    click.echo(
        f"  HMM QC: dropped {len(failing_ids)} isolate(s) "
        f"({n_dropped_segs} segment-records). Top reasons: {top}"
    )
    return kept


def _run_hmm_qc_non_segmented(sequences, cfg, qc_report):
    """Per-sequence HMM QC for non-segmented runs.

    Each dict-form spec in ``clustering.cluster_protein`` that defines
    ``hmms:`` becomes a required marker. Within one spec the tokens are
    alternative architectures (OR): the sequence passes that spec when at
    least one CDS satisfies *any* of its tokens. Across specs the rule is
    AND — every HMM-defining spec is an independent required marker, so the
    sequence must satisfy them all. Alias-only specs are ignored for QC (no
    HMM gate to enforce). Specs with no ``hmms`` at all don't trigger this
    step.
    """
    cluster_protein = cfg.get("clustering", {}).get("cluster_protein", []) or []
    hmm_specs: list[tuple[str, list[str]]] = []
    for entry in cluster_protein:
        if not isinstance(entry, dict):
            continue
        toks = list(entry.get("hmms") or [])
        if toks:
            hmm_specs.append((entry.get("name") or ",".join(toks), toks))
    if not hmm_specs:
        return sequences

    overlap_tol = int(
        cfg.get("hmm", {}).get("multidomain_overlap_tolerance", 30)
    )
    policy = ProtectionPolicy.from_cfg(cfg)
    kept: list = []
    for seq in sequences:
        failed_spec: Optional[tuple[str, str]] = None
        for spec_name, tokens in hmm_specs:
            failed_token = _segment_fails_hmm_gate(seq, tokens, overlap_tol)
            if failed_token is not None:
                failed_spec = (spec_name, failed_token)
                break
        if failed_spec is None:
            kept.append(seq)
            continue
        spec_name, failed_token = failed_spec
        reason = f"hmm_failed:{spec_name}:{failed_token}"
        if protected_keep(seq, "hmm", reason, policy, qc_report):
            kept.append(seq)
            continue
        seq.qc_passed = False
        seq.qc_fail_reason = reason
        qc_report.add_removed(seq.id, reason)
        qc_report.removed_hmm_failed += 1
        qc_report.removed_hmm_by_marker[spec_name] = (
            qc_report.removed_hmm_by_marker.get(spec_name, 0) + 1
        )

    n_dropped = len(sequences) - len(kept)
    if n_dropped:
        top = ", ".join(
            f"{k}={v}"
            for k, v in sorted(
                qc_report.removed_hmm_by_marker.items(),
                key=lambda kv: -kv[1],
            )[:5]
        )
        click.echo(f"  HMM QC: dropped {n_dropped} sequence(s). Top reasons: {top}")
    return kept


_PROTEIN_BAD_CHARS = frozenset("XBZJ")


def _protein_bad_fraction(aa: Optional[str]) -> float:
    """Fraction of ambiguous residues (X/B/Z/J) in a protein translation.

    Uses the same protein-ambiguity set as ``Sequence.ambiguous_fraction``
    (U/O are definite residues, not ambiguity codes). An empty or missing
    translation returns 1.0 — a CDS feature with no usable ``/translation``
    is not a real protein and should be treated as fully bad.
    """
    if not aa:
        return 1.0
    s = aa.upper()
    return sum(1 for c in s if c in _PROTEIN_BAD_CHARS) / len(s)


def _segment_worst_bad_protein(seq, threshold: float) -> Optional[float]:
    """Return the worst (highest) over-threshold bad-fraction among the
    segment's CDS proteins, or ``None`` if every protein is clean.

    Sequences whose ``proteins`` is ``None`` (never fetched — UniProt
    input, no accession, --no-resolve) are not assessable and return
    ``None`` (treated as clean here; they fall through to later steps).
    """
    if seq.proteins is None:
        return None
    worst: Optional[float] = None
    for prot in seq.proteins:
        frac = _protein_bad_fraction(prot.get("sequence"))
        if frac > threshold:
            worst = frac if worst is None else max(worst, frac)
    return worst


def _run_protein_quality_qc(sequences, cfg, qc_report, ncbi):
    """Protein-quality QC: drop proteins whose translation is too noisy.

    The amino-acid analogue of the nucleotide ambiguous-character filter.
    The presence-only ``protein_annotation`` count check (and segmented
    completeness) verify that the *expected number* of proteins exists,
    but never inspect the residues — a segment carrying a translation
    filled with ambiguous residues (X/B/Z/J) would pass. This step closes
    that gap: a CDS protein whose ambiguous-residue fraction exceeds
    ``qc.protein_quality.max_bad_fraction`` is considered missing, which
    fails its segment and drops the whole isolate (segmented) or the
    sequence (non-segmented).

    Pipeline placement: after ``_filter_taxonomy_consistent``, immediately
    before ``_run_hmm_qc`` (so a noisy translation is gone before the HMM
    scan, which may not be enabled). When enabled it force-fetches GenBank
    CDS translations if no earlier step did; soft-skips with a stderr line
    under ``--no-resolve``.
    """
    pq_cfg = (cfg.get("qc", {}) or {}).get("protein_quality", {}) or {}
    if not pq_cfg.get("enabled", False):
        return sequences
    threshold = pq_cfg.get("max_bad_fraction", 0.05)

    needs_proteins = any(
        seq.proteins is None
        and seq.accession
        and seq.source != SequenceSource.UNIPROT
        for seq in sequences
    )
    if needs_proteins:
        if ncbi is None:
            click.echo(
                "[protein-quality] requires GenBank CDS translations, but "
                "--no-resolve was set. Protein-quality QC skipped — set "
                "qc.protein_quality.enabled=false to silence this warning, "
                "or drop --no-resolve.",
                err=True,
            )
            return sequences
        click.echo("Fetching CDS proteins for protein-quality QC ...")
        attach_proteins(sequences, ncbi)

    policy = ProtectionPolicy.from_cfg(cfg)
    segmented = bool(cfg.get("segmented", {}).get("enabled"))
    if segmented:
        return _run_protein_quality_qc_segmented(
            sequences, qc_report, threshold, policy
        )
    return _run_protein_quality_qc_non_segmented(
        sequences, qc_report, threshold, policy
    )


def _run_protein_quality_qc_segmented(sequences, qc_report, threshold, policy=None):
    """Per-isolate protein-quality QC for segmented runs.

    Group by ``isolate_id``; if any segment carries an over-threshold
    protein, drop the whole isolate. The failing segment(s) get a
    structured ``protein_quality:<seg>:bad_fraction=<f>>thr`` reason; the
    isolate's other segments get a ``protein_quality_sibling:<seg>:...``
    reason referencing the primary failure. Sequences without an
    ``isolate_id`` are skipped (the completeness regex fallback may still
    group them).
    """
    by_isolate: dict[str, list] = {}
    for seq in sequences:
        if seq.isolate_id:
            by_isolate.setdefault(seq.isolate_id, []).append(seq)

    failures_by_isolate: dict[str, list[tuple[str, float]]] = {}
    for isolate_id, segs in by_isolate.items():
        for seq in segs:
            if not seq.segment:
                continue
            worst = _segment_worst_bad_protein(seq, threshold)
            if worst is not None:
                failures_by_isolate.setdefault(isolate_id, []).append(
                    (seq.segment, worst)
                )

    # Force-keep whitelist: drop protected isolates from the failing set
    # (naming any one segment protects the whole isolate). Record each
    # protected isolate's segments before they would have been dropped.
    if policy is not None:
        for isolate_id in list(failures_by_isolate):
            segs = by_isolate.get(isolate_id, [])
            if policy.protects_any(segs, "protein_quality"):
                primary = failures_by_isolate.pop(isolate_id)[0]
                reason = (
                    f"protein_quality:{primary[0]}:"
                    f"bad_fraction={primary[1]:.3f}>{threshold}"
                )
                for seq in segs:
                    qc_report.add_protected(seq.id, "protein_quality", reason)

    if not failures_by_isolate:
        return sequences

    failing_ids = set(failures_by_isolate)
    kept: list = []
    n_dropped_segs = 0
    for seq in sequences:
        if seq.isolate_id and seq.isolate_id in failing_ids:
            failures = failures_by_isolate[seq.isolate_id]
            own = next(((s, f) for (s, f) in failures if s == seq.segment), None)
            if own is not None:
                reason = (
                    f"protein_quality:{own[0]}:bad_fraction={own[1]:.3f}>{threshold}"
                )
            else:
                primary = failures[0]
                reason = (
                    f"protein_quality_sibling:{primary[0]}:"
                    f"bad_fraction={primary[1]:.3f}>{threshold}"
                )
            seq.qc_passed = False
            seq.qc_fail_reason = reason
            qc_report.add_removed(seq.id, reason)
            n_dropped_segs += 1
        else:
            kept.append(seq)

    qc_report.removed_protein_quality += len(failing_ids)
    click.echo(
        f"  Protein-quality QC: dropped {len(failing_ids)} isolate(s) "
        f"({n_dropped_segs} segment-records) carrying a protein with "
        f"> {threshold:.0%} ambiguous residues."
    )
    return kept


def _run_protein_quality_qc_non_segmented(sequences, qc_report, threshold, policy=None):
    """Per-sequence protein-quality QC for non-segmented runs.

    A sequence is dropped when any of its CDS proteins exceeds the
    ambiguous-residue threshold.
    """
    kept: list = []
    for seq in sequences:
        worst = _segment_worst_bad_protein(seq, threshold)
        if worst is None:
            kept.append(seq)
            continue
        reason = f"protein_quality:bad_fraction={worst:.3f}>{threshold}"
        if protected_keep(seq, "protein_quality", reason, policy, qc_report):
            kept.append(seq)
            continue
        seq.qc_passed = False
        seq.qc_fail_reason = reason
        qc_report.add_removed(seq.id, reason)
        qc_report.removed_protein_quality += 1

    n_dropped = len(sequences) - len(kept)
    if n_dropped:
        click.echo(
            f"  Protein-quality QC: dropped {n_dropped} sequence(s) "
            f"carrying a protein with > {threshold:.0%} ambiguous residues."
        )
    return kept


def _setup_protein_alphabet(sequences, cfg, qc_report, ncbi):
    """When clustering.alphabet_for_clustering=protein, set
    seq.protein_sequence on each non-segmented sequence to its marker CDS.

    In v0.14.0 the HMM scan + QC drops have already happened upstream in
    ``_run_hmm_qc``, so this step is now pure marker selection: it picks
    the longest CDS that satisfies any HMM token in the spec (HMM tier
    active), or the legacy alias/longest CDS otherwise. No-op when
    ``alphabet_for_clustering=nucleotide`` — the HMM-QC tier already
    fired upstream regardless of alphabet.

    Triggers a CDS fetch if QC didn't already do it (only needed when
    no HMM specs configured, so ``_run_hmm_qc`` didn't pre-fetch).
    """
    alphabet = cfg.get("clustering", {}).get("alphabet_for_clustering", "protein")
    if alphabet == "nucleotide":
        return sequences
    needs_proteins = any(seq.proteins is None for seq in sequences)
    if needs_proteins:
        if ncbi is None:
            click.echo(
                "[error] clustering.alphabet_for_clustering='protein' "
                "requires GenBank CDS fetch, but --no-resolve was set. "
                "Either drop --no-resolve or switch to "
                "clustering.alphabet_for_clustering: nucleotide.",
                err=True,
            )
            sys.exit(1)
        click.echo("Fetching CDS proteins for protein-alphabet clustering ...")
        attach_proteins(sequences, ncbi)
    # Non-segmented: pick the marker per-sequence now. Segmented isolates
    # get a per-isolate concat in _handle_segmented.
    if not cfg.get("segmented", {}).get("enabled"):
        cluster_protein = cfg.get("clustering", {}).get("cluster_protein", []) or []
        concatenate = bool(cfg.get("clustering", {}).get("concatenate_markers", False))
        hmm_rt = cfg.get("_hmm_runtime", {}) or {}
        before = len(sequences)
        hmm_before = qc_report.removed_hmm_failed
        if concatenate:
            names = [
                (s.get("name") if isinstance(s, dict) else s)
                for s in cluster_protein
            ]
            click.echo(
                "Clustering on concatenated markers: "
                + "+".join(str(n) for n in names if n)
            )
        sequences = populate_protein_sequences(
            sequences,
            cluster_protein,
            qc_report,
            hmm_active=hmm_rt.get("active", False),
            ga_cutoffs=hmm_rt.get("ga_cutoffs"),
            hmm_cfg=hmm_rt.get("hmm_cfg"),
            concatenate=concatenate,
        )
        dropped = before - len(sequences)
        hmm_dropped = qc_report.removed_hmm_failed - hmm_before
        if dropped:
            base = f"  Dropped {dropped} sequence(s) with no marker protein"
            if hmm_dropped:
                breakdown = ", ".join(
                    f"{k}={v}"
                    for k, v in sorted(qc_report.removed_hmm_by_marker.items())
                )
                base += f" ({hmm_dropped} via HMM gate: {breakdown})"
            click.echo(base + ".")
    return sequences


def _populate_genbank_isolate_segment(sequences, cfg, ncbi):
    """Populate seq.isolate_id / seq.segment / seq.strain from the GenBank
    source feature, when segmented mode is on and the toggle is enabled.

    The downstream segmented filter prefers these fields over its regex fallback
    (extract_isolate_id and identify_segment both short-circuit when the
    Sequence already carries the value). UniProt input and sequences without
    an accession are skipped — the regex remains the fallback for those.
    A failed/empty NCBI fetch is also silent: the regex still runs.
    """
    if ncbi is None:
        return  # --no-resolve
    seg_cfg = cfg.get("segmented", {})
    if not seg_cfg.get("enabled"):
        return
    if not seg_cfg.get("use_genbank_metadata", True):
        return
    accessions: list[str] = []
    by_acc: dict[str, list] = {}
    for seq in sequences:
        if seq.source == SequenceSource.UNIPROT:
            continue
        if not seq.accession:
            continue
        accessions.append(seq.accession)
        by_acc.setdefault(seq.accession, []).append(seq)
    if not accessions:
        return
    click.echo("Fetching GenBank source metadata for segmented mode ...")
    meta_by_acc = ncbi.fetch_source_metadata_batch(
        accessions, progress=_stderr_batch_progress,
    )
    populated = 0
    for acc, meta in meta_by_acc.items():
        # Prefer /isolate; fall back to /strain when /isolate is absent.
        # Track which qualifier supplied the value on
        # ``seq.isolate_id_source`` so downstream tooling (TSV writer,
        # collision detector, run summary) can flag the over-merge risk
        # of strain-derived ids — a single named strain is often shared
        # across distinct biological samples.
        if meta.get("isolate"):
            isolate, isolate_source = meta["isolate"], "isolate"
        elif meta.get("strain"):
            isolate, isolate_source = meta["strain"], "strain"
        else:
            isolate, isolate_source = None, None
        segment = meta.get("segment")
        strain = meta.get("strain")
        serotype = meta.get("serotype")
        for seq in by_acc.get(acc, []):
            if isolate and not seq.isolate_id:
                seq.isolate_id = isolate
                seq.isolate_id_source = isolate_source
            if segment and not seq.segment:
                seq.segment = segment
            if strain and not seq.strain:
                seq.strain = strain
            if serotype and not seq.subtype:
                seq.subtype = serotype
            if isolate or segment:
                populated += 1
    click.echo(
        f"  Populated isolate/segment fields on {populated} of "
        f"{len(sequences)} sequences from GenBank metadata."
    )
    # Provenance breakdown — surfaces the strain-as-isolate fallback
    # before clustering happens, so a bench scientist sees it in the
    # run log without having to grep the TSV.
    src_counts = {"isolate": 0, "strain": 0}
    for seq in sequences:
        if seq.isolate_id_source in src_counts:
            src_counts[seq.isolate_id_source] += 1
    if src_counts["isolate"] or src_counts["strain"]:
        click.echo(
            f"  Isolate IDs from: {src_counts['isolate']} /isolate, "
            f"{src_counts['strain']} /strain "
            f"(strain-derived ids may over-merge — see the strain-collision "
            f"check below if any)."
        )


def _filter_taxonomy_consistent(sequences, cfg, qc_report):
    """Drop isolates whose segments disagree on the configured taxonomy
    rank (default species). No-op when segmented mode is off, when the
    check is disabled, or when no segment carries a populated rank
    value. Mutates qc_report; returns the surviving sequence list.
    """
    seg_cfg = cfg.get("segmented", {}) or {}
    if not seg_cfg.get("enabled"):
        return sequences
    tc_cfg = seg_cfg.get("taxonomy_consistency", {}) or {}
    if not tc_cfg.get("enabled", True):
        return sequences
    rank = tc_cfg.get("rank", "species")
    before = len(sequences)
    policy = ProtectionPolicy.from_cfg(cfg)
    protected_out: list[tuple[str, str, str]] = []
    kept, removed = filter_taxonomy_consistent_isolates(
        sequences, rank=rank, policy=policy, protected_out=protected_out,
    )
    for seq_id, stage, would_be_reason in protected_out:
        qc_report.add_protected(seq_id, stage, would_be_reason)
    if not removed:
        return kept
    for accession, reason in removed:
        qc_report.add_removed(accession, reason)
        qc_report.removed_taxonomy_mismatch += 1
    click.echo(
        f"  Taxonomy QC: dropped {before - len(kept)} segment(s) "
        f"for {rank}-level mismatch within isolate."
    )
    return kept


def _check_strain_collisions(sequences, cfg, qc_report):
    """Run the strain-collision detector and act on it per config.

    A collision is two or more distinct accessions sharing the same
    strain-derived isolate_id AND the same segment — the over-merge
    signature of the /strain -> isolate_id fallback. ``warn`` (default)
    prints one line per collision; ``drop`` removes every accession
    involved in any collision before the completeness filter runs.

    Returns the (possibly filtered) sequence list.
    """
    seg_cfg = cfg.get("segmented", {}) or {}
    action = seg_cfg.get("strain_collision_action", "warn")
    collisions = detect_strain_collisions(sequences)
    if not collisions:
        return sequences
    click.echo(
        f"  Strain-collision check: found {len(collisions)} (isolate, segment) "
        f"pair(s) where /strain is shared across distinct accessions.",
        err=True,
    )
    for (iso, seg), accs in sorted(collisions.items()):
        click.echo(
            f"    isolate '{iso}' segment '{seg}': "
            f"{len(accs)} accessions ({', '.join(accs)})",
            err=True,
        )
    if action != "drop":
        click.echo(
            "    Action: warn (no records dropped). Set "
            "segmented.strain_collision_action: drop to remove them.",
            err=True,
        )
        return sequences
    # Drop every sequence whose (isolate_id, segment) appears in
    # collisions. Use the same (iso, seg) key as the detector so the
    # set membership is exact.
    bad_keys = set(collisions.keys())
    kept: list = []
    for seq in sequences:
        key = (seq.isolate_id, seq.segment)
        if key in bad_keys:
            reason = f"strain_collision:{seq.segment}"
            seq.qc_passed = False
            seq.qc_fail_reason = reason
            qc_report.removed_strain_collisions += 1
            qc_report.add_removed(seq.accession or seq.id, reason)
        else:
            kept.append(seq)
    click.echo(
        f"    Action: drop ({len(sequences) - len(kept)} segment(s) removed).",
        err=True,
    )
    return kept


def _handle_segmented(sequences, cfg, qc_report):
    virus_cfg = get_virus_config(cfg)
    if not virus_cfg:
        # Non-segmented mode: the input list IS what feeds the mode,
        # so this is the final survivor count. Units = "sequences".
        qc_report.final_survivors = len(sequences)
        qc_report.final_survivors_unit = "sequences"
        return sequences, None, None
    sequences = _check_strain_collisions(sequences, cfg, qc_report)
    click.echo("Applying segmented virus completeness filter ...")
    seg_cfg = cfg.get("segmented", {}) or {}
    extra_action = seg_cfg.get("extra_segments_action", "warn")
    kept, complete_isolates, extras_by_isolate = filter_complete_isolates(
        sequences, virus_cfg, qc_report, extra_segments_action=extra_action
    )
    if extras_by_isolate:
        click.echo(
            f"  Extra-segments check: found {len(extras_by_isolate)} isolate(s) "
            f"with segments outside the configured set "
            f"({', '.join(virus_cfg['segments'])}).",
            err=True,
        )
        for iso, extras in sorted(extras_by_isolate.items()):
            click.echo(
                f"    isolate '{iso}': extras = {', '.join(extras)}", err=True
            )
        if extra_action == "drop":
            click.echo(
                f"    Action: drop ({len(extras_by_isolate)} isolate(s) "
                f"removed).",
                err=True,
            )
        else:
            click.echo(
                "    Action: warn (no isolates dropped; only the expected "
                "segments enter the concat). Set "
                "segmented.extra_segments_action: drop to remove them.",
                err=True,
            )
    segment_lengths = virus_cfg.get("segment_lengths")
    if segment_lengths:
        complete_isolates = segment_length_filter(
            complete_isolates, virus_cfg["segments"], segment_lengths, qc_report
        )
        kept = [seq for segs in complete_isolates.values() for seq in segs]
        per_seg = qc_report.removed_length_by_segment
        if per_seg:
            total = sum(c["too_short"] + c["too_long"] for c in per_seg.values())
            click.echo(f"  Segment-length filter: dropped {total} isolate(s)")
            for seg_name in virus_cfg["segments"]:
                counts = per_seg.get(seg_name)
                if not counts:
                    continue
                bounds = segment_lengths.get(seg_name) or {}
                mn = bounds.get("min")
                mx = bounds.get("max")
                if counts["too_short"]:
                    bound = f"<{mn}" if mn is not None else ""
                    click.echo(
                        f"    {seg_name} too short {bound}: {counts['too_short']}"
                    )
                if counts["too_long"]:
                    bound = f">{mx}" if mx is not None else ""
                    click.echo(
                        f"    {seg_name} too long  {bound}: {counts['too_long']}"
                    )
    click.echo(f"  Complete isolates : {len(complete_isolates)}")
    click.echo(f"  Individual seqs   : {len(kept)}")
    alphabet = cfg.get("clustering", {}).get("alphabet_for_clustering", "protein")
    require_protein = alphabet == "protein"
    cluster_protein = virus_cfg.get("cluster_protein")
    segment_markers = virus_cfg.get("segment_markers")
    hmm_rt = cfg.get("_hmm_runtime", {}) or {}
    hmm_before = qc_report.removed_hmm_failed
    hmm_before_by_marker = dict(qc_report.removed_hmm_by_marker)
    concat_seqs = build_concatenated_sequences(
        complete_isolates,
        segment_names=virus_cfg["segments"],
        cluster_protein=cluster_protein,
        require_protein=require_protein,
        report=qc_report,
        segment_markers=segment_markers,
        hmm_active=hmm_rt.get("active", False),
        ga_cutoffs=hmm_rt.get("ga_cutoffs"),
        hmm_cfg=hmm_rt.get("hmm_cfg"),
        policy=ProtectionPolicy.from_cfg(cfg),
    )
    # Drop isolates that lost a marker on any segment from complete_isolates
    # too, so the per-segment FASTA writer and isolate_proteins.tsv don't
    # resurrect them.
    if require_protein:
        survivors = {s.isolate_id for s in concat_seqs}
        before = len(complete_isolates)
        complete_isolates = {
            k: v for k, v in complete_isolates.items() if k in survivors
        }
        dropped = before - len(complete_isolates)
        hmm_dropped = qc_report.removed_hmm_failed - hmm_before
        if dropped:
            base = (
                f"  Dropped {dropped} isolate(s) with no marker protein "
                f"on one or more segments"
            )
            if hmm_dropped:
                # Per-marker delta since the call started
                delta = {
                    k: qc_report.removed_hmm_by_marker[k]
                    - hmm_before_by_marker.get(k, 0)
                    for k in qc_report.removed_hmm_by_marker
                }
                delta = {k: v for k, v in delta.items() if v}
                breakdown = ", ".join(
                    f"{k}={v}" for k, v in sorted(delta.items())
                )
                base += f" ({hmm_dropped} via HMM gate: {breakdown})"
            click.echo(base + ".")

    # Exact-duplicate removal was skipped on the segment pool (a conserved
    # segment shared between distinct isolates must not knock either isolate
    # out as incomplete). Apply it now, on the concatenated isolates: two
    # isolates are true duplicates only if every segment is identical.
    if qc_report.dedup_skipped:
        before = len(concat_seqs)
        concat_seqs = remove_duplicates(concat_seqs, qc_report)
        dropped = before - len(concat_seqs)
        if dropped:
            # Drop the de-duplicated isolates from complete_isolates too, so the
            # per-segment output files don't resurrect them. concatenate_isolate
            # carries the complete_isolates key through as Sequence.isolate_id.
            survivors = {s.isolate_id for s in concat_seqs}
            complete_isolates = {
                k: v for k, v in complete_isolates.items() if k in survivors
            }
        click.echo(f"  Duplicate isolates: {dropped} removed")
    # Segmented mode: the mode consumes one CONCAT per surviving
    # isolate, so the final survivor unit is "isolates", not segments.
    qc_report.final_survivors = len(concat_seqs)
    qc_report.final_survivors_unit = "isolates"
    return concat_seqs, complete_isolates, virus_cfg.get("segments")


def _drop_unwanted_newick(out_files: list, cfg: dict) -> int:
    """Delete every ``*_tree.nwk`` and prune it from ``out_files`` unless
    ``phylo.newick`` is true.

    The Newick is an unavoidable intermediate — the phyloXML is re-parsed
    from it and the per-protein incongruence RF table reads it — so it's
    generated during the run and dropped here, AFTER every phylo step and
    the incongruence/conservation sweeps have consumed it. The annotated
    phyloXML (``*_tree.xml``) is a topological superset; the short-id
    ``*_tree_id_map.tsv`` is kept regardless because the retained
    ``*_msa.fasta`` uses the same short-id leaves. Mutates ``out_files`` in
    place; returns the number of files dropped.
    """
    if (cfg.get("phylo", {}) or {}).get("newick", False):
        return 0
    kept: list = []
    dropped = 0
    for f in out_files:
        p = Path(f)
        if p.name.endswith("_tree.nwk"):
            try:
                p.unlink()
            except OSError:
                pass
            dropped += 1
        else:
            kept.append(f)
    if dropped:
        out_files[:] = kept
    return dropped


def _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names,
                  pre_clustering_sequences=None, plot: bool = False, phylo: bool = False,
                  per_protein_phylo: bool = False, per_segment_phylo: bool = False,
                  pre_cluster_tree: bool = False):
    # Force-select (overrides.force_select): guarantee pinned sequences are
    # representatives BEFORE anything is written or treed, so they land in
    # the rep FASTA/TSV, the phylogeny, and the plot. No-op unless enabled.
    from .overrides import apply_force_select

    apply_force_select(result, pre_clustering_sequences, cfg)
    if result.force_selected:
        from collections import Counter

        counts = Counter(e["action"] for e in result.force_selected)
        click.echo(
            "  Force-select: "
            + ", ".join(f"{k}={v}" for k, v in sorted(counts.items()))
            + " (see {prefix}_force_selected.tsv).".replace(
                "{prefix}", cfg.get("output", {}).get("prefix", "repseq")
            )
        )
        unavailable = [e["id"] for e in result.force_selected
                       if e["action"] == "unavailable"]
        if unavailable:
            click.echo(
                f"  [overrides] WARNING: {len(unavailable)} force-select id(s) "
                f"did not survive QC and could not be selected: "
                f"{', '.join(unavailable[:10])}"
                f"{' …' if len(unavailable) > 10 else ''}. Add them to "
                f"--protect-ids / overrides.protect_qc to keep them through QC.",
                err=True,
            )
    out_files = write_results(result, cfg, complete_isolates, segment_names)
    # Pre-cluster overview tree (2H): a rough single-pass FastTree
    # over every post-QC sequence with [repr] prefixes on the elected
    # representatives. Honoured either by --pre-cluster-tree or by
    # phylo.pre_cluster_tree.enabled in the YAML. Runs BEFORE the
    # post-cluster phylo block since logically it depicts the
    # pre-clustering view; soft-fails like the other phylo steps.
    pc_cfg = (cfg.get("phylo", {}) or {}).get("pre_cluster_tree", {}) or {}
    if pre_cluster_tree or pc_cfg.get("enabled", False):
        if pre_clustering_sequences is None or len(pre_clustering_sequences) < 3:
            click.echo(
                "[pre-cluster tree skipped] need >= 3 post-QC sequences",
                err=True,
            )
        else:
            out_dir = Path(cfg["output"]["dir"])
            prefix = cfg["output"].get("prefix", "repseq")
            try:
                from .phylo import PhyloError, run_pre_cluster_phylogeny
                pc_files = run_pre_cluster_phylogeny(
                    pre_clustering_sequences,
                    result.representatives,
                    cfg, out_dir, prefix,
                )
                out_files.extend(pc_files)
            except PhyloError as exc:
                click.echo(f"[pre-cluster tree skipped] {exc}", err=True)
            except Exception as exc:
                click.echo(f"[pre-cluster tree failed] {exc}", err=True)
    if plot:
        out_dir = Path(cfg["output"]["dir"])
        prefix = cfg["output"].get("prefix", "repseq")
        plot_path = out_dir / f"{prefix}_clustering.png"
        try:
            from .viz.clustering_plot import write_clustering_plot
            written = write_clustering_plot(
                result, plot_path, seed=cfg.get("seed", 42),
            )
            if written:
                out_files.append(written)
        except ImportError as exc:
            click.echo(f"[plot skipped] {exc}", err=True)
        except Exception as exc:
            click.echo(f"[plot failed] {exc}", err=True)
    if phylo:
        out_dir = Path(cfg["output"]["dir"])
        prefix = cfg["output"].get("prefix", "repseq")
        try:
            from .phylo import PhyloError, run_phylogeny
            phylo_files = run_phylogeny(
                result.representatives, cfg, out_dir, prefix,
            )
            out_files.extend(phylo_files)
        except PhyloError as exc:
            click.echo(f"[phylo skipped] {exc}", err=True)
        except Exception as exc:
            click.echo(f"[phylo failed] {exc}", err=True)
        # Phylogeny-based taxonomy review: the review TSV (imputation
        # ledger) was written inside the phylo step; here we materialise the
        # corrected rep TSV + protein FASTA with high-confidence imputed
        # blanks filled (clean values — the ledger records which cells).
        # Gated on write_corrected + actual imputations; soft-fails.
        review = cfg.get("_taxonomy_review")
        rcfg = (cfg.get("phylo", {}) or {}).get("taxonomy_review", {}) or {}
        if review and review.get("imputations") and rcfg.get("write_corrected", True):
            try:
                import copy as _copy
                from .output.report import (
                    write_proteins_fasta,
                    write_representative_isolates_tsv,
                    write_representative_sequences_tsv,
                )
                from .phylo.taxonomy_review import apply_imputations
                out_dir = Path(cfg["output"]["dir"])
                prefix = cfg["output"].get("prefix", "repseq")
                segmented = bool(cfg.get("segmented", {}).get("enabled"))
                corrected_reps, corrected_ci = apply_imputations(
                    result.representatives, complete_isolates,
                    review["imputations"],
                )
                if segmented:
                    tsv_path = out_dir / f"{prefix}_representative_isolates_corrected.tsv"
                    write_representative_isolates_tsv(corrected_reps, tsv_path)
                    fasta_path = out_dir / f"{prefix}_representative_isolate_proteins_corrected.fasta"
                else:
                    tsv_path = out_dir / f"{prefix}_representative_sequences_corrected.tsv"
                    write_representative_sequences_tsv(corrected_reps, tsv_path)
                    fasta_path = out_dir / f"{prefix}_representative_sequence_proteins_corrected.fasta"
                out_files.append(tsv_path)
                corrected_result = _copy.copy(result)
                corrected_result.representatives = corrected_reps
                if write_proteins_fasta(corrected_result, corrected_ci, fasta_path):
                    out_files.append(fasta_path)
                n_imp = sum(len(v) for v in review["imputations"].values())
                click.echo(
                    f"Wrote taxonomy-corrected copies ({n_imp} imputed "
                    f"value(s) filled): {tsv_path.name}",
                )
            except Exception as exc:
                click.echo(f"[taxonomy-corrected copy skipped] {exc}", err=True)
    if per_protein_phylo:
        out_dir = Path(cfg["output"]["dir"])
        prefix = cfg["output"].get("prefix", "repseq")
        try:
            from .phylo import PhyloError, run_per_protein_phylogeny
            pp_files = run_per_protein_phylogeny(
                result.representatives, cfg, out_dir, prefix,
            )
            out_files.extend(pp_files)
        except PhyloError as exc:
            click.echo(f"[per-protein phylo skipped] {exc}", err=True)
        except Exception as exc:
            click.echo(f"[per-protein phylo failed] {exc}", err=True)
    if per_segment_phylo:
        out_dir = Path(cfg["output"]["dir"])
        prefix = cfg["output"].get("prefix", "repseq")
        try:
            from .phylo import PhyloError, run_per_segment_phylogeny
            ps_files = run_per_segment_phylogeny(
                result.representatives, cfg, out_dir, prefix,
            )
            out_files.extend(ps_files)
        except PhyloError as exc:
            click.echo(f"[per-segment phylo skipped] {exc}", err=True)
        except Exception as exc:
            click.echo(f"[per-segment phylo failed] {exc}", err=True)
    # Per-MSA conservation scoring — a post-hoc sweep over every
    # alignment any phylo step wrote this run, collected into one
    # {prefix}_msa_conservation.tsv (mean per-column JSD-to-background,
    # Henikoff-weighted, gap-penalised). Decoupled and soft: finds no
    # MSAs and writes nothing when no tree was built, and a scoring bug
    # never voids the trees already on disk.
    try:
        from .phylo.conservation import write_msa_conservation_report
        out_dir = Path(cfg["output"]["dir"])
        prefix = cfg["output"].get("prefix", "repseq")
        cons_tsv = write_msa_conservation_report(out_dir, prefix, cfg)
        if cons_tsv is not None:
            out_files.append(cons_tsv)
    except Exception as exc:
        click.echo(f"[MSA conservation report skipped] {exc}", err=True)
    # Newick retention is opt-in (phylo.newick, default false). Drop the
    # *_tree.nwk files HERE — after every phylo step + the incongruence table
    # + the conservation sweep have consumed them — and prune them from the
    # tracked list so _summary.md / the report writers never name a file
    # that is no longer on disk.
    n_dropped = _drop_unwanted_newick(out_files, cfg)
    if n_dropped:
        click.echo(
            f"Dropped {n_dropped} Newick tree file(s) (phylo.newick: false). "
            f"The annotated phyloXML (_tree.xml) is retained; pass --newick "
            f"to keep the .nwk files too.",
        )
    write_all_reports(
        result, qc_report, cfg, list(input_paths), out_files,
        complete_isolates=complete_isolates,
        pre_clustering_sequences=pre_clustering_sequences,
    )
    # One-line provenance header shared by the four plain-text taxonomic
    # reports, so each is self-describing when read in isolation (dataset
    # type, clustering substrate, tool, rep count). Computed once; soft so a
    # build failure never blocks the reports.
    try:
        from .output.summary import build_provenance_header
        _report_provenance = build_provenance_header(
            cfg, result, segmented=bool(complete_isolates),
        )
    except Exception:
        _report_provenance = None
    # Taxonomic diversity report — distinct taxa per rank before/after
    # clustering, plus a per-taxon breakdown for low-diversity ranks. The
    # "before" pool is the post-QC sequence list fed to the mode (CONCAT
    # isolates in segmented mode). Soft-fail so a render bug never voids
    # a real selection.
    if pre_clustering_sequences is not None:
        try:
            out_dir = Path(cfg["output"]["dir"])
            prefix = cfg["output"].get("prefix", "repseq")
            tax_path = out_dir / f"{prefix}_taxonomic_report.txt"
            write_taxonomic_report(
                pre_clustering_sequences,
                result.representatives,
                segmented=bool(complete_isolates),
                path=tax_path,
                provenance=_report_provenance,
            )
            out_files.append(tax_path)
        except Exception as exc:
            click.echo(f"[taxonomic report skipped] {exc}", err=True)
        try:
            out_dir = Path(cfg["output"]["dir"])
            prefix = cfg["output"].get("prefix", "repseq")
            tax_tsv = out_dir / f"{prefix}_taxonomic_report.tsv"
            if write_taxonomic_report_tsv(
                pre_clustering_sequences,
                result.representatives,
                path=tax_tsv,
            ):
                out_files.append(tax_tsv)
        except Exception as exc:
            click.echo(f"[taxonomic report TSV skipped] {exc}", err=True)
    # Per-protein coverage report — for each declared marker /
    # extra_protein, the fraction of isolates / sequences carrying it
    # per taxonomic rank, plus length statistics. Only emitted when at
    # least one protein spec is declared; soft-fail so a render bug
    # never voids a real selection.
    if pre_clustering_sequences is not None:
        try:
            out_dir = Path(cfg["output"]["dir"])
            prefix = cfg["output"].get("prefix", "repseq")
            pr_path = out_dir / f"{prefix}_protein_taxonomic_report.txt"
            if write_protein_taxonomic_report(
                pre_clustering_sequences,
                result.representatives,
                cfg,
                segmented=bool(complete_isolates),
                path=pr_path,
                provenance=_report_provenance,
            ):
                out_files.append(pr_path)
        except Exception as exc:
            click.echo(f"[protein taxonomic report skipped] {exc}", err=True)
        try:
            out_dir = Path(cfg["output"]["dir"])
            prefix = cfg["output"].get("prefix", "repseq")
            pr_tsv = out_dir / f"{prefix}_protein_taxonomic_report.tsv"
            if write_protein_taxonomic_report_tsv(
                pre_clustering_sequences,
                result.representatives,
                cfg,
                path=pr_tsv,
            ):
                out_files.append(pr_tsv)
        except Exception as exc:
            click.echo(f"[protein taxonomic report TSV skipped] {exc}", err=True)
    # Per-rank NT length statistics: per-segment lengths + a `total`
    # column (segmented) or a single `genome` column (non-segmented).
    # Always-on; soft-fails so a render bug never voids the selection.
    if pre_clustering_sequences is not None:
        try:
            out_dir = Path(cfg["output"]["dir"])
            prefix = cfg["output"].get("prefix", "repseq")
            nt_path = out_dir / f"{prefix}_nucleotide_taxonomic_report.txt"
            if write_nucleotide_taxonomic_report(
                pre_clustering_sequences,
                result.representatives,
                cfg,
                segmented=bool(complete_isolates),
                path=nt_path,
                provenance=_report_provenance,
            ):
                out_files.append(nt_path)
        except Exception as exc:
            click.echo(f"[nucleotide taxonomic report skipped] {exc}", err=True)
        try:
            out_dir = Path(cfg["output"]["dir"])
            prefix = cfg["output"].get("prefix", "repseq")
            nt_tsv = out_dir / f"{prefix}_nucleotide_taxonomic_report.tsv"
            if write_nucleotide_taxonomic_report_tsv(
                pre_clustering_sequences,
                result.representatives,
                cfg,
                segmented=bool(complete_isolates),
                path=nt_tsv,
            ):
                out_files.append(nt_tsv)
        except Exception as exc:
            click.echo(f"[nucleotide taxonomic report TSV skipped] {exc}", err=True)
    # Per-peptide coverage + length report for declared polyprotein
    # specs — the sliced-peptide analogue of the protein taxonomic
    # report. Only emitted when at least one `polyprotein:` spec is
    # declared AND the HMM tier ran; soft-fail so a render bug never
    # voids a real selection.
    if pre_clustering_sequences is not None:
        try:
            out_dir = Path(cfg["output"]["dir"])
            prefix = cfg["output"].get("prefix", "repseq")
            pp_path = out_dir / f"{prefix}_polyprotein_taxonomic_report.txt"
            if write_polyprotein_taxonomic_report(
                pre_clustering_sequences,
                result.representatives,
                cfg,
                segmented=bool(complete_isolates),
                path=pp_path,
                provenance=_report_provenance,
            ):
                out_files.append(pp_path)
        except Exception as exc:
            click.echo(f"[polyprotein taxonomic report skipped] {exc}", err=True)
        try:
            out_dir = Path(cfg["output"]["dir"])
            prefix = cfg["output"].get("prefix", "repseq")
            pp_tsv = out_dir / f"{prefix}_polyprotein_taxonomic_report.tsv"
            if write_polyprotein_taxonomic_report_tsv(
                pre_clustering_sequences,
                result.representatives,
                cfg,
                path=pp_tsv,
            ):
                out_files.append(pp_tsv)
        except Exception as exc:
            click.echo(f"[polyprotein taxonomic report TSV skipped] {exc}", err=True)
    # Methods-section starter — written after every successful run so
    # a bench scientist can copy it into a paper. Soft-fail (one stderr
    # line) so a render bug never voids a real selection.
    try:
        from .output.summary import write_summary
        summary_path = write_summary(
            cfg, qc_report, result, list(input_paths),
            complete_isolates=complete_isolates,
            segment_names=segment_names,
            phylo_ran=phylo,
            per_protein_ran=per_protein_phylo,
            per_segment_ran=per_segment_phylo,
            pre_cluster_ran=bool(
                pre_cluster_tree
                or (cfg.get("phylo", {}) or {})
                    .get("pre_cluster_tree", {}).get("enabled", False)
            ),
            command=" ".join(sys.argv),
        )
        out_files.append(summary_path)
    except Exception as exc:
        click.echo(f"[summary skipped] {exc}", err=True)
    # Lockfile — machine-readable reproducibility record. Always
    # emitted (no flag); soft-fails so a serialisation bug never
    # voids a real selection.
    try:
        from .lockfile import build_lockfile, write_lockfile
        out_dir = Path(cfg["output"]["dir"])
        prefix = cfg["output"].get("prefix", "repseq")
        lf_path = out_dir / f"{prefix}_lockfile.json"
        lockfile = build_lockfile(
            cfg, result, list(input_paths), command=" ".join(sys.argv),
        )
        write_lockfile(lockfile, lf_path)
        out_files.append(lf_path)
    except Exception as exc:
        click.echo(f"[lockfile skipped] {exc}", err=True)
    # Sanitized, fully-resolved config snapshot — every setting at the
    # value it ran with (defaults filled in), comments stripped, NCBI
    # credentials blanked. Re-runnable as a config. Always emitted;
    # soft-fails so a dump bug never voids a real selection.
    try:
        from .config import effective_config_filename, write_effective_config
        out_dir = Path(cfg["output"]["dir"])
        prefix = cfg["output"].get("prefix", "repseq")
        cfg_path = out_dir / effective_config_filename(prefix)
        write_effective_config(cfg, cfg_path)
        out_files.append(cfg_path)
    except Exception as exc:
        click.echo(f"[config snapshot skipped] {exc}", err=True)
    click.echo(f"\nOutput written to: {cfg['output']['dir']}")
    for f in out_files:
        click.echo(f"  {f.name}")
    _final_summary(result, qc_report, cfg)


def _final_summary(result, qc_report, cfg) -> None:
    """Close the run with a one-glance summary, or — if nothing came out the
    other end — a warning plus the most likely reasons why."""
    n_reps = len(result.representatives)
    n_clusters = len(result.clusters)
    segmented = bool(cfg.get("segmented", {}).get("enabled"))

    if n_reps > 0:
        cluster_note = f" across {n_clusters} cluster(s)" if n_clusters else ""
        unit = "isolate(s)" if segmented else "representative sequence(s)"
        click.echo(f"\nDone — selected {n_reps} {unit}{cluster_note}.")
        if qc_report is not None and qc_report.total_input:
            msg = (
                f"  {qc_report.passed} of {qc_report.total_input} input "
                f"sequences passed basic QC"
            )
            # When later stages (protein QC, segmented, taxonomy_consistency,
            # strain-collision) trimmed further, show the post-everything
            # count so the user isn't misled into thinking ``passed`` was
            # the survivor count fed into selection.
            if (
                qc_report.final_survivors is not None
                and qc_report.final_survivors != qc_report.passed
            ):
                msg += (
                    f"; {qc_report.final_survivors} "
                    f"{qc_report.final_survivors_unit} reached selection."
                )
            else:
                msg += "."
            click.echo(msg)
        return

    # Nothing was selected — diagnose.
    click.echo("\nWARNING: no representative sequences were selected.", err=True)
    reasons: list[str] = []
    if qc_report is None or qc_report.total_input == 0:
        reasons.append(
            "No sequences were loaded — check the input FASTA path(s) and, "
            "if headers are unusual, pass --source explicitly."
        )
    elif qc_report.passed == 0:
        bits = []
        if qc_report.removed_duplicates:
            bits.append(f"{qc_report.removed_duplicates} duplicate(s)")
        if qc_report.removed_length:
            bits.append(f"{qc_report.removed_length} on length")
        if qc_report.removed_ambiguous:
            bits.append(f"{qc_report.removed_ambiguous} on ambiguous chars")
        if qc_report.removed_annotation:
            bits.append(f"{qc_report.removed_annotation} on annotation keywords")
        if qc_report.removed_proteins:
            bits.append(f"{qc_report.removed_proteins} on protein count")
        if qc_report.removed_taxonomy_mismatch:
            bits.append(
                f"{qc_report.removed_taxonomy_mismatch} on taxonomy mismatch"
            )
        detail = f" ({', '.join(bits)})" if bits else ""
        reasons.append(
            f"QC removed all {qc_report.total_input} input sequences{detail} — "
            "loosen the relevant qc.* settings (genome_length_filter, "
            "ambiguous_threshold, annotation_filter keywords, protein_annotation)."
        )
    elif segmented and qc_report.removed_incomplete_isolates:
        reasons.append(
            f"{qc_report.passed} sequences passed basic QC, but the segmented "
            f"completeness/length filter dropped everything "
            f"({qc_report.removed_incomplete_isolates} removed) — no isolate "
            "had all expected segments. Check isolate_regex, the segment "
            "names/aliases, and any segment_lengths bounds."
        )
    else:
        reasons.append(
            f"{qc_report.passed} sequences passed basic QC but selection produced "
            "nothing — check that MMseqs2 is on PATH and that the grouping "
            "field actually has values (taxonomic/host/geographic modes fall "
            "back to a single 'Unknown' group without metadata resolution)."
        )

    for r in reasons:
        click.echo(f"  - {r}", err=True)
    click.echo("  See the run log for the full QC and selection breakdown.", err=True)


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def _as_external_tool_error(exc: BaseException) -> Optional[BaseException]:
    """Return ``exc`` when it's a known missing/failed external-tool error,
    else ``None``.

    Only the clustering-backend errors can reach the top-level boundary
    uncaught — the phylo / HMM / plot tools already soft-fail inside
    ``_write_output`` and ``_run_hmm_scan``. Imported lazily so a broken
    optional dependency never breaks plain ``repseq --help``.
    """
    try:
        from .clustering.cdhit import CDHitError
        from .clustering.mmseqs2 import MMseqs2Error
    except Exception:
        return None
    return exc if isinstance(exc, (MMseqs2Error, CDHitError)) else None


class _RepseqGroup(click.Group):
    """Click group that renders known user errors as friendly one-liners.

    A :class:`~repseq.errors.RepseqError` (config / input problems) or a
    missing-external-tool failure (MMseqs2 / cd-hit not on PATH) is printed
    as a single ``Error: ...`` line to stderr and exits 1, with **no**
    traceback. Every other exception propagates unchanged — those are
    either click's own ``UsageError`` (already friendly) or a genuine bug
    whose traceback we deliberately keep so it gets reported.
    """

    def main(self, *args, **kwargs):  # type: ignore[override]
        try:
            return super().main(*args, **kwargs)
        except RepseqError as e:
            click.echo(f"Error: {e}", err=True)
            sys.exit(1)
        except Exception as e:  # noqa: BLE001 — narrowed immediately below
            tool_err = _as_external_tool_error(e)
            if tool_err is None:
                raise  # genuine bug → full traceback (by design)
            click.echo(f"Error: {tool_err}", err=True)
            click.echo(
                "       Run 'repseq doctor' to check which external tools "
                "are installed and on your PATH.",
                err=True,
            )
            sys.exit(1)


@click.group(cls=_RepseqGroup)
@click.version_option(version=__version__, prog_name="repseq")
def main():
    """repseq — representative sequence selection for large bioinformatics datasets."""
    pass


# ---------------------------------------------------------------------------
# doctor — self-test
# ---------------------------------------------------------------------------

@main.command("doctor")
@click.option("--config", "-c", "config_path", default=None,
              help="Path to YAML config file to validate (optional).")
@click.option("--no-network", is_flag=True, default=False,
              help="Skip the NCBI and UniProt reachability checks.")
def run_doctor_cmd(config_path, no_network):
    """Self-test: check dependencies, external tools, network, and config.

    Run this whenever something stops working or after a fresh install.
    Each check is tagged [OK], [WARN] (optional piece missing), or [FAIL]
    (something required is broken). repseq exits non-zero only if any
    [FAIL] is reported.
    """
    from .doctor import run_doctor
    cfg = load_config(config_path)
    report = run_doctor(cfg, config_path=config_path, no_network=no_network)
    click.echo(report.render(__version__))
    sys.exit(1 if report.has_failures else 0)


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
               segmented, dry_run, no_resolve, overflow, plot, phylo, per_protein_phylo, per_segment_phylo, source_override, alphabet_for_clustering, concatenate_markers, fast, verbose, pre_cluster_tree, protect_ids, pin_ids, newick,threshold, n_select):
    """Global mode: cluster at a threshold or select N diverse sequences."""
    if threshold is None and n_select is None:
        raise click.UsageError("Provide --threshold or --n-select.")

    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed,
                             alphabet_for_clustering=alphabet_for_clustering, concatenate_markers=concatenate_markers, fast=fast,
                             verbose=verbose, protect_ids=protect_ids, pin_ids=pin_ids, newick=newick)
    if segmented:
        cfg["segmented"]["enabled"] = True

    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    ncbi = _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    _populate_genbank_isolate_segment(sequences, cfg, ncbi)
    sequences = _run_protein_qc(sequences, cfg, qc_report, ncbi)
    sequences = _filter_taxonomy_consistent(sequences, cfg, qc_report)
    sequences = _run_protein_quality_qc(sequences, cfg, qc_report, ncbi)
    sequences = _run_hmm_qc(sequences, cfg, qc_report, ncbi)
    sequences = _setup_protein_alphabet(sequences, cfg, qc_report, ncbi)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)
    click.echo(qc_report.summary())

    from .modes.global_mode import GlobalMode
    mode = GlobalMode(cfg, threshold=threshold, n_select=n_select)
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names, pre_clustering_sequences=sequences, plot=plot, phylo=phylo, per_protein_phylo=per_protein_phylo, per_segment_phylo=per_segment_phylo, pre_cluster_tree=pre_cluster_tree)


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
                   segmented, dry_run, no_resolve, overflow, plot, phylo, per_protein_phylo, per_segment_phylo, source_override, alphabet_for_clustering, concatenate_markers, fast, verbose, pre_cluster_tree, protect_ids, pin_ids, newick,rank, n_per_group):
    """Taxonomic mode 1: N representatives per taxonomic rank group."""
    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed,
                             alphabet_for_clustering=alphabet_for_clustering, concatenate_markers=concatenate_markers, fast=fast,
                             verbose=verbose, protect_ids=protect_ids, pin_ids=pin_ids, newick=newick)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    ncbi = _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    _populate_genbank_isolate_segment(sequences, cfg, ncbi)
    sequences = _run_protein_qc(sequences, cfg, qc_report, ncbi)
    sequences = _filter_taxonomy_consistent(sequences, cfg, qc_report)
    sequences = _run_protein_quality_qc(sequences, cfg, qc_report, ncbi)
    sequences = _run_hmm_qc(sequences, cfg, qc_report, ncbi)
    sequences = _setup_protein_alphabet(sequences, cfg, qc_report, ncbi)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)
    click.echo(qc_report.summary())

    from .modes.taxonomic1 import TaxonomicMode1
    mode = TaxonomicMode1(cfg, rank=rank, n_per_group=n_per_group, overflow=overflow)
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names, pre_clustering_sequences=sequences, plot=plot, phylo=phylo, per_protein_phylo=per_protein_phylo, per_segment_phylo=per_segment_phylo, pre_cluster_tree=pre_cluster_tree)


# ---------------------------------------------------------------------------
# run taxonomic2
# ---------------------------------------------------------------------------

@main.command("taxonomic2")
@_shared_options
@click.option("--rank-levels", "-r", required=True,
              help='JSON list of {rank, n_per_group} dicts. E.g. \'[{"rank":"family","n_per_group":20},{"rank":"genus","n_per_group":5}]\'')
def run_taxonomic2(config_path, input_paths, output_dir, prefix, threads, seed,
                   segmented, dry_run, no_resolve, overflow, plot, phylo, per_protein_phylo, per_segment_phylo, source_override, alphabet_for_clustering, concatenate_markers, fast, verbose, pre_cluster_tree, protect_ids, pin_ids, newick,rank_levels):
    """Taxonomic mode 2: hierarchical multi-rank nested clustering."""
    import json as _json
    try:
        rank_levels_parsed = _json.loads(rank_levels)
    except Exception:
        raise click.UsageError("--rank-levels must be valid JSON.")

    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed,
                             alphabet_for_clustering=alphabet_for_clustering, concatenate_markers=concatenate_markers, fast=fast,
                             verbose=verbose, protect_ids=protect_ids, pin_ids=pin_ids, newick=newick)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    ncbi = _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    _populate_genbank_isolate_segment(sequences, cfg, ncbi)
    sequences = _run_protein_qc(sequences, cfg, qc_report, ncbi)
    sequences = _filter_taxonomy_consistent(sequences, cfg, qc_report)
    sequences = _run_protein_quality_qc(sequences, cfg, qc_report, ncbi)
    sequences = _run_hmm_qc(sequences, cfg, qc_report, ncbi)
    sequences = _setup_protein_alphabet(sequences, cfg, qc_report, ncbi)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)
    click.echo(qc_report.summary())

    from .modes.taxonomic2 import TaxonomicMode2
    mode = TaxonomicMode2(cfg, rank_levels=rank_levels_parsed, overflow=overflow)
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names, pre_clustering_sequences=sequences, plot=plot, phylo=phylo, per_protein_phylo=per_protein_phylo, per_segment_phylo=per_segment_phylo, pre_cluster_tree=pre_cluster_tree)


# ---------------------------------------------------------------------------
# run host
# ---------------------------------------------------------------------------

@main.command("host")
@_shared_options
@click.option("--n-per-host", "-n", required=True, type=int,
              help="Target representatives per host organism.")
def run_host(config_path, input_paths, output_dir, prefix, threads, seed,
             segmented, dry_run, no_resolve, overflow, plot, phylo, per_protein_phylo, per_segment_phylo, source_override, alphabet_for_clustering, concatenate_markers, fast, verbose, pre_cluster_tree, protect_ids, pin_ids, newick,n_per_host):
    """Host-stratified mode: N representatives per host organism."""
    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed,
                             alphabet_for_clustering=alphabet_for_clustering, concatenate_markers=concatenate_markers, fast=fast,
                             verbose=verbose, protect_ids=protect_ids, pin_ids=pin_ids, newick=newick)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    ncbi = _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    _populate_genbank_isolate_segment(sequences, cfg, ncbi)
    sequences = _run_protein_qc(sequences, cfg, qc_report, ncbi)
    sequences = _filter_taxonomy_consistent(sequences, cfg, qc_report)
    sequences = _run_protein_quality_qc(sequences, cfg, qc_report, ncbi)
    sequences = _run_hmm_qc(sequences, cfg, qc_report, ncbi)
    sequences = _setup_protein_alphabet(sequences, cfg, qc_report, ncbi)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)
    click.echo(qc_report.summary())

    from .modes.host_mode import HostMode
    mode = HostMode(cfg, n_per_host=n_per_host, overflow=overflow)
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names, pre_clustering_sequences=sequences, plot=plot, phylo=phylo, per_protein_phylo=per_protein_phylo, per_segment_phylo=per_segment_phylo, pre_cluster_tree=pre_cluster_tree)


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
             segmented, dry_run, no_resolve, overflow, plot, phylo, per_protein_phylo, per_segment_phylo, source_override, alphabet_for_clustering, concatenate_markers, fast, verbose, pre_cluster_tree, protect_ids, pin_ids, newick,n_per_window, window):
    """Time-stratified mode: N representatives per time window."""
    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed,
                             alphabet_for_clustering=alphabet_for_clustering, concatenate_markers=concatenate_markers, fast=fast,
                             verbose=verbose, protect_ids=protect_ids, pin_ids=pin_ids, newick=newick)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    ncbi = _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    _populate_genbank_isolate_segment(sequences, cfg, ncbi)
    sequences = _run_protein_qc(sequences, cfg, qc_report, ncbi)
    sequences = _filter_taxonomy_consistent(sequences, cfg, qc_report)
    sequences = _run_protein_quality_qc(sequences, cfg, qc_report, ncbi)
    sequences = _run_hmm_qc(sequences, cfg, qc_report, ncbi)
    sequences = _setup_protein_alphabet(sequences, cfg, qc_report, ncbi)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)
    click.echo(qc_report.summary())

    from .modes.time_mode import TimeMode
    mode = TimeMode(cfg, n_per_window=n_per_window, window=window, overflow=overflow)
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names, pre_clustering_sequences=sequences, plot=plot, phylo=phylo, per_protein_phylo=per_protein_phylo, per_segment_phylo=per_segment_phylo, pre_cluster_tree=pre_cluster_tree)


# ---------------------------------------------------------------------------
# run geographic
# ---------------------------------------------------------------------------

@main.command("geographic")
@_shared_options
@click.option("--n-per-country", "-n", required=True, type=int,
              help="Target representatives per country.")
def run_geographic(config_path, input_paths, output_dir, prefix, threads, seed,
                   segmented, dry_run, no_resolve, overflow, plot, phylo, per_protein_phylo, per_segment_phylo, source_override, alphabet_for_clustering, concatenate_markers, fast, verbose, pre_cluster_tree, protect_ids, pin_ids, newick,n_per_country):
    """Geographic mode: N representatives per country."""
    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed,
                             alphabet_for_clustering=alphabet_for_clustering, concatenate_markers=concatenate_markers, fast=fast,
                             verbose=verbose, protect_ids=protect_ids, pin_ids=pin_ids, newick=newick)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    ncbi = _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    _populate_genbank_isolate_segment(sequences, cfg, ncbi)
    sequences = _run_protein_qc(sequences, cfg, qc_report, ncbi)
    sequences = _filter_taxonomy_consistent(sequences, cfg, qc_report)
    sequences = _run_protein_quality_qc(sequences, cfg, qc_report, ncbi)
    sequences = _run_hmm_qc(sequences, cfg, qc_report, ncbi)
    sequences = _setup_protein_alphabet(sequences, cfg, qc_report, ncbi)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)
    click.echo(qc_report.summary())

    from .modes.geographic_mode import GeographicMode
    mode = GeographicMode(cfg, n_per_country=n_per_country, overflow=overflow)
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names, pre_clustering_sequences=sequences, plot=plot, phylo=phylo, per_protein_phylo=per_protein_phylo, per_segment_phylo=per_segment_phylo, pre_cluster_tree=pre_cluster_tree)


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
               segmented, dry_run, no_resolve, overflow, plot, phylo, per_protein_phylo, per_segment_phylo, source_override, alphabet_for_clustering, concatenate_markers, fast, verbose, pre_cluster_tree, protect_ids, pin_ids, newick,field, n_per_group,
               metadata_table, field_regex):
    """Custom metadata mode: group by any field or metadata table column."""
    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed,
                             alphabet_for_clustering=alphabet_for_clustering, concatenate_markers=concatenate_markers, fast=fast,
                             verbose=verbose, protect_ids=protect_ids, pin_ids=pin_ids, newick=newick)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    ncbi = _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    _populate_genbank_isolate_segment(sequences, cfg, ncbi)
    sequences = _run_protein_qc(sequences, cfg, qc_report, ncbi)
    sequences = _filter_taxonomy_consistent(sequences, cfg, qc_report)
    sequences = _run_protein_quality_qc(sequences, cfg, qc_report, ncbi)
    sequences = _run_hmm_qc(sequences, cfg, qc_report, ncbi)
    sequences = _setup_protein_alphabet(sequences, cfg, qc_report, ncbi)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)
    click.echo(qc_report.summary())

    from .modes.custom_mode import CustomMode
    mode = CustomMode(
        cfg, field=field, n_per_group=n_per_group,
        metadata_table_path=metadata_table,
        field_regex=field_regex,
        overflow=overflow,
    )
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names, pre_clustering_sequences=sequences, plot=plot, phylo=phylo, per_protein_phylo=per_protein_phylo, per_segment_phylo=per_segment_phylo, pre_cluster_tree=pre_cluster_tree)


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
               segmented, dry_run, no_resolve, overflow, plot, phylo, per_protein_phylo, per_segment_phylo, source_override, alphabet_for_clustering, concatenate_markers, fast, verbose, pre_cluster_tree, protect_ids, pin_ids, newick,fields, n_per_group,
               metadata_table):
    """Hybrid mode: multi-dimensional stratification (e.g. genus × host × year)."""
    field_list = [f.strip() for f in fields.split(",")]
    cfg = _load_and_validate(config_path, output_dir, prefix, threads, seed,
                             alphabet_for_clustering=alphabet_for_clustering, concatenate_markers=concatenate_markers, fast=fast,
                             verbose=verbose, protect_ids=protect_ids, pin_ids=pin_ids, newick=newick)
    if segmented:
        cfg["segmented"]["enabled"] = True
    if dry_run:
        click.echo("Dry run — config valid. Exiting.")
        return

    sequences = _load_sequences(input_paths, source_override)
    ncbi = _resolve_metadata(sequences, cfg, no_resolve)
    sequences, qc_report = _run_qc(sequences, cfg)
    _populate_genbank_isolate_segment(sequences, cfg, ncbi)
    sequences = _run_protein_qc(sequences, cfg, qc_report, ncbi)
    sequences = _filter_taxonomy_consistent(sequences, cfg, qc_report)
    sequences = _run_protein_quality_qc(sequences, cfg, qc_report, ncbi)
    sequences = _run_hmm_qc(sequences, cfg, qc_report, ncbi)
    sequences = _setup_protein_alphabet(sequences, cfg, qc_report, ncbi)
    sequences, complete_isolates, segment_names = _handle_segmented(sequences, cfg, qc_report)
    click.echo(qc_report.summary())

    from .modes.hybrid_mode import HybridMode
    mode = HybridMode(
        cfg, fields=field_list, n_per_group=n_per_group,
        overflow=overflow, metadata_table_path=metadata_table,
    )
    result = mode.run(sequences)
    result.qc_report = qc_report

    _write_output(result, qc_report, cfg, input_paths, complete_isolates, segment_names, pre_clustering_sequences=sequences, plot=plot, phylo=phylo, per_protein_phylo=per_protein_phylo, per_segment_phylo=per_segment_phylo, pre_cluster_tree=pre_cluster_tree)


# ---------------------------------------------------------------------------
# replay — re-materialise representatives from a lockfile
# ---------------------------------------------------------------------------

from .replay import replay_command as _replay_command  # noqa: E402

main.add_command(_replay_command)


# ---------------------------------------------------------------------------
# cache management
# ---------------------------------------------------------------------------

@main.group("cache")
def cache_group():
    """Manage the taxonomy/metadata cache."""
    pass


@cache_group.command("stats")
@click.option("--config", "-c", "config_path", default=None,
              help="Path to YAML config file (locates the cache directory).")
def cache_stats(config_path):
    """Show cache statistics: total entries, size, and a per-source breakdown."""
    cfg = load_config(config_path)
    cache = TaxonomyCache(cfg["cache_dir"])
    stats = cache.stats()
    click.echo(f"Cache: {stats['db_path']}")
    click.echo(f"Total entries : {stats['total_entries']}")
    click.echo(f"DB size       : {stats['db_size_mb']} MB")
    for src, count in stats.get("by_source", {}).items():
        click.echo(f"  {src}: {count}")


@cache_group.command("clear")
@click.option("--config", "-c", "config_path", default=None,
              help="Path to YAML config file (locates the cache directory).")
@click.option("--source", default=None,
              help=(
                  "Clear only this source. Common values: 'ncbi_taxonomy' "
                  "(lineages), 'ncbi_proteins' (GenBank CDS records), "
                  "'ncbi_nuc_seq' (nucleotide bodies for `repseq replay`), "
                  "'uniprot' (UniProt entries), 'hmmscan' (HMM hit lists). "
                  "Run 'repseq cache stats' to see which sources are present."
              ))
@click.confirmation_option(prompt="This will delete cached data. Continue?")
def cache_clear(config_path, source):
    """Clear the cache (all sources, or one with --source)."""
    cfg = load_config(config_path)
    cache = TaxonomyCache(cfg["cache_dir"])
    n = cache.clear(source)
    click.echo(f"Deleted {n} cache entries.")


@cache_group.command("purge-expired")
@click.option("--config", "-c", "config_path", default=None,
              help="Path to YAML config file (locates the cache directory).")
def cache_purge(config_path):
    """Remove expired cache entries (older than taxonomy.cache_ttl_days)."""
    cfg = load_config(config_path)
    cache = TaxonomyCache(cfg["cache_dir"])
    n = cache.purge_expired()
    click.echo(f"Purged {n} expired cache entries.")


# ---------------------------------------------------------------------------
# stats — pre-flight inspection of an input FASTA
# ---------------------------------------------------------------------------

@main.command("stats")
@click.option("--input", "-i", "input_paths", multiple=True, required=True,
              help="Input FASTA file(s). Repeat for multiple files.")
@click.option("--source", "source_override",
              type=click.Choice(["auto", "uniprot", "ncbi", "ncbi_virus"]),
              default="auto",
              help="Force input source instead of auto-detecting from headers.")
@click.option("--top-n", default=10, type=int,
              help="Top N taxa per rank in the breakdown (default 10).")
@click.option("--resolve/--no-resolve", default=False,
              help=(
                  "Resolve taxonomy via NCBI/UniProt before stats — slower but "
                  "fills in missing organism/taxonomy. Default: don't resolve "
                  "(use header-derived metadata only — fast, no network)."
              ))
@click.option("--config", "-c", "config_path", default=None,
              help="Path to YAML config file (only needed with --resolve).")
def run_stats(input_paths, source_override, top_n, resolve, config_path):
    """Pre-flight inspection of an input FASTA — count, taxonomy spread,
    length distribution, missing-metadata fractions. No clustering, no
    output files. Use this to decide whether your input is shaped the way
    you expect before committing to a long run.
    """
    sequences = _load_sequences(input_paths, source_override)
    if not sequences:
        click.echo("No sequences loaded.", err=True)
        sys.exit(1)

    if resolve:
        cfg = load_config(config_path)
        _resolve_metadata(sequences, cfg, no_resolve=False)

    _print_input_stats(sequences, top_n=top_n)


def _print_input_stats(sequences, *, top_n: int = 10) -> None:
    """Render the stats panel to stdout. Pure presentation — no I/O
    beyond ``click.echo``."""
    from collections import Counter
    import statistics

    n = len(sequences)
    click.echo("")
    click.echo("=" * 60)
    click.echo(f"Input statistics — {n:,} sequence(s)")
    click.echo("=" * 60)

    # Source breakdown.
    src_counts: Counter = Counter(s.source.value for s in sequences)
    click.echo("\nSource:")
    for src, c in src_counts.most_common():
        click.echo(f"  {src:<15} {c:>8,}  ({100*c/n:.1f}%)")

    # Length distribution.
    lengths = [s.length for s in sequences if s.length]
    if lengths:
        click.echo("\nSequence length (nucleotides):")
        click.echo(f"  min     {min(lengths):>10,}")
        click.echo(f"  max     {max(lengths):>10,}")
        click.echo(f"  median  {int(statistics.median(lengths)):>10,}")
        if len(lengths) >= 4:
            q1, _q2, q3 = statistics.quantiles(lengths, n=4)
            click.echo(f"  Q3-Q1   {int(round(q3 - q1)):>10,}")

    # Missing-metadata fractions — the actionable "what's empty" report.
    fields = (
        ("organism", lambda s: s.organism),
        ("accession", lambda s: s.accession),
        ("strain", lambda s: s.strain),
        ("host", lambda s: s.host),
        ("subtype", lambda s: s.subtype),
        ("collection_date", lambda s: s.collection_date),
        ("country", lambda s: s.country),
        ("segment", lambda s: s.segment),
        ("isolate_id", lambda s: s.isolate_id),
        ("taxonomy", lambda s: s.taxonomy),
    )
    click.echo("\nMetadata coverage (populated / total):")
    name_w = max(len(name) for name, _ in fields)
    for name, getter in fields:
        c = sum(1 for s in sequences if getter(s))
        click.echo(
            f"  {name:<{name_w}}  {c:>6,} / {n:<6,}  ({100*c/n:.1f}%)"
        )

    # RefSeq / reviewed flags.
    n_refseq = sum(1 for s in sequences if s.is_refseq)
    n_reviewed = sum(1 for s in sequences if s.is_reviewed)
    click.echo("\nQuality flags:")
    click.echo(f"  is_refseq        {n_refseq:>6,}  ({100*n_refseq/n:.1f}%)")
    click.echo(f"  is_reviewed      {n_reviewed:>6,}  ({100*n_reviewed/n:.1f}%)")

    # Taxonomy breakdown — top N per rank with any populated value.
    rank_names = ("species", "genus", "family", "order", "class")
    click.echo(f"\nTaxonomy (top {top_n} per rank by count):")
    any_tax = False
    for rank in rank_names:
        c: Counter = Counter()
        for s in sequences:
            if s.taxonomy is None:
                continue
            v = s.taxonomy.get_rank(rank)
            if v:
                c[v] += 1
        if not c:
            continue
        any_tax = True
        n_distinct = len(c)
        suffix = (
            f" ({n_distinct} distinct, top {top_n} shown)"
            if n_distinct > top_n
            else f" ({n_distinct} distinct)"
        )
        click.echo(f"  {rank}{suffix}:")
        for name, count in c.most_common(top_n):
            click.echo(f"    {count:>6,}  {name}")
    if not any_tax:
        click.echo(
            "  (no taxonomy resolved on input; rerun with --resolve "
            "to fetch from NCBI/UniProt)"
        )

    click.echo("")


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

    # Whole-genome length filter — non-segmented only, absolute nt bounds.
    glf: dict = {"enabled": False, "min": None, "max": None}
    if not cfg.get("segmented", {}).get("enabled") and click.confirm(
        "Enable whole-genome length filter (non-segmented only)?", default=False
    ):
        glf["enabled"] = True
        glf["min"] = click.prompt(
            "Minimum genome length in nt (blank for none)", default=None, type=int,
        )
        glf["max"] = click.prompt(
            "Maximum genome length in nt (blank for none)", default=None, type=int,
        )
    cfg["qc"]["genome_length_filter"] = glf

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

    if click.confirm("Enable protein-annotation QC (fetches CDS counts from NCBI)?", default=False):
        cfg["qc"]["protein_annotation"] = {
            "enabled": True,
            "min_proteins": click.prompt(
                "Minimum number of annotated proteins per sequence", default=1, type=int
            ),
        }
    else:
        cfg["qc"]["protein_annotation"] = {"enabled": False}

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
        virus_def: dict = {
            "expected_segments": n_seg,
            "segments": segments,
            "isolate_regex": isolate_regex,
        }

        if click.confirm(
            "Configure expected protein counts per segment?", default=False
        ):
            expected_proteins: dict = {}
            click.echo(
                "  Enter a single integer (exact count) or comma-separated integers\n"
                "  (any of those counts accepted, e.g. '1,2' for PB1 ± PB1-F2).\n"
                "  Press Enter to skip a segment."
            )
            for seg in segments:
                raw = click.prompt(f"  {seg}", default="").strip()
                if raw:
                    parts = [p.strip() for p in raw.split(",") if p.strip()]
                    if len(parts) == 1:
                        expected_proteins[seg] = int(parts[0])
                    else:
                        expected_proteins[seg] = [int(p) for p in parts]
            if expected_proteins:
                virus_def["expected_proteins_per_segment"] = expected_proteins

        if click.confirm("Configure length bounds per segment?", default=False):
            segment_lengths: dict = {}
            click.echo(
                "  Enter min and/or max nucleotide length per segment.\n"
                "  Enter 0 (or press Enter) to leave a bound unset."
            )
            for seg in segments:
                mn = click.prompt(f"  {seg} min length", default=0, type=int)
                mx = click.prompt(f"  {seg} max length", default=0, type=int)
                bounds: dict = {}
                if mn:
                    bounds["min"] = mn
                if mx:
                    bounds["max"] = mx
                if bounds:
                    segment_lengths[seg] = bounds
            if segment_lengths:
                virus_def["segment_lengths"] = segment_lengths

        cfg["segmented"] = {
            "enabled": False,
            "virus": virus_name,
            "viruses": {virus_name: virus_def},
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
    backend = click.prompt(
        "Clustering backend",
        type=click.Choice(["mmseqs2", "cdhit"]),
        default="mmseqs2",
    )
    cfg["clustering"] = {"backend": backend}
    if backend == "mmseqs2":
        cfg["clustering"]["mmseqs2_mode"] = click.prompt(
            "MMseqs2 mode", type=click.Choice(["easy-linclust", "easy-cluster"]),
            default="easy-linclust",
        )
        cfg["clustering"]["coverage"] = click.prompt(
            "Coverage threshold", default=0.8, type=float,
        )
        cfg["clustering"]["coverage_mode"] = click.prompt(
            "Coverage mode (0-4)", default=0, type=int,
        )
    else:
        # cd-hit: defaults are reasonable; leave word_size on auto-pick.
        cfg["clustering"]["cdhit"] = {
            "binary": None,
            "word_size": None,
            "coverage": click.prompt(
                "cd-hit coverage (-aS, used only when global_alignment=false)",
                default=0.8, type=float,
            ),
            "global_alignment": click.confirm(
                "Use global identity (-G 1)?", default=True,
            ),
            "accurate": click.confirm(
                "Use accurate mode (-g 1, slower)?", default=False,
            ),
            "memory_mb": click.prompt(
                "Memory cap MB (0 = unlimited)", default=0, type=int,
            ),
            "extra_args": [],
        }

    # Output
    cfg["output"] = {
        "dir": click.prompt("Output directory", default="./repseq_output"),
        "prefix": click.prompt("Output file prefix", default="repseq"),
    }

    out_path = Path(output)
    with open(out_path, "w") as fh:
        yaml.dump(cfg, fh, default_flow_style=False, sort_keys=False)

    # Insert an explanatory comment before `enabled: false` in the segmented
    # block — yaml.dump cannot emit comments, so we post-process the file.
    text = out_path.read_text()
    text = text.replace(
        "segmented:\n  enabled: false\n",
        "segmented:\n"
        "  # Set enabled to true here, or pass --segmented on the command line.\n"
        "  # The virus definition below is stored so you can activate it per-run\n"
        "  # without re-editing this file.\n"
        "  enabled: false\n",
    )
    out_path.write_text(text)

    click.echo(f"\nConfig written to: {out_path}")


if __name__ == "__main__":
    main()
