# CLAUDE.md

Guidance for Claude Code (or any future agent) working on this repo.

## What this project is

`repseq` is a Python CLI for **selecting representative sequences** from
large FASTA datasets — protein or nucleotide, with a strong viral-sequence
focus. The user is a JCVI bioinformatics researcher; treat them as fluent in
the domain (taxonomy, clustering, FASTA conventions) but not necessarily in
every implementation detail, since most of the code was written by Claude.

Entry point: `repseq` (defined in `repseq/cli.py`).

## Top-level layout

```
repseq/
├── cli.py                  ← click commands, dispatch into modes
├── config.py               ← YAML defaults + validation
├── models.py               ← Sequence, Cluster, QCReport, RunResult dataclasses
├── io/fasta.py             ← parse 3 header formats, read/write FASTA
├── qc/
│   ├── pipeline.py         ← duplicates → length → ambiguous → annotation
│   └── protein_qc.py       ← NCBI-backed protein-count filter (opt-in)
├── segmented/completeness.py ← isolate grouping, completeness, concat
├── taxonomy/
│   ├── cache.py            ← SQLite TTL cache
│   ├── ncbi.py             ← Entrez API client
│   ├── uniprot.py          ← UniProt REST client
│   └── resolver.py         ← DB → header → strain-label fallback chain
├── clustering/
│   ├── mmseqs2.py          ← shells out to `mmseqs`
│   └── diversity.py        ← MaxMin selection (alignment-free, k-mer Jaccard)
├── representative/selector.py ← RefSeq > reviewed > longest priority
├── modes/                  ← one file per selection mode, all extend BaseMode
└── output/{writer,report}.py ← FASTA + JSON/text report writers
```

`config/default_config.yaml` is the documented config schema; the same dict
shape is hard-coded in `repseq/config.py:DEFAULTS`.

## Data flow

1. `_load_sequences` (cli.py) → `read_fasta` yields `Sequence` objects with
   parsed `accession`, `organism`, `is_refseq`, `is_reviewed`, etc.
2. `_resolve_metadata` builds a `MetadataResolver` and runs
   `resolve_batch` in a thread pool. Failures are non-fatal but **tracked**
   on `resolver.failures` and surfaced via stderr. Returns the `NCBITaxonomy`
   client so downstream steps can reuse it.
3. `_run_qc` produces a passing-list + `QCReport`.
4. `_run_protein_qc` (optional) — if `qc.protein_annotation.enabled` or
   `virus.expected_proteins_per_segment` is set, fetches CDS features via
   `NCBITaxonomy.fetch_proteins_batch` (200 accessions per request, cached
   in SQLite under source `ncbi_proteins`), populates `seq.proteins`, and
   drops sequences failing the count check.
5. If `segmented.enabled`, `filter_complete_isolates` then
   `build_concatenated_sequences` replace the sequence list with one
   concatenated sequence per complete isolate.
6. The chosen mode runs (`modes/<mode>.py`), returning a `RunResult` with
   `representatives` + `clusters`.
7. `write_results` writes the FASTA(s); `write_all_reports` writes the
   plain-text/TSV reports — including `{prefix}_isolate_proteins.tsv`
   (one row per CDS per passing isolate) when proteins were fetched
   in segmented mode.

## Invariants worth knowing

- **Representative priority** (`representative.priority`) must include
  `"longest"` — `validate_config` enforces this. The selector always
  uses `length` as final tiebreaker even when not listed.
- **Cluster representative swapping**: MMseqs2 picks the longest as rep by
  default; `apply_representative_selection` re-elects across all
  members so RefSeq/reviewed wins. When swapping, the old rep is appended
  to `cluster.members`.
- **Segmented mode** requires three fields per virus in config:
  `expected_segments`, `segments`, `isolate_regex`. The regex must capture
  either a named group `isolate` or group 1.
- **Taxonomy cache** keys are `(source, key)`; TTL eviction happens lazily
  on `get` (expired entries are deleted in place).
- **Threads**: only `MetadataResolver.resolve_batch` and the MMseqs2 binary
  consume `cfg["threads"]`. Internal Python work is single-threaded.

## Modes

Every mode subclasses `BaseMode` and returns a `RunResult`. Most modes
share the same shape:

1. Group sequences by some key (rank / host / decade / country / custom field).
2. For each group, either keep all (if `<= n_per_group`) or run a
   **binary-search over MMseqs2 thresholds** to land at `n_per_group`
   clusters (see `taxonomic1._binary_search_threshold` for the canonical
   implementation; other modes import or mirror this).
3. `--overflow keep` returns whatever cluster count the binary search
   produced; `--overflow trim` follows up with `select_diverse` to enforce
   an exact count.

`GlobalMode` is the exception: it runs *either* a single threshold pass
(`-T`) or one shot of `select_diverse` (`-n`).

## How to run things

```bash
# Tests (offline, mocked)
pytest tests/

# Smoke a single mode
repseq global -c config/default_config.yaml -i some.fasta -T 0.95 --no-resolve

# Dry-run config validation
repseq taxonomic1 -c my.yaml -i x.fasta --rank genus -n 5 --dry-run
```

`--no-resolve` skips all network lookups — useful for fast local iteration.

## Common gotchas

- **MMseqs2 not in PATH** ⇒ `MMseqs2Error` from `clustering/mmseqs2.py`.
  Cluster-based modes fail; `--no-resolve` does not help — only `global -n`
  (diversity-only) and modes where every group is small enough to skip
  clustering will run without it.
- **No taxonomy info** ⇒ taxonomic modes fall back to grouping under
  `"Unknown"`. Always run with metadata resolution unless you've
  pre-resolved.
- **Source auto-detection on weird headers**: prefer `--source <name>` for
  guaranteed parsing. NCBI Virus headers can resemble plain NCBI; the
  auto-detector checks for multiple `[...]` brackets as a tell.
- **Concatenated sequences in output**: in segmented mode, the
  representatives list contains synthetic `CONCAT|<isolate_id>` entries.
  `output/writer._write_segmented` expands these back into per-segment
  files.

## When making changes

- Keep the `Sequence` dataclass surface stable — it's threaded through
  every layer.
- Don't add network calls outside `taxonomy/`. Network-touching tests
  should `unittest.mock.MagicMock` the NCBI/UniProt clients (see
  `tests/test_resolver.py`).
- New modes: subclass `BaseMode`, register a new click subcommand in
  `cli.py`, and follow the existing `_load_sequences → _resolve_metadata
  → _run_qc → _handle_segmented → mode.run → _write_output` flow.
- Config additions: extend `DEFAULTS` and (if validation matters) add a
  check in `validate_config`. Document in `config/default_config.yaml`.

## Status

`v0.1.0` draft. All 8 modes structurally complete. 70 offline regression
tests cover IO, QC, selector, segmented logic, cache TTL, config
validation, diversity selection, resolver fallback + failure tracking,
and mode dispatch (clustering mocked). Not yet exercised against live
NCBI/UniProt endpoints.
