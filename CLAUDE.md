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
├── output/{writer,report}.py ← FASTA + TSV/text report writers
└── viz/clustering_plot.py  ← UMAP scatter (optional, [viz] extras)
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
3. `_run_qc` produces a passing-list + `QCReport`. In segmented mode it
   skips both the whole-pool length filter and exact-duplicate removal
   (`QCReport.length_filter_skipped` / `dedup_skipped`) — see step 5.
4. `_run_protein_qc` (optional) — if `qc.protein_annotation.enabled` or
   `virus.expected_proteins_per_segment` is set, fetches CDS features via
   `NCBITaxonomy.fetch_proteins_batch` (200 accessions per request, cached
   in SQLite under source `ncbi_proteins`), populates `seq.proteins`, and
   drops sequences failing the count check.
5. If `segmented.enabled`, `_handle_segmented` runs `filter_complete_isolates`
   then `build_concatenated_sequences`, replacing the sequence list with one
   concatenated sequence per complete isolate. Exact-duplicate removal is then
   applied here, on the *concatenated* sequences (not the segment pool — a
   conserved segment shared between two distinct isolates must not knock
   either out as incomplete); de-duplicated isolates are also pruned from
   `complete_isolates` so the per-segment output files stay consistent.
6. The chosen mode runs (`modes/<mode>.py`), returning a `RunResult` with
   `representatives` + `clusters`.
7. `write_results` writes the FASTA(s); `write_all_reports` writes the
   plain-text/TSV reports — including `{prefix}_group_counts.tsv`
   (one row per stratification group: `grouping, group, n_before,
   n_after, clustered, cutoff` — populated from `RunResult.group_stats`,
   which every mode fills in), `{prefix}_isolate_proteins.tsv`
   (one row per CDS per passing isolate, segmented mode — columns:
   `protein_id, product, length, isolate_id, segment, segment_length,
   accession, species, subgenus, genus, subfamily, family, suborder,
   order, subclass, class`; the four sub-ranks are populated only from
   the resolver lineage map and commonly blank for viruses) and
   `{prefix}_proteins.fasta` (amino-acid sequences for all proteins of
   the selected representatives). The protein FASTA is reconstructed
   from the same cached GenBank records — no extra network calls.
8. If `--plot` is passed, `viz.clustering_plot.write_clustering_plot`
   embeds the clustered sequences with UMAP on k-mer Jaccard distance
   and writes `{prefix}_clustering.png`. Cost-bounded by a default
   2000-point subsample (representatives always kept); skipped when
   the run produced no clusters. Requires the `[viz]` extras
   (`matplotlib` + `umap-learn`) — `ImportError` is surfaced gracefully.

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
- **Output directory must be empty**: `_load_and_validate` calls
  `_check_output_dir`, which aborts (exit 1) if `output.dir` already
  exists and is non-empty — runs never overwrite or mix into prior
  results. Tests that exercise the writers point at a fresh `tmp_path`.

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

`v0.5.6`. All 8 modes structurally complete, optional protein-annotation
QC (batched GenBank fetch + per-segment count check), a protein-FASTA
output reconstructed from cached records, per-segment nucleotide length
bounds (`segment_lengths` in virus config, applied after completeness
filter), and an optional UMAP scatter of the clustering result (`--plot`,
behind the `[viz]` extras: `matplotlib` + `umap-learn`). In segmented mode
the whole-pool QC length filter is skipped automatically (mixed segment
lengths make a pooled median meaningless and would drop short segments,
leaving every isolate incomplete) — use `segment_lengths` instead.
Exact-duplicate removal is likewise skipped on the segment pool and
applied instead to the concatenated per-isolate sequences in
`_handle_segmented` — a segment conserved across two distinct isolates
must not delete either isolate by leaving it "incomplete". Every
run ends with a one-line CLI summary (representatives/clusters selected, QC
pass rate) or, when nothing is selected, a stderr warning naming the most
likely cause (`cli._final_summary`). Every mode also records per-group
before/after counts (`RunResult.group_stats`) written to
`{prefix}_group_counts.tsv`. Taxonomic lineage is resolved via NCBI
`efetch` XML (`ncbi._parse_taxonomy_xml`) — the taxonomy *esummary*
endpoint carries no lineage for viruses, which used to group every
viral sequence under "Unknown". Offline regression tests cover IO,
QC, selector, segmented logic (including per-segment length filtering),
cache TTL, config validation, diversity selection, resolver fallback, mode
dispatch (clustering mocked), protein parser + filter, FASTA writer, segment
aliases, the closing CLI summary, per-group count reporting, and the
clustering plot (auto-skipped if `umap-learn` is missing). NCBI Entrez paths have been
live-tested against the H1N1 RefSeq genome (8 segments, 11 proteins).

A full pipeline audit (see git history around this note) fixed, with
regression tests: clustering ID round-trip for UniProt/CONCAT inputs,
RefSeq-nucleotide accession misrouting, the `MAG:` keyword filter, an
**inverted binary-search direction** (the threshold↔cluster-count
relationship was backwards, so `n-per-group` searches walked away from
the target), NCBI host/country/date now harvested from esummary
`subtype`/`subname` with the DB taking precedence over header heuristics,
a length-robust containment distance for diversity selection, thread-safe
caches and rate limiters, and assorted QC/parsing corrections.

`v0.5.5` hardens every place where header or identifier text crosses a
format boundary. The trigger was a Peribunyaviridae run where a 73-seq
genus came out of clustering as 73 reps at `cutoff = 1.0000`: isolate
names like `yaba-7 virus strain yaba 7` contain whitespace, so the
normalised id (and the derived `CONCAT|<isolate>` `seq.id`) carried
spaces; MMseqs2 truncates a FASTA header token at the first whitespace,
the parser then could not match cluster-TSV rows back to inputs, and
the binary search misread the empty cluster list as "≤ target" — falling
through to the final `best_reps = list(sequences)` fallback at `hi = 1.0`.
Fixes, with regression tests:
`_normalise_isolate_id` replaces both whitespace **and** the pipe
separator with `_` so `seq.id` survives MMseqs2 truncation and
`split("|")[1]` CONCAT parsing; `run_clustering` raises
`MMseqs2Error("Cluster round-trip mismatch...")` if the parsed cluster
membership does not account for every input; both FASTA writers
(`io/fasta.write_fasta`, `clustering/mmseqs2._write_id_fasta`)
collapse embedded `\r`/`\n` in headers to a space so a malformed
header cannot inject a phantom record; every TSV writer routes free-
text fields through a new `output/report._tsv_safe` helper that
neutralises tab and line-break characters so a stray value cannot
shift columns or rows. A 200-case brittleness probe
(`tests/test_unusual_characters.py`) exercises every (surface,
character) pair across whitespace and the punctuation classes most
common in viral metadata. 362 offline regression tests total.

`v0.5.6` enriches `{prefix}_isolate_proteins.tsv`. The previous six
columns (`isolate_id, segment, accession, protein_id, product, length`)
are reordered and extended to sixteen: `protein_id, product, length,
isolate_id, segment, segment_length, accession, species, subgenus,
genus, subfamily, family, suborder, order, subclass, class`. The
protein columns lead so the table reads as "what gene, on which
segment of which isolate, in what organism." `segment_length` is the
nucleotide length of the parent segment (`seq.length`). Taxonomy is
read via `TaxonomyInfo.get_rank(rank)`, which checks the standard
fields first and falls back to the resolver's lineage map (NCBI
`LineageEx`); the four sub-ranks (`subgenus`, `subfamily`, `suborder`,
`subclass`) have no standard field and come exclusively from the
lineage map — they are commonly blank for viruses, where ICTV often
skips intermediate ranks. Two regression tests cover the lineage-
backed sub-ranks and the missing-taxonomy fallback (all 9 rank cells
blank). 364 offline regression tests total.
