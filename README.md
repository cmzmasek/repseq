# repseq

**Representative sequence selection for large bioinformatics datasets.**

`repseq` curates a small, diverse, high-quality reference set from large FASTA
dumps (UniProt, NCBI, NCBI Virus). It handles QC, taxonomy resolution,
optional segmented-virus completeness filtering, and a choice of stratified
selection strategies — useful for phylogenetics, ML training sets, and
reference-database construction.

---

## Installation

```bash
pip install -e .
```

Requires **Python ≥ 3.10** and the external **[MMseqs2](https://github.com/soedinglab/MMseqs2)** binary in `PATH` (only needed for clustering-based modes).

```bash
# macOS
brew install mmseqs2
# Linux (conda)
conda install -c bioconda mmseqs2
```

---

## Quickstart

Generate a config:

```bash
repseq init-config -o my_config.yaml
```

Pick a mode and run:

```bash
# Cluster at 90% identity, keep representatives
repseq global -c my_config.yaml -i seqs.fasta -T 0.90

# Up to 5 representatives per genus
repseq taxonomic1 -c my_config.yaml -i seqs.fasta -r genus -n 5

# 10 representatives per host organism
repseq host -c my_config.yaml -i seqs.fasta -n 10
```

---

## Selection modes

| Command | What it does |
| --- | --- |
| `global` | One pass: cluster at threshold `-T`, or pick `-n` diverse sequences. |
| `taxonomic1` | `-n` per group at a single rank (`--rank genus`, etc.). |
| `taxonomic2` | Hierarchical multi-rank: `-r '[{"rank":"family","n_per_group":20},{"rank":"genus","n_per_group":5}]'`. |
| `host` | `-n` per host organism. |
| `time` | `-n` per time window (`--window year`, `decade`, or `5` for 5-year bins). |
| `geographic` | `-n` per country. |
| `custom` | `-n` per any field — from taxonomy, sequence attributes, or a metadata TSV. |
| `hybrid` | `-n` per multi-dimensional stratum (e.g. `--fields genus,host,decade`). |

All commands share these options: `--input/-i`, `--output-dir/-o`, `--config/-c`, `--threads`, `--seed`, `--segmented`, `--dry-run`, `--no-resolve`, `--source {auto,uniprot,ncbi,ncbi_virus}`, `--overflow {keep,trim}`.

---

## Pipeline

```
input FASTA(s)
    │
    ▼
parse + auto-detect source (UniProt / NCBI / NCBI Virus)
    │
    ▼
resolve taxonomy & metadata        (NCBI Entrez + UniProt, SQLite cached)
    │
    ▼
QC: duplicates → length → ambiguous chars → annotation keyword filter
    │
    ▼
[optional] segmented-virus completeness filter (concatenates segments per isolate)
    │
    ▼
mode-specific selection (clustering via MMseqs2 + diversity / stratification)
    │
    ▼
representative selection: RefSeq > reviewed UniProt > longest
    │
    ▼
output FASTA(s) + JSON report
```

---

## Segmented virus mode

For multi-segment viruses (e.g. influenza), `repseq` can:

1. Group sequences by isolate (regex on the header).
2. Keep only isolates that have **all** expected segments.
3. Concatenate per-isolate segments into a single sequence for clustering.
4. Write back both the concatenated FASTA and one FASTA per segment.

Configure under `segmented:` and enable with `--segmented`. See `config/examples/influenza_a.yaml`.

---

## Config

`config/default_config.yaml` documents every option. Highlights:

```yaml
qc:
  remove_duplicates: true
  length_filter:
    mode: median_percent     # or "min_max"
    min_percent: 50
  ambiguous_threshold: 0.05
  annotation_filter:
    enabled: true
    keywords: [MAG, synthetic, partial, hypothetical, ...]

clustering:
  backend: mmseqs2
  mmseqs2_mode: easy-linclust   # or easy-cluster
  coverage: 0.8

representative:
  priority: [refseq, reviewed_uniprot, longest]
```

Environment variables `REPSEQ_NCBI_EMAIL` and `REPSEQ_NCBI_API_KEY` override `taxonomy.ncbi_email` / `taxonomy.ncbi_api_key`.

---

## Cache management

Taxonomy lookups are cached in a SQLite DB (default `~/.repseq/cache/taxonomy.db`):

```bash
repseq cache stats
repseq cache purge-expired
repseq cache clear --source ncbi
```

---

## Testing

```bash
pip install pytest
pytest tests/
```

Tests run offline — all network calls (NCBI, UniProt) are mocked.

---

## Status

`v0.1.0` — feature-complete draft of all 8 modes. Network-dependent paths
(NCBI Entrez, UniProt REST) have not been exercised end-to-end against live
APIs yet; offline regression tests cover everything else.
