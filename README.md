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
# Clone and install in editable mode (recommended for development)
git clone https://github.com/cmzmasek/repseq.git
cd repseq
pip install -e .

# Or install directly from GitHub
pip install git+https://github.com/cmzmasek/repseq.git

# With the optional visualization extras (matplotlib + umap-learn)
pip install -e '.[viz]'
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

All commands share these options: `--input/-i`, `--output-dir/-o`, `--config/-c`, `--threads`, `--seed`, `--segmented`, `--dry-run`, `--no-resolve`, `--source {auto,uniprot,ncbi,ncbi_virus}`, `--overflow {keep,trim}`, `--plot`.

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
QC: duplicates → length → ambiguous chars → annotation keyword → protein-count (optional)
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
output FASTA(s) + TSV metadata + plain-text run log
```

---

## Output files

Each run writes these files to `output.dir` (default `./repseq_output/`):

| File | Contents |
| --- | --- |
| `{prefix}_representatives.fasta` | Selected representative sequences |
| `{prefix}_representatives.tsv` | Metadata for each representative (accession, organism, host, country, date, taxonomy ranks, …) |
| `{prefix}_clusters.tsv` | Per-cluster summary (cluster ID, representative, size) |
| `{prefix}_qc_removed.tsv` | Sequences removed by QC and the reason |
| `{prefix}_run.log` | Plain-text run summary: parameters (YAML), QC stats, output file list |
| `{prefix}_isolate_proteins.tsv` | (Segmented + protein QC) One row per annotated CDS per passing isolate: `isolate_id, segment, accession, protein_id, product, length` |
| `{prefix}_proteins.fasta` | (Protein QC enabled) All protein amino-acid sequences from the selected representatives, in a single FASTA. Tagged headers: `>protein_id product [isolate=…] [segment=…] [parent=accession]` |
| `{prefix}_clustering.png` | (`--plot`, `[viz]` extras) Two-panel UMAP scatter — left colored by genus, right colored by cluster with `√(cluster size)` point scaling, member→rep lines, and an inset cluster-size histogram |

Segmented-virus runs additionally produce `{prefix}_concatenated.fasta` and one `{prefix}_segment_{name}.fasta` per segment.

---

## Segmented virus mode

For multi-segment viruses (e.g. influenza), `repseq` can:

1. Group sequences by isolate (regex on the header).
2. Identify each sequence's segment by its **canonical name, numeric index,
   or any user-defined synonym** (e.g. `large segment` → `L`, `hemagglutinin` → `HA`).
3. Keep only isolates that have **all** expected segments.
4. Concatenate per-isolate segments into a single sequence for clustering.
5. Write back both the concatenated FASTA and one FASTA per segment.

Configure under `segmented:` and enable with `--segmented`:

```yaml
segmented:
  enabled: false                  # or pass --segmented per-run
  virus: influenza_a
  viruses:
    influenza_a:
      expected_segments: 8
      segments: [PB2, PB1, PA, HA, NP, NA, M, NS]
      isolate_regex: "(?P<isolate>[AB]/[^/(\\s]+/[^/(\\s]+/[^/(\\s]+/\\d{4})"
      segment_aliases:            # optional: synonyms recognised in headers
        HA: [hemagglutinin]
        NA: [neuraminidase]
        NP: [nucleoprotein, "nucleocapsid protein"]
```

See `config/examples/influenza_a.yaml` for a fully annotated example.

---

## Clustering visualization (optional)

Pass `--plot` to render a two-panel UMAP scatter alongside the standard
outputs:

```bash
repseq taxonomic1 -c my.yaml -i seqs.fasta -r genus -n 5 --plot
```

- **Left panel** — every clustered sequence embedded with UMAP on a k-mer
  Jaccard distance, colored by genus (top 10 + "Other").
- **Right panel** — same coordinates, colored by cluster, with point size
  scaling with `√(cluster size)`. Faint lines link each member to its
  representative; representatives are outlined in black. An inset histogram
  shows the cluster-size distribution.

For large runs the embedding is subsampled (default cap 2000 points;
representatives always kept) and the member→rep lines are auto-suppressed
above 500 non-rep points to avoid spaghetti.

Requires `pip install 'repseq[viz]'`. Skipped silently for diversity-only
runs (`global -n`) since those produce no clusters to visualize.

---

## Protein-annotation QC (optional)

An optional pre-clustering step fetches GenBank CDS counts from NCBI and
drops sequences with insufficient or unexpected protein annotations:

```yaml
qc:
  protein_annotation:
    enabled: true       # off by default — opt in
    min_proteins: 1     # global floor: drop any sequence with fewer CDS features

segmented:
  viruses:
    influenza_a:
      # Per-segment protein-count filter (segmented mode only).
      # Each value is either an int (exact count required) or a list
      # of ints (any of the listed counts is acceptable — useful for
      # strain variation, e.g. nonfunctional PB1-F2 in 2009 H1N1).
      expected_proteins_per_segment:
        HA: 1
        M:  2          # M1 + M2
        PB1: [1, 2]    # PB1 alone, or PB1 + PB1-F2
        NS: [1, 2]     # NS1 alone, or NS1 + NEP
        # …
```

GenBank records are fetched in **batches of 200 accessions per request**
and cached in the same SQLite store as taxonomy lookups, so subsequent runs
on the same dataset are network-free. Skipped automatically with `--no-resolve`.

Translations from each CDS are captured alongside the metadata so the
optional `{prefix}_proteins.fasta` output requires no additional network
calls — everything is reconstructed from the same cached records.

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
    keywords: ["MAG:", synthetic, partial, hypothetical, ...]

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

Taxonomy and protein-annotation lookups are cached in a SQLite DB
(default `~/.repseq/cache/taxonomy.db`):

```bash
repseq cache stats
repseq cache purge-expired
repseq cache clear --source ncbi_taxonomy
repseq cache clear --source ncbi_nuccore
repseq cache clear --source ncbi_proteins    # batched GenBank/CDS records
repseq cache clear --source uniprot
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

`v0.3.0` — all 8 selection modes implemented, optional protein-annotation
QC with batched GenBank fetching, per-segment count checks (int or list-of-int),
segment-name synonyms, a protein FASTA writer reconstructed from cached
records, and an optional UMAP visualization of the clustering result
(`--plot`, behind the `[viz]` extras). **107 offline regression tests pass.**
The NCBI-backed paths have been exercised end-to-end against live Entrez
(influenza A H1N1 RefSeq genome, 8 segments + 11 proteins).
