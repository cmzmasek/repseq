# repseq

**Pick a small, clean, representative set of sequences out of a big messy FASTA file.**

You downloaded 80,000 influenza sequences from NCBI. You need 200 good ones that
still cover the real diversity — for a tree, for a figure, for a training set, or
to seed a curated reference database. Doing that by hand is miserable. `repseq`
does it for you: it cleans up the data, looks up what each sequence actually is,
groups similar sequences together, and keeps one good example from each group.

It works on protein **or** nucleotide FASTA files, from UniProt, NCBI, or NCBI
Virus, and it has a strong focus on viral sequences (including segmented viruses
like influenza).

---

## What "representative selection" means here

If you have 5,000 nearly-identical H3N2 sequences and 3 unusual ones, a random
sample of 200 will be 200 H3N2 sequences and you'll lose the unusual ones.
`repseq` instead:

1. **Cleans** — drops duplicates, truncated junk, sequences full of `N`s, and
   records labelled "hypothetical", "synthetic", "partial", etc.
2. **Identifies** — looks up each sequence's organism, host, country, and
   collection date from NCBI/UniProt (results are cached, so it's only slow once).
3. **Groups** — buckets the sequences (by how similar they are, or by genus, host,
   year, country… your choice).
4. **Keeps the best example from each group** — preferring curated records
   (RefSeq, reviewed UniProt) over random ones, and longer over shorter.

The result is a FASTA file that's small enough to work with but still spans the
diversity that was in the original.

---

## Installation

You need **Python 3.10 or newer**. Then:

```bash
git clone https://github.com/cmzmasek/repseq.git
cd repseq
pip install -e .
```

That's it for most uses. Two optional pieces:

**A sequence-clustering program** — `repseq` calls one to group sequences
*by similarity*. You only need it for the similarity-based modes (see the
table below); if you group by genus/host/year and your groups are small, or
you just want N diverse sequences, you can skip it. Two backends are
supported; pick **one**.

*MMseqs2* (default — fast, scales to very large datasets):

```bash
brew install mmseqs2                       # macOS
conda install -c bioconda mmseqs2          # Linux (conda)
```

*cd-hit* (classic alternative — set `clustering.backend: cdhit` in your
config). Slower than MMseqs2 on big inputs and has minimum-identity floors
(0.40 for protein, 0.80 for nucleotide), but produces tight, all-vs-all
clusters that many groups prefer for reference-set work:

```bash
brew install cd-hit                        # macOS
conda install -c bioconda cd-hit           # Linux (conda)
```

`repseq` finds whichever binary it needs (`mmseqs`, or `cd-hit` /
`cd-hit-est`) on your `PATH` automatically.

**Plots** — if you want the optional diagnostic scatter plot of the clustering
result:

```bash
pip install -e '.[viz]'
```

**MSA + tree (optional, `--phylo`)** — to build a multiple-sequence alignment
and an approximate-ML phylogeny over the final representatives, install
MAFFT and FastTree:

```bash
brew install mafft fasttree                     # macOS
conda install -c bioconda mafft fasttree        # Linux (conda)
```

If either is missing the rest of the run still finishes; `repseq` just skips
the tree step with a warning.

---

## Quickstart

**Step 1 — make a config file.** This is a small text file with your settings
(cleaning thresholds, your NCBI email, etc.). A wizard asks you the questions:

```bash
repseq init-config -o my_config.yaml
```

**Step 2 — run a mode.** Point it at your FASTA file and pick how to select:

```bash
# Keep one sequence per cluster of 90%-identical sequences
repseq global -c my_config.yaml -i seqs.fasta -T 0.90

# Keep up to 5 sequences per genus
repseq taxonomic1 -c my_config.yaml -i seqs.fasta -r genus -n 5

# Keep up to 10 sequences per host species
repseq host -c my_config.yaml -i seqs.fasta -n 10
```

When it finishes you'll see a one-line summary — how many sequences passed
cleaning and how many representatives were selected. If *nothing* came out, it
tells you the most likely reason (see [Troubleshooting](#troubleshooting)).

> **Tip:** add `--no-resolve` while you're experimenting. It skips the NCBI/UniProt
> lookups, so runs are fast — at the cost of not knowing each sequence's
> organism/host/country. Drop it for the real run.

---

## Choosing a mode

A "mode" is just *how you want the sequences grouped before one is picked from
each group*. All modes take `-i` (input), `-c` (config), and most take `-n`
(how many to keep per group).

| Command | What it does | Use it when… |
| --- | --- | --- |
| `global` | One big pass over everything. Either cluster at a similarity threshold `-T` (e.g. `-T 0.95` = group sequences ≥95% identical), or just ask for `-n` maximally-different sequences. | You want a flat, even sampling of the whole dataset. |
| `taxonomic1` | Up to `-n` sequences per taxonomic rank — `--rank genus`, `family`, `species`, etc. | You want even coverage across the tree of life (or of viruses). |
| `taxonomic2` | Like `taxonomic1` but **nested**: e.g. 20 per family, then 5 per genus within each. | One rank isn't enough — you want a hierarchy. |
| `host` | Up to `-n` sequences per host organism. | You care about which host the virus came from. |
| `time` | Up to `-n` per time window — `--window year`, `decade`, or a number like `5` for 5-year bins. | You want even coverage across collection dates. |
| `geographic` | Up to `-n` per country. | You want even geographic coverage. |
| `custom` | Up to `-n` per *any* field — a sequence attribute, a taxonomy rank, or a column in a metadata spreadsheet you provide. | Your grouping isn't one of the built-ins. |
| `hybrid` | Up to `-n` per *combination* of fields, e.g. `--fields genus,host,decade`. | You want a balanced grid across several variables at once. |

**Within each group**, if the group already has `-n` or fewer sequences, all of
them are kept. If it's bigger, `repseq` clusters it down to about `-n`
representatives. Add `--overflow trim` if you need *exactly* `-n` and not "about
`-n`".

Every mode also accepts: `--input/-i`, `--output-dir/-o`, `--config/-c`,
`--threads`, `--seed`, `--segmented`, `--dry-run`, `--no-resolve`,
`--source {auto,uniprot,ncbi,ncbi_virus}`, `--overflow {keep,trim}`, `--plot`.

---

## What happens during a run

```
your FASTA file(s)
        │
        ▼
  read it, figure out the format (UniProt / NCBI / NCBI Virus)
        │
        ▼
  look up organism, host, country, date     (from NCBI & UniProt, cached locally)
        │
        ▼
  clean: drop duplicates → too-short/long → too many ambiguous chars
         → bad-keyword annotations → (optional) wrong protein count
        │
        ▼
  (optional) segmented-virus step: keep only isolates that have ALL segments
        │
        ▼
  group + pick the best representative from each group
        │
        ▼
  write: representative FASTA + metadata tables + a plain-text run log
```

The run log (`{prefix}_run.log`) records exactly what settings were used and what
got dropped at each step — keep it with your results so the selection is
reproducible.

---

## Output files

Everything is written to the output directory (`./repseq_output/` by default):

| File | What's in it |
| --- | --- |
| `{prefix}_representatives.fasta` | **The main result** — your selected sequences. |
| `{prefix}_representatives.tsv` | A spreadsheet: one row per representative, with accession, organism, host, country, date, taxonomy. Opens in Excel. |
| `{prefix}_clusters.tsv` | Which sequences ended up grouped together, and which one was picked. |
| `{prefix}_group_counts.tsv` | One row per group (genus, host, year, country, … — whatever your mode stratified on): how many sequences went *in*, how many came *out*, whether clustering ran, and the similarity cutoff it settled on. The quickest way to see where the reduction happened. |
| `{prefix}_qc_removed.tsv` | Every sequence that was dropped during cleaning, and *why*. Check this if you lost more than expected. |
| `{prefix}_run.log` | Plain-text record of the settings used and the per-step counts. |
| `{prefix}_proteins.fasta` | *(if protein QC is on)* The protein sequences of all your representatives. |
| `{prefix}_isolate_proteins.tsv` | *(segmented + protein QC)* One row per gene per kept isolate. Columns: `protein_id`, `product`, `length` (aa), `isolate_id`, `segment`, `segment_length` (nt), `accession`, and the taxonomic ranks `species`, `subgenus`, `genus`, `subfamily`, `family`, `suborder`, `order`, `subclass`, `class` (sub-ranks come from the NCBI lineage and are often blank for viruses). |
| `{prefix}_clustering.png` | *(if `--plot`)* A diagnostic scatter plot of the clustering — see below. |
| `{prefix}_msa.fasta` | *(if `--phylo`)* MAFFT alignment of the representatives, FASTA headers are short ids (`S0001`…) for compatibility with downstream tools. |
| `{prefix}_tree.nwk` | *(if `--phylo`)* FastTree Newick — leaf names are the same short ids as in the MSA. |
| `{prefix}_tree.xml` | *(if `--phylo`)* **The tree you'll usually open** — phyloXML with the original sequence names restored on every leaf. |
| `{prefix}_tree_id_map.tsv` | *(if `--phylo`)* Two columns, `short_id` ↔ `original_id`, for decoding the MSA / Newick. |

Segmented-virus runs also write `{prefix}_concatenated.fasta` (all segments of an
isolate joined head-to-tail) and one `{prefix}_segment_<name>.fasta` per segment.

---

## Cleaning (QC) — what gets dropped and how to control it

All cleaning settings live under `qc:` in your config file. Defaults are sensible;
loosen them if you're losing sequences you want to keep.

```yaml
qc:
  remove_duplicates: true        # drop byte-identical sequences (keeps the curated copy)

  length_filter:
    mode: median_percent         # judge length relative to the dataset's median…
    min_percent: 50              #   …drop anything shorter than 50% of the median
    # max_percent: 200           #   …optionally also drop anything over 200%
    # ── or ──
    # mode: min_max              # judge length against fixed numbers instead
    # min_length: 1000
    # max_length: 20000

  ambiguous_threshold: 0.05      # drop sequences that are >5% N / X / other ambiguous letters

  annotation_filter:
    enabled: true
    keywords: ["MAG:", synthetic, partial, hypothetical, fragment, uncultured, ...]
    # any sequence whose description contains one of these words is dropped
```

A few things worth knowing:

- **`median_percent` compares every sequence to the median length of the whole
  file.** That's perfect for a single gene, but **wrong for a mixed file** (e.g.
  several different genes, or a whole genome plus its individual genes) — the
  median is meaningless and you'll drop things unfairly. For mixed files, use
  `min_max` with explicit numbers instead.
- **In segmented-virus mode, the whole-file length filter is skipped
  automatically** — a file of influenza segments mixes 2,300-nt and 890-nt
  sequences, so a single median can't work. Use per-segment length bounds instead
  (see below).

---

## Segmented viruses (influenza, etc.)

Segmented viruses store their genome in several separate pieces. NCBI gives you
one FASTA record per segment, so a single isolate is spread across multiple
records. `repseq` can stitch them back together:

1. **Group records by isolate** — using a pattern (regex) that matches the strain
   name in the header.
2. **Identify each record's segment** — by its name, its number, or a synonym you
   define (e.g. `hemagglutinin` → `HA`).
3. **Keep only complete isolates** — an isolate missing any expected segment is
   dropped.
4. **(Optional) length-check each segment** — drop an isolate if, say, its HA
   segment is suspiciously short.
5. **Concatenate** the segments of each complete isolate into one sequence, so the
   normal grouping/selection can run on whole isolates.

Configure it under `segmented:` and turn it on with `--segmented` on the command
line (or `enabled: true` in the config):

```yaml
segmented:
  enabled: false
  virus: influenza_a              # which entry below to use
  viruses:
    influenza_a:
      expected_segments: 8
      segments: [PB2, PB1, PA, HA, NP, NA, M, NS]   # canonical order
      isolate_regex: "(?P<isolate>[AB]/[^/(\\s]+/[^/(\\s]+/[^/(\\s]+/\\d{4})"
      segment_aliases:            # optional: words in headers that mean a segment
        HA: [hemagglutinin]
        NA: [neuraminidase]
        NP: [nucleoprotein, "nucleocapsid protein"]
      segment_lengths:            # optional: drop an isolate if a segment is out of range
        HA: {min: 1600, max: 1800}
        NS: {min: 800,  max: 1000}
```

`config/examples/influenza_a.yaml` is a complete, commented example you can copy.

> The `isolate_regex` is the part people get wrong most often. It has to match the
> strain identifier as it appears in *your* headers, and it must capture it either
> as a group named `isolate` or as the first parenthesised group. If no isolates
> come through, this is the first thing to check.

---

## Optional: protein-annotation QC

This step asks NCBI how many protein-coding genes (CDS features) each record has,
and drops records with too few — or, for segmented viruses, the wrong number per
segment. It's off by default; turn it on in the config:

```yaml
qc:
  protein_annotation:
    enabled: true
    min_proteins: 1     # drop any record with fewer annotated proteins than this

segmented:
  viruses:
    influenza_a:
      expected_proteins_per_segment:
        HA: 1
        M:  2           # M1 + M2
        PB1: [1, 2]     # PB1 alone, or PB1 + PB1-F2 — a list means "any of these"
        NS: [1, 2]      # NS1 alone, or NS1 + NEP
```

Records are fetched from NCBI in batches and cached locally, so a second run on
the same data needs no network. Skipped automatically under `--no-resolve`.

---

## Optional: clustering plot

Pass `--plot` (and install with `'.[viz]'`) to get a two-panel scatter plot,
`{prefix}_clustering.png`, that lets you eyeball whether the clustering looks
sensible:

- **Left** — every sequence as a dot, positioned so similar sequences sit close
  together, coloured by genus.
- **Right** — the same dots, coloured by cluster, with bigger dots for bigger
  clusters and faint lines from each sequence to its chosen representative.

For big datasets the plot is drawn from a subsample (the representatives are
always included). It's skipped for `global -n` runs, which produce no clusters.

---

## The config file

`config/default_config.yaml` is fully commented and documents every option — read
that file as the reference. The `repseq init-config` wizard writes a starter
config for you. The most-changed settings:

```yaml
qc:
  remove_duplicates: true
  length_filter:
    mode: median_percent
    min_percent: 50
  ambiguous_threshold: 0.05

taxonomy:
  ncbi_email: you@institute.org   # NCBI asks for this; without it you'll be rate-limited
  ncbi_api_key: null              # optional — get one from NCBI for faster lookups

clustering:
  backend: mmseqs2                # or "cdhit"
  mmseqs2_mode: easy-linclust     # fast; use easy-cluster for tighter, slower clustering
  coverage: 0.8
  # cd-hit options (only used when backend == cdhit) live under
  # `clustering.cdhit:` — see default_config.yaml for the full block.

representative:
  priority: [refseq, reviewed_uniprot, longest]   # tie-break order for picking the "best"
```

You can also set your NCBI email/key via the environment variables
`REPSEQ_NCBI_EMAIL` and `REPSEQ_NCBI_API_KEY` instead of putting them in the file.

---

## The local cache

Every NCBI/UniProt lookup is saved to a small database (`~/.repseq/cache/` by
default) so you only pay the network cost once. Manage it with:

```bash
repseq cache stats                           # how big is it, what's in it
repseq cache purge-expired                   # remove stale entries
repseq cache clear                           # wipe everything
repseq cache clear --source ncbi_taxonomy    # wipe just one kind of lookup
```

---

## Troubleshooting

**"WARNING: no representative sequences were selected."**
The run finished but nothing came out. `repseq` prints the most likely cause; the
usual ones are:

- *No sequences were loaded* — the input path is wrong, the file is empty, or its
  header format wasn't recognised. Try `--source ncbi_virus` (or `ncbi` /
  `uniprot`) to force it.
- *QC removed everything* — your cleaning thresholds are too strict for this data.
  Look at `{prefix}_qc_removed.tsv` to see which step did it, then loosen that
  setting. A common one: `median_percent` length filtering on a mixed-gene file —
  switch to `min_max`.
- *The segmented step dropped everything* — no isolate had all its segments. Most
  often the `isolate_regex` doesn't match your headers; also check the segment
  names/aliases and any `segment_lengths` bounds.

**`MMseqs2Error` / "mmseqs not found"** (or **`CDHitError` / "cd-hit not
found"** if you've set `clustering.backend: cdhit`) — the similarity-clustering
program isn't installed or isn't on your `PATH`. Install it (see
[Installation](#installation)), or use a mode that doesn't need it (`global -n`,
or a stratified mode where every group is already small).

**`cd-hit identity threshold X is below the supported floor`** — cd-hit refuses
identities below 0.40 (protein) or 0.80 (nucleotide). Either raise your
threshold to the floor, or switch the backend to `mmseqs2`, which has no
identity floor.

**Everything is grouped under "Unknown"** in a taxonomic/host/geographic run — the
metadata lookups didn't run or didn't find anything. Don't use `--no-resolve` for
the real run, and make sure your `ncbi_email` is set in the config.

**Lookups are slow the first time** — that's expected; they're cached, so the
*second* run on the same data is fast. An NCBI API key speeds up the first run.

**`[phylo skipped]` / `[phylo failed]`** — the `--phylo` step is fail-soft: if
fewer than 3 representatives survived, or `mafft` / `FastTree` are missing or
errored, the message is printed to stderr and the rest of the run's outputs
are still written. To enable it, install MAFFT and FastTree (see
[Installation](#installation)).

---

## Testing

```bash
pip install pytest
pytest tests/
```

The tests run fully offline — all network calls are simulated — so they're safe
to run anywhere and finish in a couple of seconds.

---

## Status

**`v0.5.4`** — all 8 selection modes, optional protein-annotation QC (with
per-segment counts and per-segment length bounds), segment-name synonyms, a
protein-FASTA output, and an optional UMAP plot of the clustering.

New in `v0.5.4`:

- **Won't overwrite a previous run.** If the output directory already exists
  and is not empty, the program now stops immediately with a clear error
  instead of writing new files alongside (or on top of) the old ones. Empty
  it, delete it, or point `--output-dir` somewhere else.

New in `v0.5.3`:

- **Correct duplicate removal for segmented viruses.** Exact-duplicate
  removal used to run on the pool of individual segments *before*
  concatenation. A segment that happens to be identical between two
  otherwise distinct isolates (a conserved segment) would get one copy
  dropped — and the affected isolate was then silently discarded as
  "incomplete". Duplicate removal now runs on the *concatenated isolates*
  instead: two isolates collapse only when every segment matches. This
  only affects segmented-mode runs; you may now see slightly more isolates
  retained.

New in `v0.5.2`:

- **Taxonomy lineage fix.** Genus/family/order are now resolved from NCBI's
  `efetch` endpoint, which actually returns the ranked lineage. The previous
  code read NCBI's taxonomy *summary* endpoint, which carries no lineage for
  viruses — so taxonomic modes silently grouped every viral sequence under
  "Unknown". **If you ran an earlier version, clear the cached taxonomy first:**
  `repseq cache clear --source ncbi_taxonomy` and
  `repseq cache clear --source ncbi_nuccore`.
- `repseq --version` always reports the real version (single source of truth
  in the package, no stale-install surprises).

New in `v0.5.1`:

- **Per-group counts report.** Every run now writes `{prefix}_group_counts.tsv`
  — one row per group (genus, host, year, country, …) with how many sequences
  went in, how many came out, whether clustering ran, and the similarity cutoff
  used. The fastest way to see exactly where the reduction happened.
- **Honest plot-dependency errors.** When `--plot` is skipped, the message now
  distinguishes "the plotting extras aren't installed" from "they're installed
  but failing to import" (usually a NumPy/SciPy version clash in the
  environment) instead of always telling you to reinstall.

New in `v0.5.0`:

- **Clearer endings.** Every run finishes with a one-line summary (how many
  sequences passed cleaning, how many representatives were selected) — or, if
  nothing came out, a warning naming the most likely cause.
- **Smarter segmented-virus cleaning.** The whole-file length filter is now
  skipped automatically in segmented mode, where a single median length is
  meaningless and would wrongly discard the short segments. Use per-segment
  `segment_lengths` instead.
- A full pipeline audit — corrected sequence-ID handling through clustering,
  RefSeq accession routing, the `MAG:` keyword filter, the similarity-threshold
  search direction, NCBI host/country/date harvesting, a length-robust diversity
  metric, and thread-safe caching.

**160 offline regression tests pass.** The NCBI-backed paths have been tested
end-to-end against a live influenza A H1N1 RefSeq genome (8 segments, 11 proteins).
