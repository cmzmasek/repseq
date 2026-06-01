# Filoviridae HMM set (bundled, family-specific)

Drop one or more HMMER3 profile files (`*.hmm`) into this directory — one
protein per file is fine, or a few concatenated profiles per file. repseq
**combines every `*.hmm` in this directory into one pressed database** at
run time (cached, rebuilt only when a file changes), so you don't have to
`cat` and `hmmpress` them yourself.

Select this set from a config with the bare family name:

```yaml
hmm:
  enabled: true
  database: Filoviridae      # → repseq/data/hmms/Filoviridae/  (this directory)
```

(`database` also accepts an absolute path to a single `.hmm` file or to any
directory of `.hmm` files; `null` uses the bundled `repseq_viral_core.hmm`.)

Curated profiles that carry a `GA` gathering-threshold line are gated by it
automatically (`hmm.use_ga_when_available`); profiles without one fall back to
`hmm.default_evalue` + `hmm.relative_length_cutoff`. This README is ignored by
the loader (only `*.hmm` files are read).

## Bundled profiles (v0.40.1)

All from Pfam-A (CC0), carrying curated GA cutoffs:

| Pfam | Profile (`NAME`) | Filovirus protein |
|------|------------------|-------------------|
| PF05505 | `Ebola_NP`             | Nucleoprotein (NP) |
| PF02097 | `Filo_VP35`            | VP35 (polymerase cofactor) |
| PF07447 | `Matrix_Filo`          | VP40 (matrix) |
| PF01611 | `Filo_glycop`          | GP (glycoprotein) |
| PF22307 | `Ebola-like_HR1-HR2`   | GP2 fusion HR1–HR2 |
| PF11507 | `Transcript_VP30`      | VP30 (transcription activator) |
| PF06389 | `Filo_VP24`            | VP24 |
| PF00946 | `Mononeg_RNA_pol`      | L (RdRp) |
| PF21080 | `Methyltrans_Mon_1st`  | L — methyltransferase |
| PF14314 | `Methyltrans_Mon_2nd`  | L — 2'-O-methyltransferase |
| PF21081 | `Methyltrans_Mon_3rd`  | L — methyltransferase |
| PF14318 | `Mononeg_mRNAcap`      | L — mRNA-capping domain |

Reference a profile in config by its `NAME`, e.g. `hmms: ["Filo_glycop"]` or
a multidomain token like `"Methyltrans_Mon_1st--Methyltrans_Mon_2nd--Methyltrans_Mon_3rd"`.
