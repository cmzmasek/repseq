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

Suggested Filoviridae marker proteins to profile here: **NP** (nucleoprotein),
**VP35**, **VP40** (matrix), **GP** (glycoprotein), **VP30**, **VP24**, and
**L** (RdRp). Curated profiles that carry a `GA` gathering-threshold line are
gated by it automatically (`hmm.use_ga_when_available`); profiles without one
fall back to `hmm.default_evalue` + `hmm.relative_length_cutoff`.

This file is a placeholder so the directory is tracked and packaged; it is
ignored by the loader (only `*.hmm` files are read).
