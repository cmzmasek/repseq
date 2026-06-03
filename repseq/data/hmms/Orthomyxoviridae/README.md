# Orthomyxoviridae HMM set (bundled, family-specific)

Drop one or more HMMER3 profile files (`*.hmm`) into this directory — one
protein per file is fine, or a few concatenated profiles per file. repseq
**combines every `*.hmm` in this directory into one pressed database** at
run time (cached, rebuilt only when a file changes), so you don't have to
`cat` and `hmmpress` them yourself.

Select this set from a config with the bare family name:

```yaml
hmm:
  enabled: true
  database: Orthomyxoviridae   # → repseq/data/hmms/Orthomyxoviridae/  (this directory)
```

(`database` also accepts an absolute path to a single `.hmm` file or to any
directory of `.hmm` files; `null` uses the bundled `repseq_viral_core.hmm`.)

Curated profiles that carry a `GA` gathering-threshold line are gated by it
automatically (`hmm.use_ga_when_available`); profiles without one fall back to
`hmm.default_evalue` + `hmm.relative_length_cutoff`. This README is ignored by
the loader (only `*.hmm` files are read).

## Scope: Influenza A virus (*Alphainfluenzavirus*)

These profiles target the **influenza A** proteome. The polymerase (PB2/PB1/PA),
NP, M1/M2 and NS1/NS2(NEP) profiles are A-specific; `Hemagglutinin` (PF00509)
and `Neur` (PF00064) are the broader haemagglutinin / neuraminidase families
(canonical for influenza HA/NA, but not A-exclusive). Influenza **B** and **C**
have their own Pfam families (e.g. PF02942 NS1-B, PF03021/PF03026 M2-C/M1-C,
PF08720 HA-stalk-C) — those are **not** bundled here; this set is tuned for the
*Alphainfluenzavirus* genus.

## Bundled profiles

All from Pfam-A (CC0), carrying curated GA cutoffs. Influenza A has 8
negative-sense RNA segments; listed below by segment in canonical order. Pfam
splits the PB2 subunit into seven sequential domains (`Flu_PB2_1st`…`_7th`,
N→C); PB1 and PA are each a single large family.

| Segment | Pfam | Profile (`NAME`) | Influenza A protein / domain |
|---------|------|------------------|-------------------------------|
| 1 (PB2) | PF20947 | `Flu_PB2_1st`   | PB2 — N-terminal region |
| 1 (PB2) | PF20948 | `Flu_PB2_2nd`   | PB2 — second domain |
| 1 (PB2) | PF20949 | `Flu_PB2_3rd`   | PB2 — middle domain |
| 1 (PB2) | PF20950 | `Flu_PB2_4th`   | PB2 — helical domain |
| 1 (PB2) | PF00604 | `Flu_PB2_5th`   | PB2 — cap-binding domain |
| 1 (PB2) | PF20951 | `Flu_PB2_6th`   | PB2 — 6th domain |
| 1 (PB2) | PF20952 | `Flu_PB2_7th`   | PB2 — C-terminal domain |
| 2 (PB1) | PF00602 | `Flu_PB1`       | PB1 — RNA-dependent RNA polymerase subunit |
| 2 (PB1) | PF11986 | `PB1-F2`        | PB1-F2 proapoptotic accessory protein (+1 ORF) |
| 3 (PA)  | PF00603 | `Flu_PA`        | PA — RdRp subunit (endonuclease + C-terminal) |
| 4 (HA)  | PF00509 | `Hemagglutinin` | HA — haemagglutinin (HA1 + HA2) |
| 5 (NP)  | PF00506 | `Flu_NP`        | NP — nucleoprotein |
| 6 (NA)  | PF00064 | `Neur`          | NA — neuraminidase (sialidase) |
| 7 (M)   | PF00598 | `Flu_M1`        | M1 — matrix protein |
| 7 (M)   | PF08289 | `Flu_M1_C`      | M1 — C-terminal domain |
| 7 (M)   | PF00599 | `Flu_M2`        | M2 — ion-channel protein (spliced ORF) |
| 8 (NS)  | PF00600 | `Flu_NS1`       | NS1 — non-structural protein |
| 8 (NS)  | PF00601 | `Flu_NS2`       | NS2 / NEP — nuclear export protein (spliced ORF) |

Reference a profile in config by its `NAME`, e.g. `hmms: ["Flu_NP"]`, an
OR-list of alternatives `hmms: ["Flu_M1", "Flu_M2"]`, or a multidomain token
like `"Flu_PB2_5th--Flu_PB2_6th--Flu_PB2_7th"` (consecutive PB2 domains,
N→C). See `repseq_config_alphainfluenzavirus.yaml` for the per-segment
`segment_markers` wiring.
