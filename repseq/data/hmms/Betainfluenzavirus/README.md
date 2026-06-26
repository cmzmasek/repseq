# Betainfluenzavirus HMM set (bundled)

Influenza **B** virus (genus *Betainfluenzavirus*, family Orthomyxoviridae) —
the B-tuned companion to the *Alphainfluenzavirus*-tuned `Orthomyxoviridae` set.
Both are influenza, but B's NS1, M2 (BM2) and the segment-6 NB glycoprotein are
non-homologous to their influenza A counterparts and need B-specific profiles;
the polymerase / HA / NA / NP / M1 / NEP families are pan-influenza and shared.

Drop one or more HMMER3 profile files (`*.hmm`) into this directory — one
protein per file is fine, or a few concatenated profiles per file. repseq
**combines every `*.hmm` in this directory into one pressed database** at
run time (cached, rebuilt only when a file changes), so you don't have to
`cat` and `hmmpress` them yourself.

Select this set from a config with the bare set name:

```yaml
hmm:
  enabled: true
  database: Betainfluenzavirus   # → repseq/data/hmms/Betainfluenzavirus/  (this directory)
```

(`database` also accepts an absolute path to a single `.hmm` file or to any
directory of `.hmm` files; `null` uses the bundled `repseq_viral_core.hmm`.)

Curated profiles that carry a `GA` gathering-threshold line are gated by it
automatically (`hmm.use_ga_when_available`); profiles without one fall back to
`hmm.default_evalue` + `hmm.relative_length_cutoff`. This README is ignored by
the loader (only `*.hmm` files are read).

## Scope: Influenza B virus (*Betainfluenzavirus*)

Influenza B has 8 negative-sense RNA segments. **Note the segment numbering
differs from influenza A**: in B, segment 1 is PB1 and segment 2 is PB2 (the
reverse of A), verified against the B/Lee/1940 reference (NC_002204 = PB1,
NC_002205 = PB2). Every one of the 18 profiles below was confirmed to hit the
B/Lee/1940 reference proteome at its curated Pfam GA cutoff.

## Bundled profiles

All from Pfam-A (CC0), carrying curated GA cutoffs. Listed by segment in
canonical (B) order. Most are the pan-influenza families also used by the
`Orthomyxoviridae` set; the three B-specific profiles are marked **(B)**.

| Segment | Pfam | Profile (`NAME`) | Influenza B protein / domain |
|---------|------|------------------|-------------------------------|
| 1 (PB1) | PF00602 | `Flu_PB1`       | PB1 — RNA-dependent RNA polymerase subunit |
| 2 (PB2) | PF20947 | `Flu_PB2_1st`   | PB2 — N-terminal region |
| 2 (PB2) | PF20948 | `Flu_PB2_2nd`   | PB2 — second domain |
| 2 (PB2) | PF20949 | `Flu_PB2_3rd`   | PB2 — middle domain |
| 2 (PB2) | PF20950 | `Flu_PB2_4th`   | PB2 — helical domain |
| 2 (PB2) | PF00604 | `Flu_PB2_5th`   | PB2 — cap-binding domain |
| 2 (PB2) | PF20951 | `Flu_PB2_6th`   | PB2 — 6th domain |
| 2 (PB2) | PF20952 | `Flu_PB2_7th`   | PB2 — C-terminal domain |
| 3 (PA)  | PF00603 | `Flu_PA`        | PA — RdRp subunit (endonuclease + C-terminal) |
| 4 (HA)  | PF00509 | `Hemagglutinin` | HA — haemagglutinin (HA1 + HA2) |
| 5 (NP)  | PF00506 | `Flu_NP`        | NP — nucleoprotein |
| 6 (NA)  | PF00064 | `Neur`          | NA — neuraminidase (sialidase) |
| 6 (NB)  | PF04159 | `NB` **(B)**    | NB — segment-6 glycoprotein (+1 ORF) |
| 7 (M)   | PF00598 | `Flu_M1`        | M1 — matrix protein |
| 7 (M)   | PF08289 | `Flu_M1_C`      | M1 — C-terminal domain |
| 7 (M)   | PF04772 | `Flu_B_M2` **(B)** | BM2 — ion-channel protein (segment-7 +1 ORF) |
| 8 (NS)  | PF02942 | `Flu_B_NS1` **(B)** | NS1 — non-structural protein (B-specific; A's `Flu_NS1` does NOT match B) |
| 8 (NS)  | PF00601 | `Flu_NS2`       | NS2 / NEP — nuclear export protein (spliced ORF) |

The 15 pan-influenza profiles are byte-identical copies of the same Pfam-A
families bundled in `Orthomyxoviridae`; the three **(B)** profiles were fetched
from InterPro by accession (`https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/<ACC>?annotation=hmm`).
The A-only `Flu_M2` (A M2), `PB1-F2`, and `Flu_NS1` (A NS1) are deliberately
**not** here — they don't apply to influenza B.

Reference a profile in config by its `NAME`, e.g. `hmms: ["Flu_NP"]`, an
OR-list of alternatives `hmms: ["Flu_M1", "Flu_B_M2"]`, or a multidomain token
like `"Flu_PB2_5th--Flu_PB2_6th--Flu_PB2_7th"` (consecutive PB2 domains, N→C).
See `repseq_config_betainfluenzavirus.yaml` for the per-segment `segment_markers`
wiring.
