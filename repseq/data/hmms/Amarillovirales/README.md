# Amarillovirales HMM set (bundled)

This set covers the **order Amarillovirales** — the clade historically treated
as a single family *Flaviviridae*, which NCBI now splits into three families:

| Family | Representative genera |
|--------|-----------------------|
| Flaviviridae   | Orthoflavivirus (dengue, Zika, WNV, YFV, JEV, TBEV, …) |
| Hepaciviridae  | Hepacivirus (HCV), Pegivirus |
| Pestiviridae   | Pestivirus (BVDV, CSFV) |

The bundled profiles are Orthoflavivirus-tuned `Flavi_*` markers (which match
*only* family Flaviviridae) **plus** the pan-genus `RdRP_3` RdRp (which reaches
Hepaci-/Pegi-/Pestivirus). OR-ing `Flavi_NS5` with `RdRP_3` therefore gives a
marker that spans the whole order; using `Flavi_NS5` alone scopes a run to
family Flaviviridae (Orthoflavivirus).

Drop one or more HMMER3 profile files (`*.hmm`) into this directory — one
protein per file is fine, or a few concatenated profiles per file. repseq
**combines every `*.hmm` in this directory into one pressed database** at
run time (cached, rebuilt only when a file changes), so you don't have to
`cat` and `hmmpress` them yourself.

Select this set from a config with the bare set name:

```yaml
hmm:
  enabled: true
  database: Amarillovirales    # → repseq/data/hmms/Amarillovirales/  (this directory)
```

(`database` also accepts an absolute path to a single `.hmm` file or to any
directory of `.hmm` files; `null` uses the bundled `repseq_viral_core.hmm`.)

Curated profiles that carry a `GA` gathering-threshold line are gated by it
automatically (`hmm.use_ga_when_available`); profiles without one fall back to
`hmm.default_evalue` + `hmm.relative_length_cutoff`. This README is ignored by
the loader (only `*.hmm` files are read).

## Bundled profiles

All from Pfam-A (CC0), carrying curated GA cutoffs. The first block is the
**Orthoflavivirus** polyprotein, listed in N→C order (structural C–prM/M–E
first, then the non-structural NS1–NS5 cassette); these `Flavi_*` /
`Peptidase_S7` / `FtsJ` profiles are *Orthoflavivirus-specific* (family
Flaviviridae) and do **not** match the structurally distinct proteins of the
other two families (Hepaciviridae, Pestiviridae):

| Pfam | Profile (`NAME`) | Flaviviral protein / domain |
|------|------------------|------------------------------|
| PF01003 | `Flavi_capsid`      | Capsid protein C |
| PF01570 | `Flavi_propep`      | Polyprotein propeptide (pr of prM) |
| PF01004 | `Flavi_M`           | Membrane glycoprotein M |
| PF00869 | `Flavi_glycoprot`   | Envelope E — central + dimerisation domains |
| PF02832 | `Flavi_glycop_C`    | Envelope E — immunoglobulin-like domain |
| PF21659 | `Flavi_E_stem`      | Envelope E — stem/anchor domain |
| PF00948 | `Flavi_NS1`         | NS1 (non-structural) |
| PF01005 | `Flavi_NS2A`        | NS2A |
| PF01002 | `Flavi_NS2B`        | NS2B (protease cofactor) |
| PF00949 | `Peptidase_S7`      | NS3 — serine protease (peptidase S7) |
| PF07652 | `Flavi_DEAD`        | NS3 — DEAD-box helicase domain |
| PF20907 | `Flav_NS3-hel_C`    | NS3 — helicase C-terminal helical domain |
| PF01350 | `Flavi_NS4A`        | NS4A |
| PF01349 | `Flavi_NS4B`        | NS4B |
| PF01728 | `FtsJ`              | NS5 — FtsJ-like methyltransferase (N-terminal) |
| PF00972 | `Flavi_NS5`         | NS5 — RdRp fingers + palm domains |
| PF20483 | `Flavi_NS5_thumb`   | NS5 — RdRp thumb domain |

**Pan-genus RdRp** (added so the non-Orthoflavivirus families can be identified
at all — their RdRp is a different Pfam family, not `Flavi_NS5`):

| Pfam | Profile (`NAME`) | Flaviviral protein / domain |
|------|------------------|------------------------------|
| PF00998 | `RdRP_3`            | NS5B RdRp of **Hepacivirus / Pegivirus / Pestivirus** |

`RdRP_3` is the HCV NS5B family (representative structure 1gx6). Hepacivirus and
Pegivirus cover it near-fully; the more divergent Pestivirus NS5B clears the GA
bit-score by a wide margin but only matches its conserved core (~43 % of the
486-aa model), so admitting Pestivirus needs `hmm.relative_length_cutoff` ≤ 0.4.

## Non-Orthoflavivirus peptide profiles (Hepacivirus / Pegivirus / Pestivirus)

The `Flavi_*` profiles above identify the *Orthoflavivirus* polyprotein peptides
and match **only** that genus. The other genera have structurally distinct,
non-homologous mature peptides with their **own** Pfam families, bundled here so
each clade's polyprotein can be sliced into its real peptides (not just the
shared NS5/RdRp). See the three per-clade `polyprotein:` specs in
`repseq_config_amarillovirales.yaml`.

**Hepacivirus + Pegivirus** (`Core–E1–E2–NS2–NS3–NS4A–NS4B–NS5A–NS5B`; a
pegivirus lacks a distinct Core/E1 and uses one envelope, `GBV-C_env`):

| Pfam | Profile (`NAME`) | Peptide |
|------|------------------|---------|
| PF01542 | `HCV_core`         | Core |
| PF01543 | `HCV_capsid`       | Core (capsid domain) |
| PF01539 | `HCV_env`          | E1 |
| PF01560 | `HCV_NS1`          | E2 (HCV "E2/NS1") |
| PF12786 | `GBV-C_env`        | Pegivirus envelope (E2 slot) |
| PF01538 | `HCV_NS2`          | NS2 |
| PF02907 | `Peptidase_S29`    | NS3 protease |
| PF22027 | `NS3_helicase_C`   | NS3 helicase C-terminal |
| PF01006 | `HCV_NS4a`         | NS4A |
| PF01001 | `HCV_NS4b`         | NS4B |
| PF01506 | `HCV_NS5a`         | NS5A membrane anchor |
| PF08300 | `HCV_NS5a_1a`      | NS5A zinc-finger domain |
| PF08301 | `HCV_NS5a_1b`      | NS5A domain 1b |
| PF12941 | `HCV_NS5a_C`       | NS5A C-terminal |

**Pestivirus** (`Npro–C–Erns–E1–E2–NS2–NS3–NS4A–NS4B–NS5A–NS5B`; Erns, E1 and
NS4A/NS4B/NS5A have **no** Pfam family and are left unassigned, not mis-cut):

| Pfam | Profile (`NAME`) | Peptide |
|------|------------------|---------|
| PF05550 | `Peptidase_C53`    | Npro (leader autoprotease) |
| PF11889 | `Capsid_pestivir`  | Core (C) |
| PF16329 | `Pestivirus_E2`    | E2 |
| PF12387 | `Peptidase_C74`    | NS2 protease |
| PF05578 | `Peptidase_S31`    | NS3 protease |

NS3 helicase and NS5B for both groups reuse profiles listed above (`Flavi_DEAD`,
`Flav_NS3-hel_C`, `RdRP_3`). Coverage validated on the RefSeq set (0 N→C
ordering violations): Hepacivirus/Pegivirus reach 100 % on NS3/NS4B/NS5B and
45–80 % on the sparser structural peptides; Pestivirus reaches 82–100 % on all
six of its declared peptides. Before these profiles, the non-Orthoflavivirus
genera got **only** their NS5/RdRp sliced (via `RdRP_3`) and 0 % on everything
else.

## Scoping a run: order vs. family

- **Whole order Amarillovirales** (all genera) — OR the two RdRp profiles as the
  polyprotein marker: `hmms: ["Flavi_NS5", "RdRP_3"]` (Orthoflavivirus →
  `Flavi_NS5`, the other three genera → `RdRP_3`). See the worked example
  `repseq_config_amarillovirales.yaml`.
- **Family Flaviviridae only** (Orthoflavivirus) — use `hmms: ["Flavi_NS5"]`
  alone; the `RdRP_3`-only genera fail the gate and drop. See
  `repseq_config_flaviviridae.yaml`.

Reference a profile in config by its `NAME`, e.g. `hmms: ["Flavi_NS1"]` or a
multidomain token like
`"Peptidase_S7--Flavi_DEAD--Flav_NS3-hel_C"` (NS3 protease→helicase) or
`"FtsJ--Flavi_NS5--Flavi_NS5_thumb"` (NS5 MTase→RdRp).
