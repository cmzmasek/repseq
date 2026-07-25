# Hareavirales HMM set (bundled, order-specific)

Drop one or more HMMER3 profile files (`*.hmm`) into this directory — one
protein per file is fine, or a few concatenated profiles per file. repseq
**combines every `*.hmm` in this directory into one pressed database** at
run time (cached, rebuilt only when a file changes), so you don't have to
`cat` and `hmmpress` them yourself.

Select this set from a config with the bare order name:

```yaml
hmm:
  enabled: true
  database: Hareavirales     # → repseq/data/hmms/Hareavirales/  (this directory)
```

(`database` also accepts an absolute path to a single `.hmm` file or to any
directory of `.hmm` files; `null` uses the bundled `repseq_viral_core.hmm`.)

Curated profiles that carry a `GA` gathering-threshold line are gated by it
automatically (`hmm.use_ga_when_available`); profiles without one fall back to
`hmm.default_evalue` + `hmm.relative_length_cutoff`. This README is ignored by
the loader (only `*.hmm` files are read).

## What Hareavirales is

Hareavirales is one of the two orders the old *Bunyavirales* was split into when
that order became class **Bunyaviricetes** (the sister order is Elliovirales —
Peribunyaviridae, Hantaviridae, Tospoviridae, Phasmaviridae; see that set).
Its member families are segmented, negative-sense ssRNA viruses:

| Family | Segments | Hosts | RefSeq records |
|---|---|---|---|
| Phenuiviridae | 3 (L,M,S) mostly | vertebrates / arthropods / plants | ~473 |
| Nairoviridae | 3 (L,M,S) | ticks / vertebrates (CCHFV) | ~153 |
| Arenaviridae | **2 (L,S)** ambisense | rodents / snakes (Lassa, LCMV) | ~133 |
| Discoviridae, Leishbuviridae, Wupedeviridae, Mypoviridae, Tosoviridae | 2–3 | protists / fungi | ~27 total |

**Do not confuse it with Elliovirales.** Peribunyaviridae and Hantaviridae are
in *that* order. Phenuiviridae has its own family-level set
(`repseq/data/hmms/Phenuiviridae/`); this order set is the broader backbone.

## Two architectural facts that shape the config

**1. The RdRp is universal but split across two non-homologous Pfam families.**
Unlike Elliovirales (one `Bunya_RdRp` for the whole order), Hareavirales needs an
**OR gate**: the arenaviral L polymerase is a *distinct* Pfam
(`Arena_RNA_pol`, PF06317), not `Bunya_RdRp` (PF04196). Together they gate the
whole order:

| Family | RdRp profile | Organisms gated |
|---|---|---|
| Phenuiviridae | `Bunya_RdRp` | 150 / 154 |
| Nairoviridae | `Bunya_RdRp` | 45 / 49 |
| Arenaviridae | `Arena_RNA_pol` | 62 / 63 |
| small families | either | 6 / 6 |
| **order total** | `Bunya_RdRp OR Arena_RNA_pol` | **262 / 271 = 97 %** |

**2. Segment count differs: Arenaviridae is bi-segmented (L, S).** The
tri-segmented families (Nairoviridae, most Phenuiviridae) carry an M glycoprotein
segment; Arenaviridae carries only L and S (its S is ambisense, encoding both
nucleocapsid **and** glycoprotein). The **only two segments every family shares
are L and S**, so the bundled `repseq_config_hareavirales.yaml` declares
`segments: [L, S]` (`expected_segments: 2`), gates L on the RdRp OR-gate, and
leaves S ungated (longest-CDS = nucleocapsid). **The M glycoprotein of the
tri-segmented families is therefore excluded from the order backbone** — for full
tri-segmented genomes use the family configs.

### Two consequences worth stating in Methods

* **The order backbone is RdRp (L) + nucleocapsid (S) only.** Nairovirus /
  phenuivirus M-segment glycoprotein records are dropped as an undeclared
  segment. That is the price of spanning bi- and tri-segmented families in one
  definition.
* **Arenaviruses are non-homologous to the rest at *both* markers** — a distinct
  RdRp (`Arena_RNA_pol`) *and* a distinct nucleocapsid (`Arena_nucleocap`). So
  clustering correctly separates Arenaviridae from the bunyavirus-like families
  (identity-based, no cross-family alignment needed), but the **deep tree node
  joining arenaviruses to the rest is unreliable** — the concatenated L+S MSA
  does not align across that divide. Read the whole-genome tree as two
  well-resolved sub-backbones (arenaviral / bunyavirus-like), not one calibrated
  deep phylogeny. This is the same caveat class as the Amarillovirales order set
  (`Flavi_NS5` vs `RdRP_3`).

### `relative_length_cutoff`

`repseq_config_hareavirales.yaml` sets `hmm.relative_length_cutoff: 0.25`. The
RdRp models are long (`Bunya_RdRp` 741 aa core, `Arena_RNA_pol` 2048 aa
full-length) and cover the bulk of real records well (Nairoviridae min 0.67,
Phenuiviridae/Arenaviridae median 1.00) but have a fragment tail; 0.25 keeps
divergent full-length polymerases while dropping partial-L fragments. The GA
bit-score is the real identity gate.

## Bundled profiles

All from Pfam-A (CC0); **all 20 carry curated GA cutoffs**. Membership was taken
from InterPro's Pfam-entries-by-taxonomy query for *Hareavirales* (taxid
3151839). `Bunya_RdRp` is byte-identical to the `repseq_viral_core.hmm` copy
(bundled sets are self-contained; shared profiles duplicated deliberately).

| Pfam | Profile (`NAME`) | Segment | Role |
|------|------------------|---------|------|
| PF04196 | `Bunya_RdRp` | L | **RdRp — nairo/phenui/others (gate)** |
| PF06317 | `Arena_RNA_pol` | L | **RdRp — arenaviruses (gate)** |
| PF15518 | `L_protein_N` | L | L protein N-terminus |
| PF12603 | `L_PA-C-like` | L | L PA-C-like domain |
| PF17296 | `ArenaCapSnatch` | L | Arenavirus cap-snatching endonuclease |
| PF02338 | `OTU` | L | Nairovirus L OTU protease |
| PF05733 | `Tenui_N` | S | *Phenuiviridae* nucleocapsid |
| PF02477 | `Nairo_nucleo` | S | *Nairoviridae* nucleocapsid |
| PF00843 | `Arena_nucleocap` | S | *Arenaviridae* nucleocapsid N-terminal |
| PF17290 | `Arena_ncap_C` | S | *Arenaviridae* nucleocapsid C-terminal |
| PF11073 | `NSs` | S | Phlebovirus-like NSs |
| PF00798 | `Arena_glycoprot` | S | *Arenaviridae* GPC (ambisense on S) |
| PF07243 | `Phlebovirus_G1` | M | *Phenuiviridae* glycoprotein Gn |
| PF07245 | `Phlebovirus_G2` | M | *Phenuiviridae* glycoprotein Gc |
| PF19019 | `Phlebo_G2_C` | M | *Phenuiviridae* glycoprotein Gc C-term |
| PF07246 | `Phlebovirus_NSM` | M | *Phenuiviridae* NSm |
| PF20726 | `Nairovirus_Gn` | M | *Nairoviridae* glycoprotein Gn |
| PF07948 | `Nairovirus_GP38` | M | *Nairoviridae* GP38 |
| PF20727 | `Nairovirus_MLD` | M | *Nairoviridae* mucin-like domain |
| PF20728 | `Nairovirus_NSm` | M | *Nairoviridae* NSm |

The M-segment glycoprotein profiles don't gate anything in the order run (M is
undeclared), but they earn their place two ways: `Arena_glycoprot` annotates the
arenaviral GPC (which sits on the **S** segment), and all of them let a
family-specific config (which *does* declare M) reuse this set. In the order run,
`{prefix}_hmm_diagnostic.tsv` still annotates each L/S CDS with its
family-specific identity.

## Reassortment / full glycoprotein analysis lives in the family configs

Because M is excluded and S is ungated, the order run produces one whole-genome
(L+S concat) tree, not per-segment marker trees. Per-segment reassortment is
within-family by nature and belongs in the family configs
(`repseq_config_phenuiviridae.yaml`; a Nairoviridae or Arenaviridae config can
reuse this set with all its segments gated).
