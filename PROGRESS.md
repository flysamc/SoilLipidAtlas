# Methods index

One row per methodology step, in the order of the paper's Methods section.
Each folder's `README.md` documents the step in full.

| Step | What it does | Headline result | In this repository |
|---|---|---|---|
| 01 Taxonomy | Assigns each sample its NCBI phylum | 19 collected → **16 analysed** phyla | full (`methods/01_taxonomy/`) |
| 02 Feature lists | Cross-batch alignment → one feature table per mode | **273,248 POS / 122,571 NEG** features | scripts + schema; consensus tables on Zenodo |
| 03 Biomarker atlas | Selects phylum-enriched features | **11,371 POS / 5,697 NEG** biomarkers | full (`methods/03_biomarker_atlas/`) |
| 04 Annotation | Identifies the biomarker lipids (11-step pipeline) | 6,360 POS / 1,889 NEG annotated; SIRIUS coverage **70.4 % POS** | full (`methods/04_annotation/`) |
| 05 Public validation | Searches biomarkers in public MS/MS data | **2,287 POS / 486 NEG** soil-detected | summaries in git; per-feature match sets on Zenodo |
| 06 Fingerprints | SIMPER fingerprints, MS2LDA motifs, cross-method checks | 4 methods agree; **134/103** motif–phylum pairs | full; SIMPER atlas on Zenodo |
| 07 Decomposition | Source decomposition + pure-isolate negative control | **79.3 %** correct dominant group (n = 164) | full (`methods/07_decomposition/`) |
| 08 ClimGrass | Decodes real field soil; quantification correction | qSIP drought replication **q = 0.005** | full (`methods/08_climgrass/`) |
| Figures | Main 1–5 + Supplementary 1–8 | data figures rendered in R | `results/figures/` |
| Tables | Supplementary S1–S15 | live-formula workbooks | `results/tables/` |

## Notes on interpretation

- **Annotation scope.** Tiers come from LipidSearch, diagnostic ions,
  family propagation, RT screening, and ArchLips (steps 1–5 and 8).
  Retention-time sum-composition prediction does not change any tier.
  SIRIUS/CANOPUS/CSI report formula and class coverage; they do not
  promote or demote Gold/Silver/Bronze (`methods/04_annotation/README.md`).
- **Figure 4 panel b.** The lipid-class heatmap uses the **full atlas** as
  denominator. A classified-features-only reading is also defensible.
- **Leave-one-out accuracy** is **64.6 %** (106/164) on the 45,525-feature
  substrate (Supplementary Fig. 4, Table S13).
- **Table S15 / Supplementary Fig. 1.** Streptophyta's clade label and tree
  placement follow the five-rank framework in the Table S15 notes.
- **Figure 5 / ClimGrass.** Animalia ≈ 14.7 % lies above the literature
  range (reference-panel bias). Negative-mode archaeal annotation has no
  MS2 support and is not treated as a validated identification.
