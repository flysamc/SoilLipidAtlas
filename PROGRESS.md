# Methods index

One row per methodology step, in the order of the paper's Methods section.
Each folder's `README.md` documents the step in full.

| Step | What it does | Headline result | In this repository |
|---|---|---|---|
| 01 Taxonomy | Assigns each sample its NCBI phylum | 19 collected → **16 analysed** phyla | full (`methods/01_taxonomy/`) |
| 02 Feature lists | Cross-batch alignment → one feature table per mode | **273,248 POS / 122,571 NEG** features | scripts + schema; consensus tables on Zenodo |
| 03 Biomarker atlas | Selects phylum-enriched features | **11,371 POS / 5,697 NEG** biomarkers | full (`methods/03_biomarker_atlas/`) |
| 04 Annotation | Identifies the biomarker lipids (11-step pipeline) | 6,360 POS / 1,889 NEG annotated; SIRIUS coverage **70.4 % POS** | full; steps 6, 7 and 9 are declared exceptions |
| 05 Public validation | Searches biomarkers in public MS/MS data | **2,287 POS / 486 NEG** soil-detected | summaries in git; per-feature match sets on Zenodo |
| 06 Fingerprints | SIMPER fingerprints, MS2LDA motifs, cross-method checks | 4 methods agree; **134/103** motif–phylum pairs | full; SIMPER atlas on Zenodo |
| 07 Decomposition | Source decomposition + pure-isolate negative control | **79.3 %** correct dominant group (n = 164) | full (declared reimplementation) |
| 08 ClimGrass | Decodes real field soil; quantification correction | qSIP drought replication **q = 0.005** | full (`methods/08_climgrass/`) |
| Figures | Main 1–5 + Supplementary 1–8 | data figures rendered in R | `results/figures/` |
| Tables | Supplementary S1–S15 | live-formula workbooks | `results/tables/` |

## Declared limitations

These are documented choices and known limits, stated in full in the
relevant README:

- **Annotation steps 6, 7 and 9.** Earlier producers were not recoverable.
  Each is a declared exception and is not silently substituted
  (`methods/04_annotation/README.md`).
- **Figure 4 panel b.** The lipid-class heatmap can be read over the full
  atlas or over classified features only; both are defensible. The figure
  shows the full atlas (`results/figures/README.md`).
- **Supplementary Fig. 4 / Table S13.** Leave-one-out accuracy is **64.6 %**
  (106/164). This is a declared reimplementation: the producer of the
  earlier 52.8 % figure was not recovered.
- **Table S15 / Supplementary Fig. 1.** Streptophyta's clade label and tree
  placement follow the five-rank framework in the Table S15 notes.
