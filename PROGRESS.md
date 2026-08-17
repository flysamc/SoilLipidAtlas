# Analysis overview

One row per methodology step, in the order of the paper's Methods section.
Each folder's `README.md` documents the step in full.

| Step | What it does | Headline result | Status | In this repo |
|---|---|---|---|---|
| 01 Taxonomy | Assigns each sample its NCBI phylum | 19 collected → **16 analysed** phyla | ✅ complete (locked) | ✅ full |
| 02 Feature lists | Cross-batch alignment → one feature table per mode | **273,248 POS / 122,571 NEG** features | ✅ complete | 🗂 scripts + schema (tables → Zenodo) |
| 03 Biomarker atlas | Selects phylum-enriched features | **11,371 POS / 5,697 NEG** biomarkers | ✅ complete | ✅ full |
| 04 Annotation | Identifies the biomarker lipids (11-step pipeline) | 6,360 POS / 1,889 NEG annotated; SIRIUS coverage **70.4 % POS** | ✅ complete — steps 6/7/9 are declared exceptions | ✅ full |
| 05 Public validation | Searches biomarkers in public MS/MS data | **2,287 POS / 486 NEG** soil-detected | ✅ complete | 🗂 summaries (match sets → Zenodo) |
| 06 Fingerprints | SIMPER fingerprints, MS2LDA motifs, cross-method checks | 4 methods agree; **134/103** motif–phylum pairs | ✅ complete | ✅ full (SIMPER atlas → Zenodo) |
| 07 Decomposition | Source decomposition + pure-isolate negative control | **79.3 %** correct dominant group (n = 164) | ✅ complete (declared reimplementation) | ✅ full |
| 08 ClimGrass | Decodes real field soil; quantification correction | qSIP drought replication **q = 0.005** | ✅ complete | ✅ full |
| Figures | Main 1–5 + Supplementary 1–8 | all R-rendered | ✅ built | ✅ `results/figures/` |
| Tables | Supplementary S1–S15 workbooks | all built | ✅ built | ✅ `results/tables/` |

## Declared limitations and open items

These are documented choices and known limitations, stated in full in the
relevant step README:

- **Figure 4 panel b** — the denominator (full atlas vs classified-only
  features) admits two defensible readings; both are documented in
  `results/figures/README.md`.
- **Supplementary Fig. 3** — the v3 render is used; v2 is retained in the
  analysis archive for comparison.
- **Table S15 / Supp. Fig. 1** — the clade label and tree placement for
  Streptophyta follow the framework documented in the table notes.
- **Annotation steps 6, 7 and 9** — historical producers were not recoverable;
  each is handled as a declared exception, never silently substituted
  (`methods/04_annotation/README.md`).
