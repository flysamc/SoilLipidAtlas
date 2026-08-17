# Analysis progress

One row per methodology step, in the order of the paper's Methods section.
Details inside each folder's README.

| Step | What it does | Status | In repo? |
|---|---|---|---|
| 01 Taxonomy | Assigns each sample its NCBI phylum; 19 collected → 16 analysed | ✅ done, locked | ✅ |
| 02 Feature lists | Cross-batch alignment → one feature table per mode (273k POS / 123k NEG) | ✅ done | 🟡 pointers + schema (tables → Zenodo) |
| 03 Biomarker atlas | Selects phylum-enriched features (11,371 POS / 5,697 NEG) | ✅ done | ✅ |
| 04 Annotation | Identifies what the biomarker lipids are (11-step pipeline) | ✅ done — SIRIUS gap-fill complete 2026-08-13 (POS 70.4% coverage); steps 6/7/9 = declared exceptions | ✅ |
| 05 Public validation | Searches biomarkers in public data (fastMASST, Pan-ReDU) | ✅ done 2026-08-11 (POS 2,287 / NEG 486 soil-detected) | ✅ summaries (match sets → Zenodo) |
| 06 Fingerprints | SIMPER per-phylum fingerprints, MS2LDA motifs, cross-method validation | ✅ done | ✅ (SIMPER atlas table → Zenodo) |
| 07 Decomposition | Two-estimator source decomposition + pure-isolate negative control (79.3% / n=164) | ✅ done (declared reimplementation) | ✅ |
| 08 ClimGrass | Fingerprints applied to field soil; Rule A correction; qSIP drought replication | ✅ done (Animalia ~14.7% flagged) | ✅ |
| Figures | Main 1–5 + Supplementary 1–8, submission picks + producers | ✅ built | ✅ `results/figures/` |
| Tables | Supplementary S1–S15 workbooks + producers | ✅ built | ✅ `results/tables/` |

## Open decisions (pending author confirmation)

- Figure 4 panel-b denominator: full-atlas vs classified-only.
- Supplementary Fig 3: v3 render used — confirm vs v2.
- Table S15 / Supp Fig 1: Streptophyta clade label + tree placement.
- Step 04 wording decisions listed in `methods/04_annotation/README.md`.

Details for each row live in that step's `README.md`.
