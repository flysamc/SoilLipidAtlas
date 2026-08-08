# Analysis progress

One row per methodology step, in the order of the paper's Methods section.
Details inside each folder's README.

| Step | What it does | Status | In repo? |
|---|---|---|---|
| 01 Taxonomy | Assigns each sample its NCBI phylum; 19 collected → 16 analysed | ✅ done, locked | ✅ |
| 02 Feature lists | Cross-batch alignment → one feature table per mode (273k POS / 123k NEG) | ✅ done | 🟡 pointers + schema (tables → Zenodo) |
| 03 Biomarker atlas | Selects phylum-enriched features (11,371 POS / 5,697 NEG) | ✅ done | ✅ |
| 04 Annotation | Identifies what the biomarker lipids are (11-step pipeline) | 🟡 steps 1–5 + 8 done; 6 partial; 7 + 9 waiting on missing resources | ✅ |
| 05 Public validation | Searches biomarkers in public data (fastMASST, Pan-ReDU) | ⏳ resumed 2026-08-08, ~6,000 queries left | 🟡 skeleton |
| 06 Fingerprints | SIMPER per-phylum fingerprints, MS2LDA motifs, cross-method validation | ✅ SIMPER done; others to review | ☐ |
| 07 Decomposition | Source decomposition framework + negative control | ☐ to review | ☐ |
| 08 ClimGrass | Applies fingerprints to field soil; quantification correction (IS/RIE) | 🚧 in progress — spectral matcher being reworked | ☐ |
| Figures | Built from finished steps above | ☐ not started | ☐ |
| Figure 1 | Being updated by hand (Rahul); editable SVG sources recovered 2026-08-07 (2 versions in `figures/figure1/`) | ⏸ | 🟡 sources only |
