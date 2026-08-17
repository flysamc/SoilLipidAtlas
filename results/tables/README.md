# Supplementary tables S1–S15

Workbooks (`Table_S{N}.xlsx`) and the scripts that build them.
Producers `build_table_s{N}_strict16.py` live in this folder, except Table S3
(see below). Sums and percentages in the workbooks are live formulas, not
pasted values.

| Table | Content | Notes |
|---|---|---|
| S1 | Per-organism / per-phylum sample inventory | 19 collection phyla with a 16-analysis flag. Reports genera; a species count is not derived from the metadata |
| S2 | FBMN batch composition | 338,878 features across the twelve GNPS2 batches |
| S3 | MS2 diagnostic-ion database | **50-rule** executable classifier. No `build_table_s3` script — the workbook is the export of `../../methods/04_annotation/step03_diagnostic_ions/pos_50rule_authority/` |
| S4 | Mantel distance-metric sensitivity | Bray–Curtis, Jaccard, and cosine lipid distances versus the five-rank NCBI matrix of Supplementary Fig. 1 (16 phyla / 120 pairs). Figure 3b uses the SSU tree instead (`methods/09_phylogeny/`) |
| S5 | Per-phylum biomarker counts by annotation tier | From `tier_counts.csv`. Tiers are not updated by SIRIUS/CANOPUS (Step 10) |
| S6 | Expected composition ranges (literature) | Organism-group labels: Viridiplantae, Protists |
| S7 | Empirical response (RIE) factors | Rule A: out-of-window calibration → uncalibrated factor 1.0 |
| S8 | Per-treatment ClimGrass composition | Source for Figure 5; CLR permutation tests |
| S9 | SIRIUS / CANOPUS / CSI coverage by phylum | POS SIRIUS 8,004 formulas (70.4 %) |
| S10 | ClimGrass benchmark studies (literature) | |
| S11 | Cross-method dendrogram Mantel r | Supplementary Fig. 8 panel a |
| S12 | Cross-method consensus per phylum | Discosea and Evosea kept as separate phyla (not collapsed) |
| S13 | Leave-one-out accuracy by method × K | L1 best ≈ 74 %. Supplementary Fig. 8 panel b |
| S14 | ClimGrass overlap fractions | 696-feature verified-soil substrate. Supplementary Fig. 8 panel c |
| S15 | Five-rank taxonomic framework | Used by Supplementary Figs. 1–2. Figure 3b uses SSU patristic distance (`methods/09_phylogeny/`). Streptophyta clade label is marked in the workbook notes |
