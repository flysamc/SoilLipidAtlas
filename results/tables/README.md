# Supplementary tables S1–S15

Built workbooks (submission copies) + their producers.
Each producer `build_table_s{N}_strict16.py` writes a `RUN_SUMMARY.json` and
CSV mirrors beside the workbook in the analysis archive's release folder
(`tables/` under `ncbi-phylum-2026-08-04-v1`); the xlsx here are the
submission-named copies. Derived cells (sums, percentages) are live formulas
in the workbooks, not pasted values.

| Table | Content | Note |
|---|---|---|
| S1 | Per-organism / per-phylum sample inventory | 19 collection + 16-analysis flag; corrects two inconsistencies in the submitted version |
| S2 | FBMN batch composition | reproduce-first PASS (338,878 features) |
| S3 | MS2 diagnostic-ion database | **no `build_table_s3` producer** — copied verbatim from the 50-rule authority in `../../methods/04_annotation/step03_diagnostic_ions/pos_50rule_authority/` |
| S4 | Mantel distance-metric sensitivity | 16 phyla / 120 pairs; permutation seed handling corrected |
| S5 | Per-phylum biomarker counts by tier | ties to `tier_counts.csv`; frozen vs Step-9 guardrail |
| S6 | Expected composition ranges (literature) | Viridiplantae/Protists relabel |
| S7 | Empirical response (RIE) factors | Rule A: out-of-window → uncalibrated 1.0 |
| S8 | Per-treatment ClimGrass composition | Figure 5 v2; CLR permutation effects |
| S9 | SIRIUS/CANOPUS/CSI coverage by phylum | ✅ gap-fill complete — POS SIRIUS 8,004 (70.4 %) |
| S10 | ClimGrass benchmark studies (literature) | |
| S11 | Cross-method dendrogram Mantel r | Supp Fig 8 panel a |
| S12 | Cross-method consensus per phylum | Discosea/Evosea re-derived, not relabelled |
| S13 | LOO accuracy by method × K | declared reimplementation; L1 best ≈74% |
| S14 | ClimGrass overlap fractions | 696-feature verified-soil substrate |
| S15 | Five-rank taxonomic framework | Streptophyta clade label marked provisional in the workbook notes |
