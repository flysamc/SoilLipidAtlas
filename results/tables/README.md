# Supplementary tables S1–S15

Built workbooks (submission copies, 2026-08-12/13 builds) + their producers.
Each producer `build_table_s{N}_strict16.py` writes a `RUN_SUMMARY.json` and
CSV mirrors beside the workbook in the P2R release (`REL/tables/`); the xlsx
here are the submission-named copies. Live formulas (SUM/derived cells) are
kept live in the workbooks, not pasted values.

| Table | Content | Note |
|---|---|---|
| S1 | Per-organism / per-phylum sample inventory | 19 collection + 16-analysis flag; 2 submitted-S1 defects fixed |
| S2 | FBMN batch composition | reproduce-first PASS (338,878 features) |
| S3 | MS2 diagnostic-ion database | **no `build_table_s3` producer** — copied verbatim from the 50-rule authority in `../../methods/04_annotation/step03_diagnostic_ions/pos_50rule_authority/` |
| S4 | Mantel distance-metric sensitivity | 16 phyla / 120 pairs; seed bug fixed |
| S5 | Per-phylum biomarker counts by tier | ties to `tier_counts.csv`; frozen vs Step-9 guardrail |
| S6 | Expected composition ranges (literature) | Viridiplantae/Protists relabel |
| S7 | Empirical response (RIE) factors | Rule A: out-of-window → uncalibrated 1.0 |
| S8 | Per-treatment ClimGrass composition | Figure 5 v2; CLR permutation effects |
| S9 | SIRIUS/CANOPUS/CSI coverage by phylum | **2026-08-13 v2, gap-fill complete** — POS SIRIUS 8,004 (70.4%) |
| S10 | ClimGrass benchmark studies (literature) | |
| S11 | Cross-method dendrogram Mantel r | Supp Fig 8 panel a |
| S12 | Cross-method consensus per phylum | Discosea/Evosea re-derived, not relabelled |
| S13 | LOO accuracy by method × K | declared reimplementation; L1 best ≈74% |
| S14 | ClimGrass overlap fractions | 696-feature verified-soil substrate |
| S15 | Five-rank taxonomic framework | Streptophyta clade label pending coauthor confirmation |
