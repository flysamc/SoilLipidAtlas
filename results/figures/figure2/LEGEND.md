# Figure 2 — legend of record

**Fig. 2 | A cross-kingdom lipid atlas identifies phylum-enriched biomarkers that are widely detected in public and soil datasets.**

**a** Number of phylum-enriched biomarkers per phylum (positive ionisation mode; biomarkers are the union of composite scoring and indicator value selections; *n* = 11,371 across 16 phyla), coloured by annotation confidence tier: Gold, molecular species (headgroup and acyl chains, e.g. PE(P-16:0_16:1)); Silver, sum composition (headgroup and total carbons:double bonds, e.g. PC 34:2); Bronze, lipid class only; Unidentified, no structural assignment. Numbers at bar ends are per-phylum totals; phylum names are coloured by organism group. The count axis is broken (0–1,150 and 4,800–5,300; both segments at identical scale) to accommodate Streptophyta (*n* = 5,168), whose Bronze–Unidentified tier boundary falls within the omitted range. Silver annotations (*n* = 119) occur only in Methanobacteriota and Thermoproteota.

**b** Percentage of biomarkers detected in public MS/MS repositories (fastMASST, cosine similarity ≥ 0.7; 10,867 of 11,371 biomarkers searchable), by organism group (rows; *n* = 678–5,333 biomarkers per group) and public sample type (columns; the five organism-source categories are shown, the remaining categories in Source Data). Percentages are lower bounds because 52.4% of matched files lack sample-type annotation and count towards no category.

**c** Percentage of biomarkers detected in 149 soil and environmental datasets (Pan-ReDU), by selection method: indicator value (IndVal; features significant in ≥2 acquisition batches; 2,244 of 7,254 features, 30.9%) versus composite scoring (528 of 3,613, 14.6%). Cross-batch indicator features, by definition shared across acquisition batches and typically across organisms, transfer to soil at roughly twice the rate of composite-scored features.

Source data are provided as a Source Data file.

## Plotted data in this repository

The legend numbers are those drawn by `r/fig2_atlas.R` from:

| File | Panel |
|---|---|
| `r/data/tier_counts.csv` | a |
| `r/data/kingdom_sampletype_summary.csv` | b |
| `r/data/shared_vs_exclusive_soil.csv` | c |
