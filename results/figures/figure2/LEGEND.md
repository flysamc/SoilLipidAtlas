# Figure 2 — legend of record

**Fig. 2 | A cross-kingdom lipid atlas identifies phylum-enriched biomarkers that are widely detected in public and soil datasets.**

**a** Number of phylum-enriched biomarkers per phylum (positive ionisation mode; *n* = 11,371 across 16 phyla), coloured by annotation confidence tier: Gold, molecular species (headgroup and acyl chains, e.g. PE(P-16:0_16:1)); Silver, sum composition (headgroup and total carbons:double bonds, e.g. PC 34:2); Bronze, lipid class only; Unidentified, no structural assignment. Numbers at bar ends are per-phylum totals; phylum names are coloured by organism group. The count axis is broken (0–1,150 and 4,800–5,300; both segments at identical scale) to accommodate Streptophyta (*n* = 5,168). Silver annotations (*n* = 119) occur only in Methanobacteriota and Thermoproteota.

**b** Percentage of biomarkers detected in public MS/MS repositories (fastMASST, cosine similarity ≥ 0.7; 10,867 of 11,371 biomarkers searchable), by organism group (rows; *n* = 678–5,333 biomarkers per group) and ReDU sample type (columns; the five organism-source categories are shown, remaining categories in Source Data). Percentages are lower bounds, as 52.4% of matched files lack ReDU sample-type annotation and count towards no category.

**c** Percentage of biomarkers detected in 149 soil and environmental datasets (Pan-ReDU), by selection method: Indicator Value (features significant in ≥2 acquisition batches; 2,244 of 7,254, 30.9%) versus composite scoring (528 of 3,613, 14.6%).

Source data are provided as a Source Data file.

---

## Provenance note (not for publication)

Written from the rendered figure as the authority (strict release
`ncbi-phylum-2026-08-04-v1`, render workspace `figure2_strict16_2026-08-11_v2_wide_a`
in the analysis archive). All numbers re-derived from the three CSVs plotted by
`r/fig2_atlas.R` (`tier_counts.csv`, `kingdom_sampletype_summary.csv`,
`shared_vs_exclusive_soil.csv`). Tier definitions follow the annotation-confidence
axis of the positive-mode methods record (Gold = molecular species; Silver = sum
composition; Bronze = class only).

The panel arrangement departs from the submitted figure on author instruction
(panel **a** full height with a broken count axis, **b** top right, **c** beneath);
palette, theme, type sizes and column width follow the submitted figure. All three
panels are corrected builds, not reproductions: the submitted Figure 2a producer is
only partially recoverable, the historical fastMASST sample-type categoriser is
unrecovered, and both panel-**c** selection sets changed size under the strict
release (Indicator Value 1,647 → 7,254; composite 8,962 → 3,613), so the submitted
45.8% / 11.2% are not recomputable. Every number in the submitted caption except
the 149 soil and environmental datasets is superseded.

Notes for the manuscript pass: the main text should use a single term for the
panel-**c** Indicator Value set (currently also called "cross-batch consensus
biomarkers" in the submitted legend); the per-figure Source Data file promised
above is still to be assembled and must include the per-group *n* values and the
five ReDU sample-type categories not shown in panel **b**.
