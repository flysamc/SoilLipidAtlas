# Manuscript figures

Rendered figures and the code that draws them.

**All data figures are generated in R.** Main Figures 2–5 and Supplementary
Figures 1–8 render through the shared house style `soilmass_style.R` (Nature
Communications sizing, Wong colour-blind-safe palette). Each figure folder
carries an `r/` bundle: the render script, its `soilmass_style.R` copy as
run, and the plotted data CSVs, so the figure re-renders with only R and
that folder. Upstream data preparation lives in `methods/`.

**Figure 1 is the exception:** a workflow schematic drawn as an editable
SVG, not a data plot.

## Main figures

| Fig | Folder | What it shows | How it is drawn | Notes |
|---|---|---|---|---|
| 1 | `figure1/` | Study workflow (four panels) | SVG generator `build_figure1_concept_v1.py`; legend in `figure1/LEGEND.md` | Concept figure: six design-scale numbers only. Panel d's composition plot is schematic (quantitative results are in Fig. 5) |
| 2 | `figure2/` | Atlas biomarker counts (a) + public/soil validation (b, c) | `r/fig2_atlas.R`; plotted CSVs in `r/data/`; legend in `figure2/LEGEND.md` | Panel c reports both selection sets: Indicator Value 7,254 features (30.9 % soil-detected) and composite 3,613 (14.6 %) |
| 3 | `figure3/` | Dendrograms, lipidome–phylogeny Mantel, similarity heatmaps | `figure3_strict_composite_review_only.R` (assembles the three R-rendered panels) | 16 phyla / 120 pairs; Mantel r = 0.509 (POS) / 0.553 (NEG) against SSU v3 evolutionary distance; four multi-anchor units are non-monophyletic |
| 4 | `figure4/` | SIMPER fingerprints (a), lipid-class heatmap (b), redundancy (c) | `r/fig5_fingerprint.R` (historical filename; this is Main Figure 4) | Panel b uses the **full atlas** as denominator. A classified-features-only reading is also defensible; both are documented with the figure |
| 5 | `figure5/` | ClimGrass composition (a) + treatment effects (b) | `r/fig5_final.R`; plotted CSVs in `r/data/` and `methods/08_climgrass/results/` | Two estimators (fc-weighted bars, marker-panel diamonds). Animalia ≈ 14.7 % lies above the literature range (reference-panel bias; flagged in the manuscript). Display labels: Viridiplantae, Protists |

## Supplementary figures (`supplementary/`)

One R bundle per figure in `supplementary/r/suppfig{1..8}/` (render script +
`soilmass_style.R` + plotted data).

| SF | What it shows | R render | Notes |
|---|---|---|---|
| 1 | PCoA ordination + within-batch Mantel | `r/suppfig1/supp_fig1_submitted_layout.R` | Streptophyta tree placement follows the five-rank framework (Table S15 notes) |
| 2 | Clade-conserved vs clade-exclusive features | `r/suppfig2/supp_fig2_clade.R` | Uses the Table S15 five-rank framework |
| 3 | MS2LDA motif enrichment | `r/suppfig3/supp_fig3_ms2lda.R` | Motif–phylum pairs from `methods/06_fingerprints/` |
| 4 | Leave-one-out confusion matrix | `r/suppfig4/supp_fig4_loo.R` | **64.6 %** (106/164); declared reimplementation (earlier 52.8 % producer not recovered) |
| 5 | Negative control (pure-isolate decomposition) | `r/suppfig5/supp_fig5_negative_control.R` | Declared reimplementation; data in `methods/07_decomposition/` |
| 6 | Negative-mode sensitivity decomposition | `r/suppfig6/supp_fig6_neg_sensitivity.R` | Qualitative cross-check of the positive-mode result |
| 7 | Annotation-source contributions + evidence tiers | `r/suppfig7/supp_fig7_annotation.R` | Sampling depth over the 19 collection phyla |
| 8 | Cross-method fingerprint validation | `r/suppfig8/supp_fig8_cross_method.R` | SIMPER and SCBD recovered exactly; CAP and L1 are declared reimplementations. Producer in `methods/06_fingerprints/` |
