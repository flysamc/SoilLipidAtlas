# Manuscript figures

Final rendered figures (submission picks) + the R code that generates them.
**REL/** below = the release folder `ncbi-phylum-2026-08-04-v1` in the
analysis archive. `../sync_results.py` re-copies rendered assets from the
archive by checksum.

**All final figures are generated in R.** Every data figure (main 2–5,
supplementary 1–8) renders through the shared house style
`soilmass_style.R` (Nature Comms sizing, Wong colour-blind-safe palette).
Each figure folder carries an `r/` bundle: the exact R render script, its
`soilmass_style.R` copy as run, and the plotted data CSVs — so any figure
re-renders with only R and that folder. The upstream *data-preparation*
steps (Python) live in the `methods/` pipeline and the analysis archive;
they compute the CSVs but generate no figures. The one exception is
Figure 1: a concept schematic generated as an editable SVG (the submitted
original was hand-drawn vector art, not a data plot).

## Main figures

| Fig | Folder | What it shows | R render | Plotted data | Caveat / open decision |
|---|---|---|---|---|---|
| 1 | `figure1/` | Concept-only workflow figure (4 panels, data-driven SVG; v2.1) | — (schematic exception: SVG generator `build_figure1_concept_v1.py`; legend in `figure1/LEGEND.md`) | reads release summaries live (six design-scale numbers only) | replaces the submitted Fig 1; no results shown by design (panel-d plot is schematic) |
| 2 | `figure2/` | Atlas biomarker counts (a) + public/soil validation (b, c) | `r/fig2_atlas.R` | `figure2/r/` CSVs (tier counts, soil validation; prepared by archive scripts `build_fig2_strict16_render.py`, `build_fig2bc_strict.py`) | wide-a layout locked; panel-c selection sets inverted vs submitted (IndVal 7,254 / composite 3,613) |
| 3 | `figure3/` | Dendrograms, lipidome–phylogeny Mantel, similarity heatmaps | `figure3_strict_composite_review_only.R` | `REL/figure3/main_figure3_review_only/` + curated SSU v3 | 16 phyla / 120 pairs; Mantel r 0.509/0.553 (locked to SSU v3); 4 multi-anchor units non-monophyletic |
| 4 | `figure4/` | SIMPER fingerprints (a), lipid-class heatmap (b), redundancy (c) | `r/fig5_fingerprint.R` (the recovered original R script) | strict SIMPER + `REL/annotation/lipid_classes.csv` (prepared by archive scripts `prepare/finalize_figure4_old_r_strict.py`) | panel-b denominator (full-atlas vs classified-only) still an open decision |
| 5 | `figure5/` | ClimGrass two-estimator composition (a) + treatment effects (b) | `r/fig5_final.R` (2026-08-17, the locked final design); earlier per-treatment variant `r/fig6_climgrass_v2.R` | `figure5/r/data/` + `../../methods/08_climgrass/results/` | Animalia ≈14.7% out-of-range flag; Rule A correction; locked Viridiplantae/Protists labels |

## Supplementary figures (`supplementary/`)

One R bundle per figure in `supplementary/r/suppfig{1..8}/` (render script +
`soilmass_style.R` + plotted data). Data preparation lives in the `methods/`
pipeline and the analysis archive.

| SF | What it shows | R render | Caveat |
|---|---|---|---|
| 1 | PCoA ordination + within-batch Mantel | `r/suppfig1/supp_fig1_submitted_layout.R` | Streptophyta tree placement pending author confirmation |
| 2 | Clade-conserved vs clade-exclusive features | `r/suppfig2/supp_fig2_clade.R` | uses the S15 five-rank framework |
| 3 | MS2LDA motif enrichment | `r/suppfig3/supp_fig3_ms2lda.R` | v3 render used (confirm vs v2) |
| 4 | Leave-one-out confusion matrix | `r/suppfig4/supp_fig4_loo.R` | 64.6% = declared reimplementation (submitted 52.8% unrecovered) |
| 5 | Negative control (pure-isolate decomposition) | `r/suppfig5/supp_fig5_negative_control.R` | declared reimplementation; data in `methods/07_decomposition/` |
| 6 | Negative-mode sensitivity decomposition | `r/suppfig6/supp_fig6_neg_sensitivity.R` | qualitative cross-check only |
| 7 | Annotation-source contributions + evidence tiers | `r/suppfig7/supp_fig7_annotation.R` | sampling-depth over 19 collection phyla |
| 8 | Cross-method fingerprint validation | `r/suppfig8/supp_fig8_cross_method.R` | SIMPER/SCBD exact; CAP/L1 reimplemented; producer in `methods/06_fingerprints/` |
