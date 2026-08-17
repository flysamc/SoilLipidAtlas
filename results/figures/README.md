# Manuscript figures

Final rendered figures (submission picks) + the producer scripts that build
them. **REL/** below = the release folder `ncbi-phylum-2026-08-04-v1` in the
analysis archive. `../sync_results.py` re-copies rendered assets from the
archive by checksum.

**Render engine: every data figure (main 2–5, supplementary 1–8) is rendered
in R** through the shared house style `soilmass_style.R` (Nature Comms sizing,
Wong colour-blind-safe palette). Each figure folder carries an `r/`
subdirectory with the exact R render script, its `soilmass_style.R` copy as
run, and the plotted data CSVs — so any figure can be re-rendered with only R
and that folder. Python scripts alongside are the *data preparers* that
compute those CSVs from the release. The one exception is Figure 1: a concept
schematic generated as an SVG by a Python script (the submitted original was
hand-drawn vector art, not a data plot).

## Main figures

| Fig | Folder | What it shows | Producer(s) here | Source data | Caveat / open decision |
|---|---|---|---|---|---|
| 1 | `figure1/` | Concept schematic (4 panels, data-driven SVG) | `build_figure1_concept_v1.py` | reads release CSVs + RUN_SUMMARYs live | final skeleton; submitted Fig 1 may be reused pending author confirmation |
| 2 | `figure2/` | Atlas biomarker counts (a) + public/soil validation (b, c) | **R: `r/fig2_atlas.R`**; data prep `build_fig2_strict16_render.py`, `build_fig2bc_strict.py` | `figure2/r/` data CSVs (tier counts, soil validation) | wide-a layout locked; panel-c selection sets inverted vs submitted (IndVal 7,254 / composite 3,613) |
| 3 | `figure3/` | Dendrograms, lipidome–phylogeny Mantel, similarity heatmaps | `figure3_strict_composite_review_only.R` | `REL/figure3/main_figure3_review_only/` + curated SSU v3 | 16 phyla / 120 pairs; Mantel r 0.509/0.553 (locked to SSU v3); 4 multi-anchor units non-monophyletic |
| 4 | `figure4/` | SIMPER fingerprints (a), lipid-class heatmap (b), redundancy (c) | **R: `r/fig5_fingerprint.R`** (the recovered original R script); data prep `prepare_figure4_old_r_strict.py`, `finalize_figure4_old_r_strict.py` | strict SIMPER + `REL/annotation/lipid_classes.csv` | panel-b denominator (full-atlas vs classified-only) still an open decision |
| 5 | `figure5/` | ClimGrass two-estimator composition (a) + treatment effects (b) | **R: `r/fig5_final.R`** (2026-08-17, ports the locked final design); data prep `figure5_redesign_v2.py`; earlier per-treatment variant `r/fig6_climgrass_v2.R` | `figure5/r/data/` + `../../methods/08_climgrass/results/` | Animalia ≈14.7% out-of-range flag; Rule A correction; R render applies locked labels (Viridiplantae/Protists) |

## Supplementary figures (`supplementary/`)

R render bundles live in `supplementary/r/suppfig{1..8}/` (script +
`soilmass_style.R` + plotted data). The Python scripts below are the data
preparers.

| SF | What it shows | Producer here | Caveat |
|---|---|---|---|
| 1 | PCoA ordination + within-batch Mantel | `suppfig1_full_strict16.py`, `suppfig1_panel_d_permutations.py` | Streptophyta tree placement pending author confirmation |
| 2 | Clade-conserved vs clade-exclusive features | `suppfig2_clade_strict16.py` | uses the S15 five-rank framework |
| 3 | MS2LDA motif enrichment | `suppfig3_ms2lda_strict16.py` | v3 render used (confirm vs v2) |
| 4 | Leave-one-out confusion matrix | `suppfig4_loo_strict16.py` | 64.6% = declared reimplementation (submitted 52.8% unrecovered) |
| 5 | Negative control (pure-isolate decomposition) | in `../../methods/07_decomposition/` | declared reimplementation |
| 6 | Negative-mode sensitivity decomposition | `suppfig6_stage2_strict16.py` | qualitative cross-check only |
| 7 | Annotation-source contributions + evidence tiers | `suppfig7_evidence_tiers_release.py` | sampling-depth over 19 collection phyla |
| 8 | Cross-method fingerprint validation | in `../../methods/06_fingerprints/` | SIMPER/SCBD exact; CAP/L1 reimplemented |
