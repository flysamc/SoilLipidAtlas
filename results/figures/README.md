# Manuscript figures

Final rendered figures (submission picks) + the producer scripts that build
them. **REL/** below = the release folder `ncbi-phylum-2026-08-04-v1` in the
analysis archive. `../sync_results.py` re-copies rendered assets from the
archive by checksum.

## Main figures

| Fig | Folder | What it shows | Producer(s) here | Source data | Caveat / open decision |
|---|---|---|---|---|---|
| 1 | `figure1/` | Concept schematic (4 panels, data-driven SVG) | `build_figure1_concept_v1.py` | reads release CSVs + RUN_SUMMARYs live | final skeleton; submitted Fig 1 may be reused pending author confirmation |
| 2 | `figure2/` | Atlas biomarker counts (a) + public/soil validation (b, c) | `build_fig2_strict16_render.py`, `build_fig2bc_strict.py`, `fig2_atlas_v2_layout.R` | `REL/annotation/tier_counts.csv`; `REL/figure2_strict16_2026-08-11_v2_wide_a/` | wide-a layout locked; panel-c selection sets inverted vs submitted (IndVal 7,254 / composite 3,613) |
| 3 | `figure3/` | Dendrograms, lipidome–phylogeny Mantel, similarity heatmaps | `figure3_strict_composite_review_only.R` | `REL/figure3/main_figure3_review_only/` + curated SSU v3 | 16 phyla / 120 pairs; Mantel r 0.509/0.553 (locked to SSU v3); 4 multi-anchor units non-monophyletic |
| 4 | `figure4/` | SIMPER fingerprints (a), lipid-class heatmap (b), redundancy (c) | `prepare_figure4_old_r_strict.py`, `finalize_figure4_old_r_strict.py` (render via the recovered `fig5_fingerprint.R`, retained in the analysis archive) | strict SIMPER + `REL/annotation/lipid_classes.csv` | panel-b denominator (full-atlas vs classified-only) still an open decision |
| 5 | `figure5/` | ClimGrass corrected composition (a) + treatment effects (b) | `figure5_redesign_v2.py`, `render_figure5_final.py`, `fig6_climgrass_v2.R` + `soilmass_style.R` | `../../methods/08_climgrass/results/` | Animalia ≈14.7% out-of-range flag; Rule A correction |

## Supplementary figures (`supplementary/`)

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
