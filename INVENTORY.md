# Manuscript inventory — figures and tables

This is the master checklist for building this repository. One row = one unit of
review. A unit is added to the repository only after it has been read and
understood; its checkbox is then ticked and the row's folder gets its own README.

**Where this information comes from** (in the P2R working package):

- Script→figure mapping: `paper2_repro/portable.py` (the release build definition)
- Figure scripts + captions: `manuscript_2_clean/06_figures/figures_r/`
- Panel-level source data: `manuscript_2_clean/06_figures/source_folders/analysis18_09_figures/`
- Table builder: `manuscript_2_clean/07_tables/build_tables.py`
- As-submitted files: `submission_source/` (main article, supplementary, legends, tables workbook)
- Rendered outputs: `outputs/figures/`, `outputs/tables/`

**Numbering warning:** script filenames and caption files use an older figure
numbering than the submitted manuscript. The release build (`portable.py`) maps
old→new. During each figure's review, the final check is always against the
submitted legends file `SLA-figures-and-legends-nature-comms.docx`.

## Main figures

| ☐ | Unit | What it shows | Made by (script) | Data it reads | Status |
|---|---|---|---|---|---|
| ⏸ | **Figure 1** | Study design / overview (frozen image) | ⚠ **none — editable producer unknown**, only `Figure_1_FROZEN_SUBMISSION.png` exists | unknown | **Deferred — Rahul is updating this figure himself; added when ready** |
| ☐ | **Figure 2** | Cross-kingdom lipid atlas; phylum-enriched biomarkers; detection in public + soil data | `fig2_atlas.R` | `tier_counts.csv`, `kingdom_sampletype_summary.csv`, `shared_vs_exclusive_soil.csv` (+ `masst_by_kingdom.csv`, `tier2_by_kingdom.csv` in same folder) | Reproduced in Docker (frozen baseline); corrected 2A rerun exists but is provisional |
| ☐ | **Figure 3** | Chemotaxonomic dendrograms recapitulate phylogeny, 17 phyla, both modes | `fig4_dendrogram.R` (old numbering!) | `fig1a_dendrogram/` + `fig1b_mantel/` data (Bray-Curtis matrices, phylo distances, Mantel results) | Reproduced (frozen baseline); taxonomy-dependent → strict rerun expected |
| ☐ | **Figure 4** | Distributed lipid fingerprints (SIMPER curves, class composition, ensemble stability) | `fig5_fingerprint.R` (old numbering!) | `fig3a_simper_curves/`, `fig3b_lipid_class/`, `fig3d_ensemble_stability/`, phylum-kingdom map | Reproduced (frozen baseline); strict-release review copy exists (`figure4/` in release outputs) |
| ☐ | **Figure 5** | ClimGrass field application (soil composition, drought/climate effects) | `fig6_climgrass.R` (old numbering!) | `fig5a_climgrass_*` data dirs | Reproduced (frozen baseline); ⚠ corrected version blocked: spectral matcher must be fixed (cosine>1 bug), then strict 19-unit rerun |

## Supplementary figures

| ☐ | Unit | What it shows | Made by | Status |
|---|---|---|---|---|
| ☐ | **Supp Fig 1** | Within-batch PCoA recovers kingdom structure | `supp_fig1_within_batch.R` | Reproduced (frozen baseline) |
| ☐ | **Supp Fig 2** | Clades share feature pools, almost no exclusive features | `supp_fig2_clade.R` | Reproduced; taxonomy-dependent |
| ☐ | **Supp Fig 3** | MS2LDA motifs enriched across phyla | `supp_fig3_ms2lda.R` | Reproduced; ⚠ new strict MS2LDA rerun exists but is deliberately NOT connected to this legacy figure yet |
| ☐ | **Supp Fig 4** | NNLS leave-one-out confusion matrix (52.8% accuracy) | `supp_fig4_loo.R` | Reproduced; ⚠ the underlying 52.8% classifier producer is unrecovered |
| ☐ | **Supp Fig 5** | Negative control: pure isolates through decomposition | `supp_fig5_negative_control.R` | Reproduced (frozen baseline) |
| ☐ | **Supp Fig 6** | Negative-mode sensitivity vs positive mode | `supp_fig6_neg_sensitivity.R` | Reproduced (frozen baseline) |
| ☐ | **Supp Fig 7** | Annotation confidence pipeline + evidence tiers | `supp_fig7_annotation.R` | Reproduced; ⚠ annotation counts will change under the strict release (blocked steps documented) |
| ☐ | **Supp Fig 8** | Cross-method fingerprint validation | `supp_fig8_cross_method.R` | Reproduced (frozen baseline) |

## Supplementary tables (S1–S15)

All fifteen are built by one script: `manuscript_2_clean/07_tables/build_tables.py`.
Rebuilt copies exist in `outputs/tables/` alongside the as-submitted frozen workbook
(`SLA_Supplementary_Tables_SUBMITTED_FROZEN.xlsx`). Each table gets its own row/folder
during review; per-table inputs are documented then.

| ☐ | Unit | ☐ | Unit | ☐ | Unit |
|---|---|---|---|---|---|
| ☐ | Table S1 sample inventory | ☐ | Table S6 | ☐ | Table S11 cross-method Mantel |
| ☐ | Table S2 | ☐ | Table S7 | ☐ | Table S12 consensus per phylum |
| ☐ | Table S3 | ☐ | Table S8 | ☐ | Table S13 LOO accuracy |
| ☐ | Table S4 | ☐ | Table S9 | ☐ | Table S14 ClimGrass overlap |
| ☐ | Table S5 | ☐ | Table S10 ClimGrass benchmark | ☐ | Table S15 phylogeny cladogram |

⚠ Table S1 is where the reviewer-visible taxonomy error lives; its corrected
version depends on the strict release, so the review of S1 must present both
as-submitted and corrected states.

## Known leftovers (documented, not part of the manuscript mapping)

- `fig3_ordination.R` — PCoA ordination main-figure script, **not in the release
  mapping**; likely a demoted earlier main figure. Decide: archive with note.
- `figures_r/supp/` — an older generation of four supplementary scripts
  (`SuppFig1_nnls_confusion.R` etc.) with different numbering. Superseded; archive.
- `fig5_caption.md` is empty and no caption file exists for the ClimGrass figure —
  captions for those come from the submitted legends docx.

## Review order

Manuscript order: Figure 1 → 2 → 3 → 4 → 5 → Supp Fig 1–8 → Tables S1–S15 →
leftovers. One unit = one commit.
