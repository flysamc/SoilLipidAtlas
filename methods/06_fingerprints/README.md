# Step 06 — Distributed fingerprints (SIMPER, MS2LDA, cross-method validation)

Beyond single biomarkers, each phylum has a **distributed fingerprint**: the
weighted set of features that most explains its compositional distinctness.
This step builds those fingerprints and shows the conclusion does not depend
on the method used to build them.

## SIMPER fingerprints (`simper/`, producer `scripts/simper.py`)

Per-phylum SIMPER over the 16-phylum atlas.
`reproduction_validation.csv` records that recomputed contributions match
the stored atlas. The full per-feature table (`simper_atlas.csv`, ~7 MB)
and K-vs-recovery curves stay in the data registry.

## MS2LDA substructure motifs (`scripts/run_ms2lda_deterministic.py`)

Deterministic MS2LDA (fixed seeds, sorted inputs) on the atlas MS2 spectra,
then motif–phylum enrichment
(`../04_annotation/scripts/strict_ms2lda_phylum_enrichment.py`):
**134 POS / 103 NEG significant motif–phylum pairs**. The enrichment table
feeds Supplementary Fig. 3.

## Cross-method validation (`cross_method/`, producer `scripts/suppfig8_cross_method_rerun.py`)

Four fingerprint-construction methods (SIMPER, SCBD, CAP, L1) compared on
the same substrate (164 samples × 45,525 features). Panel data
(Supplementary Fig. 8, Tables S11–S14):

- `panel_a_mantel_curves.csv` + `panel_a_random_null.csv` — fingerprint-based
  Mantel r vs fingerprint size K, against a random-selection null (observed
  curves sit ≈33 SD above the null, P < 0.001).
- `panel_b_loo_accuracy.csv` — leave-one-out phylum assignment by method × K
  (L1 best ≈74%; CAP weakest) → Table S13. Full-substrate baseline
  **106/164 = 64.6 %**.
- `panel_c_climgrass_overlap.csv` — overlap with the 696-feature verified
  ClimGrass soil substrate (SCBD captures 96% at K = 2,500) → Table S14.
- `panel_d_consensus_per_phylum.csv` — cross-method consensus features per
  phylum → Table S12.

`RUN_SUMMARY.json` carries substrates, seeds, and method configurations.
