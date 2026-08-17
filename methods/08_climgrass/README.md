# Step 08 — ClimGrass: decoding real field soil

Applies the atlas fingerprints to the ClimGrass climate-change experiment
(2 × 2 factorial: warming × elevated CO₂, plus drought) — the test of whether
lab-derived fingerprints decode a real soil community.

## Pipeline (`scripts/`)

1. `climgrass_strict16_release.py` — fingerprints matched into ClimGrass
   soil MS/MS; 2,313 candidate matches reduced to a **736-feature
   verified-soil substrate**.
2. `climgrass_strict16_rulefix.py` — quantification correction. **Rule A**:
   features whose response-factor calibration is out-of-window get an
   explicit *uncalibrated 1.0* factor.
3. `rie_prediction.py` / `rie_wide_anchor.py` — empirical response (RIE)
   factor models behind Table S7.
4. `climgrass_strict16_archlips_extension.py` + `scan_soil_for_archlips.py` —
   archaeal coverage (ArchLips hits searched directly in soil).

## Results (`results/`)

- `kingdom_composition.csv` / `composition_fcweighted_*` — corrected
  community composition, fc-weighted primary: Bacteria 38.2%, Fungi 27.8%,
  Animalia 14.7%, Viridiplantae 9.6%, Archaea 5.5%, Protists 4.1% (bootstrap
  CIs alongside). **Animalia ≈14.7% is above the literature expectation
  range** (reference-panel bias). The marker-panel estimator disagrees
  with fc-weighted on Bacteria and Viridiplantae; both are shown in Fig. 5a.
- `phylum_effects.csv` — treatment effects per phylum (CLR + permutation
  tests). No phylum reaches omnibus FDR *q* < 0.05. Pseudomonadota drought
  *q* = 0.08 (16-phylum omnibus) / 0.005 (pre-specified qSIP directional
  family); Evosea is nominally drought-depleted (*p* = 0.050).
- `qsip_replication_test.csv` — the drought–Pseudomonadota signal
  replicates a published qSIP growth-suppression result (q = 0.005).
- `benchmark_summary.csv` + `fingerprint_set_effects.csv` — synthetic-mixture
  benchmark and fingerprint-set sensitivity.
- `RUN_SUMMARY.json` — parameters and checksums. The `kingdom_ci` block
  in that file is the **marker-panel** estimator (Bacteria 18.7%,
  Viridiplantae 28.6%), not the fc-weighted bars in Figure 5a.

Figure 5 renders `composition_fcweighted_kingdom_ci.csv` (bars),
`kingdom_ci_marker_panel.csv` (diamonds), and `phylum_effects.csv` via
`results/figures/figure5/r/fig5_final.R`. Large intermediates (soil MGFs,
benchmark mixtures) are listed in `data_registry/registry.csv`.
