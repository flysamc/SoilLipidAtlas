# Step 08 — ClimGrass: decoding real field soil

Applies the atlas fingerprints to the ClimGrass climate-change experiment
(2 × 2 factorial: warming × elevated CO₂, plus drought) — the test of whether
lab-derived fingerprints decode a real soil community.

## Pipeline (`scripts/`)

1. `climgrass_strict16_release.py` — base release: strict-16 fingerprints
   matched into ClimGrass soil MS/MS; 2,313 candidate matches reduced to a
   **736-feature verified-soil substrate**.
2. `climgrass_strict16_rulefix.py` — quantification correction. **Rule A**:
   features whose response-factor calibration is out-of-window get an
   explicit *uncalibrated 1.0* factor (replaces the earlier RIE hard floor).
3. `rie_prediction.py` / `rie_wide_anchor.py` — empirical response (RIE)
   factor models behind Table S7.
4. `climgrass_strict16_archlips_extension.py` + `scan_soil_for_archlips.py` —
   archaeal coverage extension (ArchLips hits searched directly in soil).

## Results (`results/`, from the v2 ArchLips release)

- `kingdom_composition.csv` / `composition_fcweighted_*` — corrected
  community composition, fc-weighted primary: Bacteria 38.2%, Fungi 27.8%,
  Animalia 14.7%, Viridiplantae 9.6%, Archaea 5.5%, Protists 4.1% (bootstrap
  CIs alongside). ⚠ Known caveat: **Animalia ≈14.7% is above the literature
  expectation range** — explicitly flagged in the manuscript as a
  reference-panel bias. The winner-take-all estimator disagrees with
  fc-weighted on several groups; both estimators are reported.
- `phylum_effects.csv` — treatment effects per phylum (CLR + permutation
  tests). Drought: only two FDR-supported phylum responses; Pseudomonadota
  q = 0.08 omnibus / 0.005 directional.
- `qsip_replication_test.csv` — the drought–Pseudomonadota signal replicates
  the published qSIP growth-suppression result (q = 0.005).
- `benchmark_summary.csv` + `fingerprint_set_effects.csv` — synthetic-mixture
  benchmark and fingerprint-set sensitivity.
- `RUN_SUMMARY.json` — full provenance.

Figure 5 (results section) renders `kingdom_composition.csv` +
`phylum_effects.csv` via `results/figures/figure5/r/fig5_final.R`.

Big intermediates (soil MGFs, benchmark mixtures, superseded correction-rule
workspaces) are retained in the analysis archive; they join the Zenodo
deposit via `data_registry/registry.csv` where needed.
