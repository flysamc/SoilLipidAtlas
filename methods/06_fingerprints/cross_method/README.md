# Cross-method fingerprint validation (Supplementary Figure 8), strict 16 phyla

Producer: `../scripts/suppfig8_cross_method_rerun.py` (seed 20260811, refuses
overwrite). Source workspace in the analysis archive:
`suppfig8_cross_method_strict16_2026-08-11_v1/` under release
`ncbi-phylum-2026-08-04-v1` (locked strict policy, 16 analysis phyla).

The historical producer behind the published panel was not recovered; the four
methods were reconstructed by fitting the frozen per-method top-K union sizes,
with the per-method validation status below.

## Method recovery status (validated against frozen historical outputs)

| Method | Recovered form | vs frozen union sizes (5 fingerprints) |
|---|---|---|
| SIMPER | classic pairwise Bray-Curtis contribution (vegan form), raw intensities | ✅ **exact, hard gate** — all five match (800/1695/2921/5032/10851). This is the pairwise variant, distinct from the centroid SIMPER used for the atlas; the two variants coexist in the published work |
| SCBD | per-phylum block of global SCBD sum-of-squares, Hellinger data | ✅ **exact, hard gate** (1191/2823/5294/9490/19290) |
| CAP | PCoA(BC on TSS), first 7 axes, one-vs-rest least-squares discriminant, features by abs Pearson r vs Hellinger data | 🔶 bounded best fit, max deviation 5.5 % (K=100); the exact historical configuration is not identifiable from available evidence |
| L1 stability | LARS lasso entry order over B=50 stratified half-subsamples, standardised Hellinger data; rank = frequency, then mean entry step, then abs point-biserial r | 🔶 declared reimplementation; the original was stochastic — deep tail deviates up to 24 % at K=2500 (head deviation ~1 % at K=100) |

## Strict-16 headline outputs

- Full substrate: 164 samples × 45,525 features; full-dendrogram cophenetic r
  and NNLS LOO baseline **106/164 = 64.6 %** (identical code path to the gated
  `suppfig4_loo_strict16_2026-08-11_v1` run — panels are internally consistent).
- **Panel a**: SIMPER and SCBD reconstruct the full dendrogram at Mantel
  r > 0.98 from K=100 and sit above the random-features null band; CAP
  intermediate; L1 below until large K.
- **Panel b**: L1 sets exceed the full-substrate baseline (up to ~75 %);
  SIMPER/SCBD track the baseline; CAP underperforms.
- **Panel c**: the verified-soil substrate is the corrected 5-ppm mapping
  (`strict16_verified_simper_mapping_5ppm.csv`): 696 of its 722 unique
  feature_ids sit on the strict-16 substrate. SCBD captures 96 % at K=2500.
- **Panel d**: Evosea (36) and Discosea (25) lead all-4-method consensus at
  K=500; kingdom groupings follow the locked display policy.

## Differences from the published panel

The published figure was built on the pre-correction taxonomy; these values
supersede it under the corrected labels:

| Item | Published | This release | Reason |
|---|---|---|---|
| Dendrogram substrate | 44,534 features | **45,525** features | corrected metadata maps 164 samples; quality mask recomputed |
| Analysis units | 17 phyla (implicit) | **16 phyla** | locked strict taxonomy policy |
| Panel-b baseline | 57.1 % | **64.6 %** (106/164) | the published baseline's producer is unrecovered — and is internally inconsistent with the published Supp. Fig. 4 value (52.8 %) on the same samples; the release value is a declared reimplementation |
| Panel-c feature set | 769 spectrally verified features | **696** corrected 5-ppm verified features | the published list is not reconstructible; the corrected mapping is the release-consistent replacement |
| Method provenance | (as described) | SIMPER/SCBD exact; CAP bounded fit; L1 seeded reimplementation | per-method recovery status above |

## Files

- `panel_[a-d]_*.csv` (top level and under `r_render/data/supp_cross_method/`,
  identical) — schemas match `supp_fig8_cross_method.R` unchanged
- `RUN_SUMMARY.json` — validation tables, parameters, seeds
- `r_render/out/Supplementary_Fig8_cross_method.{png,pdf,svg}` — rendered
  figure

## Downstream

Tables S11–S14 are built from these rankings (see `results/tables/`).
