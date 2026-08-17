# Cross-method fingerprint validation (Supplementary Figure 8)

Producer: `../scripts/suppfig8_cross_method_rerun.py` (seed 20260811).
Sixteen analysis phyla, locked release `ncbi-phylum-2026-08-04-v1`.

Four methods build fingerprints on the same substrate. Configurations and
top-K union sizes are recorded in `RUN_SUMMARY.json`.

| Method | Configuration | Top-K union sizes (5 fingerprints) |
|---|---|---|
| SIMPER | Pairwise Bray-Curtis contribution (vegan form), raw intensities | 800 / 1,695 / 2,921 / 5,032 / 10,851 |
| SCBD | Per-phylum block of global SCBD sum-of-squares, Hellinger data | 1,191 / 2,823 / 5,294 / 9,490 / 19,290 |
| CAP | PCoA (Bray-Curtis on TSS), first 7 axes, one-vs-rest least-squares discriminant; features by abs Pearson r vs Hellinger data | see `RUN_SUMMARY.json` |
| L1 stability | LARS lasso entry order over B = 50 stratified half-subsamples, standardised Hellinger data; rank = frequency, then mean entry step, then abs point-biserial r | see `RUN_SUMMARY.json` |

The atlas SIMPER (centroid form, step 06) is a different variant from the
pairwise SIMPER used here.

## Results

- Substrate: 164 samples × **45,525** features. Leave-one-out baseline
  **106/164 = 64.6 %** (same code path as Supplementary Fig. 4).
- **Panel a.** SIMPER and SCBD reconstruct the full dendrogram at Mantel
  r > 0.98 from K = 100 and sit above the random-features null; CAP
  intermediate; L1 below until large K.
- **Panel b.** L1 sets exceed the full-substrate baseline (up to ~75 %);
  SIMPER/SCBD track the baseline; CAP is weakest.
- **Panel c.** 696 of 722 unique verified-soil feature IDs sit on this
  substrate (`strict16_verified_simper_mapping_5ppm.csv`). SCBD captures
  96 % at K = 2,500.
- **Panel d.** Evosea (36) and Discosea (25) lead all-4-method consensus at
  K = 500. Organism-group labels follow the locked display policy
  (Viridiplantae, Protists).

## Files

- `panel_[a-d]_*.csv` (top level and under `r_render/data/supp_cross_method/`,
  identical)
- `RUN_SUMMARY.json` — parameters, seeds, union sizes
- `r_render/out/Supplementary_Fig8_cross_method.{png,pdf,svg}` — rendered
  figure

Tables S11–S14 are built from these rankings (`results/tables/`).
