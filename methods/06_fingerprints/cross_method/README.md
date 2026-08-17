# Supplementary Figure 8 rerun — cross-method fingerprint validation, strict 16 phyla

Workspace: `outputs/analysis/ncbi-phylum-2026-08-04-v1/suppfig8_cross_method_strict16_2026-08-11_v1/`
Producer: `paper2_repro/scripts/suppfig8_cross_method_rerun.py` (seed 20260811, refuses overwrite)
Taxonomy: `ncbi-phylum-2026-08-04-v1` locked strict policy (16 analysis phyla)
Status: review-only. The published producer ("analysis-20", deleted worktree) is
unrecovered; methods were recovered by fitting the frozen per-method top-K union
sizes, with per-method validation modes below.

## Method recovery status (old-label validation, hard evidence)

| method | recovered form | vs frozen union sizes (5 fingerprints) |
|---|---|---|
| SIMPER | classic pairwise Bray-Curtis contribution (vegan form), raw intensities | **EXACT, hard gate** — all five match (800/1695/2921/5032/10851). NOT the centroid SIMPER of the atlas; the two variants coexist in the published work. |
| SCBD | per-phylum block of global SCBD sum-of-squares, Hellinger data | **EXACT, hard gate** (1191/2823/5294/9490/19290) |
| CAP | PCoA(BC on TSS), first 7 axes, one-vs-rest least-squares discriminant, features by abs Pearson r vs Hellinger data | bounded best fit, max deviation 5.5% (K=100); exact historical config not identifiable from available evidence |
| L1 stability | LARS lasso entry order over B=50 stratified half-subsamples, standardised Hellinger data, rank = frequency, then mean entry step, then abs point-biserial r | declared reimplementation; original was stochastic, deep tail deviates up to 24% at K=2500 (head deviation ~1% at K=100) |

## Corrected (strict-16) headline outputs

- Full substrate: 164 samples x 45,525 features; full-dendrogram cophenetic r
  and NNLS LOO baseline **106/164 = 64.6%** (identical code path to the gated
  `suppfig4_loo_strict16_2026-08-11_v1` run — panels are internally consistent).
- Panel a: SIMPER and SCBD reconstruct the full dendrogram at Mantel r > 0.98
  from K=100 and sit above the random-features null band; CAP intermediate;
  L1 below until large K.
- Panel b: L1 sets exceed the full-substrate baseline (up to ~75%); SIMPER/SCBD
  track the baseline; CAP underperforms.
- Panel c: verified-soil substrate is the corrected 5-ppm mapping
  (`strict16_verified_simper_mapping_5ppm.csv`): 696 of its 722 unique
  feature_ids sit on the strict-16 substrate (n=696 replaces the published 769,
  which used old labels and is not reconstructible). SCBD captures 96% at K=2500.
- Panel d: Evosea (36) and Discosea (25) lead all-4-method consensus at K=500;
  kingdoms follow the locked policy display groups.

## Caption changes required (vs submitted legend)

| item | published | revised | reason |
|---|---|---|---|
| dendrogram | full 44,534-feature dendrogram | full 45,525-feature dendrogram | corrected metadata maps 164 samples; quality mask recomputed |
| units | (17 phyla implicit) | 16 phyla | locked strict policy |
| panel b baseline | 57.1% | 64.6% (106/164), reimplementation | the published 57.1% (92/161) producer is unrecovered and is inconsistent with Supp Fig 4's published 52.8% on the same samples — flag in the response letter; the R script's "Full 44k-substrate baseline" legend label needs the one-word 44k->45k edit if adopted |
| panel c set | 769 ClimGrass spectrally verified features | 696 corrected 5-ppm verified features | published list not reconstructible; corrected mapping is the release-consistent replacement |
| methods | (as described) | CAP and L1 reconstructed | Methods must state SIMPER/SCBD reproduce the historical selection exactly; CAP is a bounded best-fit; L1 is a seeded reimplementation of a stochastic original |

## Files

- `panel_[a-d]_*.csv` (top level and under `r_render/data/supp_cross_method/`,
  identical) — schemas match `supp_fig8_cross_method.R` unchanged
- `RUN_SUMMARY.json` — validation tables, parameters, seeds
- `r_render/out/Supplementary_Fig8_cross_method.{png,pdf,svg}` — rendered figure
  (R script byte-identical to the repo copy; no edits were needed)

## Downstream

Table S12 (cross-method consensus per unit) can now be rebuilt from these
rankings — it was one of the three tables blocked on this producer.
