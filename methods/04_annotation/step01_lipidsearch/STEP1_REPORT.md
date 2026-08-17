# Annotation Step 1 (LipidSearch) — release `ncbi-phylum-2026-08-04-v1`

Supplementary Method 3, Step 1: LipidSearch database matching.
Producer: `scripts/step1_lipidsearch_release.py` (this step's `scripts/`
folder). Checksums and parameters: `STEP1_MANIFEST.json`.

Matching parameters: 5 ppm m/z, 0.3 min RT, RT window 1.5–25 min.
Tiers: Gold = Grade A/B, IDScore ≥ 0.20, ≥ 2 samples; Silver = IDScore ≥ 0.10
with Grade A/B or ≥ 2 samples; Bronze = IDScore ≥ 0.01.

## How each mode is mapped

| Mode | Method |
|---|---|
| POSITIVE | Direct mapping from raw LipidSearch exports onto all 273,248 POS consensus features (`consensus_lipidsearch_direct.csv`) |
| NEGATIVE | Same matching rules on 177 `*NEG.raw.txt` exports → 74,106 rows → 18,179 clusters (17,214 after reject) → 122,571 NEG consensus features |

## Coverage

| Feature set | Features | Annotated | % | Grade A/B | % | Gold | Silver | Bronze |
|---|---|---|---|---|---|---|---|---|
| POS strict atlas | 11,371 | 3,619 | 31.83 | 1,304 | 11.47 | 979 | 1,334 | 1,306 |
| POS newly selected | 6,636 | 2,275 | 34.28 | 880 | 13.26 | 636 | 752 | 887 |
| NEG strict atlas | 5,697 | 961 | 16.87 | 512 | 8.99 | 394 | 316 | 251 |

Grade A/B rates are the Gold-eligible floor from this step alone, not the
final atlas tier rate.

A network-routed POS alternative yields 592 Grade A/B. This analysis uses
the **direct** mapping (1,304). The choice is recorded in
`STEP1_MANIFEST.json`.

## NEG mapping checks

All matches respect the tolerances (max 4.998 ppm, max 0.300 min; medians
0.926 ppm and 0.113 min). Adducts are negative-mode with no positive-mode
leakage: `M-H` 5,481, `M+HCOO` 2,844, `M-CH3` 312, `M-2H` 12. Classes are
NEG-appropriate (PC via formate adduct, FA, Cer, PG, PE, PS, CL, DGDG,
MGDG, PI).

## Archaeal coverage is low — expected

LipidSearch does not cover archaeal ether lipids (GDGT, archaeol,
hydroxy-archaeol, caldarchaeol). Archaeal annotation in this analysis comes
from ArchLips (step 8).

| Phylum | POS features | POS annotated | POS Grade A/B | NEG features | NEG annotated |
|---|---|---|---|---|---|
| Methanobacteriota | 1,017 | 58 (5.70%) | 14 (1.38%) | 18 | 1 (5.56%) |
| Thermoproteota | 316 | 22 (6.96%) | 1 (0.32%) | 3 | 0 (0.00%) |

Heterolobosea is also sparsely annotated here: 31 POS features at 12.90%
with **0** Grade A/B, and 561 NEG features at **0.53%** with 0 Grade A/B.

## Files

| File | Contents |
|---|---|
| `pos_strict_atlas_lipidsearch.csv` | 11,371 strict POS features + LipidSearch fields |
| `pos_newly_selected_lipidsearch.csv` | 6,636 newly selected POS features + LipidSearch fields |
| `neg_consensus_lipidsearch_direct.csv` | all 122,571 NEG consensus features, 8,649 matched |
| `neg_strict_atlas_lipidsearch.csv` | 5,697 strict NEG features + LipidSearch fields |
| `*_by_phylum.csv` | Per-phylum coverage for each set |
| `step1_coverage_summary.csv` | The coverage table above |
| `STEP1_MANIFEST.json` | Parameters and input/output sha256 |
