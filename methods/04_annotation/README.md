# Step 04 — Annotation

Steps 01–03 found *which* features mark each phylum. This step asks: **what
lipid is each feature actually?** Each answer gets a confidence tier:
**Gold** (molecular species) > **Silver** (class + partial structure) >
**Bronze** (class only) > Unidentified.

## Result

**6,360 annotated POS features (1,668 Gold+Silver) · 1,889 NEG (512 Gold+Silver)**
— assembled by `scripts/annotation_summaries.py`. Headline tables in
`summaries/` (`tier_counts.csv`, `waterfall.csv`); full per-feature evidence
in `evidence/`.

## The eleven sub-steps + three methods

| # | What it does | Result | Folder |
|---|---|---|---|
| 1 | LipidSearch vendor IDs (434 exported reports) | POS 3,619 annotated (1,304 grade A/B); NEG 961. POS uses the **direct** mapping from raw exports | `step01_lipidsearch/` |
| 2 | Harmonise naming across sources | POS 6,291 (55.3 %), 1,706 G+S; NEG 1,791, 512 G+S | `step02_harmonization/` |
| 3 | MS2 diagnostic-ion rules | 2,075 Bronze class candidates (POS); **50-rule** executable classifier; 2,168 non-matches stay explicit | `step03_diagnostic_ions/` |
| 4 | Molecular-family propagation | +155 POS, +98 NEG upgrades | `step04_family_propagation/` |
| 5 | RT validation (impossible retention times) | 86 archaeal ArchLips calls elute at 1.8–9.3 min vs ~18 min expected — screened out | `step05_rt_validation/` |
| 6 | RT → sum-composition prediction | Screened; **1 candidate (PE 32:1), 0 tier changes**. Not applied in NEG | `step06_rt_prediction_pos_v1/` |
| 7 | Custom archaeal compound list | Not used. Archaeal spectral annotation is ArchLips (step 8) | — |
| 8 | ArchLips archaeal library | POS 546 validated after RT screen (**27.6 %** Gold+Silver); NEG 10 Bronze, 0 Gold+Silver, zero diagnostic ions | `step08_archlips*/` |
| 9 | Evidence-depth tier transitions | Not applied. Tiers remain those from steps 1–5 and 8 | in `evidence/` |
| 10 | SIRIUS / CANOPUS / CSI:FingerID | POS 8,004 formulas (**70.4 %**) / 7,336 classes / 7,199 structures; NEG 3,914 / 3,531 / 3,352 | in `evidence/` + `step10_sirius_gapfill/` |
| 11 | DreaMS deep-learning similarity | 7,095 POS + all 5,695 usable NEG features | in `evidence/` |
| M5 | MS2LDA substructure motifs | Deterministic; 134 POS / 103 NEG significant motif–phylum pairs | in `evidence/` |
| M6 | fastMASST + Pan-ReDU public search | 13,103 queries | step 05 of this repo |

`scripts/` holds the producers. SIRIUS input MGFs and primary result TSVs
are in `step10_sirius_gapfill/` (large files in the data registry).

## SIRIUS / CANOPUS / CSI coverage

Every eligible atlas biomarker (usable MS2, m/z ≤ 850) was searched:
formula, CANOPUS class, and CSI:FingerID structure against the database
union BIO ∪ PubChem ∪ HMDB ∪ GNPS ∪ YMDB ∪ PLANTCYC ∪ KNAPSACK.

POS: **8,004** formulas (70.4 %), 7,336 CANOPUS classes, 7,199 structures.
NEG: 3,914 / 3,531 / 3,352. Per-phylum Bacteria and Fungi SIRIUS coverage
is 74–88 %. SIRIUS evidence does **not** change Gold/Silver/Bronze tiers,
so the annotated totals (6,360 POS / 1,889 NEG) and every tier-based
figure stay as assigned by steps 1–5 and 8.

## Notes

1. **LipidSearch POS mapping.** The direct mapping from raw vendor exports
   is used (1,304 Grade A/B on the atlas). A network-routed alternative
   (592 Grade A/B) is not used; the choice is recorded in
   `step01_lipidsearch/STEP1_MANIFEST.json`.
2. **Negative-mode archaeal annotation has no MS2 support** (0 diagnostic
   ions and 0 neutral losses across all 21 spectra). It is reported as
   mass-accuracy-only, not as a validated identification
   (`step08_archlips/STEP8_REPORT.md`).
