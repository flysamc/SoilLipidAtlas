# Annotation Steps 2 and 4 — release `ncbi-phylum-2026-08-04-v1`

Step 2 (harmonisation): `scripts/step2_harmonization_release.py`
Step 4 (family propagation): `scripts/step4_family_propagation_release.py`

## Step 2: harmonisation

Priority hierarchy, with ArchLips as an archaeal-specific structural source
directly after LipidSearch Grade A/B (LipidSearch has no archaeal ether
lipids):

1. LipidSearch molecular species (Grade A/B)
2. ArchLips Gold/Silver spectral match
3. LipidSearch sum composition / lower grade
4. ArchLips Bronze
5. MS2 diagnostic class
6. Molecular family propagation → Step 4

| Mode | Features | Annotated | % | Gold+Silver | % | Classes | Superclasses |
|---|---|---|---|---|---|---|---|
| POS | 11,371 | 6,291 | 55.32 | 1,706 | 15.00 | 49 | 9 |
| NEG | 5,697 | 1,791 | 31.44 | 512 | 8.99 | 34 | 6 |

Sources — POS: LipidSearch 3,599, MS2 diagnostic 2,075, ArchLips 617, none 5,080.
NEG: LipidSearch 961, MS2 diagnostic 820, ArchLips 10, none 3,906.

Superclass distribution (POS): Glycerolipids 2,086; Glycerophospholipids 1,309;
Glycerophospholipids|Sphingolipids 948; Prenol lipids 823; Archaeal lipids 647;
Sphingolipids 243; Sterol lipids 108; Fatty acyls 101; Betaine lipids 26.

### Vocabulary choices

Recorded in `step2_harmonization/STEP2_MANIFEST.json`:

- **`Archaeal lipids` and `Betaine lipids`** extend beyond the six LIPID MAPS
  superclasses because this dataset requires them.
- **`PC/SM` → `Glycerophospholipids|Sphingolipids`.** The diagnostic engine
  calls this from m/z 184.07, the phosphocholine headgroup shared by PC and
  SM. The joint label keeps 948 POS features comparable to each other without
  implying a pure PC or pure SM assignment.
- **`No_MS2` and `Unknown` are excluded** as status flags, not lipid classes
  (removed 2 spurious NEG annotations).
- Four NEG classes added: `OAHFA` → Fatty acyls; `MGMG`, `DGMG`, `SQMG` →
  Glycerolipids. Unmapped count is zero in both modes.

## Step 4: molecular family propagation

Map `feature_id → (reference batch, cluster index) → component`, exclude
singleton families, require ≥1 annotated member and ≥50% label agreement,
upgrade Unidentified → Bronze with source `Network_propagation`, two rounds.

| Mode | Unidentified before | In families | Families | Upgrades (normalised) | Upgrades (verbatim) | Delta |
|---|---|---|---|---|---|---|
| POS | 5,080 | 1,441 | 1,081 | **155** | 154 | +1 |
| NEG | 3,906 | 1,651 | 1,187 | **98** | 98 | 0 |

Both modes converged after round 1; round 2 produced no further upgrades.
Normalising class-name spelling changes +1 POS upgrade and 0 NEG.

### Family coverage caps this step

Only **1,441 of 11,371 POS features (12.7%)** and **1,651 of 5,697 NEG
features (29.0%)** land in a non-singleton molecular family. Everything
else is a singleton or has no cluster-summary entry, so propagation can
reach only a small fraction of the Unidentified pool.

## Files

`annotation/step2_harmonization/` — `harmonised_annotations_{pos,neg}.csv`,
`step2_summary.csv`, `STEP2_MANIFEST.json`

`annotation/step4_family_propagation/` — `propagated_{pos,neg}_{normalised,verbatim}.csv`,
`step4_summary.csv`, `STEP4_MANIFEST.json`
