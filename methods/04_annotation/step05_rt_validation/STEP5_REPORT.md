# Annotation Step 5 (RT validation) — release `ncbi-phylum-2026-08-04-v1`

Producer: `paper2_repro/scripts/step5_rt_validation_release.py`
Parameters and checksums: `STEP5_MANIFEST.json`

Tukey fences at 1.5 × IQR, severity from the group z-score (|z| > 3 severe, |z| > 2
moderate, else mild). Groups with fewer than 4 members are skipped: an IQR computed
over 2–3 points is meaningless.

## Results

| Mode | Features | In families | Families | **By family (documented)** | By superclass (comparison) |
|---|---|---|---|---|---|
| POS | 11,371 | 1,441 | 1,081 | **11** | 97 |
| NEG | 5,697 | 1,651 | 1,187 | **23** | 11 |

Severity — POS by family: 9 moderate, 2 mild. POS by superclass: 68 moderate, 16 severe,
13 mild. NEG by family: 13 moderate, 9 mild, 1 severe. NEG by superclass: 10 moderate,
1 severe.

## The submitted description of this step is probably wrong

Supplementary Method 3 Step 5 states:

> "Retention time coherence **within molecular families** was assessed using IQR-based
> outlier detection. **104 features** flagged as RT-uncertain."

Implemented exactly as written, that method **cannot produce a number near 104**. The
reason is structural: only **29 of 1,081 POS families (2.7%)** and **38 of 1,187 NEG
families (3.2%)** contain 4 or more members. 1,052 POS and 1,149 NEG families were
skipped for being too small to support an IQR at all. The documented method yields
**11 (POS) and 23 (NEG)**.

The class-grouped variant gives **97 in POS** — close to the reported 104. That, plus
the fact that the only recovered RT producer (`step12_rt_validation.py`) groups by
`h_class` and reports `class_median_rt`, points to the historical analysis having been
**grouped by lipid class, not by molecular family**, with the supplementary text
describing it incorrectly.

## Producer status: unrecovered, and 104 is unverified

No located producer reproduces 104:

| Candidate | Grouping | Result |
|---|---|---|
| `step12_rt_validation.py` | lipid class | 577 anomalies |
| `platinum_unified_annotations_rt_validated.csv` | expected RT range per class | 209 violations, 243 `rt_valid=False` |
| documented family-grouped method | molecular family | not located |

**The 104 figure must not be restated in the revision until its producer is recovered
or the number is regenerated under this release.** Both variants computed here are
emitted side by side, labelled by their `grouping` column, so neither is mistaken for
a reproduction of it.

## Recommended action for the revision

Either recover the Step 5 producer, or correct the Methods sentence to state the
grouping actually used and report the regenerated count. This is the third unrecovered
producer in the package, alongside the editable Figure 1 source and the 52.8%
classifier.

## Files

| File | Contents |
|---|---|
| `rt_uncertain_pos.csv` | POS flags, both groupings, with severity, z-score, fences, phylum |
| `rt_uncertain_neg.csv` | NEG equivalent |
| `step5_summary.csv` | The results table above |
| `STEP5_MANIFEST.json` | Parameters, producer-status warning, candidate comparison |
