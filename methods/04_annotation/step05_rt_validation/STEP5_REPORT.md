# Annotation Step 5 (RT validation) — release `ncbi-phylum-2026-08-04-v1`

Producer: `scripts/step5_rt_validation_release.py`. Parameters and
checksums: `STEP5_MANIFEST.json`.

Tukey fences at 1.5 × IQR; severity from the group z-score (|z| > 3 severe,
|z| > 2 moderate, else mild). Groups with fewer than 4 members are skipped:
an IQR over 2–3 points is not meaningful.

## Results

This analysis flags RT-uncertain features **within molecular families**.
A class/superclass grouping is also emitted as a comparison (same files,
`grouping` column).

| Mode | Features | In families | Families | By family (used) | By superclass (comparison) |
|---|---|---|---|---|---|
| POS | 11,371 | 1,441 | 1,081 | **11** | 97 |
| NEG | 5,697 | 1,651 | 1,187 | **23** | 11 |

Severity — POS by family: 9 moderate, 2 mild. POS by superclass: 68
moderate, 16 severe, 13 mild. NEG by family: 13 moderate, 9 mild, 1 severe.
NEG by superclass: 10 moderate, 1 severe.

Most families are too small for an IQR: only 29 of 1,081 POS families
(2.7%) and 38 of 1,187 NEG families (3.2%) have 4 or more members. That is
why the family-grouped counts are low.

Separately, 86 POS archaeal ArchLips annotations elute at 1.8–9.3 min
versus ~18 min expected for intact ether lipids; those calls are screened
out (see step 8).

## Files

| File | Contents |
|---|---|
| `rt_uncertain_pos.csv` | POS flags, both groupings, with severity, z-score, fences, phylum |
| `rt_uncertain_neg.csv` | NEG equivalent |
| `step5_summary.csv` | The results table above |
| `STEP5_MANIFEST.json` | Parameters and checksums |
