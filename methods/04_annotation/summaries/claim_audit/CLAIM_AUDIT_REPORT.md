# Annotation numbers — release `ncbi-phylum-2026-08-04-v1`

Current values from the step outputs. Machine-readable table:
`annotation_claim_audit.csv`. Re-running the producer refreshes the file.

| Quantity | This analysis |
|---|---|
| LipidSearch Grade A/B | POS 11.47 %; NEG 8.99 % |
| Diagnostic-ion rule set | **50** executable rules (24 fragment, 5 neutral-loss, 21 terpenoid) |
| Family-propagation upgrades | POS 155 + NEG 98 |
| RT-uncertain (family-grouped IQR) | POS 11; NEG 23 |
| Archaeal Gold+Silver (ArchLips, after RT screen) | **27.6 %** POS (30.5 % before screen); NEG 0 |
| NEG archaeal MS2 support | 0 diagnostic ions, 0 neutral losses (21 spectra) |
| SIRIUS / CANOPUS / CSI | POS 8,004 / 7,336 / 7,199; NEG 3,914 / 3,531 / 3,352 |
| DreaMS | POS 7,095; NEG all 5,695 usable |
| MS2LDA documents | POS 10,611; NEG 5,466 |
| RT sum-composition → tier change | 1 candidate reviewed (PE 32:1); **0 tier changes** |

Archaeal Gold+Silver is reported **after** the step-5 retention-time
screen: 86 POS ArchLips annotations elute 8–16 minutes earlier than the
archaeal median of 17.92 min, which is chromatographically impossible for
intact ether lipids.

Analysis phyla follow `taxonomy_policy.json` (Methanobacteriota,
Thermoproteota; not Euryarchaeota/Crenarchaeota). ArchLips assignments:
879 Methanobacteriota + 216 Thermoproteota.

Step reports with the underlying tables:
`../step01_lipidsearch/STEP1_REPORT.md`,
`../step05_rt_validation/STEP5_REPORT.md`,
`../step08_archlips/STEP8_REPORT.md`,
and the step-04 README.
