# Shadow PC handoff

This package is the reproducible correction for SoilMass Paper 2 Supplementary Table S3.

## Use these files

- `table/Table_S3_diagnostic_ion_database_50_rules.xlsx` — corrected supplementary workbook for review/use.
- `table/Table_S3_diagnostic_ion_database_50_rules.csv` — machine-readable version.
- `documentation/CORRECTION_NOTE.md` — evidence, limitation, and recommended manuscript wording.
- `source/ms2_classify_indval.py` — authoritative recovered executable classifier.
- `source/build_supp_table_S3.py` — recovered table exporter.

## Instruction for an agent on the Shadow PC

Treat the 50-rule classifier as the current reproducible source of truth. Replace manuscript references to a “66-entry diagnostic fragment-ion database” with “a reproducible 50-rule diagnostic MS/MS classifier comprising 24 diagnostic fragment rules, five neutral-loss rules, and 21 terpenoid-fragment rules across 13 classifier labels.” Do not invent or infer the missing 16 rules. Preserve the note that the historical 66-entry source was not recovered.

The older 33-entry quick-reference table is literature-oriented background and must not be substituted for the executable classifier.
