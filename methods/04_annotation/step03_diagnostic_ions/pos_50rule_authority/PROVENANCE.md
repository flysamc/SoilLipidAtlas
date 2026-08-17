# Provenance — the 50-rule diagnostic classifier (Table S3 authority)

This package is the reproducible correction for Supplementary Table S3.

## Files

- `table/Table_S3_diagnostic_ion_database_50_rules.xlsx` — corrected
  supplementary workbook.
- `table/Table_S3_diagnostic_ion_database_50_rules.csv` — machine-readable
  version.
- `documentation/CORRECTION_NOTE.md` — evidence, limitation, and recommended
  manuscript wording.
- `source/ms2_classify_indval.py` — the authoritative recovered executable
  classifier.
- `source/build_supp_table_S3.py` — the recovered table exporter.

## Authority statement

The 50-rule executable classifier is the reproducible source of truth for
Table S3. The manuscript's original claim of a "66-entry diagnostic
fragment-ion database" is not supported by any recoverable source: no
standalone 66-entry artifact was found in any archive, and the missing 16
entries cannot be reconstructed reliably from feature-level outputs.
Manuscript wording therefore references *"a reproducible 50-rule diagnostic
MS/MS classifier comprising 24 diagnostic fragment rules, five neutral-loss
rules, and 21 terpenoid-fragment rules across 13 classifier labels"*, with
the historical 66-entry source explicitly noted as unrecovered. The missing
rules were **not** invented or inferred.

A separate, older 33-entry quick-reference table is literature-oriented
background and must not be substituted for the executable classifier.
