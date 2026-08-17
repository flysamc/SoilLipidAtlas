# Diagnostic-ion classifier (Table S3)

Supplementary Table S3 is the **50-rule positive-mode diagnostic
classifier**:

- 24 diagnostic fragment rules
- 5 diagnostic neutral-loss rules
- 21 terpenoid-fragment rules
- 50 rules in total, spanning 13 classifier labels when `PS` and
  `Terpenoid` are included

A 33-entry literature quick-reference table is background material, not
the executable classifier. The feature-level file of 8,362 classifications
is an output of the classifier, not a definition of it.

## Methods wording

> a reproducible 50-rule diagnostic MS/MS classifier comprising 24
> diagnostic fragment rules, five neutral-loss rules, and 21
> terpenoid-fragment rules across 13 classifier labels (Supplementary
> Table S3)

## Files

- `source/ms2_classify_indval.py` — executable definition
- `source/build_supp_table_S3.py` — exporter that flattens the classifier dictionaries
- `table/Table_S3_diagnostic_ion_database_50_rules.csv` — flattened rules
- `table/Table_S3_diagnostic_ion_database_50_rules.xlsx` — supplementary workbook
