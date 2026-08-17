# SoilMass Paper 2 — diagnostic-ion correction note

## The correction

Use the recovered, executable **50-rule positive-mode diagnostic classifier** as the reproducible source for Supplementary Table S3. Do not retain the unsupported statement that the database contains 66 entries.

The recovered classifier contains:

- 24 diagnostic fragment rules;
- 5 diagnostic neutral-loss rules;
- 21 terpenoid-fragment rules;
- 50 rules in total, spanning 13 distinct classifier labels when `PS` and `Terpenoid` are included.

The previous curated quick-reference table contains only 33 literature-reference entries. It is useful background material, but it is not the complete executable classifier. The feature-level file containing 8,362 classifications is an output and is not a definition of the diagnostic database.

## What remains unavailable

No standalone source containing the manuscript's claimed 66 entries was found in the project repository, coauthor packages, submission bundle, or producer-recovery archive. The missing 16 entries cannot be reconstructed reliably from feature-level outputs. The 66-entry statement is therefore unsupported by the recoverable evidence.

## Wording adopted in the revised manuscript

Wording equivalent to:

> a 66-entry diagnostic fragment-ion database

is replaced with:

> a reproducible 50-rule diagnostic MS/MS classifier comprising 24 diagnostic fragment rules, five neutral-loss rules, and 21 terpenoid-fragment rules across 13 classifier labels (Supplementary Table S3)

Where the Methods summary said “66 entries covering 12 lipid classes,” the revised text uses:

> 50 diagnostic rules across 13 classifier labels

## Files and authority

- `source/ms2_classify_indval.py` — authoritative recovered executable definition.
- `source/build_supp_table_S3.py` — recovered exporter that flattens the classifier dictionaries.
- `table/Table_S3_diagnostic_ion_database_50_rules.csv` — machine-readable flattened rules.
- `table/Table_S3_diagnostic_ion_database_50_rules.xlsx` — corrected supplementary workbook.

This correction does not invent the missing rules and does not imply that the historical 66-entry set has been recovered.
