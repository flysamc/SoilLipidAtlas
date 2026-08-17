# Supplementary Methods — index

The manuscript's Supplementary Methods (SM1–SM9) and Supplementary Note (SN1)
are each implemented by one of the numbered method steps in this repository.
This index maps every SM section to the folder that holds its code, data, and
plain-language description, with the release (strict-16) headline values.

| SM | Topic | Repository folder | Strict-16 release values |
|---|---|---|---|
| SM1 | Cross-batch feature alignment | [`../02_features/`](../02_features/) | 273,248 POS / 122,571 NEG consensus features; >15,000 anchors; LOESS residual 0.35–0.52 min |
| SM2 | Biomarker selection (composite score + IndVal) | [`../03_biomarker_atlas/`](../03_biomarker_atlas/) | 11,371 POS / 5,697 NEG; composite Platinum 2,480 / Silver 222 (selected); IndVal 7,751 unique (Streptophyta 5,167 dominant) |
| SM3 | Annotation pipeline (11 steps) | [`../04_annotation/`](../04_annotation/) | 6,360 POS / 1,889 NEG annotated; Step 10 SIRIUS 8,004 formulas (70.4%) / 7,336 CANOPUS / 7,199 CSI after the gap-fill; Steps 6 & 9 = declared exceptions |
| SM4 | Fingerprint extraction & cross-method validation | [`../06_fingerprints/`](../06_fingerprints/) (`cross_method/`) | 45,525-feature substrate, 164 samples, 16 phyla; SIMPER/SCBD exact, CAP bounded fit, L1 seeded reimplementation; LOO baseline 64.6% |
| SM5 | MS2LDA substructure-motif discovery | [`../06_fingerprints/`](../06_fingerprints/) | deterministic pipeline; 134 POS / 103 NEG significant motif–phylum pairs |
| SM6 | Public-data validation (fastMASST + Pan-ReDU) | [`../05_public_validation/`](../05_public_validation/) | 7,406 POS / 5,695 NEG queried; POS 3,978 public / 2,287 soil; IndVal 30.9% vs composite 14.6% soil-detected |
| SM7 | Archaeal lipid annotation (ArchLips) | [`../04_annotation/`](../04_annotation/) (`step08_archlips*/`) | 546 validated POS after RT screening (368 Gold+Silver); impossible-RT calls screened out |
| SM8 | Quantification correction & expected composition | [`../08_climgrass/`](../08_climgrass/) | Rule A (out-of-window RIE → uncalibrated 1.0) replaces the RIE floor; literature expectation ranges in Table S6 |
| SM9 | Negative control (pure-isolate decomposition) | [`../07_decomposition/`](../07_decomposition/) | n = 164; dominant group 79.3% corrected; archaeal self-recovery 93.1% |
| SN1 | Feature-count substrate reconciliation | cross-cutting (see below) | one table, single source of truth for every substrate size cited |

## SN1 — the substrate numbers, reconciled

Different analyses legitimately run on different views of the data. These are
the strict-16 values; any other number appearing in a legacy document is
superseded:

| Substrate | Size | Where used |
|---|---|---|
| Consensus aligned features | 273,248 POS / 122,571 NEG | SM1, feature tables |
| Quality-filtered analysis substrate | 45,525 features × 164 samples | dendrograms, SM4, Supp Fig 8 |
| Biomarker atlas | 11,371 POS / 5,697 NEG | SM2, Figure 2 |
| SIRIUS-eligible / queried | 7,406 POS queryable of 7,751 strict validation set | SM6 |
| ClimGrass verified-soil substrate | **736** features (Fig 5 v2; 722 SIMPER + 14 archaeol); 696 on the strict cross-method substrate (S14) | SM4, SM8, Fig 5 |

Legacy counts (44,534 / 830 / 769 / 12,710 …) belong to the superseded
pre-correction analysis and must not be cited.
