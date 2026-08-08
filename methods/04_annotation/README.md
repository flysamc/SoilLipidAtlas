# Step 04 — Annotation

Steps 01–03 found *which* features mark each phylum. This step asks: **what
lipid is each feature actually?** Each answer gets a confidence tier:
**Gold** (molecular species) > **Silver** (class + partial structure) >
**Bronze** (class only) > Unidentified.

Because the 16-phylum atlas re-selected biomarkers from corrected labels,
6,636 of 11,371 POS and 3,521 of 5,697 NEG features had never been annotated
before — old caches only partially transfer.

## Current release-facing result

**6,360 annotated POS features (1,668 Gold+Silver) · 1,889 NEG (512 Gold+Silver)**
— assembled by `scripts/annotation_summaries.py` with blocked steps declared as
exceptions, never guessed. Headline tables in `summaries/`
(`tier_counts.csv`, `waterfall.csv`), full per-feature evidence in `evidence/`.

## The eleven sub-steps + three methods

| # | What it does | Result | Folder |
|---|---|---|---|
| 1 | LipidSearch vendor IDs (434 exported reports) | POS 3,619 annotated (1,304 grade A/B); NEG 961 — up from 74 cached (new producer) | `step01_lipidsearch/` |
| 2 | Harmonise naming across sources | POS 6,291 (55.3%), 1,706 G+S; NEG 1,791, 512 G+S | `step02_harmonization/` |
| 3 | MS2 diagnostic-ion rules | 2,075 Bronze class candidates (POS). Manuscript claimed 66 rules; **50 recovered and used**, the missing 16 not invented; 2,168 non-matches stay explicit | `step03_diagnostic_ions/` |
| 4 | Molecular-family propagation | +155 POS, +98 NEG upgrades | `step04_family_propagation/` |
| 5 | RT validation (impossible retention times) | 86 archaeal calls elute at 1.8–9.3 min vs ~18 min real — all from ArchLips, all screened out | `step05_rt_validation/` |
| 6 | RT → sum-composition prediction | Historical producer lost (1,764 predictions unreproducible). Conservative replacement v1: **1 candidate (PE 32:1), 0 tier changes**. NEG blocked | `step06_rt_prediction_pos_v1/` |
| 7 | Custom 124-entry archaeal DB | ❌ database absent from every archive | — |
| 8 | ArchLips archaeal library | POS 546 validated after RT screen (27.6% G+S — use this, not unscreened 30.5%); NEG 10 Bronze, 0 G+S, zero diagnostic ions | `step08_archlips*/` |
| 9 | Evidence-depth → tier corrections | Audit done; the producer behind the submitted 564+640+12 tier changes is lost — **not reproduced, not faked** | in `evidence/` |
| 10 | SIRIUS / CANOPUS / CSI:FingerID (LISC cluster) | POS 5,865 formulas / 5,325 classes / 5,243 structures; NEG 3,862 / 3,495 / 3,325 — hash-validated, exact-ID integrated | in `evidence/` |
| 11 | DreaMS deep-learning similarity | 7,095 POS + all 5,695 usable NEG features | in `evidence/` |
| M5 | MS2LDA substructure motifs | Deterministic (bit-for-bit across reruns); 134 POS / 103 NEG significant motif–phylum pairs | in `evidence/` |
| M6 | fastMASST + Pan-ReDU public search | ⏳ running (13,103 queries; resumed 2026-08-08) | step 05 of this repo |

`scripts/` holds all producers written for this release. Step-6 recovery
evidence (35 MB) stays in P2R (`biomarker_discovery/rt_prediction_recovery_pos/`).

## Open items

1. ⏳ fastMASST/Pan-ReDU completion (running).
2. Lost producers blocking restated numbers: Step 5's "104 RT-uncertain",
   Step 6 POS caller, Step 7 database, Step 9 tier-action producer.
3. Decisions needed (Rahul): Step 1 mapping choice (direct 1,304 A/B —
   recommended — vs stale-network 592); Step 5 grouping wording (text says
   family, recovered code groups by class); how to phrase unvalidated NEG
   archaeal annotation.
4. Manuscript numbers update only after the release completes.
