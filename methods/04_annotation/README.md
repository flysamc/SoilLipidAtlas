# Step 04 — Annotation

Steps 01–03 found *which* features mark each phylum. This step asks: **what
lipid is each feature actually?** Each answer gets a confidence tier:
**Gold** (molecular species) > **Silver** (class + partial structure) >
**Bronze** (class only) > Unidentified.

Because the 16-phylum atlas re-selected biomarkers from corrected labels,
6,636 of 11,371 POS and 3,521 of 5,697 NEG features had never been annotated
before — historical annotation caches only partially transfer.

## Release-facing result

**6,360 annotated POS features (1,668 Gold+Silver) · 1,889 NEG (512 Gold+Silver)**
— assembled by `scripts/annotation_summaries.py` with unrecoverable steps
declared as exceptions, never guessed. Headline tables in `summaries/`
(`tier_counts.csv`, `waterfall.csv`), full per-feature evidence in `evidence/`.

## The eleven sub-steps + three methods

| # | What it does | Result | Folder |
|---|---|---|---|
| 1 | LipidSearch vendor IDs (434 exported reports) | ✅ POS 3,619 annotated (1,304 grade A/B); NEG 961 — up from 74 cached (new producer) | `step01_lipidsearch/` |
| 2 | Harmonise naming across sources | ✅ POS 6,291 (55.3 %), 1,706 G+S; NEG 1,791, 512 G+S | `step02_harmonization/` |
| 3 | MS2 diagnostic-ion rules | ✅ 2,075 Bronze class candidates (POS); the executable classifier holds **50 rules** (the manuscript's 66-entry count is unsupported — see `step03_diagnostic_ions/`); 2,168 non-matches stay explicit | `step03_diagnostic_ions/` |
| 4 | Molecular-family propagation | ✅ +155 POS, +98 NEG upgrades | `step04_family_propagation/` |
| 5 | RT validation (impossible retention times) | ✅ 86 archaeal calls elute at 1.8–9.3 min vs ~18 min expected — all from ArchLips, all screened out | `step05_rt_validation/` |
| 6 | RT → sum-composition prediction | ⚠ declared exception — historical producer lost (1,764 predictions not reproducible). Conservative replacement v1: **1 candidate (PE 32:1), 0 tier changes**. NEG not applicable | `step06_rt_prediction_pos_v1/` |
| 7 | Custom 124-entry archaeal DB | ⚠ declared exception — database absent from every archive | — |
| 8 | ArchLips archaeal library | ✅ POS 546 validated after RT screen (27.6 % G+S — the screened value, not the unscreened 30.5 %); NEG 10 Bronze, 0 G+S, zero diagnostic ions | `step08_archlips*/` |
| 9 | Evidence-depth → tier corrections | ⚠ declared exception — audit complete; the producer behind the submitted 564+640+12 tier changes is lost, so those changes are **not reproduced and not faked** | in `evidence/` |
| 10 | SIRIUS / CANOPUS / CSI:FingerID | ✅ **complete after gap-fill:** POS 8,004 formulas (70.4 %) / 7,336 classes / 7,199 structures; NEG 3,914 / 3,531 / 3,352 — hash-validated, exact-ID integrated | in `evidence/` + `step10_sirius_gapfill/` |
| 11 | DreaMS deep-learning similarity | ✅ 7,095 POS + all 5,695 usable NEG features | in `evidence/` |
| M5 | MS2LDA substructure motifs | ✅ deterministic (bit-for-bit across reruns); 134 POS / 103 NEG significant motif–phylum pairs | in `evidence/` |
| M6 | fastMASST + Pan-ReDU public search | ✅ complete (13,103 queries) | step 05 of this repo |

`scripts/` holds all producers written for this release. Step-6 recovery
evidence (35 MB) is retained in the analysis archive
(`biomarker_discovery/rt_prediction_recovery_pos/`).

## The SIRIUS submission-gap fill

An audit found only ~52 % of POS strict biomarkers had ever been submitted to
SIRIUS (Bacteria 3 %, Fungi 10 % of their biomarkers), while SIRIUS succeeds on
98 % of what it receives — the low per-phylum coverage was a **submission gap,
not an instrument limit**. The 2,209 SIRIUS-eligible unsubmitted POS features
(usable MS2, m/z ≤ 850) + 80 NEG were exported
(`scripts/build_sirius_gapfill_package.py`), run on the LISC compute cluster
(SLURM scripts + sha256 manifest in `step10_sirius_gapfill/`), and integrated
by exact feature ID. CSI:FingerID initially returned 0 structures in both modes
(a structure-DB cache fault during the run); a corrected structure-only
re-search over the BIO ∪ PubChem ∪ HMDB ∪ GNPS ∪ YMDB ∪ PLANTCYC ∪ KNAPSACK
union recovered POS 3,911 / NEG 27 structures. Primary per-adduct result TSVs
are in `step10_sirius_gapfill/results/` (top-5 variants and input MGFs are in
the data registry).

Net effect: POS SIRIUS coverage 5,865 → **8,004** biomarkers (70.4 %), CANOPUS
5,325 → 7,336, CSI 5,243 → 7,199; per-phylum Bacteria/Fungi SIRIUS coverage
rose from near-zero to 74–88 %. Validations: gap-fill features are disjoint
from previously-submitted and a subset of the strict atlas. Confidence-**tier**
counts are unchanged (no tier transition is inferred from Step-10 evidence
alone), so the release-eligible annotated totals (6,360 POS / 1,889 NEG) and
every tier-based figure are unaffected.

## Declared limitations

1. **Unrecoverable historical producers** (steps 5 "104 RT-uncertain" count,
   6 POS caller, 7 database, 9 tier-action producer): each is a declared
   exception. The affected historical numbers are documented, not restated as
   reproduced — see the step reports for the evidence in each case.
2. **Two LipidSearch POS mappings exist** and disagree (direct 1,304 Grade A/B
   vs network-routed 592); the direct mapping is used, with the choice and its
   rationale recorded in `step01_lipidsearch/STEP1_MANIFEST.json`.
3. **Negative-mode archaeal annotation has no MS2 support** (0 diagnostic ions
   and 0 neutral losses across all 21 spectra); it is therefore reported as
   mass-accuracy-only and never as validated identification
   (`step08_archlips/STEP8_REPORT.md`).
