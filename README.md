# SoilLipidAtlas

Code, methods documentation, and result assets for **SoilMass Paper 2** —
*Distributed lipid fingerprints encode phylum-level identity across the soil
food web*.

Soil organisms from across the food web — bacteria, archaea, fungi, plants,
animals, protists — were cultured, their lipidomes measured by untargeted
LC-MS/MS, and their phylum-level lipid signatures extracted, annotated, and
validated. The resulting reference atlas is then applied to real field soil
(the ClimGrass climate-change experiment) to decode community composition
from lipids alone.

**The study in numbers**

| | |
|---|---|
| Organisms cultured | 19 phyla collected → **16 phyla analysed** (≥2 samples in both ionisation modes) |
| Aligned features | 273,248 (positive mode) / 122,571 (negative mode) |
| Biomarker atlas | **11,371 POS / 5,697 NEG** phylum-enriched features |
| Annotation | 6,360 POS / 1,889 NEG annotated; SIRIUS formula coverage 70.4% POS |
| Public validation | 2,287 POS / 486 NEG biomarkers found in public soil datasets |
| Negative control | 79.3% correct dominant-group recovery over 164 pure isolates |
| Field application | 736 verified soil features; drought response replicates published qSIP results |

Analysis release: `ncbi-phylum-2026-08-04-v1` (locked — any label change
requires a new release ID).

## Repository map

```
README.md            you are here
PROGRESS.md          one-row-per-step status + open author decisions
environment.yml      pinned conda environment (Python 3.12, R 4.4)
methods/
  01_taxonomy/           sample → phylum labels (the locked foundation)
  02_features/           cross-batch alignment → one feature table per mode
  03_biomarker_atlas/    phylum-enriched feature selection
  04_annotation/         what each biomarker lipid is (11-step pipeline)
  05_public_validation/  biomarkers searched in the world's public MS/MS data
  06_fingerprints/       SIMPER fingerprints, MS2LDA motifs, cross-method checks
  07_decomposition/      source decomposition + pure-isolate negative control
  08_climgrass/          decoding real field soil; quantification correction
results/
  figures/               main figures 1–5 + supplementary 1–8, with producers
  tables/                supplementary tables S1–S15, with producers
data_registry/       checksummed index of files too large for git
```

## How to follow the analysis (recommended route)

1. **`PROGRESS.md`** — the one-page status: what each step does and where it
   stands.
2. **`methods/01_…` through `methods/08_…`**, in order. Every step folder
   contains:
   - a `README.md` that explains in plain language what the step does, what
     it reads, what it produces, and every number it contributes to the paper;
   - the producer scripts that built the released artifacts;
   - the step's key result tables (small enough for git), so claims can be
     checked without downloading anything.
3. **`results/`** — each figure and table sits next to the script that builds
   it; the READMEs there map every asset to its producer, its source data,
   and any open decision.
4. **`data_registry/registry.csv`** — everything too large for git, with
   SHA-256 checksums and its public home:
   - Raw LC-MS/MS data: **MassIVE MSV000102115**
   - Derived large tables: versioned **Zenodo** deposit (DOI
     10.5281/zenodo.20811187, reserved)

## Reproducibility

- **Environment.** `environment.yml` pins the Python/R environment the
  release was built with.
- **The analysis archive.** Producer scripts are committed exactly as they
  ran; they reference the authors' local analysis workspace (the *analysis
  archive*), whose released outputs live under the release ID above. To
  re-run a producer, obtain its registry-listed inputs (verify by SHA-256)
  and point the path constants at the top of the script to them. Generated
  run reports (`STEP*_REPORT.md`, manifests) cite producer paths as they
  were inside the archive (e.g. `paper2_repro/scripts/…`); the same scripts
  are committed in each step's `scripts/` folder here.
- **Reproduce-first policy.** Every rebuilt artifact was gated against the
  release before acceptance; gates are recorded in per-step
  `RUN_SUMMARY.json` / manifest files. Where a historical producer could not
  be recovered, its replacement is a **declared reimplementation**, labelled
  as such in the step README and never silently substituted.
- **Fail-closed producers.** Multi-day external searches (SIRIUS, fastMASST)
  run behind checkpointed, checksum-gated runners that refuse to summarise
  incomplete inputs.
- **Taxonomy labels.** Display labels follow the locked `ecological_group`
  policy (Viridiplantae, Protists). Some producer-emitted CSVs retain legacy
  internal keys (Plantae, Protozoa); this is noted wherever it occurs.

## Status

All methods steps (01–08), all figures, and all supplementary tables are
built. Open items awaiting author decisions are listed at the bottom of
`PROGRESS.md`.
