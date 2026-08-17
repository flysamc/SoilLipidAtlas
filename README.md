# SoilLipidAtlas

Analysis repository for **SoilMass Paper 2** — *Distributed lipid fingerprints
encode phylum-level identity across the soil food web*.

The analysis uses the **16-phylum NCBI taxonomy set** (release
`ncbi-phylum-2026-08-04-v1`, defined and locked in `methods/01_taxonomy/`):
19 phyla were collected, 16 enter statistics (≥2 samples in both ionisation
modes).

## How to follow this repository (reviewer's route)

1. **`PROGRESS.md`** — one row per methodology step: what it does, its status,
   and what of it lives in this repo. Start here.
2. **`methods/01_…` → `methods/08_…`** — the pipeline in the order of the
   paper's Methods section. Each folder has a plain-language `README.md`
   explaining what the step does, what it reads, what it produces, and every
   number it contributes to the paper — plus the producer scripts and the
   small key result tables themselves.
3. **`results/figures/`** and **`results/tables/`** — the rendered manuscript
   figures (main 1–5, supplementary 1–8) and supplementary tables (S1–S15),
   each next to the script that builds it. The READMEs there map every asset
   to its producer, its source data, and any open decision.
4. **`data_registry/registry.csv`** — files too big for git (raw spectra,
   full feature tables, bulk match sets) with sha256 checksums and their
   public home. Raw LC-MS/MS data: **MassIVE MSV000102115**. Derived big
   tables: versioned **Zenodo** deposit (DOI reserved; populated from this
   registry).

## Reproducibility notes

- **Environment:** `environment.yml` (conda; Python 3.12 + R 4.4 with pinned
  package versions) is the environment the release was built with.
- **Paths:** producer scripts record the exact code that built each release
  artifact and reference the analysis working tree they ran in (a local
  `P2R/` directory). To re-run one, place the registry-listed inputs at the
  paths named at the top of the script (or edit those constants) — inputs are
  verifiable by sha256 against `data_registry/registry.csv`.
- **Reproduce-first policy:** every rebuilt artifact was gated against the
  release before being accepted (gates recorded in each step's
  `RUN_SUMMARY.json` / manifest files). Where a historical producer could not
  be recovered, the replacement is a **declared reimplementation** and is
  labelled as such in the step README — never silently swapped.
- **Taxonomy policy:** organism-group display labels are Viridiplantae and
  Protists (`ecological_group` policy); some producer-emitted CSVs keep
  legacy internal keys (Plantae/Protozoa) — noted wherever they appear.

## Status

All methods steps (01–08), all figures, and all supplementary tables are
built. Remaining open decisions (coauthor input) are listed at the bottom of
`PROGRESS.md`.
