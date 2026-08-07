# SoilLipidAtlas

Analysis repository for **SoilMass Paper 2** — *Distributed lipid fingerprints
encode phylum-level identity across the soil food web* (coauthor phase).

The analysis uses the **16-phylum NCBI taxonomy set** (release
`ncbi-phylum-2026-08-04-v1`, defined in `methods/01_taxonomy/`).

## How this repository is organized

- **`PROGRESS.md`** — where the analysis stands, one row per step.
- **`methods/`** — the methodology in order (01, 02, …). Each folder explains in
  plain language what the step does, what it reads, and what it produces.
- **`figures/`** — manuscript figures, added as the steps behind them finish.
- **`data_registry/`** — big files (raw spectra, bulk outputs) don't live in
  git; they live on MassIVE / Zenodo and are listed here with checksums so any
  copy can be verified.

Content is added one reviewed step at a time.
