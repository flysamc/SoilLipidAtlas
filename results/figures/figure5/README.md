# Figure 5 — ClimGrass composition and treatment responses

Corrected organism-group composition of ClimGrass soils, and phylum-level
responses to drought and warming.

| Panel | What it shows | Data in `r/data/` |
|---|---|---|
| a | Two-estimator community composition (fc-weighted bars + 95% CI; marker-panel diamonds; literature ranges) | `composition_fcweighted_kingdom_ci.csv`, `kingdom_ci_marker_panel.csv` |
| b | CLR fingerprint-set drought and warming effects for 16 phyla | `phylum_effects.csv` |

`r/fig5_final.R` draws the two-panel layout (183 × 150 mm). House style:
`r/soilmass_style.R`. Organism-group colours and legend names are set in
the R script (Viridiplantae, Protists). CSV keys remain Plantae / Protozoa.

The fc-weighted estimator is primary. Animalia ≈ 14.7% lies above the
literature range (Table S6). No phylum reaches omnibus FDR *q* < 0.05;
Pseudomonadota drought *q* = 0.08. Archaeal calibration is noted in the
Results, not on the figure.

Legend of record: `LEGEND.md`. Upstream methods: `methods/08_climgrass/`.
`r/data/kingdom_composition.csv` is per-sample composition for Table S8
and is not plotted here.
