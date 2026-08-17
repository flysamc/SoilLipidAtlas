# Figure 5 — legend of record

**Fig. 5 | Lipid fingerprints recover community composition in ClimGrass soils and a drought-associated bacterial decline.**

**a**, Mean fraction of matched soil lipid signal attributed to each organism group in 12 ClimGrass soils (positive ionisation; 736 spectrally verified features across 16 phyla). Bars are the fold-change-weighted (fc-weighted) estimate after quantification correction, with 95% bootstrap confidence intervals; diamonds are a marker-panel estimate restricted to phylum-specific markers. Grey segments are literature biomass ranges (Table S6). Bacteria (38.2%) and fungi (27.8%) fall within those ranges; Animalia (14.7%) lies above the 1–5% expectation. The two estimators disagree most for Bacteria (38.2% versus 18.7%) and Viridiplantae (9.6% versus 28.6%).

**b**, Drought (x) and warming (y) log₂ fold changes for the 16 analysis phyla (*n* = 12 soils). Point size is mean attributed abundance; a heavier outline marks unadjusted *p* < 0.05. No phylum reaches false-discovery-rate *q* < 0.05. Pseudomonadota declines under drought (log₂ fold change −0.60, *p* = 0.005, *q* = 0.08); Evosea is nominally drought-depleted (*p* = 0.050). Warming effects are not significant.

Source data are provided as a Source Data file.

## Plotted data in this repository

The legend numbers are those drawn by `r/fig5_final.R` from:

| File | Panel |
|---|---|
| `r/data/composition_fcweighted_kingdom_ci.csv` | a (bars, 95% CI) |
| `r/data/kingdom_ci_marker_panel.csv` | a (diamonds) |
| `r/data/phylum_effects.csv` | b |

`r/data/kingdom_composition.csv` is the per-sample fc-weighted composition used by Table S8; it is not drawn by this figure. CSV organism-group keys remain Plantae / Protozoa; the script remaps them to Viridiplantae / Protists for display.
