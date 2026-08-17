# Figure 4 — distributed fingerprints

SIMPER features needed for 50% phylum separation, lipid-class composition
of those fingerprints, and how much of the distance matrix is recovered
from a small feature subset.

| Panel | What it shows | Data in `r/data/` |
|---|---|---|
| a | SIMPER features to 50% separation, 16 analysis phyla | `simper_curves.csv` |
| b | Lipid-class heatmap in positive and negative mode | `class_composition.csv`, `class_order.json` |
| c | Distance-matrix correlation versus features retained | `subsampling_curve.csv` |

`r/fig4.R` draws the three-panel layout (183 × 160 mm). House style:
`r/soilmass_style.R`. Organism-group colours and legend names are set in
the R script (Viridiplantae, Protists). CSV kingdom columns are not used.

Panel b uses the **full atlas** as the denominator, so the Unclassified
row is part of the plotted fraction. The heatmap is ordered by organism
group; cell gutters and group separators mark the grid.
