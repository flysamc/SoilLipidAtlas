# Figure 3 — chemotaxonomic structure

Phylum-level lipid dendrograms, their correspondence to SSU rRNA
evolutionary distance, and pairwise lipid similarity.

Methods, accessions, sequences, and the tree:
[`methods/09_phylogeny/`](../../../methods/09_phylogeny/).

| Panel | What it shows | Data in `r/data/` |
|---|---|---|
| a | UPGMA dendrograms of phylum lipidomes (Bray–Curtis; cophenetic r = 0.880 POS / 0.867 NEG) | `lipid_braycurtis.csv`, `phyla.csv` |
| b | Mantel scatter of SSU distance against lipid Bray–Curtis (r = 0.509 POS / 0.553 NEG). Triangles = cross-group pairs; circles = within-group | `panel_b_pairs.csv`, `mantel.csv` |
| c | POS and NEG lipid similarity heatmaps (1 − Bray–Curtis) | `panel_c_heatmap.csv`, `phyla.csv` |

Panels a and c use lipids only. Panel b’s x-axis is SSU patristic
distance, not the five-rank NCBI framework in Table S15.

Panel a dendrogram branches are rotated toward the fixed biological
display order of `r/data/phyla.csv` (Bacteria → Archaea → Protists →
Fungi → Viridiplantae → Animalia) so the POS and NEG trees and the
panel c heatmaps read in a comparable order. Rotation is
topology-preserving: UPGMA merge heights and cophenetic correlations
are unchanged.

`figure3_strict_composite_review_only.R` draws the three-panel layout
(183 × 195 mm). House style: `r/soilmass_style.R`.
