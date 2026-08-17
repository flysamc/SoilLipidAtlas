# Figure 3 — chemotaxonomic structure

Dendrograms, lipidome–phylogeny Mantel scatter, and POS/NEG similarity
heatmaps.

**How the SSU tree was built** (SILVA + NCBI accessions, alignment,
GTR+Γ, Mantel): [`methods/09_phylogeny/`](../../../methods/09_phylogeny/).
That folder holds one accession table, one FASTA, and one Newick.
**Plotted numbers** are here, combined into five CSVs.

The assembler `figure3_strict_composite_review_only.R` pastes three
already-rendered R panels onto one 183 × 195 mm canvas. It requires
panel b to be the SSU v3 freeze
(`anchor_set_mean__inventory_weighted`).

| Panel | What is plotted | File in `r/data/` |
|---|---|---|
| a | POS/NEG UPGMA dendrograms from phylum-centroid lipid Bray–Curtis (cophenetic r = 0.880 / 0.867) | `lipid_braycurtis.csv`, `phyla.csv` |
| b | Mantel scatter: SSU distance vs lipid Bray–Curtis. r = 0.509 (POS) / 0.553 (NEG). Triangles = cross-group; circles = within-group | `panel_b_pairs.csv`, `mantel.csv` |
| c | POS/NEG lipid similarity heatmaps (1 − Bray–Curtis). No point markers | `panel_c_heatmap.csv`, `phyla.csv` |

Table S15 (five-rank NCBI framework) is **not** the x-axis of panel b.
