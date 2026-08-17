# Step 09 — Chemotaxonomic structure

Tests whether phylum-level lipid distances recover evolutionary
relationships across the soil food web (Main Figure 3).

Lipid distances are computed from the quality-filtered feature table
(45,525 features × 164 samples in positive mode; 14,896 × 192 in
negative mode). Evolutionary distances come from a small-subunit (SSU)
rRNA tree of the same isolates. The two matrices are compared by Mantel
test: **r = 0.509 (POS) / 0.553 (NEG)** over 16 phyla (120 pairs;
9,999 permutations, P < 0.001). Partial Mantel, controlling for shared
LC-MS batches, is **0.510 / 0.559**.

## Figure 3

| Panel | Content |
|---|---|
| **a** | UPGMA dendrograms of phylum-centroid lipidomes (Bray–Curtis). Cophenetic r = 0.880 (POS) / 0.867 (NEG). Lipids only. |
| **b** | Mantel scatter: SSU evolutionary distance (x) versus lipid Bray–Curtis (y). Hollow triangles are cross-group pairs; filled circles are within-group pairs. Groups are ecological display strata (Bacteria, Archaea, Fungi, Viridiplantae, Animalia, Protists), not a phylogenetic rank. |
| **c** | Pairwise lipid similarity heatmaps (1 − Bray–Curtis). Lipids only; no point markers. |

Table S15 and Supplementary Figures 1–2 use a five-rank NCBI framework
(rank counts to the last common ancestor). That matrix is **not** the
x-axis of panel b.

Plotted tables: `results/figures/figure3/r/data/`.

## Evolutionary tree

Each cultured isolate is represented by a public SSU record of the same
NCBI taxon. The lipid measurements are from this study; the sequences
are not.

**Sources.** SILVA 138.2 SSU Ref NR99 (CC BY 4.0; clipped 16S/18S) is
used when the verified taxon is present. NCBI Nucleotide supplies 16S/18S
records for taxa missing from SILVA (organelle genomes excluded; RefSeq
`NR_` preferred). Sequences shorter than 900 nt are not used (observed
minimum 939 nt). Uncultured, unidentified, and `sp.`/`cf.` names are
rejected.

**One tip per taxon, except genera.** LC-MS replicates are not tree
tips: one NCBI TaxID is one analysis unit. A unit identified only to
genus is represented by up to three named descendant species; its
distance to another unit is the mean of all pairwise patristic distances
among those anchors. Six units whose exact records failed the length
gate use a named same-genus SILVA sequence instead (*Auricularia*,
*Cephalotrichum*, *Peziza*, *Amanita*, *Taraxacum*, *Armillaria*).
*Heydenia* and *Warcupia* had no usable SSU or proxy and are omitted
(105 taxa → **103 units**, **148 sequences**, 16 phyla).

**Inference.** Structural SSU alignment (DECIPHER); columns with
occupancy ≥ 0.50 that vary among A/C/G/T are retained (1,593 of 4,996).
Maximum-likelihood tree: GTR+Γ (four rate categories), NNI search,
100 bootstrap replicates (Felsenstein support on `ssu_tree.nwk`).
Random seed 20260804. Four multi-anchor units are non-monophyletic on
this tree; they are retained rather than forced into a clade.

**Phylum distances.** Unit distances are averaged to 16 phyla, weighted
by each unit’s sample count in the inventory. That inventory-weighted
anchor-set mean is the Figure 3b x-axis.

## Lipid distances and Mantel

Phylum centroids on the aligned feature table → Bray–Curtis → UPGMA
(panel a) and similarity 1 − Bray–Curtis (panel c, colour scale 0–0.55).

Mantel: Pearson correlation of the two 16 × 16 distance matrices, with
9,999 permutations of the evolutionary labels. Partial Mantel
residualizes both vectors against batch-set Jaccard overlap
(|batches in common| / |batches in either phylum|) before correlating.

Producers: `scripts/figure3a_strict.R`, `figure3c_strict_heatmap.R`,
`apply_figure3_ssu_distance_v3.R`. Sequence and tree producers:
`freeze_figure3_ssu_accessions.R`, `collect_figure3_ssu_sequences_rough.R`,
`curate_figure3_ssu_freeze_v3.R`, `align_figure3_ssu_full.R`,
`fit_figure3_ssu_tree_full.R`, `bootstrap_figure3_ssu_chunk.R`,
`finalize_figure3_ssu_tree_v3.R`.

## Files

| File | What it is |
|---|---|
| `ssu_accessions.csv` | Accession, SILVA tip, length, and selection rule for each of the 148 sequences |
| `ssu_sequences.fasta` | The 148 SSU sequences |
| `ssu_tree.nwk` | ML tree with bootstrap support |
| `SILVA_LICENSE.txt` | SILVA 138.2 licence (CC BY 4.0) |
| `scripts/` | R producers |
| `../../results/figures/figure3/r/data/phyla.csv` | Display order, organism group, sample counts |
| `../../results/figures/figure3/r/data/lipid_braycurtis.csv` | Lipid distances (panels a, c) |
| `../../results/figures/figure3/r/data/panel_b_pairs.csv` | SSU vs lipid distance for each phylum pair |
| `../../results/figures/figure3/r/data/mantel.csv` | Mantel and partial Mantel (both modes) |
| `../../results/figures/figure3/r/data/panel_c_heatmap.csv` | Lipid similarities in heatmap order |

The full SILVA 138.2 NR99 collection (~510k sequences) is not stored
here; only the 148 sequences used in the tree are. Feature tables for
the lipid distances are the step-02 consensus tables (Zenodo).

Software: R 4.6.1; ape 5.8.1; phangorn 2.12.1; DECIPHER; Biostrings.
