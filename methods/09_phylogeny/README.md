# Step 09 — SSU phylogeny and Figure 3 pairing

This step builds the organism small-subunit (SSU) rRNA tree used as the
**x-axis of Figure 3b**, then pairs it with phylum-level lipid distances
from the quality-filtered feature table (not the biomarker atlas).

**Headline result:** Mantel r = **0.509 (POS)** / **0.553 (NEG)** between
lipidome Bray–Curtis and SSU evolutionary distance (16 phyla, 120 pairs,
9,999 permutations, P < 0.001). Partial Mantel controlling for batch-set
overlap: **0.510 / 0.559**.

Freeze: `figure3-ssu-curated-freeze-2026-08-04-v3`
(148 sequences, 103 analysis units, 16 phyla).

Reads: ← `01_taxonomy/` organism list; lipid distances from 45,525 POS
features × 164 samples and 14,896 NEG × 192.
Used by: → Main Figure 3.

## Files in this folder (one copy of each essential object)

| File | What it is |
|---|---|
| `ssu_accessions.csv` | Locked sequence list: accession, SILVA tip, length, proxy flag, selection rule (148 rows) |
| `ssu_sequences.fasta` | Those 148 SSU sequences |
| `ssu_tree.nwk` | GTR+Γ ML tree with Felsenstein bootstrap % on nodes |
| `SILVA_LICENSE.txt` | SILVA 138.2 is CC BY 4.0 |
| `scripts/` | The R producers, in pipeline order |
| `../results/figures/figure3/r/data/` | Combined plotted tables for panels a–c (see below) |

Plotted data (not duplicated here):

| File | Panel |
|---|---|
| `phyla.csv` | 16 phyla: display order, group, sample counts |
| `lipid_braycurtis.csv` | POS+NEG lipid distances (120 pairs × 2 modes) — dendrograms 3a and heatmaps 3c |
| `panel_b_pairs.csv` | SSU distance vs lipid distance, with within/cross-group class — scatter 3b |
| `mantel.csv` | Mantel and partial Mantel, both modes, primary metric only |
| `panel_c_heatmap.csv` | Lipid similarity = 1 − Bray–Curtis, in heatmap order |

Intermediates (v2 freeze, coverage-audit CSVs, six sensitivity distance
matrices, unit-level 103×103 matrices, unmasked alignment) stay in the
analysis archive. They are not needed to read the figure.

## What Figure 3 is, and what it is not

| Panel | What it plots | Distance source |
|---|---|---|
| **3a** | Two UPGMA dendrograms (POS, NEG) | **Lipids only.** Phylum-centroid Bray–Curtis. Cophenetic r = 0.880 / 0.867 |
| **3b** | Mantel scatter: evolution vs lipids | **x = SSU patristic distance** (this step). y = the same lipid Bray–Curtis as 3a. Hollow triangles = cross-group; filled circles = within-group (ecological display groups, not a phylogenetic rank) |
| **3c** | POS / NEG similarity heatmaps | **Lipids only.** Similarity = 1 − Bray–Curtis from 3a. No point markers |

Two other distance objects in this repository are **not** the Figure 3b
x-axis:

- **Table S15** is the five-rank NCBI framework (Domain > Supergroup >
  organism group > clade > phylum). It orders Supplementary Figures 1–2
  and Table S4.
- An NCBI hierarchy edge-count script exists in the analysis archive as
  a diagnostic. It was not used for Figure 3.

The primary 3b aggregation is
`anchor_set_mean__inventory_weighted`: for every pair of analysis units,
mean patristic distance across all SSU anchor sequences, then roll up
to phyla weighted by each unit’s sample count.

---

## Data sources (where the sequences came from)

The lipid table is measured in this study. The **tree is not**: each
cultured isolate is represented by a public SSU record of the same
NCBI taxon (or an explicit same-genus proxy — step 6).

| Source | What was taken | How it was obtained | In git? |
|---|---|---|---|
| **This study** | Verified NCBI TaxIDs of the isolates (105 unique taxa before curation) | `../01_taxonomy/` | taxon list yes |
| **SILVA 138.2 SSU Ref NR99** | Clipped 16S/18S already truncated to the SSU gene | SILVA archive, CC BY 4.0 (`SILVA_LICENSE.txt`) | full 510,495-tip FASTA **no**; the 148 used sequences **yes** (`ssu_sequences.fasta`) |
| **NCBI Nucleotide** | Records for taxa absent from SILVA; annotated 16S/18S feature extracted | E-utilities `efetch` (`db=nuccore&rettype=fasta`, tool `SoilMassFigure3`), batches of 15, 0.4 s pause | selected sequences **yes** (same FASTA). Raw EFetch cache **no** |

SILVA files that were downloaded but not stored here (too large):

- `SILVA_138.2_SSURef_NR99_tax_silva_trunc.fasta.gz`
- `SILVA_138.2_SSURef_NR99.accessions.ntree`
- `taxmap_ncbi_ssu_ref_nr99_138.2.txt.gz` and the two tax tables

NCBI URL pattern:

```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text&tool=SoilMassFigure3&id=<accession1,accession2,...>
```

---

## Steps taken

Scripts in `scripts/` were run from the analysis archive; path constants
inside them still point there. Random seed for alignment, ML fit, and
the first bootstrap chunk: **20260804**.

### 1. Organism list (105 unique TaxIDs)

Locked taxonomy `ncbi-phylum-2026-08-04-v1`. LC-MS replicates are not
tree tips: one verified NCBI TaxID is one analysis unit. 19 collection
phyla → 16 analysis phyla. Figure 3 lipid matrices use 164 POS / 192 NEG
samples after dropping Cercozoa, Cyanobacteriota, and Mortierellomycota.

### 2–3. Coverage audits (not stored; numbers only)

SILVA 138.2 NR99 vs the 105 TaxIDs: 55 species-or-better, 19 genus/proxy
pools, 74 with any SILVA candidate, **31 with none**. NCBI Nucleotide
then found a full-length SSU candidate for 88 / 105 taxa and rescued
20 of those 31. Audits did not pick accessions.

### 4. Accession candidate freeze v2 (145 records / 105 units)

Script: `scripts/freeze_figure3_ssu_accessions.R`.

Rules:

- NCBI taxonomy is the identity of the cultured isolate.
- Prefer a SILVA 138.2 NR99 clipped SSU when the verified taxon is
  present (target ~1,500 nt bacteria/archaea, ~1,800 nt eukaryotes).
- NCBI Nucleotide only for exact taxa missing from SILVA; must be a
  16S/18S marker, not an organelle genome. RefSeq `NR_` preferred.
- A **genus-labelled** unit is represented by up to three named
  descendant-species anchors. Later distance is the mean of those
  anchors, not one arbitrary tip.
- A genus with fewer than two named descendant SSU records stays a
  single explicit exemplar.
- Uncultured / unidentified / `sp.` / `cf.` names are rejected.

v2: 115 SILVA + 30 NCBI; 81 exact tips; 20 genus anchor-sets; 2
single-exemplar genera; 2 approved genus proxies.

### 5. Collect sequences

Script: `scripts/collect_figure3_ssu_sequences_rough.R`
(filename says “rough”; this is the collector that was run).

- SILVA rows: look up `silva_tree_tip` in the truncated NR99 FASTA,
  convert U → T.
- NCBI rows: EFetch accession.version in two cached batches of 15,
  four retries with exponential backoff.

### 6. Curated freeze v3 — the locked set (148 sequences / 103 units)

Script: `scripts/curate_figure3_ssu_freeze_v3.R`.
Length gate: SSU **≥ 900 nt** (observed minimum 939 nt).

| Action | What happened |
|---|---|
| Exact replacements | *Hyaloscypha finlandica* (2482753) → `L76625.1`; *Petroselinum crispum* (4043) → `AH001742.2` (annotated 18S feature joined) |
| Same-genus SILVA proxy (exact record < 900 nt) | *Auricularia*, *Cephalotrichum*, *Peziza*, *Amanita*, *Taraxacum*, *Armillaria* |
| Excluded | *Heydenia* (931642), *Warcupia* (352928) — no defensible SSU/proxy. 105 → **103** units |
| Distance rule | Former `genus_mrca` labels renamed `genus_anchor_set_mean`. Primary distance averages all anchors |

Locked files: `ssu_accessions.csv`, `ssu_sequences.fasta`.

### 7. Structural alignment

Script: `scripts/align_figure3_ssu_full.R` (DECIPHER `AlignSeqs`,
structures on, 2 iterations, 1 refinement). 4,996 columns → **1,593**
kept (occupancy ≥ 0.50 and variable). Alignment file not stored; the
tree below is the product.

### 8. Maximum-likelihood tree

Script: `scripts/fit_figure3_ssu_tree_full.R`. F81 distances → NJ start
→ GTR+Γ (4 categories), NNI, log-likelihood −58,918.2. A stochastic
ratchet was tried and not kept (did not bound). **NNI is the declared
search.**

### 9. 100 bootstrap replicates

Scripts: `bootstrap_figure3_ssu_chunk.R` then
`finalize_figure3_ssu_tree_v3.R`. Four chunks × 25; seeds 20261804,
20262804, 20263804, 20264804. Felsenstein bootstrap % labelled on
`ssu_tree.nwk`. Four multi-anchor units are non-monophyletic; they are
kept and flagged, not forced into a clade.

### 10. Phylum evolutionary distances

Between analysis units: **primary** = mean patristic distance over all
anchor-sequence pairs; medoid and MRCA are sensitivities (not plotted).
Unit distances are then inventory-weighted to 16 phyla. Those 16×16
values are the `evolutionary_distance` column of `panel_b_pairs.csv`.

### 11. Lipid distances (Figure 3a / 3c; no SSU)

Scripts: `figure3a_strict.R`, `figure3c_strict_heatmap.R`. Phylum
centroids → Bray–Curtis → UPGMA. Heatmap similarity = 1 − Bray–Curtis,
scale 0–0.55. Tables: `lipid_braycurtis.csv`, `panel_c_heatmap.csv`.

### 12. Mantel (Figure 3b)

Script: `apply_figure3_ssu_distance_v3.R`. 9,999 permutations. Seeds:
Mantel POS 20720247, NEG 21211386; partial POS 20720248, NEG 21211387.
Partial covariate: batch-set Jaccard `|B_i ∩ B_j| / |B_i ∪ B_j|`.

| | Mantel r | P | Partial r | P |
|---|---:|---:|---:|---:|
| POS (n = 164) | 0.509 | < 0.001 | 0.510 | < 0.001 |
| NEG (n = 192) | 0.553 | < 0.001 | 0.559 | < 0.001 |

### 13. Composite figure

`results/figures/figure3/figure3_strict_composite_review_only.R`
pastes the three R panels onto one 183 × 195 mm canvas.

## Software

R 4.6.1; ape 5.8.1; phangorn 2.12.1; DECIPHER; Biostrings; ggplot2 4.0.3.

## Not in git

- SILVA 138.2 NR99 truncated FASTA and guide tree (~510k tips)
- NCBI EFetch cache
- ML fit `.rds` and bootstrap-chunk `.rds`
- Full aligned feature tables (step 02 / Zenodo)
- Sensitivity distance matrices (medoid, MRCA, unweighted)
