# Step 03 — Biomarker atlas

Finds the lipid features that mark each phylum. Two independent selection
methods are run on the 16-phylum sample set and their union is the atlas:

- **Composite scoring** — ranks features by how strongly they are enriched in
  one phylum (tiers: Platinum = strongest, then Silver).
- **IndVal (Indicator Value)** — classic indicator-species statistics, run for
  single phyla and phylum pairs.

**Result:** 11,371 biomarker features in positive mode, 5,697 in negative mode
(row counts in the atlas files; one header line each).

Reads: ← `01_taxonomy/` sample metadata.
Used by: → 04 annotation, 05 public validation, 06 fingerprints.

## Files

| File | What it is |
|---|---|
| `pos/atlas_pos_strict.csv`, `neg/strict_atlas_NEG.csv` | **The atlases** — one row per biomarker feature, with its phylum and how it was selected |
| `*/…composite_platinum…`, `*/…composite_silver…` | Composite-score selections by tier |
| `*/…indval_unique…`, `*/…indval_pairs…` | IndVal selections (single phylum / phylum pairs) |
| `*/…counts_by_phylum…` | Biomarkers per phylum, the headline table |
| `pos/discovery_method_by_phylum.csv` | Which method found what, per phylum |
| `*/…summary.json` | Run summary with selection counts |
| `pos/reproducibility_check.json`, `pos/freeze_manifest.json` | Verification + freeze record |
| `scripts/figure2a_strict.py` | Producer, positive mode |
| `scripts/negative_biomarker_pipeline.py` | Producer, negative mode |
| `scripts/freeze_biomarker_atlas.py` | Freezes the atlas so later steps use one fixed version |
