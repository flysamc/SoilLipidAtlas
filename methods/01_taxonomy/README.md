# Step 01 — Taxonomy

Every sample gets one phylum label under **NCBI Taxonomy**. Everything downstream
(atlas, annotation, fingerprints, ClimGrass) uses these labels, so this step
comes first.

**The rule:** a phylum enters the analysis only if it has ≥2 samples in *both*
ionisation modes. 19 phyla were collected → **16 phyla are analysed**
(Cercozoa, Cyanobacteriota, Mortierellomycota have too few samples).
Core samples: 168 positive mode, 195 negative mode.

Release ID: `ncbi-phylum-2026-08-04-v1` (locked — label changes require a new
release ID, never an edit).

## Files

| File | What it is |
|---|---|
| `taxonomy_policy.json` | The full rulebook: NCBI as primary, PR2 for protists, name mappings (e.g. Euryarchaeota → Methanobacteriota), exclusion rules |
| `taxonomy_summary.json` | The result in numbers: phyla lists, per-phylum sample counts |
| `sample_metadata_POS_ncbi_phylum.csv` | Sample → phylum table, positive mode |
| `sample_metadata_NEG_ncbi_phylum.csv` | Sample → phylum table, negative mode |
| `taxonomy_audit.csv` | Per-sample audit: how each label was verified |

Used by: → every later step.
