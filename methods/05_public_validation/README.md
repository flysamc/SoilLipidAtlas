# Step 05 — Public validation (fastMASST + Pan-ReDU)

Searches every atlas biomarker against public MS/MS data.
Final summaries in `results/`; per-feature match sets (GBs) in the data
registry.

## What this step does

1. **T1 — fastMASST**: each biomarker spectrum is searched against the GNPS2
   Pan-Repo index (`metabolomicspanrepo_index_latest`, cosine ≥ 0.7).
   13,103 queries total: 7,408 POS + 5,695 NEG.
2. **T2 — Pan-ReDU**: maps those matches to sample-type metadata (which
   kingdoms/sample types contain the biomarkers, including 149 soil
   datasets). Summaries are written only when both T1 queues are complete,
   checksummed, and non-empty.

## Results (database `metabolomicspanrepo_index_latest`)

| Mode | Strict features | Queryable | Any public match | Any soil match |
|---|---|---|---|---|
| POS | 7,751 | 7,406 | 3,978 (54 %) | **2,287 (31 %)** |
| NEG | 5,697 | 5,695 | 2,283 (40 %) | **486 (8.5 %)** |

Per-kingdom and per-phylum breakdowns: `results/panredu_summary.csv`
(NEG rows still carry Plantae/Protozoa keys; display labels are
Viridiplantae/Protists). Checksums:
`results/PUBLIC_VALIDATION_MANIFEST.json`.

A second queue, restricted to composite-selected features, gives the
Figure 2c comparison: Indicator Value 7,254 features / 30.9 % soil-detected
versus composite 3,613 / 14.6 %.

## Negative-mode charge

The FASST API expects precursor charge as a **magnitude**. A signed
`precursor_charge: -1` returns empty matches. This analysis queries NEG
features with charge `1`. Completeness checks verify result *content*, not
only that a file exists.

## Runner

Each query writes its own status file and a compact gzip match file
(USI, Dataset, Cosine, Matching Peaks, Delta Mass). A supervisor retries
durable errors until every expected success is present. Parameters and
input checksums: `run_manifest_POS.json` / `run_manifest_NEG.json`.

## Files

| File | What it is |
|---|---|
| `run_manifest_POS.json`, `run_manifest_NEG.json` | Query-queue provenance |
| `scripts/build_fastmasst_pending_queue.py` | Builds the query queue from the atlas MGFs |
| `scripts/run_fastmasst_async.py` + `supervise_fastmasst_async.py` | Checkpointed runner + retry supervisor |
| `scripts/compact_fastmasst_matches.py`, `reconcile_fastmasst_durable_matches.py` | Match compaction + checkpoint reconciliation |
| `scripts/validate_fastmasst_pilots.py` | Pilot validation |
| `scripts/finalize_public_validation_when_ready.py` | Triggers T2 when both queues complete |
| `scripts/public_validation.py` | T2 — Pan-ReDU combined producer |
| `scripts/build_composite_fastmasst_queue.py` | Figure 2c composite-selection queue |

Per-feature match files (~1.4 GB) stay outside git — see
`data_registry/registry.csv`.
