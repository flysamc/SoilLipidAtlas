# Step 05 — Public validation (fastMASST + Pan-ReDU)

Searches every atlas biomarker against the world's public MS/MS data.
Final summaries in `results/`; per-feature match sets (GBs) in the data
registry.

## What this step does

1. **T1 — fastMASST**: each biomarker spectrum is searched against the GNPS2
   Pan-Repo index (`metabolomicspanrepo_index_latest`, cosine ≥ 0.7).
   13,103 queries total: 7,408 POS + 5,695 NEG.
2. **T2 — Pan-ReDU**: maps those matches to sample-type metadata (which
   kingdoms/sample types contain the biomarkers, incl. 149 soil datasets).
   Fail-closed: it refuses to run until **both** T1 queues are 100 % complete
   with exact checksums and denominators.

## Results (database `metabolomicspanrepo_index_latest`)

| Mode | Strict features | Queryable | Any public match | Any soil match |
|---|---|---|---|---|
| POS | 7,751 | 7,406 | 3,978 (54 %) | ✅ **2,287 (31 %)** |
| NEG | 5,697 | 5,695 | 2,283 (40 %) | ✅ **486 (8.5 %)** |

Per-kingdom and per-phylum breakdowns: `results/panredu_summary.csv`
(NEG rows carry the producer's legacy Plantae/Protozoa keys; display labels
are Viridiplantae/Protists per the locked taxonomy policy). Producer +
output checksums: `results/PUBLIC_VALIDATION_MANIFEST.json`.

A second, composite-selection queue (`scripts/build_composite_fastmasst_queue.py`)
validated the alternative Figure-2c selection: IndVal 7,254 features / 30.9 %
soil-detected vs composite 3,613 / 14.6 % — the numbers printed in
Figure 2b/c.

## Data-integrity note: the negative-mode charge convention

> ⚠ The FASST API expects the precursor charge **magnitude**; a signed
> `precursor_charge: -1` silently matches nothing. A first NEG run completed
> its 5,695/5,695 queue gate with every result file empty for exactly this
> reason. The fault was proven with known-positive features (features with
> 12/17/20 historical matches returned 0 with charge −1 and exactly their
> historical matches with +1). Positive mode is unaffected (+1 either way).

The released NEG match set is the **corrected v2 run**
(`fastmasst_async_neg_v2_20260808`), which completed the full 5,695/5,695
queue with verified real content; the void first run is quarantined in the
analysis archive as the fault record. The wrapper's charge handling is fixed
and commented at its charge-parsing block. Consequence for reuse: bulk-search
completeness gates must verify result *content*, not only queue status — an
all-empty result set passes format checks.

## Engineering: how the run survives outages

- **Per-feature checkpointing**: every query writes its own status file and a
  compact gzip match file (fields: USI, Dataset, Cosine, Matching Peaks,
  Delta Mass; raw response hashed then discarded). A crash loses nothing.
- **Supervisors** retry durable errors in rounds until the completeness gate
  (all expected successes present) passes. The POS queue survived a multi-day
  GNPS2 backend outage mid-run this way.
- Full run provenance in `run_manifest_POS.json` / `run_manifest_NEG.json`:
  input MGF checksums, every parameter, wrapper + recovered-producer hashes.

## Files

| File | What it is |
|---|---|
| `run_manifest_POS.json`, `run_manifest_NEG.json` | Complete provenance of the query queues |
| `scripts/build_fastmasst_pending_queue.py` | Builds the query queue from the atlas MGFs |
| `scripts/run_fastmasst_async.py` + `supervise_fastmasst_async.py` | The checkpointed runner + its retry supervisor |
| `scripts/compact_fastmasst_matches.py`, `reconcile_fastmasst_durable_matches.py` | Match compaction + checkpoint reconciliation |
| `scripts/validate_fastmasst_pilots.py` | Pilot validation gates |
| `scripts/finalize_public_validation_when_ready.py` | Watcher that triggers T2 when both queues complete |
| `scripts/public_validation.py` | T2 — the fail-closed Pan-ReDU combined producer |

The per-feature match files (~1.4 GB POS + NEG v2) stay outside git — see
`data_registry/registry.csv`.

Completion path as executed: both denominators verified (POS queryable
7,406 of 7,751 strict; NEG 5,695), match sets checksummed, then the
fail-closed T2 producer wrote the two summary CSVs now in `results/`.
