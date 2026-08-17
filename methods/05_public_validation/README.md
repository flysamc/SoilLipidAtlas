# Step 05 — Public validation (fastMASST + Pan-ReDU)

**Status: ✅ complete (2026-08-11), both polarities, one current index.**
Final summaries in `results/`; per-feature match sets (GBs) in the data
registry.

## What this step does

Searches every atlas biomarker against the world's public MS/MS data:

1. **T1 — fastMASST**: each biomarker spectrum is searched against the GNPS2
   Pan-Repo index (`metabolomicspanrepo_index_latest`, cosine ≥ 0.7).
   13,103 queries total: 7,408 POS + 5,695 NEG.
2. **T2 — Pan-ReDU**: maps those matches to sample-type metadata (which
   kingdoms/sample types contain our biomarkers, incl. 149 soil datasets).
   Fail-closed: it refuses to run until **both** T1 queues are 100% complete
   with exact checksums and denominators.

## Results (T2 run 2026-08-11, database `metabolomicspanrepo_index_latest`)

| Mode | Strict features | Queryable | Any public match | Any soil match |
|---|---|---|---|---|
| POS | 7,751 | 7,406 | 3,978 (54%) | 2,287 (31%) |
| NEG | 5,697 | 5,695 | 2,283 (40%) | 486 (8.5%) |

Per-kingdom and per-phylum breakdowns: `results/panredu_summary.csv`
(NEG rows carry the producer's legacy Plantae/Protozoa keys; display labels
are Viridiplantae/Protists per the locked taxonomy policy). Producer +
output checksums: `results/PUBLIC_VALIDATION_MANIFEST.json`.

A second, composite-selection queue (`scripts/build_composite_fastmasst_queue.py`,
run `fastmasst_composite_20260811`) validated the alternative Figure-2c
selection: IndVal 7,254 features / 30.9% soil-detected vs composite 3,613 /
14.6% — the number printed in Figure 2b/c.

## History: the outage and the signed-charge bug

Paused 2026-08-06 when the GNPS2 backend died mid-run; resumed 2026-08-08
after verified recovery.

⚠ **NEG restarted from zero the same day.** The first NEG run completed its
5,695/5,695 gate but every result file was empty: the wrapper sent
`precursor_charge: -1`, and the FASST API expects the charge **magnitude** —
negative values silently match nothing. Proven with known-positive features
(12/17/20 historical matches: 0 with −1, exact historical reproduction with
+1). POS is unaffected (+1 either way). The wrapper is fixed (see the comment
at its charge-parsing block), the void run is quarantined in the analysis
archive as the bug
record, and the corrected NEG run (`fastmasst_async_neg_v2_20260808`)
completed its full 5,695/5,695 gate with real content. Lesson recorded:
verify result *content* after bulk runs — an all-empty set passes format
checks.

## How it survives outages (the engineering)

- **Per-feature checkpointing**: every query writes its own status file and a
  compact gzip match file (fields: USI, Dataset, Cosine, Matching Peaks,
  Delta Mass; raw response hashed then discarded). A crash loses nothing.
- **Supervisors** retry durable errors in rounds until the completeness gate
  (all expected successes present) passes.
- Full run provenance in `run_manifest_POS.json` / `run_manifest_NEG.json`:
  input MGF checksums, every parameter, wrapper + recovered-producer hashes.

## Files

| File | What it is |
|---|---|
| `run_manifest_POS.json`, `run_manifest_NEG.json` | Complete provenance of the running queues |
| `scripts/build_fastmasst_pending_queue.py` | Builds the query queue from the atlas MGFs |
| `scripts/run_fastmasst_async.py` + `supervise_fastmasst_async.py` | The checkpointed runner + its retry supervisor |
| `scripts/compact_fastmasst_matches.py`, `reconcile_fastmasst_durable_matches.py` | Match compaction + checkpoint reconciliation (used after the outage) |
| `scripts/validate_fastmasst_pilots.py` | Pilot validation gates |
| `scripts/finalize_public_validation_when_ready.py` | Watcher that triggers T2 when both queues complete |
| `scripts/public_validation.py` | T2 — the fail-closed Pan-ReDU combined producer |

The per-feature match files (~1.4 GB POS + NEG v2) stay outside git — see
`data_registry/registry.csv`.

Completion path as executed: both denominators verified (POS queryable
7,406 of 7,751 strict; NEG 5,695), match sets checksummed, then the
fail-closed T2 producer wrote the two summary CSVs now in `results/`.
