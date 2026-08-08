# Step 05 — Public validation (fastMASST + Pan-ReDU)

**Status: ⏳ running.** This is the skeleton — scripts and provenance now,
result tables when the searches finish.

## What this step does

Searches every atlas biomarker against the world's public MS/MS data:

1. **T1 — fastMASST**: each biomarker spectrum is searched against the GNPS2
   Pan-Repo index (`metabolomicspanrepo_index_latest`, cosine ≥ 0.7).
   13,103 queries total: 7,408 POS + 5,695 NEG.
2. **T2 — Pan-ReDU**: maps those matches to sample-type metadata (which
   kingdoms/sample types contain our biomarkers, incl. 149 soil datasets).
   Fail-closed: it refuses to run until **both** T1 queues are 100% complete
   with exact checksums and denominators.

## Where it stands (2026-08-08)

Paused 2026-08-06 when the GNPS2 backend died mid-run; resumed 2026-08-08
after verified recovery.

⚠ **NEG restarted from zero the same day.** The first NEG run completed its
5,695/5,695 gate but every result file was empty: the wrapper sent
`precursor_charge: -1`, and the FASST API expects the charge **magnitude** —
negative values silently match nothing. Proven with known-positive features
(12/17/20 historical matches: 0 with −1, exact historical reproduction with
+1). POS is unaffected (+1 either way). The wrapper is fixed (see the comment
at its charge-parsing block), the void run is quarantined in P2R as the bug
record, and the corrected NEG run (`fastmasst_async_neg_v2_20260808`) is in
progress. Lesson recorded: verify result *content* after bulk runs — an
all-empty set passes format checks.

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

The per-feature match files (~1.4 GB POS and growing) stay outside git — see
`data_registry/registry.csv`; final checksums land there at completion.

## When the queues finish

1. Verify both success denominators (7,408 / 5,695) and checksum the match sets.
2. Run `public_validation.py` (T2) with verified paths.
3. Add the combined validation summaries here, update the registry, and the
   step flips to ✅.
