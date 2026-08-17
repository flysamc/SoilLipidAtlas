# Annotation Step 8 (ArchLips) — release `ncbi-phylum-2026-08-04-v1`

Supplementary Method 3 Step 8 / Supplementary Method 7, targeted search of
archaeal-assigned features. Producer: `scripts/step8_archlips_release.py`.
Parameters and library checksums: `STEP8_MANIFEST.json`.

## The two polarities are different procedures

The ArchLips libraries hold **positive-mode reference spectra only**:
29,446 `[M+H]+` and 29,446 `[M+NH4]+` in the high-confidence library, zero
negative entries. Consequently:

- **POS** — direct spectral matching (`archlips_indval_search.py`): precursor
  ±10 ppm, modified cosine at 0.02 Da, requiring cosine ≥ 0.3 and ≥ 2 matched
  peaks. Gold cos ≥ 0.7 & ≥ 4 peaks; Silver cos ≥ 0.5 & ≥ 3; Bronze cos ≥ 0.3
  & ≥ 2.
- **NEG** — no cosine (`archlips_search_neg.py`): negative adducts are
  computed from the neutral MW (`[M-H]⁻` −1.00728, `[M+HCOO]⁻` +44.99820,
  `[M+CH₃COO]⁻` +59.01385), matched at ±10 ppm, then checked against 12
  archaeal diagnostic ions and 6 neutral losses. Cosine against an `[M+H]+`
  reference would be meaningless.

Libraries were streamed (1.17 GB + 154 MB); a compound was retained only if
it could match some query precursor. 9,531 of 348,430 entries retained.

## Results

| Mode | Archaeal strict features | With match | Validated | % | Gold+Silver | % |
|---|---|---|---|---|---|---|
| POS | 1,333 | 632 | 632 | 47.41 | 406 | **30.46** |
| NEG | 21 | 17 | 10 | 47.62 | 0 | **0.00** |

After the step-5 RT screen (86 impossible-RT calls removed), POS
Gold+Silver is **27.6 %**. That screened value is the one used downstream.

POS tiers by phylum (before the RT screen):

| Phylum | Gold | Silver | Bronze |
|---|---|---|---|
| Methanobacteriota | 225 | 106 | 178 |
| Thermoproteota | 46 | 29 | 48 |

6 of the 1,333 POS archaeal features have no usable MS2 and cannot be matched.

## POS matches are chemically coherent

Cosine median 0.867 (max 1.000), maximum mass error 9.95 ppm, reference
adducts `[M+H]+` 384 and `[M+NH4]+` 248. Frequent hits are canonical
archaeal ether lipids (archaeol and related). LipidSearch (step 1) barely
covers these classes; ArchLips is the archaeal annotation source in this
analysis.

## Negative mode has no MS2 support

All 21 NEG archaeal spectra are of good quality (12–160 peaks, median 49).
Against 12 archaeal diagnostic ions and 6 neutral losses:

```
total diagnostic ions across all 21 spectra: 0
total neutral losses across all 21 spectra:  0
```

The 10 Bronze calls (102 rows) rest on `|error| ≤ 5 ppm` and
`detection_rate ≥ 0.3` only. Derived columns make that explicit:

- `spectral_support` — `False` for all 102
- `evidence_basis` — `mass_and_detection_rate_only` for all 102

The NEG mass match is also unspecific: 918 candidate compounds match 17
features (~54 candidates per feature). Negative-mode archaeal ArchLips
annotation is therefore reported as mass-accuracy-only, not as a validated
identification.

This step is targeted (archaeal-assigned features only).

## Files

| File | Contents |
|---|---|
| `archlips_pos_strict_matches.csv` | 632 POS matches with cosine, matched peaks, tier, phylum |
| `archlips_neg_strict_matches.csv` | 918 NEG mass matches with diagnostic counts, tier, `spectral_support`, `evidence_basis` |
| `step8_coverage_summary.csv` | The results table above |
| `STEP8_MANIFEST.json` | Parameters, library and producer sha256 |
