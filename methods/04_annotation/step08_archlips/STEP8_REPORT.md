# Annotation Step 8 (ArchLips) — release `ncbi-phylum-2026-08-04-v1`

Supplementary Method 3 Step 8 / Supplementary Method 7, Round 1 (targeted).
Producer: `paper2_repro/scripts/step8_archlips_release.py`.
Parameters, checksums and library provenance: `STEP8_MANIFEST.json`.

## The two polarities are different procedures

The ArchLips libraries hold **positive-mode reference spectra only** — verified:
29,446 `[M+H]+` and 29,446 `[M+NH4]+` in the high-confidence library, zero negative
entries. Consequently:

- **POS** — direct spectral matching, per `archlips_indval_search.py`: precursor
  ±10 ppm, modified cosine at 0.02 Da, requiring cosine ≥ 0.3 and ≥ 2 matched peaks.
  Gold cos ≥ 0.7 & ≥ 4 peaks; Silver cos ≥ 0.5 & ≥ 3; Bronze cos ≥ 0.3 & ≥ 2.
- **NEG** — no cosine, per `archlips_search_neg.py`: negative adducts are *computed*
  from the neutral MW (`[M-H]⁻` −1.00728, `[M+HCOO]⁻` +44.99820, `[M+CH₃COO]⁻`
  +59.01385), matched at ±10 ppm, then validated against 12 archaeal diagnostic ions
  and 6 neutral losses. Cosine against an `[M+H]+` reference would be meaningless.

Both libraries were streamed rather than loaded (1.17 GB + 154 MB); a compound was
retained only if it could match some query precursor. 9,531 of 348,430 entries retained.

## Results

| Mode | Archaeal strict features | With match | Validated | % | Gold+Silver | % |
|---|---|---|---|---|---|---|
| POS | 1,333 | 632 | 632 | 47.41 | 406 | **30.46** |
| NEG | 21 | 17 | 10 | 47.62 | 0 | **0.00** |

POS tiers by phylum:

| Phylum | Gold | Silver | Bronze |
|---|---|---|---|
| Methanobacteriota | 225 | 106 | 178 |
| Thermoproteota | 46 | 29 | 48 |

6 of the 1,333 POS archaeal features have no usable MS2 and can never be matched.

## POS validation

Chemically coherent. Cosine median 0.867 (max 1.000), maximum mass error 9.95 ppm,
reference adducts `[M+H]+` 384 and `[M+NH4]+` 248. The most frequently matched
compounds are canonical archaeal ether lipids:

```
Archaeol(20:0_20:1)  17    PA-Archaeol(20:0_20:1)       7
Archaeol(20:0_20:0)  15    MP(5:2(OH))                  6
Archaeol(20:1_20:1)  10    Archaeol(20:0(OH)_20:0)      6
phytoene              7    Gly-Archaeol(20:0(OH)_20:0)  6
```

This is what confirms Step 8 is doing its job: Step 1 (LipidSearch) reached only
1.38% Grade A/B for Methanobacteriota and 0.32% for Thermoproteota, because
LipidSearch has no archaeal ether lipids. Step 8 lifts POS archaeal Gold+Silver to
30.46%, roughly a twentyfold improvement.

For context, Supplementary Method 7 reports the historical Round 1 improving
archaeal Gold+Silver from 3.9% to 21.3% over 895 features. This release reaches
30.46% over 1,333 features. The numbers are not comparable without re-deriving the
historical denominator, so **treat 30.46% as a new release number, not as a
confirmation of 21.3%**.

## Negative mode has no spectral support — never presented as validated

This is the step's most consequential finding.

All 21 NEG archaeal spectra are of good quality (12–160 peaks, median 49). Scanning
them against the 12 archaeal diagnostic ions and 6 neutral losses returns:

```
total diagnostic ions across all 21 spectra: 0
total neutral losses across all 21 spectra:  0
```

Yet 102 rows across 10 features are labelled Bronze. Those calls come entirely from
the fallback rung of the recovered tier ladder, `|error| ≤ 5 ppm and detection_rate
≥ 0.3` — and every NEG archaeal feature has a detection rate of 0.412–0.500, so that
rung fires on mass accuracy alone. The tier column is left exactly as the recovered
producer assigns it; two derived columns make the basis explicit:

- `spectral_support` — `False` for all 102
- `evidence_basis` — `mass_and_detection_rate_only` for all 102

Compounding this, the NEG mass match is not specific: 918 candidate compounds match
just 17 features, about **54 candidates per feature**.

**Conclusion: negative-mode archaeal ArchLips annotation in this release has no MS2
evidence.** It is reported as mass-accuracy-only, never as validated archaeal
identification. With 0 Gold+Silver across 21 features, the manuscript makes no
archaeal negative-mode annotation claim from this step.

## Files

| File | Contents |
|---|---|
| `archlips_pos_strict_matches.csv` | 632 POS matches with cosine, matched peaks, tier, phylum |
| `archlips_neg_strict_matches.csv` | 918 NEG mass matches with diagnostic counts, tier, `spectral_support`, `evidence_basis` |
| `step8_coverage_summary.csv` | The results table above |
| `STEP8_MANIFEST.json` | Parameters, library and producer sha256, polarity note |

## Scope note

This release runs Round 1 (targeted, archaeal-assigned features only).
Supplementary Method 7 Round 2 is a full-batch untargeted search across all
six FBMN batches (historically 42,983 raw mass matches narrowed to 7,328
archaeal-enriched candidates). Round 2 was not rerun for this release, so the
historical Round-2 impact figures in Supplementary Method 3 are not restated.
