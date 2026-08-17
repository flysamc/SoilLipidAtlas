# Figure 1 — Conceptual workflow figure

Data-driven, fully editable vector figure: a Python generator builds the SVG
directly from the analysis release, so every printed number is live, not
hand-typed. This is the **documented exception** to the R-only figure policy —
Figure 1 is a workflow schematic (programmatic vector drawing), not a data
plot.

## Design contract (v2, 2026-08-17)

Figure 1 is a **workflow figure**: it carries the study design and its logic,
not results. Everything that pre-reported later figures was removed —
annotation-tier percentages and validation rates (Fig. 2), the live
composition forest plot (Fig. 5), negative-control statistics (Supp. Fig. 5),
composite weights and alignment tolerances (Supplementary Methods). Panel d's
composition plot is an **explicitly schematic sketch** labelled "see Fig. 5".
The organism count is intentionally absent (the species count is not reliably
derivable; see Table S1 notes), as are per-mode sample counts (Table S1).

Six design-scale numbers remain, all read live from the release or verified
against it: 19 collection phyla (16 in analysis), 6 batches, 273,248 POS /
122,571 NEG consensus features, 11,371 biomarkers across 16 phyla, 736
verified soil fingerprint features.

## Files

| File | Purpose |
|---|---|
| `build_figure1_concept_v1.py` | Generator — builds the SVG programmatically from release CSVs and run summaries |
| `Figure1_concept_skeleton.svg` | The figure: editable vector (1080×720 units; 180×120 mm, Nature double-column) |
| `Figure1_concept_skeleton_preview.png` | 2160×1440 raster preview of the SVG |
| `render_preview.html` | Wrapper used to rasterise the preview (headless Chrome; see comment inside) |
| `LEGEND.md` | The figure legend of record, in the manuscript's legend style |

## The four panels

- **a** — Reference atlas: six organism groups (vector icons), collection
  scale, wet-lab/measurement pipeline, cross-batch alignment into the
  per-mode consensus atlases.
- **b** — Two analytical layers, each phrased as the question it answers:
  the phylum-enriched biomarker atlas (with annotation-tier concept chips)
  and the distributed-fingerprint analysis (similarity → fingerprints →
  independent-method confirmation).
- **c** — Soil decoding: ClimGrass 2×2 design, MS/MS spectral matching
  sketch, the 736-feature verified soil substrate, decomposition.
- **d** — From lipid signal to community composition: plain-language
  correction steps (Calibrate / Share / Archaea), validation pointer
  (Supp. Fig. 5), schematic output sketch (Fig. 5).

## Inputs read live by the generator

| Release file | Values used |
|---|---|
| `biomarker_discovery/summary.json` | 11,371 biomarkers, 16 phyla |
| `climgrass/strict16_archlips_extended_*/RUN_SUMMARY.json` | 736-feature soil substrate (722 SIMPER + 14 ArchLips) |
| `figure3/substrate/substrate_summary.json` | core sample counts (kept as an internal cross-check; not printed) |

## Regenerating

Run the generator inside the analysis archive (it resolves the release by
relative path), then re-copy the SVG here with the checksum sync and re-raster
the preview:

```bash
python results/sync_results.py --apply
```

```bash
chrome --headless --disable-gpu --window-size=2160,1440 --screenshot=Figure1_concept_skeleton_preview.png render_preview.html
```

## Notes

- Fully editable in any vector editor (Illustrator, Inkscape, Affinity) or a
  text editor; each panel and each organism icon is an independent SVG group.
- Group colour palette: Bacteria `#7B52AB`, Archaea `#1B9E8F`, Fungi
  `#D9A420`, Viridiplantae `#3E9C35`, Animalia `#D64541`, Protists `#3B6FD4`
  (locked `ecological_group` display labels; some internal CSV keys remain
  Plantae/Protozoa).
- The two hand-drawn 2026-06 drafts formerly kept here as a visual record are
  superseded and removed; they remain in the analysis archive's recovery
  bundle.
