# Figure 1 — Conceptual study-design schematic

Data-driven, fully editable vector figure: a Python generator builds the SVG
directly from the analysis release, so every printed number is live, not
hand-typed.

## Files

| File | Purpose |
|---|---|
| `build_figure1_concept_v1.py` | Generator — builds the SVG programmatically from release CSVs and run summaries |
| `Figure1_concept_skeleton.svg` | Rendered editable vector (1080×720 px; 180×120 mm, Nature double-column) |
| `versionA_*/versionB_*` renders | Historical hand-drawn drafts (superseded; kept as visual record) |

## The four panels

- **a** — Reference lipidome atlas: organisms from six ecological groups
  (19 collection phyla, 16 analysed), with vector organism icons.
- **b** — Two complementary analytical layers: the biomarker atlas
  (annotation-tier bar) and distributed fingerprints; external validation
  line reports the fastMASST public-data result.
- **c** — Soil community decoding: ClimGrass 2×2 factorial design, MS/MS
  matching, source decomposition.
- **d** — Framework validation: quantification correction and the
  pure-isolate negative control, read live from
  `methods/07_decomposition/negative_control/RUN_SUMMARY.json`.

## Regenerating

Run the generator inside the analysis archive (it reads release CSVs and run
summaries by relative path), then re-copy the SVG here with the checksum
sync:

```bash
python results/sync_results.py --apply
```

## Notes

- The SVG is directly editable in any vector editor (Illustrator, Inkscape,
  Affinity) or a text editor.
- Group colour palette: Bacteria `#7B52AB`, Archaea `#1B9E8F`, Fungi
  `#D9A420`, Viridiplantae `#3E9C35`, Animalia `#D64541`, Protists `#3B6FD4`
  (locked `ecological_group` display labels; some internal CSV keys remain
  Plantae/Protozoa).
- Annotation tiers: Gold (molecular species) > Silver (partial structure) >
  Bronze (lipid class) > Unidentified.
- All data-driven values come from release `ncbi-phylum-2026-08-04-v1`.
