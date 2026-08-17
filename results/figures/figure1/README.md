# Figure 1 — Conceptual study-design schematic

**Status: updated.** Editable Python generator and rendered SVG.

## Files

| File | Purpose |
|---|---|
| `build_figure1_concept_v1.py` | Python script that programmatically builds the SVG from analysis outputs |
| `Figure1_concept_skeleton.svg` | Rendered editable vector (1080×720, 180×120 mm Nature Comms double-column) |

## Description

4-panel conceptual overview:

- **Panel a** — Reference lipidome atlas from 6 organism groups (19 collection phyla, 16 analysed)
- **Panel b** — Two complementary analytical layers (biomarker atlas + distributed fingerprints)
- **Panel c** — Soil community decoding (ClimGrass 2×2 factorial design, MS/MS matching, decomposition)
- **Panel d** — Lipid framework (quantification correction, validation, forest plot of lipid-derived proportions)

## Updating

To regenerate with current data:

```bash
cd C:\Users\Shadow\Desktop\P2R
python paper2_repro\scripts\build_figure1_concept_v1.py
```

Then copy the output SVG back here (checksum-synced):

```bash
python results/sync_results.py --apply
```

## Notes

- The SVG is fully editable in any vector editor (Illustrator, Inkscape, Affinity, or text editor).
- Six group colour palette: Bacteria `#7B52AB`, Archaea `#1B9E8F`, Fungi `#D9A420`, Viridiplantae `#3E9C35`, Animalia `#D64541`, Protists `#3B6FD4` (locked `ecological_group` labels; legacy CSV keys remain Plantae/Protozoa internally).
- Annotation tiers: Gold (molecular species), Silver (partial), Bronze (lipid class), Unidentified.
- All data-driven values read from the ncbi-phylum-2026-08-04-v1 analysis release.
