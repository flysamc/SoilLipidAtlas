# Analysis progress

One row per methodology step. A step appears in `methods/` once its files are
reviewed and added. Details inside each folder's README.

| Step | What it does | Status | In repo? |
|---|---|---|---|
| 01 Taxonomy | Assigns each sample its NCBI phylum; 19 collected → 16 analysed | ✅ done, locked | ✅ |
| 02 Biomarker atlas | Selects phylum-enriched features (11,371 POS / 5,697 NEG) | ✅ done | ☐ |
| 03 Annotation | Identifies what the biomarker lipids are (11-step pipeline) | 🟡 steps 1–5 + 8 done; 6 partial; 7 + 9 waiting on missing resources | ☐ |
| 04 SIMPER fingerprints | Per-phylum lipid fingerprint atlas | ✅ done | ☐ |
| 05 Public validation | Searches biomarkers in public data (fastMASST, Pan-ReDU) | ⏸ ~half done, paused — GNPS2 server down | ☐ |
| 06 ClimGrass | Applies fingerprints to field soil samples | 🚧 in progress — spectral matcher being reworked | ☐ |
| Figures | Built from finished steps above | ☐ not started | ☐ |
| Figure 1 | Being updated by hand (Rahul) | ⏸ | ☐ |
