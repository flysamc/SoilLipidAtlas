# Step 07 — Source decomposition + negative control

The decomposition framework estimates what fraction of a mixed lipidome comes
from each organism group, using the reference fingerprints from step 06. Two
estimators are carried through the paper:

- **fc-weighted** (primary) — fold-change-weighted contribution of every
  fingerprint feature.
- **marker-panel** (secondary) — restricted to top discriminating markers.

They fail on different groups by design, which is why both are reported.

## Negative control: pure-isolate leave-one-out (`negative_control/`)

Producer `scripts/suppfig5_negative_control_strict16.py` (a **declared
reimplementation** — the submitted producer was never recovered, so this
panel is never claimed to reproduce the published one). Each of **n = 164**
pure isolates is decomposed against phylum centroids rebuilt *without* that
isolate:

| Metric | Uncorrected | Corrected |
|---|---|---|
| Dominant group correct (fc-weighted) | 75.0% | **79.3%** (130/164) |
| Archaeal self-recovery | 68.2% | **93.1%** |
| Non-archaeal samples leaking into Archaea | — | ≤0.8% |

Marker-panel secondary: 83.5% dominant-group correct overall but archaeal
self-recovery drops to 57.9% — the estimator trade-off in one number.

Files: per-sample decompositions (uncorrected / corrected / marker-panel),
group confusion matrices, self-recovery summary, `RUN_SUMMARY.json` (read
live by Figure 1 panel d), and the Supplementary Fig. 5 caption.

Supporting diagnostics in `scripts/`: `test_capped_marker_panel.py`,
`diagnose_marker_panel.py`.
