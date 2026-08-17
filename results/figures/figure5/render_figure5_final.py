#!/usr/bin/env python3
"""Render the final Figure 5 proposal from the saved v2 outputs.

Panel a primary: fc-weighted + rules ("lipid-derived community composition",
the estimand matched to the biomass literature ranges and the treatment
analysis). Diamonds: marker-panel estimator ("provenance of matched lipid
signal", where plant necromass legitimately dominates). Panel b unchanged.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
WS = (PROJECT_ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1" / "climgrass"
      / "figure5_redesign_2026-08-08_v2_archlips")
FIG5_V1 = PROJECT_ROOT / "paper2_repro" / "scripts" / "figure5_redesign.py"


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    f5 = load_module("f5v1", FIG5_V1)
    comp = pd.read_csv(WS / "composition_fcweighted_by_sample.csv", index_col=0)
    boot_kingdom = pd.read_csv(WS / "composition_fcweighted_kingdom_ci.csv")
    marker_kingdom = pd.read_csv(WS / "capped_marker_panel_test" / "kingdom_ci_marker_panel.csv")
    effects = pd.read_csv(WS / "fingerprint_set_effects.csv")
    replication = pd.read_csv(WS / "qsip_replication_test.csv")
    unit_groups = dict(zip(effects["unit"], effects["kingdom"]))

    note = ("bars: lipid-derived community composition\n"
            "(fc-weighted, corrected; 95% CI, reference-sample\n"
            "bootstrap); diamonds: provenance of matched lipid\n"
            "signal (specific-marker panel; incl. plant necromass);\n"
            "grey bars: literature biomass ranges. † Archaea from\n"
            "ArchLips-validated ether lipids only (14 diagnostic\n"
            "markers); scale uncertain, no ether-lipid RIE standards")
    f5.render_figure(comp, boot_kingdom, marker_kingdom, effects, replication,
                     unit_groups, "community estimate",
                     WS / "Figure5_final",
                     note_text=note, note_pos=(0.98, 0.05))
    print(f"Wrote {WS / 'Figure5_final.png'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
