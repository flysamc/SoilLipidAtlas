#!/usr/bin/env python3
"""Test a dominance-capped marker panel on the extended substrate.

Motivation: the marker-quality diagnostic showed that some phyla's specific
markers are concentrated in one chemically generic feature (Mollusca: one LPC
at 50.6% of its evidence; Discosea: 48.6%), while others are distributed
(Streptophyta: 92 markers, top 29%). Rule: no single marker may contribute more
than CAP of a phylum's summed marker overlap (per sample, waterfill-free simple
clip + renormalise within the overlap vector). Direction-blind, all phyla.

Runs the capped variant through the SAME held-out mixture benchmark and reports
its real-soil composition next to the uncapped winner and fc-weighted.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
V2 = PROJECT_ROOT / "paper2_repro" / "scripts" / "figure5_redesign_v2.py"
FIG5_V1 = PROJECT_ROOT / "paper2_repro" / "scripts" / "figure5_redesign.py"
OUT = (PROJECT_ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1" / "climgrass"
       / "figure5_redesign_2026-08-08_v2_archlips" / "capped_marker_panel_test")
CAP = 0.25


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_capped(cap=CAP):
    def est_marker_capped(ref, y, ctx):
        phyla = ctx["phyla"]
        masks = ctx["enriched_masks"]
        specific = ctx["n_enriched"] == 1
        scores = np.zeros(len(phyla))
        for i, p in enumerate(phyla):
            m = masks[p] & specific
            if m.sum() == 0:
                continue
            overlap = np.minimum(y[m], ref[i][m])
            total = overlap.sum()
            if total <= 0:
                continue
            # clip each marker's contribution to cap * (renormalised total),
            # iterated a few times so the cap holds after renormalisation
            ov = overlap.copy()
            for _ in range(6):
                total = ov.sum()
                if total <= 0:
                    break
                over = ov > cap * total
                if not over.any():
                    break
                ov[over] = cap * total
            scores[i] = ov.sum()
        return scores / scores.sum() if scores.sum() > 0 else scores
    return est_marker_capped


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    v2 = load_module("f5v2", V2)
    f5 = load_module("f5v1", FIG5_V1)
    inputs = v2.build_inputs_extended(f5)

    f5.ESTIMATORS["marker_capped"] = make_capped()

    print("Benchmark including capped marker panel...")
    bench, summary = f5.mixture_benchmark(inputs)
    summary.to_csv(OUT / "benchmark_summary_with_capped.csv", index=False)
    print(summary.round(4).to_string(index=False))

    print("\nReal-soil kingdom means (%):")
    rows = {}
    for est, sqrt_t in (("marker_panel", False), ("marker_capped", False),
                        ("excess_bc", False), ("fc_weighted_bc", False)):
        comp, _, king = f5.soil_composition(inputs, est, sqrt_t)
        rows[est] = (king.set_index("kingdom")["mean"] * 100).round(1)
        king.to_csv(OUT / f"kingdom_ci_{est}.csv", index=False)
    table = pd.DataFrame(rows)
    table.to_csv(OUT / "kingdom_means_comparison.csv")
    print(table.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
