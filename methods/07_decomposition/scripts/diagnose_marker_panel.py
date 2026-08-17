#!/usr/bin/env python3
"""What carries each phylum's marker-panel score on the extended substrate?

For every phylum: number of phylum-specific enriched markers (k=1 within the
substrate), their summed soil overlap, the top marker's share of that overlap,
and its annotation — to test whether the high Plantae/Protozoa fractions are
distributed marker evidence or single-feature artifacts.
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


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    v2 = load_module("f5v2", V2)
    f5 = load_module("f5v1", FIG5_V1)
    inputs = v2.build_inputs_extended(f5)

    atlas = inputs["atlas"].to_numpy(dtype=float)
    soil = inputs["soil"].to_numpy(dtype=float)
    phyla = inputs["phyla"]
    masks = inputs["enriched_masks"]
    n_enriched = inputs["n_enriched"]
    feature_ids = inputs["feature_ids"]
    sample_phylum = inputs["sample_phylum"]
    cols = list(inputs["atlas"].columns)

    # phylum mean reference (same as build_reference, all samples)
    ref = np.zeros((len(phyla), atlas.shape[0]))
    for i, p in enumerate(phyla):
        idx = [cols.index(s) for s, pp in sample_phylum.items() if pp == p]
        ref[i] = atlas[:, idx].mean(axis=1)

    # annotations
    annot = pd.read_csv(
        PROJECT_ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1" / "climgrass"
        / "strict16_archlips_extended_2026-08-08_v1" / "prepared_inputs" / "annotations_union.csv",
        low_memory=False).set_index("feature_id")

    specific = n_enriched == 1
    rows = []
    total_scores = np.zeros(len(phyla))
    per_phylum_shares = {}
    for i, p in enumerate(phyla):
        m = masks[p] & specific
        if m.sum() == 0:
            continue
        overlap = np.zeros(int(m.sum()))
        for s in range(soil.shape[1]):
            overlap += np.minimum(soil[m, s], ref[i][m])
        total = overlap.sum()
        total_scores[i] = total
        order = np.argsort(overlap)[::-1]
        m_ids = np.array(feature_ids)[m]
        top_id = m_ids[order[0]]
        cls = annot["ls_ClassKey"].get(top_id) if top_id in annot.index else None
        lipid = annot["ls_LipidIon"].get(top_id) if (top_id in annot.index and "ls_LipidIon" in annot.columns) else None
        shares = overlap / max(total, 1e-12)
        per_phylum_shares[p] = np.sort(shares)[::-1][:3]
        rows.append({
            "phylum": p, "kingdom": inputs["unit_groups"][p],
            "n_specific_markers": int(m.sum()),
            "score_share_pct": 0.0,  # filled below
            "top_marker": top_id,
            "top_marker_share_pct": round(float(shares[order[0]]) * 100, 1),
            "top3_share_pct": round(float(np.sort(shares)[::-1][:3].sum()) * 100, 1),
            "top_annotation": (str(lipid)[:30] if pd.notna(lipid) else
                               (str(cls) if pd.notna(cls) else "unannotated")),
        })
    grand = total_scores.sum()
    frame = pd.DataFrame(rows)
    frame["score_share_pct"] = [round(float(total_scores[phyla.index(p)] / grand) * 100, 2)
                                for p in frame["phylum"]]
    frame = frame.sort_values("score_share_pct", ascending=False)
    print(frame.to_string(index=False))
    print("\n(score_share_pct approximates the marker-panel composition summed over the 12 soils)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
