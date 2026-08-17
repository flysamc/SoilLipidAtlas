#!/usr/bin/env python3
"""Screen the full ClimGrass soil quant table (all features, not just the
722-feature SIMPER substrate) for ArchLips-annotated archaeal lipid masses.

Match criteria: precursor m/z within PPM tolerance; RT compared where both
sides carry it (atlas and soil share the same LC method family, but RT is
reported descriptively, not used as a gate). Also screens the 2,313 corrected
soil-to-atlas spectral matches for atlas m/z values that coincide with
ArchLips features.
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(r"C:\Users\Shadow\Desktop\P2R")
HANDOFF = ROOT / "cg" / "SHADOW_PC_CLIMGRASS_STRICT_PHYLA_HANDOFF_2026-08-06"
PREPARED = (ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1" / "climgrass"
            / "strict16_all_simper_sensitivity_2026-08-06" / "prepared_inputs")

PPM = 5.0


def ppm_match(soil_mz: np.ndarray, targets: pd.DataFrame, mz_col: str, ppm: float):
    """Return list of (soil_idx, target_idx) pairs within ppm."""
    t_mz = targets[mz_col].astype(float).to_numpy()
    order = np.argsort(t_mz)
    t_sorted = t_mz[order]
    pairs = []
    for i, mz in enumerate(soil_mz):
        tol = mz * ppm / 1e6
        lo = np.searchsorted(t_sorted, mz - tol, "left")
        hi = np.searchsorted(t_sorted, mz + tol, "right")
        for j in range(lo, hi):
            pairs.append((i, order[j]))
    return pairs


def main():
    quant = pd.read_csv(HANDOFF / "data_contract" / "soil-data-without-background_quant.csv",
                        low_memory=False)
    area_cols = [c for c in quant.columns if "Peak area" in c]
    soil_mz = quant["row m/z"].astype(float).to_numpy()
    soil_rt = quant["row retention time"].astype(float).to_numpy()
    soil_max_area = quant[area_cols].astype(float).max(axis=1).to_numpy()
    soil_n_detected = (quant[area_cols].astype(float) > 0).sum(axis=1).to_numpy()

    lists = {
        "release_rt_screened": (
            pd.read_csv(PREPARED / "archlips_pos_release_eligible_rt_screened.csv", low_memory=False),
            "query_mz", "archlips_name", None,
        ),
        "combined": (
            pd.read_csv(HANDOFF / "data_contract" / "archaeal_identification_combined.csv", low_memory=False),
            "consensus_mz", "archlips_name", "consensus_rt",
        ),
        "fullbatch_ms2_validated": (
            pd.read_csv(HANDOFF / "data_contract" / "archlips_fullbatch_ms2_validated.csv", low_memory=False),
            "precursor_mz", "archlips_name", "rt_minutes",
        ),
    }

    all_hits = []
    for list_name, (table, mz_col, name_col, rt_col) in lists.items():
        table = table.dropna(subset=[mz_col]).reset_index(drop=True)
        pairs = ppm_match(soil_mz, table, mz_col, PPM)
        for soil_idx, t_idx in pairs:
            row = table.iloc[t_idx]
            rt_ref = float(row[rt_col]) if rt_col and pd.notna(row.get(rt_col)) else np.nan
            hit = {
                "list": list_name,
                "soil_row_id": int(quant.iloc[soil_idx]["row ID"]),
                "soil_mz": round(float(soil_mz[soil_idx]), 5),
                "soil_rt": round(float(soil_rt[soil_idx]), 2),
                "ref_mz": round(float(row[mz_col]), 5),
                "ppm": round(abs(soil_mz[soil_idx] - row[mz_col]) / row[mz_col] * 1e6, 2),
                "ref_rt": round(rt_ref, 2) if np.isfinite(rt_ref) else None,
                "rt_delta": round(abs(soil_rt[soil_idx] - rt_ref), 2) if np.isfinite(rt_ref) else None,
                "archlips_name": str(row.get(name_col))[:45],
                "ref_feature_id": row.get("feature_id", row.get("aligned_feature_id")),
                "soil_n_samples_detected": int(soil_n_detected[soil_idx]),
                "soil_max_area": float(soil_max_area[soil_idx]),
            }
            for extra in ("archlips_confidence", "overall_confidence", "confidence",
                          "archlips_tier", "biomarker_quality", "exclusive_to_archaea",
                          "cosine", "cosine_score", "phylum"):
                if extra in row.index and pd.notna(row[extra]):
                    hit[extra] = row[extra]
            all_hits.append(hit)

    hits = pd.DataFrame(all_hits)
    print(f"Soil quant rows: {len(quant)}; PPM tolerance: {PPM}")
    if hits.empty:
        print("NO precursor-level ArchLips matches in the soil quant table.")
    else:
        print(f"\nTotal precursor-level hits: {len(hits)} "
              f"({hits['soil_row_id'].nunique()} distinct soil features)")
        print(hits.groupby("list")["soil_row_id"].nunique().to_string())
        out = ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1" / "climgrass" / "archlips_soil_precursor_screen_5ppm.csv"
        hits.sort_values(["list", "ppm"]).to_csv(out, index=False)
        print(f"\nSaved: {out}")
        # Show most credible: distinct soil features, best ppm per feature, prefer RT-consistent
        best = hits.sort_values("ppm").drop_duplicates(subset=["soil_row_id"])
        best = best.sort_values(["rt_delta"], na_position="last")
        cols = ["list", "soil_row_id", "soil_mz", "soil_rt", "ref_rt", "rt_delta", "ppm",
                "archlips_name", "soil_n_samples_detected", "soil_max_area"]
        print("\nDistinct soil features (best hit each, sorted by RT agreement):")
        with pd.option_context("display.width", 250):
            print(best[cols].head(40).to_string(index=False))

    # Spectral-level: corrected match atlas m/z vs ArchLips masses
    matches = pd.read_csv(HANDOFF / "results" / "strict_19unit_corrected_2026-08-06" / "corrected_spectral_matches.csv")
    amz = matches["atlas_mz"].astype(float).to_numpy()
    print("\n--- Spectral-level: 2,313 corrected soil-to-atlas matches vs ArchLips masses ---")
    for list_name, (table, mz_col, name_col, rt_col) in lists.items():
        table = table.dropna(subset=[mz_col]).reset_index(drop=True)
        pairs = ppm_match(amz, table, mz_col, PPM)
        n = len({p[0] for p in pairs})
        print(f"  {list_name}: {n} of 2,313 matched atlas spectra coincide with an ArchLips mass at {PPM} ppm")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())
