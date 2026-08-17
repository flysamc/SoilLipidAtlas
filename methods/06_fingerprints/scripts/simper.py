#!/usr/bin/env python3
"""Reproduce the published centroid-SIMPER method, then run final NCBI phyla."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

TOP_N = 3000
MIN_BATCHES = 2
MIN_REPRODUCTION_R = 0.995
MAX_FROZEN_N50_RELATIVE_DELTA = 0.05

LEGACY_UNITS = sorted([
    "Actinomycetota", "Amoebozoa", "Arthropoda", "Ascomycota", "Bacillota",
    "Basidiomycota", "Bryophyta", "Chlorophyta", "Euryarchaeota",
    "Marchantiophyta", "Methanobacteriota", "Mollusca", "Mucoromycota",
    "Nematoda", "Pseudomonadota", "Thermoproteota", "Trachaeophyta", "Virus",
])


def gini(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[values >= 0]
    if values.size == 0 or values.sum() == 0:
        return 0.0
    ordered = np.sort(values)
    n = ordered.size
    return float((2 * np.dot(np.arange(1, n + 1), ordered) / (n * ordered.sum())) - (n + 1) / n)


def mapped_matrix(table: pd.DataFrame, metadata: pd.DataFrame, units: list[str]) -> tuple[np.ndarray, np.ndarray, list[str]]:
    sample_columns = [column for column in table.columns if str(column).startswith("sample:")]
    sample_to_unit = {
        row["original_header"]: row["phylum"]
        for _, row in metadata.iterrows()
        if row.get("original_header") in table.columns and row.get("phylum") in units
    }
    mapped = [column for column in sample_columns if column in sample_to_unit]
    labels = np.asarray([sample_to_unit[column] for column in mapped])
    values = table[mapped].fillna(0.0).to_numpy(dtype=np.float64).T
    return values, labels, mapped


def build(
    values: np.ndarray,
    labels: np.ndarray,
    feature_ids: np.ndarray,
    feature_meta: pd.DataFrame,
    units: list[str],
    groups: dict[str, str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    atlas_parts: list[pd.DataFrame] = []
    summary: list[dict] = []
    for unit in units:
        target = labels == unit
        if target.sum() == 0:
            continue
        mean_target = values[target].mean(axis=0)
        mean_background = values[~target].mean(axis=0)
        contribution = np.abs(mean_target - mean_background)
        total = contribution.sum()
        contribution_pct = contribution / total * 100 if total > 0 else np.zeros_like(contribution)
        order = np.argsort(-contribution_pct)
        cumulative_full = np.cumsum(contribution_pct[order])
        thresholds = {
            pct: int(np.searchsorted(cumulative_full, pct, side="left") + 1)
            for pct in (50, 80, 90, 95)
        }
        keep = order[: min(TOP_N, len(order))]
        atlas_parts.append(pd.DataFrame({
            "feature_id": feature_ids[keep],
            "phylum": unit,
            "ecological_group": groups.get(unit, ""),
            "consensus_mz": feature_meta["consensus_mz"].to_numpy()[keep],
            "consensus_rt": feature_meta["consensus_rt"].to_numpy()[keep],
            "n_batches": feature_meta["n_batches"].to_numpy()[keep],
            "simper_rank": np.arange(1, len(keep) + 1),
            "contribution_pct": contribution_pct[keep],
            "cumulative_pct": np.cumsum(contribution_pct[keep]),
            "direction": np.where(mean_target[keep] >= mean_background[keep], "enriched", "depleted"),
            "mean_target": mean_target[keep],
            "mean_background": mean_background[keep],
            "fold_change": (mean_target[keep] + 1e-9) / (mean_background[keep] + 1e-9),
        }))
        summary.append({
            "phylum": unit,
            "ecological_group": groups.get(unit, ""),
            "n_samples": int(target.sum()),
            "n_features": int(values.shape[1]),
            "n_50pct": thresholds[50],
            "n_80pct": thresholds[80],
            "n_90pct": thresholds[90],
            "n_95pct": thresholds[95],
            "max_single_pct": float(contribution_pct[order[0]]),
            "gini": gini(contribution_pct),
        })
    return pd.concat(atlas_parts, ignore_index=True), pd.DataFrame(summary)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--taxonomy-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--frozen-curves",
        type=Path,
        default=None,
        help=(
            "Frozen SIMPER curves used for the bounded historical n50 gate when "
            "the historical per-feature atlas is unavailable."
        ),
    )
    args = parser.parse_args()

    table_path = args.source_root / "analysis/analysis-15/03_alignment/consensus_aligned_table.csv"
    table = pd.read_csv(table_path, index_col=0, low_memory=False)
    table = table.loc[table["n_batches"] >= MIN_BATCHES].copy()
    feature_meta = table[["consensus_mz", "consensus_rt", "n_batches"]].copy()
    feature_ids = np.asarray(table.index)

    legacy_meta = pd.read_csv(
        args.source_root / "analysis/analysis-15/04_biomarker_discovery/01_composite_scoring/sample_metadata.csv"
    )
    legacy_values, legacy_labels, _ = mapped_matrix(table, legacy_meta, LEGACY_UNITS)
    legacy_atlas, legacy_summary = build(
        legacy_values, legacy_labels, feature_ids, feature_meta, LEGACY_UNITS, {}
    )
    historical_atlas_path = (
        args.source_root
        / "analysis/analysis-15/04_biomarker_discovery/fingerprint_atlas/"
        / "simper_fingerprint_atlas_enriched.csv"
    )
    frozen_curves_path = args.frozen_curves or (
        args.source_root
        / "manuscript_2_clean/06_figures/source_folders/analysis18_09_figures/"
        / "fig3a_simper_curves/data/simper_curves.csv"
    )
    if historical_atlas_path.exists():
        published_atlas = pd.read_csv(
            historical_atlas_path,
            usecols=["feature_id", "phylum", "contribution_pct"],
            low_memory=False,
        )
        joined = published_atlas.merge(
            legacy_atlas[["feature_id", "phylum", "contribution_pct"]],
            on=["feature_id", "phylum"], suffixes=("_published", "_reproduced")
        )
        reproduction_r = float(np.corrcoef(
            joined["contribution_pct_published"], joined["contribution_pct_reproduced"]
        )[0, 1])
        published_summary = pd.read_csv(
            args.source_root / "analysis/analysis-15/04_biomarker_discovery/fingerprint_atlas/simper_summary.csv"
        )
        validation = published_summary[["phylum", "n_50pct"]].rename(
            columns={"n_50pct": "published_n50"}
        ).merge(
            legacy_summary[["phylum", "n_50pct"]].rename(columns={"n_50pct": "reproduced_n50"}),
            on="phylum",
        )
        validation["delta"] = validation["reproduced_n50"] - validation["published_n50"]
        validation["per_feature_pearson_r"] = reproduction_r
        validation["validation_mode"] = "per_feature_exact"
        if reproduction_r < MIN_REPRODUCTION_R:
            raise AssertionError(
                f"published SIMPER method reproduction r={reproduction_r:.6f} < {MIN_REPRODUCTION_R}"
            )
        validation_status = "pass"
        validation_description = "per-feature Pearson reproduction against recovered historical atlas"
    else:
        if not frozen_curves_path.exists():
            raise FileNotFoundError(
                "Neither the historical SIMPER atlas nor frozen curves are available for validation."
            )
        frozen_curves = pd.read_csv(frozen_curves_path)
        frozen_n50 = (
            frozen_curves.loc[frozen_curves["cumulative_pct"] >= 50]
            .sort_values(["phylum", "simper_rank"])
            .groupby("phylum", as_index=False)
            .first()[["phylum", "simper_rank"]]
            .rename(columns={"simper_rank": "frozen_n50"})
        )
        validation = frozen_n50.merge(
            legacy_summary[["phylum", "n_50pct"]].rename(columns={"n_50pct": "reproduced_n50"}),
            on="phylum",
            how="inner",
        )
        if sorted(validation["phylum"].tolist()) != sorted(frozen_n50["phylum"].tolist()):
            raise AssertionError("Frozen SIMPER curve gate does not cover every historical phylum.")
        validation["delta"] = validation["reproduced_n50"] - validation["frozen_n50"]
        validation["absolute_relative_delta"] = validation["delta"].abs() / validation["frozen_n50"]
        reproduction_r = float(np.corrcoef(validation["frozen_n50"], validation["reproduced_n50"])[0, 1])
        validation["per_feature_pearson_r"] = np.nan
        validation["n50_pearson_r"] = reproduction_r
        validation["validation_mode"] = "bounded_n50_frozen_curve"
        if reproduction_r < MIN_REPRODUCTION_R or validation["absolute_relative_delta"].max() > MAX_FROZEN_N50_RELATIVE_DELTA:
            raise AssertionError(
                "frozen SIMPER curve reproduction failed: "
                f"r={reproduction_r:.6f}, max relative delta="
                f"{validation['absolute_relative_delta'].max():.6f}"
            )
        validation_status = "pass_bounded_not_exact"
        validation_description = "bounded per-phylum n50 reproduction against frozen curves; historical per-feature atlas absent"

    taxonomy = json.loads((args.taxonomy_dir / "taxonomy_summary.json").read_text(encoding="utf-8"))
    policy = json.loads((Path(__file__).resolve().parent / "config/taxonomy_policy.json").read_text(encoding="utf-8"))
    units = taxonomy["analysis_phyla"]
    final_meta = pd.read_csv(args.taxonomy_dir / "sample_metadata_POS_ncbi_phylum.csv")
    final_values, final_labels, mapped = mapped_matrix(table, final_meta, units)
    final_atlas, final_summary = build(
        final_values, final_labels, feature_ids, feature_meta, units, policy["ecological_group"]
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    final_atlas.to_csv(args.output_dir / "simper_atlas.csv", index=False)
    final_summary.to_csv(args.output_dir / "simper_summary.csv", index=False)
    final_atlas[["phylum", "simper_rank", "cumulative_pct"]].to_csv(
        args.output_dir / "simper_curves.csv", index=False
    )
    validation.to_csv(args.output_dir / "reproduction_validation.csv", index=False)
    result = {
        "taxonomy_release": taxonomy["taxonomy_release"],
        "method": "centroid absolute-difference SIMPER on raw intensities; n_batches >= 2; top 3000 per phylum",
        "published_method_reproduction": {
            "validation_mode": validation_status,
            "description": validation_description,
            "reproduction_correlation": reproduction_r,
            "minimum_required_r": MIN_REPRODUCTION_R,
            "status": "pass",
            "shared_rows": int(len(validation)),
            "mean_absolute_n50_delta": float(validation["delta"].abs().mean()),
        },
        "n_features": int(final_values.shape[1]),
        "n_samples": int(len(mapped)),
        "n_phyla": len(units),
        "n50_min": int(final_summary["n_50pct"].min()),
        "n50_max": int(final_summary["n_50pct"].max()),
        "max_single_feature_pct": float(final_summary["max_single_pct"].max()),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    print("REPRODUCTION: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
