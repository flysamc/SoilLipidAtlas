#!/usr/bin/env python3
"""Versioned replacement POS RT sum-composition method (review-only v1).

The historical POS caller is unrecovered.  This replacement trains transparent
class-specific RT regressions from strict Grade A/B LipidSearch identities and
learns robust class/adduct mass intercepts.  It never mutates annotation tiers.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RELEASE = "ncbi-phylum-2026-08-04-v1"
METHOD = "pos-rt-sumcomp-v1-2026-08-05"
CH2_MASS = 14.01565
DB_MASS = -2.01565
MIN_CLASS_N = 20
MIN_ADDUCT_N = 8
MIN_CV_R2 = 0.20
MAX_CV_MAE = 2.0
MAX_MASS_ERROR_PPM = 10.0
MAX_RT_ERROR_MIN = 2.0
MIN_RT_MARGIN_MIN = 0.25
CLASS_ALIASES = {"Ceramide": "Cer"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_cdb(value: object) -> tuple[float, float]:
    pairs = re.findall(r"(\d+):(\d+)", str(value))
    if not pairs:
        return np.nan, np.nan
    return float(sum(int(c) for c, _ in pairs)), float(sum(int(db) for _, db in pairs))


def fit_ols(frame: pd.DataFrame) -> np.ndarray:
    x = np.column_stack(
        [np.ones(len(frame)), frame["total_c"].to_numpy(float), frame["total_db"].to_numpy(float)]
    )
    return np.linalg.lstsq(x, frame["consensus_rt"].to_numpy(float), rcond=None)[0]


def predict_ols(beta: np.ndarray, total_c: float, total_db: float) -> float:
    return float(np.array([1.0, total_c, total_db]) @ beta)


def cross_validate(frame: pd.DataFrame) -> dict[str, float]:
    ordered = frame.sort_values("feature_id").reset_index(drop=True)
    observed = ordered["consensus_rt"].to_numpy(float)
    predicted = np.full(len(ordered), np.nan)
    folds = np.arange(len(ordered)) % 5
    for fold in range(5):
        train = ordered.loc[folds != fold]
        test = ordered.loc[folds == fold]
        beta = fit_ols(train)
        predicted[folds == fold] = [
            predict_ols(beta, c, db) for c, db in zip(test["total_c"], test["total_db"])
        ]
    residual = predicted - observed
    denom = float(np.sum((observed - observed.mean()) ** 2))
    return {
        "cv_mae_min": float(np.mean(np.abs(residual))),
        "cv_rmse_min": float(np.sqrt(np.mean(residual**2))),
        "cv_r2": float(1.0 - np.sum(residual**2) / denom) if denom else float("nan"),
    }


def parse_args() -> argparse.Namespace:
    base = ROOT / "outputs/analysis" / RELEASE
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=base / "biomarker_discovery/atlas_pos_strict.csv")
    parser.add_argument(
        "--lipidsearch", type=Path,
        default=base / "annotation/step1_lipidsearch/pos_strict_atlas_lipidsearch.csv",
    )
    parser.add_argument(
        "--harmonised", type=Path,
        default=base / "annotation/step2_harmonization/harmonised_annotations_pos.csv",
    )
    parser.add_argument(
        "--diagnostic", type=Path,
        default=base / "biomarker_discovery/annotation_local_diagnostic/diagnostic_ms2_classification.csv",
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=base / "annotation/step6_rt_prediction_pos_v1_review_only",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    atlas_path = args.atlas.resolve()
    lipid_path = args.lipidsearch.resolve()
    harmonised_path = args.harmonised.resolve()
    diagnostic_path = args.diagnostic.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    atlas = pd.read_csv(atlas_path, usecols=["feature_id", "consensus_mz", "consensus_rt"])
    lipid = pd.read_csv(lipid_path, low_memory=False)
    harmonised = pd.read_csv(harmonised_path, low_memory=False)
    diagnostic = pd.read_csv(
        diagnostic_path,
        usecols=["feature_id", "has_usable_ms2", "assigned_subclass", "classification_confidence"],
        low_memory=False,
    ).rename(
        columns={
            "assigned_subclass": "diagnostic_subclass",
            "classification_confidence": "diagnostic_confidence",
        }
    )

    training = lipid.loc[
        lipid["Grade"].isin(["A", "B"])
        & lipid["LipidMolec"].notna()
        & lipid["ClassKey"].notna()
        & lipid["Adduct"].notna()
    ].merge(atlas, on="feature_id", how="inner", validate="one_to_one")
    parsed = training["LipidMolec"].map(parse_cdb)
    training[["total_c", "total_db"]] = pd.DataFrame(parsed.tolist(), index=training.index)
    training = training.dropna(
        subset=["total_c", "total_db", "consensus_mz", "consensus_rt"]
    ).copy()
    training["total_c"] = training["total_c"].astype(int)
    training["total_db"] = training["total_db"].astype(int)
    training["mass_intercept"] = (
        training["consensus_mz"] - CH2_MASS * training["total_c"] - DB_MASS * training["total_db"]
    )

    mass_models: list[dict[str, object]] = []
    mass_inlier_ids: set[str] = set()
    for (lipid_class, adduct), group in training.groupby(["ClassKey", "Adduct"], sort=True):
        if len(group) < MIN_ADDUCT_N:
            continue
        median_intercept = float(group["mass_intercept"].median())
        residual_ppm = (
            (group["mass_intercept"] - median_intercept).abs()
            / group["consensus_mz"] * 1e6
        )
        inliers = group.loc[residual_ppm <= 20.0].copy()
        if len(inliers) < MIN_ADDUCT_N:
            continue
        median_intercept = float(inliers["mass_intercept"].median())
        residual_ppm = (
            (inliers["mass_intercept"] - median_intercept).abs()
            / inliers["consensus_mz"] * 1e6
        )
        mass_inlier_ids.update(inliers["feature_id"])
        mass_models.append(
            {
                "model_class": str(lipid_class),
                "adduct": str(adduct),
                "n_raw": int(len(group)),
                "n_inlier": int(len(inliers)),
                "mass_intercept": median_intercept,
                "residual_ppm_median": float(residual_ppm.median()),
                "residual_ppm_p95": float(residual_ppm.quantile(0.95)),
                "carbon_min": int(inliers["total_c"].min()),
                "carbon_max": int(inliers["total_c"].max()),
                "db_min": int(inliers["total_db"].min()),
                "db_max": int(inliers["total_db"].max()),
            }
        )
    mass_df = pd.DataFrame(mass_models).sort_values(["model_class", "adduct"])
    mass_df.to_csv(output / "mass_models.csv", index=False, lineterminator="\n")

    rt_models: dict[str, dict[str, object]] = {}
    rt_rows: list[dict[str, object]] = []
    rt_training = training.loc[training["feature_id"].isin(mass_inlier_ids)].copy()
    for lipid_class, group in rt_training.groupby("ClassKey", sort=True):
        if len(group) < MIN_CLASS_N:
            continue
        metrics = cross_validate(group)
        accepted = metrics["cv_r2"] >= MIN_CV_R2 and metrics["cv_mae_min"] <= MAX_CV_MAE
        beta = fit_ols(group)
        row = {
            "model_class": str(lipid_class),
            "n_train": int(len(group)),
            **metrics,
            "intercept": float(beta[0]),
            "coef_total_c": float(beta[1]),
            "coef_total_db": float(beta[2]),
            "carbon_min": int(group["total_c"].min()),
            "carbon_max": int(group["total_c"].max()),
            "db_min": int(group["total_db"].min()),
            "db_max": int(group["total_db"].max()),
            "accepted": bool(accepted),
        }
        rt_rows.append(row)
        if accepted:
            rt_models[str(lipid_class)] = {"beta": beta, **row}
    rt_df = pd.DataFrame(rt_rows).sort_values("model_class")
    rt_df.to_csv(output / "rt_model_validation.csv", index=False, lineterminator="\n")

    targets = harmonised.loc[
        harmonised["annotation_tier"].eq("Bronze")
        & harmonised["annotation_molecular_species"].isna()
        & harmonised["annotation_class_normalised"].notna()
    ].copy()
    targets["model_class"] = targets["annotation_class_normalised"].replace(CLASS_ALIASES)
    targets = targets.merge(atlas, on="feature_id", how="left", validate="one_to_one")
    targets = targets.merge(diagnostic, on="feature_id", how="left", validate="one_to_one")
    targets["diagnostic_model_class"] = targets["diagnostic_subclass"].replace(CLASS_ALIASES)

    candidates: list[dict[str, object]] = []
    target_audit: list[dict[str, object]] = []
    for _, target in targets.sort_values("feature_id").iterrows():
        feature_id = str(target["feature_id"])
        lipid_class = str(target["model_class"])
        rt_model = rt_models.get(lipid_class)
        class_mass = mass_df.loc[mass_df["model_class"].eq(lipid_class)]
        reason = ""
        diagnostic_compatible = (
            bool(target.get("has_usable_ms2", False))
            and str(target.get("diagnostic_model_class", "")) == lipid_class
        )
        if not diagnostic_compatible:
            reason = "current_diagnostic_subclass_not_compatible"
        elif rt_model is None:
            reason = "no_accepted_class_rt_model"
        elif class_mass.empty:
            reason = "no_accepted_class_adduct_mass_model"
        elif pd.isna(target["consensus_mz"]) or pd.isna(target["consensus_rt"]):
            reason = "missing_mz_or_rt"
        else:
            possible: list[dict[str, object]] = []
            for _, mass_model in class_mass.iterrows():
                for total_c in range(int(mass_model["carbon_min"]), int(mass_model["carbon_max"]) + 1):
                    for total_db in range(int(mass_model["db_min"]), int(mass_model["db_max"]) + 1):
                        if total_db > total_c * 0.45:
                            continue
                        predicted_mz = (
                            float(mass_model["mass_intercept"])
                            + CH2_MASS * total_c + DB_MASS * total_db
                        )
                        mass_error_ppm = abs(predicted_mz - float(target["consensus_mz"])) / float(target["consensus_mz"]) * 1e6
                        if mass_error_ppm > MAX_MASS_ERROR_PPM:
                            continue
                        predicted_rt = predict_ols(rt_model["beta"], total_c, total_db)
                        rt_error = abs(predicted_rt - float(target["consensus_rt"]))
                        possible.append(
                            {
                                "feature_id": feature_id,
                                "source_class": target["annotation_class_normalised"],
                                "model_class": lipid_class,
                                "diagnostic_subclass": target["diagnostic_subclass"],
                                "diagnostic_confidence": target["diagnostic_confidence"],
                                "adduct": mass_model["adduct"],
                                "predicted_total_carbon": total_c,
                                "predicted_total_db": total_db,
                                "predicted_sum_species": f"{lipid_class} {total_c}:{total_db}",
                                "observed_mz": float(target["consensus_mz"]),
                                "predicted_mz": predicted_mz,
                                "mass_error_ppm": mass_error_ppm,
                                "observed_rt": float(target["consensus_rt"]),
                                "predicted_rt": predicted_rt,
                                "rt_error_min": rt_error,
                                "rt_model_cv_mae_min": rt_model["cv_mae_min"],
                                "rt_model_cv_r2": rt_model["cv_r2"],
                            }
                        )
            possible.sort(key=lambda row: (row["rt_error_min"], row["mass_error_ppm"], row["adduct"], row["predicted_total_carbon"], row["predicted_total_db"]))
            if not possible:
                reason = "no_candidate_within_mass_tolerance"
            else:
                best = possible[0]
                alternatives = [p for p in possible[1:] if p["predicted_sum_species"] != best["predicted_sum_species"]]
                second_rt = alternatives[0]["rt_error_min"] if alternatives else np.inf
                margin = float(second_rt - best["rt_error_min"]) if np.isfinite(second_rt) else np.inf
                best["alternative_compositions"] = len({p["predicted_sum_species"] for p in possible}) - 1
                best["rt_margin_min"] = margin
                if best["rt_error_min"] > MAX_RT_ERROR_MIN:
                    reason = "best_candidate_exceeds_rt_tolerance"
                elif margin < MIN_RT_MARGIN_MIN:
                    reason = "ambiguous_rt_margin"
                else:
                    if best["mass_error_ppm"] <= 5.0 and best["rt_error_min"] <= 1.0 and margin >= 0.5:
                        best["candidate_confidence"] = "high"
                    else:
                        best["candidate_confidence"] = "medium"
                    best["status"] = "review_only_candidate_no_tier_change"
                    candidates.append(best)
                    reason = "candidate_staged_review_only"
        target_audit.append(
            {
                "feature_id": feature_id,
                "source_class": target["annotation_class_normalised"],
                "model_class": lipid_class,
                "diagnostic_subclass": target.get("diagnostic_subclass", ""),
                "diagnostic_confidence": target.get("diagnostic_confidence", ""),
                "decision": reason,
            }
        )

    candidate_df = pd.DataFrame(candidates)
    candidate_df.to_csv(output / "rt_sum_composition_candidates.csv", index=False, lineterminator="\n")
    audit_df = pd.DataFrame(target_audit)
    audit_df.to_csv(output / "target_audit.csv", index=False, lineterminator="\n")
    decision_counts = audit_df["decision"].value_counts().sort_index().to_dict()
    diagnostic_compatible_targets = int(
        (
            targets["has_usable_ms2"].fillna(False).astype(bool)
            & targets["diagnostic_model_class"].eq(targets["model_class"])
        ).sum()
    )
    confidence_counts = (
        candidate_df["candidate_confidence"].value_counts().sort_index().to_dict()
        if not candidate_df.empty else {}
    )
    checks = {
        "atlas_feature_ids_unique": bool(atlas["feature_id"].is_unique),
        "lipidsearch_feature_ids_unique": bool(lipid["feature_id"].is_unique),
        "harmonised_feature_ids_unique": bool(harmonised["feature_id"].is_unique),
        "diagnostic_feature_ids_unique": bool(diagnostic["feature_id"].is_unique),
        "training_is_grade_a_or_b_only": bool(training["Grade"].isin(["A", "B"]).all()),
        "accepted_rt_models_pass_declared_gates": bool(
            rt_df.loc[rt_df["accepted"], "cv_r2"].ge(MIN_CV_R2).all()
            and rt_df.loc[rt_df["accepted"], "cv_mae_min"].le(MAX_CV_MAE).all()
        ),
        "candidate_feature_ids_unique": bool(candidate_df.empty or candidate_df["feature_id"].is_unique),
        "candidates_have_current_diagnostic_class_agreement": bool(
            candidate_df.empty
            or candidate_df["diagnostic_subclass"].replace(CLASS_ALIASES).eq(candidate_df["model_class"]).all()
        ),
        "candidates_pass_mass_gate": bool(
            candidate_df.empty or candidate_df["mass_error_ppm"].le(MAX_MASS_ERROR_PPM).all()
        ),
        "candidates_pass_rt_gate": bool(
            candidate_df.empty or candidate_df["rt_error_min"].le(MAX_RT_ERROR_MIN).all()
        ),
        "candidates_pass_ambiguity_gate": bool(
            candidate_df.empty or candidate_df["rt_margin_min"].ge(MIN_RT_MARGIN_MIN).all()
        ),
        "no_annotation_atlas_output_written": not (output / "atlas_unified_annotations.csv").exists(),
    }

    script_path = Path(__file__).resolve()
    outputs = {}
    for path in sorted(output.glob("*.csv")):
        outputs[path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    manifest = {
        "schema_version": 1,
        "method_id": METHOD,
        "taxonomy_release": RELEASE,
        "polarity": "POS",
        "status": "pass_review_only_no_tier_mutation" if all(checks.values()) else "fail",
        "method_date": "2026-08-05",
        "producer": str(script_path.relative_to(ROOT)).replace("\\", "/"),
        "producer_sha256": sha256(script_path),
        "historical_reproduction": False,
        "method_change_authority": "User explicitly approved a new versioned POS RT method on 2026-08-05.",
        "inputs": {
            "strict_atlas": {"path": str(atlas_path.relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(atlas_path)},
            "strict_lipidsearch": {"path": str(lipid_path.relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(lipid_path)},
            "harmonised_annotations": {"path": str(harmonised_path.relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(harmonised_path)},
            "current_diagnostic_classification": {"path": str(diagnostic_path.relative_to(ROOT)).replace("\\", "/"), "sha256": sha256(diagnostic_path)},
        },
        "contract": {
            "training": "strict POS Grade A/B LipidSearch identities with parseable C:DB, exact class and adduct",
            "rt_model": "per-class OLS: RT ~ total_carbon + total_double_bonds; deterministic five-fold CV",
            "mass_model": "robust median class/adduct intercept for mz = intercept + 14.01565*C - 2.01565*DB; 20 ppm training-inlier filter",
            "model_acceptance": {"min_n": MIN_CLASS_N, "min_cv_r2": MIN_CV_R2, "max_cv_mae_min": MAX_CV_MAE},
            "candidate_acceptance": {"max_mass_error_ppm": MAX_MASS_ERROR_PPM, "max_rt_error_min": MAX_RT_ERROR_MIN, "min_rt_margin_min": MIN_RT_MARGIN_MIN},
            "target": "strict Bronze class-level feature without an existing molecular/sum species; exact current usable-MS2 diagnostic subclass agreement required",
            "tier_action": "none; candidates are review-only",
        },
        "counts": {
            "grade_ab_parseable_training_rows": int(len(training)),
            "mass_inlier_training_rows": int(len(rt_training)),
            "mass_models": int(len(mass_df)),
            "rt_models_evaluated": int(len(rt_df)),
            "rt_models_accepted": int(rt_df["accepted"].sum()),
            "audited_bronze_class_level_targets": int(len(targets)),
            "current_diagnostic_compatible_targets": diagnostic_compatible_targets,
            "review_only_candidates": int(len(candidate_df)),
            "candidate_confidence": confidence_counts,
            "target_decisions": {str(k): int(v) for k, v in decision_counts.items()},
        },
        "guardrails": [
            "No historical pickle was loaded or reverse-engineered into invented feature names.",
            "No mass-only reuse or prediction is allowed.",
            "Ambiguous PC/SM and generic classes are excluded because no exact class model exists.",
            "No annotation tier or manuscript-facing file is modified.",
            "This is a versioned method replacement, not reproduction of the submitted Step 6.",
        ],
        "checks": checks,
        "outputs": outputs,
    }
    (output / "stage_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if manifest["status"] == "fail":
        raise SystemExit("POS RT replacement v1 validation failed")
    print(json.dumps({"status": manifest["status"], "method_id": METHOD, "counts": manifest["counts"]}, indent=2))


if __name__ == "__main__":
    main()
