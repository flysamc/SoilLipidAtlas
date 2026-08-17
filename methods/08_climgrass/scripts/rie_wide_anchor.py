#!/usr/bin/env python3
"""Wide-domain logIE model (Kruve-lab unified dataset) anchored to the
platform's Table S10 standards, applied to archaeal ether lipids.

Strategy (the route the platform-only QSPR feasibility test pointed to):
  1  Train a structure->logIE regressor on the Kruve-lab unified ionisation
     efficiency collection (1,421 unique compounds, 13 instruments, unified
     scale; downloaded from github.com/kruvelab/MS2Quant, development/data,
     file all_datasets_unified_IEs_20221121.csv, SHA256 f159d8c2...).
  2  Anchor the wide-scale predictions to the platform's own logRIE_LPE scale
     with a linear map fitted on the 57 positive-mode standards (slope,
     intercept, adduct offset). Anchoring is in-domain interpolation on the
     wide scale, unlike the failed out-of-class training on 57 rows.
  3  Severity test: ANCHORED leave-one-class-out on the standards. For each
     held-out lipid class the anchor is refitted without it and the class is
     predicted through the wide model + anchor. Comparison baseline: the
     platform-only QSPR (LOGO RMSE 1.26 log10).
  4  Predict the 546 release-eligible ArchLips structures, with applicability
     domain measured against BOTH the wide training cloud and the standards.

Review-only; predictions feed at most a clearly-labelled sensitivity arm.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.decomposition import PCA

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RELEASE_ROOT = PROJECT_ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1"
OUTPUT = RELEASE_ROOT / "climgrass" / "rie_prediction_2026-08-08_v2_wide_anchor"
WIDE_CSV = PROJECT_ROOT / "external" / "kruvelab_logIE_2026-08-08" / "all_datasets_unified_IEs_20221121.csv"
STANDARDS = RELEASE_ROOT / "climgrass" / "rie_prediction_2026-08-08_v1" / "standards_smiles.csv"
ARCHLIPS = (RELEASE_ROOT / "climgrass" / "strict16_all_simper_sensitivity_2026-08-06"
            / "prepared_inputs" / "archlips_pos_release_eligible_rt_screened.csv")

RNG_SEED = 20260808
BASELINE_LOGO_RMSE = 1.264  # platform-only QSPR, rie_prediction_2026-08-08_v1


def descriptor_frame(smiles_list):
    rows = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        rows.append(Descriptors.CalcMolDescriptors(mol) if mol is not None else None)
    template = next(r for r in rows if r is not None)
    return pd.DataFrame([r if r is not None else {k: np.nan for k in template} for r in rows])


def main() -> int:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite: {OUTPUT}")
    staging = OUTPUT.with_name(OUTPUT.name + ".incomplete")
    staging.mkdir(parents=True)

    # --- wide dataset: aggregate per canonical structure ---
    wide = pd.read_csv(WIDE_CSV, low_memory=False, encoding="latin-1")
    wide = wide.dropna(subset=["SMILES", "unified_IEs"])
    canon = []
    for smi in wide["SMILES"]:
        mol = Chem.MolFromSmiles(str(smi))
        canon.append(Chem.MolToSmiles(mol) if mol is not None else None)
    wide["canonical"] = canon
    wide = wide.dropna(subset=["canonical"])
    agg = (wide.groupby("canonical")["unified_IEs"]
           .agg(["median", "std", "count"]).reset_index()
           .rename(columns={"median": "logIE_unified"}))
    print(f"Wide training compounds: {len(agg)} "
          f"(replicate sd median {agg['std'].dropna().median():.2f})")

    stds = pd.read_csv(STANDARDS)
    stds = stds.dropna(subset=["smiles"]).reset_index(drop=True)
    arch = pd.read_csv(ARCHLIPS, low_memory=False).dropna(subset=["archlips_smiles"]).reset_index(drop=True)

    # --- descriptors, stability-screened across all three sets ---
    d_wide = descriptor_frame(agg["canonical"].tolist())
    d_std = descriptor_frame(stds["smiles"].tolist())
    d_arch = descriptor_frame(arch["archlips_smiles"].tolist())
    stable = []
    for c in d_wide.columns:
        cols = (d_wide[c], d_std[c], d_arch[c])
        if all(v.notna().all() and np.isfinite(v).all() for v in cols) \
                and d_wide[c].std() > 0 \
                and max(v.abs().max() for v in cols) < 1e8:
            stable.append(c)
    print(f"Stable descriptors: {len(stable)}")
    Xw = d_wide[stable].to_numpy()
    Xs = d_std[stable].to_numpy()
    Xa = d_arch[stable].to_numpy()
    yw = agg["logIE_unified"].to_numpy()

    # --- wide model + internal CV ---
    cv = KFold(n_splits=5, shuffle=True, random_state=RNG_SEED)
    cv_pred = np.full(len(yw), np.nan)
    for tr, te in cv.split(Xw):
        m = RandomForestRegressor(n_estimators=400, min_samples_leaf=2,
                                  random_state=RNG_SEED, n_jobs=-1)
        m.fit(Xw[tr], yw[tr])
        cv_pred[te] = m.predict(Xw[te])
    wide_rmse = float(np.sqrt(np.mean((cv_pred - yw) ** 2)))
    wide_r2 = float(np.corrcoef(cv_pred, yw)[0, 1] ** 2)
    print(f"Wide model 5-fold CV: RMSE {wide_rmse:.2f} log, R2 {wide_r2:.2f}")

    model = RandomForestRegressor(n_estimators=400, min_samples_leaf=2,
                                  random_state=RNG_SEED, n_jobs=-1)
    model.fit(Xw, yw)
    pred_std_wide = model.predict(Xs)
    pred_arch_wide = model.predict(Xa)

    # --- platform anchor: obs logRIE ~ wide prediction + adduct offset ---
    y_std = stds["logRIE_LPE"].astype(float).to_numpy()
    is_nh4 = (stds["adduct"].astype(str).str.contains("NH4")).astype(float).to_numpy()

    def fit_anchor(pred, nh4, obs):
        A = np.column_stack([np.ones(len(pred)), pred, nh4])
        coef, *_ = np.linalg.lstsq(A, obs, rcond=None)
        return coef

    def apply_anchor(coef, pred, nh4):
        return coef[0] + coef[1] * pred + coef[2] * nh4

    coef_all = fit_anchor(pred_std_wide, is_nh4, y_std)
    anchor_fit = apply_anchor(coef_all, pred_std_wide, is_nh4)
    anchor_r2 = float(np.corrcoef(anchor_fit, y_std)[0, 1] ** 2)
    print(f"Anchor (all standards): intercept {coef_all[0]:.2f}, slope {coef_all[1]:.2f}, "
          f"NH4 offset {coef_all[2]:.2f}, in-sample R2 {anchor_r2:.2f}")

    # --- anchored leave-one-class-out ---
    classes = stds["class"].to_numpy()
    logo_pred = np.full(len(y_std), np.nan)
    for held in sorted(set(classes)):
        mask = classes == held
        coef = fit_anchor(pred_std_wide[~mask], is_nh4[~mask], y_std[~mask])
        logo_pred[mask] = apply_anchor(coef, pred_std_wide[mask], is_nh4[mask])
    logo_err = logo_pred - y_std
    logo_rmse = float(np.sqrt(np.mean(logo_err ** 2)))
    logo_mae = float(np.mean(np.abs(logo_err)))
    logo_r2 = float(np.corrcoef(logo_pred, y_std)[0, 1] ** 2)
    print(f"ANCHORED LOGO on standards: RMSE {logo_rmse:.3f} log10 (x{10**logo_rmse:.1f}), "
          f"MAE {logo_mae:.3f}, R2 {logo_r2:.2f}  "
          f"[platform-only baseline: {BASELINE_LOGO_RMSE:.3f} (x{10**BASELINE_LOGO_RMSE:.1f})]")

    cv_rows = stds[["class", "species", "adduct"]].copy()
    cv_rows["obs_logRIE"] = y_std
    cv_rows["pred_wide_scale"] = pred_std_wide
    cv_rows["pred_anchored_logo"] = logo_pred
    cv_rows.to_csv(staging / "standards_anchored_logo.csv", index=False)
    per_class = (cv_rows.assign(abs_err=lambda f: (f["pred_anchored_logo"] - f["obs_logRIE"]).abs())
                 .groupby("class")["abs_err"].agg(["mean", "count"]).reset_index()
                 .rename(columns={"mean": "anchored_logo_mae", "count": "n"})
                 .sort_values("anchored_logo_mae", ascending=False))
    per_class.to_csv(staging / "standards_per_class.csv", index=False)
    print("\nWorst classes (anchored LOGO):")
    print(per_class.head(6).round(2).to_string(index=False))

    # --- archaeol predictions ---
    arch_nh4 = (arch["archlips_adduct"].astype(str).str.contains("NH4")).astype(float).to_numpy()
    pred_arch = apply_anchor(coef_all, pred_arch_wide, arch_nh4)

    # applicability domain vs the wide cloud (PCA + 5-NN percentile)
    mu, sd = Xw.mean(axis=0), Xw.std(axis=0)
    sd[sd == 0] = 1.0
    Zw = (Xw - mu) / sd
    Za = (Xa - mu) / sd
    Zs = (Xs - mu) / sd
    pca = PCA(n_components=10, random_state=RNG_SEED)
    Pw = pca.fit_transform(Zw)
    Pa = pca.transform(Za)
    Ps = pca.transform(Zs)

    def knn(points, ref, k=5, exclude_self=False):
        d = np.sqrt(((points[:, None, :] - ref[None, :, :]) ** 2).sum(axis=2))
        if exclude_self:
            np.fill_diagonal(d, np.inf)
        return np.sort(d, axis=1)[:, :k].mean(axis=1)

    self_d = knn(Pw, Pw, exclude_self=True)
    arch_d = knn(Pa, Pw)
    std_d = knn(Ps, Pw)
    arch_pct = np.array([float((self_d <= v).mean() * 100) for v in arch_d])
    std_pct = np.array([float((self_d <= v).mean() * 100) for v in std_d])
    print(f"\nAD vs wide cloud: standards median percentile {np.median(std_pct):.0f}, "
          f"archaeols median percentile {np.median(arch_pct):.0f}")

    out = arch[["feature_id", "archlips_name", "archlips_adduct", "archlips_mz",
                "archlips_confidence", "archlips_tier", "phylum"]].copy()
    out["pred_logIE_wide_scale"] = pred_arch_wide
    out["pred_logRIE_LPE_anchored"] = pred_arch
    out["pred_RIE_LPE"] = 10 ** pred_arch
    out["pred_RIE_lo"] = 10 ** (pred_arch - logo_rmse)
    out["pred_RIE_hi"] = 10 ** (pred_arch + logo_rmse)
    out["ad_percentile_vs_wide"] = arch_pct
    out.to_csv(staging / "archlips_predicted_rie_anchored.csv", index=False)
    print("\nAnchored ArchLips RIE (vs LPE):")
    print(out["pred_RIE_LPE"].describe(percentiles=[0.05, 0.5, 0.95]).round(3).to_string())

    # --- figure ---
    plt.rcParams.update({"font.size": 7, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(9.8, 3.2))
    ax = axes[0]
    for cls in sorted(set(classes)):
        m = classes == cls
        ax.scatter(y_std[m], logo_pred[m], s=13, alpha=0.85, label=cls)
    lims = [min(y_std.min(), np.nanmin(logo_pred)) - 0.3,
            max(y_std.max(), np.nanmax(logo_pred)) + 0.3]
    ax.plot(lims, lims, color="#999999", linewidth=0.8)
    ax.set_xlabel("observed logRIE (LPE = 0)")
    ax.set_ylabel("anchored prediction (class held out)")
    ax.set_title(f"a  Anchored LOGO: RMSE {logo_rmse:.2f} (x{10**logo_rmse:.1f})\n"
                 f"platform-only baseline {BASELINE_LOGO_RMSE:.2f} (x{10**BASELINE_LOGO_RMSE:.1f})",
                 loc="left")

    ax = axes[1]
    ax.hist(self_d, bins=30, density=True, alpha=0.6, color="#0072B2", label="wide cloud self")
    ax.hist(std_d, bins=15, density=True, alpha=0.6, color="#009E73", label="platform standards")
    ax.hist(arch_d, bins=20, density=True, alpha=0.6, color="#CC79A7", label="ArchLips")
    ax.set_xlabel("5-NN distance in wide descriptor PCA space")
    ax.set_title("b  Applicability domain", loc="left")
    ax.legend(frameon=False, fontsize=6)

    ax = axes[2]
    class_means = (pd.DataFrame({"cls": classes, "y": y_std})
                   .groupby("cls")["y"].mean().sort_values())
    ax.scatter(class_means.to_numpy(), np.arange(len(class_means)), s=12, c="#0072B2")
    ax.set_yticks(np.arange(len(class_means)))
    ax.set_yticklabels(class_means.index, fontsize=5.5)
    med = float(np.nanmedian(pred_arch))
    ax.axvline(med, color="#CC79A7", linewidth=1.4, label=f"ArchLips median ({med:.2f})")
    ax.axvspan(float(np.nanpercentile(pred_arch, 5)), float(np.nanpercentile(pred_arch, 95)),
               color="#CC79A7", alpha=0.18, label="ArchLips 5-95%")
    ax.set_xlabel("logRIE (LPE = 0)")
    ax.set_title("c  Anchored archaeol RIE vs standard classes", loc="left")
    ax.legend(frameon=False, fontsize=5.5)
    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(staging / "rie_wide_anchor_overview.png", dpi=350, facecolor="white")
    fig.savefig(staging / "rie_wide_anchor_overview.pdf", facecolor="white")
    plt.close(fig)

    summary = {
        "wide_dataset": {
            "source": "github.com/kruvelab/MS2Quant development/data/all_datasets_unified_IEs_20221121.csv",
            "sha256": "F159D8C20099D2646C5409C765722CE52A22EB21F1E9C139B0313486D593297F".lower(),
            "unique_compounds": int(len(agg)),
            "replicate_sd_median": float(agg["std"].dropna().median()),
            "cv_rmse_log10": wide_rmse, "cv_r2": wide_r2,
        },
        "anchor": {
            "n_standards": int(len(stds)),
            "intercept": float(coef_all[0]), "slope": float(coef_all[1]),
            "nh4_offset": float(coef_all[2]), "in_sample_r2": anchor_r2,
        },
        "anchored_logo": {
            "rmse_log10": logo_rmse, "fold_error": float(10 ** logo_rmse),
            "mae_log10": logo_mae, "r2": logo_r2,
            "baseline_platform_only_rmse": BASELINE_LOGO_RMSE,
        },
        "archlips": {
            "n": int(len(out)),
            "pred_RIE_median": float(out["pred_RIE_LPE"].median()),
            "pred_RIE_p5_p95": [float(out["pred_RIE_LPE"].quantile(0.05)),
                                float(out["pred_RIE_LPE"].quantile(0.95))],
            "ad_median_percentile_vs_wide": float(np.median(arch_pct)),
            "standards_ad_median_percentile_vs_wide": float(np.median(std_pct)),
        },
    }
    (staging / "RUN_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    staging.rename(OUTPUT)
    print(f"\nOutput: {OUTPUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
