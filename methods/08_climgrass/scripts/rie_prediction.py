#!/usr/bin/env python3
"""Predict relative ionisation efficiency (RIE) for archaeal ether lipids from
molecular descriptors, trained on the platform's own Table S10 standards.

Rationale: archaeal ether lipids (archaeol family) have no authentic standards,
so their response factors are unknown and the archaeal fraction of the ClimGrass
decomposition is an upper bound. Ionisation efficiency is predictable from
molecular structure (Kruve and co-workers: descriptor-based logIE regression;
MS2-based prediction in MS2Quant). Here we mirror that logic on-platform:

  1  Build structures (SMILES) for the positive-mode Table S10 standards from
     class + chain templates (deuterium labels ignored; double bonds placed
     methylene-interrupted from delta-9, cis).
  2  Compute RDKit 2D descriptors, append adduct one-hots.
  3  Train logRIE_LPE regressors (random forest + ridge). Severity check:
     LEAVE-ONE-CLASS-OUT cross-validation, because predicting archaeol is by
     construction an out-of-class prediction.
  4  Predict the 546 release-eligible ArchLips structures with an
     applicability-domain (AD) diagnostic: distance to the training cloud in
     standardised descriptor space, expressed as a percentile of the training
     set's own nearest-neighbour distances.

Outputs land in outputs/analysis/<release>/climgrass/rie_prediction_2026-08-08_v1/.
Review-only; nothing here feeds the strict DAG.
"""

from __future__ import annotations

import json
import re
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
from sklearn.linear_model import Ridge
from sklearn.decomposition import PCA

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RELEASE_ROOT = PROJECT_ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1"
OUTPUT = RELEASE_ROOT / "climgrass" / "rie_prediction_2026-08-08_v1"
HANDOFF = PROJECT_ROOT / "cg" / "SHADOW_PC_CLIMGRASS_STRICT_PHYLA_HANDOFF_2026-08-06"
RIE_TABLE = HANDOFF / "data_contract" / "rie_table_s10.csv"
ARCHLIPS = (RELEASE_ROOT / "climgrass" / "strict16_all_simper_sensitivity_2026-08-06"
            / "prepared_inputs" / "archlips_pos_release_eligible_rt_screened.csv")

POS_ADDUCTS = ("M+H", "M+NH4")
RNG_SEED = 20260808


# ---------------------------------------------------------------------------
# Chain and class SMILES templates
# ---------------------------------------------------------------------------
def acyl_tail(n_carbons: int, n_db: int) -> str:
    """Tail after the carbonyl carbon of an N:M acyl (carbonyl = carbon 1).
    Cis double bonds, methylene-interrupted, starting at delta-9 (shifted down
    if the chain is too short)."""
    if n_carbons < 2:
        return ""
    start = 9
    while n_db and start + 3 * (n_db - 1) + 1 > n_carbons:
        start -= 1
    starts = {start + 3 * k for k in range(n_db)} if n_db else set()
    out = []
    i = 2
    while i <= n_carbons:
        if i in starts and i + 1 <= n_carbons:
            # terminal double bond cannot carry a trailing stereo marker
            out.append("C=C" if i + 1 == n_carbons else "/C=C\\")
            i += 2
        else:
            out.append("C")
            i += 1
    return "".join(out)


SPHINGO_TAIL = "CCCCCCCCCCCCC"  # d18:1(4E) carbons 6-18


def gpl(head: str, t1: str, t2: str | None) -> str:
    sn2 = f"(OC(=O){t2})" if t2 is not None else "(O)"
    return f"C(COC(=O){t1}){sn2}COP(=O)(O){head}"


def build_smiles(cls: str, chains: list[tuple[int, int]]) -> str:
    tails = [acyl_tail(n, d) for n, d in chains]
    t1 = tails[0] if tails else ""
    t2 = tails[1] if len(tails) > 1 else None
    t3 = tails[2] if len(tails) > 2 else None
    heads = {
        "PA": "O", "PE": "OCCN", "PG": "OCC(O)CO",
        "PI": "OC1C(O)C(O)C(O)C(O)C1O", "PS": "OCC(N)C(=O)O",
    }
    if cls in ("PA", "PE", "PG", "PI", "PS"):
        return gpl(heads[cls], t1, t2)
    if cls in ("LPA", "LPE", "LPG", "LPI", "LPS"):
        return gpl(heads[cls[1:]], t1, None)
    if cls == "PC":
        return f"C(COC(=O){t1})(OC(=O){t2})COP(=O)([O-])OCC[N+](C)(C)C"
    if cls == "LPC":
        return f"C(COC(=O){t1})(O)COP(=O)([O-])OCC[N+](C)(C)C"
    if cls == "MG":
        return f"C(COC(=O){t1})(O)CO"
    if cls == "DG":
        return f"C(COC(=O){t1})(OC(=O){t2})CO"
    if cls == "TG":
        return f"C(COC(=O){t1})(OC(=O){t2})COC(=O){t3}"
    if cls == "MGDG":
        return f"C(COC(=O){t1})(OC(=O){t2})COC1OC(CO)C(O)C(O)C1O"
    if cls == "DGDG":
        return (f"C(COC(=O){t1})(OC(=O){t2})"
                f"COC1OC(COC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O")
    if cls == "SQDG":
        return f"C(COC(=O){t1})(OC(=O){t2})COC1OC(CS(=O)(=O)O)C(O)C(O)C1O"
    if cls == "DGTS":
        return f"C(COC(=O){t1})(OC(=O){t2})COCCC(C(=O)[O-])[N+](C)(C)C"
    if cls == "CE":
        return f"CC(C)CCCC(C)C1CCC2C3CC=C4CC(OC(=O){t1})CCC4(C)C3CCC12C"
    if cls == "Cer":
        return f"OCC(NC(=O){t1})C(O)/C=C/{SPHINGO_TAIL}"
    if cls == "CerP":
        return f"O=P(O)(O)OCC(NC(=O){t1})C(O)/C=C/{SPHINGO_TAIL}"
    if cls == "HexCer":
        return f"OCC1OC(OCC(NC(=O){t1})C(O)/C=C/{SPHINGO_TAIL})C(O)C(O)C1O"
    if cls == "SM":
        return f"C[N+](C)(C)CCOP(=O)([O-])OCC(NC(=O){t1})C(O)/C=C/{SPHINGO_TAIL}"
    raise ValueError(f"No template for class {cls}")


def parse_species(cls: str, species: str) -> list[tuple[int, int]]:
    """Return acyl chains (carbons, double bonds), sphingoid base excluded
    (the base is part of the class template)."""
    if species.startswith("C15 Ceramide"):
        return [(15, 0)]
    tokens = re.findall(r"(d?)(\d+):(\d+)", species)
    chains = []
    for d_flag, n, db in tokens:
        if d_flag == "d":       # sphingoid base, in template
            continue
        n, db = int(n), int(db)
        if n == 0:              # lyso 0:0 placeholder
            continue
        chains.append((n, db))
    expected = {"TG": 3, "DG": 2, "PC": 2, "PE": 2, "PG": 2, "PI": 2, "PS": 2,
                "PA": 2, "MGDG": 2, "DGDG": 2, "SQDG": 2, "DGTS": 2}
    if cls in expected and len(chains) != expected[cls]:
        raise ValueError(f"{cls} {species}: expected {expected[cls]} chains, got {chains}")
    if cls in ("MG", "CE", "Cer", "CerP", "HexCer", "SM",
               "LPA", "LPC", "LPE", "LPG", "LPI", "LPS") and len(chains) != 1:
        raise ValueError(f"{cls} {species}: expected 1 chain, got {chains}")
    return chains


# ---------------------------------------------------------------------------
# Descriptors
# ---------------------------------------------------------------------------
def descriptor_frame(smiles_list: list[str]) -> pd.DataFrame:
    rows = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            rows.append(None)
            continue
        rows.append(Descriptors.CalcMolDescriptors(mol))
    ok = [r for r in rows if r is not None]
    if not ok:
        raise ValueError("No molecule parsed")
    frame = pd.DataFrame([r if r is not None else {k: np.nan for k in ok[0]} for r in rows])
    return frame


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite: {OUTPUT}")
    staging = OUTPUT.with_name(OUTPUT.name + ".incomplete")
    staging.mkdir(parents=True)

    # --- training table (positive mode only) ---
    rie = pd.read_csv(RIE_TABLE)
    rie = rie[rie["adduct"].isin(POS_ADDUCTS)].reset_index(drop=True)
    smiles, fails = [], []
    for _, row in rie.iterrows():
        try:
            chains = parse_species(row["class"], str(row["species"]))
            smi = build_smiles(row["class"], chains)
            if Chem.MolFromSmiles(smi) is None:
                raise ValueError("RDKit parse failed")
            smiles.append(smi)
        except Exception as exc:
            smiles.append(None)
            fails.append(f"{row['class']} {row['species']}: {exc}")
    rie["smiles"] = smiles
    if fails:
        print("UNPARSED STANDARDS:")
        for f in fails:
            print("  ", f)
    train = rie.dropna(subset=["smiles"]).reset_index(drop=True)
    print(f"Training standards (POS): {len(train)}/{len(rie)} rows, "
          f"{train['class'].nunique()} classes")
    train.to_csv(staging / "standards_smiles.csv", index=False)

    # ArchLips structures are loaded up front so the descriptor screen can
    # reject columns that are numerically unstable on EITHER set (e.g. the
    # information-content descriptors reach ~1e21 on archaeol-sized molecules).
    arch = pd.read_csv(ARCHLIPS, low_memory=False)
    arch = arch.dropna(subset=["archlips_smiles"]).reset_index(drop=True)
    arch_desc_full = descriptor_frame(arch["archlips_smiles"].tolist())

    desc = descriptor_frame(train["smiles"].tolist())
    stable = []
    for c in desc.columns:
        tv, av = desc[c], arch_desc_full[c]
        if (tv.notna().all() and av.notna().all()
                and np.isfinite(tv).all() and np.isfinite(av).all()
                and tv.std() > 0
                and max(tv.abs().max(), av.abs().max()) < 1e8):
            stable.append(c)
    keep = pd.Index(stable)
    desc = desc[keep]
    print(f"Descriptors retained (stable on both sets): {desc.shape[1]}")

    def add_adduct(frame: pd.DataFrame, adducts: pd.Series) -> pd.DataFrame:
        out = frame.copy()
        for a in POS_ADDUCTS:
            out[f"adduct_{a}"] = (adducts.str.replace("[", "", regex=False)
                                  .str.replace("]", "", regex=False)
                                  .str.rstrip("+").eq(a)).astype(float)
        return out

    MECHANISTIC = ["MolLogP", "TPSA", "NumHAcceptors", "NumHDonors",
                   "NumRotatableBonds", "FractionCSP3", "MolWt", "MolMR",
                   "NOCount", "NHOHCount", "LabuteASA", "BalabanJ"]
    mech_cols = [c for c in MECHANISTIC if c in desc.columns]

    def quaternary_n(smiles_series):
        return smiles_series.str.count(r"\[N\+\]").astype(float)

    X_raw = add_adduct(desc, train["adduct"])
    X_raw["quaternary_N"] = quaternary_n(train["smiles"]).to_numpy()
    Xm_raw = add_adduct(desc[mech_cols], train["adduct"])
    Xm_raw["quaternary_N"] = quaternary_n(train["smiles"]).to_numpy()
    y = train["logRIE_LPE"].astype(float).to_numpy()
    mu, sd = X_raw.mean(), X_raw.std().replace(0, 1.0)
    X = ((X_raw - mu) / sd).to_numpy()
    mu_m, sd_m = Xm_raw.mean(), Xm_raw.std().replace(0, 1.0)
    Xm = ((Xm_raw - mu_m) / sd_m).to_numpy()

    # --- leave-one-class-out CV ---
    classes = train["class"].to_numpy()
    models = {
        "random_forest": (lambda: RandomForestRegressor(
            n_estimators=500, min_samples_leaf=2, random_state=RNG_SEED), "full"),
        "ridge": (lambda: Ridge(alpha=10.0), "full"),
        "rf_mechanistic": (lambda: RandomForestRegressor(
            n_estimators=500, min_samples_leaf=2, random_state=RNG_SEED), "mech"),
        "ridge_mechanistic": (lambda: Ridge(alpha=1.0), "mech"),
    }
    feature_sets = {"full": X, "mech": Xm}
    cv_rows = []
    preds_cv = {name: np.full(len(y), np.nan) for name in models}
    for held in sorted(set(classes)):
        mask = classes == held
        for name, (make, fs) in models.items():
            model = make()
            XX = feature_sets[fs]
            model.fit(XX[~mask], y[~mask])
            preds_cv[name][mask] = model.predict(XX[mask])
    for name in models:
        err = preds_cv[name] - y
        cv_rows.append({
            "model": name,
            "logo_rmse_log10": float(np.sqrt(np.mean(err ** 2))),
            "logo_mae_log10": float(np.mean(np.abs(err))),
            "r2_vs_obs": float(np.corrcoef(preds_cv[name], y)[0, 1] ** 2),
        })
    cv_summary = pd.DataFrame(cv_rows)
    print(cv_summary.round(3).to_string(index=False))
    best = cv_summary.sort_values("logo_rmse_log10").iloc[0]["model"]
    logo_rmse = float(cv_summary.set_index("model").loc[best, "logo_rmse_log10"])
    best_fs = models[best][1]
    print(f"Selected model: {best} [{best_fs} features] "
          f"(LOGO RMSE = {logo_rmse:.3f} log10 units, x{10**logo_rmse:.1f} fold)")

    per_row = train[["class", "species", "adduct"]].copy()
    for name in models:
        per_row[f"pred_{name}"] = preds_cv[name]
    per_row["obs_logRIE"] = y
    per_row.to_csv(staging / "cv_predictions.csv", index=False)
    per_class = (per_row.assign(abs_err=lambda f: (f[f"pred_{best}"] - f["obs_logRIE"]).abs())
                 .groupby("class")["abs_err"].agg(["mean", "count"]).reset_index()
                 .rename(columns={"mean": "logo_mae_log10", "count": "n"})
                 .sort_values("logo_mae_log10", ascending=False))
    per_class.to_csv(staging / "cv_per_class.csv", index=False)
    print("\nWorst held-out classes:")
    print(per_class.head(8).round(2).to_string(index=False))

    # --- fit final model on all standards ---
    final_model = models[best][0]()
    X_best = feature_sets[best_fs]
    final_model.fit(X_best, y)

    # --- ArchLips prediction ---
    arch_desc = arch_desc_full.reindex(columns=keep)
    if best_fs == "mech":
        Xa_raw = add_adduct(arch_desc[mech_cols], arch["archlips_adduct"].astype(str))
        Xa_raw["quaternary_N"] = quaternary_n(arch["archlips_smiles"]).to_numpy()
        Xa = ((Xa_raw - mu_m) / sd_m).to_numpy()
        X_ad = Xm
    else:
        Xa_raw = add_adduct(arch_desc, arch["archlips_adduct"].astype(str))
        Xa_raw["quaternary_N"] = quaternary_n(arch["archlips_smiles"]).to_numpy()
        Xa = ((Xa_raw - mu) / sd).to_numpy()
        X_ad = X
    valid = np.isfinite(Xa).all(axis=1)
    print(f"ArchLips structures: {len(arch)} rows, parsed+finite: {int(valid.sum())}")

    pred = np.full(len(arch), np.nan)
    pred[valid] = final_model.predict(Xa[valid])
    if best == "random_forest":
        tree_preds = np.stack([t.predict(Xa[valid]) for t in final_model.estimators_])
        spread = np.full(len(arch), np.nan)
        spread[valid] = tree_preds.std(axis=0)
    else:
        spread = np.full(len(arch), np.nan)

    # --- applicability domain: 5-NN distance in PCA space ---
    pca = PCA(n_components=min(10, X_ad.shape[1], X_ad.shape[0] - 1), random_state=RNG_SEED)
    Z = pca.fit_transform(X_ad)
    Za = pca.transform(np.nan_to_num(Xa, nan=0.0))

    def knn_dist(points, ref, k=5, exclude_self=False):
        d = np.sqrt(((points[:, None, :] - ref[None, :, :]) ** 2).sum(axis=2))
        if exclude_self:
            np.fill_diagonal(d, np.inf)
        return np.sort(d, axis=1)[:, :k].mean(axis=1)

    train_self = knn_dist(Z, Z, exclude_self=True)
    arch_dist = knn_dist(Za, Z)
    ad_percentile = np.array([float((train_self <= dv).mean() * 100) for dv in arch_dist])

    out = arch[["feature_id", "archlips_name", "archlips_adduct", "archlips_mz",
                "archlips_confidence", "archlips_tier", "phylum"]].copy()
    out["pred_logRIE_LPE"] = pred
    out["pred_RIE_LPE"] = 10 ** pred
    out["pred_RIE_lo"] = 10 ** (pred - logo_rmse)
    out["pred_RIE_hi"] = 10 ** (pred + logo_rmse)
    out["rf_tree_sd_log10"] = spread
    out["ad_knn_distance"] = arch_dist
    out["ad_percentile_vs_train"] = ad_percentile
    out["ad_extrapolating"] = ad_percentile >= 100.0
    out.to_csv(staging / "archlips_predicted_rie.csv", index=False)

    q = out["pred_RIE_LPE"].describe(percentiles=[0.05, 0.25, 0.5, 0.75, 0.95])
    print("\nPredicted ArchLips RIE (relative to LPE):")
    print(q.round(3).to_string())
    print(f"AD: median training-percentile {np.nanmedian(ad_percentile):.0f}; "
          f"{(ad_percentile >= 100).mean() * 100:.0f}% beyond every training distance")

    # --- figure ---
    plt.rcParams.update({"font.size": 7, "pdf.fonttype": 42})
    fig, axes = plt.subplots(1, 3, figsize=(9.5, 3.1))
    ax = axes[0]
    ax.scatter(y, preds_cv[best], s=14, c="#0072B2", alpha=0.8)
    lims = [min(y.min(), np.nanmin(preds_cv[best])) - 0.3,
            max(y.max(), np.nanmax(preds_cv[best])) + 0.3]
    ax.plot(lims, lims, color="#999999", linewidth=0.8)
    ax.set_xlabel("observed logRIE (LPE = 0)")
    ax.set_ylabel("predicted logRIE (class held out)")
    ax.set_title(f"a  LOGO CV, {best}\nRMSE = {logo_rmse:.2f} log10 (x{10**logo_rmse:.1f})", loc="left")

    ax = axes[1]
    ax.hist(train_self, bins=20, alpha=0.65, label="training self-distance", color="#0072B2", density=True)
    ax.hist(arch_dist[valid], bins=20, alpha=0.65, label="ArchLips", color="#CC79A7", density=True)
    ax.set_xlabel("5-NN distance in descriptor PCA space")
    ax.set_title("b  Applicability domain", loc="left")
    ax.legend(frameon=False)

    ax = axes[2]
    class_means = train.assign(logRIE=y).groupby("class")["logRIE"].mean().sort_values()
    ax.scatter(class_means.to_numpy(), np.arange(len(class_means)), s=12, c="#0072B2")
    ax.set_yticks(np.arange(len(class_means)))
    ax.set_yticklabels(class_means.index, fontsize=5.5)
    med = float(np.nanmedian(pred))
    ax.axvline(med, color="#CC79A7", linewidth=1.4,
               label=f"ArchLips median ({med:.2f})")
    ax.axvspan(np.nanpercentile(pred, 5), np.nanpercentile(pred, 95),
               color="#CC79A7", alpha=0.18, label="ArchLips 5-95%")
    ax.set_xlabel("logRIE (LPE = 0)")
    ax.set_title("c  Predicted archaeol RIE vs standard classes", loc="left")
    ax.legend(frameon=False, fontsize=5.5)
    for a in axes:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(staging / "rie_prediction_overview.png", dpi=350, facecolor="white")
    fig.savefig(staging / "rie_prediction_overview.pdf", facecolor="white")
    plt.close(fig)

    summary = {
        "training_rows_pos": int(len(train)),
        "training_classes": int(train["class"].nunique()),
        "descriptors": int(desc.shape[1]),
        "unparsed_standards": fails,
        "model_selected": best,
        "logo_rmse_log10": logo_rmse,
        "logo_fold_error": float(10 ** logo_rmse),
        "cv_summary": cv_summary.to_dict("records"),
        "archlips_rows": int(len(arch)),
        "archlips_predicted": int(valid.sum()),
        "pred_RIE_median": float(np.nanmedian(out["pred_RIE_LPE"])),
        "pred_RIE_p5_p95": [float(np.nanpercentile(out["pred_RIE_LPE"], 5)),
                            float(np.nanpercentile(out["pred_RIE_LPE"], 95))],
        "ad_median_percentile": float(np.nanmedian(ad_percentile)),
        "ad_fraction_beyond_training": float((ad_percentile >= 100).mean()),
        "note": ("Deuterium labels ignored in structure building; double bonds "
                 "placed methylene-interrupted from delta-9 (cis). Predictions are "
                 "an out-of-class extrapolation; one authentic archaeol standard "
                 "would anchor or refute the model."),
    }
    (staging / "RUN_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    staging.rename(OUTPUT)
    print(f"\nOutput: {OUTPUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
