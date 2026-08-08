#!/usr/bin/env python3
"""Assemble release-facing annotation summaries without inventing blocked tier actions."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

RELEASE = "ncbi-phylum-2026-08-04-v1"
TIER_ORDER = ["Gold", "Silver", "Bronze", "Unidentified"]
EXPECTED = {"POS": 11371, "NEG": 5697}
EXPECTED_PROPAGATION = {"POS": 155, "NEG": 98}
EXPECTED_ARCHLIPS_RT_EXCLUSIONS = {"POS": 86, "NEG": 0}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_csv_atomic(frame: pd.DataFrame, path: Path) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    os.replace(temporary, path)


def assemble_mode(run_root: Path, mode: str) -> tuple[pd.DataFrame, dict]:
    lower = mode.lower()
    annotation = run_root / "annotation"
    harmonised = pd.read_csv(
        annotation / "step2_harmonization" / f"harmonised_annotations_{lower}.csv",
        low_memory=False,
    )
    propagated = pd.read_csv(
        annotation / "step4_family_propagation" / f"propagated_{lower}_normalised.csv",
        low_memory=False,
    )
    rt = pd.read_csv(
        annotation / "step5_rt_validation" / f"rt_uncertain_{lower}.csv",
        low_memory=False,
    )
    evidence_path = (
        run_root
        / (
            "biomarker_discovery/external_annotation_results/integration_review_only/"
            "strict_pos_annotation_evidence.csv"
            if mode == "POS"
            else "annotation_recovery_neg/integration_review_only/strict_neg_annotation_evidence.csv"
        )
    )
    evidence = pd.read_csv(evidence_path, low_memory=False)

    for label, frame in (("harmonised", harmonised), ("evidence", evidence)):
        if "feature_id" not in frame or frame["feature_id"].duplicated().any():
            raise ValueError(f"{mode} {label} must contain unique feature_id values")
        frame["feature_id"] = frame["feature_id"].astype(str)
    propagated["feature_id"] = propagated["feature_id"].astype(str)
    rt["feature_id"] = rt["feature_id"].astype(str)
    if len(harmonised) != EXPECTED[mode] or set(harmonised.feature_id) != set(evidence.feature_id):
        raise ValueError(f"{mode} strict annotation/evidence denominator mismatch")
    if propagated.feature_id.nunique() != EXPECTED_PROPAGATION[mode]:
        raise ValueError(f"{mode} propagation count differs from the declared release result")

    out = harmonised.copy()
    out = out.rename(
        columns={
            "annotation_tier": "annotation_tier_pre_release",
            "annotation_source": "annotation_source_pre_release",
            "annotation_level": "annotation_level_pre_release",
        }
    )
    out["annotation_tier"] = out["annotation_tier_pre_release"]
    out["annotation_source"] = out["annotation_source_pre_release"]
    out["annotation_level"] = out["annotation_level_pre_release"]
    out["family_propagation_applied"] = False
    out["rt_uncertain_any"] = out.feature_id.isin(set(rt.feature_id))
    out["archlips_rt_excluded"] = False
    out["release_exclusion_reason"] = ""

    prop = propagated.sort_values(["round", "feature_id"]).drop_duplicates("feature_id")
    prop_map = prop.set_index("feature_id")["propagated_label"].astype(str).to_dict()
    prop_mask = out.feature_id.isin(prop_map) & out.annotation_tier.eq("Unidentified")
    if int(prop_mask.sum()) != EXPECTED_PROPAGATION[mode]:
        raise ValueError(f"{mode} propagation is not a pure Unidentified-to-Bronze transition")
    out.loc[prop_mask, "annotation_tier"] = "Bronze"
    out.loc[prop_mask, "annotation_source"] = "Network_propagation"
    out.loc[prop_mask, "annotation_level"] = 1
    out.loc[prop_mask, "annotation_superclass"] = out.loc[prop_mask, "feature_id"].map(prop_map)
    out.loc[prop_mask, "annotation_class_normalised"] = out.loc[prop_mask, "feature_id"].map(prop_map)
    out.loc[prop_mask, "family_propagation_applied"] = True

    arch_excluded = out.annotation_source_pre_release.eq("ArchLips") & out.rt_uncertain_any
    if int(arch_excluded.sum()) != EXPECTED_ARCHLIPS_RT_EXCLUSIONS[mode]:
        raise ValueError(f"{mode} RT-screened ArchLips exclusion count changed")
    out.loc[arch_excluded, "archlips_rt_excluded"] = True
    out.loc[arch_excluded, "release_exclusion_reason"] = "archlips_rt_uncertain"
    out.loc[arch_excluded, "annotation_tier"] = "Unidentified"
    out.loc[arch_excluded, "annotation_source"] = "none_after_rt_screen"
    out.loc[arch_excluded, "annotation_level"] = 0
    for column in (
        "annotation_molecular_species",
        "annotation_class_verbatim",
        "annotation_class_normalised",
        "annotation_superclass",
    ):
        out.loc[arch_excluded, column] = ""

    duplicate_evidence = {
        "phylum",
        "kingdom",
        "annotation_tier",
        "annotation_source",
        "annotation_level",
        "has_family_propagation",
        "has_rt_screened_archlips",
    }
    evidence_columns = [
        column for column in evidence.columns if column != "feature_id" and column not in duplicate_evidence
    ]
    out = out.merge(
        evidence[["feature_id", *evidence_columns]], on="feature_id", how="left", validate="one_to_one"
    )
    out["annotation_release_status"] = "eligible"
    out.loc[out.annotation_tier.eq("Unidentified"), "annotation_release_status"] = "unidentified"
    out.loc[out.archlips_rt_excluded, "annotation_release_status"] = "rt_excluded_archlips"

    summary = {
        "mode": mode,
        "strict_features": len(out),
        "annotated_pre_release": int(out.annotation_tier_pre_release.ne("Unidentified").sum()),
        "family_propagation_upgrades": int(out.family_propagation_applied.sum()),
        "archlips_rt_exclusions": int(out.archlips_rt_excluded.sum()),
        "release_eligible_annotated": int(out.annotation_tier.ne("Unidentified").sum()),
        "release_gold_silver": int(out.annotation_tier.isin(["Gold", "Silver"]).sum()),
        "sirius_formula_evidence": int(out.get("has_sirius_formula", False).fillna(False).sum()),
        "canopus_class_evidence": int(out.get("has_canopus_class", False).fillna(False).sum()),
        "dreams_evidence": int(out.get("has_dreams_result", False).fillna(False).sum()),
        "ms2lda_membership": int(out.get("in_strict_ms2lda_model", False).fillna(False).sum()),
    }
    return out, summary


def build_waterfall(summaries: list[dict]) -> pd.DataFrame:
    rows = []
    for summary in summaries:
        mode = summary["mode"]
        n = summary["strict_features"]
        entries = [
            (1, "Step 1-2", "harmonised primary annotation", "complete", summary["annotated_pre_release"], 0, 0,
             "Direct LipidSearch mapping plus declared harmonisation hierarchy."),
            (2, "Step 3", "diagnostic-ion evidence", "complete_recovered_50_rule_scope", None, 0, 0,
             "Diagnostic calls are already represented in the harmonised primary annotation."),
            (3, "Step 4", "molecular-family propagation", "complete", None,
             summary["family_propagation_upgrades"], 0, "Declared Unidentified-to-Bronze transitions applied."),
            (4, "Step 5", "RT validation", "complete_original_104_unverified", None, 0,
             summary["archlips_rt_exclusions"], "RT-inconsistent ArchLips calls excluded from release-eligible tiers."),
            (5, "Step 6", "RT sum-composition prediction", "review_only_pos_blocked_neg", None, 0, 0,
             "POS replacement candidate retained review-only; no tier mutation. NEG historical producer missing."),
            (6, "Step 7", "custom 124-entry archaeal database", "blocked_asset_missing", None, 0, 0,
             "No substitute database used."),
            (7, "Step 8", "ArchLips", "complete_rt_screened", None, 0, 0,
             "Eligible ArchLips calls are included through the harmonisation source; screened calls remain auditable."),
            (8, "Step 9", "annotation-depth tier actions", "evidence_complete_actions_blocked", None, 0, 0,
             "SIRIUS/DreaMS/MS2LDA evidence attached without recreating missing 564/640/12 tier rules."),
            (9, "Step 10", "SIRIUS/CANOPUS/CSI", "complete_review_only_evidence",
             summary["sirius_formula_evidence"], 0, 0, "Exact-ID evidence attached; no tier mutation."),
            (10, "Step 11", "DreaMS", "complete_review_only_evidence", summary["dreams_evidence"], 0, 0,
             "Exact-ID top-1 evidence attached; no tier mutation."),
            (11, "Method 5", "MS2LDA", "complete_review_only_evidence", summary["ms2lda_membership"], 0, 0,
             "Strict model membership attached; no tier mutation."),
        ]
        cumulative = summary["annotated_pre_release"]
        for order, step, method, status, evidence_features, upgrades, exclusions, note in entries:
            cumulative += upgrades - exclusions
            rows.append(
                {
                    "mode": mode,
                    "stage_order": order,
                    "method_step": step,
                    "method": method,
                    "stage_status": status,
                    "strict_features": n,
                    "evidence_features": evidence_features,
                    "tier_upgrades_applied": upgrades,
                    "release_exclusions": exclusions,
                    "release_annotation_count_after": cumulative,
                    "note": note,
                }
            )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    run_root = args.run_root.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    assembled = []
    summaries = []
    inputs = []
    for mode in ("POS", "NEG"):
        frame, summary = assemble_mode(run_root, mode)
        assembled.append(frame.assign(mode=mode))
        summaries.append(summary)
        path = output_dir / f"annotation_evidence_{mode.lower()}.csv"
        write_csv_atomic(frame.sort_values("feature_id"), path)

    combined = pd.concat(assembled, ignore_index=True)
    tiers = (
        combined.groupby(["mode", "phylum", "kingdom", "annotation_tier"], dropna=False)
        .size()
        .rename("n_features")
        .reset_index()
    )
    tiers["annotation_tier"] = pd.Categorical(tiers.annotation_tier, TIER_ORDER, ordered=True)
    tiers = tiers.sort_values(["mode", "phylum", "annotation_tier"])
    write_csv_atomic(tiers, output_dir / "tier_counts.csv")

    annotated = combined[
        combined.annotation_tier.ne("Unidentified") & combined.annotation_superclass.fillna("").ne("")
    ].copy()
    classes = (
        annotated.groupby(["mode", "phylum", "kingdom", "annotation_superclass"])
        .size()
        .rename("n_features")
        .reset_index()
    )
    classes["phylum_annotated_total"] = classes.groupby(["mode", "phylum"])["n_features"].transform("sum")
    classes["fraction"] = classes.n_features / classes.phylum_annotated_total
    write_csv_atomic(classes, output_dir / "lipid_classes.csv")
    write_csv_atomic(build_waterfall(summaries), output_dir / "waterfall.csv")
    write_csv_atomic(pd.DataFrame(summaries), output_dir / "annotation_release_summary.csv")

    output_paths = sorted(output_dir.glob("*.csv"))
    manifest = {
        "schema_version": 1,
        "taxonomy_release": RELEASE,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "producer": str(Path(__file__).resolve()),
        "producer_sha256": sha256(Path(__file__).resolve()),
        "status": "complete_with_declared_blocked_steps_no_invented_tier_actions",
        "summary": summaries,
        "blocked_steps_quarantined": [
            "Step 6 NEG historical RT producer",
            "Step 7 custom 124-entry archaeal database",
            "Step 9 historical 564/640/12 tier-action producer",
        ],
        "guardrails": [
            "SIRIUS, CANOPUS, CSI:FingerID, DreaMS and MS2LDA are evidence only.",
            "No Step 9 tier transition is inferred from evidence availability.",
            "RT-inconsistent ArchLips evidence is retained in pre-release audit columns but excluded from release tiers.",
            "No manuscript-facing or submission_source file is modified.",
        ],
        "outputs": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in output_paths
        },
    }
    temporary = output_dir / "ANNOTATION_SUMMARIES_MANIFEST.json.tmp"
    temporary.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, output_dir / "ANNOTATION_SUMMARIES_MANIFEST.json")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
