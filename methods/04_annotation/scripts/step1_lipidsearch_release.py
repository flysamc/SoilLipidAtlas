#!/usr/bin/env python3
"""Annotation pipeline Step 1 (LipidSearch) for taxonomy release ncbi-phylum-2026-08-04-v1.

Two provenance classes, kept strictly separate and labelled in the manifest:

POSITIVE  exact reuse. The recovered producer step6b_direct_lipidsearch_mapping.py
          already mapped all 273,248 POS consensus features. Nothing is recomputed;
          its output is joined onto the strict POS feature sets. The older
          JSON-routed table (step6_lipidsearch_mapping.py) is carried as a
          cross-check only, never as the primary source.

NEGATIVE  new producer. No NEG equivalent of the POS consensus LipidSearch mapping
          exists anywhere in the package. This script ports the step6b semantics
          (identical tolerances, clustering rule and gold/silver/bronze thresholds)
          to negative mode and records that it is a NEW producer, not a recovery.

Nothing here depends on the GNPS2 fastMASST service.
"""

from __future__ import annotations

import bisect
import csv
import hashlib
import json
import os
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------
# Paths and constants. Tolerances are those documented in Supplementary Method 3
# Step 1 and implemented in the recovered step6b producer.
# ----------------------------------------------------------------------------
RELEASE = "ncbi-phylum-2026-08-04-v1"

ROOT = Path(__file__).resolve().parents[2]
RECOVERY = (ROOT / "external" / "SOILMASS_PRODUCER_RECOVERY_2026-08-04" / "payload"
            / "core" / "workspace")
LS_RAW_DIR = RECOVERY / "soilmass-viewer" / "data" / "LIPIDSEARCH"
POS_LS_DIR = (RECOVERY / "analysis" / "analysis-15" / "04_biomarker_discovery"
              / "02_lipidsearch_annotations")
STEP6B = RECOVERY / "analysis" / "analysis-15" / "scripts" / "step6b_direct_lipidsearch_mapping.py"

REL_DIR = ROOT / "outputs" / "analysis" / RELEASE / "biomarker_discovery"
REL_NEG_DIR = ROOT / "outputs" / "analysis" / RELEASE / "biomarker_discovery_neg"
NEG_CONSENSUS = (ROOT / "analysis" / "analysis-16" / "negative_mode" / "03_alignment"
                 / "consensus_aligned_table.csv")

OUT_DIR = ROOT / "outputs" / "analysis" / RELEASE / "annotation" / "step1_lipidsearch"

MZ_TOL_PPM = 5
RT_TOL_MIN = 0.3
RT_MIN = 1.5
RT_MAX = 25.0

GOLD_MIN_SCORE = 0.2
SILVER_MIN_SCORE = 0.1
BRONZE_MIN_SCORE = 0.01


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def nonempty(series: pd.Series) -> pd.Series:
    text = series.astype(str).str.strip()
    return series.notna() & (text != "") & (
        ~text.str.lower().isin({"nan", "none", "unknown", "unclassified"})
    )


# ----------------------------------------------------------------------------
# NEG producer: ported step6b
# ----------------------------------------------------------------------------
def parse_ls_files(polarity: str) -> tuple[list[dict], list[tuple[str, int]]]:
    """Parse raw LipidSearch exports for one polarity (step6b parsing rules)."""
    suffix = f"{polarity}.raw.txt"
    files = sorted(f for f in os.listdir(LS_RAW_DIR) if f.endswith(suffix))
    print(f"  found {len(files)} {polarity} files ending '{suffix}'")

    annotations: list[dict] = []
    stats: list[tuple[str, int]] = []

    for name in files:
        sample = name.replace(".raw.txt", "")
        with (LS_RAW_DIR / name).open(encoding="utf-8", errors="replace") as handle:
            lines = [ln for ln in handle if not ln.startswith("#") and ln.strip()]
        if len(lines) < 2:
            stats.append((sample, 0))
            continue

        header = lines[0].rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        count = 0
        for line in lines[1:]:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < len(header):
                continue
            try:
                obs_mz = float(cols[idx["ObsMz"]])
                ms2rt = float(cols[idx["MS2RT"]])
                calc_mz = float(cols[idx["CalcMz"]])
                idscore = float(cols[idx["IDScore"]])
            except (ValueError, KeyError):
                continue
            if ms2rt < RT_MIN or ms2rt > RT_MAX:
                continue
            annotations.append({
                "source_file": sample,
                "ObsMz": obs_mz,
                "CalcMz": calc_mz,
                "MS2RT": ms2rt,
                "LipidMolec": cols[idx["LipidMolec"]] if "LipidMolec" in idx else "",
                "LipidID": cols[idx["LipidID"]] if "LipidID" in idx else "",
                "ClassKey": cols[idx["ClassKey"]] if "ClassKey" in idx else "",
                "SubClassKey": cols[idx["SubClassKey"]] if "SubClassKey" in idx else "",
                "Adduct": cols[idx["Adduct"]] if "Adduct" in idx else "",
                "Grade": cols[idx["Grade"]] if "Grade" in idx else "",
                "IDScore": idscore,
                "IonFormula": cols[idx["IonFormula"]] if "IonFormula" in idx else "",
            })
            count += 1
        stats.append((sample, count))

    print(f"  parsed {len(annotations):,} annotation rows")
    return annotations, stats


def cluster_annotations(annotations: list[dict]) -> list[dict]:
    """Greedy cluster on (LipidMolec, Adduct) within 5 ppm / 0.3 min (step6b rule)."""
    by_lipid: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for ann in annotations:
        by_lipid[(ann["LipidMolec"], ann["Adduct"])].append(ann)

    summaries: list[dict] = []
    for (_, _), anns in by_lipid.items():
        anns.sort(key=lambda x: x["ObsMz"])
        used = [False] * len(anns)
        for i in range(len(anns)):
            if used[i]:
                continue
            cluster = [anns[i]]
            used[i] = True
            ref_mz, ref_rt = anns[i]["ObsMz"], anns[i]["MS2RT"]
            for j in range(i + 1, len(anns)):
                if used[j]:
                    continue
                ppm = abs(anns[j]["ObsMz"] - ref_mz) / ref_mz * 1e6
                if ppm > MZ_TOL_PPM:
                    break
                if abs(anns[j]["MS2RT"] - ref_rt) <= RT_TOL_MIN:
                    cluster.append(anns[j])
                    used[j] = True

            grades = [a["Grade"] for a in cluster]
            scores = [a["IDScore"] for a in cluster]
            best = max(cluster, key=lambda x: x["IDScore"])
            n_a = sum(1 for g in grades if g == "A")
            n_b = sum(1 for g in grades if g == "B")
            summaries.append({
                "LipidMolec": best["LipidMolec"],
                "LipidID": best["LipidID"],
                "ClassKey": best["ClassKey"],
                "SubClassKey": best["SubClassKey"],
                "Adduct": best["Adduct"],
                "ObsMz": float(np.median([a["ObsMz"] for a in cluster])),
                "CalcMz": best["CalcMz"],
                "MS2RT": float(np.median([a["MS2RT"] for a in cluster])),
                "best_grade": "A" if n_a else ("B" if n_b else "C"),
                "best_score": max(scores),
                "median_score": float(np.median(scores)),
                "n_samples": len({a["source_file"] for a in cluster}),
                "n_matches": len(cluster),
                "IonFormula": best["IonFormula"],
            })

    print(f"  built {len(summaries):,} annotation clusters")
    return summaries


def apply_tiers(clusters: list[dict]) -> list[dict]:
    """step6b tier thresholds, verbatim."""
    for cl in clusters:
        if (cl["best_grade"] in ("A", "B") and cl["best_score"] >= GOLD_MIN_SCORE
                and cl["n_samples"] >= 2):
            cl["tier"] = "gold"
        elif cl["best_score"] >= SILVER_MIN_SCORE and (
                cl["best_grade"] in ("A", "B") or cl["n_samples"] >= 2):
            cl["tier"] = "silver"
        elif cl["best_score"] >= BRONZE_MIN_SCORE:
            cl["tier"] = "bronze"
        else:
            cl["tier"] = "reject"

    counts = Counter(cl["tier"] for cl in clusters)
    for tier in ("gold", "silver", "bronze", "reject"):
        print(f"    {tier:<7} {counts.get(tier, 0):>8,}")
    kept = [cl for cl in clusters if cl["tier"] != "reject"]
    print(f"  kept after removing rejects: {len(kept):,}")
    return kept


def load_consensus_features(path: Path) -> list[dict]:
    """Read only the four columns needed; the NEG table is 111 MB with sample columns."""
    features: list[dict] = []
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            try:
                features.append({
                    "feature_id": row["feature_id"],
                    "consensus_mz": float(row["consensus_mz"]),
                    "consensus_rt": float(row["consensus_rt"]),
                    "n_batches": int(float(row.get("n_batches", 0) or 0)),
                })
            except (ValueError, KeyError):
                continue
    print(f"  loaded {len(features):,} consensus features")
    return features


def match_to_consensus(clusters: list[dict], features: list[dict]) -> pd.DataFrame:
    """step6b matching and candidate-preference rule, verbatim."""
    clusters.sort(key=lambda x: x["ObsMz"])
    cluster_mzs = [cl["ObsMz"] for cl in clusters]
    tier_rank = {"gold": 0, "silver": 1, "bronze": 2}

    rows: list[dict] = []
    matched = 0
    multi = 0

    for feat in features:
        fmz, frt = feat["consensus_mz"], feat["consensus_rt"]
        lo = bisect.bisect_left(cluster_mzs, fmz * (1 - MZ_TOL_PPM / 1e6))
        hi = bisect.bisect_right(cluster_mzs, fmz * (1 + MZ_TOL_PPM / 1e6))

        candidates = []
        for i in range(lo, hi):
            cl = clusters[i]
            rt_diff = abs(cl["MS2RT"] - frt)
            if rt_diff <= RT_TOL_MIN:
                ppm = abs(cl["ObsMz"] - fmz) / fmz * 1e6
                candidates.append((cl, ppm, rt_diff))

        base = {
            "feature_id": feat["feature_id"],
            "consensus_mz": fmz,
            "consensus_rt": frt,
            "n_batches": feat["n_batches"],
        }
        if not candidates:
            base["matched"] = False
            rows.append(base)
            continue

        if len(candidates) > 1:
            multi += 1
        candidates.sort(key=lambda x: (
            tier_rank.get(x[0]["tier"], 3), -x[0]["n_samples"], x[2], -x[0]["best_score"]))
        best, ppm, rt_diff = candidates[0]
        matched += 1
        base.update({
            "matched": True,
            "LipidMolec": best["LipidMolec"],
            "LipidID": best["LipidID"],
            "ClassKey": best["ClassKey"],
            "SubClassKey": best["SubClassKey"],
            "Adduct": best["Adduct"],
            "Grade": best["best_grade"],
            "IDScore": best["best_score"],
            "MedianScore": best["median_score"],
            "CalcMz": best["CalcMz"],
            "ls_mz": best["ObsMz"],
            "ls_rt": best["MS2RT"],
            "ls_ppm_diff": ppm,
            "ls_rt_diff": rt_diff,
            "n_samples_detected": best["n_samples"],
            "n_matches_total": best["n_matches"],
            "tier": best["tier"],
            "IonFormula": best["IonFormula"],
            "n_candidates": len(candidates),
        })
        rows.append(base)

    print(f"  matched {matched:,}/{len(features):,} "
          f"({100 * matched / max(len(features), 1):.1f}%), multi-candidate {multi:,}")
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------
# Coverage reporting
# ----------------------------------------------------------------------------
def coverage(target: pd.DataFrame, mapping: pd.DataFrame, label: str) -> dict:
    """Coverage of a strict feature set against a LipidSearch mapping table."""
    merged = target.merge(mapping, on="feature_id", how="left", suffixes=("", "_ls"))
    n = len(target)
    has_row = int(merged["matched"].notna().sum())
    annotated = merged[nonempty(merged.get("LipidMolec", pd.Series(dtype=object)))]
    gold_ab = annotated[annotated["Grade"].astype(str).str.upper().isin({"A", "B"})]

    tiers = annotated["tier"].astype(str).value_counts().to_dict() if "tier" in annotated else {}
    result = {
        "set": label,
        "strict_features": n,
        "present_in_mapping": has_row,
        "with_lipid_annotation": len(annotated),
        "annotation_pct": round(100 * len(annotated) / max(n, 1), 2),
        "grade_A_or_B": len(gold_ab),
        "grade_AB_pct": round(100 * len(gold_ab) / max(n, 1), 2),
        "tier_gold": int(tiers.get("gold", 0)),
        "tier_silver": int(tiers.get("silver", 0)),
        "tier_bronze": int(tiers.get("bronze", 0)),
    }
    return result, merged


def by_phylum(merged: pd.DataFrame) -> pd.DataFrame:
    if "phylum" not in merged.columns:
        return pd.DataFrame()
    merged = merged.copy()
    merged["_annot"] = nonempty(merged.get("LipidMolec", pd.Series(dtype=object)))
    merged["_ab"] = merged["_annot"] & merged["Grade"].astype(str).str.upper().isin({"A", "B"})
    grouped = merged.groupby("phylum").agg(
        features=("feature_id", "count"),
        annotated=("_annot", "sum"),
        grade_AB=("_ab", "sum"),
    ).reset_index()
    grouped["annotation_pct"] = (100 * grouped.annotated / grouped.features).round(2)
    grouped["grade_AB_pct"] = (100 * grouped.grade_AB / grouped.features).round(2)
    return grouped.sort_values("features", ascending=False)


# ----------------------------------------------------------------------------
def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest: dict = {
        "taxonomy_release": RELEASE,
        "step": "Supplementary Method 3, Step 1: LipidSearch database matching",
        "generated_utc": utc_now(),
        "producer": str(Path(__file__).resolve()),
        "producer_sha256": sha256(Path(__file__).resolve()),
        "parameters": {
            "mz_tolerance_ppm": MZ_TOL_PPM,
            "rt_tolerance_min": RT_TOL_MIN,
            "rt_window_min": [RT_MIN, RT_MAX],
            "gold_min_idscore": GOLD_MIN_SCORE,
            "silver_min_idscore": SILVER_MIN_SCORE,
            "bronze_min_idscore": BRONZE_MIN_SCORE,
            "note": "Tolerances match Supplementary Method 3 Step 1 and the recovered step6b producer.",
        },
        "positive": {},
        "negative": {},
        "outputs": {},
        "decisions_required": [],
    }
    summaries: list[dict] = []

    # ------------------------------------------------------------------ POS
    print("\n=== POSITIVE: exact reuse of recovered step6b output ===")
    direct = pd.read_csv(POS_LS_DIR / "consensus_lipidsearch_direct.csv", low_memory=False)
    legacy = pd.read_csv(POS_LS_DIR / "consensus_lipidsearch_annotations.csv", low_memory=False)
    for frame in (direct, legacy):
        frame["feature_id"] = frame["feature_id"].astype(str)
        if "matched" not in frame.columns:
            frame["matched"] = nonempty(frame["LipidMolec"])

    pos_atlas = pd.read_csv(REL_DIR / "atlas_pos_strict.csv", low_memory=False)
    pos_atlas["feature_id"] = pos_atlas["feature_id"].astype(str)
    pos_new = pd.read_csv(REL_DIR / "annotation_queue_new_features.csv", low_memory=False)
    pos_new["feature_id"] = pos_new["feature_id"].astype(str)

    keep = ["feature_id", "matched", "LipidMolec", "LipidID", "ClassKey", "SubClassKey",
            "Adduct", "Grade", "IDScore", "tier", "ls_ppm_diff", "ls_rt_diff",
            "n_samples_detected", "n_candidates", "IonFormula"]
    direct_keep = direct[[c for c in keep if c in direct.columns]]

    for frame, label, cols in (
        (pos_atlas, "POS strict atlas", ["feature_id", "phylum"]),
        (pos_new, "POS newly selected", ["feature_id"] + (["phylum"] if "phylum" in pos_new.columns else [])),
    ):
        target = frame[[c for c in cols if c in frame.columns]]
        stats, merged = coverage(target, direct_keep, label)
        summaries.append(stats)
        print(f"  {label}: {stats['with_lipid_annotation']:,}/{stats['strict_features']:,} "
              f"({stats['annotation_pct']}%) annotated, A/B {stats['grade_AB_pct']}%")
        slug = label.lower().replace(" ", "_")
        merged.to_csv(OUT_DIR / f"{slug}_lipidsearch.csv", index=False)
        ph = by_phylum(merged)
        if not ph.empty:
            ph.to_csv(OUT_DIR / f"{slug}_lipidsearch_by_phylum.csv", index=False)

    # cross-check against the superseded JSON-routed mapping
    legacy_keep = legacy[[c for c in keep if c in legacy.columns]]
    xstats, _ = coverage(pos_atlas[["feature_id", "phylum"]], legacy_keep,
                         "POS strict atlas (superseded JSON-routed mapping)")
    summaries.append(xstats)

    manifest["positive"] = {
        "provenance": "exact_reuse_of_recovered_producer_output",
        "primary_source": str(POS_LS_DIR / "consensus_lipidsearch_direct.csv"),
        "primary_source_sha256": sha256(POS_LS_DIR / "consensus_lipidsearch_direct.csv"),
        "primary_source_rows": len(direct),
        "crosscheck_source": str(POS_LS_DIR / "consensus_lipidsearch_annotations.csv"),
        "crosscheck_source_sha256": sha256(POS_LS_DIR / "consensus_lipidsearch_annotations.csv"),
        "crosscheck_source_rows": len(legacy),
        "recovered_producer": str(STEP6B),
        "recovered_producer_sha256": sha256(STEP6B) if STEP6B.exists() else None,
        "note": ("step6b parses the raw LipidSearch exports directly and supersedes "
                 "step6, which routed through a stale GNPS classical-network JSON."),
    }

    # ------------------------------------------------------------------ NEG
    print("\n=== NEGATIVE: NEW producer, step6b semantics ported ===")
    neg_ann, neg_stats = parse_ls_files("NEG")
    neg_clusters = apply_tiers(cluster_annotations(neg_ann))
    neg_features = load_consensus_features(NEG_CONSENSUS)
    neg_map = match_to_consensus(neg_clusters, neg_features)
    neg_map["feature_id"] = neg_map["feature_id"].astype(str)

    neg_out = OUT_DIR / "neg_consensus_lipidsearch_direct.csv"
    neg_map.to_csv(neg_out, index=False)

    neg_atlas = pd.read_csv(REL_NEG_DIR / "strict_atlas_NEG.csv", low_memory=False)
    neg_atlas["feature_id"] = neg_atlas["feature_id"].astype(str)
    neg_target = neg_atlas[[c for c in ("feature_id", "phylum") if c in neg_atlas.columns]]
    neg_cov, neg_merged = coverage(neg_target, neg_map[[c for c in keep if c in neg_map.columns]],
                                   "NEG strict atlas")
    summaries.append(neg_cov)
    print(f"  NEG strict atlas: {neg_cov['with_lipid_annotation']:,}/"
          f"{neg_cov['strict_features']:,} ({neg_cov['annotation_pct']}%) annotated, "
          f"A/B {neg_cov['grade_AB_pct']}%")
    neg_merged.to_csv(OUT_DIR / "neg_strict_atlas_lipidsearch.csv", index=False)
    neg_ph = by_phylum(neg_merged)
    if not neg_ph.empty:
        neg_ph.to_csv(OUT_DIR / "neg_strict_atlas_lipidsearch_by_phylum.csv", index=False)

    manifest["negative"] = {
        "provenance": "NEW_producer_ported_from_recovered_step6b",
        "semantic_source": str(STEP6B),
        "semantic_source_sha256": sha256(STEP6B) if STEP6B.exists() else None,
        "raw_export_dir": str(LS_RAW_DIR),
        "raw_files_parsed": len(neg_stats),
        "raw_annotation_rows": len(neg_ann),
        "clusters_retained": len(neg_clusters),
        "consensus_table": str(NEG_CONSENSUS),
        "consensus_table_sha256": sha256(NEG_CONSENSUS),
        "consensus_features": len(neg_features),
        "matched_features": int(neg_map["matched"].sum()),
        "warning": ("This mapping did not previously exist. It is a new computation, not a "
                    "recovered result, and must be declared as such in any manuscript text."),
    }

    # ------------------------------------------------------------------ outputs
    summary_df = pd.DataFrame(summaries)
    summary_df.to_csv(OUT_DIR / "step1_coverage_summary.csv", index=False)

    for path in sorted(OUT_DIR.glob("*.csv")):
        manifest["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}

    ab_direct = next(s for s in summaries if s["set"] == "POS strict atlas")["grade_A_or_B"]
    ab_legacy = next(s for s in summaries
                     if s["set"].startswith("POS strict atlas (superseded"))["grade_A_or_B"]
    manifest["decisions_required"].append({
        "id": "pos_lipidsearch_mapping_choice",
        "question": ("Which POS LipidSearch mapping defines the Gold tier for the revised "
                     "manuscript: the direct raw-file mapping (step6b) or the superseded "
                     "JSON-routed mapping (step6)?"),
        "grade_AB_direct": ab_direct,
        "grade_AB_superseded": ab_legacy,
        "recommendation": ("Use the direct mapping. step6b was written to replace step6 and "
                           "parses the raw exports instead of a stale network JSON."),
        "status": "OPEN - user decision required before tier rates are published",
    })

    (OUT_DIR / "STEP1_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print("\n=== coverage summary ===")
    print(summary_df.to_string(index=False))
    print(f"\nwrote {OUT_DIR}")


if __name__ == "__main__":
    sys.exit(main())
