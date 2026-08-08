#!/usr/bin/env python3
"""Annotation pipeline Step 5 (RT validation), release ncbi-phylum-2026-08-04-v1.

Supplementary Method 3 Step 5, as submitted:

    "Retention time coherence within molecular families was assessed using
     IQR-based outlier detection. 104 features flagged as RT-uncertain."

PRODUCER STATUS: the exact producer of the reported 104 is NOT recovered. Three
candidates exist in the package and none yields 104:

  step12_rt_validation.py                        IQR within LIPID CLASS  -> 577 anomalies
  platinum_unified_annotations_rt_validated.csv  expected RT range/class -> 209 violations
                                                                            (243 rt_valid=False)
  (the described family-grouped method)          not located

This producer therefore implements the method AS DOCUMENTED (IQR within molecular
family) and additionally computes the class-grouped variant for comparison. Neither
is presented as a reproduction of 104. That number must be treated as unverified
until its producer is recovered.

Outlier rule, standard Tukey fences as used by both candidates:
    outlier if rt < Q1 - 1.5*IQR  or  rt > Q3 + 1.5*IQR
Severity from the group z-score: |z|>3 severe, |z|>2 moderate, else mild.
Groups smaller than MIN_GROUP are skipped: an IQR over 2-3 points is meaningless.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

RELEASE = "ncbi-phylum-2026-08-04-v1"
ROOT = Path(__file__).resolve().parents[2]
WS = (ROOT / "external" / "SOILMASS_PRODUCER_RECOVERY_2026-08-04" / "payload"
      / "core" / "workspace" / "analysis")

REL = ROOT / "outputs" / "analysis" / RELEASE
STEP2 = REL / "annotation" / "step2_harmonization"
OUT_DIR = REL / "annotation" / "step5_rt_validation"

POS_ATLAS = REL / "biomarker_discovery" / "atlas_pos_strict.csv"
NEG_ATLAS = REL / "biomarker_discovery_neg" / "strict_atlas_NEG.csv"
POS_CONSENSUS = WS / "analysis-15" / "03_alignment" / "consensus_aligned_table.csv"
NEG_CONSENSUS = (ROOT / "analysis" / "analysis-16" / "negative_mode" / "03_alignment"
                 / "consensus_aligned_table.csv")
POS_FBMN = WS / "FBMN_all_batches_POS"
NEG_FBMN = WS / "FBMN_all_batches_NEG"

STEP12 = WS / "analysis-15" / "scripts" / "step12_rt_validation.py"

IQR_MULTIPLIER = 1.5
MIN_GROUP = 4


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def discover_cluster_summaries(fbmn_root: Path) -> dict[str, Path]:
    found: dict[str, Path] = {}
    for path in sorted(fbmn_root.glob("batch_*/results/*_cluster_summary.tsv")):
        token = path.name.replace("_cluster_summary.tsv", "")
        if "_" in token and token.split("_", 1)[0].isdigit():
            token = token.split("_", 1)[1]
        found[token] = path
    return found


def build_scan_to_component(summaries: dict[str, Path]) -> dict[str, dict[int, int]]:
    out: dict[str, dict[int, int]] = {}
    for batch, path in summaries.items():
        mapping: dict[int, int] = {}
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                try:
                    idx = int(float(row["cluster index"]))
                    comp = int(float(row["component"]))
                except (KeyError, TypeError, ValueError):
                    continue
                if comp != -1:
                    mapping[idx] = comp
        out[batch] = mapping
    return out


def feature_to_family(consensus: Path, wanted: set[str],
                      scan_to_component: dict[str, dict[int, int]]) -> dict[str, tuple[str, int]]:
    out: dict[str, tuple[str, int]] = {}
    with consensus.open(newline="", encoding="utf-8", errors="replace") as handle:
        for row in csv.DictReader(handle):
            fid = row.get("feature_id")
            if fid not in wanted:
                continue
            batch = str(row.get("batches", "")).split(",")[0].strip()
            try:
                ref_id = int(float(row.get("ref_batch_id", "")))
            except (TypeError, ValueError):
                continue
            table = scan_to_component.get(batch)
            if table is None:
                continue
            comp = table.get(ref_id)
            if comp is not None:
                out[fid] = (batch, comp)
    return out


def detect_outliers(groups: dict, rt: dict[str, float], grouping: str) -> pd.DataFrame:
    """Tukey-fence outlier detection within each group."""
    rows = []
    skipped_small = 0
    for key, fids in groups.items():
        values = [rt[f] for f in fids if f in rt and np.isfinite(rt[f])]
        if len(values) < MIN_GROUP:
            skipped_small += 1
            continue
        arr = np.array(values, dtype=float)
        q1, q3 = np.percentile(arr, [25, 75])
        iqr = q3 - q1
        if iqr <= 0:
            continue
        lower = q1 - IQR_MULTIPLIER * iqr
        upper = q3 + IQR_MULTIPLIER * iqr
        median = float(np.median(arr))
        std = float(np.std(arr, ddof=0))
        for fid in fids:
            value = rt.get(fid)
            if value is None or not np.isfinite(value):
                continue
            if lower <= value <= upper:
                continue
            z = (value - median) / std if std > 0 else 0.0
            severity = "severe" if abs(z) > 3 else "moderate" if abs(z) > 2 else "mild"
            rows.append({
                "feature_id": fid, "grouping": grouping,
                "group_key": f"{key[0]}:{key[1]}" if isinstance(key, tuple) else str(key),
                "group_size": len(values),
                "consensus_rt": round(value, 4),
                "group_median_rt": round(median, 4),
                "q1": round(q1, 4), "q3": round(q3, 4), "iqr": round(iqr, 4),
                "lower_fence": round(lower, 4), "upper_fence": round(upper, 4),
                "rt_deviation_min": round(value - median, 4),
                "z_score": round(z, 3), "severity": severity,
                "rt_status": "RT_uncertain",
            })
    print(f"      groups skipped (< {MIN_GROUP} members): {skipped_small:,}")
    return pd.DataFrame(rows)


def run_mode(mode: str) -> dict:
    print(f"\n=== {mode} ===")
    atlas = pd.read_csv(POS_ATLAS if mode == "POS" else NEG_ATLAS, low_memory=False)
    atlas["feature_id"] = atlas["feature_id"].astype(str)
    harm = pd.read_csv(STEP2 / f"harmonised_annotations_{mode.lower()}.csv", low_memory=False)
    harm["feature_id"] = harm["feature_id"].astype(str)

    rt = {f: float(v) for f, v in zip(atlas.feature_id, atlas.consensus_rt)
          if pd.notna(v)}
    print(f"  features with RT: {len(rt):,}/{len(atlas):,}")

    fbmn = POS_FBMN if mode == "POS" else NEG_FBMN
    consensus = POS_CONSENSUS if mode == "POS" else NEG_CONSENSUS
    families = feature_to_family(consensus, set(atlas.feature_id),
                                 build_scan_to_component(discover_cluster_summaries(fbmn)))
    fam_groups: dict[tuple[str, int], list[str]] = defaultdict(list)
    for fid, key in families.items():
        fam_groups[key].append(fid)
    print(f"  features in families: {len(families):,} across {len(fam_groups):,} families")

    print("    [documented method] IQR within molecular family:")
    fam_df = detect_outliers(fam_groups, rt, "molecular_family")
    print(f"      RT-uncertain: {len(fam_df):,}")

    class_groups: dict[str, list[str]] = defaultdict(list)
    for fid, superclass in zip(harm.feature_id, harm.annotation_superclass):
        text = "" if pd.isna(superclass) else str(superclass).strip()
        if text:
            class_groups[text].append(fid)
    print(f"    [comparison] IQR within lipid superclass ({len(class_groups)} groups):")
    cls_df = detect_outliers(class_groups, rt, "lipid_superclass")
    print(f"      RT-uncertain: {len(cls_df):,}")

    combined = pd.concat([fam_df, cls_df], ignore_index=True)
    if not combined.empty:
        combined = combined.merge(harm[["feature_id", "phylum", "kingdom",
                                        "annotation_superclass", "annotation_tier"]],
                                  on="feature_id", how="left")
    combined.to_csv(OUT_DIR / f"rt_uncertain_{mode.lower()}.csv", index=False)

    def sev(df: pd.DataFrame) -> dict:
        return df.severity.value_counts().to_dict() if not df.empty else {}

    return {
        "mode": mode,
        "features_with_rt": len(rt),
        "features_in_families": len(families),
        "families_assessed": len(fam_groups),
        "rt_uncertain_by_family": int(len(fam_df)),
        "rt_uncertain_by_family_features": int(fam_df.feature_id.nunique()) if not fam_df.empty else 0,
        "family_severity": sev(fam_df),
        "rt_uncertain_by_superclass": int(len(cls_df)),
        "rt_uncertain_by_superclass_features": int(cls_df.feature_id.nunique()) if not cls_df.empty else 0,
        "superclass_severity": sev(cls_df),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries = [run_mode("POS"), run_mode("NEG")]

    flat = pd.DataFrame([{k: (json.dumps(v) if isinstance(v, dict) else v)
                          for k, v in s.items()} for s in summaries])
    flat.to_csv(OUT_DIR / "step5_summary.csv", index=False)

    manifest = {
        "taxonomy_release": RELEASE,
        "step": "Supplementary Method 3, Step 5: RT validation",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "producer": str(Path(__file__).resolve()),
        "producer_sha256": sha256(Path(__file__).resolve()),
        "parameters": {
            "iqr_multiplier": IQR_MULTIPLIER,
            "min_group_size": MIN_GROUP,
            "severity": "|z|>3 severe, |z|>2 moderate, else mild",
        },
        "producer_status": {
            "state": "DOCUMENTED_METHOD_REIMPLEMENTED_ORIGINAL_PRODUCER_UNRECOVERED",
            "submitted_claim": ("Supplementary Method 3 Step 5: 'Retention time coherence "
                                "within molecular families was assessed using IQR-based "
                                "outlier detection. 104 features flagged as RT-uncertain.'"),
            "candidates_checked": {
                "step12_rt_validation.py": {
                    "path": str(STEP12), "sha256": sha256(STEP12) if STEP12.exists() else None,
                    "grouping": "lipid class", "historical_output_rows": 577,
                },
                "platinum_unified_annotations_rt_validated.csv": {
                    "grouping": "expected RT range per class",
                    "historical_violations": 209, "historical_rt_valid_false": 243,
                },
            },
            "warning": ("No located producer yields 104. The submitted figure is therefore "
                        "UNVERIFIED. Nothing here reproduces it and it must not be restated "
                        "in the revision until its producer is recovered."),
        },
        "summary": summaries,
        "outputs": {},
    }
    for path in sorted(OUT_DIR.glob("*.csv")):
        manifest["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    (OUT_DIR / "STEP5_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print("\n=== summary ===")
    for s in summaries:
        print(f"  {s['mode']}: family-grouped {s['rt_uncertain_by_family']:,} "
              f"({s['family_severity']}), superclass-grouped "
              f"{s['rt_uncertain_by_superclass']:,} ({s['superclass_severity']})")
    print(f"\nwrote {OUT_DIR}")


if __name__ == "__main__":
    sys.exit(main())
