#!/usr/bin/env python3
"""Apply the Step 5 RT screen to the Step 8 ArchLips results and report both views.

Two questions answered:

  1. How many strict features carry an ArchLips annotation, before and after the
     retention-time consistency screen.
  2. How many DISTINCT compounds from the ArchLips spectral library (Zheng et al.,
     2026) those annotations draw on. This counts entries in the external ArchLips
     database only. The project's own 124-entry custom archaeal reference
     (Supplementary Method 3 Step 7) is a different resource and is NOT included:
     its producer has not been located and that step has not been run.

RT-flagged rows are retained and marked, never deleted: the mass and spectral
evidence still exists, it simply conflicts with chromatography.
"""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

RELEASE = "ncbi-phylum-2026-08-04-v1"
ROOT = Path(__file__).resolve().parents[2]
REL = ROOT / "outputs" / "analysis" / RELEASE
STEP5 = REL / "annotation" / "step5_rt_validation"
STEP8 = REL / "annotation" / "step8_archlips"
OUT_DIR = REL / "annotation" / "step8_archlips_rt_filtered"

VALID_TIERS = ["Gold", "Silver", "Bronze"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rt_flagged(mode: str) -> set[str]:
    """Feature ids flagged RT-uncertain by either grouping."""
    path = STEP5 / f"rt_uncertain_{mode.lower()}.csv"
    if not path.exists():
        return set()
    df = pd.read_csv(path)
    return set(df.feature_id.astype(str))


def run(mode: str) -> dict:
    matches = pd.read_csv(STEP8 / f"archlips_{mode.lower()}_strict_matches.csv", low_memory=False)
    if matches.empty:
        return {"mode": mode, "note": "no ArchLips matches"}
    matches["feature_id"] = matches["feature_id"].astype(str)

    flagged = rt_flagged(mode)
    matches["rt_uncertain"] = matches.feature_id.isin(flagged)

    validated = matches[matches.archlips_tier.isin(VALID_TIERS)].copy()
    kept = validated[~validated.rt_uncertain]
    dropped = validated[validated.rt_uncertain]

    def compounds(df: pd.DataFrame) -> int:
        return int(df.archlips_name.astype(str).str.strip().replace("", pd.NA).dropna().nunique())

    def tiers(df: pd.DataFrame) -> dict:
        return (df.drop_duplicates("feature_id").archlips_tier.value_counts().to_dict()
                if not df.empty else {})

    matches.to_csv(OUT_DIR / f"archlips_{mode.lower()}_rt_screened.csv", index=False)

    gs_before = validated[validated.archlips_tier.isin(["Gold", "Silver"])].feature_id.nunique()
    gs_after = kept[kept.archlips_tier.isin(["Gold", "Silver"])].feature_id.nunique()

    return {
        "mode": mode,
        "features_validated_before_rt": int(validated.feature_id.nunique()),
        "features_validated_after_rt": int(kept.feature_id.nunique()),
        "features_rt_excluded": int(dropped.feature_id.nunique()),
        "gold_silver_before_rt": int(gs_before),
        "gold_silver_after_rt": int(gs_after),
        "tiers_before_rt": tiers(validated),
        "tiers_after_rt": tiers(kept),
        "tiers_rt_excluded": tiers(dropped),
        "archlips_compounds_all_matches": compounds(matches),
        "archlips_compounds_validated": compounds(validated),
        "archlips_compounds_after_rt": compounds(kept),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries = [run("POS"), run("NEG")]

    rows = []
    for s in summaries:
        rows.append({k: (json.dumps(v) if isinstance(v, dict) else v) for k, v in s.items()})
    pd.DataFrame(rows).to_csv(OUT_DIR / "archlips_rt_screened_summary.csv", index=False)

    manifest = {
        "taxonomy_release": RELEASE,
        "purpose": ("Step 8 ArchLips results screened against Step 5 retention-time "
                    "consistency; ArchLips reference-compound counts."),
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "producer": str(Path(__file__).resolve()),
        "producer_sha256": sha256(Path(__file__).resolve()),
        "database_scope": ("Counts refer to the external ArchLips spectral library "
                           "(Zheng et al., 2026) ONLY. The project's own 124-entry custom "
                           "archaeal reference database (Supplementary Method 3 Step 7) is a "
                           "separate resource, its producer is unrecovered, and Step 7 has "
                           "not been run for this release."),
        "rt_screen": ("Features flagged RT-uncertain in Step 5 by either grouping are marked "
                      "rt_uncertain=True and excluded from the post-screen counts. Rows are "
                      "retained, never deleted."),
        "summary": summaries,
        "outputs": {},
    }
    for path in sorted(OUT_DIR.glob("*.csv")):
        manifest["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    (OUT_DIR / "ARCHLIPS_RT_SCREENED_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    for s in summaries:
        print(f"\n=== {s['mode']} ===")
        for key in ("features_validated_before_rt", "features_validated_after_rt",
                    "features_rt_excluded", "gold_silver_before_rt", "gold_silver_after_rt",
                    "archlips_compounds_all_matches", "archlips_compounds_validated",
                    "archlips_compounds_after_rt"):
            print(f"  {key:34} {s.get(key)}")
        print(f"  tiers before RT: {s.get('tiers_before_rt')}")
        print(f"  tiers after  RT: {s.get('tiers_after_rt')}")
    print(f"\nwrote {OUT_DIR}")


if __name__ == "__main__":
    sys.exit(main())
