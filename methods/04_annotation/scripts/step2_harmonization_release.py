#!/usr/bin/env python3
"""Annotation pipeline Step 2 (harmonisation) for release ncbi-phylum-2026-08-04-v1.

Supplementary Method 3 Step 2: merge annotation sources under the documented
priority hierarchy:

  (1) LipidSearch molecular species (Grade A/B)
  (2) harmonised molecular species (cross-validated)
  (3) LipidSearch sum composition
  (4) MS2 diagnostic class assignment
  (5) molecular family propagated annotation   <- Step 4, not yet available

ArchLips (Step 8) is an archaeal-specific structural source. It is inserted
immediately after (1), because an ArchLips Gold/Silver call is a spectral match
to an authentic archaeal reference spectrum, which LipidSearch cannot provide at
all (its database holds no archaeal ether lipids).

TWO class fields are emitted:

  annotation_class_verbatim    the winning source's raw string, reproducing the
                               historical behaviour exactly
  annotation_class_normalised  mapped onto one controlled vocabulary

The historical vocabulary was never normalised: it carries both
'Glycerophospholipid' and 'Glycerophospholipids' as distinct values plus a
'Spingolipids' misspelling. Because Step 4 propagates by string agreement at a
50% threshold, those variants counted as disagreement and suppressed real
upgrades. Emitting both fields lets the propagation delta be quantified rather
than asserted.

LS_CLASS_TO_STANDARD, CLASSYFIRE_TO_SUPERCLASS and CLASSYFIRE_SUBCLASS_TO_CLASS
are read straight out of the recovered step11_harmonization.py, not retyped.
The MS2-diagnostic vocabulary has no mapping in that producer, so MS2_TO_SUPERCLASS
below is a NEW declared decision and is recorded as such in the manifest.
"""

from __future__ import annotations

import ast
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

RELEASE = "ncbi-phylum-2026-08-04-v1"
ROOT = Path(__file__).resolve().parents[2]
RECOVERY = (ROOT / "external" / "SOILMASS_PRODUCER_RECOVERY_2026-08-04" / "payload"
            / "core" / "workspace" / "analysis")
STEP11 = RECOVERY / "analysis-15" / "scripts" / "step11_harmonization.py"

REL = ROOT / "outputs" / "analysis" / RELEASE
STEP1 = REL / "annotation" / "step1_lipidsearch"
STEP8 = REL / "annotation" / "step8_archlips"
POS_DIAG = REL / "biomarker_discovery" / "annotation_local_diagnostic" / "diagnostic_ms2_classification.csv"
NEG_DIAG = REL / "annotation_recovery_neg" / "diagnostic_annotation" / "strict_diagnostic_annotation_NEG.csv"
POS_ATLAS = REL / "biomarker_discovery" / "atlas_pos_strict.csv"
NEG_ATLAS = REL / "biomarker_discovery_neg" / "strict_atlas_NEG.csv"
OUT_DIR = REL / "annotation" / "step2_harmonization"

# NEW declared decision: MS2 diagnostic vocabulary -> controlled superclass.
# 'Archaeal lipids' and 'Betaine lipids' are outside the six LIPID MAPS superclasses
# used by LS_CLASS_TO_STANDARD but are required by this dataset.
MS2_TO_SUPERCLASS = {
    "Phospholipid": "Glycerophospholipids",
    "Glycerophospholipid": "Glycerophospholipids",
    "PC": "Glycerophospholipids", "PE": "Glycerophospholipids",
    "PG": "Glycerophospholipids", "PI": "Glycerophospholipids",
    "PS": "Glycerophospholipids", "PA": "Glycerophospholipids",
    "Sphingolipid": "Sphingolipids", "GlcCer": "Sphingolipids",
    "Ceramide": "Sphingolipids", "SM": "Sphingolipids",
    "Glycolipid": "Glycerolipids", "MGDG": "Glycerolipids",
    "DGDG": "Glycerolipids", "SQDG": "Glycerolipids",
    "Sterol": "Sterol lipids",
    "Terpenoid": "Prenol lipids",
    "Pentacyclic_Triterpenoid": "Prenol lipids",
    "Hopanoid": "Prenol lipids",
    "Quinone": "Prenol lipids", "Ubiquinone": "Prenol lipids",
    "Archaeal_Lipid": "Archaeal lipids",
    "Betaine_Lipid": "Betaine lipids",
    "FA": "Fatty acyls", "Fatty_Acid": "Fatty acyls",
    # negative-mode classes absent from the recovered mappings
    "OAHFA": "Fatty acyls",          # (O-acyl)-omega-hydroxy fatty acid
    "MGMG": "Glycerolipids",         # monogalactosyl monoacylglycerol
    "DGMG": "Glycerolipids",         # digalactosyl monoacylglycerol
    "SQMG": "Glycerolipids",         # sulfoquinovosyl monoacylglycerol
}

# Not lipid classes: status flags that must never become an annotation.
MS2_NON_CLASS = {"No_MS2", "No_MS2_spectrum", "Unknown"}

# Combined headgroup calls the diagnostic engine cannot resolve to one superclass.
# m/z 184.07 (phosphocholine) is shared by PC and SM, which sit in different
# superclasses. These are given an explicit joint label rather than being dropped:
# they can still agree with each other during Step 4 propagation, but can never be
# mistaken for a pure PC or pure SM assignment.
MS2_COMBINED = {"PC/SM": "Glycerophospholipids|Sphingolipids"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_recovered_maps() -> dict:
    """Read the mapping dicts out of the recovered producer without executing it."""
    tree = ast.parse(STEP11.read_text(encoding="utf-8", errors="replace"))
    want = {"LS_CLASS_TO_STANDARD", "CLASSYFIRE_TO_SUPERCLASS", "CLASSYFIRE_SUBCLASS_TO_CLASS"}
    out = {}
    for node in tree.body:
        if isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in want:
                    out[target.id] = ast.literal_eval(node.value)
    missing = want - set(out)
    if missing:
        raise SystemExit(f"could not extract {missing} from {STEP11}")
    return out


def clean(value) -> str:
    if value is None:
        return ""
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "unknown", "unclassified"}:
        return ""
    return text


def harmonise(mode: str, maps: dict) -> pd.DataFrame:
    ls_map = maps["LS_CLASS_TO_STANDARD"]

    if mode == "POS":
        atlas = pd.read_csv(POS_ATLAS, low_memory=False)
        step1 = pd.read_csv(STEP1 / "pos_strict_atlas_lipidsearch.csv", low_memory=False)
        diag = pd.read_csv(POS_DIAG, low_memory=False)
        diag_class_col, diag_sub_col = "assigned_class", "assigned_subclass"
        arch = pd.read_csv(STEP8 / "archlips_pos_strict_matches.csv", low_memory=False)
    else:
        atlas = pd.read_csv(NEG_ATLAS, low_memory=False)
        step1 = pd.read_csv(STEP1 / "neg_strict_atlas_lipidsearch.csv", low_memory=False)
        diag = pd.read_csv(NEG_DIAG, low_memory=False)
        diag_class_col, diag_sub_col = "ms2_class", None
        arch = pd.read_csv(STEP8 / "archlips_neg_strict_matches.csv", low_memory=False)

    for frame in (atlas, step1, diag, arch):
        if "feature_id" in frame.columns:
            frame["feature_id"] = frame["feature_id"].astype(str)

    ls = {r.feature_id: r for r in step1.itertuples()}
    dg = {r.feature_id: r for r in diag.itertuples()}

    # ArchLips: keep the best validated call per feature.
    order = {"Gold": 0, "Silver": 1, "Bronze": 2, "Unvalidated": 3}
    arch = arch[arch.get("archlips_tier").notna()].copy()
    arch["_rank"] = arch["archlips_tier"].map(order).fillna(9)
    arch = arch.sort_values("_rank").drop_duplicates("feature_id", keep="first")
    ar = {r.feature_id: r for r in arch.itertuples()}

    rows = []
    for fid, phylum, kingdom in zip(atlas.feature_id.astype(str), atlas.phylum, atlas.kingdom):
        verbatim = normalised = superclass = source = tier = ""
        level = 0
        molec = ""

        l = ls.get(fid)
        a = ar.get(fid)
        d = dg.get(fid)

        ls_class = clean(getattr(l, "ClassKey", "")) if l is not None else ""
        ls_grade = clean(getattr(l, "Grade", "")).upper() if l is not None else ""
        ls_molec = clean(getattr(l, "LipidMolec", "")) if l is not None else ""
        arch_tier = clean(getattr(a, "archlips_tier", "")) if a is not None else ""
        arch_name = clean(getattr(a, "archlips_name", "")) if a is not None else ""

        # (1) LipidSearch molecular species, Grade A/B
        if ls_molec and ls_grade in {"A", "B"}:
            verbatim, molec, source, tier, level = ls_class or ls_molec, ls_molec, "LipidSearch", "Gold", 4
        # (1b) ArchLips Gold/Silver structural match
        elif arch_tier in {"Gold", "Silver"}:
            verbatim, molec, source, tier, level = "Archaeal_Lipid", arch_name, "ArchLips", arch_tier, 4
        # (3) LipidSearch sum composition / lower grade
        elif ls_molec:
            verbatim, molec, source, tier, level = ls_class or ls_molec, ls_molec, "LipidSearch", "Bronze", 3
        # (1c) ArchLips Bronze
        elif arch_tier == "Bronze":
            verbatim, molec, source, tier, level = "Archaeal_Lipid", arch_name, "ArchLips", "Bronze", 2
        # (4) MS2 diagnostic class
        elif d is not None:
            dc = clean(getattr(d, diag_class_col, ""))
            ds = clean(getattr(d, diag_sub_col, "")) if diag_sub_col else ""
            if dc and dc not in MS2_NON_CLASS:
                candidate = ds or dc
                if candidate not in MS2_NON_CLASS:
                    verbatim, source, tier, level = candidate, "MS2_diagnostic", "Bronze", 2

        # normalise
        if verbatim:
            if ls_class and ls_class in ls_map:
                normalised = ls_map[ls_class]["std_class"]
                superclass = ls_map[ls_class]["superclass"]
            elif verbatim in MS2_COMBINED:
                normalised, superclass = verbatim, MS2_COMBINED[verbatim]
            else:
                normalised = verbatim
                superclass = MS2_TO_SUPERCLASS.get(verbatim, "")
                if not superclass:
                    base = clean(getattr(d, diag_class_col, "")) if d is not None else ""
                    superclass = MS2_TO_SUPERCLASS.get(base, "")

        rows.append({
            "feature_id": fid, "phylum": phylum, "kingdom": kingdom,
            "annotation_molecular_species": molec,
            "annotation_class_verbatim": verbatim,
            "annotation_class_normalised": normalised,
            "annotation_superclass": superclass,
            "annotation_source": source or "none",
            "annotation_tier": tier or "Unidentified",
            "annotation_level": level,
        })

    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame, mode: str) -> dict:
    n = len(df)
    annotated = int((df.annotation_tier != "Unidentified").sum())
    gs = int(df.annotation_tier.isin(["Gold", "Silver"]).sum())
    return {
        "mode": mode, "strict_features": n,
        "annotated": annotated, "annotated_pct": round(100 * annotated / max(n, 1), 2),
        "gold_silver": gs, "gold_silver_pct": round(100 * gs / max(n, 1), 2),
        "distinct_class_verbatim": int(df.loc[df.annotation_class_verbatim != "",
                                              "annotation_class_verbatim"].nunique()),
        "distinct_class_normalised": int(df.loc[df.annotation_class_normalised != "",
                                                "annotation_class_normalised"].nunique()),
        "distinct_superclass": int(df.loc[df.annotation_superclass != "",
                                          "annotation_superclass"].nunique()),
        "superclass_unmapped": int(((df.annotation_class_verbatim != "") &
                                    (df.annotation_superclass == "")).sum()),
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    maps = load_recovered_maps()
    print(f"extracted recovered maps: "
          f"{ {k: len(v) for k, v in maps.items()} }")

    summaries = []
    for mode in ("POS", "NEG"):
        df = harmonise(mode, maps)
        path = OUT_DIR / f"harmonised_annotations_{mode.lower()}.csv"
        df.to_csv(path, index=False)
        stats = summarize(df, mode)
        summaries.append(stats)
        print(f"\n=== {mode} ===")
        print(f"  features {stats['strict_features']:,}, annotated {stats['annotated']:,} "
              f"({stats['annotated_pct']}%), Gold+Silver {stats['gold_silver']:,} "
              f"({stats['gold_silver_pct']}%)")
        print(f"  distinct class: verbatim {stats['distinct_class_verbatim']} -> "
              f"normalised {stats['distinct_class_normalised']}, "
              f"superclasses {stats['distinct_superclass']}, "
              f"unmapped {stats['superclass_unmapped']}")
        print(f"  sources: {df.annotation_source.value_counts().to_dict()}")
        sc = df.loc[df.annotation_superclass != "", "annotation_superclass"].value_counts()
        if not sc.empty:
            print(f"  superclass: {sc.to_dict()}")

    pd.DataFrame(summaries).to_csv(OUT_DIR / "step2_summary.csv", index=False)

    manifest = {
        "taxonomy_release": RELEASE,
        "step": "Supplementary Method 3, Step 2: annotation harmonisation",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "producer": str(Path(__file__).resolve()),
        "producer_sha256": sha256(Path(__file__).resolve()),
        "recovered_mapping_source": {"path": str(STEP11), "sha256": sha256(STEP11)},
        "priority_hierarchy": [
            "1. LipidSearch molecular species (Grade A/B)",
            "1b. ArchLips Gold/Silver structural match (archaeal-specific)",
            "3. LipidSearch sum composition / lower grade",
            "1c. ArchLips Bronze",
            "4. MS2 diagnostic class assignment",
            "5. molecular family propagation (Step 4, not yet run)",
        ],
        "decisions_required": [{
            "id": "ms2_vocabulary_normalisation",
            "status": "APPLIED - user approved normalisation",
            "note": ("MS2_TO_SUPERCLASS is a NEW mapping declared by this producer. The "
                     "recovered step11_harmonization.py maps only LipidSearch and ClassyFire "
                     "vocabularies. 'Archaeal lipids' and 'Betaine lipids' extend beyond the "
                     "six LIPID MAPS superclasses because this dataset requires them."),
            "ms2_to_superclass": MS2_TO_SUPERCLASS,
            "combined_headgroup_labels": MS2_COMBINED,
            "non_class_flags_excluded": sorted(MS2_NON_CLASS),
        }],
        "summary": summaries,
        "outputs": {},
    }
    for path in sorted(OUT_DIR.glob("*.csv")):
        manifest["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    (OUT_DIR / "STEP2_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"\nwrote {OUT_DIR}")


if __name__ == "__main__":
    sys.exit(main())
