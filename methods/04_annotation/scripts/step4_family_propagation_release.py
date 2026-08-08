#!/usr/bin/env python3
"""Annotation pipeline Step 4 (molecular family propagation), release ncbi-phylum-2026-08-04-v1.

Supplementary Method 3 Step 4: propagate class-level annotations from annotated to
unannotated neighbours within GNPS molecular families.

Rules, from the recovered network_propagation_neg.py:
  - map feature -> (reference batch, batch cluster index) -> molecular family component
  - singleton families (component -1) are excluded
  - for each Unidentified feature in a non-singleton family, look at annotated members
  - require >= 1 annotated member and >= 50% agreement on the class label
  - upgrade Unidentified -> Bronze, source Network_propagation
  - two rounds; round 2 may use round 1's upgrades as sources

Propagation is run twice over, on two different label columns, so the effect of
vocabulary normalisation is measured rather than assumed:

  verbatim    the winning source's raw string, as the historical pipeline used
  normalised  the controlled superclass from Step 2

The historical vocabulary held 'Glycerophospholipid' and 'Glycerophospholipids'
as separate strings plus a 'Spingolipids' misspelling, so string agreement at the
50% threshold was suppressed. The normalised run is the release default; the
verbatim run is retained as the comparison.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

RELEASE = "ncbi-phylum-2026-08-04-v1"
ROOT = Path(__file__).resolve().parents[2]
WS = (ROOT / "external" / "SOILMASS_PRODUCER_RECOVERY_2026-08-04" / "payload"
      / "core" / "workspace" / "analysis")

REL = ROOT / "outputs" / "analysis" / RELEASE
STEP2 = REL / "annotation" / "step2_harmonization"
OUT_DIR = REL / "annotation" / "step4_family_propagation"

POS_CONSENSUS = WS / "analysis-15" / "03_alignment" / "consensus_aligned_table.csv"
NEG_CONSENSUS = (ROOT / "analysis" / "analysis-16" / "negative_mode" / "03_alignment"
                 / "consensus_aligned_table.csv")
POS_FBMN = WS / "FBMN_all_batches_POS"
NEG_FBMN = WS / "FBMN_all_batches_NEG"

MIN_AGREEMENT = 0.50
MIN_ANNOTATED = 1
N_ROUNDS = 2


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def discover_cluster_summaries(fbmn_root: Path) -> dict[str, Path]:
    """Map the batch token embedded in each cluster-summary filename to its path."""
    found: dict[str, Path] = {}
    for path in sorted(fbmn_root.glob("batch_*/results/*_cluster_summary.tsv")):
        token = path.name.replace("_cluster_summary.tsv", "")
        # filenames may carry a numeric prefix, e.g. 01_OE11-3-NEG
        if "_" in token and token.split("_", 1)[0].isdigit():
            token = token.split("_", 1)[1]
        found[token] = path
    return found


def build_scan_to_component(summaries: dict[str, Path]) -> dict[str, dict[int, int]]:
    """batch token -> {cluster index -> component}, excluding singletons."""
    out: dict[str, dict[int, int]] = {}
    for batch, path in summaries.items():
        mapping: dict[int, int] = {}
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                try:
                    idx = int(float(row["cluster index"]))
                    comp = int(float(row["component"]))
                except (KeyError, TypeError, ValueError):
                    continue
                if comp == -1:
                    continue
                mapping[idx] = comp
        out[batch] = mapping
        print(f"    {batch}: {len(mapping):,} clustered scans")
    return out


def feature_to_family(consensus: Path, wanted: set[str],
                      scan_to_component: dict[str, dict[int, int]]) -> dict[str, tuple[str, int]]:
    """feature_id -> (batch, component) using the reference batch and its cluster index."""
    out: dict[str, tuple[str, int]] = {}
    unmatched_batches: Counter = Counter()
    with consensus.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
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
                unmatched_batches[batch] += 1
                continue
            comp = table.get(ref_id)
            if comp is not None:
                out[fid] = (batch, comp)
    if unmatched_batches:
        print(f"    batches without a cluster summary: {dict(unmatched_batches)}")
    return out


def propagate(labels: dict[str, str], families: dict[str, tuple[str, int]]) -> list[dict]:
    """Two-round majority propagation within molecular families."""
    members: dict[tuple[str, int], list[str]] = defaultdict(list)
    for fid, key in families.items():
        members[key].append(fid)

    current = dict(labels)
    upgrades: list[dict] = []

    for rnd in range(1, N_ROUNDS + 1):
        round_upgrades = []
        for key, fids in members.items():
            if len(fids) < 2:
                continue
            annotated = [current[f] for f in fids if current.get(f)]
            if len(annotated) < MIN_ANNOTATED:
                continue
            counts = Counter(annotated)
            top_label, top_n = counts.most_common(1)[0]
            if top_n / len(annotated) < MIN_AGREEMENT:
                continue
            for fid in fids:
                if current.get(fid):
                    continue
                round_upgrades.append({
                    "feature_id": fid, "batch": key[0], "component": key[1],
                    "propagated_label": top_label,
                    "family_size": len(fids), "annotated_members": len(annotated),
                    "agreement": round(top_n / len(annotated), 3),
                    "round": rnd,
                })
        for up in round_upgrades:
            current[up["feature_id"]] = up["propagated_label"]
        upgrades.extend(round_upgrades)
        print(f"      round {rnd}: {len(round_upgrades):,} upgrades")
        if not round_upgrades:
            break
    return upgrades


def run_mode(mode: str) -> dict:
    print(f"\n=== {mode} ===")
    harm = pd.read_csv(STEP2 / f"harmonised_annotations_{mode.lower()}.csv", low_memory=False)
    harm["feature_id"] = harm["feature_id"].astype(str)
    wanted = set(harm.feature_id)

    fbmn = POS_FBMN if mode == "POS" else NEG_FBMN
    consensus = POS_CONSENSUS if mode == "POS" else NEG_CONSENSUS

    print("  cluster summaries:")
    summaries = discover_cluster_summaries(fbmn)
    scan_to_component = build_scan_to_component(summaries)

    print("  mapping features to families...")
    families = feature_to_family(consensus, wanted, scan_to_component)
    fam_keys = {k for k in families.values()}
    print(f"    {len(families):,}/{len(wanted):,} features in {len(fam_keys):,} families")

    results = {}
    for variant, column in (("normalised", "annotation_superclass"),
                            ("verbatim", "annotation_class_verbatim")):
        labels = {}
        for fid, val in zip(harm.feature_id, harm[column]):
            text = "" if pd.isna(val) else str(val).strip()
            labels[fid] = text if text and text.lower() not in {"nan", "unclassified"} else ""
        print(f"    propagating on {variant} ({column}):")
        ups = propagate(labels, families)
        df = pd.DataFrame(ups)
        if not df.empty:
            df = df.merge(harm[["feature_id", "phylum", "kingdom"]], on="feature_id", how="left")
        df.to_csv(OUT_DIR / f"propagated_{mode.lower()}_{variant}.csv", index=False)
        results[variant] = {
            "upgrades": int(len(df)),
            "features_upgraded": int(df.feature_id.nunique()) if not df.empty else 0,
            "round1": int((df["round"] == 1).sum()) if not df.empty else 0,
            "round2": int((df["round"] == 2).sum()) if not df.empty else 0,
        }

    unannotated = int((harm.annotation_tier == "Unidentified").sum())
    return {
        "mode": mode,
        "strict_features": int(len(harm)),
        "unidentified_before": unannotated,
        "features_in_families": len(families),
        "families": len(fam_keys),
        "normalised_upgrades": results["normalised"]["features_upgraded"],
        "verbatim_upgrades": results["verbatim"]["features_upgraded"],
        "delta_from_normalisation": (results["normalised"]["features_upgraded"]
                                     - results["verbatim"]["features_upgraded"]),
        "detail": results,
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    summaries = [run_mode("POS"), run_mode("NEG")]

    flat = pd.DataFrame([{k: v for k, v in s.items() if k != "detail"} for s in summaries])
    flat.to_csv(OUT_DIR / "step4_summary.csv", index=False)

    manifest = {
        "taxonomy_release": RELEASE,
        "step": "Supplementary Method 3, Step 4: molecular family propagation",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "producer": str(Path(__file__).resolve()),
        "producer_sha256": sha256(Path(__file__).resolve()),
        "semantic_source": {
            "path": str(WS / "analysis-16" / "negative_mode" / "04_biomarker_discovery"
                        / "08_network_propagation" / "network_propagation_neg.py"),
        },
        "parameters": {
            "min_agreement": MIN_AGREEMENT, "min_annotated": MIN_ANNOTATED,
            "rounds": N_ROUNDS, "singletons_excluded": True,
            "upgrade_tier": "Bronze", "upgrade_source": "Network_propagation",
        },
        "vocabulary_decision": {
            "release_default": "normalised (annotation_superclass)",
            "comparison": "verbatim (annotation_class_verbatim)",
            "why": ("The historical vocabulary was unnormalised, so string agreement at the "
                    "50% threshold suppressed genuine upgrades. Both runs are emitted so the "
                    "difference is measured, not asserted."),
        },
        "summary": summaries,
        "outputs": {},
    }
    for path in sorted(OUT_DIR.glob("*.csv")):
        manifest["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    (OUT_DIR / "STEP4_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print("\n=== summary ===")
    print(flat.to_string(index=False))
    print(f"\nwrote {OUT_DIR}")


if __name__ == "__main__":
    sys.exit(main())
