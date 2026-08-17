#!/usr/bin/env python3
"""Build corrected Figure 2b/2c data tables under the strict 16-phylum release.

Panel b: kingdom x ReDU-sample-type detection of searched atlas biomarkers
         (IndVal + composite fastMASST runs, one current index), using the
         declared sampletype_category_map_v1.csv (the historical categoriser
         is unrecovered; this mapping is a new versioned artifact).
Panel c: soil detection rate, cross-batch consensus (IndVal) selection vs
         composite-scoring selection.

Match->sample-type resolution is FILE-level: each match USI's (dataset,
filename) is joined to the full ReDU metadata; unmatched filenames fall back
to dataset-level majority sample type, and unresolvable matches are counted
and reported, never silently dropped.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def norm_file(name: str) -> str:
    base = str(name).replace("\\", "/").split("/")[-1].lower()
    return re.sub(r"\.(mzml|mzxml|raw|cdf|d|wiff|mgf)$", "", base)


def iter_feature_pairs(per_feature_dir: Path):
    """Yield (feature_id, set of unique (dataset, filename)) per feature.

    Streams one file at a time so the ~100M total match rows never sit in
    memory together; per-feature unique pairs are small.
    """
    for path in sorted(per_feature_dir.glob("*_matches.tsv.gz")):
        feature_id = path.name[: -len("_matches.tsv.gz")]
        pairs: set[tuple[str, str]] = set()
        with gzip.open(path, "rt", encoding="utf-8", errors="ignore") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            try:
                usi_i = header.index("USI")
            except ValueError:
                usi_i = 0
            try:
                ds_i = header.index("Dataset")
            except ValueError:
                ds_i = 1
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) <= max(usi_i, ds_i):
                    continue
                usi = parts[usi_i]
                dataset = parts[ds_i]
                fname = ""
                bits = usi.split(":")
                if len(bits) >= 3:
                    if not dataset:
                        dataset = bits[1]
                    fname = bits[2]
                pairs.add((dataset, norm_file(fname)))
        yield feature_id, pairs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", required=True, type=Path)
    parser.add_argument("--indval-run", required=True, type=Path)
    parser.add_argument("--composite-run", required=True, type=Path)
    parser.add_argument("--redu-all", required=True, type=Path)
    parser.add_argument("--soil-metadata", required=True, type=Path)
    parser.add_argument("--category-map", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    atlas = pd.read_csv(args.atlas, low_memory=False)
    atlas["feature_id"] = atlas["feature_id"].astype(str)
    kingdom_of = dict(zip(atlas.feature_id, atlas.kingdom))
    method_of = dict(zip(atlas.feature_id, atlas.discovery_method))

    cat_map = pd.read_csv(args.category_map)
    type_to_cat = dict(zip(cat_map.redu_sample_type, cat_map.category))

    redu = pd.read_csv(args.redu_all, sep="\t", low_memory=False)
    redu["ds"] = redu["ATTRIBUTE_DatasetAccession"].astype(str)
    redu["nf"] = redu["filename"].map(norm_file)
    file_type = dict(zip(zip(redu.ds, redu.nf), redu.SampleType))
    ds_majority = (
        redu.groupby("ds")["SampleType"]
        .agg(lambda s: s.value_counts().index[0])
        .to_dict()
    )

    soil = pd.read_csv(args.soil_metadata, sep="\t", low_memory=False)
    soil_datasets = set(soil["dataset_id"].astype(str)) if "dataset_id" in soil else set(soil["ATTRIBUTE_DatasetAccession"].astype(str))

    counters = {"file_level": 0, "dataset_fallback": 0, "unresolved": 0}
    feature_cats: dict[str, set[str]] = defaultdict(set)
    feature_soil: dict[str, bool] = {}

    for run_dir in (args.indval_run, args.composite_run):
        for fid, pairs in iter_feature_pairs(run_dir / "per_feature"):
            in_soil = False
            for ds, nf in pairs:
                if ds in soil_datasets:
                    in_soil = True
                stype = file_type.get((ds, nf))
                if stype is not None:
                    counters["file_level"] += 1
                elif ds in ds_majority:
                    stype = ds_majority[ds]
                    counters["dataset_fallback"] += 1
                else:
                    counters["unresolved"] += 1
                    continue
                cat = type_to_cat.get(stype)
                if cat and cat != "EXCLUDED":
                    feature_cats[fid].add(cat)
            feature_soil[fid] = feature_soil.get(fid, False) or in_soil

    searched = sorted(feature_soil)
    categories = [c for c in cat_map.category.unique() if c != "EXCLUDED"]

    # Panel b table (kingdom x category), frozen schema
    rows_b = []
    frame = pd.DataFrame({"feature_id": searched})
    frame["kingdom"] = frame.feature_id.map(kingdom_of)
    for kingdom, g in frame.groupby("kingdom"):
        row = {"kingdom": kingdom, "n_biomarkers": len(g)}
        any_ids = 0
        for cat in categories:
            n = sum(1 for f in g.feature_id if cat in feature_cats.get(f, ()))
            row[f"n_{cat}"] = n
            row[f"pct_{cat}"] = round(100 * n / len(g), 1)
        n_any = sum(1 for f in g.feature_id if feature_cats.get(f))
        row["n_any"] = n_any
        row["pct_any"] = round(100 * n_any / len(g), 1)
        rows_b.append(row)
    panel_b = pd.DataFrame(rows_b)
    panel_b.to_csv(args.output_dir / "kingdom_sampletype_summary_strict.csv", index=False)

    # Panel c table, frozen schema
    rows_c = []
    # Labels renamed 2026-08-11: a coauthor found "cross-batch consensus" vs
    # "composite atlas" confusing; the selection-method names are used instead.
    for label, method in (
        ("Indicator Value (IndVal)", "indval_consensus"),
        ("Composite scoring", "composite"),
    ):
        ids = [f for f in searched if method_of.get(f) == method]
        n_soil = sum(1 for f in ids if feature_soil.get(f))
        rows_c.append(
            {
                "category": label,
                "n_total": len(ids),
                "n_soil": n_soil,
                "pct_soil": round(100 * n_soil / len(ids), 1) if ids else 0.0,
            }
        )
    panel_c = pd.DataFrame(rows_c)
    panel_c.to_csv(args.output_dir / "shared_vs_exclusive_soil_strict.csv", index=False)

    manifest = {
        "schema_version": 1,
        "taxonomy_release": "ncbi-phylum-2026-08-04-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "searched_features": len(searched),
        "by_method": dict(Counter(method_of.get(f) for f in searched)),
        "match_resolution": counters,
        "panel_b": panel_b.to_dict("records"),
        "panel_c": panel_c.to_dict("records"),
    }
    (args.output_dir / "FIG2BC_MANIFEST.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({k: manifest[k] for k in ("searched_features", "by_method", "match_resolution")}, indent=1))
    print(panel_c.to_string(index=False))


if __name__ == "__main__":
    main()
