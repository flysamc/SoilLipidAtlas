#!/usr/bin/env python3
"""Build the uncached strict-IndVal fastMASST query package without submitting it."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from paper2_repro.scripts.run_local_diagnostic_annotation import sha256


TAXONOMY_RELEASE = "ncbi-phylum-2026-08-04-v1"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--indval", required=True, type=Path)
    parser.add_argument("--atlas", required=True, type=Path)
    parser.add_argument("--full-mgf", required=True, type=Path)
    parser.add_argument("--cache-dir", required=True, type=Path)
    parser.add_argument("--recovered-producer", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--refresh-cached",
        action="store_true",
        help="Stage every strict queryable feature so one current library index is used.",
    )
    return parser.parse_args()


def read_mgf_blocks(path: Path) -> dict[str, tuple[str, int]]:
    blocks: dict[str, tuple[str, int]] = {}
    current: list[str] = []
    feature_id: str | None = None
    peak_count = 0
    in_block = False
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.rstrip("\r\n") == "BEGIN IONS":
                if in_block:
                    raise ValueError("Nested BEGIN IONS record")
                current = [line]
                feature_id = None
                peak_count = 0
                in_block = True
                continue
            if not in_block:
                continue
            current.append(line)
            if line.startswith("FEATURE_ID="):
                feature_id = line.rstrip("\r\n").split("=", 1)[1]
            stripped = line.strip()
            if stripped and stripped[0] in "-0123456789" and "=" not in stripped:
                fields = stripped.split()
                if len(fields) >= 2:
                    try:
                        float(fields[0])
                        float(fields[1])
                        peak_count += 1
                    except ValueError:
                        pass
            if line.rstrip("\r\n") == "END IONS":
                if not feature_id:
                    raise ValueError("MGF record lacks FEATURE_ID")
                if feature_id in blocks:
                    raise ValueError(f"Duplicate MGF FEATURE_ID: {feature_id}")
                blocks[feature_id] = ("".join(current), peak_count)
                current = []
                feature_id = None
                in_block = False
    if in_block:
        raise ValueError("Unclosed MGF record")
    return blocks


def main() -> None:
    args = parse_args()
    indval_path = args.indval.resolve()
    atlas_path = args.atlas.resolve()
    full_mgf_path = args.full_mgf.resolve()
    cache_dir = args.cache_dir.resolve()
    recovered_producer = args.recovered_producer.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    indval = pd.read_csv(indval_path, low_memory=False)
    atlas = pd.read_csv(atlas_path, low_memory=False)
    for label, frame in (("IndVal", indval), ("atlas", atlas)):
        if "feature_id" not in frame:
            raise ValueError(f"{label} table lacks feature_id")
        if frame["feature_id"].duplicated().any():
            raise ValueError(f"{label} table contains duplicate feature IDs")

    strict_ids = set(indval["feature_id"].astype(str))
    cached_ids = {
        path.name.removesuffix("_matches.tsv")
        for path in cache_dir.glob("*_matches.tsv")
    }
    blocks = read_mgf_blocks(full_mgf_path)
    uncached_ids = strict_ids - cached_ids
    target_ids = strict_ids if args.refresh_cached else uncached_ids
    queryable_ids = {
        feature_id
        for feature_id in target_ids & set(blocks)
        if blocks[feature_id][1] >= 3
    }
    fewer_than_three_ids = {
        feature_id
        for feature_id in target_ids & set(blocks)
        if blocks[feature_id][1] < 3
    }
    no_ms2_ids = target_ids - set(blocks)
    unqueryable_ids = fewer_than_three_ids | no_ms2_ids

    atlas_index = atlas.set_index("feature_id", drop=False)
    missing_metadata = sorted((queryable_ids | unqueryable_ids) - set(atlas_index.index))
    if missing_metadata:
        raise KeyError(f"Missing atlas metadata for {len(missing_metadata)} features")

    query_mgf = output_dir / "figure2a_strict_indval_fastmasst_pending.mgf"
    query_csv = output_dir / "figure2a_strict_indval_fastmasst_pending.csv"
    no_ms2_csv = output_dir / "figure2a_strict_indval_fastmasst_unqueryable.csv"
    manifest_path = output_dir / "manifest.json"
    ordered_query_ids = sorted(queryable_ids)
    with query_mgf.open("w", encoding="utf-8", newline="\n") as handle:
        for feature_id in ordered_query_ids:
            handle.write(blocks[feature_id][0].rstrip("\r\n") + "\n\n")
    atlas_index.loc[ordered_query_ids].to_csv(query_csv, index=False)
    unqueryable_frame = atlas_index.loc[sorted(unqueryable_ids)].copy()
    unqueryable_frame["fastmasst_exclusion_reason"] = [
        "no_usable_ms2" if feature_id in no_ms2_ids else "fewer_than_3_fragment_peaks"
        for feature_id in unqueryable_frame.index
    ]
    unqueryable_frame.to_csv(no_ms2_csv, index=False)

    manifest = {
        "taxonomy_release": TAXONOMY_RELEASE,
        "scope": (
            "All strict queryable IndVal features staged for one current fastMASST index"
            if args.refresh_cached
            else "Uncached strict IndVal features staged for fastMASST; no queries submitted"
        ),
        "strict_indval_features": len(strict_ids),
        "cached_files_total": len(cached_ids),
        "cached_strict_indval_features": len(strict_ids & cached_ids),
        "uncached_strict_indval_features": len(uncached_ids),
        "refresh_cached": args.refresh_cached,
        "target_strict_indval_features": len(target_ids),
        "queryable_at_least_3_fragment_peaks": len(queryable_ids),
        "unqueryable_without_usable_ms2": len(no_ms2_ids),
        "unqueryable_fewer_than_3_fragment_peaks": len(fewer_than_three_ids),
        "submitted": False,
        "producer_status": "recovered_requires_resume_hardening",
        "inputs": {
            "indval": {"path": str(indval_path), "sha256": sha256(indval_path)},
            "atlas": {"path": str(atlas_path), "sha256": sha256(atlas_path)},
            "full_mgf": {"path": str(full_mgf_path), "sha256": sha256(full_mgf_path)},
            "cache_manifest": {
                "path": str(cache_dir),
                "cached_feature_ids_sha256": __import__("hashlib").sha256(
                    ("\n".join(sorted(cached_ids)) + "\n").encode("utf-8")
                ).hexdigest(),
            },
            "producer": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__).resolve())},
            "recovered_fastmasst_producer": {
                "path": str(recovered_producer),
                "sha256": sha256(recovered_producer),
            },
        },
        "outputs": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in (query_mgf, query_csv, no_ms2_csv)
        },
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in manifest.items() if key != "inputs"}, indent=2))


if __name__ == "__main__":
    main()
