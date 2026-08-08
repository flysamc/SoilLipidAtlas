#!/usr/bin/env python3
"""Finalize strict fastMASST and Pan-ReDU summaries after both search gates pass."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

RELEASE = "ncbi-phylum-2026-08-04-v1"


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


def run_pan_redu(
    run_root: Path, polarity: str, soil_metadata: Path
) -> Path:
    if polarity == "POS":
        feature_metadata = run_root / "biomarker_discovery/indval_unique_strict.csv"
        status_dir = run_root / (
            "biomarker_discovery/external_annotation_package/"
            "fastmasst_async_full_20260806/statuses"
        )
        pan_dir = run_root / "biomarker_discovery/external_annotation_package/pan_redu_current_index_pos"
        extra = [
            "--exclusions",
            str(
                run_root
                / (
                    "biomarker_discovery/external_annotation_package/fastmasst_all_current_index/"
                    "figure2a_strict_indval_fastmasst_unqueryable.csv"
                )
            ),
        ]
        expected_total, expected_queryable = 7751, 7408
    else:
        feature_metadata = run_root / "annotation_recovery_neg/strict_feature_spectrum_ledger.csv"
        status_dir = run_root / "annotation_recovery_neg/fastmasst_async_full_20260806/statuses"
        pan_dir = run_root / "annotation_recovery_neg/pan_redu_current_index"
        extra = ["--usable-column", "has_usable_ms2"]
        expected_total, expected_queryable = 5697, 5695

    command = [
        sys.executable,
        "-m",
        "paper2_repro.scripts.build_pan_redu_strict",
        "--polarity",
        polarity,
        "--feature-metadata",
        str(feature_metadata),
        "--status-dir",
        str(status_dir),
        "--soil-metadata",
        str(soil_metadata),
        "--output-dir",
        str(pan_dir),
        "--expected-total",
        str(expected_total),
        "--expected-queryable",
        str(expected_queryable),
        *extra,
    ]
    completed = subprocess.run(command)
    if completed.returncode != 0:
        raise RuntimeError(f"{polarity} Pan-ReDU readiness/finalization gate failed")
    return pan_dir


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--soil-metadata", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    run_root = args.run_root.resolve()
    soil_metadata = args.soil_metadata.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    pan_dirs = {
        polarity: run_pan_redu(run_root, polarity, soil_metadata)
        for polarity in ("POS", "NEG")
    }
    fast_rows = []
    pan_rows = []
    for polarity, directory in pan_dirs.items():
        features = pd.read_csv(directory / "pan_redu_feature_summary.csv", low_memory=False)
        groups = pd.read_csv(directory / "pan_redu_group_summary.csv", low_memory=False)
        fast_rows.append(
            {
                "mode": polarity,
                "strict_features": len(features),
                "queryable_features": int(features.queryable.sum()),
                "excluded_no_usable_ms2": int((~features.queryable).sum()),
                "features_with_public_match": int((features.n_total_matches > 0).sum()),
                "features_with_soil_match": int((features.n_soil_matches > 0).sum()),
                "total_public_matches": int(features.n_total_matches.sum()),
                "total_soil_matches": int(features.n_soil_matches.sum()),
            }
        )
        pan_rows.append(groups.assign(mode=polarity))

    fast_path = output_dir / "fastmasst_summary.csv"
    pan_path = output_dir / "panredu_summary.csv"
    write_csv_atomic(pd.DataFrame(fast_rows), fast_path)
    write_csv_atomic(pd.concat(pan_rows, ignore_index=True), pan_path)
    manifest = {
        "schema_version": 1,
        "taxonomy_release": RELEASE,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "producer": str(Path(__file__).resolve()),
        "producer_sha256": sha256(Path(__file__).resolve()),
        "database": "metabolomicspanrepo_index_latest",
        "status": "complete_both_polarities_one_current_index",
        "outputs": {
            path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in (fast_path, pan_path)
        },
    }
    temporary = output_dir / "PUBLIC_VALIDATION_MANIFEST.json.tmp"
    temporary.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, output_dir / "PUBLIC_VALIDATION_MANIFEST.json")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
