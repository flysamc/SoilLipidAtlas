#!/usr/bin/env python3
"""Validate and freeze the strict Figure 2A feature-selection release."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    atlas_dir = args.atlas_dir.resolve()
    summary = json.loads((atlas_dir / "summary.json").read_text(encoding="utf-8"))
    reproducibility = json.loads(
        (atlas_dir / "reproducibility_check.json").read_text(encoding="utf-8")
    )
    annotation = json.loads(
        (atlas_dir / "annotation_recovery" / "audit_report.json").read_text(
            encoding="utf-8"
        )
    )

    checks = {
        "historical_reproduction_pass": bool(summary["reproduction"]["pass"]),
        "archlips_reproduction_pass": summary["archlips_assignment_reproduction_accuracy"] == 1.0,
        "deterministic_scientific_outputs": bool(reproducibility["scientific_outputs_pass"]),
        "taxonomy_release_match": summary["taxonomy_release"] == annotation["taxonomy_release"],
        "sixteen_analysis_phyla": summary["validation"]["expected_phyla"] == summary["validation"]["observed_phyla"] == 16,
        "unique_feature_ids": summary["validation"]["duplicate_feature_ids"] == 0,
        "queue_partition_matches": (
            summary["annotation_coverage"]["new_need_full_annotation_pipeline"]
            == annotation["queue_unique_features"]
            == annotation["any_recovered_annotation_evidence"] + annotation["needs_annotation_pipeline"]
        ),
        "queue_ms2_partition_matches": (
            annotation["queue_unique_features"]
            == annotation["queue_with_recovered_consensus_ms2_spectrum"]
            + annotation["queue_without_recovered_consensus_ms2_spectrum"]
        ),
    }
    if not all(checks.values()):
        raise RuntimeError(f"Biomarker freeze checks failed: {checks}")

    feature_files = [
        "atlas_pos_strict.csv",
        "composite_platinum_strict.csv",
        "composite_silver_strict.csv",
        "indval_pairs_strict.csv",
        "indval_unique_strict.csv",
        "archlips_assignments_strict.csv",
        "figure2a_tier_counts.csv",
        "biomarker_counts_by_phylum.csv",
        "discovery_method_by_phylum.csv",
        "annotation_queue_new_features.csv",
        "historical_control/composite_platinum_reproduced.csv",
        "historical_control/composite_silver_reproduced.csv",
        "historical_control/indval_pairs_reproduced.csv",
        "historical_control/indval_unique_reproduced.csv",
        "historical_control/reproduction_gate.json",
        "summary.json",
        "reproducibility_check.json",
        "annotation_recovery/audit_report.json",
        "annotation_recovery/annotation_queue_recovered_evidence.csv",
        "annotation_recovery/annotation_queue_unresolved.csv",
    ]
    outputs = []
    for relative in feature_files:
        path = atlas_dir / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        outputs.append(
            {
                "path": relative,
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )

    manifest = {
        "taxonomy_release": summary["taxonomy_release"],
        "status": "feature_set_frozen_annotation_pending",
        "checks": checks,
        "counts": {
            "biomarkers": summary["validation"]["atlas_rows"],
            "analysis_phyla": summary["validation"]["observed_phyla"],
            "composite_platinum": summary["strict_selection"]["composite_platinum"],
            "composite_silver": summary["strict_selection"]["composite_silver"],
            "indval_pair_rows": summary["strict_selection"]["indval_pair_rows"],
            "indval_unique": summary["strict_selection"]["indval_unique"],
            "new_feature_queue": annotation["queue_unique_features"],
            "recovered_annotation_evidence": annotation["any_recovered_annotation_evidence"],
            "unresolved_annotation_queue": annotation["needs_annotation_pipeline"],
            "queue_with_recovered_ms2": annotation["queue_with_recovered_consensus_ms2_spectrum"],
            "queue_without_recovered_ms2": annotation["queue_without_recovered_consensus_ms2_spectrum"],
        },
        "outputs": outputs,
    }
    output_path = atlas_dir / "freeze_manifest.json"
    output_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: value for key, value in manifest.items() if key != "outputs"}, indent=2))


if __name__ == "__main__":
    main()
