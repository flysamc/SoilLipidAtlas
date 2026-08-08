#!/usr/bin/env python3
"""Validate and stage the 2026-08-05 positive-mode RT recovery archive.

This importer does not execute the recovered pickle.  It preserves the archive as
historical/method evidence and audits exact-feature compatibility with the strict
positive atlas.  Predictions are never transferred by mass or retention time alone.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[2]
PREFIX = "DIAGNOSTIC_RT_POS_RECOVERY_2026-08-05/"


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256(path: Path) -> str:
    return digest(path.read_bytes())


def read_csv_bytes(data: bytes) -> list[dict[str, str]]:
    return list(csv.DictReader(io.StringIO(data.decode("utf-8-sig"))))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument(
        "--strict-atlas",
        type=Path,
        default=ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/biomarker_discovery/atlas_pos_strict.csv",
    )
    parser.add_argument(
        "--harmonised",
        type=Path,
        default=ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/annotation/step2_harmonization/harmonised_annotations_pos.csv",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    archive = args.archive.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    with ZipFile(archive) as bundle:
        names = set(bundle.namelist())
        sums_member = PREFIX + "SHA256SUMS.txt"
        if sums_member not in names:
            raise ValueError("SHA256SUMS.txt is missing")
        expected: dict[str, str] = {}
        for line in bundle.read(sums_member).decode("utf-8").splitlines():
            value, relative = line.split("  ", 1)
            expected[relative] = value

        member_records = []
        checksums_ok = True
        for relative, expected_hash in sorted(expected.items()):
            member = PREFIX + relative
            if member not in names:
                raise ValueError(f"Checksummed archive member missing: {member}")
            safe = PurePosixPath(relative)
            if safe.is_absolute() or ".." in safe.parts:
                raise ValueError(f"Unsafe archive member path: {relative}")
            data = bundle.read(member)
            actual_hash = digest(data)
            checksums_ok &= actual_hash == expected_hash
            destination = output.joinpath(*safe.parts)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(data)
            member_records.append(
                {
                    "archive_member": member,
                    "path": str(destination.relative_to(ROOT)).replace("\\", "/"),
                    "bytes": len(data),
                    "sha256": actual_hash,
                    "expected_sha256": expected_hash,
                    "checksum_pass": actual_hash == expected_hash,
                }
            )
        # Preserve the checksum manifest itself; it intentionally cannot hash itself.
        sums_data = bundle.read(sums_member)
        (output / "SHA256SUMS.txt").write_bytes(sums_data)

    primary = read_csv(output / "historical_outputs/platinum_rt_predicted_sum_species.csv")
    hybrid = read_csv(output / "historical_outputs/platinum_rt_predicted_sum_species_hybrid.csv")
    v2 = read_csv(output / "historical_outputs/platinum_unified_annotations_v2.csv")
    final = read_csv(output / "historical_outputs/platinum_unified_annotations_final.csv")
    atlas = {row["feature_id"]: row for row in read_csv(args.strict_atlas.resolve())}
    harmonised = {row["feature_id"]: row for row in read_csv(args.harmonised.resolve())}

    primary_ids = {row["feature_id"] for row in primary}
    hybrid_ids = {row["feature_id"] for row in hybrid}
    prediction_ids = primary_ids | hybrid_ids
    prediction_by_id = {row["feature_id"]: row["predicted_sum_species"] for row in primary + hybrid}
    v2_by_id = {row["feature_id"]: row for row in v2}
    final_by_id = {row["feature_id"]: row for row in final}
    high = sum(row["upgrade_confidence"] == "high" for row in primary + hybrid)
    medium = sum(row["upgrade_confidence"] == "medium" for row in primary + hybrid)
    hybrid_over_60 = sum(float(row["predicted_total_carbon"]) > 60 for row in hybrid)

    ledger: list[dict[str, object]] = []
    for source, rows in (("primary", primary), ("hybrid", hybrid)):
        for row in rows:
            feature_id = row["feature_id"]
            strict = atlas.get(feature_id)
            if strict is None:
                continue
            current = harmonised.get(feature_id, {})
            current_class = (current.get("annotation_class_normalised") or "").strip()
            model_class = row["model_class"].strip()
            rt_delta = abs(float(row["observed_rt"]) - float(strict["consensus_rt"]))
            rt_exact = rt_delta <= 1e-9
            class_match = bool(current_class) and current_class.casefold() == model_class.casefold()
            if source == "hybrid":
                decision = "quarantine_hybrid_removed_from_historical_final"
            elif not rt_exact:
                decision = "reject_incompatible_rt_provenance"
            elif class_match:
                decision = "candidate_exact_id_rt_class_cache_reuse_review_only"
            elif current_class:
                decision = "reject_current_class_mismatch"
            else:
                decision = "review_missing_current_class_compatibility"
            ledger.append(
                {
                    "feature_id": feature_id,
                    "historical_output": source,
                    "model_class": model_class,
                    "current_harmonised_class": current_class,
                    "observed_rt_historical": row["observed_rt"],
                    "strict_consensus_rt": strict["consensus_rt"],
                    "rt_abs_delta_min": f"{rt_delta:.12g}",
                    "rt_exact": str(rt_exact).lower(),
                    "class_exact": str(class_match).lower(),
                    "predicted_sum_species": row["predicted_sum_species"],
                    "upgrade_confidence": row["upgrade_confidence"],
                    "decision": decision,
                }
            )
    fields = [
        "feature_id", "historical_output", "model_class", "current_harmonised_class",
        "observed_rt_historical", "strict_consensus_rt", "rt_abs_delta_min", "rt_exact",
        "class_exact", "predicted_sum_species", "upgrade_confidence", "decision",
    ]
    write_csv(output / "exact_id_cache_compatibility.csv", ledger, fields)
    decision_counts: dict[str, int] = {}
    source_counts: dict[str, int] = {}
    for row in ledger:
        decision_counts[str(row["decision"])] = decision_counts.get(str(row["decision"]), 0) + 1
        source_counts[str(row["historical_output"])] = source_counts.get(str(row["historical_output"]), 0) + 1

    checks = {
        "all_internal_checksums_pass": checksums_ok,
        "checksummed_member_count_is_23": len(member_records) == 23,
        "primary_unique_rows_are_980": len(primary) == len(primary_ids) == 980,
        "primary_high_medium_are_760_220": (
            sum(r["upgrade_confidence"] == "high" for r in primary) == 760
            and sum(r["upgrade_confidence"] == "medium" for r in primary) == 220
        ),
        "hybrid_unique_rows_are_784": len(hybrid) == len(hybrid_ids) == 784,
        "hybrid_high_medium_are_215_569": (
            sum(r["upgrade_confidence"] == "high" for r in hybrid) == 215
            and sum(r["upgrade_confidence"] == "medium" for r in hybrid) == 569
        ),
        "primary_hybrid_disjoint_union_is_1764": not (primary_ids & hybrid_ids) and len(prediction_ids) == 1764,
        "combined_high_medium_are_975_789": high == 975 and medium == 789,
        "hybrid_carbon_over_60_is_403": hybrid_over_60 == 403,
        "v2_applies_all_1764_predictions": all(
            v2_by_id.get(feature_id, {}).get("rt_predicted_species", "") == predicted
            for feature_id, predicted in prediction_by_id.items()
        ),
        "final_retains_all_980_primary_predictions": all(
            final_by_id.get(feature_id, {}).get("rt_predicted_species", "") == prediction_by_id[feature_id]
            for feature_id in primary_ids
        ),
        "final_removes_all_784_hybrid_predictions": all(
            not final_by_id.get(feature_id, {}).get("rt_predicted_species", "").strip()
            for feature_id in hybrid_ids
        ),
        "strict_overlap_is_276": len(ledger) == 276,
        "strict_primary_overlap_is_167": source_counts.get("primary") == 167,
        "strict_hybrid_overlap_is_109": source_counts.get("hybrid") == 109,
        "all_strict_overlap_observed_rt_exact": all(r["rt_exact"] == "true" for r in ledger),
    }
    manifest = {
        "schema_version": 1,
        "stage_id": "positive_rt_sum_composition_recovery_import",
        "taxonomy_release": "ncbi-phylum-2026-08-04-v1",
        "polarity": "POS",
        "status": "pass_recovery_evidence_only_producer_blocked" if all(checks.values()) else "fail",
        "source_archive": {
            "original_path": str(archive),
            "bytes": archive.stat().st_size,
            "sha256": sha256(archive),
        },
        "strict_inputs": {
            "atlas": {
                "path": str(args.strict_atlas.resolve().relative_to(ROOT)).replace("\\", "/"),
                "bytes": args.strict_atlas.resolve().stat().st_size,
                "sha256": sha256(args.strict_atlas.resolve()),
            },
            "harmonised_annotations": {
                "path": str(args.harmonised.resolve().relative_to(ROOT)).replace("\\", "/"),
                "bytes": args.harmonised.resolve().stat().st_size,
                "sha256": sha256(args.harmonised.resolve()),
            },
        },
        "historical_reconciliation": {
            "primary_rows": len(primary),
            "hybrid_rows": len(hybrid),
            "disjoint_union_rows": len(prediction_ids),
            "high": high,
            "medium": medium,
            "hybrid_carbon_over_60": hybrid_over_60,
            "v2_rows": len(v2),
            "final_rows": len(final),
            "historical_final_retains_primary_predictions": 980,
            "historical_final_removes_hybrid_predictions": 784,
            "hybrid_carbon_at_most_60_also_removed": 784 - hybrid_over_60,
        },
        "strict_exact_id_audit": {
            "strict_atlas_features": len(atlas),
            "exact_prediction_id_overlap": len(ledger),
            "primary_overlap": source_counts.get("primary", 0),
            "hybrid_overlap": source_counts.get("hybrid", 0),
            "decision_counts": decision_counts,
            "ledger": str((output / "exact_id_cache_compatibility.csv").relative_to(ROOT)).replace("\\", "/"),
            "ledger_sha256": sha256(output / "exact_id_cache_compatibility.csv"),
        },
        "checks": checks,
        "members": member_records,
        "producer_gate": {
            "exact_pos_producer_recovered": False,
            "blocking_missing_contract": [
                "five model feature names, order, and transformations",
                "exact training dataframe and filters",
                "adduct correction and candidate generation code",
                "primary versus hybrid orchestration",
                "atlas application and downgrade code",
                "original command and matching Python environment",
            ],
        },
        "decision": [
            "Treat 1,764 and 403 as reconciled historical frozen-output counts, not a strict rerun.",
            "Do not identify step13_rt_family_modeling.py or the NEG analogue as the missing POS producer.",
            "Do not execute or reuse a prediction unless feature ID, observed RT, class, and provenance are compatible.",
            "Quarantine all 784 hybrid predictions because the historical final table removed them.",
            "Keep any exact-compatible primary overlap review-only until the missing producer is recovered or a new versioned method is approved.",
            "Do not update manuscript, figures, submitted documents, or supplementary tables from this stage.",
        ],
    }
    (output / "stage_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if manifest["status"] == "fail":
        raise SystemExit("POS RT recovery validation failed")
    print(json.dumps({"status": manifest["status"], "historical": manifest["historical_reconciliation"], "strict": manifest["strict_exact_id_audit"]}, indent=2))


if __name__ == "__main__":
    main()
