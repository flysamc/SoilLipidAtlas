#!/usr/bin/env python3
"""Summarize fastMASST pilot outcomes and gate bulk spectrum submission."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--taxonomy-release", required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    pilots = []
    input_hashes: set[str] = set()
    for pilot_dir_arg in args.pilot:
        pilot_dir = pilot_dir_arg.resolve()
        manifest_path = pilot_dir / "run_manifest.json"
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        status_paths = sorted((pilot_dir / "statuses").glob("*.json"))
        statuses = [json.loads(path.read_text(encoding="utf-8")) for path in status_paths]
        input_hashes.add(str(manifest.get("input_sha256")))
        pilots.append(
            {
                "directory": str(pilot_dir),
                "database": manifest.get("database"),
                "recovered_database": manifest.get("recovered_database", manifest.get("database")),
                "database_overridden": manifest.get("database_overridden", False),
                "endpoint": manifest.get("endpoint"),
                "api_root": manifest.get("api_root"),
                "protocol": manifest.get("protocol"),
                "scheduled_queries": manifest.get("scheduled_queries"),
                "success_records": manifest.get("success_records"),
                "error_records": manifest.get("error_records"),
                "complete": manifest.get("complete"),
                "errors": sorted({str(status.get("error")) for status in statuses if status.get("error")}),
                "input_sha256": manifest.get("input_sha256"),
                "manifest_sha256": sha256(manifest_path),
                "summary_sha256": sha256(pilot_dir / "masst_summary.csv"),
            }
        )

    for pilot in pilots:
        pilot["pilot_passed"] = bool(
            (pilot.get("scheduled_queries") or 0) > 0
            and pilot.get("success_records") == pilot.get("scheduled_queries")
            and (pilot.get("error_records") or 0) == 0
            and pilot.get("complete") is True
        )
    any_success = any((pilot.get("success_records") or 0) > 0 for pilot in pilots)
    all_pilots_passed = all(pilot["pilot_passed"] for pilot in pilots)
    all_same_input = len(input_hashes) == 1
    report = {
        "taxonomy_release": args.taxonomy_release,
        "scope": "Strict Figure 2A fastMASST external-service pilot gate",
        "pilot_count": len(pilots),
        "all_same_input": all_same_input,
        "service_query_succeeded": any_success,
        "all_pilots_passed": all_pilots_passed,
        "bulk_submission_safe": any_success and all_same_input and all_pilots_passed,
        "decision": (
            "Pilot succeeded; bulk submission may be considered."
            if any_success and all_same_input and all_pilots_passed
            else "Do not submit the remaining spectra until a pilot succeeds."
        ),
        "pilots": pilots,
    }
    args.output.resolve().write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
