#!/usr/bin/env python3
"""Restore success statuses from validated durable compact fastMASST files."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
from datetime import datetime, timezone
from pathlib import Path


SCHEMA = "soilmasst-fastmasst-minimal-v1"
FIELDS = ("USI", "Dataset", "Cosine", "Matching Peaks", "Delta Mass")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def validate_matches(path: Path) -> tuple[int, int, float]:
    rows = 0
    datasets: set[str] = set()
    top_cosine = 0.0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != FIELDS:
            raise ValueError(f"Unexpected compact schema in {path}: {reader.fieldnames}")
        for row in reader:
            rows += 1
            dataset = str(row.get("Dataset", ""))
            if dataset:
                datasets.add(dataset)
            try:
                top_cosine = max(top_cosine, float(row.get("Cosine", "")))
            except (TypeError, ValueError):
                pass
    return rows, len(datasets), top_cosine


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--label", required=True)
    args = parser.parse_args()

    run_dir = args.run_dir.resolve()
    statuses = run_dir / "statuses"
    matches = run_dir / "per_feature"
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    backup = run_dir / f"status_reconciliation_backup_{stamp}"
    recovered: list[dict] = []
    already_success = 0
    invalid: list[dict] = []

    for match_path in sorted(matches.glob("*_matches.tsv.gz")):
        feature_id = match_path.name[: -len("_matches.tsv.gz")]
        status_path = statuses / f"{feature_id}.json"
        existing = read_json(status_path) if status_path.exists() else {}
        if existing.get("state") == "success":
            if existing.get("matches_sha256") != sha256(match_path):
                invalid.append({"feature_id": feature_id, "reason": "success_checksum_mismatch"})
            else:
                already_success += 1
            continue
        try:
            n_matches, n_datasets, top_cosine = validate_matches(match_path)
        except Exception as exc:
            invalid.append({"feature_id": feature_id, "reason": f"{type(exc).__name__}:{exc}"})
            continue
        backup.mkdir(parents=True, exist_ok=True)
        if status_path.exists():
            shutil.copy2(status_path, backup / status_path.name)
        payload = {
            "feature_id": feature_id,
            "state": "success",
            "task_id": existing.get("task_id"),
            "result_task": existing.get("result_task"),
            "submitted_utc": existing.get("submitted_utc"),
            "completed_utc": datetime.fromtimestamp(
                match_path.stat().st_mtime, timezone.utc
            ).isoformat(),
            "request_sha256": existing.get("request_sha256"),
            "n_matches": n_matches,
            "n_datasets": n_datasets,
            "top_cosine": round(top_cosine, 6),
            "matches_file": str(match_path.resolve()),
            "matches_sha256": sha256(match_path),
            "match_storage_schema": SCHEMA,
            "retained_match_fields": list(FIELDS),
            "reconciliation": "restored from validated durable compact match file",
        }
        atomic_json(status_path, payload)
        recovered.append(
            {
                "feature_id": feature_id,
                "prior_state": existing.get("state", "missing"),
                "matches_sha256": payload["matches_sha256"],
                "n_matches": n_matches,
            }
        )

    state_counts: dict[str, int] = {}
    for status_path in statuses.glob("*.json"):
        state = str(read_json(status_path).get("state", "missing"))
        state_counts[state] = state_counts.get(state, 0) + 1
    manifest = {
        "schema_version": 1,
        "label": args.label,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "run_dir": str(run_dir),
        "durable_match_files": len(list(matches.glob("*_matches.tsv.gz"))),
        "already_valid_success": already_success,
        "restored_success": len(recovered),
        "invalid_match_files": invalid,
        "post_reconciliation_states": state_counts,
        "backup_dir": str(backup) if backup.exists() else None,
        "recovered": recovered,
    }
    manifest_path = run_dir / f"status_reconciliation_{stamp}.json"
    atomic_json(manifest_path, manifest)
    print(json.dumps({k: v for k, v in manifest.items() if k != "recovered"}, indent=2))
    if invalid:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
