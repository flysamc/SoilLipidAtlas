#!/usr/bin/env python3
"""Compact completed async fastMASST match files to the release storage contract."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from pathlib import Path

from paper2_repro.scripts.run_fastmasst_async import (
    MATCH_STORAGE_SCHEMA,
    RETAINED_MATCH_FIELDS,
    atomic_json,
    compact_match_row,
    read_json,
    sha256,
)


def iter_stored_rows(path: Path):
    if path.name.endswith(".jsonl.gz"):
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                if line.strip():
                    yield json.loads(line)
        return
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def compact_file(source: Path, destination: Path) -> int:
    return compact_file_from_rows(iter_stored_rows(source), destination)


def compact_file_from_rows(rows, destination: Path) -> int:
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    count = 0
    with gzip.open(temporary, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RETAINED_MATCH_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(compact_match_row(row))
            count += 1
    os.replace(temporary, destination)
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    run_dir = args.run_dir.resolve()
    statuses_dir = run_dir / "statuses"
    per_feature_dir = run_dir / "per_feature"
    migrated = 0
    removed_bytes = 0
    for status_path in sorted(statuses_dir.glob("*.json")):
        status = read_json(status_path)
        if not status or status.get("state") != "success":
            continue
        source_value = status.get("matches_file")
        source = Path(source_value) if source_value else None
        if source is not None and not source.is_absolute():
            source = (Path.cwd() / source).resolve()
        if source is None or not source.is_file():
            candidates = sorted(per_feature_dir.glob(f"{status_path.stem}_matches.*"))
            if not candidates:
                if int(status.get("n_matches", -1)) == 0:
                    destination = per_feature_dir / f"{status_path.stem}_matches.tsv.gz"
                    compact_file_from_rows([], destination)
                    status.update(
                        {
                            "matches_file": str(destination),
                            "matches_sha256": sha256(destination),
                            "match_storage_schema": MATCH_STORAGE_SCHEMA,
                            "retained_match_fields": list(RETAINED_MATCH_FIELDS),
                        }
                    )
                    atomic_json(status_path, status)
                    migrated += 1
                    continue
                raise FileNotFoundError(f"No match file for {status_path.stem}")
            source = candidates[0]
        if status.get("match_storage_schema") == MATCH_STORAGE_SCHEMA:
            continue

        original_bytes = source.stat().st_size
        original_sha256 = sha256(source)
        destination = per_feature_dir / f"{status_path.stem}_matches.tsv.gz"
        row_count = compact_file(source, destination)
        if row_count != int(status.get("n_matches", row_count)):
            destination.unlink(missing_ok=True)
            raise ValueError(
                f"Row-count mismatch for {status_path.stem}: {row_count} != {status.get('n_matches')}"
            )
        compact_sha256 = sha256(destination)
        status.update(
            {
                "matches_file": str(destination),
                "matches_sha256": compact_sha256,
                "match_storage_schema": MATCH_STORAGE_SCHEMA,
                "retained_match_fields": list(RETAINED_MATCH_FIELDS),
                "pre_compaction_matches_sha256": original_sha256,
                "pre_compaction_matches_bytes": original_bytes,
            }
        )
        atomic_json(status_path, status)
        if source != destination and source.exists():
            removed_bytes += source.stat().st_size
            source.unlink()
        migrated += 1
    print(json.dumps({"migrated": migrated, "removed_redundant_bytes": removed_bytes}, indent=2))


if __name__ == "__main__":
    main()
