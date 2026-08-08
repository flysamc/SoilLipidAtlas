#!/usr/bin/env python3
"""Run the recovered fastMASST query semantics with durable per-feature checkpoints."""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import re
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mgf", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--recovered-producer", required=True, type=Path)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--pm-tolerance", type=float, default=0.05)
    parser.add_argument("--fragment-tolerance", type=float, default=0.05)
    parser.add_argument("--min-cos", type=float, default=0.7)
    parser.add_argument(
        "--database",
        help="Override the recovered producer database while retaining its request semantics",
    )
    parser.add_argument("--retry-errors", action="store_true")
    parser.add_argument("--limit", type=int)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_producer(path: Path):
    spec = importlib.util.spec_from_file_location("recovered_fastmasst", path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, data: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def atomic_tsv(path: Path, rows: list[dict]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    pd.DataFrame(rows).to_csv(temporary, sep="\t", index=False)
    os.replace(temporary, path)


def read_status(path: Path) -> dict | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def parse_matches(payload) -> list[dict]:
    if isinstance(payload, dict) and "results" in payload:
        matches = payload["results"]
    elif isinstance(payload, list):
        matches = payload
    else:
        matches = []
    return [row for row in matches if isinstance(row, dict)]


def summarize_matches(matches: list[dict]) -> tuple[int, int, float]:
    datasets: set[str] = set()
    top_cosine = 0.0
    for match in matches:
        usi = str(match.get("usi", ""))
        dataset = str(match.get("dataset_id", ""))
        if not dataset and ":" in usi:
            dataset = usi.split(":", 2)[1]
        if dataset:
            datasets.add(dataset)
        try:
            top_cosine = max(top_cosine, float(match.get("cosine", match.get("score", 0))))
        except (TypeError, ValueError):
            pass
    return len(matches), len(datasets), top_cosine


def main() -> None:
    args = parse_args()
    mgf_path = args.mgf.resolve()
    producer_path = args.recovered_producer.resolve()
    out_dir = args.out_dir.resolve()
    statuses_dir = out_dir / "statuses"
    matches_dir = out_dir / "per_feature"
    statuses_dir.mkdir(parents=True, exist_ok=True)
    matches_dir.mkdir(parents=True, exist_ok=True)

    if args.workers < 1 or args.workers > 10:
        raise ValueError("workers must be between 1 and 10")
    producer = load_producer(producer_path)
    recovered_database = producer.DATABASE
    if args.database:
        producer.DATABASE = args.database
    spectra = producer.parse_mgf(str(mgf_path), min_peaks=3)
    tasks = []
    skipped_success = 0
    skipped_error = 0
    for spectrum in spectra:
        feature_id = str(spectrum.get("FEATURE_ID", spectrum.get("SCANS", "unknown")))
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", feature_id):
            raise ValueError(f"Unsafe feature ID: {feature_id}")
        existing = read_status(statuses_dir / f"{feature_id}.json")
        if existing and existing.get("state") == "success":
            skipped_success += 1
            continue
        if existing and existing.get("state") == "error" and not args.retry_errors:
            skipped_error += 1
            continue
        precursor_mz = float(str(spectrum.get("PEPMASS", "0")).split()[0])
        charge = int(str(spectrum.get("CHARGE", "1")).replace("+", "").replace("-", ""))
        tasks.append((feature_id, precursor_mz, charge, spectrum["peaks"]))
    if args.limit is not None:
        if args.limit < 1:
            raise ValueError("limit must be positive")
        tasks = tasks[: args.limit]

    run_manifest = {
        "started_utc": utc_now(),
        "input_mgf": str(mgf_path),
        "input_sha256": sha256(mgf_path),
        "input_eligible_spectra": len(spectra),
        "recovered_producer": str(producer_path),
        "recovered_producer_sha256": sha256(producer_path),
        "wrapper": str(Path(__file__).resolve()),
        "wrapper_sha256": sha256(Path(__file__).resolve()),
        "endpoint": producer.FASST_URL,
        "database": producer.DATABASE,
        "recovered_database": recovered_database,
        "database_overridden": producer.DATABASE != recovered_database,
        "parameters": {
            "workers": args.workers,
            "pm_tolerance": args.pm_tolerance,
            "fragment_tolerance": args.fragment_tolerance,
            "cosine_threshold": args.min_cos,
            "retry_errors": args.retry_errors,
            "limit": args.limit,
        },
        "skipped_existing_success": skipped_success,
        "skipped_existing_error": skipped_error,
        "scheduled_queries": len(tasks),
        "complete": False,
    }
    atomic_json(out_dir / "run_manifest.json", run_manifest)

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(
                producer.search_spectrum,
                feature_id,
                precursor_mz,
                charge,
                peaks,
                args.pm_tolerance,
                args.fragment_tolerance,
                args.min_cos,
            ): feature_id
            for feature_id, precursor_mz, charge, peaks in tasks
        }
        progress = tqdm(as_completed(futures), total=len(futures), desc="fastMASST")
        for future in progress:
            feature_id = futures[future]
            try:
                returned_id, payload, error = future.result()
                if returned_id != feature_id:
                    raise RuntimeError(f"Feature mismatch: {feature_id} != {returned_id}")
            except Exception as exc:  # preserve a durable error checkpoint
                payload = None
                error = f"worker_exception:{type(exc).__name__}:{str(exc)[:200]}"

            status_path = statuses_dir / f"{feature_id}.json"
            if error:
                atomic_json(
                    status_path,
                    {
                        "feature_id": feature_id,
                        "state": "error",
                        "error": str(error),
                        "completed_utc": utc_now(),
                    },
                )
                continue

            matches = parse_matches(payload)
            n_matches, n_datasets, top_cosine = summarize_matches(matches)
            matches_path = matches_dir / f"{feature_id}_matches.tsv"
            if matches:
                atomic_tsv(matches_path, matches)
            atomic_json(
                status_path,
                {
                    "feature_id": feature_id,
                    "state": "success",
                    "n_matches": n_matches,
                    "n_datasets": n_datasets,
                    "top_cosine": round(top_cosine, 6),
                    "matches_file": str(matches_path) if matches else None,
                    "completed_utc": utc_now(),
                },
            )

    statuses = []
    for path in sorted(statuses_dir.glob("*.json")):
        status = read_status(path)
        if status:
            statuses.append(status)
    summary_path = out_dir / "masst_summary.csv"
    temporary_summary = summary_path.with_suffix(".csv.tmp")
    pd.DataFrame(statuses).sort_values("feature_id").to_csv(temporary_summary, index=False)
    os.replace(temporary_summary, summary_path)

    errors = sum(status.get("state") == "error" for status in statuses)
    run_manifest.update(
        {
            "finished_utc": utc_now(),
            "status_records": len(statuses),
            "success_records": sum(status.get("state") == "success" for status in statuses),
            "error_records": errors,
            "summary_sha256": sha256(summary_path),
            "complete": errors == 0,
        }
    )
    atomic_json(out_dir / "run_manifest.json", run_manifest)
    print(json.dumps(run_manifest, indent=2))
    if errors:
        sys.exit(2)


if __name__ == "__main__":
    main()
