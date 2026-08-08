#!/usr/bin/env python3
"""Run strict-release fastMASST queries through the asynchronous GNPS API.

The recovered SoilMass producer used the former synchronous form endpoint.  This
wrapper preserves its spectrum normalisation and search parameters, while adding
durable submit/status/result checkpoints for ``https://api.fasst.gnps2.org``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import os
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd
import requests
from tqdm import tqdm


TAXONOMY_RELEASE = "ncbi-phylum-2026-08-04-v1"
DEFAULT_API_ROOT = "https://api.fasst.gnps2.org"
DEFAULT_DATABASE = "metabolomicspanrepo_index_latest"
TERMINAL_FAILURE_STATES = {"FAILED", "FAILURE", "ERROR", "REVOKED"}
RETAINED_MATCH_FIELDS = ("USI", "Dataset", "Cosine", "Matching Peaks", "Delta Mass")
MATCH_STORAGE_SCHEMA = "soilmasst-fastmasst-minimal-v1"


class TaskNotFoundError(RuntimeError):
    """The API no longer recognises a previously checkpointed task ID."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mgf", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--recovered-producer", required=True, type=Path)
    parser.add_argument("--api-root", default=DEFAULT_API_ROOT)
    parser.add_argument("--database", default=DEFAULT_DATABASE)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--pm-tolerance", type=float, default=0.05)
    parser.add_argument("--fragment-tolerance", type=float, default=0.05)
    parser.add_argument("--min-cos", type=float, default=0.7)
    parser.add_argument("--poll-interval", type=float, default=2.0)
    parser.add_argument("--task-timeout", type=float, default=300.0)
    parser.add_argument("--http-timeout", type=float, default=30.0)
    parser.add_argument("--max-submit-retries", type=int, default=3)
    parser.add_argument("--limit", type=int)
    parser.add_argument(
        "--selection",
        choices=("first", "evenly-spaced"),
        default="first",
        help="How --limit selects spectra; evenly-spaced is preferred for pilots.",
    )
    parser.add_argument("--retry-errors", action="store_true")
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_sha256(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def load_producer(path: Path):
    spec = importlib.util.spec_from_file_location("recovered_fastmasst_async", path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def atomic_json(path: Path, data: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def read_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None


def select_spectra(spectra: list[dict], limit: int | None, strategy: str) -> list[dict]:
    if limit is None or limit >= len(spectra):
        return list(spectra)
    if limit < 1:
        raise ValueError("limit must be positive")
    if strategy == "first" or limit == 1:
        return list(spectra[:limit])
    indices = [(position * (len(spectra) - 1)) // (limit - 1) for position in range(limit)]
    if len(set(indices)) != limit:
        raise RuntimeError("evenly-spaced selection produced duplicate indices")
    return [spectra[index] for index in indices]


def spectrum_feature_id(spectrum: dict) -> str:
    feature_id = str(spectrum.get("FEATURE_ID", spectrum.get("SCANS", "unknown")))
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", feature_id):
        raise ValueError(f"Unsafe feature ID: {feature_id}")
    return feature_id


def build_query_payload(
    producer,
    spectrum: dict,
    database: str,
    pm_tolerance: float,
    fragment_tolerance: float,
    cosine_threshold: float,
) -> dict[str, Any]:
    filtered_peaks = producer.normalize_and_filter_peaks(spectrum["peaks"])
    if len(filtered_peaks) < 3:
        raise ValueError("too_few_peaks_after_filter")
    precursor_mz = float(str(spectrum.get("PEPMASS", "0")).split()[0])
    charge_text = str(spectrum.get("CHARGE", "1")).strip()
    charge_magnitude = int(charge_text.replace("+", "").replace("-", ""))
    charge = -charge_magnitude if "-" in charge_text else charge_magnitude
    query_spectrum = {
        "n_peaks": len(filtered_peaks),
        "peaks": filtered_peaks,
        "precursor_mz": precursor_mz,
        "precursor_charge": charge,
    }
    return {
        "library": database,
        "query_spectrum": json.dumps(query_spectrum, separators=(",", ":")),
        "analog": "No",
        "cache": "Yes",
        "pm_tolerance": pm_tolerance,
        "fragment_tolerance": fragment_tolerance,
        "cosine_threshold": cosine_threshold,
    }


def parse_matches(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, dict) and isinstance(payload.get("results"), list):
        return [row for row in payload["results"] if isinstance(row, dict)]
    if isinstance(payload, list):
        return [row for row in payload if isinstance(row, dict)]
    return []


def summarize_matches(matches: list[dict[str, Any]]) -> tuple[int, int, float]:
    datasets: set[str] = set()
    top_cosine = 0.0
    for match in matches:
        usi = str(match.get("USI", match.get("usi", "")))
        dataset = str(
            match.get("Dataset", match.get("dataset_id", match.get("dataset", "")))
        )
        if not dataset and ":" in usi:
            dataset = usi.split(":", 2)[1]
        if dataset:
            datasets.add(dataset)
        value = match.get("Cosine", match.get("cosine", match.get("score", 0)))
        try:
            top_cosine = max(top_cosine, float(value))
        except (TypeError, ValueError):
            pass
    return len(matches), len(datasets), top_cosine


def request_json(
    session: requests.Session,
    method: str,
    url: str,
    http_timeout: float,
    **kwargs,
) -> dict[str, Any]:
    response = session.request(method, url, timeout=http_timeout, **kwargs)
    response.raise_for_status()
    payload = response.json()
    if not isinstance(payload, dict):
        raise ValueError(f"Expected JSON object from {url}")
    return payload


def submit_query(
    session: requests.Session,
    api_root: str,
    payload: dict[str, Any],
    http_timeout: float,
    max_retries: int,
) -> dict[str, Any]:
    last_error: Exception | None = None
    for attempt in range(max_retries):
        try:
            receipt = request_json(
                session,
                "POST",
                f"{api_root}/search",
                http_timeout,
                json=payload,
            )
            if not receipt.get("id"):
                raise ValueError("FASST submission response lacks task id")
            return receipt
        except (requests.RequestException, ValueError) as exc:
            last_error = exc
            if attempt + 1 < max_retries:
                time.sleep(2**attempt)
    assert last_error is not None
    raise last_error


def poll_task(
    session: requests.Session,
    api_root: str,
    task_id: str,
    poll_interval: float,
    task_timeout: float,
    http_timeout: float,
) -> dict[str, Any]:
    deadline = time.monotonic() + task_timeout
    last_status: dict[str, Any] = {}
    transient_failures = 0
    while time.monotonic() < deadline:
        try:
            last_status = request_json(
                session,
                "GET",
                f"{api_root}/search/status/{task_id}",
                http_timeout,
            )
            transient_failures = 0
        except (requests.RequestException, ValueError):
            transient_failures += 1
            time.sleep(min(poll_interval * (2 ** min(transient_failures, 4)), 15.0))
            continue
        state = str(last_status.get("status", "")).upper()
        if state == "COMPLETED":
            return last_status
        if state == "NOT_FOUND":
            raise TaskNotFoundError(
                f"FASST task {task_id} is no longer present in the API task store"
            )
        if state in TERMINAL_FAILURE_STATES:
            raise RuntimeError(f"FASST task {task_id} ended as {state}: {last_status}")
        time.sleep(poll_interval)
    raise TimeoutError(f"FASST task {task_id} did not complete within {task_timeout}s")


def download_result_ranged(
    session: requests.Session,
    api_root: str,
    task_id: str,
    destination: Path,
    http_timeout: float,
    max_retries: int,
    chunk_bytes: int = 8 * 1024 * 1024,
) -> tuple[int, str | None]:
    """Download a result in validated byte ranges and resume a partial file."""
    url = f"{api_root}/search/result/{task_id}"
    probe = session.get(
        url,
        headers={"Range": "bytes=0-0", "Accept-Encoding": "identity"},
        timeout=http_timeout,
    )
    probe.raise_for_status()
    content_range = probe.headers.get("Content-Range", "")
    match = re.fullmatch(r"bytes 0-0/(\d+)", content_range)
    if probe.status_code != 206 or not match or len(probe.content) != 1:
        raise RuntimeError(
            f"Result endpoint lacks validated range support: status={probe.status_code}, "
            f"Content-Range={content_range!r}"
        )
    total_bytes = int(match.group(1))
    etag = probe.headers.get("ETag")
    current_bytes = destination.stat().st_size if destination.exists() else 0
    if current_bytes > total_bytes:
        # A partial created before identity-encoding was enforced is not byte-range
        # compatible with the current representation. It is temporary evidence only;
        # restart the same completed task rather than submitting a new query.
        destination.unlink()
        current_bytes = 0
    if current_bytes == 0:
        destination.write_bytes(probe.content)
        current_bytes = 1

    with destination.open("ab") as handle:
        while current_bytes < total_bytes:
            end = min(current_bytes + chunk_bytes, total_bytes) - 1
            expected = end - current_bytes + 1
            last_error: Exception | None = None
            for attempt in range(max_retries):
                try:
                    response = session.get(
                        url,
                        headers={
                            "Range": f"bytes={current_bytes}-{end}",
                            "Accept-Encoding": "identity",
                        },
                        timeout=http_timeout,
                    )
                    response.raise_for_status()
                    expected_range = f"bytes {current_bytes}-{end}/{total_bytes}"
                    if response.status_code != 206:
                        raise RuntimeError(f"Expected HTTP 206, received {response.status_code}")
                    if response.headers.get("Content-Range") != expected_range:
                        raise RuntimeError(
                            f"Unexpected Content-Range: {response.headers.get('Content-Range')!r}"
                        )
                    if len(response.content) != expected:
                        raise RuntimeError(
                            f"Incomplete range {current_bytes}-{end}: "
                            f"received {len(response.content)} of {expected} bytes"
                        )
                    handle.write(response.content)
                    handle.flush()
                    current_bytes += expected
                    last_error = None
                    break
                except (requests.RequestException, RuntimeError) as exc:
                    last_error = exc
                    if attempt + 1 < max_retries:
                        time.sleep(2**attempt)
            if last_error is not None:
                raise last_error
    return total_bytes, etag


def iter_result_rows(path: Path):
    """Yield objects from the top-level ``results`` array without loading all JSON."""
    decoder = json.JSONDecoder()
    buffer = ""
    position = 0
    eof = False
    marker = re.compile(r'"results"\s*:\s*\[')
    with path.open("r", encoding="utf-8") as handle:
        while True:
            match = marker.search(buffer, position)
            if match:
                position = match.end()
                break
            chunk = handle.read(1024 * 1024)
            if not chunk:
                raise ValueError("FASST result lacks a top-level results array")
            buffer += chunk
            if len(buffer) > 2 * 1024 * 1024:
                buffer = buffer[-1024:]
                position = 0

        while True:
            while position < len(buffer) and buffer[position] in " \t\r\n,":
                position += 1
            if position < len(buffer) and buffer[position] == "]":
                return
            try:
                row, end = decoder.raw_decode(buffer, position)
            except json.JSONDecodeError:
                if eof:
                    raise
                buffer = buffer[position:]
                position = 0
                chunk = handle.read(1024 * 1024)
                if chunk:
                    buffer += chunk
                else:
                    eof = True
                continue
            if not isinstance(row, dict):
                raise ValueError("FASST results array contains a non-object row")
            yield row
            position = end
            if position > 1024 * 1024:
                buffer = buffer[position:]
                position = 0


def compact_match_row(row: dict[str, Any]) -> dict[str, Any]:
    usi = str(row.get("USI", row.get("usi", "")))
    dataset = str(row.get("Dataset", row.get("dataset_id", row.get("dataset", ""))))
    if not dataset and ":" in usi:
        dataset = usi.split(":", 2)[1]
    return {
        "USI": usi,
        "Dataset": dataset,
        "Cosine": row.get("Cosine", row.get("cosine", row.get("score", ""))),
        "Matching Peaks": row.get("Matching Peaks", row.get("matching_peaks", "")),
        "Delta Mass": row.get("Delta Mass", row.get("delta_mass", "")),
    }


def convert_result_to_compact_tsv_gzip(
    raw_path: Path, output_path: Path
) -> tuple[int, int, float]:
    temporary = output_path.with_suffix(output_path.suffix + ".tmp")
    n_matches = 0
    datasets: set[str] = set()
    top_cosine = 0.0
    with gzip.open(temporary, "wt", encoding="utf-8", newline="\n") as handle:
        writer = csv.DictWriter(handle, fieldnames=RETAINED_MATCH_FIELDS, delimiter="\t")
        writer.writeheader()
        for row in iter_result_rows(raw_path):
            retained = compact_match_row(row)
            writer.writerow(retained)
            n_matches += 1
            dataset = str(retained["Dataset"])
            if dataset:
                datasets.add(dataset)
            value = retained["Cosine"]
            try:
                top_cosine = max(top_cosine, float(value))
            except (TypeError, ValueError):
                pass
    os.replace(temporary, output_path)
    return n_matches, len(datasets), top_cosine


def run_one(
    producer,
    spectrum: dict,
    status_path: Path,
    matches_dir: Path,
    args: argparse.Namespace,
) -> dict[str, Any]:
    feature_id = spectrum_feature_id(spectrum)
    payload = build_query_payload(
        producer,
        spectrum,
        args.database,
        args.pm_tolerance,
        args.fragment_tolerance,
        args.min_cos,
    )
    existing = read_json(status_path)
    task_id = None
    result_task = None
    submitted_utc = None
    receipt = None
    if (
        existing
        and existing.get("state") in {"pending", "error"}
        and existing.get("task_id")
        and existing.get("request_sha256") == canonical_sha256(payload)
    ):
        task_id = str(existing["task_id"])
        result_task = existing.get("result_task")
        submitted_utc = existing.get("submitted_utc")

    session = requests.Session()
    try:
        if not task_id:
            receipt = submit_query(
                session,
                args.api_root,
                payload,
                args.http_timeout,
                args.max_submit_retries,
            )
            task_id = str(receipt["id"])
            result_task = receipt.get("result_task")
            submitted_utc = utc_now()
            atomic_json(
                status_path,
                {
                    "feature_id": feature_id,
                    "state": "pending",
                    "task_id": task_id,
                    "result_task": result_task,
                    "submitted_utc": submitted_utc,
                    "request_sha256": canonical_sha256(payload),
                },
            )

        try:
            service_status = poll_task(
                session,
                args.api_root,
                task_id,
                args.poll_interval,
                args.task_timeout,
                args.http_timeout,
            )
        except TaskNotFoundError:
            # A checkpointed task may be evicted from the API task store before a
            # later retry round reaches it. Resubmit only after the API explicitly
            # reports NOT_FOUND; a merely slow PENDING task is never duplicated.
            if not existing or str(existing.get("task_id", "")) != task_id:
                raise
            receipt = submit_query(
                session,
                args.api_root,
                payload,
                args.http_timeout,
                args.max_submit_retries,
            )
            task_id = str(receipt["id"])
            result_task = receipt.get("result_task")
            submitted_utc = utc_now()
            atomic_json(
                status_path,
                {
                    "feature_id": feature_id,
                    "state": "pending",
                    "task_id": task_id,
                    "result_task": result_task,
                    "submitted_utc": submitted_utc,
                    "request_sha256": canonical_sha256(payload),
                    "resubmission_reason": "checkpointed_task_not_found",
                },
            )
            service_status = poll_task(
                session,
                args.api_root,
                task_id,
                args.poll_interval,
                args.task_timeout,
                args.http_timeout,
            )
        raw_result_path = matches_dir / f"{feature_id}_result.json.partial"
        result_bytes, result_etag = download_result_ranged(
            session,
            args.api_root,
            task_id,
            raw_result_path,
            args.http_timeout,
            args.max_submit_retries,
        )
        result_sha256 = sha256(raw_result_path)
        matches_path = matches_dir / f"{feature_id}_matches.tsv.gz"
        n_matches, n_datasets, top_cosine = convert_result_to_compact_tsv_gzip(
            raw_result_path, matches_path
        )
        raw_result_path.unlink()
        status = {
            "feature_id": feature_id,
            "state": "success",
            "task_id": task_id,
            "result_task": result_task,
            "submitted_utc": submitted_utc,
            "completed_utc": utc_now(),
            "service_status": service_status,
            "request_sha256": canonical_sha256(payload),
            "result_sha256": result_sha256,
            "result_bytes": result_bytes,
            "result_etag": result_etag,
            "n_matches": n_matches,
            "n_datasets": n_datasets,
            "top_cosine": round(top_cosine, 6),
            "matches_file": str(matches_path),
            "matches_sha256": sha256(matches_path),
            "match_storage_schema": MATCH_STORAGE_SCHEMA,
            "retained_match_fields": list(RETAINED_MATCH_FIELDS),
        }
        atomic_json(status_path, status)
        return status
    except Exception as exc:
        status = {
            "feature_id": feature_id,
            "state": "error",
            "task_id": task_id,
            "result_task": result_task,
            "submitted_utc": submitted_utc,
            "completed_utc": utc_now(),
            "request_sha256": canonical_sha256(payload),
            "error": f"{type(exc).__name__}:{str(exc)[:500]}",
        }
        atomic_json(status_path, status)
        return status
    finally:
        session.close()


def main() -> None:
    args = parse_args()
    args.api_root = args.api_root.rstrip("/")
    if args.workers < 1 or args.workers > 10:
        raise ValueError("workers must be between 1 and 10")
    if args.poll_interval <= 0 or args.task_timeout <= 0 or args.http_timeout <= 0:
        raise ValueError("timeouts and poll interval must be positive")
    if args.max_submit_retries < 1:
        raise ValueError("max-submit-retries must be positive")

    mgf_path = args.mgf.resolve()
    producer_path = args.recovered_producer.resolve()
    out_dir = args.out_dir.resolve()
    statuses_dir = out_dir / "statuses"
    matches_dir = out_dir / "per_feature"
    statuses_dir.mkdir(parents=True, exist_ok=True)
    matches_dir.mkdir(parents=True, exist_ok=True)

    producer = load_producer(producer_path)
    all_spectra = producer.parse_mgf(str(mgf_path), min_peaks=3)
    selected = select_spectra(all_spectra, args.limit, args.selection)
    selected_ids = [spectrum_feature_id(spectrum) for spectrum in selected]
    if len(set(selected_ids)) != len(selected_ids):
        raise ValueError("Selected MGF records contain duplicate feature IDs")

    tasks: list[dict] = []
    skipped_success = 0
    skipped_error = 0
    for spectrum, feature_id in zip(selected, selected_ids):
        existing = read_json(statuses_dir / f"{feature_id}.json")
        if existing and existing.get("state") == "success":
            skipped_success += 1
            continue
        if existing and existing.get("state") == "error" and not args.retry_errors:
            skipped_error += 1
            continue
        tasks.append(spectrum)

    run_manifest: dict[str, Any] = {
        "schema_version": 1,
        "taxonomy_release": TAXONOMY_RELEASE,
        "scope": "Bounded strict POS fastMASST API pilot" if args.limit else "Strict POS fastMASST run",
        "protocol": "GNPS Fast Search Library API v1 asynchronous JSON task API",
        "match_storage": {
            "schema": MATCH_STORAGE_SCHEMA,
            "format": "per-feature gzip-compressed TSV",
            "retained_fields": list(RETAINED_MATCH_FIELDS),
            "raw_response_policy": (
                "hash and byte count retained; raw response removed after compact output verifies"
            ),
        },
        "started_utc": utc_now(),
        "input_mgf": str(mgf_path),
        "input_sha256": sha256(mgf_path),
        "input_eligible_spectra": len(all_spectra),
        "selected_feature_ids": selected_ids,
        "recovered_producer": str(producer_path),
        "recovered_producer_sha256": sha256(producer_path),
        "wrapper": str(Path(__file__).resolve()),
        "wrapper_sha256": sha256(Path(__file__).resolve()),
        "api_root": args.api_root,
        "database": args.database,
        "parameters": {
            "workers": args.workers,
            "pm_tolerance": args.pm_tolerance,
            "fragment_tolerance": args.fragment_tolerance,
            "cosine_threshold": args.min_cos,
            "poll_interval_seconds": args.poll_interval,
            "task_timeout_seconds": args.task_timeout,
            "http_timeout_seconds": args.http_timeout,
            "max_submit_retries": args.max_submit_retries,
            "retry_errors": args.retry_errors,
            "limit": args.limit,
            "selection": args.selection,
        },
        "skipped_existing_success": skipped_success,
        "skipped_existing_error": skipped_error,
        "scheduled_queries": len(tasks),
        "complete": False,
    }
    atomic_json(out_dir / "run_manifest.json", run_manifest)

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {}
        for spectrum in tasks:
            feature_id = spectrum_feature_id(spectrum)
            future = executor.submit(
                run_one,
                producer,
                spectrum,
                statuses_dir / f"{feature_id}.json",
                matches_dir,
                args,
            )
            futures[future] = feature_id
        for future in tqdm(as_completed(futures), total=len(futures), desc="fastMASST async"):
            future.result()

    statuses = []
    for feature_id in selected_ids:
        status = read_json(statuses_dir / f"{feature_id}.json")
        if status:
            statuses.append(status)
    summary_path = out_dir / "masst_summary.csv"
    temporary_summary = summary_path.with_suffix(".csv.tmp")
    pd.DataFrame(statuses).sort_values("feature_id").to_csv(temporary_summary, index=False)
    os.replace(temporary_summary, summary_path)

    successes = sum(status.get("state") == "success" for status in statuses)
    errors = sum(status.get("state") == "error" for status in statuses)
    pending = sum(status.get("state") == "pending" for status in statuses)
    run_manifest.update(
        {
            "finished_utc": utc_now(),
            "status_records": len(statuses),
            "success_records": successes,
            "error_records": errors,
            "pending_records": pending,
            "features_with_matches": sum(
                status.get("state") == "success" and int(status.get("n_matches", 0)) > 0
                for status in statuses
            ),
            "total_matches": sum(
                int(status.get("n_matches", 0))
                for status in statuses
                if status.get("state") == "success"
            ),
            "summary_sha256": sha256(summary_path),
            "complete": len(statuses) == len(selected_ids) and errors == 0 and pending == 0,
        }
    )
    atomic_json(out_dir / "run_manifest.json", run_manifest)
    print(json.dumps(run_manifest, indent=2))
    if not run_manifest["complete"]:
        sys.exit(2)


if __name__ == "__main__":
    main()
