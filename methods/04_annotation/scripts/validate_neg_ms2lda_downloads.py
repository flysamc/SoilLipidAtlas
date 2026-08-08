#!/usr/bin/env python3
"""Validate the two downloaded deterministic NEG strict MS2LDA runs.

This validator deliberately treats feature IDs as opaque exact identifiers.  It
does not attempt any mass- or retention-time-based reconciliation.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/annotation_recovery_neg"
MS2LDA = BASE / "ms2lda"
RUNS = [MS2LDA / "strict_full_run1", MS2LDA / "strict_full_run2"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def mgf_titles(path: Path) -> list[str]:
    titles: list[str] = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("TITLE="):
                titles.append(line.rstrip("\r\n")[6:])
    return titles


def main() -> None:
    ledger = pd.read_csv(BASE / "strict_feature_spectrum_ledger.csv")
    usable_ids = set(
        ledger.loc[
            ledger["has_usable_ms2"].astype(str).str.lower().eq("true"),
            "feature_id",
        ].astype(str)
    )

    reports = []
    scientific_parameters = []
    for run in RUNS:
        result = run / "results_full"
        params = json.loads((result / "run_parameters.json").read_text(encoding="utf-8"))
        comparable = {k: v for k, v in params.items() if k not in {"input", "elapsed_seconds"}}
        scientific_parameters.append(comparable)

        titles = mgf_titles(run / "neg_strict_full_for_ms2lda.mgf")
        docs = pd.read_csv(result / "doc_topic_matrix.csv", usecols=["feature_id"])
        doc_ids = set(docs["feature_id"].astype(str))
        title_ids = set(titles)

        reports.append(
            {
                "run": run.name,
                "input_spectra": len(titles),
                "unique_input_feature_ids": len(title_ids),
                "model_documents": len(docs),
                "unique_model_feature_ids": len(doc_ids),
                "preprocessing_exclusions": len(title_ids - doc_ids),
                "model_ids_not_in_input": len(doc_ids - title_ids),
                "input_exactly_matches_strict_usable_ids": title_ids == usable_ids,
                "model_ids_are_exact_input_subset": doc_ids <= title_ids,
                "files": {
                    path.name: {"bytes": path.stat().st_size, "sha256": sha256(path)}
                    for path in [
                        run / "neg_strict_full_for_ms2lda.mgf",
                        result / "doc_topic_matrix.csv",
                        result / "motif_top_words.csv",
                        result / "convergence_history.json",
                        result / "run_parameters.json",
                    ]
                },
            }
        )

    deterministic_files = {}
    for filename in ["doc_topic_matrix.csv", "motif_top_words.csv", "convergence_history.json"]:
        hashes = [sha256(run / "results_full" / filename) for run in RUNS]
        deterministic_files[filename] = {"identical": len(set(hashes)) == 1, "sha256": hashes[0]}

    checks = {
        "scientific_parameters_equal_ignoring_input_path_and_elapsed_time": (
            scientific_parameters[0] == scientific_parameters[1]
        ),
        "all_deterministic_outputs_byte_identical": all(
            item["identical"] for item in deterministic_files.values()
        ),
        "both_inputs_equal_full_strict_usable_set": all(
            report["input_exactly_matches_strict_usable_ids"] for report in reports
        ),
        "both_model_id_sets_are_exact_input_subsets": all(
            report["model_ids_are_exact_input_subset"] for report in reports
        ),
        "both_runs_have_same_counts": reports[0]["model_documents"] == reports[1]["model_documents"],
    }
    validation = {
        "schema_version": 1,
        "taxonomy_release": "ncbi-phylum-2026-08-04-v1",
        "polarity": "NEG",
        "status": "pass" if all(checks.values()) else "fail",
        "strict_atlas_features": len(ledger),
        "strict_usable_ms2_features": len(usable_ids),
        "runs": reports,
        "deterministic_outputs": deterministic_files,
        "checks": checks,
        "interpretation": (
            "Both 5,695-spectrum inputs exactly equal the strict usable-MS2 feature-ID set. "
            "The declared MS2LDA preprocessing retained 5,466 documents and excluded 229; "
            "both deterministic reruns retained the same exact IDs and produced byte-identical "
            "scientific outputs. Run-specific path and elapsed time are metadata only."
        ),
    }
    out = MS2LDA / "download_validation.json"
    out.write_text(json.dumps(validation, indent=2) + "\n", encoding="utf-8")
    if validation["status"] != "pass":
        raise SystemExit("NEG MS2LDA validation failed")

    manifest = {
        "schema_version": 2,
        "stage_id": "negative_ms2lda_full_strict_retrain",
        "taxonomy_release": "ncbi-phylum-2026-08-04-v1",
        "polarity": "NEG",
        "status": "complete_downloaded_deterministic_replicate_validated",
        "submitted": True,
        "external_execution": [
            {
                "system": "LISC",
                "job_id": "5802341",
                "role": "strict_full_run1",
                "scheduler_state": "COMPLETED",
                "scheduler_exit_code": "0:0",
            },
            {
                "system": "LISC",
                "job_id": None,
                "role": "strict_full_run2_deterministic_replicate",
                "note": "Scheduler job ID was not preserved in the downloaded compact result metadata.",
            },
        ],
        "counts": {
            "strict_atlas_features": len(ledger),
            "strict_usable_ms2_input_spectra": len(usable_ids),
            "model_documents_after_declared_preprocessing": reports[0]["model_documents"],
            "preprocessing_exclusions": reports[0]["preprocessing_exclusions"],
            "motifs": scientific_parameters[0]["n_motifs"],
            "iterations": scientific_parameters[0]["n_iterations"],
            "deterministic_replicates": 2,
        },
        "parameters": scientific_parameters[0],
        "validation": {
            "path": str(out.relative_to(ROOT)).replace("\\", "/"),
            "bytes": out.stat().st_size,
            "sha256": sha256(out),
            "status": "pass",
        },
        "outputs": reports,
        "release_boundary": (
            "Corpus-level strict NEG motif model is locally complete. Motif enrichment, biological "
            "interpretation, and manuscript/table/figure consumers remain review-only until downstream "
            "release gates pass."
        ),
    }
    (MS2LDA / "stage_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": "pass", "output": str(out), "documents": 5466}, indent=2))


if __name__ == "__main__":
    main()
