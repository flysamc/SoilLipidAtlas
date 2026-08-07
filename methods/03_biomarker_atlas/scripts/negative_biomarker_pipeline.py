#!/usr/bin/env python3
"""Audit and gated strict negative-mode biomarker recovery.

The command first audits the recovered NEG substrate, historical Composite and
IndVal producers, exact cached annotations, and consensus-MGF usability.  It
only runs the locked strict selector after the historical Composite gate
passes.  A failed gate writes a complete evidence package and exits non-zero;
it never fabricates a strict atlas or submits external annotation jobs.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from paper2_repro.strict_biomarker_core import (
    audit_mgf_reference_ids,
    exact_feature_ids,
    fingerprint,
    load_locked_taxonomy,
    normalize_mode_metadata,
    pair_report,
    sha256,
)
from working_corrections.figure_2a import rerun_fig2a as engine


RELEASE_ID = "ncbi-phylum-2026-08-04-v1"
NEG_ALIGNMENT_BATCH_ORDER = (
    "OE23-NEG",
    "OE21-4-NEG",
    "ALL-25-2-NEG",
    "OE26-1-NEG",
    "OE11-3-NEG",
    "OE25-1-NEG",
)


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")


def configure_engine(consensus: Path) -> None:
    # Reuse the recovered POS implementation without changing its POS globals
    # on disk or duplicating selection thresholds in the NEG adapter.
    engine.CONSENSUS = consensus


def read_metadata(path: Path) -> pd.DataFrame:
    return normalize_mode_metadata(pd.read_csv(path, low_memory=False), sample_column="sample_col")


def read_recovery_metadata(path: Path) -> pd.DataFrame:
    """Adapt the recovered historical producer's sample map to the shared engine contract."""
    frame = pd.read_csv(path, low_memory=False)
    frame = frame[frame["include_in_analysis"].eq("Yes")].copy()
    frame["sample_col"] = "sample:" + frame["sample_name"].astype(str).str.strip()
    return normalize_mode_metadata(frame, sample_column="sample_col")


def pair_counts(expected: pd.DataFrame, observed: pd.DataFrame, phyla: list[str]) -> pd.DataFrame:
    rows = []
    for phylum in phyla:
        left = expected.loc[expected["phylum"].eq(phylum), "feature_id"].nunique()
        right = observed.loc[observed["phylum"].eq(phylum), "feature_id"].nunique()
        overlap = len(
            set(expected.loc[expected["phylum"].eq(phylum), "feature_id"])
            & set(observed.loc[observed["phylum"].eq(phylum), "feature_id"])
        )
        rows.append(
            {
                "phylum": phylum,
                "expected": int(left),
                "observed": int(right),
                "overlap": int(overlap),
                "retained_from_historical": int(overlap),
                "new_relative_to_historical": int(max(right - overlap, 0)),
            }
        )
    return pd.DataFrame(rows)


def source_paths(source_root: Path, recovery_root: Path | None = None) -> dict[str, Path]:
    external = source_root / "external/SOILMASS_PRODUCER_RECOVERY_2026-08-04/payload/core/workspace"
    paths = {
        "consensus": source_root / "analysis/analysis-16/negative_mode/03_alignment/consensus_aligned_table.csv",
        "historical_metadata": source_root / "analysis/analysis-16/negative_mode/04_biomarker_discovery/01_composite_scoring/sample_metadata.csv",
        "historical_atlas": external / "analysis/analysis-16/negative_mode/04_biomarker_discovery/06_unified_annotations/atlas_unified_annotations.csv",
        "historical_indval": external / "analysis/analysis-16/negative_mode/04_biomarker_discovery/12_indval/multi_batch_indicators.csv",
        "historical_indval_script": external / "analysis/analysis-16/negative_mode/scripts/indval_within_batch_neg.py",
        "external_inputs": source_root / "config/external_inputs.csv",
        "taxonomy_policy": source_root / "paper2_repro/config/taxonomy_policy.json",
        "diagnostic_ms2": external / "analysis/analysis-16/negative_mode/04_biomarker_discovery/05_ms2_structural/atlas_ms2_diagnostic_classification.csv",
        "fastmasst": external / "analysis/analysis-16/negative_mode/04_biomarker_discovery/17_masst_validation/masst_results/masst_summary.csv",
        "dreams": external / "analysis/analysis-16/negative_mode/dreams_results/neg_atlas_dreams_search_results.tsv",
        "ms2lda": external / "analysis/analysis-16/negative_mode/ms2lda_results/results_neg/doc_topic_matrix.csv",
        "ms2lda_annotated": external / "analysis/analysis-16/negative_mode/ms2lda_results/results_neg_annotated/doc_topic_matrix.csv",
        "sirius": external / "Dreams/results/sirius_atlas_neg/formula_identifications.tsv",
        "canopus": external / "Dreams/results/sirius_atlas_neg/canopus_formula_summary.tsv",
        "archlips": external / "analysis/analysis-16/negative_mode/04_biomarker_discovery/07_archlips/archlips_atlas_matches.csv",
        "mgf_dir": external / "analysis/FBMN_all_batches_NEG",
    }
    if recovery_root is not None:
        paths.update(
            {
                "recovery_root": recovery_root,
                "recovery_archive": recovery_root.parent.parent / "NEG_REC.zip",
                "consensus": recovery_root / "negative_mode/03_alignment/consensus_aligned_table.csv",
                "historical_metadata": recovery_root / "negative_mode/00_sample_mapping/neg_sample_metadata.csv",
                "recovery_expected_atlas": recovery_root / "reference/expected_atlas_combined_biomarkers.csv",
                "recovery_atlas": recovery_root / "negative_mode/04_biomarker_discovery/02_platinum_diamond/atlas_combined_biomarkers.csv",
                "recovery_composite_script": recovery_root / "negative_mode/04_biomarker_discovery/01_composite_scoring/composite_scoring_neg.py",
                "recovery_selector_script": recovery_root / "negative_mode/04_biomarker_discovery/02_platinum_diamond/platinum_diamond_neg.py",
            }
        )
    return paths


def semantic_atlas_report(expected: pd.DataFrame, observed: pd.DataFrame) -> dict[str, object]:
    """Compare historical outputs by identity and values, tolerating platform serialization."""
    expected = expected.copy()
    observed = observed.copy()
    expected_ids = set(expected["feature_id"].astype(str))
    observed_ids = set(observed["feature_id"].astype(str))
    merged = expected.merge(
        observed,
        on="feature_id",
        how="outer",
        suffixes=("_expected", "_observed"),
        indicator=True,
    )
    categorical_mismatches: dict[str, int] = {}
    numeric_max_abs_diff: dict[str, float] = {}
    for column in ["phylum", "kingdom", "biomarker_tier", "tier"]:
        left = f"{column}_expected"
        right = f"{column}_observed"
        if left not in merged or right not in merged:
            continue
        mismatches = int(
            (
                merged[left].fillna("<NA>").astype(str)
                != merged[right].fillna("<NA>").astype(str)
            ).sum()
        )
        categorical_mismatches[column] = mismatches
    for column in [
        "consensus_mz",
        "consensus_rt",
        "composite_score",
        "detection_rate",
        "specificity",
        "exclusivity",
        "log2fc",
        "mean_target",
        "mean_background",
        "batch_rate",
    ]:
        left = f"{column}_expected"
        right = f"{column}_observed"
        if left not in merged or right not in merged:
            continue
        values = pd.to_numeric(merged[left], errors="coerce") - pd.to_numeric(
            merged[right], errors="coerce"
        )
        numeric_max_abs_diff[column] = float(values.abs().max(skipna=True))
    max_numeric_diff = max(numeric_max_abs_diff.values(), default=0.0)
    return {
        "expected_rows": int(len(expected)),
        "observed_rows": int(len(observed)),
        "expected_unique_feature_ids": int(len(expected_ids)),
        "observed_unique_feature_ids": int(len(observed_ids)),
        "feature_id_overlap": int(len(expected_ids & observed_ids)),
        "rerun_only_ids": int(len(observed_ids - expected_ids)),
        "historical_only_ids": int(len(expected_ids - observed_ids)),
        "categorical_mismatches": categorical_mismatches,
        "numeric_max_abs_diff": numeric_max_abs_diff,
        "numeric_tolerance": 1e-9,
        "byte_exact": False,
        "exact": (
            expected_ids == observed_ids
            and all(value == 0 for value in categorical_mismatches.values())
            and max_numeric_diff <= 1e-9
        ),
        "note": "Semantic identity/value comparison; CSV byte hashes can differ across macOS and Windows line endings and stable tie ordering.",
    }


def historical_gate(
    paths: dict[str, Path],
    taxonomy: dict,
    historical_meta: pd.DataFrame,
    features: pd.DataFrame,
    samples: list[str],
    matrix,
    output_dir: Path,
    recovery_root: Path | None = None,
) -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    atlas = pd.read_csv(paths["historical_atlas"], usecols=["feature_id", "phylum", "composite_score", "biomarker_tier", "tier"], low_memory=False)
    indval_expected = pd.read_csv(paths["historical_indval"], usecols=["feature_id", "phylum"], low_memory=False)
    recovery_report = None
    if recovery_root is not None and paths.get("recovery_atlas", Path()).is_file():
        observed_composite = pd.read_csv(paths["recovery_atlas"], low_memory=False)
        observed_platinum = observed_composite[
            observed_composite["biomarker_tier"].ne("Silver")
        ].copy()
        observed_silver = observed_composite[
            observed_composite["biomarker_tier"].eq("Silver")
        ].copy()
        recovery_expected = pd.read_csv(paths["recovery_expected_atlas"], low_memory=False)
        recovery_report = semantic_atlas_report(recovery_expected, observed_composite)
        selection_engine = "recovered historical NEG Composite producer outputs"
    else:
        observed_platinum, observed_silver, observed_composite = engine.composite_select(
            features, samples, matrix, historical_meta, "historical NEG labels"
        )
        selection_engine = "working_corrections.figure_2a.rerun_fig2a diagnostic fallback"
    observed_pairs, observed_indval = engine.run_indval(
        features, samples, matrix, historical_meta, "historical NEG labels"
    )
    expected_pairs = atlas[["feature_id", "phylum"]].drop_duplicates()
    observed_composite_pairs = observed_composite[["feature_id", "phylum"]]
    composite_report = pair_report(expected_pairs, observed_composite_pairs)
    indval_report = pair_report(indval_expected, observed_pairs[["feature_id", "phylum"]])

    composite_dir = paths["historical_metadata"].parent
    if recovery_root is not None:
        producer_files = [
            path
            for path in [
                paths.get("recovery_composite_script", Path()),
                paths.get("recovery_selector_script", Path()),
                paths.get("recovery_atlas", Path()),
            ]
            if path.is_file()
        ]
    else:
        producer_files = [
            path for path in composite_dir.rglob("*")
            if path.is_file() and path.name != paths["historical_metadata"].name
        ]
    producer_missing = not producer_files
    source_control = {
        "composite_scoring_directory": str(composite_dir),
        "producer_files_found": [str(path) for path in producer_files],
        "producer_missing": producer_missing,
        "interpretation": (
            "Recovered historical producer and selected atlas are present and used for the gate."
            if recovery_root is not None and not producer_missing
            else "Only sample_metadata.csv is present in the recovered NEG Composite directory; the original score/selection producer and candidate table are absent."
            if producer_missing
            else "Recovered Composite producer artifacts are present and require source-level verification."
        ),
    }
    gate = {
        "release_id": RELEASE_ID,
        "polarity": "NEG",
        "status": "pass"
        if composite_report["exact"]
        and indval_report["exact"]
        and not producer_missing
        and (recovery_report is None or recovery_report["exact"])
        else "failed",
        "historical_atlas_rows": int(len(atlas)),
        "historical_atlas_distinct_features": int(atlas["feature_id"].nunique()),
        "historical_atlas_nonviral_rows": int((atlas["phylum"] != "Virus").sum()),
        "historical_atlas_unique_phyla": sorted(atlas["phylum"].dropna().unique()),
        "historical_composite": composite_report,
        "historical_recovery_output": recovery_report,
        "historical_recovery_reference_hash": (
            sha256(paths["recovery_expected_atlas"])
            if paths.get("recovery_expected_atlas", Path()).is_file()
            else None
        ),
        "historical_recovery_observed_hash": (
            sha256(paths["recovery_atlas"])
            if paths.get("recovery_atlas", Path()).is_file()
            else None
        ),
        "historical_indval": indval_report,
        "historical_selector_diagnostic_only": {
            "composite_observed_rows": int(len(observed_composite)),
            "composite_retained_from_historical": int(composite_report["overlap"]),
            "composite_new_relative_to_historical": int(
                composite_report["observed"] - composite_report["overlap"]
            ),
            "note": "Historical gate evidence; strict counts are generated only from locked 16-phylum metadata after this gate passes.",
        },
        "indval_reproduction_exact": bool(indval_report["exact"]),
        "composite_reproduction_exact": bool(composite_report["exact"]),
        "producer_control": source_control,
        "selection_engine": selection_engine,
        "gate_rule": "Composite and IndVal exact pair reproduction plus recovered Composite producer evidence",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    atlas_pairs = atlas[["feature_id", "phylum"]]
    observed_composite_pairs = observed_composite[["feature_id", "phylum"]]
    observed_indval_pairs = observed_pairs[["feature_id", "phylum"]]
    diagnostic_phyla = sorted(set(atlas["phylum"].dropna()) | set(observed_composite["phylum"].dropna()))
    pair_counts(atlas_pairs, observed_composite_pairs, diagnostic_phyla).to_csv(
        output_dir / "historical_composite_counts_by_analysis_phylum.csv", index=False
    )
    indval_phyla = sorted(set(indval_expected["phylum"].dropna()) | set(observed_indval["phylum"].dropna()))
    pair_counts(indval_expected, observed_indval_pairs, indval_phyla).to_csv(
        output_dir / "historical_indval_counts_by_analysis_phylum.csv", index=False
    )
    observed_composite.to_csv(output_dir / "historical_composite_reproduced_diagnostic.csv", index=False)
    observed_platinum.to_csv(output_dir / "historical_composite_platinum_reproduced_diagnostic.csv", index=False)
    observed_silver.to_csv(output_dir / "historical_composite_silver_reproduced_diagnostic.csv", index=False)
    observed_pairs.to_csv(output_dir / "historical_indval_pairs_reproduced.csv", index=False)
    observed_indval.to_csv(output_dir / "historical_indval_unique_reproduced.csv", index=False)
    write_json(output_dir / "historical_neg_reproduction_gate.json", gate)
    return gate, atlas, observed_composite


def read_ids(path: Path, *, sep: str = ",", column: str = "feature_id", valid: set[str]) -> dict:
    if not path.is_file():
        return {"status": "missing", "source": fingerprint(path), "source_rows": 0, "raw_ids": 0, "normalized_known_ids": 0, "duplicate_normalized_ids": 0, "exact_ids": set(), "note": "No recovered cache file at the configured path."}
    frame = pd.read_csv(path, sep=sep, usecols=[column], low_memory=False)
    exact, raw_count, duplicates = exact_feature_ids(frame[column].tolist(), valid)
    return {"status": "cached", "source": fingerprint(path), "source_rows": int(len(frame)), "raw_ids": int(raw_count), "normalized_known_ids": int(len(exact)), "duplicate_normalized_ids": int(duplicates), "exact_ids": exact, "note": "Exact feature-ID reuse only; numeric IDs were normalized only against the known NEG feature whitelist."}


def _batch_token(value: object) -> str:
    return "".join(character for character in str(value).casefold() if character.isalnum())


def _reference_batch(batches: object) -> str:
    """Return the alignment-order batch that owns ``ref_batch_id``."""
    available = {value.strip() for value in str(batches).split(",") if value.strip()}
    for batch in NEG_ALIGNMENT_BATCH_ORDER:
        if batch in available:
            return batch
    raise ValueError(f"No known NEG alignment batch in {batches!r}")


def audited_mgf_feature_sets(paths: dict[str, Path], atlas: pd.DataFrame) -> dict[str, object]:
    """Map MGF spectra to features only through one-to-one batch-local references."""
    consensus = pd.read_csv(
        paths["consensus"], usecols=["feature_id", "ref_batch_id", "batches"], low_memory=False
    )
    target = atlas[["feature_id"]].merge(consensus, on="feature_id", how="left")
    pair_to_features: dict[tuple[str, str], set[str]] = {}
    for row in target.itertuples(index=False):
        if pd.isna(row.ref_batch_id) or pd.isna(row.batches):
            continue
        ref_id = str(int(row.ref_batch_id))
        # ``ref_batch_id`` is local to the first batch in the recovered
        # alignment order.  Testing that scan against every listed batch can
        # create false-positive spectra with unrelated local scan IDs.
        batch = _reference_batch(row.batches)
        pair_to_features.setdefault((_batch_token(batch), ref_id), set()).add(str(row.feature_id))

    known_batches = sorted({pair[0] for pair in pair_to_features})
    paths_by_batch: dict[str, list[Path]] = {batch: [] for batch in known_batches}
    unmatched_paths: list[str] = []
    for path in sorted(paths["mgf_dir"].rglob("*_consensus_spectra.mgf")):
        path_token = _batch_token(path.parent.parent.name)
        matches = [batch for batch in known_batches if batch in path_token]
        if len(matches) == 1:
            paths_by_batch[matches[0]].append(path)
        else:
            unmatched_paths.append(str(path))

    found_by_feature: dict[str, int] = {}
    peak_by_feature: dict[str, int] = {}
    file_reports: list[dict[str, object]] = []
    ambiguous_pairs = {
        pair for pair, feature_ids in pair_to_features.items() if len(feature_ids) != 1
    }
    for batch, mgf_paths in paths_by_batch.items():
        valid_refs = {
            "NEG_" + ref_id
            for pair, feature_ids in pair_to_features.items()
            if pair[0] == batch and len(feature_ids) == 1
            for ref_id in [pair[1]]
        }
        if not valid_refs:
            continue
        report = audit_mgf_reference_ids(mgf_paths, valid_refs)
        file_reports.extend(report["files"])
        for ref_key in report["target_feature_ids_found"]:
            ref_id = ref_key.removeprefix("NEG_")
            feature_ids = pair_to_features.get((batch, ref_id), set())
            if len(feature_ids) == 1:
                feature_id = next(iter(feature_ids))
                found_by_feature[feature_id] = found_by_feature.get(feature_id, 0) + 1
        for ref_key in report["target_peak_bearing_feature_ids"]:
            ref_id = ref_key.removeprefix("NEG_")
            feature_ids = pair_to_features.get((batch, ref_id), set())
            if len(feature_ids) == 1:
                feature_id = next(iter(feature_ids))
                peak_by_feature[feature_id] = peak_by_feature.get(feature_id, 0) + 1
    found = set(found_by_feature)
    peak = set(peak_by_feature)
    return {
        "files": file_reports,
        "n_files": len(file_reports),
        "feature_ids_found": len(found),
        "peak_bearing_feature_ids": len(peak),
        "feature_ids_without_peak_bearing_spectrum": len(found - peak),
        "feature_ids_in_multiple_files": sum(count > 1 for count in found_by_feature.values()),
        "ambiguous_reference_pairs_excluded": len(ambiguous_pairs),
        "features_with_ambiguous_reference_excluded": len(
            set().union(*(pair_to_features[pair] for pair in ambiguous_pairs))
            if ambiguous_pairs
            else set()
        ),
        "unmatched_mgf_files": unmatched_paths,
        "target_feature_ids_found": found,
        "target_peak_bearing_feature_ids": peak,
        "provenance_note": "MGF reuse requires an exact consensus reference plus one-to-one batch-local provenance; ambiguous ref_batch_id mappings are excluded.",
    }


def cache_coverage(paths: dict[str, Path], atlas: pd.DataFrame, output_dir: Path) -> dict:
    all_ids = set(atlas["feature_id"].astype(str))
    nonviral_ids = set(atlas.loc[atlas["phylum"].ne("Virus"), "feature_id"].astype(str))
    stages: dict[str, dict] = {}
    stages["diagnostic_annotation"] = read_ids(paths["diagnostic_ms2"], valid=all_ids)
    stages["fastMASST"] = read_ids(paths["fastmasst"], valid=all_ids)
    stages["DreaMS"] = read_ids(paths["dreams"], sep="\t", valid=all_ids)
    stages["MS2LDA"] = read_ids(paths["ms2lda"], valid=all_ids)
    stages["MS2LDA_annotated"] = read_ids(paths["ms2lda_annotated"], valid=all_ids)
    sirius = read_ids(paths["sirius"], sep="\t", column="mappingFeatureId", valid=all_ids)
    canopus = read_ids(paths["canopus"], sep="\t", column="mappingFeatureId", valid=all_ids)
    stages["SIRIUS"] = sirius
    stages["CANOPUS"] = canopus
    sirius_canopus = {
        "status": "cached" if sirius["status"] == "cached" or canopus["status"] == "cached" else "missing",
        "source": {"formula_identifications": sirius["source"], "canopus_formula_summary": canopus["source"]},
        "source_rows": int(sirius["source_rows"] + canopus["source_rows"]),
        "raw_ids": int(sirius["raw_ids"] + canopus["raw_ids"]),
        "normalized_known_ids": len(sirius["exact_ids"] | canopus["exact_ids"]),
        "duplicate_normalized_ids": int(sirius["duplicate_normalized_ids"] + canopus["duplicate_normalized_ids"]),
        "exact_ids": sirius["exact_ids"] | canopus["exact_ids"],
        "note": "Union of exact mappingFeatureId values from SIRIUS formula and CANOPUS caches.",
    }
    stages["SIRIUS_CANOPUS"] = sirius_canopus
    stages["Pan-ReDU"] = {
        "status": "not_located",
        "source": {"configured_status": "partly_cached", "path": "external/panredu_full_hits/"},
        "source_rows": 0,
        "raw_ids": 0,
        "normalized_known_ids": 0,
        "duplicate_normalized_ids": 0,
        "exact_ids": set(),
        "note": "Configured as partly cached, but no complete per-feature NEG hit table was located; no reuse claimed."
    }

    mgf = audited_mgf_feature_sets(paths, atlas)
    found_features = mgf["target_feature_ids_found"]
    peak_features = mgf["target_peak_bearing_feature_ids"]
    rows = []
    for stage, report in stages.items():
        ids = report.pop("exact_ids")
        report["exact_all_atlas"] = len(ids)
        report["exact_nonviral_atlas"] = len(ids & nonviral_ids)
        report["coverage_nonviral_pct"] = round(100 * report["exact_nonviral_atlas"] / max(len(nonviral_ids), 1), 4)
        rows.append({"stage": stage, **{key: value for key, value in report.items() if not isinstance(value, dict)}})
    coverage = {
        "release_id": RELEASE_ID,
        "polarity": "NEG",
        "scope": "historical unified NEG atlas; nonviral rows are the recoverable core comparison",
        "atlas_rows": int(len(atlas)),
        "nonviral_atlas_rows": int(len(nonviral_ids)),
        "stages": stages,
        "mgf": {key: value for key, value in mgf.items() if key not in {"target_feature_ids_found", "target_peak_bearing_feature_ids"}},
        "mgf_reference_ids_expected": None,
        "mgf_reference_ids_found": None,
        "mgf_peak_bearing_reference_ids": None,
        "mgf_exact_atlas_feature_ids_found": len(found_features & all_ids),
        "mgf_exact_nonviral_feature_ids_found": len(found_features & nonviral_ids),
        "mgf_peak_bearing_atlas_feature_ids": len(peak_features & all_ids),
        "mgf_peak_bearing_nonviral_feature_ids": len(peak_features & nonviral_ids),
        "id_reuse_rule": "Exact feature ID plus source/provenance compatibility; no mass-only reuse.",
    }
    pd.DataFrame(rows).to_csv(output_dir / "historical_neg_annotation_cache_coverage.csv", index=False)
    write_json(output_dir / "historical_neg_annotation_cache_coverage.json", coverage)
    return coverage


def strict_annotation_queue(
    paths: dict[str, Path],
    strict_atlas: pd.DataFrame,
    historical_atlas: pd.DataFrame,
    output_dir: Path,
) -> dict:
    """Partition strict features into exact cache reuse and genuinely new queues."""
    strict_ids = set(strict_atlas["feature_id"].astype(str))
    historical_nonviral = historical_atlas[historical_atlas["phylum"].ne("Virus")]
    historical_ids = set(historical_nonviral["feature_id"].astype(str))
    new_ids = strict_ids - historical_ids
    retained_ids = strict_ids & historical_ids
    old_labels = historical_nonviral.groupby("feature_id")["phylum"].apply(set).to_dict()
    new_labels = strict_atlas.groupby("feature_id")["phylum"].apply(set).to_dict()
    reassigned_ids = {
        feature_id
        for feature_id in retained_ids
        if old_labels.get(feature_id, set()) != new_labels.get(feature_id, set())
    }

    mgf = audited_mgf_feature_sets(paths, strict_atlas)
    mgf_found_ids = mgf["target_feature_ids_found"]
    usable_ms2_ids = mgf["target_peak_bearing_feature_ids"]
    new_with_ms2 = new_ids & usable_ms2_ids
    new_without_ms2 = new_ids - usable_ms2_ids

    stage_paths = {
        "diagnostic_annotation": (paths["diagnostic_ms2"], ",", "feature_id"),
        "SIRIUS": (paths["sirius"], "\t", "mappingFeatureId"),
        "CANOPUS": (paths["canopus"], "\t", "mappingFeatureId"),
        "DreaMS": (paths["dreams"], "\t", "feature_id"),
        "MS2LDA": (paths["ms2lda"], ",", "feature_id"),
        "MS2LDA_annotated": (paths["ms2lda_annotated"], ",", "feature_id"),
        "fastMASST": (paths["fastmasst"], ",", "feature_id"),
    }
    stage_reports = {}
    stage_exact_sets: dict[str, set[str]] = {}
    queue_rows = []
    for stage, (path, sep, column) in stage_paths.items():
        report = read_ids(path, sep=sep, column=column, valid=strict_ids)
        exact_ids = report.pop("exact_ids")
        stage_exact_sets[stage] = exact_ids
        cached_new_ids = exact_ids & new_ids
        external_ids = new_ids - exact_ids
        stage_reports[stage] = {
            **report,
            "exact_strict_ids": int(len(exact_ids)),
            "exact_retained_ids": int(len(exact_ids & retained_ids)),
            "exact_reassigned_ids": int(len(exact_ids & reassigned_ids)),
            "exact_new_ids": int(len(cached_new_ids)),
            "external_required_new_ids": int(len(external_ids)),
            "external_required_new_with_usable_ms2": int(len(external_ids & usable_ms2_ids)),
            "external_required_new_without_usable_ms2": int(len(external_ids - usable_ms2_ids)),
        }
        queue_rows.append(
            {
                "stage": stage,
                "strict_ids": len(strict_ids),
                "exact_cache_reuse": len(exact_ids),
                "retained_ids": len(retained_ids),
                "reassigned_ids": len(reassigned_ids),
                "new_ids": len(new_ids),
                "new_ids_with_cached_evidence": len(cached_new_ids),
                "external_required_new_ids": len(external_ids),
                "external_required_new_with_usable_ms2": len(external_ids & usable_ms2_ids),
                "external_required_new_without_usable_ms2": len(external_ids - usable_ms2_ids),
            }
        )

    sirius_canopus_ids = stage_exact_sets["SIRIUS"] | stage_exact_sets["CANOPUS"]
    sirius_canopus_new = new_ids - sirius_canopus_ids
    stage_reports["SIRIUS_CANOPUS"] = {
        "status": "cached",
        "exact_strict_ids": len(sirius_canopus_ids),
        "exact_retained_ids": len(sirius_canopus_ids & retained_ids),
        "exact_reassigned_ids": len(sirius_canopus_ids & reassigned_ids),
        "exact_new_ids": len(sirius_canopus_ids & new_ids),
        "external_required_new_ids": len(sirius_canopus_new),
        "external_required_new_with_usable_ms2": len(sirius_canopus_new & usable_ms2_ids),
        "external_required_new_without_usable_ms2": len(sirius_canopus_new - usable_ms2_ids),
        "note": "Union of exact feature IDs from SIRIUS formula and CANOPUS caches.",
    }
    queue_rows.append(
        {
            "stage": "SIRIUS_CANOPUS",
            "strict_ids": len(strict_ids),
            "exact_cache_reuse": len(sirius_canopus_ids),
            "retained_ids": len(retained_ids),
            "reassigned_ids": len(reassigned_ids),
            "new_ids": len(new_ids),
            "new_ids_with_cached_evidence": len(sirius_canopus_ids & new_ids),
            "external_required_new_ids": len(sirius_canopus_new),
            "external_required_new_with_usable_ms2": len(sirius_canopus_new & usable_ms2_ids),
            "external_required_new_without_usable_ms2": len(sirius_canopus_new - usable_ms2_ids),
        }
    )

    panredu_external = new_ids
    stage_reports["Pan-ReDU"] = {
        "status": "not_located",
        "exact_strict_ids": 0,
        "exact_retained_ids": 0,
        "exact_reassigned_ids": 0,
        "exact_new_ids": 0,
        "external_required_new_ids": len(panredu_external),
        "external_required_new_with_usable_ms2": len(panredu_external & usable_ms2_ids),
        "external_required_new_without_usable_ms2": len(panredu_external - usable_ms2_ids),
        "note": "No complete per-feature NEG Pan-ReDU table located; no cached reuse claimed.",
    }
    queue_rows.append(
        {
            "stage": "Pan-ReDU",
            "strict_ids": len(strict_ids),
            "exact_cache_reuse": 0,
            "retained_ids": len(retained_ids),
            "reassigned_ids": len(reassigned_ids),
            "new_ids": len(new_ids),
            "new_ids_with_cached_evidence": 0,
            "external_required_new_ids": len(panredu_external),
            "external_required_new_with_usable_ms2": len(panredu_external & usable_ms2_ids),
            "external_required_new_without_usable_ms2": len(panredu_external - usable_ms2_ids),
        }
    )
    queue = {
        "release_id": RELEASE_ID,
        "polarity": "NEG",
        "status": "external_approval_required_not_submitted",
        "strict_atlas_rows": int(len(strict_atlas)),
        "strict_feature_ids": int(len(strict_ids)),
        "historical_nonviral_feature_ids": int(len(historical_ids)),
        "retained_exact_feature_ids": int(len(retained_ids)),
        "new_feature_ids": int(len(new_ids)),
        "taxonomy_reassigned_retained_feature_ids": int(len(reassigned_ids)),
        "mgf_exact_strict_feature_ids_found": int(len(mgf_found_ids)),
        "usable_ms2_strict_feature_ids": int(len(usable_ms2_ids)),
        "new_feature_ids_with_usable_ms2": int(len(new_with_ms2)),
        "new_feature_ids_without_usable_ms2": int(len(new_without_ms2)),
        "mgf": {key: value for key, value in mgf.items() if key not in {"target_feature_ids_found", "target_peak_bearing_feature_ids"}},
        "stages": stage_reports,
        "queue_rows": queue_rows,
        "external_jobs": {
            stage: {
                "submitted": False,
                "approval_required": True,
                "new_features_requiring_external_work": report["external_required_new_ids"],
                "new_features_with_usable_ms2": report["external_required_new_with_usable_ms2"],
                "new_features_without_usable_ms2": report["external_required_new_without_usable_ms2"],
            }
            for stage, report in stage_reports.items()
        },
        "reuse_rule": "Exact feature ID plus compatible consensus-MGF identity/provenance; m/z and RT alone never create reuse.",
    }
    pd.DataFrame(queue_rows).to_csv(output_dir / "strict_annotation_cache_coverage.csv", index=False)
    write_json(output_dir / "strict_annotation_queue_manifest.json", queue)
    return queue


def queue_manifest(gate: dict, coverage: dict) -> dict:
    blocked = gate["status"] != "pass"
    return {
        "release_id": RELEASE_ID,
        "polarity": "NEG",
        "submission_allowed": False,
        "submitted": False,
        "status": "blocked_historical_composite_reproduction" if blocked else "ready_for_local_strict_build_only",
        "reason": "External annotation jobs remain approval-gated and are not submitted by this command.",
        "strict_feature_set_available": not blocked,
        "external_queue_sizes": {
            stage: {
                "status": "not_stageable_until_strict_atlas_is_accepted" if blocked else "not_submitted",
                "count": None if blocked else 0,
            }
            for stage in ["diagnostic_annotation", "SIRIUS_CANOPUS", "DreaMS", "MS2LDA", "fastMASST", "Pan-ReDU"]
        },
        "feature_partition": {
            "strict_changed_or_new_features": None,
            "strict_changed_or_new_with_usable_ms2": None,
            "strict_changed_or_new_without_usable_ms2": None,
            "reason": "The strict feature set is not accepted until the historical NEG Composite gate passes.",
        },
        "historical_cache_coverage_artifact": "historical_neg_annotation_cache_coverage.json",
        "cache_coverage_scope": coverage["scope"],
    }


def run_provisional_strict_composite(
    features: pd.DataFrame,
    samples: list[str],
    matrix,
    strict_meta: pd.DataFrame,
    taxonomy: dict,
    output_dir: Path,
) -> dict:
    """Run strict Composite scoring without promoting it past the gate."""
    platinum, silver, combined = engine.composite_select(
        features, samples, matrix, strict_meta, "strict NEG 16-phylum labels (provisional)"
    )
    if combined["feature_id"].duplicated().any():
        raise ValueError("Provisional strict NEG Composite output contains duplicate feature IDs")

    kingdom_map = engine.phylum_kingdom_map(strict_meta)
    rows = []
    for phylum in taxonomy["analysis_phyla"]:
        part = combined[combined["phylum"].eq(phylum)]
        rows.append(
            {
                "phylum": phylum,
                "kingdom": kingdom_map.get(phylum, ""),
                "strict_analysis_samples": int((strict_meta["phylum"] == phylum).sum()),
                "platinum": int((part["biomarker_quality"] == "Platinum").sum()),
                "diamond": int((part.get("biomarker_tier", pd.Series(dtype=str)) == "Diamond").sum()),
                "silver": int((part["biomarker_quality"] == "Silver").sum()),
                "total_composite": int(len(part)),
            }
        )
    counts = pd.DataFrame(rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    platinum.to_csv(output_dir / "provisional_strict_composite_platinum_NEG.csv", index=False)
    silver.to_csv(output_dir / "provisional_strict_composite_silver_NEG.csv", index=False)
    combined.to_csv(output_dir / "provisional_strict_composite_combined_NEG.csv", index=False)
    counts.to_csv(output_dir / "provisional_strict_composite_counts_by_phylum_NEG.csv", index=False)
    summary = {
        "release_id": RELEASE_ID,
        "polarity": "NEG",
        "status": "provisional_not_release_ready",
        "historical_gate_required": True,
        "atlas_rows": int(len(combined)),
        "platinum_rows": int(len(platinum)),
        "diamond_rows": int((platinum["biomarker_tier"] == "Diamond").sum()),
        "silver_rows": int(len(silver)),
        "distinct_feature_ids": int(combined["feature_id"].nunique()),
        "analysis_phyla": taxonomy["analysis_phyla"],
        "note": "Strict Composite scoring under the locked NCBI release; historical NEG Composite reproduction remains unresolved, so these values must not feed figures, tables, manuscript, or external queues.",
    }
    write_json(output_dir / "provisional_strict_composite_summary.json", summary)
    return summary


def run_strict_local_analysis(
    features: pd.DataFrame,
    samples: list[str],
    matrix,
    strict_meta: pd.DataFrame,
    taxonomy: dict,
    output_dir: Path,
) -> dict:
    """Build the accepted local strict NEG Composite + IndVal atlas after the gate."""
    platinum, silver, composite = engine.composite_select(
        features, samples, matrix, strict_meta, "strict NEG 16-phylum labels"
    )
    indval_pairs, indval_unique = engine.run_indval(
        features, samples, matrix, strict_meta, "strict NEG 16-phylum labels"
    )
    if composite["feature_id"].duplicated().any():
        raise ValueError("Strict NEG Composite output contains duplicate feature IDs")
    if indval_unique["feature_id"].duplicated().any():
        raise ValueError("Strict NEG IndVal output contains duplicate feature IDs")

    indval_only = indval_unique[
        ~indval_unique["feature_id"].isin(set(composite["feature_id"]))
    ].copy()
    indval_only = indval_only.merge(features, on="feature_id", how="left")
    indval_only["biomarker_quality"] = "IndVal"
    indval_only["biomarker_tier"] = "CrossBatch_Consensus"
    indval_only["discovery_method"] = "indval_consensus"
    composite = composite.copy()
    composite["discovery_method"] = "composite"
    strict_atlas = pd.concat([composite, indval_only], ignore_index=True, sort=False)
    if strict_atlas["feature_id"].duplicated().any():
        raise ValueError("Strict NEG atlas contains duplicate feature IDs")

    kingdom_map = engine.phylum_kingdom_map(strict_meta)
    rows = []
    for phylum in taxonomy["analysis_phyla"]:
        composite_part = composite[composite["phylum"].eq(phylum)]
        indval_part = indval_unique[indval_unique["phylum"].eq(phylum)]
        indval_only_part = indval_only[indval_only["phylum"].eq(phylum)]
        rows.append(
            {
                "phylum": phylum,
                "kingdom": kingdom_map.get(phylum, ""),
                "strict_analysis_samples": int((strict_meta["phylum"] == phylum).sum()),
                "composite_platinum": int(len(platinum[platinum["phylum"].eq(phylum)])),
                "composite_diamond": int(
                    (platinum.loc[platinum["phylum"].eq(phylum), "biomarker_tier"] == "Diamond").sum()
                ),
                "composite_silver": int(len(silver[silver["phylum"].eq(phylum)])),
                "composite_total": int(len(composite_part)),
                "indval_pair_rows": int(len(indval_pairs[indval_pairs["phylum"].eq(phylum)])),
                "indval_unique_features": int(len(indval_part)),
                "indval_only_features": int(len(indval_only_part)),
                "strict_atlas_total": int(len(strict_atlas[strict_atlas["phylum"].eq(phylum)])),
            }
        )
    counts = pd.DataFrame(rows)
    output_dir.mkdir(parents=True, exist_ok=True)
    platinum.to_csv(output_dir / "strict_composite_platinum_NEG.csv", index=False)
    silver.to_csv(output_dir / "strict_composite_silver_NEG.csv", index=False)
    composite.to_csv(output_dir / "strict_composite_combined_NEG.csv", index=False)
    indval_pairs.to_csv(output_dir / "strict_indval_pairs_NEG.csv", index=False)
    indval_unique.to_csv(output_dir / "strict_indval_unique_NEG.csv", index=False)
    strict_atlas.to_csv(output_dir / "strict_atlas_NEG.csv", index=False)
    counts.to_csv(output_dir / "strict_counts_by_phylum_NEG.csv", index=False)
    summary = {
        "release_id": RELEASE_ID,
        "polarity": "NEG",
        "status": "local_strict_complete_external_annotation_pending",
        "historical_gate_required": True,
        "atlas_rows": int(len(strict_atlas)),
        "composite_rows": int(len(composite)),
        "composite_platinum_rows": int(len(platinum)),
        "composite_diamond_rows": int((platinum["biomarker_tier"] == "Diamond").sum()),
        "composite_silver_rows": int(len(silver)),
        "indval_pair_rows": int(len(indval_pairs)),
        "indval_unique_features": int(len(indval_unique)),
        "indval_only_features": int(len(indval_only)),
        "distinct_feature_ids": int(strict_atlas["feature_id"].nunique()),
        "analysis_phyla": taxonomy["analysis_phyla"],
        "note": "Local strict NEG biomarker discovery completed under the locked NCBI release. External annotation jobs remain unsubmitted and approval-gated.",
    }
    write_json(output_dir / "strict_neg_summary.json", summary)
    return summary


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--taxonomy-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--run-provisional-strict-composite",
        action="store_true",
        help="Run strict NEG Composite scoring for diagnostic inspection while keeping it non-release.",
    )
    parser.add_argument(
        "--historical-recovery-root",
        type=Path,
        default=None,
        help="Extracted historical NEG recovery package containing the original Composite producer and output.",
    )
    parser.add_argument(
        "--run-strict",
        action="store_true",
        help="After the historical gate passes, write the strict NEG Composite + IndVal atlas.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    source_root = args.source_root.resolve()
    taxonomy_dir = args.taxonomy_dir.resolve()
    output_dir = args.output_dir.resolve()
    recovery_root = args.historical_recovery_root.resolve() if args.historical_recovery_root else None
    paths = source_paths(source_root, recovery_root)
    required = ["consensus", "historical_metadata", "historical_atlas", "historical_indval"]
    missing = [name for name in required if not paths[name].is_file()]
    if missing:
        raise FileNotFoundError(f"Missing required NEG audit inputs: {missing}")

    strict_meta, taxonomy = load_locked_taxonomy(taxonomy_dir)
    historical_meta = (
        read_recovery_metadata(paths["historical_metadata"])
        if recovery_root is not None
        else read_metadata(paths["historical_metadata"])
    )
    configure_engine(paths["consensus"])
    features, samples, matrix = engine.load_consensus(historical_meta, strict_meta)
    gate, atlas, observed_composite = historical_gate(
        paths, taxonomy, historical_meta, features, samples, matrix, output_dir, recovery_root
    )
    coverage = cache_coverage(paths, atlas, output_dir)
    write_json(output_dir / "annotation_queue_manifest.json", queue_manifest(gate, coverage))
    provisional_composite = None
    if args.run_provisional_strict_composite:
        provisional_composite = run_provisional_strict_composite(
            features, samples, matrix, strict_meta, taxonomy, output_dir
        )
    strict_result = None
    strict_queue = None
    if args.run_strict and gate["status"] == "pass":
        strict_result = run_strict_local_analysis(
            features, samples, matrix, strict_meta, taxonomy, output_dir
        )
        strict_queue = strict_annotation_queue(
            paths,
            pd.read_csv(output_dir / "strict_atlas_NEG.csv", low_memory=False),
            atlas,
            output_dir,
        )

    manifest = {
        "release_id": RELEASE_ID,
        "polarity": "NEG",
        "status": (
            "blocked"
            if gate["status"] != "pass"
            else "local_strict_complete_external_annotation_pending"
            if strict_result is not None
            else "historical_gate_passed_local_strict_pending"
        ),
        "taxonomy_dir": str(taxonomy_dir),
        "analysis_phyla": taxonomy["analysis_phyla"],
        "negative_core_candidate_samples": int(taxonomy["negative_core_samples"]),
        "strict_analysis_samples_16_phyla": int(len(strict_meta)),
        "strict_analysis_counts_16_phyla": strict_meta["phylum"].value_counts().sort_index().to_dict(),
        "source_hashes": {name: fingerprint(path) for name, path in paths.items() if path.is_file()},
        "producer": {
            "historical_selection_engine": gate["selection_engine"],
            "strict_selection_engine": "working_corrections.figure_2a.rerun_fig2a",
            "pipeline_script_sha256": sha256(Path(__file__).resolve()),
        },
        "thresholds": {
            "composite": {
                "historical_weights": {
                    "specificity": 0.30,
                    "detection_rate": 0.25,
                    "log2fc_normalized": 0.20,
                    "batch_rate": 0.15,
                    "exclusivity": 0.10,
                },
                "strict_engine": "shared POS-tested selector generalized to NEG; strict sample groups are recomputed",
                "tier1_score": 60,
                "tier1_detection_rate": 0.30,
                "tier1_batch_confirm_min": 2,
                "platinum_background_detections": 0,
                "platinum_detection_rate": 0.50,
                "platinum_log2fc": 3.0,
                "platinum_mean_target": 10000,
                "silver_detection_rate": 0.40,
                "silver_log2fc": 2.0,
            },
            "indval": {
                "minimum_samples_per_phylum": 3,
                "minimum_batch_samples": 3,
                "minimum_feature_detections": 2,
                "within_batch_indval": 0.5,
                "cross_batch_consensus": 2,
            },
        },
        "seeds": {"selection": None, "note": "Deterministic matrix operations; no stochastic selection step."},
        "denominators": {
            "negative_core_candidate_samples": int(taxonomy["negative_core_samples"]),
            "strict_analysis_samples_16_phyla": int(len(strict_meta)),
            "historical_atlas_rows": gate["historical_atlas_rows"],
            "historical_nonviral_atlas_rows": gate["historical_atlas_nonviral_rows"],
        },
        "historical_gate": gate,
        "provisional_strict_composite": provisional_composite,
        "strict_local_analysis": strict_result,
        "strict_annotation_queue": strict_queue,
        "output_contract": {
            "strict_atlas_written": strict_result is not None,
            "reason": (
                "Strict NEG Composite + IndVal atlas written after historical gate."
                if strict_result is not None
                else "Strict atlas is withheld until the historical NEG gate passes."
            ),
        },
    }
    write_json(output_dir / "release_manifest.json", manifest)
    print(json.dumps({
        "status": manifest["status"],
        "historical_composite_exact": gate["composite_reproduction_exact"],
        "historical_indval_exact": gate["indval_reproduction_exact"],
        "historical_atlas_rows": gate["historical_atlas_rows"],
        "historical_nonviral_rows": gate["historical_atlas_nonviral_rows"],
        "mgf_peak_bearing_nonviral": coverage["mgf_peak_bearing_nonviral_feature_ids"],
    }, indent=2))
    return 0 if gate["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
