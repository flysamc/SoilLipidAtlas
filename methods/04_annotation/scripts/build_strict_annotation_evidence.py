#!/usr/bin/env python3
"""Build exact-feature-ID review ledgers for current POS/NEG annotation evidence.

The output is intentionally an evidence audit, not a confidence-tier upgrader.
SIRIUS formulas, CANOPUS classes, DreaMS neighbours, family propagation, and
ArchLips matches remain separate evidence types with explicit provenance.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
RELEASE = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_ids(path: Path, column: str, sep: str = ",") -> set[str]:
    return set(pd.read_csv(path, sep=sep, usecols=[column])[column].dropna().astype(str))


def mgf_titles(path: Path) -> set[str]:
    result: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("TITLE="):
                result.add(line.rstrip("\r\n")[6:])
    return result


def bool_series(series: pd.Series) -> pd.Series:
    return series.astype(str).str.lower().isin({"true", "1", "1.0", "yes"})


def load_rank1(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df["mappingFeatureId"] = df["mappingFeatureId"].astype(str)
    if "formulaRank" in df:
        df = df[pd.to_numeric(df["formulaRank"], errors="coerce").eq(1)]
    return df.drop_duplicates("mappingFeatureId", keep="first")


def load_top_structure(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df["mappingFeatureId"] = df["mappingFeatureId"].astype(str)
    if "structurePerIdRank" in df:
        df = df[pd.to_numeric(df["structurePerIdRank"], errors="coerce").eq(1)]
    return df.drop_duplicates("mappingFeatureId", keep="first")


def best_dreams(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    df["feature_id"] = df["feature_id"].astype(str)
    df["DreaMS_similarity"] = pd.to_numeric(df["DreaMS_similarity"], errors="coerce")
    if "topk" in df:
        top1 = df[pd.to_numeric(df["topk"], errors="coerce").eq(1)].copy()
        if not top1.empty:
            df = top1
    return df.sort_values("DreaMS_similarity", ascending=False).drop_duplicates("feature_id")


def source_record(path: Path) -> dict:
    return {
        "path": str(path.relative_to(ROOT)).replace("\\", "/"),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def build(polarity: str) -> None:
    ann = RELEASE / "annotation"
    if polarity == "POS":
        atlas_path = RELEASE / "biomarker_discovery/atlas_pos_strict.csv"
        atlas = pd.read_csv(atlas_path, low_memory=False)
        coverage_path = RELEASE / "biomarker_discovery/annotation_coverage_by_feature.csv"
        coverage = pd.read_csv(coverage_path, low_memory=False)
        new_ids = set(
            coverage.loc[bool_series(coverage["new_biomarker_needing_full_pipeline"]), "feature_id"].astype(str)
        )
        atlas["partition"] = np.where(atlas["feature_id"].astype(str).isin(new_ids), "new_or_changed", "retained_exact")
        strict_ms2_path = RELEASE / "biomarker_discovery/external_annotation_package/figure2a_strict_atlas_with_usable_ms2.mgf"
        usable_ms2_ids = mgf_titles(strict_ms2_path)
        sirius_new_dir = RELEASE / "biomarker_discovery/external_annotation_results/sirius_canopus/new_features_mzle850_run1_summaries"
        sirius_cache_dir = ROOT / "external/SOILMASS_PRODUCER_RECOVERY_2026-08-04/payload/core/workspace/Dreams/results/sirius_atlas_full_v2"
        dreams_new_path = RELEASE / "biomarker_discovery/external_annotation_package/dreams_new_features_results/figure2a_new_features_dreams_top10.tsv"
        dreams_cache_path = ROOT / "external/SOILMASS_PRODUCER_RECOVERY_2026-08-04/payload/core/workspace/Dreams/results/dreams_gpu_search/atlas_dreams_search_results.tsv"
        ms2lda_path = RELEASE / "biomarker_discovery/external_annotation_package/ms2lda_strict_atlas_results/doc_topic_matrix.csv"
        out_dir = RELEASE / "biomarker_discovery/external_annotation_results/integration_review_only"
    else:
        atlas_path = RELEASE / "biomarker_discovery_neg/strict_atlas_NEG.csv"
        atlas = pd.read_csv(atlas_path, low_memory=False)
        coverage_path = RELEASE / "annotation_recovery_neg/strict_feature_spectrum_ledger.csv"
        coverage = pd.read_csv(coverage_path, low_memory=False)
        partition = coverage[["feature_id", "partition", "has_usable_ms2"]].copy()
        atlas = atlas.merge(partition[["feature_id", "partition"]], on="feature_id", how="left", validate="one_to_one")
        usable_ms2_ids = set(coverage.loc[bool_series(coverage["has_usable_ms2"]), "feature_id"].astype(str))
        strict_ms2_path = coverage_path
        sirius_new_dir = RELEASE / "annotation_recovery_neg/sirius_canopus/neg_new_changed_mzle850_run1_summaries"
        sirius_cache_dir = ROOT / "external/SOILMASS_PRODUCER_RECOVERY_2026-08-04/payload/core/workspace/Dreams/results/sirius_atlas_neg"
        dreams_new_path = RELEASE / "annotation_recovery_neg/dreams/neg_new_changed_dreams_top10.tsv"
        dreams_cache_path = RELEASE / "annotation_recovery_neg/dreams/neg_retained_dreams_top10.tsv"
        ms2lda_path = RELEASE / "annotation_recovery_neg/ms2lda/strict_full_run1/results_full/doc_topic_matrix.csv"
        out_dir = RELEASE / "annotation_recovery_neg/integration_review_only"

    out_dir.mkdir(parents=True, exist_ok=True)
    atlas["feature_id"] = atlas["feature_id"].astype(str)
    strict_ids = set(atlas["feature_id"])
    if len(atlas) != len(strict_ids):
        raise ValueError(f"{polarity}: strict atlas feature IDs are not unique")

    current_path = ann / f"step2_harmonization/harmonised_annotations_{polarity.lower()}.csv"
    current = pd.read_csv(current_path, low_memory=False)
    current["feature_id"] = current["feature_id"].astype(str)
    current = current.drop_duplicates("feature_id")
    ledger = atlas[[c for c in ["feature_id", "phylum", "kingdom", "discovery_method", "partition"] if c in atlas]].copy()
    ledger = ledger.merge(
        current[["feature_id", "annotation_tier", "annotation_level", "annotation_source"]],
        on="feature_id", how="left", validate="one_to_one",
    )
    ledger["has_usable_ms2"] = ledger["feature_id"].isin(usable_ms2_ids)

    retained_ids = set(ledger.loc[ledger["partition"].eq("retained_exact"), "feature_id"])
    new_ids = strict_ids - retained_ids

    formula_cache_path = sirius_cache_dir / "formula_identifications.tsv"
    formula_new_path = sirius_new_dir / "formula_identifications.tsv"
    canopus_cache_path = sirius_cache_dir / "canopus_formula_summary.tsv"
    canopus_new_path = sirius_new_dir / "canopus_formula_summary.tsv"
    csi_cache_path = sirius_cache_dir / "structure_identifications.tsv"
    csi_new_path = sirius_new_dir / "structure_identifications.tsv"
    denovo_cache_path = sirius_cache_dir / "denovo_structure_identifications.tsv"
    denovo_new_path = sirius_new_dir / "denovo_structure_identifications.tsv"
    formula_cache = load_rank1(formula_cache_path)
    formula_new = load_rank1(formula_new_path)
    canopus_cache = load_rank1(canopus_cache_path)
    canopus_new = load_rank1(canopus_new_path)
    csi_cache = load_top_structure(csi_cache_path)
    csi_new = load_top_structure(csi_new_path)
    denovo_cache = load_top_structure(denovo_cache_path)
    denovo_new = load_top_structure(denovo_new_path)

    # --- Gap-fill sources: SIRIUS/CANOPUS/CSI/de-novo for the eligible-unsubmitted
    # strict biomarkers (POS 2,209 / NEG 80). LISC jobs 5862224 (formula/CANOPUS/
    # de-novo) + 5869613 (CSI structures over BIO,PUBCHEM,HMDB,... union), 2026-08-12/13.
    # Disjoint from the submitted (cache union new) universe by construction;
    # mappingFeatureId == strict feature_id.
    gapfill_dir = RELEASE / f"biomarker_discovery/external_annotation_package/gapfill_2026-08-12/results/gapfill_{polarity}"
    formula_gapfill_path = gapfill_dir / "formula_identifications.tsv"
    canopus_gapfill_path = gapfill_dir / "canopus_formula_summary.tsv"
    csi_gapfill_path = gapfill_dir / "structure_identifications.tsv"
    denovo_gapfill_path = gapfill_dir / "denovo_structure_identifications.tsv"
    formula_gapfill = load_rank1(formula_gapfill_path)
    canopus_gapfill = load_rank1(canopus_gapfill_path)
    csi_gapfill = load_top_structure(csi_gapfill_path)
    denovo_gapfill = load_top_structure(denovo_gapfill_path)

    formula_cache_ids = set(formula_cache["mappingFeatureId"]) & retained_ids
    formula_new_ids = set(formula_new["mappingFeatureId"]) & new_ids
    canopus_cache_ids = set(canopus_cache["mappingFeatureId"]) & retained_ids
    canopus_new_ids = set(canopus_new["mappingFeatureId"]) & new_ids
    csi_cache_ids = set(csi_cache["mappingFeatureId"]) & retained_ids
    csi_new_ids = set(csi_new["mappingFeatureId"]) & new_ids
    denovo_cache_ids = set(denovo_cache["mappingFeatureId"]) & retained_ids
    denovo_new_ids = set(denovo_new["mappingFeatureId"]) & new_ids
    formula_gapfill_ids = set(formula_gapfill["mappingFeatureId"]) & strict_ids
    canopus_gapfill_ids = set(canopus_gapfill["mappingFeatureId"]) & strict_ids
    csi_gapfill_ids = set(csi_gapfill["mappingFeatureId"]) & strict_ids
    denovo_gapfill_ids = set(denovo_gapfill["mappingFeatureId"]) & strict_ids

    ledger["has_sirius_formula"] = ledger["feature_id"].isin(formula_cache_ids | formula_new_ids | formula_gapfill_ids)
    ledger["sirius_formula_provenance"] = np.select(
        [ledger["feature_id"].isin(formula_new_ids), ledger["feature_id"].isin(formula_gapfill_ids), ledger["feature_id"].isin(formula_cache_ids)],
        ["fresh_new_or_changed", "gapfill_2026_08_12", "exact_retained_cache"], default="none",
    )
    formulas = pd.concat(
        [
            formula_cache[formula_cache["mappingFeatureId"].isin(formula_cache_ids)],
            formula_new[formula_new["mappingFeatureId"].isin(formula_new_ids)],
            formula_gapfill[formula_gapfill["mappingFeatureId"].isin(formula_gapfill_ids)],
        ], ignore_index=True,
    ).drop_duplicates("mappingFeatureId")
    ledger = ledger.merge(
        formulas[["mappingFeatureId", "molecularFormula", "adduct", "SiriusScoreNormalized"]].rename(
            columns={"mappingFeatureId": "feature_id", "molecularFormula": "sirius_molecular_formula", "adduct": "sirius_adduct", "SiriusScoreNormalized": "sirius_score_normalized"}
        ), on="feature_id", how="left", validate="one_to_one",
    )

    ledger["has_canopus_class"] = ledger["feature_id"].isin(canopus_cache_ids | canopus_new_ids | canopus_gapfill_ids)
    ledger["canopus_provenance"] = np.select(
        [ledger["feature_id"].isin(canopus_new_ids), ledger["feature_id"].isin(canopus_gapfill_ids), ledger["feature_id"].isin(canopus_cache_ids)],
        ["fresh_new_or_changed", "gapfill_2026_08_12", "exact_retained_cache"], default="none",
    )
    canopus = pd.concat(
        [
            canopus_cache[canopus_cache["mappingFeatureId"].isin(canopus_cache_ids)],
            canopus_new[canopus_new["mappingFeatureId"].isin(canopus_new_ids)],
            canopus_gapfill[canopus_gapfill["mappingFeatureId"].isin(canopus_gapfill_ids)],
        ], ignore_index=True,
    ).drop_duplicates("mappingFeatureId")
    canopus_cols = ["mappingFeatureId", "NPC#pathway", "NPC#superclass", "NPC#class", "ClassyFire#most specific class"]
    ledger = ledger.merge(
        canopus[canopus_cols].rename(columns={"mappingFeatureId": "feature_id"}),
        on="feature_id", how="left", validate="one_to_one",
    )

    ledger["has_csi_fingerid_structure"] = ledger["feature_id"].isin(csi_cache_ids | csi_new_ids | csi_gapfill_ids)
    ledger["has_denovo_structure"] = ledger["feature_id"].isin(denovo_cache_ids | denovo_new_ids | denovo_gapfill_ids)

    dreams_cache = best_dreams(dreams_cache_path)
    dreams_new = best_dreams(dreams_new_path)
    dreams_cache = dreams_cache[dreams_cache["feature_id"].isin(retained_ids)]
    dreams_new = dreams_new[dreams_new["feature_id"].isin(new_ids)]
    dreams = pd.concat([dreams_cache, dreams_new], ignore_index=True).sort_values(
        "DreaMS_similarity", ascending=False
    ).drop_duplicates("feature_id")
    dreams_cache_ids = set(dreams_cache["feature_id"])
    dreams_new_ids = set(dreams_new["feature_id"])
    ledger["has_dreams_result"] = ledger["feature_id"].isin(dreams_cache_ids | dreams_new_ids)
    ledger["dreams_provenance"] = np.select(
        [ledger["feature_id"].isin(dreams_new_ids), ledger["feature_id"].isin(dreams_cache_ids)],
        ["fresh_new_or_changed", "exact_retained_cache"], default="none",
    )
    dream_cols = [c for c in ["feature_id", "DreaMS_similarity", "ref_msv_id", "ref_msv_species_resolved", "ref_name", "ref_smiles"] if c in dreams]
    ledger = ledger.merge(dreams[dream_cols], on="feature_id", how="left", validate="one_to_one")

    ms2lda_ids = exact_ids(ms2lda_path, "feature_id") & strict_ids
    ledger["in_strict_ms2lda_model"] = ledger["feature_id"].isin(ms2lda_ids)

    propagation_path = ann / f"step4_family_propagation/propagated_{polarity.lower()}_normalised.csv"
    archlips_path = ann / f"step8_archlips_rt_filtered/archlips_{polarity.lower()}_rt_screened.csv"
    propagation_ids = exact_ids(propagation_path, "feature_id") & strict_ids
    archlips = pd.read_csv(archlips_path, low_memory=False)
    archlips["feature_id"] = archlips["feature_id"].astype(str)
    archlips_ids = set(archlips.loc[~bool_series(archlips["rt_uncertain"]), "feature_id"]) & strict_ids
    ledger["has_family_propagation"] = ledger["feature_id"].isin(propagation_ids)
    ledger["has_rt_screened_archlips"] = ledger["feature_id"].isin(archlips_ids)

    existing = ledger["annotation_tier"].fillna("Unidentified").ne("Unidentified")
    ledger["external_evidence_depth"] = np.select(
        [
            ledger["has_csi_fingerid_structure"],
            ledger["has_denovo_structure"],
            ledger["has_sirius_formula"] & ledger["has_canopus_class"],
            ledger["has_sirius_formula"],
            ledger["has_canopus_class"],
            ledger["has_dreams_result"],
            ledger["in_strict_ms2lda_model"],
            ledger["has_usable_ms2"],
            existing,
        ],
        [
            "csi_fingerid_structure_candidate",
            "denovo_structure_candidate",
            "sirius_formula_plus_canopus_class",
            "sirius_formula_only",
            "canopus_class_only",
            "dreams_neighbour_only",
            "ms2lda_model_only",
            "usable_ms2_only",
            "existing_annotation_only",
        ],
        default="no_current_evidence",
    )

    evidence_path = out_dir / f"strict_{polarity.lower()}_annotation_evidence.csv"
    ledger.to_csv(evidence_path, index=False)

    depth = (
        ledger.groupby(["partition", "external_evidence_depth"], dropna=False)
        .size().rename("features").reset_index()
    )
    depth["polarity"] = polarity
    depth.to_csv(out_dir / "external_evidence_depth_counts.csv", index=False)

    rows = []
    for partition_name, subset in [("all", ledger), *list(ledger.groupby("partition"))]:
        rows.append(
            {
                "polarity": polarity,
                "partition": partition_name,
                "strict_features": len(subset),
                "usable_ms2": int(subset["has_usable_ms2"].sum()),
                "sirius_formula": int(subset["has_sirius_formula"].sum()),
                "canopus_class": int(subset["has_canopus_class"].sum()),
                "csi_fingerid_structure": int(subset["has_csi_fingerid_structure"].sum()),
                "denovo_structure": int(subset["has_denovo_structure"].sum()),
                "dreams_result": int(subset["has_dreams_result"].sum()),
                "ms2lda_model_document": int(subset["in_strict_ms2lda_model"].sum()),
                "family_propagation": int(subset["has_family_propagation"].sum()),
                "rt_screened_archlips": int(subset["has_rt_screened_archlips"].sum()),
            }
        )
    coverage_summary = pd.DataFrame(rows)
    coverage_summary.to_csv(out_dir / "tool_coverage_by_partition.csv", index=False)

    dream_rows = []
    for partition_name, subset in [("all", ledger), *list(ledger.groupby("partition"))]:
        scores = pd.to_numeric(subset["DreaMS_similarity"], errors="coerce").dropna()
        dream_rows.append(
            {
                "polarity": polarity,
                "partition": partition_name,
                "features_with_result": len(scores),
                "median_similarity": float(scores.median()) if len(scores) else None,
                "mean_similarity": float(scores.mean()) if len(scores) else None,
                "similarity_ge_0_3": int((scores >= 0.3).sum()),
                "similarity_ge_0_5": int((scores >= 0.5).sum()),
                "similarity_ge_0_7": int((scores >= 0.7).sum()),
                "similarity_ge_0_9": int((scores >= 0.9).sum()),
            }
        )
    pd.DataFrame(dream_rows).to_csv(out_dir / "dreams_similarity_summary.csv", index=False)

    validations = {
        "atlas_feature_ids_unique": len(atlas) == len(strict_ids),
        "ledger_rows_equal_atlas_rows": len(ledger) == len(atlas),
        "ledger_feature_ids_equal_atlas_ids": set(ledger["feature_id"]) == strict_ids,
        "all_sirius_formula_ids_are_exact_strict_ids": set(formulas["mappingFeatureId"]) <= strict_ids,
        "all_canopus_ids_are_exact_strict_ids": set(canopus["mappingFeatureId"]) <= strict_ids,
        "all_csi_structure_ids_are_exact_strict_ids": (csi_cache_ids | csi_new_ids) <= strict_ids,
        "all_denovo_structure_ids_are_exact_strict_ids": (denovo_cache_ids | denovo_new_ids) <= strict_ids,
        "all_dreams_ids_are_exact_strict_ids": set(dreams["feature_id"]) <= strict_ids,
        "all_ms2lda_ids_are_exact_strict_ids": ms2lda_ids <= strict_ids,
        "fresh_and_cached_formula_sets_disjoint": not (formula_cache_ids & formula_new_ids),
        "fresh_and_cached_canopus_sets_disjoint": not (canopus_cache_ids & canopus_new_ids),
        "fresh_and_cached_csi_sets_disjoint": not (csi_cache_ids & csi_new_ids),
        "fresh_and_cached_denovo_sets_disjoint": not (denovo_cache_ids & denovo_new_ids),
        "fresh_and_cached_dreams_sets_disjoint": not (dreams_cache_ids & dreams_new_ids),
        "all_gapfill_formula_ids_are_exact_strict_ids": set(formula_gapfill["mappingFeatureId"]) <= strict_ids,
        "all_gapfill_canopus_ids_are_exact_strict_ids": set(canopus_gapfill["mappingFeatureId"]) <= strict_ids,
        "all_gapfill_csi_ids_are_exact_strict_ids": csi_gapfill_ids <= strict_ids,
        "all_gapfill_denovo_ids_are_exact_strict_ids": denovo_gapfill_ids <= strict_ids,
        "gapfill_formula_disjoint_from_submitted": not (formula_gapfill_ids & (formula_cache_ids | formula_new_ids)),
        "gapfill_canopus_disjoint_from_submitted": not (canopus_gapfill_ids & (canopus_cache_ids | canopus_new_ids)),
        "gapfill_csi_disjoint_from_submitted": not (csi_gapfill_ids & (csi_cache_ids | csi_new_ids)),
        "gapfill_denovo_disjoint_from_submitted": not (denovo_gapfill_ids & (denovo_cache_ids | denovo_new_ids)),
    }
    manifest = {
        "schema_version": 1,
        "taxonomy_release": "ncbi-phylum-2026-08-04-v1",
        "polarity": polarity,
        "status": "review_only_local_evidence_integration_complete" if all(validations.values()) else "failed_validation",
        "matching_rule": "Exact feature ID only. No mass-only or retention-time-only reuse is permitted.",
        "counts": coverage_summary.iloc[0].to_dict(),
        "validations": validations,
        "inputs": [
            source_record(path) for path in [
                atlas_path, coverage_path, strict_ms2_path, current_path,
                formula_cache_path, formula_new_path, canopus_cache_path, canopus_new_path,
                csi_cache_path, csi_new_path, denovo_cache_path, denovo_new_path,
                dreams_cache_path, dreams_new_path, ms2lda_path, propagation_path, archlips_path,
                formula_gapfill_path, canopus_gapfill_path, csi_gapfill_path, denovo_gapfill_path,
            ]
        ],
        "outputs": [
            source_record(path) for path in [
                evidence_path,
                out_dir / "external_evidence_depth_counts.csv",
                out_dir / "tool_coverage_by_partition.csv",
                out_dir / "dreams_similarity_summary.csv",
            ]
        ],
        "release_boundary": [
            "This ledger records evidence availability; it does not convert SIRIUS, CANOPUS, DreaMS, or MS2LDA evidence into Gold/Silver/Bronze lipid-identification tiers.",
            "Historical Step 9 tier-transition rules and the exact submitted 564/640/12 action producer remain unrecovered, so those transitions are not claimed as reproduced.",
            "No manuscript, figure, table, or submission-source consumer is connected to these review-only outputs.",
        ],
    }
    manifest_path = out_dir / "stage_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    if manifest["status"].startswith("failed"):
        raise SystemExit(f"{polarity}: evidence integration validation failed")
    print(json.dumps({"polarity": polarity, "status": manifest["status"], "counts": manifest["counts"]}, indent=2))


def main() -> None:
    for polarity in ["POS", "NEG"]:
        build(polarity)


if __name__ == "__main__":
    main()
