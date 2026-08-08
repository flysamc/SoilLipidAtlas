#!/usr/bin/env python3
"""Run review-only strict POS/NEG MS2LDA phylum enrichment by exact feature ID."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[2]
REL = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1"
THRESHOLD = 0.05


def bh_adjust(pvalues: pd.Series) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted_ranked = np.minimum.accumulate(
        (ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1]
    )[::-1]
    adjusted = np.empty_like(adjusted_ranked)
    adjusted[order] = np.clip(adjusted_ranked, 0, 1)
    return adjusted


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def record(path: Path) -> dict:
    return {
        "path": str(path.relative_to(ROOT)).replace("\\", "/"),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def mgf_titles(path: Path) -> set[str]:
    titles = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("TITLE="):
                titles.add(line.rstrip("\r\n")[6:])
    return titles


def run(polarity: str) -> None:
    if polarity == "POS":
        model = REL / "biomarker_discovery/external_annotation_package/ms2lda_strict_atlas_results"
        atlas_path = REL / "biomarker_discovery/atlas_pos_strict.csv"
        usable_path = REL / "biomarker_discovery/external_annotation_package/figure2a_strict_atlas_with_usable_ms2.mgf"
        usable_ids = mgf_titles(usable_path)
        out = REL / "biomarker_discovery/ms2lda_enrichment_review_only"
    else:
        model = REL / "annotation_recovery_neg/ms2lda/strict_full_run1/results_full"
        atlas_path = REL / "biomarker_discovery_neg/strict_atlas_NEG.csv"
        usable_path = REL / "annotation_recovery_neg/strict_feature_spectrum_ledger.csv"
        usable = pd.read_csv(usable_path, usecols=["feature_id", "has_usable_ms2"])
        usable_ids = set(usable.loc[usable["has_usable_ms2"].astype(str).str.lower().eq("true"), "feature_id"].astype(str))
        out = REL / "annotation_recovery_neg/ms2lda_enrichment_review_only"
    out.mkdir(parents=True, exist_ok=True)

    doc_path = model / "doc_topic_matrix.csv"
    words_path = model / "motif_top_words.csv"
    docs = pd.read_csv(doc_path, low_memory=False)
    atlas = pd.read_csv(atlas_path, usecols=["feature_id", "phylum", "kingdom"], low_memory=False)
    docs["feature_id"] = docs["feature_id"].astype(str)
    atlas["feature_id"] = atlas["feature_id"].astype(str)
    if docs["feature_id"].duplicated().any() or atlas["feature_id"].duplicated().any():
        raise ValueError(f"{polarity}: duplicate feature IDs")
    strict_ids = set(atlas["feature_id"])
    model_ids = set(docs["feature_id"])
    if not model_ids <= strict_ids:
        raise ValueError(f"{polarity}: model contains non-strict feature IDs")
    if not model_ids <= usable_ids:
        raise ValueError(f"{polarity}: model contains a feature outside the strict usable-MS2 set")
    docs = docs.merge(atlas, on="feature_id", how="left", validate="one_to_one")
    if docs["phylum"].isna().any():
        raise ValueError(f"{polarity}: exact-ID phylum mapping incomplete")

    motifs = [column for column in docs if column.startswith("motif_")]
    phylum_counts = docs["phylum"].value_counts()
    phyla = sorted(phylum_counts[phylum_counts >= 3].index)
    total = len(docs)
    results = []
    for motif in motifs:
        present = pd.to_numeric(docs[motif], errors="coerce").fillna(0).ge(THRESHOLD)
        n_motif = int(present.sum())
        if n_motif < 3:
            continue
        for phylum in phyla:
            in_phylum = docs["phylum"].eq(phylum)
            a = int((present & in_phylum).sum())
            if a < 1:
                continue
            b = int(in_phylum.sum()) - a
            c = n_motif - a
            d = total - a - b - c
            odds, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
            observed = a / n_motif
            expected = int(in_phylum.sum()) / total
            results.append(
                {
                    "motif": motif,
                    "phylum": phylum,
                    "kingdom": docs.loc[in_phylum, "kingdom"].iloc[0],
                    "n_motif_in_phylum": a,
                    "n_motif_total": n_motif,
                    "n_phylum_total": int(in_phylum.sum()),
                    "observed_fraction": observed,
                    "expected_fraction": expected,
                    "enrichment_ratio": observed / expected,
                    "odds_ratio": odds,
                    "pvalue": pvalue,
                }
            )
    enrichment = pd.DataFrame(results)
    if enrichment.empty:
        raise ValueError(f"{polarity}: no enrichment tests generated")
    enrichment["padj_bh_all_tests"] = bh_adjust(enrichment["pvalue"])
    enrichment["significant"] = enrichment["padj_bh_all_tests"].lt(0.05) & enrichment["enrichment_ratio"].ge(2)
    enrichment = enrichment.sort_values(["padj_bh_all_tests", "pvalue", "motif", "phylum"])

    words = pd.read_csv(words_path)
    character = []
    for motif in motifs:
        tests = enrichment[enrichment["motif"].eq(motif)].sort_values("pvalue")
        significant = tests[tests["significant"]]
        primary = significant.iloc[0] if len(significant) else (tests.iloc[0] if len(tests) else None)
        top = words[words["motif"].eq(motif)].sort_values("probability", ascending=False).head(5)
        character.append(
            {
                "motif": motif,
                "primary_phylum": "" if primary is None else primary["phylum"],
                "primary_kingdom": "" if primary is None else primary["kingdom"],
                "enrichment_ratio": 0 if primary is None else primary["enrichment_ratio"],
                "fdr_bh_all_tests": 1 if primary is None else primary["padj_bh_all_tests"],
                "n_enriched_phyla": len(significant),
                "other_enriched_phyla": ";".join(significant.iloc[1:]["phylum"].astype(str)),
                "top_words": ";".join(top["word"].astype(str)),
            }
        )
    character = pd.DataFrame(character)
    matrix = docs.groupby("phylum", sort=True)[motifs].mean()

    enrichment_path = out / "motif_phylum_enrichment.csv"
    character_path = out / "ms2lda_motif_characterization.csv"
    matrix_path = out / "phylum_motif_probability_matrix.csv"
    enrichment.to_csv(enrichment_path, index=False)
    character.to_csv(character_path, index=False)
    matrix.to_csv(matrix_path)

    manifest = {
        "schema_version": 1,
        "stage_id": f"strict_{polarity.lower()}_ms2lda_phylum_enrichment",
        "taxonomy_release": "ncbi-phylum-2026-08-04-v1",
        "polarity": polarity,
        "status": "review_only_complete",
        "method": {
            "motif_probability_threshold": THRESHOLD,
            "minimum_model_documents_per_phylum": 3,
            "test": "one-sided Fisher exact",
            "multiple_testing": "Benjamini-Hochberg across all emitted motif-phylum tests",
            "significance": "adjusted p < 0.05 and enrichment ratio >= 2",
            "mapping": "exact feature_id only",
        },
        "counts": {
            "strict_atlas_features": len(atlas),
            "strict_usable_ms2_inputs": len(usable_ids),
            "strict_without_usable_ms2": len(atlas) - len(usable_ids),
            "model_documents": len(docs),
            "model_preprocessing_exclusions": len(usable_ids) - len(docs),
            "motifs": len(motifs),
            "tested_phyla": len(phyla),
            "tests": len(enrichment),
            "significant_tests": int(enrichment["significant"].sum()),
            "motifs_with_significant_phylum": int(character["n_enriched_phyla"].gt(0).sum()),
        },
        "inputs": [record(doc_path), record(words_path), record(atlas_path), record(usable_path)],
        "outputs": [record(enrichment_path), record(character_path), record(matrix_path)],
        "release_boundary": [
            "This is a strict candidate analysis and is not connected to the legacy supplementary MS2LDA figure or table.",
            "Motif biological labels remain unreviewed; no manuscript or submitted-document update is authorized.",
        ],
    }
    (out / "stage_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"polarity": polarity, **manifest["counts"]}, indent=2))


def main() -> None:
    for polarity in ["POS", "NEG"]:
        run(polarity)


if __name__ == "__main__":
    main()
