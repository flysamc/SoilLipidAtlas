#!/usr/bin/env python3
"""Rerun strict positive-mode Figure 2A and audit annotation coverage.

This adapter uses the locked taxonomy freeze and the recovered historical
producer inputs.  It must reproduce the archived Composite and IndVal outputs
before it is allowed to generate the corrected strict atlas.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import posixpath
import re
import time
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")

from working_corrections.figure_2a import rerun_fig2a as recovered


XLSX_MAIN_NS = "http://schemas.openxmlformats.org/spreadsheetml/2006/main"
XLSX_DOC_REL_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
XLSX_PACKAGE_REL_NS = "http://schemas.openxmlformats.org/package/2006/relationships"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _xlsx_column_index(reference: str) -> int:
    match = re.match(r"([A-Z]+)", reference)
    if not match:
        raise ValueError(f"Invalid XLSX cell reference: {reference}")
    value = 0
    for character in match.group(1):
        value = value * 26 + ord(character) - ord("A") + 1
    return value - 1


def read_xlsx_table(path: Path, sheet_name: str) -> pd.DataFrame:
    """Read one value-only XLSX table using the Python standard library.

    The historical scoring workbook has no formulas in the four fields used by
    the producer. Avoiding an optional Excel engine keeps this reproduction
    gate portable while the workbook SHA-256 remains the authoritative input.
    """
    with zipfile.ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        sheet = next(
            (
                node
                for node in workbook.findall(f".//{{{XLSX_MAIN_NS}}}sheet")
                if node.attrib.get("name") == sheet_name
            ),
            None,
        )
        if sheet is None:
            raise ValueError(f"Sheet {sheet_name!r} not found in {path}")
        relationship_id = sheet.attrib[f"{{{XLSX_DOC_REL_NS}}}id"]
        relationships = ET.fromstring(
            archive.read("xl/_rels/workbook.xml.rels")
        )
        target = next(
            node.attrib["Target"]
            for node in relationships.findall(
                f".//{{{XLSX_PACKAGE_REL_NS}}}Relationship"
            )
            if node.attrib.get("Id") == relationship_id
        )
        worksheet_path = (
            target.lstrip("/")
            if target.startswith("/")
            else posixpath.normpath(posixpath.join("xl", target))
        )

        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in shared.findall(f"{{{XLSX_MAIN_NS}}}si"):
                shared_strings.append(
                    "".join(
                        node.text or ""
                        for node in item.iter(f"{{{XLSX_MAIN_NS}}}t")
                    )
                )

        worksheet = ET.fromstring(archive.read(worksheet_path))
        rows: list[dict[int, object]] = []
        for row in worksheet.findall(f".//{{{XLSX_MAIN_NS}}}row"):
            values: dict[int, object] = {}
            for cell in row.findall(f"{{{XLSX_MAIN_NS}}}c"):
                index = _xlsx_column_index(cell.attrib["r"])
                cell_type = cell.attrib.get("t")
                value_node = cell.find(f"{{{XLSX_MAIN_NS}}}v")
                if cell_type == "inlineStr":
                    value = "".join(
                        node.text or ""
                        for node in cell.iter(f"{{{XLSX_MAIN_NS}}}t")
                    )
                elif value_node is None:
                    value = None
                elif cell_type == "s":
                    value = shared_strings[int(value_node.text or "0")]
                else:
                    value = value_node.text
                values[index] = value
            rows.append(values)

    if not rows:
        raise ValueError(f"No rows found in {sheet_name!r} of {path}")
    max_column = max(max(row.keys(), default=-1) for row in rows)
    header = [rows[0].get(index) for index in range(max_column + 1)]
    if any(value is None for value in header):
        raise ValueError(f"Blank header cells in {sheet_name!r} of {path}")
    data = [
        [row.get(index) for index in range(max_column + 1)]
        for row in rows[1:]
    ]
    return pd.DataFrame(data, columns=[str(value) for value in header])


def prepare_original_composite_meta(path: Path) -> pd.DataFrame:
    meta = read_xlsx_table(path, "Master_GNPS_Mapping")
    consensus_columns = [
        column
        for column in pd.read_csv(recovered.CONSENSUS, nrows=0).columns
        if column.startswith("sample:")
    ]
    lookup = {column.lower(): column for column in consensus_columns}
    key = (
        meta["Consensus Column"]
        .fillna("")
        .astype(str)
        .str.replace(":area", "", regex=False)
        .str.lower()
    )
    meta["original_header"] = key.map(lookup)
    meta = meta[
        meta["original_header"].notna() & meta["Phylum"].notna()
    ].copy()
    meta = meta.drop_duplicates("original_header", keep="last")
    return meta.rename(
        columns={
            "Phylum": "phylum",
            "Kingdom": "kingdom",
            "Batch Group": "batch",
        }
    )[["original_header", "phylum", "kingdom", "batch"]]


def configure(source_root: Path, output_dir: Path) -> dict[str, Path]:
    paths = {
        "consensus": source_root / "analysis/analysis-15/03_alignment/consensus_aligned_table.csv",
        "old_metadata": source_root / "analysis/analysis-15/04_biomarker_discovery/01_composite_scoring/sample_metadata.csv",
        "old_composite_workbook": source_root / "analysis/analysis-13/repository_bundle/01_input_data/GNPS_Master_Sample_Mapping_Positive.xlsx",
        "old_platinum_raw": source_root / "analysis/analysis-15/04_biomarker_discovery/04_platinum_diamond/phylum_biomarkers_platinum.csv",
        "old_silver_raw": source_root / "analysis/analysis-15/04_biomarker_discovery/04_platinum_diamond/phylum_biomarkers_silver.csv",
        "old_platinum_pre_archlips": source_root / "analysis/analysis-15/04_biomarker_discovery/04_platinum_diamond/atlas_platinum_biomarkers_final_backup_pre_fullbatch.csv",
        "old_platinum_final": source_root / "analysis/analysis-15/04_biomarker_discovery/04_platinum_diamond/atlas_platinum_biomarkers_final.csv",
        "old_indval_pairs": source_root / "analysis/analysis-17/positive/03_within_batch/multi_batch_indicators.csv",
        "old_atlas": source_root / "analysis/analysis-15/04_biomarker_discovery/04_platinum_diamond/atlas_expanded_final.csv",
        "harmonized_annotations": source_root / "analysis/analysis-15/04_biomarker_discovery/10_annotation_harmonization/harmonized_annotations.csv",
        "ms2_cache": source_root / "analysis/analysis-15/04_biomarker_discovery/06_ms2_structural/ms2_diagnostic_classification.csv",
    }
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing Figure 2A inputs: " + ", ".join(missing))

    recovered.SOURCE_ROOT = source_root
    recovered.CONSENSUS = paths["consensus"]
    recovered.META_OLD = paths["old_metadata"]
    recovered.META_ORIGINAL_COMPOSITE = paths["old_composite_workbook"]
    recovered.OLD_PLAT_RAW = paths["old_platinum_raw"]
    recovered.OLD_SILVER_RAW = paths["old_silver_raw"]
    recovered.OLD_PLAT_PRE_ARCH = paths["old_platinum_pre_archlips"]
    recovered.OLD_PLAT_FINAL = paths["old_platinum_final"]
    recovered.OLD_INDVAL = paths["old_indval_pairs"]
    recovered.OLD_EXPANDED = paths["old_atlas"]
    recovered.HARMONIZED = paths["harmonized_annotations"]
    recovered.MS2_CLASS = paths["ms2_cache"]
    recovered.OUT = output_dir
    recovered.KINGDOM_ORDER = [
        "Bacteria", "Archaea", "Fungi", "Viridiplantae", "Animalia", "Protists"
    ]
    recovered.KINGDOM_COLOURS = {
        "Bacteria": "#0072B2",
        "Archaea": "#E69F00",
        "Fungi": "#009E73",
        "Viridiplantae": "#56B4E9",
        "Animalia": "#D55E00",
        "Protists": "#CC79A7",
    }
    return paths


def load_strict_metadata(taxonomy_dir: Path) -> tuple[pd.DataFrame, dict]:
    summary_path = taxonomy_dir / "taxonomy_summary.json"
    metadata_path = taxonomy_dir / "sample_metadata_POS_ncbi_phylum.csv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("assertions", {}).get("status") != "pass":
        raise ValueError("Taxonomy freeze assertions did not pass")
    analysis_phyla = set(summary["analysis_phyla"])
    metadata = pd.read_csv(metadata_path, low_memory=False)
    metadata = metadata[
        metadata["taxonomy_scope"].eq("core_candidate")
        & metadata["phylum"].isin(analysis_phyla)
        & metadata["original_header"].notna()
    ].copy()
    metadata["kingdom"] = metadata["ecological_group"]
    metadata["phylum"] = metadata["phylum"].astype(str).str.strip()
    metadata["batch"] = metadata["batch"].astype(str).str.strip()
    if metadata["original_header"].duplicated().any():
        duplicates = metadata.loc[
            metadata["original_header"].duplicated(keep=False), "original_header"
        ].tolist()
        raise ValueError(f"Duplicated strict POS sample IDs: {duplicates[:10]}")
    return metadata, summary


def pair_overlap(expected: pd.DataFrame, observed: pd.DataFrame) -> dict:
    keys = ["feature_id", "phylum"]
    expected_pairs = expected[keys].drop_duplicates()
    observed_pairs = observed[keys].drop_duplicates()
    overlap = expected_pairs.merge(observed_pairs, on=keys)
    report = {
        "expected": int(len(expected_pairs)),
        "observed": int(len(observed_pairs)),
        "overlap": int(len(overlap)),
        "expected_coverage_pct": round(100 * len(overlap) / max(len(expected_pairs), 1), 4),
        "observed_precision_pct": round(100 * len(overlap) / max(len(observed_pairs), 1), 4),
    }
    report["exact"] = (
        report["expected"] == report["observed"] == report["overlap"]
    )
    return report


def reproduction_audit(
    old_atlas: pd.DataFrame,
    observed_platinum: pd.DataFrame,
    observed_silver: pd.DataFrame,
    observed_indval: pd.DataFrame,
) -> dict:
    expected_platinum = old_atlas[
        old_atlas["biomarker_quality"].eq("Platinum")
        & old_atlas["composite_score"].notna()
    ]
    expected_silver = old_atlas[old_atlas["biomarker_quality"].eq("Silver")]
    expected_indval = old_atlas[old_atlas["biomarker_quality"].eq("IndVal")]
    report = {
        "source": "reconstructed from packaged atlas_expanded_final.csv",
        "platinum": pair_overlap(expected_platinum, observed_platinum),
        "silver": pair_overlap(expected_silver, observed_silver),
        "indval_resolved": pair_overlap(expected_indval, observed_indval),
        "limitation": (
            "The original Composite mapping workbook and unresolved 2,906-row IndVal pair "
            "table are absent from the portable package; the packaged 2,698 resolved IndVal "
            "features are used for the available gate."
        ),
    }
    report["exact_available_gate"] = all(
        report[name]["exact"] for name in ("platinum", "silver", "indval_resolved")
    )
    return report


def annotation_catalog(old_atlas: pd.DataFrame, ms2_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    previous = old_atlas.sort_values("feature_id").drop_duplicates("feature_id").copy()
    previous["previous_annotation_tier"] = (
        previous["annot_confidence"].fillna("Unidentified").astype(str).str.title()
    )
    previous_catalog = previous[[
        "feature_id", "previous_annotation_tier", "annotation_source",
        "clean_annotation_source", "has_ms2", "has_formula", "has_canopus",
        "has_denovo", "archlips_source",
    ]].copy()

    ms2 = pd.read_csv(ms2_path, low_memory=False)
    class_field = (
        "ms2_assigned_class" if "ms2_assigned_class" in ms2.columns else "assigned_class"
    )
    ms2["ms2_cache_assigned"] = (
        ms2[class_field].notna()
        & ~ms2[class_field].astype(str).str.casefold().isin({"unknown", "nan", ""})
    )
    ms2_catalog = (
        ms2.sort_values("feature_id")
        .drop_duplicates("feature_id")[["feature_id", class_field, "ms2_cache_assigned"]]
        .rename(columns={class_field: "ms2_cached_class"})
    )

    tiers = previous_catalog[["feature_id", "previous_annotation_tier"]].rename(
        columns={"previous_annotation_tier": "annotation_tier"}
    )
    new_ms2 = ms2_catalog[
        ms2_catalog["ms2_cache_assigned"]
        & ~ms2_catalog["feature_id"].isin(set(tiers["feature_id"]))
    ][["feature_id"]].copy()
    new_ms2["annotation_tier"] = "Bronze"
    tiers = pd.concat([tiers, new_ms2], ignore_index=True)
    return tiers, previous_catalog.merge(ms2_catalog, on="feature_id", how="outer")


def annotation_audit(
    expanded: pd.DataFrame,
    catalog: pd.DataFrame,
    phyla: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    fields = [
        "feature_id", "phylum", "kingdom", "discovery_method", "biomarker_quality",
        "biomarker_tier", "annotation_tier",
    ]
    audit = expanded[fields].copy().merge(catalog, on="feature_id", how="left")
    audit["in_previous_full_annotation_pipeline"] = audit["previous_annotation_tier"].notna()
    audit["previously_annotated"] = (
        audit["in_previous_full_annotation_pipeline"]
        & ~audit["previous_annotation_tier"].eq("Unidentified")
    )
    audit["previously_processed_unidentified"] = (
        audit["in_previous_full_annotation_pipeline"]
        & audit["previous_annotation_tier"].eq("Unidentified")
    )
    audit["new_biomarker_needing_full_pipeline"] = ~audit[
        "in_previous_full_annotation_pipeline"
    ]
    audit["new_with_assigned_ms2_cache"] = (
        audit["new_biomarker_needing_full_pipeline"]
        & audit["ms2_cache_assigned"].eq(True)
    )
    audit["new_without_current_annotation"] = (
        audit["new_biomarker_needing_full_pipeline"]
        & audit["annotation_tier"].eq("Unidentified")
    )
    audit["currently_annotated"] = ~audit["annotation_tier"].eq("Unidentified")

    summary_rows = []
    for phylum in phyla:
        part = audit[audit["phylum"].eq(phylum)]
        summary_rows.append({
            "phylum": phylum,
            "total_biomarkers": int(len(part)),
            "previous_full_pipeline": int(part["in_previous_full_annotation_pipeline"].sum()),
            "previously_annotated": int(part["previously_annotated"].sum()),
            "previously_processed_unidentified": int(part["previously_processed_unidentified"].sum()),
            "new_need_full_annotation_pipeline": int(part["new_biomarker_needing_full_pipeline"].sum()),
            "new_with_assigned_ms2_cache": int(part["new_with_assigned_ms2_cache"].sum()),
            "new_without_current_annotation": int(part["new_without_current_annotation"].sum()),
            "currently_annotated_after_cache_reuse": int(part["currently_annotated"].sum()),
            "currently_unidentified": int(part["annotation_tier"].eq("Unidentified").sum()),
        })
    return audit, pd.DataFrame(summary_rows)


def add_zero_rows(counts: pd.DataFrame, phyla: list[str], policy: dict) -> pd.DataFrame:
    missing = [phylum for phylum in phyla if phylum not in set(counts["phylum"])]
    if not missing:
        return counts
    rows = [{
        "phylum": phylum,
        "kingdom": policy["ecological_group"].get(phylum, ""),
        "annotation_tier": "Unidentified",
        "n_features": 0,
    } for phylum in missing]
    return pd.concat([counts, pd.DataFrame(rows)], ignore_index=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--taxonomy-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    started = time.time()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    paths = configure(args.source_root.resolve(), args.output_dir.resolve())
    strict_meta, taxonomy = load_strict_metadata(args.taxonomy_dir.resolve())
    policy = json.loads(
        (Path(__file__).resolve().parent / "config/taxonomy_policy.json").read_text(
            encoding="utf-8"
        )
    )
    phyla = taxonomy["analysis_phyla"]
    old_meta = recovered.prepare_meta(paths["old_metadata"])
    old_composite_meta = prepare_original_composite_meta(
        paths["old_composite_workbook"]
    )
    old_atlas = pd.read_csv(paths["old_atlas"], low_memory=False)

    historical_load_meta = pd.concat(
        [old_meta, old_composite_meta], ignore_index=True, sort=False
    )
    features, samples, matrix = recovered.load_consensus(
        historical_load_meta, strict_meta
    )

    old_platinum, old_silver, _ = recovered.composite_select(
        features, samples, matrix, old_composite_meta, "archived Composite metadata"
    )
    old_pairs, old_indval = recovered.run_indval(
        features, samples, matrix, old_meta, "packaged old labels"
    )
    reproduction = recovered.assert_old_reproduction(
        old_composite_meta, old_platinum, old_silver, old_pairs, old_indval
    )
    historical_dir = args.output_dir / "historical_control"
    historical_dir.mkdir(parents=True, exist_ok=True)
    old_platinum.to_csv(
        historical_dir / "composite_platinum_reproduced.csv", index=False
    )
    old_silver.to_csv(
        historical_dir / "composite_silver_reproduced.csv", index=False
    )
    old_pairs.to_csv(historical_dir / "indval_pairs_reproduced.csv", index=False)
    old_indval.to_csv(
        historical_dir / "indval_unique_reproduced.csv", index=False
    )
    (historical_dir / "reproduction_gate.json").write_text(
        json.dumps(reproduction, indent=2) + "\n", encoding="utf-8"
    )
    if not reproduction["pass"]:
        raise RuntimeError(
            "Historical Composite/IndVal reproduction gate failed; strict rerun aborted"
        )

    new_platinum, new_silver, new_base = recovered.composite_select(
        features, samples, matrix, strict_meta, "strict 16-phylum labels"
    )
    new_pairs, new_indval = recovered.run_indval(
        features, samples, matrix, strict_meta, "strict 16-phylum labels"
    )

    pre_arch = pd.read_csv(paths["old_platinum_pre_archlips"], usecols=["feature_id"])
    final_arch_source = pd.read_csv(paths["old_platinum_final"], low_memory=False)
    arch_source = final_arch_source[
        ~final_arch_source["feature_id"].isin(set(pre_arch["feature_id"]))
    ].copy()
    if len(arch_source) != 1095:
        raise ValueError(f"Expected 1,095 packaged ArchLips additions, found {len(arch_source)}")
    old_arch = recovered.assign_archlips(
        features, samples, matrix, old_meta, arch_source["feature_id"].tolist()
    )
    arch_accuracy = float(
        old_arch.merge(
            arch_source[["feature_id", "phylum"]], on="feature_id", suffixes=("_new", "_old")
        ).eval("phylum_new == phylum_old").mean()
    )
    new_arch = recovered.assign_archlips(
        features, samples, matrix, strict_meta, arch_source["feature_id"].tolist()
    )

    tiers, catalog = annotation_catalog(old_atlas, paths["ms2_cache"])
    composite, expanded = recovered.build_expanded(
        features, new_base, new_arch, arch_source, new_indval, tiers
    )
    expanded = expanded[expanded["phylum"].isin(phyla)].copy()
    if expanded["feature_id"].duplicated().any():
        raise ValueError("Final strict atlas contains duplicated feature IDs")

    counts = recovered.tier_counts(expanded)
    counts = add_zero_rows(counts, phyla, policy)
    audit, coverage = annotation_audit(expanded, catalog, phyla)

    discovery_summary = (
        expanded.groupby(["phylum", "discovery_method"], as_index=False)
        .size()
        .rename(columns={"size": "n_biomarkers"})
        .sort_values(["phylum", "discovery_method"])
    )
    queue_fields = [
        "feature_id", "consensus_mz", "consensus_rt", "phylum", "kingdom",
        "discovery_method", "biomarker_quality", "biomarker_tier", "annotation_tier",
    ]
    annotation_queue = expanded[queue_fields].merge(
        audit[[
            "feature_id", "ms2_cached_class", "ms2_cache_assigned",
            "new_with_assigned_ms2_cache", "new_without_current_annotation",
        ]],
        on="feature_id", how="left", validate="one_to_one",
    )
    annotation_queue = annotation_queue[
        ~annotation_queue["feature_id"].isin(
            set(audit.loc[
                audit["in_previous_full_annotation_pipeline"], "feature_id"
            ])
        )
    ].sort_values(["phylum", "feature_id"])

    sample_counts = strict_meta["phylum"].value_counts().rename("pos_samples")
    tier_wide = (
        counts.pivot_table(
            index="phylum", columns="annotation_tier", values="n_features", fill_value=0
        )
        .reindex(index=phyla, columns=recovered.TIER_ORDER, fill_value=0)
        .reset_index()
    )
    biomarker_summary = tier_wide.merge(
        sample_counts.rename_axis("phylum").reset_index(), on="phylum", how="left"
    )
    biomarker_summary["ecological_group"] = biomarker_summary["phylum"].map(
        policy["ecological_group"]
    )
    biomarker_summary["total_biomarkers"] = biomarker_summary[
        recovered.TIER_ORDER
    ].sum(axis=1)
    coverage = coverage.merge(
        biomarker_summary[["phylum", "pos_samples", "ecological_group"]],
        on="phylum", how="left",
    )

    expanded.to_csv(args.output_dir / "atlas_pos_strict.csv", index=False)
    new_platinum.to_csv(args.output_dir / "composite_platinum_strict.csv", index=False)
    new_silver.to_csv(args.output_dir / "composite_silver_strict.csv", index=False)
    new_pairs.to_csv(args.output_dir / "indval_pairs_strict.csv", index=False)
    new_indval.to_csv(args.output_dir / "indval_unique_strict.csv", index=False)
    new_arch.to_csv(args.output_dir / "archlips_assignments_strict.csv", index=False)
    counts.to_csv(args.output_dir / "figure2a_tier_counts.csv", index=False)
    biomarker_summary.to_csv(args.output_dir / "biomarker_counts_by_phylum.csv", index=False)
    audit.to_csv(args.output_dir / "annotation_coverage_by_feature.csv", index=False)
    coverage.to_csv(args.output_dir / "annotation_coverage_by_phylum.csv", index=False)
    discovery_summary.to_csv(
        args.output_dir / "discovery_method_by_phylum.csv", index=False
    )
    annotation_queue.to_csv(
        args.output_dir / "annotation_queue_new_features.csv", index=False
    )
    recovered.render_figure(counts, args.output_dir / "Figure_2a_strict_16phyla")

    totals = {
        "biomarkers": int(len(expanded)),
        "phyla": int(len(phyla)),
        "previous_full_pipeline": int(audit["in_previous_full_annotation_pipeline"].sum()),
        "previously_annotated": int(audit["previously_annotated"].sum()),
        "previously_processed_unidentified": int(audit["previously_processed_unidentified"].sum()),
        "new_need_full_annotation_pipeline": int(audit["new_biomarker_needing_full_pipeline"].sum()),
        "new_with_assigned_ms2_cache": int(audit["new_with_assigned_ms2_cache"].sum()),
        "new_without_current_annotation": int(audit["new_without_current_annotation"].sum()),
        "currently_annotated_after_cache_reuse": int(audit["currently_annotated"].sum()),
        "currently_unidentified": int(audit["annotation_tier"].eq("Unidentified").sum()),
    }
    release_status = (
        "candidate"
        if reproduction["pass"] and arch_accuracy == 1.0
        else "diagnostic_reproduction_gap"
    )
    validation = {
        "atlas_rows": int(len(expanded)),
        "distinct_feature_ids": int(expanded["feature_id"].nunique()),
        "duplicate_feature_ids": int(expanded["feature_id"].duplicated().sum()),
        "observed_phyla": int(expanded["phylum"].nunique()),
        "expected_phyla": int(len(phyla)),
        "unexpected_phyla": sorted(set(expanded["phylum"]) - set(phyla)),
        "tier_total": int(counts["n_features"].sum()),
        "annotation_partition_total": int(
            totals["previous_full_pipeline"] + totals["new_need_full_annotation_pipeline"]
        ),
        "annotation_queue_rows": int(len(annotation_queue)),
    }
    if not (
        validation["atlas_rows"] == validation["distinct_feature_ids"]
        == validation["tier_total"] == validation["annotation_partition_total"]
        and validation["duplicate_feature_ids"] == 0
        and validation["observed_phyla"] == validation["expected_phyla"]
        and not validation["unexpected_phyla"]
        and validation["annotation_queue_rows"] == totals["new_need_full_annotation_pipeline"]
    ):
        raise ValueError(f"Strict Figure 2A validation failed: {validation}")
    summary = {
        "taxonomy_release": taxonomy["taxonomy_release"],
        "scope": "Figure 2A POS, strict shared 16-phylum set",
        "status": release_status,
        "reproduction": reproduction,
        "archlips_assignment_reproduction_accuracy": arch_accuracy,
        "strict_selection": {
            "composite_platinum": int(len(new_platinum)),
            "composite_silver": int(len(new_silver)),
            "composite_total": int(len(new_base)),
            "indval_pair_rows": int(len(new_pairs)),
            "indval_unique": int(len(new_indval)),
            "archlips_candidates": int(len(arch_source)),
            "composite_with_archlips": int(len(composite)),
        },
        "annotation_coverage": totals,
        "validation": validation,
        "remaining_external_producers": [
            "fresh SIRIUS/CANOPUS, DreaMS, fastMASST and Pan-ReDU results for newly selected feature IDs",
            "exact historical 66-entry/52.8-percent classifier producer",
        ],
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
        },
        "inputs": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in paths.items()
        },
        "seconds": round(time.time() - started, 3),
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    output_records = []
    for path in sorted(args.output_dir.iterdir()):
        if path.is_file() and path.name != "output_manifest.json":
            output_records.append({
                "path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)
            })
    (args.output_dir / "output_manifest.json").write_text(
        json.dumps({
            "taxonomy_release": taxonomy["taxonomy_release"],
            "status": release_status,
            "outputs": output_records,
        }, indent=2) + "\n",
        encoding="utf-8",
    )

    print(json.dumps({
        "status": release_status,
        "historical_reproduction_gate": reproduction["pass"],
        "archlips_reproduction_accuracy": arch_accuracy,
        "annotation_coverage": totals,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
