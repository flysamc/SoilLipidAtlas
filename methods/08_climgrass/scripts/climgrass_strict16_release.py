#!/usr/bin/env python3
"""Run review-only ClimGrass decomposition for the locked strict 16-phylum release.

This producer reuses the corrected taxonomy-independent soil-to-atlas spectral
matches, remaps them at 5 ppm to the current release SIMPER atlas, prepares an
exact feature subset from the full POS consensus table, and reruns the PC-primary,
PE-sensitivity, and PC+PE-sensitivity decomposition arms. The primary ArchLips
mode fails closed when no release-eligible archaeal feature reaches the mapped
soil substrate; an explicitly named all-SIMPER sensitivity mode is also available.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RELEASE_ID = "ncbi-phylum-2026-08-04-v1"
HANDOFF_ROOT = (
    PROJECT_ROOT / "cg" / "SHADOW_PC_CLIMGRASS_STRICT_PHYLA_HANDOFF_2026-08-06"
)
RECOVERY_ROOT = (
    PROJECT_ROOT / "cg" / "POS_CLIMGRASS_DECOMPOSITION_RECOVERY-2026-08-06"
)
LEGACY_RUNNER = (
    HANDOFF_ROOT / "source" / "strict_taxonomy" / "run_strict_19unit_analysis.py"
)

PC_NAME = "PC_15-18d7"
PE_NAME = "PE_15-18d7_precursor"
COMBINED_NAME = "PC_PE_two_standard"
ARMS = [
    ("strict16_PC_primary", PC_NAME, "PC-d7 primary"),
    ("strict16_PE_sensitivity", PE_NAME, "PE-d7 precursor sensitivity"),
    ("strict16_PC_PE_sensitivity", COMBINED_NAME, "PC+PE sensitivity"),
]
DETECTION_THRESHOLD = 0.005


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--spectral-matches",
        type=Path,
        default=(
            HANDOFF_ROOT
            / "results"
            / "strict_19unit_corrected_2026-08-06"
            / "corrected_spectral_matches.csv"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            PROJECT_ROOT
            / "outputs"
            / "analysis"
            / RELEASE_ID
            / "climgrass"
            / "strict16_soil_mgf_2026-08-06"
        ),
    )
    parser.add_argument(
        "--archaea-mode",
        choices=("archlips", "all-simper-sensitivity"),
        default="archlips",
        help="Primary ArchLips restriction or explicitly non-primary all-SIMPER sensitivity.",
    )
    return parser.parse_args()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_legacy_runner():
    if not LEGACY_RUNNER.is_file():
        raise FileNotFoundError(LEGACY_RUNNER)
    spec = importlib.util.spec_from_file_location("climgrass_recovered_runner", LEGACY_RUNNER)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load recovered runner: {LEGACY_RUNNER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.ARMS = ARMS
    return module


def required_sources(spectral_matches: Path) -> dict[str, Path]:
    release_root = PROJECT_ROOT / "outputs" / "analysis" / RELEASE_ID
    return {
        "taxonomy_policy": PROJECT_ROOT / "paper2_repro" / "config" / "taxonomy_policy.json",
        "taxonomy_summary": release_root / "taxonomy" / "taxonomy_summary.json",
        "metadata_pos": release_root / "taxonomy" / "sample_metadata_POS_ncbi_phylum.csv",
        "simper_atlas": release_root / "simper" / "simper_atlas.csv",
        "spectral_matches": spectral_matches.resolve(),
        "consensus_full": PROJECT_ROOT / "analysis" / "analysis-15" / "03_alignment" / "consensus_aligned_table.csv",
        "annotation_full": (
            RECOVERY_ROOT
            / "analysis"
            / "analysis-15"
            / "04_biomarker_discovery"
            / "02_lipidsearch_annotations"
            / "consensus_unified_annotations.csv"
        ),
        "climgrass_quant": HANDOFF_ROOT / "data_contract" / "soil-data-without-background_quant.csv",
        "climgrass_meta": HANDOFF_ROOT / "data_contract" / "climgrass_sample_metadata.csv",
        "rie_table": HANDOFF_ROOT / "data_contract" / "rie_table_s10.csv",
        "expected_ref": HANDOFF_ROOT / "data_contract" / "expected_kingdom_composition.csv",
        "archlips_pos_rt_screened": (
            release_root
            / "annotation"
            / "step8_archlips_rt_filtered"
            / "archlips_pos_rt_screened.csv"
        ),
    }


def map_spectral_matches(atlas: pd.DataFrame, matches: pd.DataFrame) -> pd.DataFrame:
    unique = atlas.drop_duplicates(subset="feature_id", keep="first")
    mzs = unique["consensus_mz"].to_numpy(dtype=float)
    feature_ids = unique["feature_id"].astype(str).to_numpy()
    order = np.argsort(mzs)
    mzs = mzs[order]
    feature_ids = feature_ids[order]
    mapped = []
    for row in matches.itertuples(index=False):
        atlas_mz = float(row.atlas_mz)
        tolerance = atlas_mz * 5.0 / 1e6
        lower = np.searchsorted(mzs, atlas_mz - tolerance, side="left")
        upper = np.searchsorted(mzs, atlas_mz + tolerance, side="right")
        if lower >= upper:
            continue
        candidates = np.arange(lower, upper)
        ppm = np.abs(mzs[candidates] - atlas_mz) / atlas_mz * 1e6
        best = candidates[np.argmin(ppm)]
        mapped.append(
            {
                "soil_scan": row.soil_scan,
                "atlas_mz": atlas_mz,
                "cosine": float(row.cosine),
                "feature_id": feature_ids[best],
                "ppm_diff": float(ppm[np.argmin(ppm)]),
            }
        )
    result = pd.DataFrame(mapped)
    if result.empty:
        raise ValueError("No corrected spectral matches map to the strict SIMPER atlas")
    return (
        result.sort_values("cosine", ascending=False)
        .drop_duplicates(subset="feature_id", keep="first")
        .reset_index(drop=True)
    )


def write_feature_subset(source: Path, destination: Path, feature_ids: set[str]) -> dict:
    chunks = []
    for chunk in pd.read_csv(source, chunksize=5000, low_memory=False):
        selected = chunk[chunk["feature_id"].astype(str).isin(feature_ids)]
        if not selected.empty:
            chunks.append(selected)
    if not chunks:
        raise ValueError(f"No requested features found in {source}")
    subset = pd.concat(chunks, ignore_index=True)
    observed = set(subset["feature_id"].astype(str))
    missing = sorted(feature_ids - observed)
    duplicated = int(subset["feature_id"].astype(str).duplicated().sum())
    if missing or duplicated:
        raise ValueError(
            f"Invalid feature subset from {source}: missing={missing[:10]}, duplicates={duplicated}"
        )
    subset.to_csv(destination, index=False)
    return {
        "rows": int(len(subset)),
        "columns": int(len(subset.columns)),
        "feature_ids_exact": observed == feature_ids,
        "sha256": sha256_file(destination),
    }


def configure_units(legacy, cache: dict, sources: dict[str, Path]) -> tuple[list[str], dict[str, str], dict]:
    policy = json.loads(sources["taxonomy_policy"].read_text(encoding="utf-8"))
    summary = json.loads(sources["taxonomy_summary"].read_text(encoding="utf-8"))
    if policy["release_id"] != RELEASE_ID or policy["status"] != "locked":
        raise ValueError("Taxonomy policy is not the expected locked release")
    strict_units = list(summary["analysis_phyla"])
    if len(strict_units) != policy["expected_after_freeze"]["analysis_phyla"]:
        raise ValueError("Strict analysis-phylum count disagrees with the locked policy")
    atlas_units = sorted(cache["simper"]["phylum"].astype(str).unique())
    if sorted(strict_units) != atlas_units:
        raise ValueError("SIMPER atlas phyla do not exactly match taxonomy_summary analysis_phyla")

    group_map = {}
    for unit in strict_units:
        group = policy["ecological_group"][unit]
        group_map[unit] = "Plantae" if group == "Viridiplantae" else group
    expected_groups = {"Animalia", "Archaea", "Bacteria", "Fungi", "Plantae", "Protists"}
    if set(group_map.values()) != expected_groups:
        raise ValueError(f"Unexpected policy display groups: {sorted(set(group_map.values()))}")
    # The recovered framework calls the protist display stratum Protozoa.
    framework_group_map = {
        unit: ("Protozoa" if group == "Protists" else group)
        for unit, group in group_map.items()
    }
    legacy.decomposition.PHYLUM_KINGDOM.clear()
    legacy.decomposition.PHYLUM_KINGDOM.update(framework_group_map)

    cache["sample_phylum"] = {
        sample: phylum
        for sample, phylum in cache["sample_phylum"].items()
        if phylum in set(strict_units)
    }
    counts = Counter(cache["sample_phylum"].values())
    expected_counts = {unit: int(summary["positive_counts"][unit]) for unit in strict_units}
    if dict(counts) != expected_counts:
        raise ValueError(f"Strict metadata counts disagree with taxonomy summary: {counts}")
    if min(counts.values()) < policy["analysis_rule"]["minimum_samples_per_polarity"]:
        raise ValueError("A strict analysis phylum is below the locked POS replication threshold")
    diagnostics = {
        "taxonomy_release": RELEASE_ID,
        "n_phyla": len(strict_units),
        "analysis_phyla": strict_units,
        "positive_analysis_samples": int(sum(counts.values())),
        "positive_collection_samples": int(summary["positive_core_samples"]),
        "excluded_below_threshold": summary["below_threshold"],
        "unit_sample_counts": dict(sorted(counts.items())),
        "display_group_normalization": {
            "Viridiplantae": "Plantae",
            "Protists": "Protozoa",
            "note": "display-only compatibility with the recovered six-group evaluation; analysis units remain NCBI phyla",
        },
    }
    return strict_units, framework_group_map, diagnostics


def validate_compositions(arm_dir: Path, label: str, strict_units: list[str]) -> dict:
    expected_groups = {"Animalia", "Archaea", "Bacteria", "Fungi", "Plantae", "Protozoa"}
    checks = {}
    for method in ("nnls", "standard_bc", "enriched_only_bc", "fc_weighted_bc"):
        phylum = pd.read_csv(arm_dir / f"phylum_composition_{label}_{method}.csv", index_col=0)
        group = pd.read_csv(arm_dir / f"kingdom_composition_{label}_{method}.csv", index_col=0)
        if phylum.shape != (12, len(strict_units)) or set(phylum.columns) != set(strict_units):
            raise ValueError(f"Invalid strict-16 composition: {label}/{method}")
        if group.shape != (12, 6) or set(group.columns) != expected_groups:
            raise ValueError(f"Invalid six-group display composition: {label}/{method}")
        for name, frame in (("phylum", phylum), ("display_group", group)):
            values = frame.to_numpy(dtype=float)
            if not np.isfinite(values).all() or (values < -1e-12).any():
                raise ValueError(f"Invalid values: {label}/{method}/{name}")
            error = float(np.max(np.abs(values.sum(axis=1) - 1.0)))
            if error > 1e-9:
                raise ValueError(f"Rows do not sum to one: {label}/{method}/{name}")
            checks[f"{method}_{name}_max_sum_error"] = error
    return checks


def plot_effects(effects: pd.DataFrame, destination: Path) -> None:
    colors = {
        "Bacteria": "#0072B2",
        "Archaea": "#CC79A7",
        "Fungi": "#D55E00",
        "Plantae": "#009E73",
        "Animalia": "#E69F00",
        "Protozoa": "#56B4E9",
    }
    figure, axes = plt.subplots(1, 3, figsize=(15, 5.8), sharex=True, sharey=True)
    for axis, (arm, group) in zip(axes, effects.groupby("arm", sort=False)):
        axis.axhline(0, color="#cccccc", linewidth=0.8)
        axis.axvline(0, color="#cccccc", linewidth=0.8)
        for row in group.itertuples(index=False):
            significant = row.drought_fdr_0_05 or row.climate_fdr_0_05
            axis.scatter(
                row.drought_log2FC,
                row.climate_log2FC,
                s=35 + 1200 * row.mean_fraction,
                color=colors[row.kingdom],
                edgecolor="black",
                linewidth=1.8 if significant else 0.4,
            )
            if significant:
                axis.annotate(row.unit, (row.drought_log2FC, row.climate_log2FC), fontsize=7)
        axis.set_title(arm)
        axis.spines[["top", "right"]].set_visible(False)
    figure.supxlabel("Drought effect: log2(drought / no drought)", y=0.18)
    figure.supylabel("Climate effect: log2(future / ambient)")
    figure.suptitle("Corrected ClimGrass effects: locked strict 16-phylum release", y=0.98)
    legend = [
        Line2D(
            [0], [0], marker="o", linestyle="", markersize=7,
            markerfacecolor=color, markeredgecolor="black", markeredgewidth=0.4,
            label=group,
        )
        for group, color in colors.items()
    ]
    figure.legend(
        handles=legend, loc="lower center", ncol=6, frameon=False,
        bbox_to_anchor=(0.5, 0.005), title="Display group",
    )
    n_significant = int(effects["drought_fdr_0_05"].sum() + effects["climate_fdr_0_05"].sum())
    figure.text(
        0.5, 0.115,
        f"FC-weighted composition; Mann-Whitney U, n=6 vs 6; BH-FDR significant tests: {n_significant}",
        ha="center", fontsize=9, color="#444444",
    )
    figure.tight_layout(rect=(0.02, 0.24, 1, 0.91))
    figure.savefig(destination.with_suffix(".png"), dpi=300, facecolor="white")
    figure.savefig(destination.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)


def write_readme(
    output: Path,
    spectral_qc: dict,
    mapping_contract: dict,
    taxonomy: dict,
    significance: pd.DataFrame,
    archaea_mode: str,
) -> None:
    significant = significance[
        (significance["drought_q_lt_0_05"] > 0)
        | (significance["climate_q_lt_0_05"] > 0)
    ]
    lines = [
        "# ClimGrass decomposition: locked strict 16-phylum release",
        "",
        f"Taxonomy release: `{RELEASE_ID}` (locked). Status: **review-only complete**.",
        "",
        "The corrected taxonomy-independent spectral matches were reused and remapped to the current strict SIMPER atlas at exactly 5 ppm.",
        f"Archaeal reference mode: `{archaea_mode}`.",
        "",
        "## Verified run facts",
        "",
        f"- Corrected spectral matches: {spectral_qc['rows']:,}",
        f"- Corrected cosine maximum: {spectral_qc['cosine_max']:.6f}",
        f"- Strict mapped features: {mapping_contract['unique_simper_features']:,}",
        f"- Strict analysis phyla: {taxonomy['n_phyla']}",
        f"- Strict POS reference samples: {taxonomy['positive_analysis_samples']}",
        "- PC-d7 and PE-d7 precursor soil detection: 12/12 each",
        "- Treatment design: balanced 2x2, n=3 per cell",
        "",
        "## FDR result",
        "",
        (
            "No PC-primary or sensitivity arm has an FDR-significant drought or climate effect."
            if significant.empty
            else "At least one arm has an FDR-significant effect; inspect `phylum_effects_all_arms.csv`."
        ),
        "",
        "## Boundary",
        "",
        "These outputs replace the diagnostic 19-unit downstream run only as a strict-16 review artifact. They are not connected to manuscript, table, or Figure 5 consumers and are not a final paper release.",
        "",
        (
            "The primary ArchLips-restricted contract passed."
            if archaea_mode == "archlips"
            else "This is a non-primary all-SIMPER archaeal-reference sensitivity. The primary ArchLips-restricted run is blocked because zero release-eligible ArchLips features intersect the mapped strict-16 soil substrate."
        ),
        "",
        "PE-only and PC+PE remain sensitivity analyses because the PE DDA spectrum is nondiagnostic. PC-d7 is primary.",
        "",
    ]
    (output / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    args = parse_args()
    output = args.output.resolve()
    staging = output.with_name(output.name + ".incomplete")
    if output.exists() or staging.exists():
        raise FileExistsError(f"Refusing to overwrite existing output or staging directory: {output}")

    sources = required_sources(args.spectral_matches.resolve())
    missing = [path for path in sources.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required input(s):\n" + "\n".join(map(str, missing)))
    legacy = load_legacy_runner()
    matches, spectral_qc = legacy.validate_spectral_matches(sources["spectral_matches"])
    atlas = pd.read_csv(sources["simper_atlas"], low_memory=False)
    seed_mapping = map_spectral_matches(atlas, matches)
    if float(seed_mapping["ppm_diff"].max()) > 5.0 + 1e-12:
        raise ValueError("Strict-16 mapping exceeded 5 ppm")

    staging.mkdir(parents=True)
    input_dir = staging / "prepared_inputs"
    input_dir.mkdir()
    seed_mapping.to_csv(input_dir / "strict16_seed_mapping_5ppm.csv", index=False)
    feature_ids = set(seed_mapping["feature_id"].astype(str))
    anchor_feature_ids = {"A15_179732", "A15_168325"}
    prepared_feature_ids = feature_ids | anchor_feature_ids
    subset_diagnostics = {
        "consensus": write_feature_subset(
            sources["consensus_full"],
            input_dir / "consensus_aligned_strict16_mapped.csv",
            prepared_feature_ids,
        ),
        "annotations": write_feature_subset(
            sources["annotation_full"],
            input_dir / "unified_annotations_strict16_mapped.csv",
            prepared_feature_ids,
        ),
    }

    archlips_raw = pd.read_csv(sources["archlips_pos_rt_screened"], low_memory=False)
    if "rt_uncertain" not in archlips_raw.columns:
        raise ValueError("Release ArchLips table lacks the RT-screen flag")
    rt_uncertain = archlips_raw["rt_uncertain"].astype(str).str.lower().isin({"true", "1"})
    archlips_eligible = archlips_raw.loc[~rt_uncertain].copy()
    if len(archlips_eligible) != 546:
        raise ValueError(f"Expected 546 RT-screened POS ArchLips features, found {len(archlips_eligible)}")
    archlips_eligible_path = input_dir / "archlips_pos_release_eligible_rt_screened.csv"
    archlips_eligible.to_csv(archlips_eligible_path, index=False)

    pipeline_inputs = {
        "CONSENSUS_POS": input_dir / "consensus_aligned_strict16_mapped.csv",
        "METADATA_POS": sources["metadata_pos"],
        "SIMPER_ATLAS": sources["simper_atlas"],
        "CLIMGRASS_QUANT": sources["climgrass_quant"],
        "CLIMGRASS_META": sources["climgrass_meta"],
        "SPECTRAL_MATCHES": sources["spectral_matches"],
        "RIE_TABLE": sources["rie_table"],
        "EXPECTED_REF": sources["expected_ref"],
        "UNIFIED_ANNOT": input_dir / "unified_annotations_strict16_mapped.csv",
        "ARCHAEAL_IDENTIFICATION": archlips_eligible_path,
        "ARCHLIPS_VALIDATED": archlips_eligible_path,
    }
    for name, path in pipeline_inputs.items():
        setattr(legacy.pipeline, name, path.resolve())
    legacy.pipeline.PPM_TOL = 5.0
    legacy.pipeline.REQUIRE_ARCHLIPS_INPUTS = True
    legacy.pipeline._cache.clear()
    cache = legacy.pipeline.load_inputs()
    strict_units, unit_groups, taxonomy = configure_units(legacy, cache, sources)

    mapping = legacy.pipeline.get_verified_mapping(cache)
    if set(mapping["feature_id"].astype(str)) != feature_ids:
        raise ValueError("Prepared feature set disagrees with framework 5-ppm mapping")
    if not set(mapping["feature_id"].astype(str)).issubset(set(cache["cons"]["feature_id"].astype(str))):
        raise ValueError("Mapped strict features are missing from the prepared consensus subset")
    archlips = legacy.pipeline.get_archlips_feature_ids()
    archlips_intersection = len(feature_ids & archlips)
    atlas_archlips_intersection = len(set(atlas["feature_id"].astype(str)) & archlips)
    if args.archaea_mode == "archlips" and archlips_intersection == 0:
        raise ValueError(
            "Primary ArchLips mode is blocked: zero release-eligible RT-screened "
            "ArchLips features intersect the mapped strict-16 soil substrate"
        )

    soil_meta = cache["soil_meta"]
    cell_counts = soil_meta.groupby(["climate", "drought"]).size().to_dict()
    if len(soil_meta) != 12 or len(cell_counts) != 4 or set(cell_counts.values()) != {3}:
        raise ValueError(f"Expected balanced 2x2 design with n=3 per cell: {cell_counts}")

    original_atlas, original_soil, combined = legacy.install_anchor_accessors()
    anchors = legacy.anchor_diagnostics(cache, original_atlas, original_soil, combined)
    mapping.to_csv(staging / "strict16_verified_simper_mapping_5ppm.csv", index=False)
    anchors.to_csv(staging / "anchor_diagnostics.csv", index=False)
    mapping_contract = {
        "taxonomy_release": RELEASE_ID,
        "simper_mapping_ppm": 5.0,
        "corrected_spectral_match_rows": int(len(matches)),
        "unique_simper_features": int(mapping["feature_id"].nunique()),
        "ppm_max": float(mapping["ppm_diff"].max()),
        "archlips_mapped_feature_intersection": archlips_intersection,
        "archlips_atlas_feature_intersection": atlas_archlips_intersection,
        "archaea_mode": args.archaea_mode,
        "prepared_feature_subsets": subset_diagnostics,
        "prepared_internal_standard_feature_ids": sorted(anchor_feature_ids),
    }
    (staging / "mapping_contract.json").write_text(
        json.dumps(mapping_contract, indent=2) + "\n", encoding="utf-8"
    )

    summaries = []
    composition_validation = {}
    for label, anchor, display_name in ARMS:
        config = legacy.CorrectionConfig(
            IS_normalization=True,
            IS_reference_compound=anchor,
            IS_spiked_pmol=100.0,
            RIE_correction=True,
            RIE_floor=0.20,
            restrict_archaea_to_archlips=args.archaea_mode == "archlips",
        )
        arm_dir = staging / label
        summary = legacy.pipeline.run_iteration(config, arm_dir, label, include_bayesian=False)
        if summary["n_phyla"] != len(strict_units) or set(summary["phyla"]) != set(strict_units):
            raise ValueError(f"Arm did not produce the exact strict 16 phyla: {label}")
        composition_validation[label] = validate_compositions(arm_dir, label, strict_units)
        summaries.append(
            {
                "label": label,
                "display_name": display_name,
                "anchor": anchor,
                "config": asdict(config),
                "summary": summary,
            }
        )

    effects = legacy.calculate_effects(staging, sources["climgrass_meta"], unit_groups)
    effects.to_csv(staging / "phylum_effects_all_arms.csv", index=False)
    significance = (
        effects.groupby("arm", sort=False)
        .agg(
            n_units_tested=("unit", "size"),
            drought_p_lt_0_05=("drought_p", lambda values: int((values < 0.05).sum())),
            drought_q_lt_0_05=("drought_q_bh", lambda values: int((values < 0.05).sum())),
            climate_p_lt_0_05=("climate_p", lambda values: int((values < 0.05).sum())),
            climate_q_lt_0_05=("climate_q_bh", lambda values: int((values < 0.05).sum())),
        )
        .reset_index()
    )
    significance.to_csv(staging / "phylum_significance_summary.csv", index=False)
    plot_effects(effects, staging / "phylum_effect_map_strict16")

    acceptance = {
        "status": (
            "PASS_REVIEW_ONLY_STRICT16"
            if args.archaea_mode == "archlips"
            else "PASS_REVIEW_ONLY_STRICT16_ALL_SIMPER_SENSITIVITY"
        ),
        "taxonomy_release": RELEASE_ID,
        "spectral_qc": spectral_qc,
        "mapping_contract": mapping_contract,
        "taxonomy": taxonomy,
        "design_cell_counts": {f"{key[0]}|{key[1]}": int(value) for key, value in cell_counts.items()},
        "anchors": anchors.to_dict("records"),
        "composition_validation": composition_validation,
        "significance_summary": significance.to_dict("records"),
        "consumer_boundary": "review-only; not connected to manuscript, tables, figures, or response tracker",
        "primary_archlips_gate": {
            "status": "pass" if archlips_intersection else "blocked",
            "mapped_release_eligible_features": archlips_intersection,
            "atlas_release_eligible_features": atlas_archlips_intersection,
        },
    }
    (staging / "ACCEPTANCE_REPORT.json").write_text(
        json.dumps(acceptance, indent=2, default=str) + "\n", encoding="utf-8"
    )
    (staging / "run_summaries.json").write_text(
        json.dumps(summaries, indent=2, default=str) + "\n", encoding="utf-8"
    )
    write_readme(
        staging, spectral_qc, mapping_contract, taxonomy, significance, args.archaea_mode
    )

    provenance_inputs = dict(sources)
    provenance_inputs.update(
        {
            "producer": Path(__file__).resolve(),
            "recovered_runner": LEGACY_RUNNER,
            "framework_pipeline": HANDOFF_ROOT / "source" / "framework" / "pipeline.py",
            "framework_decomposition": HANDOFF_ROOT / "source" / "framework" / "decomposition.py",
            "framework_corrections": HANDOFF_ROOT / "source" / "framework" / "corrections.py",
            "framework_evaluation": HANDOFF_ROOT / "source" / "framework" / "evaluation.py",
        }
    )
    manifest = {
        "taxonomy_release": RELEASE_ID,
        "status": (
            "review_only"
            if args.archaea_mode == "archlips"
            else "review_only_all_simper_sensitivity"
        ),
        "inputs": {
            name: {"path": str(path), "sha256": sha256_file(path)}
            for name, path in sorted(provenance_inputs.items())
        },
        "outputs": {
            str(path.relative_to(staging)): sha256_file(path)
            for path in sorted(staging.rglob("*"))
            if path.is_file() and path.name != "RUN_MANIFEST.json"
        },
    }
    (staging / "RUN_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(staging), str(output))

    print("Strict 16-phylum ClimGrass review-only run complete")
    print(significance.to_string(index=False))
    print(f"Output: {output}")


if __name__ == "__main__":
    main()
