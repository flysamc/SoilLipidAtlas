#!/usr/bin/env python3
"""ArchLips-extended ClimGrass substrate: fingerprint = SIMPER union validated
diagnostic markers.

Finding that motivates this run: of the 546 MS2-validated, RT-screened,
release-eligible ArchLips compounds, exactly one survives SIMPER selection into
the 12,020-feature atlas and zero survive to the 722-feature soil substrate.
Contribution-ranked (SIMPER) fingerprints systematically under-select
low-abundance diagnostic chemistry; archaea are the extreme case. Meanwhile the
corrected soil-to-atlas spectral matches contain hits at ArchLips masses.

Design (additive union; evidence standard unchanged):
  - The SIMPER 5-ppm mapping is reproduced EXACTLY (722 features); the accepted
    strict16 run is re-derived and diffed cell-for-cell before anything else.
  - The same 2,313 corrected matches are additionally mapped against the 546
    release-eligible ArchLips feature masses at the same 5 ppm rule, deduped
    per feature by best cosine. These mappings are ADDED; nothing is replaced.
  - Mapped ArchLips features receive enriched status for their validated
    phylum with the substrate's median enriched SIMPER weight (they carry no
    SIMPER fold change by construction).
  - Arms: E0 baseline reproduction (722, all-SIMPER archaea reference);
    E1 extended substrate, all-SIMPER archaea reference;
    E2 extended substrate, ArchLips-RESTRICTED archaea reference (the
    originally intended primary contract, blocked until now because the
    substrate contained zero ArchLips features).

Review-only; not connected to the strict DAG.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RELEASE_ROOT = PROJECT_ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1"
OUTPUT = RELEASE_ROOT / "climgrass" / "strict16_archlips_extended_2026-08-08_v1"
RULEFIX = PROJECT_ROOT / "paper2_repro" / "scripts" / "climgrass_strict16_rulefix.py"
HANDOFF_ROOT = PROJECT_ROOT / "cg" / "SHADOW_PC_CLIMGRASS_STRICT_PHYLA_HANDOFF_2026-08-06"

PC_NAME = "PC_15-18d7"
METHODS = ("nnls", "standard_bc", "enriched_only_bc", "fc_weighted_bc")


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def map_matches_to_targets(matches: pd.DataFrame, target_mz: np.ndarray,
                           target_ids: np.ndarray, ppm: float = 5.0) -> pd.DataFrame:
    order = np.argsort(target_mz)
    mz_s, id_s = target_mz[order], target_ids[order]
    rows = []
    for r in matches.itertuples(index=False):
        amz = float(r.atlas_mz)
        tol = amz * ppm / 1e6
        lo = np.searchsorted(mz_s, amz - tol, "left")
        hi = np.searchsorted(mz_s, amz + tol, "right")
        if lo >= hi:
            continue
        cand = np.arange(lo, hi)
        ppm_diffs = np.abs(mz_s[cand] - amz) / amz * 1e6
        best = cand[np.argmin(ppm_diffs)]
        rows.append({"soil_scan": r.soil_scan, "atlas_mz": amz, "cosine": float(r.cosine),
                     "feature_id": str(id_s[best]), "ppm_diff": float(ppm_diffs.min())})
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    return (frame.sort_values("cosine", ascending=False)
            .drop_duplicates(subset="feature_id", keep="first").reset_index(drop=True))


def main() -> int:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite: {OUTPUT}")
    staging = OUTPUT.with_name(OUTPUT.name + ".incomplete")
    staging.mkdir(parents=True)

    rf = load_module("rulefix", RULEFIX)
    producer = rf.load_module("strict16_producer", rf.PRODUCER_PATH)
    legacy = producer.load_legacy_runner()
    prepared = rf.BASELINE_RUN / "prepared_inputs"
    sources = producer.required_sources(
        HANDOFF_ROOT / "results" / "strict_19unit_corrected_2026-08-06" / "corrected_spectral_matches.csv"
    )
    archlips_path = prepared / "archlips_pos_release_eligible_rt_screened.csv"

    # ------------------------------------------------------------------
    # Build the extended mapping and prepared inputs
    # ------------------------------------------------------------------
    matches, spectral_qc = legacy.validate_spectral_matches(sources["spectral_matches"])
    atlas = pd.read_csv(sources["simper_atlas"], low_memory=False)
    simper_mapping = producer.map_spectral_matches(atlas, matches)
    print(f"SIMPER mapping: {simper_mapping['feature_id'].nunique()} features (expect 722)")

    archlips = pd.read_csv(archlips_path, low_memory=False)
    arch_map = map_matches_to_targets(
        matches,
        archlips["query_mz"].astype(float).to_numpy(),
        archlips["feature_id"].astype(str).to_numpy(),
    )
    simper_ids = set(simper_mapping["feature_id"].astype(str))
    arch_map["also_in_simper_mapping"] = arch_map["feature_id"].isin(simper_ids)
    arch_new = arch_map[~arch_map["also_in_simper_mapping"]].copy()
    arch_meta = archlips.set_index(archlips["feature_id"].astype(str))
    arch_new["archlips_name"] = arch_new["feature_id"].map(arch_meta["archlips_name"])
    arch_new["archlips_phylum"] = arch_new["feature_id"].map(arch_meta["phylum"])
    print(f"ArchLips mapping: {len(arch_map)} features hit at 5 ppm, "
          f"{len(arch_new)} new beyond the SIMPER set")
    print(arch_new[["feature_id", "archlips_name", "archlips_phylum", "cosine", "ppm_diff"]]
          .round(3).to_string(index=False))
    if arch_new.empty:
        print("No new ArchLips features map — nothing to extend. Stopping.")
        return 1
    arch_new.to_csv(staging / "archlips_added_mapping.csv", index=False)
    simper_mapping.to_csv(staging / "simper_mapping_722.csv", index=False)

    union_mapping = pd.concat([
        simper_mapping[["soil_scan", "atlas_mz", "cosine", "feature_id", "ppm_diff"]],
        arch_new[["soil_scan", "atlas_mz", "cosine", "feature_id", "ppm_diff"]],
    ], ignore_index=True)

    input_dir = staging / "prepared_inputs"
    input_dir.mkdir()
    anchor_ids = {"A15_179732", "A15_168325"}
    union_ids = set(union_mapping["feature_id"].astype(str))
    producer.write_feature_subset(sources["consensus_full"],
                                  input_dir / "consensus_union.csv",
                                  union_ids | anchor_ids)
    producer.write_feature_subset(sources["annotation_full"],
                                  input_dir / "annotations_union.csv",
                                  union_ids | anchor_ids)

    # ------------------------------------------------------------------
    # Configure the framework on the union subset
    # ------------------------------------------------------------------
    pipeline_inputs = {
        "CONSENSUS_POS": input_dir / "consensus_union.csv",
        "METADATA_POS": sources["metadata_pos"],
        "SIMPER_ATLAS": sources["simper_atlas"],
        "CLIMGRASS_QUANT": sources["climgrass_quant"],
        "CLIMGRASS_META": sources["climgrass_meta"],
        "SPECTRAL_MATCHES": sources["spectral_matches"],
        "RIE_TABLE": sources["rie_table"],
        "EXPECTED_REF": sources["expected_ref"],
        "UNIFIED_ANNOT": input_dir / "annotations_union.csv",
        "ARCHAEAL_IDENTIFICATION": archlips_path,
        "ARCHLIPS_VALIDATED": archlips_path,
    }
    for name, path in pipeline_inputs.items():
        setattr(legacy.pipeline, name, Path(path).resolve())
    legacy.pipeline.PPM_TOL = 5.0
    legacy.pipeline.REQUIRE_ARCHLIPS_INPUTS = True
    legacy.pipeline._cache.clear()
    cache = legacy.pipeline.load_inputs()
    strict_units, unit_groups, taxonomy = producer.configure_units(legacy, cache, sources)
    legacy.install_anchor_accessors()

    original_get_mapping = legacy.pipeline.get_verified_mapping

    def set_mapping(frame):
        legacy.pipeline.get_verified_mapping = lambda cache_: frame.copy()

    # ------------------------------------------------------------------
    # E0: reproduction gate on the unchanged 722 substrate
    # ------------------------------------------------------------------
    def make_config(restrict):
        return legacy.CorrectionConfig(
            IS_normalization=True, IS_reference_compound=PC_NAME, IS_spiked_pmol=100.0,
            RIE_correction=True, RIE_floor=0.20,
            restrict_archaea_to_archlips=restrict,
        )

    arms_dir = staging / "arms"
    arms_dir.mkdir()
    set_mapping(simper_mapping)
    print("\n=== E0 baseline reproduction (722, all-SIMPER archaea) ===")
    legacy.pipeline.run_iteration(make_config(False), arms_dir / "E0_baseline722", "E0_baseline722",
                                  include_bayesian=False)
    reproduction = {}
    for method in METHODS:
        mine = pd.read_csv(arms_dir / "E0_baseline722" / f"phylum_composition_E0_baseline722_{method}.csv",
                           index_col=0)
        theirs = pd.read_csv(rf.BASELINE_RUN / "strict16_PC_primary"
                             / f"phylum_composition_strict16_PC_primary_{method}.csv", index_col=0)
        reproduction[method] = float((mine[theirs.columns] - theirs).abs().to_numpy().max())
    max_delta = max(reproduction.values())
    print(f"REPRODUCTION: max |delta| = {max_delta:.3e} -> {'PASS' if max_delta < 1e-12 else 'FAIL'}")
    (staging / "reproduction_check.json").write_text(
        json.dumps({"per_method": reproduction, "max": max_delta, "pass": max_delta < 1e-12},
                   indent=2) + "\n", encoding="utf-8")
    if max_delta >= 1e-12:
        print("Stopping: extension effects would not be attributable.")
        return 1

    # ------------------------------------------------------------------
    # Extended substrate: augment SIMPER table with diagnostic-marker rows
    # ------------------------------------------------------------------
    median_fc = float(cache["simper"].loc[
        (cache["simper"]["direction"] == "enriched")
        & cache["simper"]["feature_id"].isin(simper_ids), "fold_change"].median())
    synthetic = pd.DataFrame({
        "feature_id": arch_new["feature_id"].astype(str),
        "phylum": arch_new["archlips_phylum"].astype(str),
        "direction": "enriched",
        "fold_change": median_fc,
    })
    # Only phyla present in the strict scheme can carry weight
    synthetic = synthetic[synthetic["phylum"].isin(strict_units)]
    print(f"\nDiagnostic-marker enrichment rows added: {len(synthetic)} "
          f"(median SIMPER enriched fold change {median_fc:.1f})")
    cache["simper"] = pd.concat([cache["simper"], synthetic], ignore_index=True)

    set_mapping(union_mapping)
    print("\n=== E1 extended substrate, all-SIMPER archaea reference ===")
    legacy.pipeline.run_iteration(make_config(False), arms_dir / "E1_extended_allsimper",
                                  "E1_extended_allsimper", include_bayesian=False)
    print("\n=== E2 extended substrate, ArchLips-RESTRICTED archaea reference (primary) ===")
    legacy.pipeline.run_iteration(make_config(True), arms_dir / "E2_extended_archlips_restricted",
                                  "E2_extended_archlips_restricted", include_bayesian=False)

    # ------------------------------------------------------------------
    # Comparison table
    # ------------------------------------------------------------------
    rows = []
    for arm in ("E0_baseline722", "E1_extended_allsimper", "E2_extended_archlips_restricted"):
        for method in METHODS:
            king = pd.read_csv(arms_dir / arm / f"kingdom_composition_{arm}_{method}.csv", index_col=0)
            entry = {"arm": arm, "method": method}
            entry.update({g: float(v) * 100 for g, v in king.mean().items()})
            rows.append(entry)
    comparison = pd.DataFrame(rows)
    comparison.to_csv(staging / "kingdom_comparison.csv", index=False)
    print("\n=== kingdom means (%) ===")
    print(comparison.round(2).to_string(index=False))

    summary = {
        "spectral_qc": spectral_qc,
        "simper_features": int(simper_mapping["feature_id"].nunique()),
        "archlips_hit_at_5ppm": int(len(arch_map)),
        "archlips_added": int(len(arch_new)),
        "added_features": arch_new[["feature_id", "archlips_name", "archlips_phylum",
                                    "cosine", "ppm_diff"]].to_dict("records"),
        "reproduction_max_delta": max_delta,
        "median_simper_enriched_fold_change_assigned": median_fc,
        "arms": {
            "E0_baseline722": "722 SIMPER features, all-SIMPER archaea reference (accepted run)",
            "E1_extended_allsimper": "union substrate, archaea reference unchanged",
            "E2_extended_archlips_restricted": "union substrate, archaea reference = ArchLips features only (primary contract, previously blocked)",
        },
    }
    (staging / "RUN_SUMMARY.json").write_text(json.dumps(summary, indent=2, default=str) + "\n",
                                              encoding="utf-8")
    staging.rename(OUTPUT)
    print(f"\nOutput: {OUTPUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
