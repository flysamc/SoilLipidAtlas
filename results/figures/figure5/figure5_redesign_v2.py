#!/usr/bin/env python3
"""Figure 5 redesign v2: the full estimator/statistics machinery of
figure5_redesign.py applied to the ArchLips-EXTENDED substrate (736 features)
with the ArchLips-restricted archaeal reference (the unblocked primary
contract from strict16_archlips_extended_2026-08-08_v1).

Deltas vs v1 (everything else identical, same seed):
  - substrate: 722 SIMPER features + 14 spectrally matched archaeol-family
    features (additive union mapping, same 2,313 matches, same 5 ppm rule);
  - archaeal reference columns are masked to ArchLips-validated features only
    (diagnostic chemistry), mirroring pipeline restrict_archaea_to_archlips;
  - the 14 diagnostic markers carry enriched status for their validated phylum
    at the substrate's median enriched SIMPER weight (rule C still applies).
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
OUTPUT = RELEASE_ROOT / "climgrass" / "figure5_redesign_2026-08-08_v2_archlips"
EXTENSION = RELEASE_ROOT / "climgrass" / "strict16_archlips_extended_2026-08-08_v1"
FIG5_V1 = PROJECT_ROOT / "paper2_repro" / "scripts" / "figure5_redesign.py"
HANDOFF_ROOT = PROJECT_ROOT / "cg" / "SHADOW_PC_CLIMGRASS_STRICT_PHYLA_HANDOFF_2026-08-06"

ARCHAEA_PHYLA = {"Methanobacteriota", "Thermoproteota", "Halobacteriota"}


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def build_inputs_extended(f5):
    rf = f5.load_module("rulefix", f5.RULEFIX)
    producer = rf.load_module("strict16_producer", rf.PRODUCER_PATH)
    legacy = producer.load_legacy_runner()
    prepared = EXTENSION / "prepared_inputs"
    sources = producer.required_sources(
        HANDOFF_ROOT / "results" / "strict_19unit_corrected_2026-08-06" / "corrected_spectral_matches.csv"
    )
    archlips_path = (rf.BASELINE_RUN / "prepared_inputs"
                     / "archlips_pos_release_eligible_rt_screened.csv")
    pipeline_inputs = {
        "CONSENSUS_POS": prepared / "consensus_union.csv",
        "METADATA_POS": sources["metadata_pos"],
        "SIMPER_ATLAS": sources["simper_atlas"],
        "CLIMGRASS_QUANT": sources["climgrass_quant"],
        "CLIMGRASS_META": sources["climgrass_meta"],
        "SPECTRAL_MATCHES": sources["spectral_matches"],
        "RIE_TABLE": sources["rie_table"],
        "EXPECTED_REF": sources["expected_ref"],
        "UNIFIED_ANNOT": prepared / "annotations_union.csv",
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

    # Union mapping = frozen SIMPER mapping + the 14 archaeol additions
    simper_map = pd.read_csv(EXTENSION / "simper_mapping_722.csv")
    arch_map = pd.read_csv(EXTENSION / "archlips_added_mapping.csv")
    union = pd.concat([
        simper_map[["soil_scan", "atlas_mz", "cosine", "feature_id", "ppm_diff"]],
        arch_map[["soil_scan", "atlas_mz", "cosine", "feature_id", "ppm_diff"]],
    ], ignore_index=True)

    pl = legacy.pipeline
    atlas_df, soil_df = pl.build_matrices(cache, union)
    if len(atlas_df) != 736:
        raise ValueError(f"Expected 736 union features, got {len(atlas_df)}")

    config = rf.make_config(legacy, rf.PC_NAME)
    atlas_is = pl.get_atlas_IS_intensity(cache, rf.PC_NAME)
    soil_is = pl.get_soil_IS_intensity(cache, rf.PC_NAME)
    both = pd.concat([atlas_is[atlas_is > 0], soil_is[soil_is > 0]])
    global_median = float(both.median())
    atlas_corr = pl.apply_IS_normalization(atlas_df, atlas_is, config.IS_spiked_pmol, global_median)
    soil_corr = pl.apply_IS_normalization(soil_df, soil_is, config.IS_spiked_pmol, global_median)

    feature_ids = list(atlas_corr.index)
    classes, adducts = pl.assign_class_adduct(feature_ids, cache["annot"])
    rule_a = rf.make_rule_a(pl.apply_RIE_correction)
    atlas_corr, _ = rule_a(atlas_corr, classes, adducts, cache["rie"], 1.0,
                           rie_floor=rf.RIE_FLOOR, rie_ceiling=rf.RIE_CEILING)
    soil_corr, _ = rule_a(soil_corr, classes, adducts, cache["rie"], 1.0,
                          rie_floor=rf.RIE_FLOOR, rie_ceiling=rf.RIE_CEILING)

    sample_phylum = {s: p for s, p in cache["sample_phylum"].items() if s in atlas_corr.columns}
    phyla = sorted(set(sample_phylum.values()))

    # ArchLips-restricted archaeal reference: in archaea reference columns,
    # non-ArchLips features carry no weight (diagnostic chemistry only).
    archlips_ids = pl.get_archlips_feature_ids()
    archaea_cols = [s for s, p in sample_phylum.items() if p in ARCHAEA_PHYLA]
    non_archlips = ~atlas_corr.index.astype(str).isin(archlips_ids)
    n_arch_features = int((~non_archlips).sum())
    atlas_corr = atlas_corr.copy()
    atlas_corr.loc[non_archlips, archaea_cols] = 0.0
    print(f"ArchLips-restricted archaea reference: {n_arch_features} diagnostic features "
          f"across {len(archaea_cols)} archaeal reference samples")

    # Diagnostic-marker enrichment rows (median enriched SIMPER weight)
    simper_ids = set(simper_map["feature_id"].astype(str))
    median_fc = float(cache["simper"].loc[
        (cache["simper"]["direction"] == "enriched")
        & cache["simper"]["feature_id"].isin(simper_ids), "fold_change"].median())
    synthetic = pd.DataFrame({
        "feature_id": arch_map["feature_id"].astype(str),
        "phylum": arch_map["archlips_phylum"].astype(str),
        "direction": "enriched",
        "fold_change": median_fc,
    })
    synthetic = synthetic[synthetic["phylum"].isin(phyla)]
    simper_aug = pd.concat([cache["simper"], synthetic], ignore_index=True)

    from decomposition import build_direction_masks
    enriched_masks, fc_weights = build_direction_masks(simper_aug, phyla, feature_ids)
    n_enriched = np.zeros(len(feature_ids))
    for p in phyla:
        n_enriched += enriched_masks[p].astype(float)
    divisor = np.maximum(n_enriched, 1.0)
    for p in phyla:
        w = fc_weights[p].copy()
        m = enriched_masks[p]
        w[m] = w[m] / divisor[m]
        fc_weights[p] = w

    meta = pd.read_csv(sources["climgrass_meta"]).set_index("column_name")
    return dict(
        atlas=atlas_corr, soil=soil_corr, sample_phylum=sample_phylum, phyla=phyla,
        feature_ids=feature_ids, enriched_masks=enriched_masks, fc_weights=fc_weights,
        unit_groups=unit_groups, soil_meta=meta, expected=cache["expected"],
        n_enriched=n_enriched,
    )


def main() -> int:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite: {OUTPUT}")
    staging = OUTPUT.with_name(OUTPUT.name + ".incomplete")
    staging.mkdir(parents=True)

    f5 = load_module("figure5_v1", FIG5_V1)

    print("Step 1: extended corrected matrices...")
    inputs = build_inputs_extended(f5)
    print(f"  atlas {inputs['atlas'].shape}, soil {inputs['soil'].shape}, phyla {len(inputs['phyla'])}")

    print("Step 2: held-out mixture benchmark on extended substrate...")
    bench, bench_summary = f5.mixture_benchmark(inputs)
    bench.to_csv(staging / "benchmark_mixtures.csv", index=False)
    bench_summary.to_csv(staging / "benchmark_summary.csv", index=False)
    print(bench_summary.round(4).to_string(index=False))
    winner = bench_summary.iloc[0]
    est_name, sqrt_t = str(winner["estimator"]), bool(winner["sqrt"])
    winner_label = f"{est_name}{' + sqrt' if sqrt_t else ''}"
    print(f"  WINNER: {winner_label}  (score={winner['score']:.4f})")

    print("Step 3: real-soil composition + bootstrap...")
    comp_df, boot, boot_kingdom = f5.soil_composition(inputs, est_name, sqrt_t)
    comp_df.to_csv(staging / "composition_winner_by_sample.csv")
    boot.to_csv(staging / "composition_winner_bootstrap_ci.csv", index=False)
    boot_kingdom.to_csv(staging / "composition_winner_kingdom_ci.csv", index=False)
    comp_pub, boot_pub, boot_pub_king = f5.soil_composition(inputs, "fc_weighted_bc", False)
    comp_pub.to_csv(staging / "composition_fcweighted_by_sample.csv")
    boot_pub.to_csv(staging / "composition_fcweighted_bootstrap_ci.csv", index=False)
    boot_pub_king.to_csv(staging / "composition_fcweighted_kingdom_ci.csv", index=False)
    print("Winner kingdom means (%):")
    print((boot_kingdom.set_index("kingdom") * 100).round(1).to_string())

    print("Step 4: fingerprint-set treatment tests...")
    effects = f5.set_tests(inputs)
    effects.to_csv(staging / "fingerprint_set_effects.csv", index=False)
    print(effects.sort_values(["factor", "p_perm"]).head(10).round(4).to_string(index=False))
    replication = f5.replication_test(effects)
    replication.to_csv(staging / "qsip_replication_test.csv", index=False)
    print("\nqSIP replication (family of 2):")
    print(replication.round(4).to_string(index=False))

    print("Step 5: figure...")
    note = ("bars: benchmark-selected estimator, 95% CI\n"
            "(reference-sample bootstrap); diamonds:\n"
            "fc-weighted estimator with rules; grey bars:\n"
            "literature ranges. † Archaea measured from\n"
            "ArchLips-validated ether lipids only (14 diagnostic\n"
            "markers); scale uncertain, no ether-lipid RIE standards")
    f5.render_figure(comp_df, boot_kingdom, boot_pub_king, effects, replication,
                     inputs["unit_groups"], winner_label, staging / "Figure5_proposal_v2",
                     note_text=note, note_pos=(0.98, 0.05))

    summary = {
        "substrate": "736 features: 722 SIMPER + 14 spectrally matched archaeol markers; "
                     "ArchLips-restricted archaeal reference",
        "upstream": str(EXTENSION),
        "benchmark_winner": winner_label,
        "winner_score": float(winner["score"]),
        "kingdom_ci": boot_kingdom.to_dict("records"),
        "top_drought_hits": effects[effects["factor"] == "drought"].nsmallest(4, "p_perm")[
            ["unit", "set_mean_log2fc", "p_perm", "q_bh"]].to_dict("records"),
        "qsip_replication": replication.to_dict("records"),
    }
    (staging / "RUN_SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    staging.rename(OUTPUT)
    print(f"\nOutput: {OUTPUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
