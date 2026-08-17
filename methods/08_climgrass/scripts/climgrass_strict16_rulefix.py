#!/usr/bin/env python3
"""Rule-based correction arms for the strict-16 ClimGrass decomposition.

Reproduce-first: the baseline arm re-runs the exact configuration of the
accepted `strict16_all_simper_sensitivity_2026-08-06` run from its own prepared
inputs and must match the stored composition CSVs cell-for-cell before any
rule arm is interpreted.

Rule arms (direction-blind, applied identically to every phylum):
  ruleA  exact-class RIE outside [floor, ceiling] is treated as UNCALIBRATED
         (RIE = 1.0) instead of being clipped to the floor/ceiling. The floor
         clip amplifies floored features 1/floor = 5x, which concentrates
         leverage in a handful of features; out-of-range calibration is not
         evidence of a 5x response factor.
  ruleB  no single feature may carry more than DOMINANCE_CAP of a phylum
         reference profile's total mass (waterfill cap). RETIRED after v1:
         overlap dominance is co-driven by soil intensity, so the cap only
         moves dominance to the next feature.
  ruleC  specificity weighting: a feature SIMPER-enriched in k phyla has its
         enriched evidence weight divided by k (ambiguous evidence is split,
         not multiply counted). Implements the intervention-audit
         recommendation of ambiguity down-weighting rather than deletion.
  ruleS  Hellinger-style square-root transform of the corrected reference and
         soil profiles before decomposition; compresses dynamic range so no
         single feature can dominate an overlap.
  ruleD  chemistry discordance gate: archaeal phyla may not draw enriched
         evidence weight from features carrying a confident LipidSearch
         ester-lipid class annotation (TG, PC, PI, WE, ...). Archaeal membranes
         are isoprenoid ethers (archaeol/GDGT); ester-glycerolipid annotations
         on archaeal-enriched culture features are media/carryover or
         misannotation, not diagnostic archaeal chemistry. Unannotated features
         are retained (they may be unrecognized ether lipids). Zero ArchLips
         features intersect the soil substrate, so this is the only available
         archaeal chemistry gate.

For each arm x estimator (nnls, standard_bc, enriched_only_bc, fc_weighted_bc)
the script records kingdom means, the framework's own plausibility diagnostics
against expected_kingdom_composition.csv, drought/climate Mann-Whitney + BH
statistics, and per-phylum top-feature leverage (share of weighted overlap,
and leave-top-feature-out delta) so single-feature-conditional units stay
visible.
"""

from __future__ import annotations

import importlib.util
import json
from dataclasses import asdict
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RELEASE_ID = "ncbi-phylum-2026-08-04-v1"
RELEASE_ROOT = PROJECT_ROOT / "outputs" / "analysis" / RELEASE_ID
BASELINE_RUN = RELEASE_ROOT / "climgrass" / "strict16_all_simper_sensitivity_2026-08-06"
HANDOFF_ROOT = PROJECT_ROOT / "cg" / "SHADOW_PC_CLIMGRASS_STRICT_PHYLA_HANDOFF_2026-08-06"
PRODUCER_PATH = PROJECT_ROOT / "paper2_repro" / "scripts" / "climgrass_strict16_release.py"

OUTPUT = RELEASE_ROOT / "climgrass" / "strict16_rulefix_2026-08-08_v3"

PC_NAME = "PC_15-18d7"
PE_NAME = "PE_15-18d7_precursor"
COMBINED_NAME = "PC_PE_two_standard"

METHODS = ("nnls", "standard_bc", "enriched_only_bc", "fc_weighted_bc")
DETECTION_THRESHOLD = 0.005
EPSILON = 1e-4
DOMINANCE_CAP = 0.25
RIE_FLOOR = 0.20
RIE_CEILING = 100.0

# (arm_id, anchor, rules) with rules a subset of {"A", "C", "S", "D"}
ARM_SPECS = [
    ("baseline_PC", PC_NAME, set()),
    ("ruleACS_PC", PC_NAME, {"A", "C", "S"}),
    ("ruleD_PC", PC_NAME, {"D"}),
    ("ruleACD_PC", PC_NAME, {"A", "C", "D"}),
    ("ruleACSD_PC", PC_NAME, {"A", "C", "S", "D"}),
    ("ruleACSD_PE", PE_NAME, {"A", "C", "S", "D"}),
    ("ruleACSD_PCPE", COMBINED_NAME, {"A", "C", "S", "D"}),
]
PRIMARY_ARM = "ruleACSD_PC"

ARCHAEA_PHYLA = {"Methanobacteriota", "Thermoproteota", "Halobacteriota"}


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make_config(legacy, anchor: str):
    return legacy.CorrectionConfig(
        IS_normalization=True,
        IS_reference_compound=anchor,
        IS_spiked_pmol=100.0,
        RIE_correction=True,
        RIE_floor=RIE_FLOOR,
        restrict_archaea_to_archlips=False,  # matches the accepted all-SIMPER sensitivity run
    )


# ---------------------------------------------------------------------------
# Rule A: out-of-range exact RIE -> uncalibrated (1.0)
# ---------------------------------------------------------------------------
def make_rule_a(original_apply_rie):
    def apply_rie_rule_a(intensity, feature_class, feature_adduct, rie_lookup,
                         fallback_RIE=1.0, rie_floor=0.05, rie_ceiling=100.0):
        # Reuse the framework's exact lookup by disabling the clip, then apply
        # the uncalibrated rule ourselves.
        _, rie_raw = original_apply_rie(
            intensity, feature_class, feature_adduct, rie_lookup,
            fallback_RIE, rie_floor=1e-12, rie_ceiling=1e12,
        )
        out_of_range = (rie_raw < rie_floor) | (rie_raw > rie_ceiling)
        rie_final = rie_raw.where(~out_of_range, fallback_RIE)
        corrected = intensity.div(rie_final.values, axis=0)
        return corrected, rie_final
    return apply_rie_rule_a


# ---------------------------------------------------------------------------
# Rule C: specificity weighting (enriched weight split across k phyla)
# ---------------------------------------------------------------------------
def make_rule_c(original_masks):
    def build_direction_masks_rule_c(simper_full, phyla, feature_ids):
        enriched_masks, fc_weights = original_masks(simper_full, phyla, feature_ids)
        n_features = len(feature_ids)
        n_enriched = np.zeros(n_features)
        for p in phyla:
            n_enriched += enriched_masks[p].astype(float)
        divisor = np.maximum(n_enriched, 1.0)
        for p in phyla:
            mask = enriched_masks[p]
            fc_weights[p] = fc_weights[p].copy()
            fc_weights[p][mask] = fc_weights[p][mask] / divisor[mask]
        return enriched_masks, fc_weights
    return build_direction_masks_rule_c


# ---------------------------------------------------------------------------
# Rule D: archaeal ester-lipid discordance gate (weight level)
# ---------------------------------------------------------------------------
def make_rule_d(original_masks, feature_class_map):
    def build_direction_masks_rule_d(simper_full, phyla, feature_ids):
        enriched_masks, fc_weights = original_masks(simper_full, phyla, feature_ids)
        annotated = np.array([
            pd.notna(feature_class_map.get(f)) for f in feature_ids
        ])
        for p in phyla:
            if p not in ARCHAEA_PHYLA:
                continue
            mask = enriched_masks[p].copy()
            discordant = mask & annotated
            mask[discordant] = False
            weights = fc_weights[p].copy()
            weights[discordant] = 0.01  # same baseline as non-enriched features
            enriched_masks[p] = mask
            fc_weights[p] = weights
        return enriched_masks, fc_weights
    return build_direction_masks_rule_d


# ---------------------------------------------------------------------------
# Rule S: Hellinger-style sqrt transform of corrected profiles
# ---------------------------------------------------------------------------
def make_rule_s_ref(original_builder):
    def build_reference_rule_s(atlas_intensity_df, sample_phylum, min_samples=2):
        ref, phyla, feature_ids = original_builder(
            atlas_intensity_df, sample_phylum, min_samples=min_samples
        )
        return np.sqrt(np.clip(ref, 0, None)), phyla, feature_ids
    return build_reference_rule_s


def make_rule_s_soil(original_builder):
    def build_soil_rule_s(soil_intensity_df, feature_ids):
        soil, cols = original_builder(soil_intensity_df, feature_ids)
        return np.sqrt(np.clip(soil, 0, None)), cols
    return build_soil_rule_s


class ArmPatches:
    """Apply/remove rule patches on the loaded framework modules."""

    def __init__(self, legacy, rules: set):
        self.legacy = legacy
        self.rules = set(rules)
        self._orig_rie = legacy.pipeline.apply_RIE_correction
        self._orig_ref_builder = legacy.pipeline.build_phylum_reference_array
        self._orig_soil_builder = legacy.pipeline.build_soil_matrix_array
        self._orig_masks = legacy.decomposition.build_direction_masks

    def __enter__(self):
        if "A" in self.rules:
            self.legacy.pipeline.apply_RIE_correction = make_rule_a(self._orig_rie)
        if "C" in self.rules:
            self.legacy.decomposition.build_direction_masks = make_rule_c(self._orig_masks)
        if "D" in self.rules:
            annot = self.legacy.pipeline._cache["annot"]
            class_map = dict(zip(annot["feature_id"].astype(str), annot["ls_ClassKey"]))
            current = self.legacy.decomposition.build_direction_masks
            self.legacy.decomposition.build_direction_masks = make_rule_d(current, class_map)
        if "S" in self.rules:
            self.legacy.pipeline.build_phylum_reference_array = make_rule_s_ref(self._orig_ref_builder)
            self.legacy.pipeline.build_soil_matrix_array = make_rule_s_soil(self._orig_soil_builder)
        return self

    def __exit__(self, *exc):
        self.legacy.pipeline.apply_RIE_correction = self._orig_rie
        self.legacy.pipeline.build_phylum_reference_array = self._orig_ref_builder
        self.legacy.pipeline.build_soil_matrix_array = self._orig_soil_builder
        self.legacy.decomposition.build_direction_masks = self._orig_masks
        return False


# ---------------------------------------------------------------------------
# Statistics across all estimators
# ---------------------------------------------------------------------------
def bh_adjust(values: pd.Series) -> pd.Series:
    p_values = values.astype(float).to_numpy()
    order = np.argsort(p_values)
    ranked = p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty(len(p_values), dtype=float)
    adjusted[order] = np.clip(ranked, 0, 1)
    return pd.Series(adjusted, index=values.index)


def effects_for(composition: pd.DataFrame, metadata: pd.DataFrame,
                unit_kingdom: dict, arm: str, method: str) -> pd.DataFrame:
    common = composition.index.intersection(metadata.index)
    if len(common) != 12:
        raise ValueError(f"Treatment metadata join covered {len(common)}/12 samples")
    composition = composition.loc[common]
    meta = metadata.loc[common]
    units = [c for c in composition if composition[c].max() >= DETECTION_THRESHOLD]
    rows = []
    for unit in units:
        values = composition[unit].astype(float)
        drought = values[meta["drought"] == "Drought"]
        control = values[meta["drought"] == "No_drought"]
        future = values[meta["climate"] == "Future"]
        ambient = values[meta["climate"] == "Ambient"]
        _, drought_p = mannwhitneyu(drought, control, alternative="two-sided")
        _, climate_p = mannwhitneyu(future, ambient, alternative="two-sided")
        rows.append({
            "arm": arm, "method": method, "unit": unit,
            "kingdom": unit_kingdom[unit],
            "mean_fraction": float(values.mean()),
            "drought_log2FC": float(np.log2((drought.mean() + EPSILON) / (control.mean() + EPSILON))),
            "climate_log2FC": float(np.log2((future.mean() + EPSILON) / (ambient.mean() + EPSILON))),
            "drought_p": float(drought_p),
            "climate_p": float(climate_p),
        })
    frame = pd.DataFrame(rows)
    frame["drought_q_bh"] = bh_adjust(frame["drought_p"])
    frame["climate_q_bh"] = bh_adjust(frame["climate_p"])
    return frame


# ---------------------------------------------------------------------------
# Leverage diagnostics (weighted-overlap share + leave-top-feature-out)
# ---------------------------------------------------------------------------
def build_corrected_matrices(legacy, config):
    """Mirror pipeline.run_iteration steps 1-4 exactly, using whatever patches
    are active, and return the decomposition inputs."""
    pl = legacy.pipeline
    cache = pl.load_inputs()
    mapped_df = pl.get_verified_mapping(cache)
    atlas_df, soil_df = pl.build_matrices(cache, mapped_df)
    atlas_corr, soil_corr = atlas_df.copy(), soil_df.copy()

    if config.IS_normalization:
        is_name = config.IS_reference_compound
        atlas_is = pl.get_atlas_IS_intensity(cache, is_name)
        soil_is = pl.get_soil_IS_intensity(cache, is_name)
        both = pd.concat([atlas_is[atlas_is > 0], soil_is[soil_is > 0]])
        global_median = float(both.median()) if len(both) else 1.0
        atlas_corr = pl.apply_IS_normalization(atlas_corr, atlas_is, config.IS_spiked_pmol, global_median)
        soil_corr = pl.apply_IS_normalization(soil_corr, soil_is, config.IS_spiked_pmol, global_median)

    if config.RIE_correction:
        feature_ids = list(atlas_corr.index)
        classes, adducts = pl.assign_class_adduct(feature_ids, cache["annot"])
        atlas_corr, rie_vec = pl.apply_RIE_correction(
            atlas_corr, classes, adducts, cache["rie"], config.fallback_RIE,
            rie_floor=config.RIE_floor, rie_ceiling=config.RIE_ceiling,
        )
        soil_corr, _ = pl.apply_RIE_correction(
            soil_corr, classes, adducts, cache["rie"], config.fallback_RIE,
            rie_floor=config.RIE_floor, rie_ceiling=config.RIE_ceiling,
        )
    else:
        rie_vec = None

    ref_arr, phyla, feature_ids = pl.build_phylum_reference_array(
        atlas_corr, cache["sample_phylum"]
    )
    soil_arr, soil_cols = pl.build_soil_matrix_array(soil_corr, feature_ids)
    return cache, ref_arr, soil_arr, phyla, feature_ids, soil_cols, rie_vec


def leverage_table(legacy, config, unit_kingdom) -> pd.DataFrame:
    """Per phylum: top feature by mean weighted-overlap share (fc-weighted BC
    contribution proxy) and the kingdom-level delta when it is removed."""
    from decomposition import build_direction_masks, decompose_fc_weighted_bc, phylum_to_kingdom_array

    cache, ref, soil, phyla, feature_ids, soil_cols, _ = build_corrected_matrices(legacy, config)
    _, fc_weights = build_direction_masks(cache["simper"], phyla, feature_ids)

    baseline = decompose_fc_weighted_bc(ref, soil, fc_weights, phyla)
    base_kingdom, kingdoms = phylum_to_kingdom_array(baseline, phyla)
    base_mean = pd.Series(base_kingdom.mean(axis=0), index=kingdoms)

    rows = []
    for p_idx, phylum in enumerate(phyla):
        w = fc_weights.get(phylum, np.ones(ref.shape[1]))
        wn = w / (w.sum() + 1e-10)
        shares = np.zeros(ref.shape[1])
        for i in range(soil.shape[0]):
            overlap = wn * np.minimum(soil[i], ref[p_idx])
            total = overlap.sum()
            if total > 0:
                shares += overlap / total
        shares /= soil.shape[0]
        top_idx = int(np.argmax(shares))
        top_share = float(shares[top_idx])

        keep = np.ones(ref.shape[1], dtype=bool)
        keep[top_idx] = False
        fc_weights_loo = {p: v[keep] for p, v in fc_weights.items()}
        loo = decompose_fc_weighted_bc(ref[:, keep], soil[:, keep], fc_weights_loo, phyla)
        loo_kingdom, _ = phylum_to_kingdom_array(loo, phyla)
        loo_mean = pd.Series(loo_kingdom.mean(axis=0), index=kingdoms)
        kingdom = unit_kingdom[phylum]
        rows.append({
            "phylum": phylum,
            "kingdom": kingdom,
            "top_feature": feature_ids[top_idx],
            "top_overlap_share": top_share,
            "phylum_mean_fraction": float(baseline[:, p_idx].mean()),
            "phylum_mean_fraction_loo": float(loo[:, p_idx].mean()),
            "kingdom_delta_pp_on_removal": float((loo_mean[kingdom] - base_mean[kingdom]) * 100),
            "single_feature_conditional": top_share > DOMINANCE_CAP,
        })
    return pd.DataFrame(rows).sort_values("top_overlap_share", ascending=False)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite: {OUTPUT}")

    producer = load_module("strict16_producer", PRODUCER_PATH)
    legacy = producer.load_legacy_runner()

    # Point the framework at the SAME inputs as the accepted baseline run,
    # including its prepared feature subsets.
    prepared = BASELINE_RUN / "prepared_inputs"
    sources = producer.required_sources(
        HANDOFF_ROOT / "results" / "strict_19unit_corrected_2026-08-06" / "corrected_spectral_matches.csv"
    )
    pipeline_inputs = {
        "CONSENSUS_POS": prepared / "consensus_aligned_strict16_mapped.csv",
        "METADATA_POS": sources["metadata_pos"],
        "SIMPER_ATLAS": sources["simper_atlas"],
        "CLIMGRASS_QUANT": sources["climgrass_quant"],
        "CLIMGRASS_META": sources["climgrass_meta"],
        "SPECTRAL_MATCHES": sources["spectral_matches"],
        "RIE_TABLE": sources["rie_table"],
        "EXPECTED_REF": sources["expected_ref"],
        "UNIFIED_ANNOT": prepared / "unified_annotations_strict16_mapped.csv",
        "ARCHAEAL_IDENTIFICATION": prepared / "archlips_pos_release_eligible_rt_screened.csv",
        "ARCHLIPS_VALIDATED": prepared / "archlips_pos_release_eligible_rt_screened.csv",
    }
    for name, path in pipeline_inputs.items():
        if not Path(path).is_file():
            raise FileNotFoundError(path)
        setattr(legacy.pipeline, name, Path(path).resolve())
    legacy.pipeline.PPM_TOL = 5.0
    legacy.pipeline.REQUIRE_ARCHLIPS_INPUTS = True
    legacy.pipeline._cache.clear()
    cache = legacy.pipeline.load_inputs()
    strict_units, unit_groups, taxonomy = producer.configure_units(legacy, cache, sources)
    legacy.install_anchor_accessors()

    staging = OUTPUT.with_name(OUTPUT.name + ".incomplete")
    staging.mkdir(parents=True)
    arms_dir = staging / "arms"
    arms_dir.mkdir()

    metadata = pd.read_csv(sources["climgrass_meta"]).set_index("column_name")

    summaries = {}
    effects_frames = []
    kingdom_rows = []
    for arm_id, anchor, rules in ARM_SPECS:
        print(f"\n{'=' * 80}\nARM {arm_id} (anchor={anchor}, rules={sorted(rules) or 'none'})\n{'=' * 80}")
        config = make_config(legacy, anchor)
        arm_dir = arms_dir / arm_id
        with ArmPatches(legacy, rules):
            summary = legacy.pipeline.run_iteration(config, arm_dir, arm_id, include_bayesian=False)
        if set(summary["phyla"]) != set(strict_units):
            raise ValueError(f"Arm {arm_id} did not produce the strict 16 phyla")
        summaries[arm_id] = {"config": asdict(config), "summary": summary}

        for method in METHODS:
            comp = pd.read_csv(arm_dir / f"phylum_composition_{arm_id}_{method}.csv", index_col=0)
            effects_frames.append(effects_for(comp, metadata, unit_groups, arm_id, method))
            eval_entry = summary["methods"][method]
            row = {
                "arm": arm_id, "method": method,
                "bc_from_expected": eval_entry["bc_from_expected"],
                "in_range_fraction": eval_entry["in_range_fraction"],
                "inflation_score": eval_entry["inflation_score"],
            }
            row.update({f"pct_{k}": v for k, v in eval_entry["mean_kingdom_pct"].items()})
            kingdom_rows.append(row)

    # Reproduction gate: baseline arm must match the accepted run exactly.
    reproduction = {}
    for method in METHODS:
        mine = pd.read_csv(
            arms_dir / "baseline_PC" / f"phylum_composition_baseline_PC_{method}.csv", index_col=0
        )
        theirs = pd.read_csv(
            BASELINE_RUN / "strict16_PC_primary" / f"phylum_composition_strict16_PC_primary_{method}.csv",
            index_col=0,
        )
        delta = float((mine[theirs.columns] - theirs).abs().to_numpy().max())
        reproduction[method] = delta
    max_delta = max(reproduction.values())
    reproduction_pass = max_delta < 1e-12
    print(f"\nREPRODUCTION: max |delta| = {max_delta:.3e} -> {'PASS' if reproduction_pass else 'FAIL'}")
    (staging / "reproduction_check.json").write_text(
        json.dumps({"per_method_max_abs_delta": reproduction,
                    "max_abs_delta": max_delta, "pass": reproduction_pass}, indent=2) + "\n",
        encoding="utf-8",
    )
    if not reproduction_pass:
        print("Stopping: rule arms are not attributable to the rules alone.")
        return 1

    effects = pd.concat(effects_frames, ignore_index=True)
    effects.to_csv(staging / "effects_all_arms_methods.csv", index=False)
    comparison = pd.DataFrame(kingdom_rows)
    comparison.to_csv(staging / "method_arm_comparison.csv", index=False)

    significance = (
        effects.groupby(["arm", "method"], sort=False)
        .agg(
            n_units_tested=("unit", "size"),
            drought_p_lt_0_05=("drought_p", lambda v: int((v < 0.05).sum())),
            drought_q_lt_0_05=("drought_q_bh", lambda v: int((v < 0.05).sum())),
            climate_p_lt_0_05=("climate_p", lambda v: int((v < 0.05).sum())),
            climate_q_lt_0_05=("climate_q_bh", lambda v: int((v < 0.05).sum())),
            min_drought_q=("drought_q_bh", "min"),
            min_climate_q=("climate_q_bh", "min"),
        )
        .reset_index()
    )
    significance.to_csv(staging / "significance_summary.csv", index=False)

    # Leverage diagnostics: baseline vs primary rule arm.
    for arm_id, anchor, rules in ARM_SPECS:
        if arm_id not in ("baseline_PC", PRIMARY_ARM):
            continue
        config = make_config(legacy, anchor)
        with ArmPatches(legacy, rules):
            table = leverage_table(legacy, config, unit_groups)
        table.to_csv(staging / f"top_feature_leverage_{arm_id}.csv", index=False)

    (staging / "run_summaries.json").write_text(
        json.dumps(summaries, indent=2, default=str) + "\n", encoding="utf-8"
    )
    meta = {
        "taxonomy_release": RELEASE_ID,
        "baseline_run": str(BASELINE_RUN),
        "dominance_cap": DOMINANCE_CAP,
        "rie_floor": RIE_FLOOR,
        "rie_ceiling": RIE_CEILING,
        "arms": [
            {"arm": a, "anchor": anchor, "rules": sorted(rules)}
            for a, anchor, rules in ARM_SPECS
        ],
        "primary_arm": PRIMARY_ARM,
        "note": (
            "Rule arms are review-only. Rules are general and direction-blind; "
            "expected_kingdom_composition.csv is a held-out development diagnostic, "
            "not an optimization target."
        ),
    }
    (staging / "RUN_META.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")

    staging.rename(OUTPUT)
    print("\n=== method x arm comparison (kingdom means and plausibility) ===")
    print(comparison.round(3).to_string(index=False))
    print(f"\nOutput: {OUTPUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
