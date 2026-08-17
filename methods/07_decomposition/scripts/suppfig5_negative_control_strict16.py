#!/usr/bin/env python3
"""Supplementary Figure 5 - negative control, rebuilt on the Figure 5 v2
(ArchLips-extended) machinery under the LOCKED strict release
ncbi-phylum-2026-08-04-v1.

Design
------
The negative control feeds pure isolates back through the SAME decomposition
that Figure 5 applies to soil, and asks whether each isolate is returned to its
own organism group.  Every isolate is decomposed leave-one-out: the phylum
centroid reference is rebuilt with that sample removed, so nothing is ever
matched against itself.

Panels (submitted layout of supp_fig5_negative_control.R):
  a  uncorrected pipeline  - raw peak areas, full archaeal reference, SIMPER
                             fold-change weights as published (no rule C)
  b  corrected pipeline v2 - IS normalisation, rule A RIE correction with the
                             0.20 floor / uncalibrated 1.0 fallback, the
                             ArchLips-restricted archaeal reference, and rule C
                             (enriched weight split across k enriched phyla)
  c  per-group self-recovery, uncorrected vs corrected (diagonals of a vs b)

Both panels use the SAME 736-feature substrate and the SAME 164 isolates, so
the a -> b delta is the correction stack alone and nothing else.  This is a
tighter contrast than the legacy figure, which also moved the substrate.

Estimator
---------
fc_weighted_bc is primary, matching the Figure 5 v2 final framing (similarity
family = living-community estimate; bars in Fig 5).  marker_panel (the mixture
benchmark winner, = provenance of matched lipid signal) is computed on the
corrected pipeline as well and written out for the caption and for Table S5.

Status
------
The submitted producer (analysis-19/16_negative_control/build_figure.py) is not
in the repository, so the published panel values are NOT reproducible here; this
is a documented reimplementation, comparable only within this recovered
workflow.  Stage 1 runs the legacy-label configuration and reports its distance
from the published matrices as a code-path diagnostic.

Units are the 16 strict analysis phyla.  Group labels (Bacteria, Archaea, Fungi,
Plantae, Animalia, Protozoa) are the pipeline's display summaries and match
Figure 5 v2's kingdom_composition.csv exactly; per policy, broad ecological
groups are display summaries only and are never analysis units.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RELEASE_ROOT = PROJECT_ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1"
OUTPUT = RELEASE_ROOT / "suppfig5_negative_control_strict16_2026-08-11_v1"
EXTENSION = RELEASE_ROOT / "climgrass" / "strict16_archlips_extended_2026-08-08_v1"
FIG5_V1 = PROJECT_ROOT / "paper2_repro" / "scripts" / "figure5_redesign.py"
FIG5_V2 = PROJECT_ROOT / "paper2_repro" / "scripts" / "figure5_redesign_v2.py"
HANDOFF_ROOT = PROJECT_ROOT / "cg" / "SHADOW_PC_CLIMGRASS_STRICT_PHYLA_HANDOFF_2026-08-06"
TAX_SUMMARY = RELEASE_ROOT / "taxonomy" / "taxonomy_summary.json"
POLICY = PROJECT_ROOT / "paper2_repro" / "config" / "taxonomy_policy.json"
PUBLISHED = (PROJECT_ROOT / "manuscript_2_clean" / "06_figures" / "figures_r"
             / "data" / "supp_negative_control")

ARCHAEA_PHYLA = {"Methanobacteriota", "Thermoproteota", "Halobacteriota"}
GROUP_ORDER = ["Bacteria", "Archaea", "Fungi", "Plantae", "Animalia", "Protozoa"]


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# ---------------------------------------------------------------------------
# Inputs - mirrors figure5_redesign_v2.build_inputs_extended with a toggle for
# the correction stack.  Kept as an explicit copy rather than a refactor so the
# adopted Figure 5 v2 producer is not touched.
# ---------------------------------------------------------------------------
def build_inputs(f5, v2, corrected: bool):
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

    simper_map = pd.read_csv(EXTENSION / "simper_mapping_722.csv")
    arch_map = pd.read_csv(EXTENSION / "archlips_added_mapping.csv")
    cols = ["soil_scan", "atlas_mz", "cosine", "feature_id", "ppm_diff"]
    union = pd.concat([simper_map[cols], arch_map[cols]], ignore_index=True)

    pl = legacy.pipeline
    atlas_df, _soil_df = pl.build_matrices(cache, union)
    if len(atlas_df) != 736:
        raise ValueError(f"Expected 736 union features, got {len(atlas_df)}")

    if corrected:
        config = rf.make_config(legacy, rf.PC_NAME)
        atlas_is = pl.get_atlas_IS_intensity(cache, rf.PC_NAME)
        soil_is = pl.get_soil_IS_intensity(cache, rf.PC_NAME)
        both = pd.concat([atlas_is[atlas_is > 0], soil_is[soil_is > 0]])
        global_median = float(both.median())
        atlas_corr = pl.apply_IS_normalization(
            atlas_df, atlas_is, config.IS_spiked_pmol, global_median)
        feature_ids = list(atlas_corr.index)
        classes, adducts = pl.assign_class_adduct(feature_ids, cache["annot"])
        rule_a = rf.make_rule_a(pl.apply_RIE_correction)
        atlas_corr, _ = rule_a(atlas_corr, classes, adducts, cache["rie"], 1.0,
                               rie_floor=rf.RIE_FLOOR, rie_ceiling=rf.RIE_CEILING)
    else:
        atlas_corr = atlas_df.copy()
        feature_ids = list(atlas_corr.index)

    sample_phylum = {s: p for s, p in cache["sample_phylum"].items()
                     if s in atlas_corr.columns}
    phyla = sorted(set(sample_phylum.values()))

    n_arch_features = None
    if corrected:
        archlips_ids = pl.get_archlips_feature_ids()
        archaea_cols = [s for s, p in sample_phylum.items() if p in ARCHAEA_PHYLA]
        non_archlips = ~atlas_corr.index.astype(str).isin(archlips_ids)
        n_arch_features = int((~non_archlips).sum())
        atlas_corr = atlas_corr.copy()
        atlas_corr.loc[non_archlips, archaea_cols] = 0.0

    simper_tab = cache["simper"]
    if corrected:
        simper_ids = set(simper_map["feature_id"].astype(str))
        median_fc = float(simper_tab.loc[
            (simper_tab["direction"] == "enriched")
            & simper_tab["feature_id"].isin(simper_ids), "fold_change"].median())
        synthetic = pd.DataFrame({
            "feature_id": arch_map["feature_id"].astype(str),
            "phylum": arch_map["archlips_phylum"].astype(str),
            "direction": "enriched",
            "fold_change": median_fc,
        })
        synthetic = synthetic[synthetic["phylum"].isin(phyla)]
        simper_tab = pd.concat([simper_tab, synthetic], ignore_index=True)

    from decomposition import build_direction_masks
    enriched_masks, fc_weights = build_direction_masks(simper_tab, phyla, feature_ids)
    n_enriched = np.zeros(len(feature_ids))
    for p in phyla:
        n_enriched += enriched_masks[p].astype(float)
    if corrected:                                    # rule C
        divisor = np.maximum(n_enriched, 1.0)
        for p in phyla:
            w = fc_weights[p].copy()
            m = enriched_masks[p]
            w[m] = w[m] / divisor[m]
            fc_weights[p] = w

    return dict(
        atlas=atlas_corr, sample_phylum=sample_phylum, phyla=phyla,
        feature_ids=feature_ids, enriched_masks=enriched_masks,
        fc_weights=fc_weights, unit_groups=unit_groups, n_enriched=n_enriched,
        n_archlips_features=n_arch_features,
    )


# ---------------------------------------------------------------------------
# Leave-one-out negative control
# ---------------------------------------------------------------------------
def negative_control(f5, inputs, est_name, sqrt_transform=False):
    atlas = inputs["atlas"]
    phyla = inputs["phyla"]
    groups = inputs["unit_groups"]
    values = atlas.to_numpy(dtype=float)
    col_index = {c: i for i, c in enumerate(atlas.columns)}
    samples = [s for s in atlas.columns if s in inputs["sample_phylum"]]
    cols_by_phylum = {
        p: [col_index[s] for s, pp in inputs["sample_phylum"].items() if pp == p]
        for p in phyla
    }
    ctx = {"phyla": phyla, "fc_weights": inputs["fc_weights"],
           "enriched_masks": inputs["enriched_masks"],
           "n_enriched": inputs["n_enriched"]}

    rows = []
    for s in samples:
        true_phylum = inputs["sample_phylum"][s]
        i = col_index[s]
        ref = f5.build_reference(values, cols_by_phylum, phyla,
                                 exclude={true_phylum: {i}})
        prop = f5.run_estimator(est_name, ref, values[:, i], ctx, sqrt_transform)
        by_group = {g: 0.0 for g in GROUP_ORDER}
        for p, v in zip(phyla, prop):
            by_group[groups[p]] += float(v)
        true_group = groups[true_phylum]
        row = {"sample": s, "true_phylum": true_phylum, "true_group": true_group}
        row.update({g: round(100 * by_group[g], 2) for g in GROUP_ORDER})
        row["correct_group_pct"] = row[true_group]
        row["dominant_group"] = max(GROUP_ORDER, key=lambda g: by_group[g])
        row["dominant_correct"] = int(row["dominant_group"] == true_group)
        rows.append(row)

    per_sample = pd.DataFrame(rows)
    conf = (per_sample.groupby("true_group")[GROUP_ORDER].mean()
            .reindex(GROUP_ORDER).round(1))
    conf.index.name = "true_group"
    return per_sample, conf


def summarise(per_sample, conf, label):
    diag = pd.Series({g: conf.loc[g, g] for g in GROUP_ORDER})
    off_arch = pd.Series({g: conf.loc[g, "Archaea"]
                          for g in GROUP_ORDER if g != "Archaea"})
    dom = 100 * per_sample["dominant_correct"].mean()
    print(f"  [{label}] dominant-group correct {dom:.1f}% "
          f"({int(per_sample['dominant_correct'].sum())}/{len(per_sample)}); "
          f"mean self-recovery {diag.mean():.1f}%; "
          f"max non-archaeal->Archaea {off_arch.max():.1f}%")
    return {"dominant_correct_pct": round(dom, 1),
            "dominant_correct_n": int(per_sample["dominant_correct"].sum()),
            "n_samples": int(len(per_sample)),
            "self_recovery_pct": {g: float(diag[g]) for g in GROUP_ORDER},
            "mean_self_recovery_pct": round(float(diag.mean()), 1),
            "max_nonarchaeal_to_archaea_pct": round(float(off_arch.max()), 1),
            "worst_nonarchaeal_sink_group": str(off_arch.idxmax())}


def stage1_diagnostic(conf_raw, conf_cor):
    """Distance from the published panels. The submitted producer is not in the
    repository, so this cannot be a PASS/FAIL gate - it is reported as-is."""
    out = {}
    for tag, fn, mine in (("uncorrected", "negative_control_kingdom_confusion.csv", conf_raw),
                          ("corrected", "negative_control_corrected_kingdom_confusion.csv", conf_cor)):
        path = PUBLISHED / fn
        if not path.exists():
            out[tag] = "published matrix not found"
            continue
        pub = pd.read_csv(path).set_index("true_kingdom")
        shared = [g for g in GROUP_ORDER if g in pub.index]
        delta = (mine.loc[shared, GROUP_ORDER] - pub.loc[shared, GROUP_ORDER]).abs()
        out[tag] = {"published_rows": list(pub.index),
                    "compared_rows": shared,
                    "max_abs_delta_pct": round(float(delta.to_numpy().max()), 1),
                    "mean_abs_delta_pct": round(float(delta.to_numpy().mean()), 1)}
        print(f"  [{tag}] vs published: max |delta| = "
              f"{out[tag]['max_abs_delta_pct']} pp over rows {shared}")
    return out


def main() -> int:
    if OUTPUT.exists():
        sys.exit(f"Refusing to overwrite {OUTPUT} - delete or rename it first.")

    tax = json.loads(TAX_SUMMARY.read_text())
    units16 = sorted(tax["analysis_phyla"])
    assert len(units16) == 16 and tax["taxonomy_release"] == "ncbi-phylum-2026-08-04-v1"
    policy = json.loads(POLICY.read_text(encoding="utf-8"))
    assert policy["status"] == "locked"

    sys.path.insert(0, str(HANDOFF_ROOT / "source" / "framework"))
    f5 = load_module("figure5_v1", FIG5_V1)
    v2 = load_module("figure5_v2", FIG5_V2)

    print("Building uncorrected pipeline (panel a)...")
    inp_raw = build_inputs(f5, v2, corrected=False)
    print("Building corrected v2 pipeline (panel b)...")
    inp_cor = build_inputs(f5, v2, corrected=True)

    for inp in (inp_raw, inp_cor):
        assert sorted(inp["phyla"]) == units16, "units are not the strict 16"
    retired = {"Amoebozoa", "Euryarchaeota", "Crenarchaeota", "Halobacteriota",
               "Bryophyta", "Marchantiophyta", "Tracheophyta", "Trachaeophyta",
               "Bicosoecida", "Rootnodules", "Mixed"}
    assert not (set(units16) & retired), "retired label present in analysis units"

    print("Leave-one-out negative control...")
    ps_raw, conf_raw = negative_control(f5, inp_raw, "fc_weighted_bc")
    s_raw = summarise(ps_raw, conf_raw, "uncorrected / fc_weighted_bc")
    ps_cor, conf_cor = negative_control(f5, inp_cor, "fc_weighted_bc")
    s_cor = summarise(ps_cor, conf_cor, "corrected v2 / fc_weighted_bc")
    ps_mrk, conf_mrk = negative_control(f5, inp_cor, "marker_panel")
    s_mrk = summarise(ps_mrk, conf_mrk, "corrected v2 / marker_panel")

    print("Stage 1 diagnostic against published panels...")
    diag = stage1_diagnostic(conf_raw, conf_cor)

    recovery = pd.DataFrame({
        "group": GROUP_ORDER,
        "n": [int((ps_cor["true_group"] == g).sum()) for g in GROUP_ORDER],
        "uncorrected_pct": [conf_raw.loc[g, g] for g in GROUP_ORDER],
        "corrected_pct": [conf_cor.loc[g, g] for g in GROUP_ORDER],
        "marker_panel_pct": [conf_mrk.loc[g, g] for g in GROUP_ORDER],
    })
    recovery["delta_pp"] = (recovery["corrected_pct"] - recovery["uncorrected_pct"]).round(1)

    OUTPUT.mkdir(parents=True)
    conf_raw.to_csv(OUTPUT / "negative_control_group_confusion_uncorrected.csv")
    conf_cor.to_csv(OUTPUT / "negative_control_group_confusion_corrected.csv")
    conf_mrk.to_csv(OUTPUT / "negative_control_group_confusion_marker_panel.csv")
    ps_raw.to_csv(OUTPUT / "negative_control_persample_uncorrected.csv", index=False)
    ps_cor.to_csv(OUTPUT / "negative_control_persample_corrected.csv", index=False)
    ps_mrk.to_csv(OUTPUT / "negative_control_persample_marker_panel.csv", index=False)
    recovery.to_csv(OUTPUT / "negative_control_self_recovery.csv", index=False)

    summary = {
        "status": "reimplementation - submitted producer "
                  "(analysis-19/16_negative_control/build_figure.py) not in repository; "
                  "published panel values are not reproducible and are not interchangeable",
        "taxonomy_release": "ncbi-phylum-2026-08-04-v1 (locked strict policy)",
        "supersedes": "manuscript_2_clean/06_figures/figures_r/data/supp_negative_control "
                      "(163 isolates, 18 legacy phyla incl. Amoebozoa/Euryarchaeota/"
                      "Bryophyta/Marchantiophyta/Trachaeophyta)",
        "built_on": "Figure 5 v2 ArchLips machinery "
                    "(figure5_redesign_v2.build_inputs_extended)",
        "substrate": "736 features: 722 SIMPER + 14 spectrally matched archaeol markers",
        "analysis_units": units16,
        "group_labels": GROUP_ORDER,
        "group_label_note": "display summaries only, identical to Figure 5 v2 "
                            "kingdom_composition.csv; never analysis units",
        "n_isolates": int(len(ps_cor)),
        "isolate_count_reconciliation":
            "164 strict isolates vs 163 in the legacy figure: the legacy set used the "
            "2026-08-03 partition and included 2 Cyanobacteriota samples, which the "
            "locked release excludes as below-threshold (n>=2 in both polarities); "
            "164 matches the Supplementary Fig 4 denominator",
        "design": "leave-one-out: phylum centroid reference rebuilt with the test "
                  "isolate removed; same substrate and same isolates in both panels, "
                  "so the a->b delta is the correction stack alone",
        "corrections_in_panel_b": ["IS normalisation",
                                   "rule A RIE correction (0.20 floor, uncalibrated 1.0)",
                                   "ArchLips-restricted archaeal reference",
                                   "rule C (enriched weight split across k enriched phyla)"],
        "n_archlips_diagnostic_features": inp_cor["n_archlips_features"],
        "estimators": {
            "primary": "fc_weighted_bc (Figure 5 v2 bars; living-community estimate)",
            "secondary": "marker_panel (Figure 5 v2 diamonds; mixture-benchmark winner)",
        },
        "uncorrected": s_raw,
        "corrected": s_cor,
        "corrected_marker_panel": s_mrk,
        "published_panel_distance": diag,
    }
    (OUTPUT / "RUN_SUMMARY.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print("\nSelf-recovery by group (%):")
    print(recovery.to_string(index=False))
    print(f"\nWrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
