#!/usr/bin/env python3
"""Write the R render folder for the strict-16 Figure 2 (panels a, b, c).

Follows the house pattern of build_suppfig5_v2_render.py: emit the data CSVs in
the schema the submitted renderer already expects, copy the style file, and copy
fig2_atlas.R with lookup-level edits only.  Layout, geometry, panel composition
and the save block are untouched.

Panel a  annotation/tier_counts.csv (POS)   -> data/tier_counts.csv
         The Aug-06 annotation release supersedes the Aug-04 figure2a tier
         split; per-phylum totals are identical (gate below), only the
         Gold/Silver/Bronze/Unidentified partition moved.
Panel b  figure2_panels_20260811/kingdom_sampletype_summary_strict.csv
Panel c  figure2_panels_20260811/shared_vs_exclusive_soil_strict.csv

Declared edits to the submitted renderer and style file:
  1. soilmass_style.R  Plantae -> Viridiplantae, Protozoa -> Protists in
     KINGDOM_COLOURS / KINGDOM_ORDER.  Same Wong hex codes, same order length;
     the strict atlas carries these labels natively and the locked taxonomy
     policy forbids the retired ones.
  2. fig2_atlas.R      panel b fill limits c(0, 40) -> c(0, 60).  The strict
     data maxes at 55.5% (Bacteria x Animal/Clinical); the submitted limit
     would push every tile above 40% out of bounds and render it grey.
  3. fig2_atlas.R      panel c factor levels and bar_colors names follow the
     selection-method rename recorded in build_fig2bc_strict.py.
  4. fig2_atlas.R      header provenance block.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RELEASE_ROOT = PROJECT_ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1"

# Two layouts share one gate implementation so the numbers can never diverge.
#   submitted -- the submitted fig2_atlas.R composition, lookup edits only
#   wide_a    -- author-directed 2026-08-11: panel a full height on the left
#                with a broken count axis, b top-right, c compact beneath
LAYOUTS = {
    "submitted": RELEASE_ROOT / "figure2_strict16_2026-08-11_v1",
    "wide_a": RELEASE_ROOT / "figure2_strict16_2026-08-11_v2_wide_a",
}
V2_LAYOUT_R = Path(__file__).resolve().parent / "fig2_atlas_v2_layout.R"

TIERS_SRC = RELEASE_ROOT / "annotation" / "tier_counts.csv"
FIG2A_TIERS = RELEASE_ROOT / "figure2a" / "figure2a_tier_counts.csv"
PANELS_DIR = RELEASE_ROOT / "figure2_panels_20260811"
PANEL_B_SRC = PANELS_DIR / "kingdom_sampletype_summary_strict.csv"
PANEL_C_SRC = PANELS_DIR / "shared_vs_exclusive_soil_strict.csv"
FIG2BC_MANIFEST = PANELS_DIR / "FIG2BC_MANIFEST.json"
PANREDU_MANIFEST = (
    RELEASE_ROOT
    / "biomarker_discovery"
    / "external_annotation_package"
    / "pan_redu_current_index_pos"
    / "manifest.json"
)

FIGURES_R = PROJECT_ROOT / "manuscript_2_clean" / "06_figures" / "figures_r"
SUBMITTED_R = FIGURES_R / "fig2_atlas.R"
STYLE_R = FIGURES_R / "soilmass_style.R"

EXPECTED_FEATURES = 11371
EXPECTED_PHYLA = 16

# -- lookup edits -----------------------------------------------------------
STYLE_EDITS = [
    ('Plantae  = "#56B4E9"', 'Viridiplantae = "#56B4E9"'),
    ('Protozoa = "#CC79A7"', 'Protists = "#CC79A7"'),
    (
        'KINGDOM_ORDER <- c("Bacteria", "Archaea", "Fungi", "Plantae", "Animalia", "Protozoa")',
        'KINGDOM_ORDER <- c("Bacteria", "Archaea", "Fungi", "Viridiplantae", "Animalia", "Protists")',
    ),
]

R_EDITS = [
    # panel b fill scale: strict data reaches 55.5%
    ("    limits  = c(0, 40),", "    limits  = c(0, 60),"),
    # panel c selection-method labels
    (
        'levels = c("Cross-batch consensus", "Composite atlas"))',
        'levels = c("Indicator Value (IndVal)", "Composite scoring"))',
    ),
    (
        'bar_colors <- c("Cross-batch consensus" = "#2C5F8A", "Composite atlas" = "#B0B0B0")',
        'bar_colors <- c("Indicator Value (IndVal)" = "#2C5F8A", "Composite scoring" = "#B0B0B0")',
    ),
]

HEADER_OLD = """# Fig 2 -- Cross-kingdom lipid atlas overview (3-panel layout)
#   Panel a: Biomarker counts per phylum (POS only) with annotation tiers
#   Panel b: MASST kingdom x sample-type heatmap
#   Panel c: Shared vs exclusive soil detection (45.8% vs 11.2%)"""

HEADER_NEW = """# Fig 2 -- Cross-kingdom lipid atlas overview (3-panel layout)
# Rebuilt under the locked strict release ncbi-phylum-2026-08-04-v1:
# 16 analysis phyla, 11,371 POS biomarkers, 10,867 searched against the
# current metabolomicspanrepo index.
#   Panel a: Biomarker counts per phylum (POS only) with annotation tiers,
#            annotation release 2026-08-06 (6,360 annotated, 55.9%)
#   Panel b: fastMASST kingdom x sample-type heatmap -- REIMPLEMENTATION.
#            The historical sample-type categoriser is unrecovered; this uses
#            the declared sampletype_category_map_v1.csv against a different
#            (current) index, so the submitted panel is NOT reproduced.
#   Panel c: Soil detection by selection method, 149 Pan-ReDU soil datasets.
#            Submitted 45.8% / 11.2% are NOT reproduced: the two selection
#            sets themselves changed size under the strict release.
# Producers: paper2_repro/scripts/build_fig2bc_strict.py (b, c),
#            paper2_repro/annotation_summaries.py (a),
#            paper2_repro/scripts/build_fig2_strict16_render.py (this folder)."""


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def apply_edits(text: str, edits, label: str) -> str:
    for old, new in edits:
        if old not in text:
            sys.exit(f"{label}: anchor not found, refusing to guess -> {old!r}")
        text = text.replace(old, new, 1)
    return text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--layout", choices=sorted(LAYOUTS), default="submitted")
    args = parser.parse_args()
    workspace = LAYOUTS[args.layout]

    for path in (TIERS_SRC, PANEL_B_SRC, PANEL_C_SRC, SUBMITTED_R, STYLE_R):
        if not path.exists():
            sys.exit(f"missing input: {path}")
    if args.layout == "wide_a" and not V2_LAYOUT_R.exists():
        sys.exit(f"missing layout script: {V2_LAYOUT_R}")

    WORKSPACE = workspace
    render = WORKSPACE / "r_render"
    data_dir = render / "data"
    if render.exists():
        shutil.rmtree(render)
    data_dir.mkdir(parents=True)

    # -- panel a ------------------------------------------------------------
    tiers = pd.read_csv(TIERS_SRC)
    tiers = tiers[tiers["mode"] == "POS"].copy()

    gate = {}
    gate["tier_total_pos"] = int(tiers.n_features.sum())
    gate["phyla_pos"] = int(tiers.phylum.nunique())
    if gate["tier_total_pos"] != EXPECTED_FEATURES or gate["phyla_pos"] != EXPECTED_PHYLA:
        sys.exit(
            f"GATE FAILED: POS tier counts {gate['tier_total_pos']} features / "
            f"{gate['phyla_pos']} phyla, expected {EXPECTED_FEATURES}/{EXPECTED_PHYLA}"
        )

    # Cross-producer membership gate: the Aug-04 Figure 2a producer and the
    # Aug-06 annotation release must agree per phylum.  They disagree on the
    # tier split by design (a real annotation change); a per-phylum total
    # mismatch would mean the atlas membership moved and must stop the build.
    fig2a = pd.read_csv(FIG2A_TIERS)
    left = tiers.groupby("phylum").n_features.sum()
    right = fig2a.groupby("phylum").n_features.sum()
    diff = (left - right).fillna(-1)
    if not (diff == 0).all():
        sys.exit(f"GATE FAILED: per-phylum totals differ from figure2a producer:\n{diff[diff != 0]}")
    gate["per_phylum_totals_match_figure2a_producer"] = True
    gate["tier_split_release_2026_08_06"] = (
        tiers.groupby("annotation_tier").n_features.sum().astype(int).to_dict()
    )
    gate["tier_split_figure2a_2026_08_04_superseded"] = (
        fig2a.groupby("annotation_tier").n_features.sum().astype(int).to_dict()
    )
    annotated = int(tiers.loc[tiers.annotation_tier != "Unidentified", "n_features"].sum())
    gate["annotated_pos"] = annotated
    gate["annotated_pct"] = round(100 * annotated / gate["tier_total_pos"], 1)

    # Emitted in the submitted schema: the release column is annotation_tier,
    # the renderer reads `tier`.
    tiers_out = tiers.rename(columns={"annotation_tier": "tier"})[
        ["mode", "phylum", "kingdom", "tier", "n_features"]
    ]
    tiers_out.to_csv(data_dir / "tier_counts.csv", index=False)

    # -- panels b and c -----------------------------------------------------
    panel_b = pd.read_csv(PANEL_B_SRC, dtype=str)
    panel_b.to_csv(data_dir / "kingdom_sampletype_summary.csv", index=False)
    panel_c = pd.read_csv(PANEL_C_SRC)
    panel_c.to_csv(data_dir / "shared_vs_exclusive_soil.csv", index=False)

    kingdoms_a = set(tiers.kingdom)
    kingdoms_b = set(panel_b.kingdom)
    if kingdoms_a != kingdoms_b:
        sys.exit(f"GATE FAILED: panel a kingdoms {sorted(kingdoms_a)} != panel b {sorted(kingdoms_b)}")
    retired = {"Plantae", "Protozoa", "Amoebozoa", "Euryarchaeota", "Crenarchaeota",
               "Bryophyta", "Marchantiophyta", "Tracheophyta", "Trachaeophyta",
               "Halobacteriota", "Bicosoecida"}
    shown = set(tiers.phylum) | kingdoms_a
    if shown & retired:
        sys.exit(f"GATE FAILED: retired labels would be displayed: {sorted(shown & retired)}")
    gate["displayed_labels_policy_clean"] = True
    gate["kingdoms"] = sorted(kingdoms_a)

    # -- renderer + style ---------------------------------------------------
    style = STYLE_R.read_text(encoding="utf-8")
    style = apply_edits(style, STYLE_EDITS, "soilmass_style.R")
    (render / "soilmass_style.R").write_text(style, encoding="utf-8")

    if args.layout == "submitted":
        rscript = SUBMITTED_R.read_text(encoding="utf-8")
        rscript = apply_edits(rscript, [(HEADER_OLD, HEADER_NEW)] + R_EDITS, "fig2_atlas.R")
        renderer_note = "manuscript_2_clean/06_figures/figures_r/fig2_atlas.R (lookup edits only)"
    else:
        # The variant already carries the strict fill limit and the panel c
        # labels; it is a layout change, not a lookup edit, and says so.
        rscript = V2_LAYOUT_R.read_text(encoding="utf-8")
        renderer_note = "paper2_repro/scripts/fig2_atlas_v2_layout.R (author-directed layout change)"
    (render / "fig2_atlas.R").write_text(rscript, encoding="utf-8")

    # -- manifest -----------------------------------------------------------
    bc = json.loads(FIG2BC_MANIFEST.read_text())
    panredu = json.loads(PANREDU_MANIFEST.read_text())
    manifest = {
        "schema_version": 1,
        "taxonomy_release": "ncbi-phylum-2026-08-04-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "producer": "paper2_repro/scripts/build_fig2_strict16_render.py",
        "layout": args.layout,
        "renderer": renderer_note,
        "gates": gate,
        "panel_a": {
            "source": str(TIERS_SRC.relative_to(PROJECT_ROOT)),
            "biomarkers_pos": gate["tier_total_pos"],
            "phyla": gate["phyla_pos"],
            "annotated": annotated,
            "annotated_pct": gate["annotated_pct"],
            "reproduction_status": (
                "corrected build; the submitted Figure 2a producer is only partially "
                "reproducible (figure2a/summary.json status=diagnostic_reproduction_gap) "
                "and the submitted panel used the retired 21-label scheme"
            ),
        },
        "panel_b": {
            "source": str(PANEL_B_SRC.relative_to(PROJECT_ROOT)),
            "searched_features": bc["searched_features"],
            "match_resolution": bc["match_resolution"],
            "reproduction_status": "reimplementation - historical categoriser unrecovered, current index",
        },
        "panel_c": {
            "source": str(PANEL_C_SRC.relative_to(PROJECT_ROOT)),
            "rows": panel_c.to_dict("records"),
            "soil_dataset_ids": panredu["soil_dataset_ids"],
            "soil_file_paths": panredu["soil_file_paths"],
            "submitted_values_not_reproduced": {
                "Cross-batch consensus": {"n_total": 1647, "pct_soil": 45.8},
                "Composite atlas": {"n_total": 8962, "pct_soil": 11.2},
            },
            "reproduction_status": "reimplementation - selection sets changed size under strict release",
        },
        "declared_edits": {
            "soilmass_style.R": [f"{o} -> {n}" for o, n in STYLE_EDITS],
            "fig2_atlas.R": (
                ["header provenance block"] + [f"{o.strip()} -> {n.strip()}" for o, n in R_EDITS]
                if args.layout == "submitted"
                else [
                    "LAYOUT CHANGE (author-directed, not a lookup edit): panel a spans "
                    "the full figure height on the left with a broken count axis "
                    "(0-1,150 | 4,800-5,300 at one shared linear scale); panel b "
                    "top-right; panel c compact beneath it; height 150 -> 160 mm",
                    "carries the same strict fill limit c(0, 60) and panel c labels "
                    "as the submitted-layout render",
                    "theme_nature, Wong palette, tier colours, heatmap ramp, bar "
                    "colours and 5-7 pt text unchanged",
                ]
            ),
        },
        "inputs": {
            str(p.relative_to(PROJECT_ROOT)): sha256(p)
            for p in (TIERS_SRC, PANEL_B_SRC, PANEL_C_SRC, SUBMITTED_R, STYLE_R)
        },
    }
    (WORKSPACE / "FIGURE2_RENDER_MANIFEST.json").write_text(json.dumps(manifest, indent=2))

    print(json.dumps(gate, indent=1))
    print(f"render folder ready: {render}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
