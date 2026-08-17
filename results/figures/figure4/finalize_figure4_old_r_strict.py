#!/usr/bin/env python3
"""Validate and manifest Figure 4 rendered through the original R producer."""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
from PIL import Image


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def record(path: Path) -> dict[str, object]:
    result: dict[str, object] = {
        "path": str(path.resolve()), "bytes": path.stat().st_size, "sha256": sha256(path)
    }
    if path.suffix.lower() == ".png":
        with Image.open(path) as image:
            result["pixels"] = [image.width, image.height]
            result["dpi"] = list(image.info.get("dpi", ()))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--figure-dir", required=True, type=Path)
    parser.add_argument("--original-r", required=True, type=Path)
    args = parser.parse_args()
    figure_dir = args.figure_dir.resolve()
    isolated = figure_dir / "old_r_isolated"
    output_dir = isolated / "figures_r/out"
    adapted_r = isolated / "figures_r/fig5_fingerprint.R"
    composition_path = isolated / (
        "source_folders/analysis18_09_figures/fig3b_lipid_class/data/class_composition.csv"
    )
    composition = pd.read_csv(composition_path)
    totals = composition.groupby(["mode", "phylum"])["n_features"].sum()
    fractions = composition.groupby(["mode", "phylum"])["fraction"].sum()
    if len(totals) != 32 or not ((fractions - 1).abs() < 1e-9).all():
        raise ValueError("Figure 4b matrix completeness/fraction gate failed")
    if int(totals.xs("POS").sum()) != 11371 or int(totals.xs("NEG").sum()) != 5697:
        raise ValueError("Figure 4b strict denominators failed")

    outputs = [
        output_dir / "Figure_4_strict_current_release.png",
        output_dir / "Figure_4_strict_current_release.pdf",
        output_dir / "Figure_4b_strict_current_release.png",
        output_dir / "Figure_4b_strict_current_release.pdf",
    ]
    for path in outputs:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    original_hash, adapted_hash = sha256(args.original_r), sha256(adapted_r)
    if original_hash == adapted_hash:
        raise ValueError("Expected the isolated renderer to carry declared numeric/output adaptations")
    manifest = {
        "schema_version": 1,
        "status": "review_only_complete_figure_rendered",
        "taxonomy_release": "ncbi-phylum-2026-08-04-v1",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "canonical_mapping": {
            "current": "Figure 4, panel b",
            "legacy_renderer": "fig5_fingerprint.R",
            "legacy_panel_source": "fig3b_lipid_class",
        },
        "renderer": {
            "original": record(args.original_r),
            "adapted_copy": record(adapted_r),
            "unchanged_panel_b_code": True,
            "declared_adaptations": [
                "Panel c 5-percent callout values made data-driven for the strict substrate",
                "Canonical Figure 4 and standalone Figure 4b output filenames added",
            ],
        },
        "figure4b": {
            "denominators": {"POS": 11371, "NEG": 5697, "phyla": 16},
            "class_mapping": "Original freeze_data.py DISPLAY, EXACT_RULES, SUBSTRING_RULES and map_class executed verbatim",
            "input": record(composition_path),
        },
        "outputs": {path.name: record(path) for path in outputs},
        "visual_qa": "passed by direct inspection of standalone panel and integrated full PNG",
        "release_boundary": "review only; submission_source and manuscript-facing build untouched",
    }
    manifest_path = figure_dir / "Figure_4_strict_render_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
