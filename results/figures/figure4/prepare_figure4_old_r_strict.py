#!/usr/bin/env python3
"""Prepare strict Figure 4 inputs for the original, stale-named R renderer."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path

import pandas as pd


RELEASE = "ncbi-phylum-2026-08-04-v1"
EXPECTED = {"POS": 11371, "NEG": 5697}
EXPECTED_ANNOTATED = {"POS": 6360, "NEG": 1889}
SUPERCLASS_TO_LEGACY_ALIAS = {
    "Archaeal lipids": "archaeal_lipid",
    "Betaine lipids": "betaine_lipid",
    "Fatty acyls": "fatty acids",
    "Glycerolipids": "glycerolipid",
    "Glycerophospholipids": "glycerophospholipid",
    "Glycerophospholipids|Sphingolipids": "phosphatidylcholine/sphingomyelin",
    "Prenol lipids": "terpenoid",
    "Sphingolipids": "sphingolipids",
    "Sterol lipids": "sterol lipid",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_author_mapping(freeze_path: Path):
    source = freeze_path.read_text(encoding="utf-8")
    definitions = source[: source.index("# Process POS atlas")]
    namespace = {
        "pd": pd,
        "re": re,
        "json": json,
        "__name__": "_freeze_definitions",
        "__file__": str(freeze_path),
    }
    exec(compile(definitions, str(freeze_path), "exec"), namespace)
    return namespace["map_class"], namespace["DISPLAY"]


def normalize_kingdom(value: str) -> str:
    return {"Viridiplantae": "Plantae", "Protists": "Protozoa"}.get(value, value)


def class_table(path: Path, mode: str, map_class, display: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    evidence = pd.read_csv(path, low_memory=False)
    if len(evidence) != EXPECTED[mode] or evidence["feature_id"].duplicated().any():
        raise ValueError(f"{mode} strict evidence denominator or ID uniqueness failed")
    eligible = evidence["annotation_release_status"].eq("eligible")
    if int(eligible.sum()) != EXPECTED_ANNOTATED[mode]:
        raise ValueError(f"{mode} release-eligible annotation denominator failed")
    buckets = []
    for is_eligible, verbatim, normalized, superclass in zip(
        eligible,
        evidence["annotation_class_verbatim"],
        evidence["annotation_class_normalised"],
        evidence["annotation_superclass"],
    ):
        alias = SUPERCLASS_TO_LEGACY_ALIAS.get(str(superclass), "")
        bucket = map_class(verbatim, normalized, alias) if is_eligible else None
        buckets.append(bucket or "Unclassified")
    evidence = evidence.assign(
        display_class=buckets,
        kingdom=evidence["kingdom"].astype(str).map(normalize_kingdom),
    )
    phylum_meta = evidence[["phylum", "kingdom"]].drop_duplicates()
    if phylum_meta["phylum"].duplicated().any() or len(phylum_meta) != 16:
        raise ValueError(f"{mode} strict phylum mapping is not one-to-one across 16 phyla")
    observed = evidence.groupby(["phylum", "kingdom", "display_class"]).size().rename("n_features")
    rows = []
    for phylum, kingdom in phylum_meta.itertuples(index=False):
        total = int((evidence["phylum"] == phylum).sum())
        for lipid_class in display:
            count = int(observed.get((phylum, kingdom, lipid_class), 0))
            rows.append({
                "mode": mode,
                "phylum": phylum,
                "kingdom": kingdom,
                "class": lipid_class,
                "n_features": count,
                "fraction": count / total,
                "phylum_total": total,
            })
    return pd.DataFrame(rows), phylum_meta


def copy_checked(source: Path, target: Path) -> dict[str, str]:
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    source_hash, target_hash = sha256(source), sha256(target)
    if source_hash != target_hash:
        raise ValueError(f"Copy hash mismatch: {source} -> {target}")
    return {"source": str(source.resolve()), "source_sha256": source_hash,
            "copy": str(target.resolve()), "copy_sha256": target_hash}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    source_root, run_root, output_dir = args.source_root.resolve(), args.run_root.resolve(), args.output_dir.resolve()
    legacy_root = source_root / "manuscript_2_clean/06_figures"
    freeze = legacy_root / "source_folders/analysis18_09_figures/fig3b_lipid_class/data/freeze_data.py"
    original_r = legacy_root / "figures_r/fig5_fingerprint.R"
    original_style = legacy_root / "figures_r/soilmass_style.R"
    map_class, display = load_author_mapping(freeze)

    frames, metadata = [], []
    for mode, name in (("POS", "annotation_evidence_pos.csv"), ("NEG", "annotation_evidence_neg.csv")):
        frame, meta = class_table(run_root / "annotation" / name, mode, map_class, display)
        frames.append(frame)
        metadata.append(meta)
    composition = pd.concat(frames, ignore_index=True)
    if "Cyanobacteriota" in set(composition["phylum"]):
        raise ValueError("Cyanobacteriota leaked into the locked strict phylum set")
    if len(composition) != 16 * 2 * len(display):
        raise ValueError("Figure 4b matrix is incomplete")

    isolated = output_dir / "old_r_isolated"
    data_root = isolated / "source_folders/analysis18_09_figures"
    figures_r = isolated / "figures_r"
    copies = {
        "renderer": copy_checked(original_r, figures_r / "fig5_fingerprint.R"),
        "style": copy_checked(original_style, figures_r / "soilmass_style.R"),
        "simper": copy_checked(run_root / "figure4/panels_ac/simper_curves_strict.csv",
                               data_root / "fig3a_simper_curves/data/simper_curves.csv"),
        "subsampling": copy_checked(run_root / "figure4/panels_ac/subsampling_curve_strict.csv",
                                    data_root / "fig3d_ensemble_stability/data/subsampling_curve.csv"),
    }
    class_dir = data_root / "fig3b_lipid_class/data"
    class_dir.mkdir(parents=True, exist_ok=True)
    class_path = class_dir / "class_composition.csv"
    order_path = class_dir / "class_order.json"
    composition.to_csv(class_path, index=False)
    order_path.write_text(json.dumps(display, indent=2) + "\n", encoding="utf-8")

    pk = pd.concat(metadata, ignore_index=True).drop_duplicates()
    pk = pk.drop_duplicates("phylum").sort_values(["kingdom", "phylum"])
    if len(pk) != 16:
        raise ValueError("Expected 16 unique phylum/kingdom rows")
    pk_path = data_root / "fig1a_dendrogram/data/phylum_kingdom_map.csv"
    pk_path.parent.mkdir(parents=True, exist_ok=True)
    pk.to_csv(pk_path, index=False)

    manifest = {
        "schema_version": 1,
        "taxonomy_release": RELEASE,
        "status": "strict_inputs_ready_for_original_r_renderer",
        "canonical_figure": "Figure 4",
        "canonical_panel": "Figure 4b",
        "legacy_renderer_filename": "fig5_fingerprint.R",
        "legacy_source_folder": "fig3b_lipid_class",
        "renderer_is_byte_identical_to_original": True,
        "copies": copies,
        "generated_inputs": {
            str(path.resolve()): {"bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in (class_path, order_path, pk_path)
        },
        "denominators": {"POS_all": 11371, "NEG_all": 5697,
                         "POS_annotated": 6360, "NEG_annotated": 1889, "phyla": 16},
        "run_command": ["Rscript", str((figures_r / "fig5_fingerprint.R").resolve())],
    }
    manifest_path = output_dir / "Figure_4_old_r_input_manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
