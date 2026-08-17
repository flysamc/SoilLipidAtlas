#!/usr/bin/env python3
"""Build the composite-scoring fastMASST query package (Figure 2c completion).

The 2026-08 strict fastMASST run covered only the IndVal consensus selection
(see fastmasst_all_current_index/manifest.json, scope field). Figure 2c
compares soil detection between the cross-batch consensus selection and the
composite-scoring selection, so the composite features must be searched under
the same current index. This builder stages exactly the atlas features with
discovery_method == 'composite' whose spectra are present in the strict
usable-MS2 MGF. It does not submit anything.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

TAXONOMY_RELEASE = "ncbi-phylum-2026-08-04-v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", required=True, type=Path)
    parser.add_argument("--full-mgf", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    atlas = pd.read_csv(args.atlas, low_memory=False)
    composite = atlas[atlas["discovery_method"] == "composite"]
    composite_ids = set(composite["feature_id"].astype(str))

    kept_blocks: list[str] = []
    kept_ids: set[str] = set()
    current: list[str] = []
    current_id = None
    with args.full_mgf.open("r", encoding="utf-8", errors="strict") as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if stripped == "BEGIN IONS":
                current = [stripped]
                current_id = None
            elif stripped == "END IONS":
                current.append(stripped)
                if current_id in composite_ids:
                    kept_blocks.append("\n".join(current))
                    kept_ids.add(current_id)
                current = []
            elif current:
                current.append(stripped)
                if stripped.startswith("FEATURE_ID="):
                    current_id = stripped.split("=", 1)[1]

    out_mgf = args.output_dir / "composite_strict_fastmasst_pending.mgf"
    out_mgf.write_text("\n".join(kept_blocks) + "\n", encoding="utf-8")

    missing = sorted(composite_ids - kept_ids)
    listing = composite[["feature_id", "phylum", "biomarker_tier"]].copy()
    listing["in_queue"] = listing["feature_id"].astype(str).isin(kept_ids)
    listing.to_csv(args.output_dir / "composite_strict_fastmasst_pending.csv", index=False)

    manifest = {
        "schema_version": 1,
        "taxonomy_release": TAXONOMY_RELEASE,
        "scope": "All strict composite-scoring features with usable MS2, staged for one current fastMASST index (Figure 2c completion)",
        "built_utc": datetime.now(timezone.utc).isoformat(),
        "atlas": {"path": str(args.atlas), "sha256": sha256(args.atlas)},
        "full_mgf": {"path": str(args.full_mgf), "sha256": sha256(args.full_mgf)},
        "composite_features_total": int(len(composite)),
        "staged_with_usable_ms2": len(kept_ids),
        "not_in_usable_ms2_mgf": len(missing),
        "submitted": False,
        "outputs": {
            out_mgf.name: {"bytes": out_mgf.stat().st_size, "sha256": sha256(out_mgf)},
        },
    }
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps({k: v for k, v in manifest.items() if k not in ("atlas", "full_mgf", "outputs")}, indent=1))


if __name__ == "__main__":
    main()
