#!/usr/bin/env python3
"""Export the SIRIUS-eligible biomarker features that were never submitted, so a
LISC SIRIUS/CANOPUS run can close the annotation submission gap (see
ANNOTATION_PROGRESS.md 2026-08-12). POS ~2,209 features (bacteria/fungi heavy);
NEG the smaller equivalent.

Eligible-unsubmitted = has usable MS2 AND consensus m/z <= 850 AND not already in
the SIRIUS submission universe (the prior new-features export UNION the historical
atlas SIRIUS cache). Spectra are subset VERBATIM from the strict full-atlas
usable-MS2 MGF so the exact CHARGE/peak formatting SIRIUS already accepts is kept.

Outputs per mode: <mgf>, <features.csv>, and a shared MANIFEST.json with counts
and sha256. A gate asserts every gap feature has a spectrum in the source MGF.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RE = ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1"
OUT = RE / "biomarker_discovery" / "external_annotation_package" / "gapfill_2026-08-12"
MZ_MAX = 850.0

CFG = {
    "POS": {
        "evidence": RE / "annotation/annotation_evidence_pos.csv",
        "atlas_mz": RE / "biomarker_discovery/atlas_pos_strict.csv",
        "source_mgf": RE / "biomarker_discovery/external_annotation_package/figure2a_strict_atlas_with_usable_ms2.mgf",
        "submitted": [
            (RE / "biomarker_discovery/external_annotation_package/figure2a_new_features_sirius_eligible_mz_le_850.csv", "feature_id", ","),
            (ROOT / "external/SOILMASS_PRODUCER_RECOVERY_2026-08-04/payload/core/workspace/Dreams/results/sirius_atlas_full_v2/formula_identifications.tsv", "mappingFeatureId", "\t"),
        ],
    },
    "NEG": {
        "evidence": RE / "annotation/annotation_evidence_neg.csv",
        "atlas_mz": RE / "biomarker_discovery_neg/strict_atlas_NEG.csv",
        "source_mgf": RE / "annotation_recovery_neg/spectra/strict_full_usable_ms2_NEG.mgf",
        "submitted": [
            (RE / "annotation_recovery_neg/sirius_canopus/new_changed_mz_le_850_features.csv", None, ","),
            (ROOT / "external/SOILMASS_PRODUCER_RECOVERY_2026-08-04/payload/core/workspace/Dreams/results/sirius_atlas_neg/formula_identifications.tsv", "mappingFeatureId", "\t"),
        ],
    },
}


def as_bool(s: pd.Series) -> pd.Series:
    if s.dtype == bool:
        return s.fillna(False)
    if s.dtype == object:
        return s.astype(str).str.strip().str.lower().isin(["true", "1", "1.0", "yes"])
    return s.fillna(0).astype(float) > 0


def id_set(path: Path, col: str | None, sep: str) -> set[str]:
    df = pd.read_csv(path, sep=sep, low_memory=False)
    if col is None:  # pick the feature-id-like column
        cands = [c for c in df.columns if "feature" in c.lower()]
        col = cands[0] if cands else df.columns[0]
    return set(df[col].astype(str))


def index_mgf(path: Path) -> dict[str, str]:
    """Return {feature_id: full BEGIN..END block text}."""
    blocks = {}
    cur, fid = [], None
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("BEGIN IONS"):
                cur, fid = [line], None
            elif line.startswith("FEATURE_ID="):
                fid = line.split("=", 1)[1].strip()
                cur.append(line)
            elif line.startswith("END IONS"):
                cur.append(line)
                if fid is None:  # fall back to TITLE
                    for l in cur:
                        if l.startswith("TITLE="):
                            fid = l.split("=", 1)[1].strip()
                            break
                if fid is not None:
                    blocks[fid] = "".join(cur)
                cur, fid = [], None
            else:
                cur.append(line)
    return blocks


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    if OUT.exists():
        sys.exit(f"Refusing to overwrite {OUT} - delete or rename it first.")
    OUT.mkdir(parents=True)
    manifest = {"created": "2026-08-12", "mz_max_inclusive": MZ_MAX, "modes": {}}

    for mode, c in CFG.items():
        ev = pd.read_csv(c["evidence"], low_memory=False)
        mz = pd.read_csv(c["atlas_mz"], low_memory=False)[["feature_id", "consensus_mz"]]
        ev = ev.merge(mz, on="feature_id", how="left")
        submitted = set().union(*[id_set(p, col, sep) for p, col, sep in c["submitted"]])
        ev["submitted"] = ev["feature_id"].astype(str).isin(submitted)
        ev["eligible"] = as_bool(ev["has_usable_ms2"]) & (ev["consensus_mz"] <= MZ_MAX)
        gap = ev[ev["eligible"] & ~ev["submitted"]].copy()
        gap_ids = set(gap["feature_id"].astype(str))
        print(f"[{mode}] biomarkers={len(ev)} submitted={ev.submitted.sum()} "
              f"eligible-unsubmitted={len(gap)}")

        blocks = index_mgf(c["source_mgf"])
        present = [f for f in gap_ids if f in blocks]
        missing = sorted(gap_ids - set(present))
        print(f"[{mode}] source MGF spectra: {len(present)}/{len(gap_ids)} present"
              + (f"  MISSING {len(missing)}: {missing[:5]}..." if missing else ""))
        assert not missing, f"{mode}: {len(missing)} gap features lack a spectrum in {c['source_mgf'].name}"

        # write MGF (verbatim blocks) in a stable feature_id order
        order = gap.sort_values("feature_id")["feature_id"].astype(str).tolist()
        mgf_path = OUT / f"sirius_gapfill_{mode}_mzle850.mgf"
        with mgf_path.open("w", encoding="utf-8") as fh:
            for fid in order:
                fh.write(blocks[fid])
        csv_path = OUT / f"sirius_gapfill_{mode}_mzle850_features.csv"
        gap[["feature_id", "consensus_mz", "phylum", "kingdom"]].sort_values("feature_id").to_csv(csv_path, index=False)

        by_king = gap.groupby("kingdom").size().sort_values(ascending=False).to_dict()
        manifest["modes"][mode] = {
            "eligible_unsubmitted": int(len(gap)),
            "spectra_written": len(order),
            "by_kingdom": {k: int(v) for k, v in by_king.items()},
            "mgf": {"path": str(mgf_path.relative_to(ROOT)).replace("\\", "/"),
                    "sha256": sha256(mgf_path), "bytes": mgf_path.stat().st_size},
            "features_csv": {"path": str(csv_path.relative_to(ROOT)).replace("\\", "/"),
                             "sha256": sha256(csv_path)},
            "source_mgf": str(c["source_mgf"].relative_to(ROOT)).replace("\\", "/"),
        }
        print(f"[{mode}] wrote {mgf_path.name} ({len(order)} spectra) + features.csv; by kingdom: {by_king}\n")

    (OUT / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
