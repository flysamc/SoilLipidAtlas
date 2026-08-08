#!/usr/bin/env python3
"""Annotation pipeline Step 8 (ArchLips) for taxonomy release ncbi-phylum-2026-08-04-v1.

Supplementary Method 3 Step 8 / Supplementary Method 7, Round 1 (targeted):
search the archaeal-assigned strict biomarkers against the ArchLips spectral
libraries (Zheng et al., 2026).

The two polarities use DIFFERENT procedures, because the ArchLips library holds
positive-mode reference spectra only ([M+H]+ and [M+NH4]+, verified: 29,446 each
in the high-confidence library, zero negative entries).

POS  direct spectral matching, per analysis-17/positive/scripts/archlips_indval_search.py
     precursor +/-10 ppm, modified cosine at 0.02 Da, require cosine >= 0.3 and
     >= 2 matched peaks; Gold cos>=0.7 & >=4 peaks, Silver cos>=0.5 & >=3,
     Bronze cos>=0.3 & >=2.

NEG  no cosine. Per analysis-16/.../07_archlips/archlips_search_neg.py, the neutral
     MW is taken from the library and negative adducts are COMPUTED
     ([M-H]- = MW-1.00728, [M+HCOO]- = MW+44.99820, [M+CH3COO]- = MW+59.01385),
     matched at +/-10 ppm, then validated against archaeal diagnostic ions and
     neutral losses. Scoring a NEG spectrum against an [M+H]+ reference would be
     chemically meaningless, so it is deliberately not done.

Memory note: the Full library is 1.17 GB. It is streamed, and a compound is kept
only if it can match some query precursor. Results are identical to a full load.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from bisect import bisect_left, bisect_right
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

RELEASE = "ncbi-phylum-2026-08-04-v1"
ROOT = Path(__file__).resolve().parents[2]
RECOVERY = (ROOT / "external" / "SOILMASS_PRODUCER_RECOVERY_2026-08-04" / "payload"
            / "core" / "workspace" / "analysis")

ARCHLIPS_HIGH = RECOVERY / "Archlips_High_confidence_spectral_library.msp"
ARCHLIPS_FULL = RECOVERY / "Archlips_Full_spectral_library.msp"
POS_PRODUCER = RECOVERY / "analysis-17" / "positive" / "scripts" / "archlips_indval_search.py"
NEG_PRODUCER = (RECOVERY / "analysis-16" / "negative_mode" / "04_biomarker_discovery"
                / "07_archlips" / "archlips_search_neg.py")

REL = ROOT / "outputs" / "analysis" / RELEASE
POS_ATLAS = REL / "biomarker_discovery" / "atlas_pos_strict.csv"
NEG_ATLAS = REL / "biomarker_discovery_neg" / "strict_atlas_NEG.csv"
POS_MGF = (REL / "biomarker_discovery" / "external_annotation_package"
           / "figure2a_strict_atlas_with_usable_ms2.mgf")
NEG_MGF = REL / "annotation_recovery_neg" / "spectra" / "strict_full_usable_ms2_NEG.mgf"
OUT_DIR = REL / "annotation" / "step8_archlips"

PPM_TOL = 10.0
FRAG_TOL = 0.02
COSINE_MIN = 0.3
NOISE_THRESHOLD = 0.05

PROTON = 1.00728
FORMATE = 44.99820
ACETATE = 59.01385

ARCHAEAL_DIAG_IONS = [
    ("isoprenoid_C5", 69.070, 0.02), ("isoprenoid_C10", 137.133, 0.02),
    ("glycerol_ether", 89.060, 0.02), ("phytanyl_chain", 339.362, 0.03),
    ("PG_head_153", 152.996, 0.02), ("PE_head_140", 140.012, 0.02),
    ("PI_head_241", 241.012, 0.02), ("phosphate_79", 78.959, 0.02),
    ("phosphate_97", 96.970, 0.02), ("GDGT_biphytanyl", 561.550, 0.04),
    ("menaquinone_187", 187.076, 0.02), ("menaquinone_225", 225.091, 0.02),
]
ARCHAEAL_NLS = [
    ("NL_phytanol", 296.344, 0.04), ("NL_hexose", 162.053, 0.03),
    ("NL_dihexose", 324.106, 0.04), ("NL_serine", 87.032, 0.03),
    ("NL_PE_head", 141.019, 0.03), ("NL_PI_head", 260.053, 0.04),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_mgf(path: Path, wanted: set[str]) -> dict[str, dict]:
    """Read only the spectra whose FEATURE_ID is in `wanted`."""
    out: dict[str, dict] = {}
    fid = None
    pmz = 0.0
    peaks: list[list[float]] = []
    inside = False
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line == "BEGIN IONS":
                inside, fid, pmz, peaks = True, None, 0.0, []
            elif line == "END IONS":
                if inside and fid in wanted and peaks:
                    out[fid] = {"precursor_mz": pmz, "peaks": peaks}
                inside = False
            elif inside:
                if line.startswith("FEATURE_ID="):
                    fid = line[11:]
                elif line.startswith("PEPMASS="):
                    try:
                        pmz = float(line[8:].split()[0])
                    except (ValueError, IndexError):
                        pmz = 0.0
                elif line and (line[0].isdigit() or line[0] == "-"):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            peaks.append([float(parts[0]), float(parts[1])])
                        except ValueError:
                            pass
    return out


def within(sorted_vals: list[float], lo: float, hi: float) -> bool:
    i = bisect_left(sorted_vals, lo)
    return i < len(sorted_vals) and sorted_vals[i] <= hi


def stream_library(path: Path, pos_mz: list[float], neg_mz: list[float]) -> list[dict]:
    """Stream an MSP, keeping only compounds able to match a POS or NEG query.

    Peaks are retained only for POS-matchable compounds; NEG never uses them.
    """
    kept: list[dict] = []
    cur: dict | None = None
    n_seen = 0

    def flush(compound: dict | None) -> None:
        nonlocal kept
        if not compound or "precursor_mz" not in compound or "mw" not in compound:
            return
        pmz = compound["precursor_mz"]
        mw = compound["mw"]
        tol = pmz * PPM_TOL / 1e6
        pos_hit = within(pos_mz, pmz - tol, pmz + tol) if pos_mz else False
        neg_hit = False
        if neg_mz:
            for target in (mw - PROTON, mw + FORMATE, mw + ACETATE):
                t = target * PPM_TOL / 1e6
                if within(neg_mz, target - t, target + t):
                    neg_hit = True
                    break
        if not (pos_hit or neg_hit):
            return
        compound["pos_candidate"] = pos_hit
        compound["neg_candidate"] = neg_hit
        if not pos_hit:
            compound["peaks"] = []
        kept.append(compound)

    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith("Name: "):
                flush(cur)
                n_seen += 1
                cur = {"name": line[6:], "peaks": []}
            elif cur is not None:
                if line.startswith("MW: "):
                    try:
                        cur["mw"] = float(line[4:])
                    except ValueError:
                        pass
                elif line.startswith("PRECURSORMZ: "):
                    try:
                        cur["precursor_mz"] = float(line[13:])
                    except ValueError:
                        pass
                elif line.startswith("Formula: "):
                    cur["formula"] = line[9:]
                elif line.startswith("PRECURSORTYPE: "):
                    cur["adduct"] = line[15:]
                elif line.startswith("SMILES: "):
                    cur["smiles"] = line[8:]
                elif line.startswith("INCHIKEY: "):
                    cur["inchikey"] = line[10:]
                elif line.startswith("Confidencelevel: "):
                    cur["confidence"] = line[17:]
                elif line and (line[0].isdigit() or line[0] == "-"):
                    parts = line.split()
                    if len(parts) >= 2:
                        try:
                            cur["peaks"].append([float(parts[0]), float(parts[1])])
                        except ValueError:
                            pass
    flush(cur)
    print(f"    {path.name}: {n_seen:,} entries scanned, {len(kept):,} retained")
    return kept


def cosine_similarity(peaks1, peaks2, tol=FRAG_TOL) -> tuple[float, int]:
    """Modified cosine, verbatim from the recovered POS producer."""
    if not peaks1 or not peaks2:
        return 0.0, 0
    max1 = max(p[1] for p in peaks1)
    max2 = max(p[1] for p in peaks2)
    if max1 == 0 or max2 == 0:
        return 0.0, 0
    norm1 = [(mz, i / max1) for mz, i in peaks1]
    norm2 = [(mz, i / max2) for mz, i in peaks2]
    used: set[int] = set()
    matched: list[tuple[float, float]] = []
    for mz1, i1 in norm1:
        best_j, best_diff = None, tol + 1
        for j, (mz2, _) in enumerate(norm2):
            if j in used:
                continue
            d = abs(mz1 - mz2)
            if d <= tol and d < best_diff:
                best_diff, best_j = d, j
        if best_j is not None:
            matched.append((i1, norm2[best_j][1]))
            used.add(best_j)
    if not matched:
        return 0.0, 0
    dot = sum(a * b for a, b in matched)
    na = math.sqrt(sum(a * a for a, _ in matched))
    nb = math.sqrt(sum(b * b for _, b in matched))
    if na == 0 or nb == 0:
        return 0.0, 0
    return dot / (na * nb), len(matched)


def pos_tier(cos: float, mp: int) -> str:
    if cos >= 0.7 and mp >= 4:
        return "Gold"
    if cos >= 0.5 and mp >= 3:
        return "Silver"
    if cos >= 0.3 and mp >= 2:
        return "Bronze"
    return "Unvalidated"


def scan_archaeal_ms2(precursor_mz: float, peaks) -> tuple[int, int, str, str]:
    """Diagnostic ion / neutral loss scan, verbatim from the recovered NEG producer."""
    if not peaks:
        return 0, 0, "", ""
    max_int = max(p[1] for p in peaks)
    if max_int <= 0:
        return 0, 0, "", ""
    clean = [(mz, i) for mz, i in peaks if i >= max_int * NOISE_THRESHOLD]
    if not clean:
        return 0, 0, "", ""
    mzs = np.array([p[0] for p in clean])
    ints = np.array([p[1] for p in clean])

    ion_hits = []
    for name, target, tol in ARCHAEAL_DIAG_IONS:
        mask = np.abs(mzs - target) <= tol
        if mask.any():
            rel = ints[mask][int(np.argmax(ints[mask]))] / max_int * 100
            ion_hits.append(f"{name} ({rel:.0f}%)")
    nl_hits = []
    for name, mass, tol in ARCHAEAL_NLS:
        target = precursor_mz - mass
        if target < 50:
            continue
        mask = np.abs(mzs - target) <= tol
        if mask.any():
            rel = ints[mask][int(np.argmax(ints[mask]))] / max_int * 100
            nl_hits.append(f"{name} ({rel:.0f}%)")
    return len(ion_hits), len(nl_hits), "; ".join(ion_hits), "; ".join(nl_hits)


def neg_tier(n_ions: int, n_nls: int, error: float, det_rate: float, conf: str) -> str:
    """Tier ladder, verbatim from the recovered NEG producer."""
    total = n_ions + n_nls
    error = abs(error)
    if error <= 5 and total >= 3 and det_rate >= 0.2:
        return "Gold"
    if error <= 5 and total >= 2:
        return "Gold"
    if error <= 7 and total >= 2:
        return "Silver"
    if error <= 10 and total >= 1 and det_rate >= 0.1:
        return "Silver"
    if error <= 10 and total >= 1:
        return "Bronze"
    if error <= 5 and det_rate >= 0.3:
        return "Bronze"
    if error <= 5 and str(conf) in ("1", "2"):
        return "Bronze"
    return "Unvalidated"


def archaeal_subset(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, low_memory=False)
    arc = df[df["kingdom"].astype(str).str.contains("Archae", case=False, na=False)].copy()
    arc["feature_id"] = arc["feature_id"].astype(str)
    return arc


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"=== Step 8 ArchLips, release {RELEASE} ===")

    pos_arc = archaeal_subset(POS_ATLAS)
    neg_arc = archaeal_subset(NEG_ATLAS)
    print(f"POS archaeal features: {len(pos_arc):,}  {pos_arc.phylum.value_counts().to_dict()}")
    print(f"NEG archaeal features: {len(neg_arc):,}  {neg_arc.phylum.value_counts().to_dict()}")

    print("\nreading query spectra...")
    pos_spec = parse_mgf(POS_MGF, set(pos_arc.feature_id))
    neg_spec = parse_mgf(NEG_MGF, set(neg_arc.feature_id))
    print(f"  POS spectra found: {len(pos_spec):,} / {len(pos_arc):,}")
    print(f"  NEG spectra found: {len(neg_spec):,} / {len(neg_arc):,}")

    pos_mz = sorted(s["precursor_mz"] for s in pos_spec.values())
    neg_mz = sorted(s["precursor_mz"] for s in neg_spec.values())

    print("\nstreaming ArchLips libraries (retaining only matchable compounds)...")
    compounds: list[dict] = []
    for lib in (ARCHLIPS_HIGH, ARCHLIPS_FULL):
        compounds.extend(stream_library(lib, pos_mz, neg_mz))
    print(f"  total retained: {len(compounds):,}")

    # ------------------------------------------------------------------ POS
    print("\n--- POS: spectral cosine matching ---")
    pos_pool = [c for c in compounds if c.get("pos_candidate") and c["peaks"]]
    pos_pool.sort(key=lambda c: c["precursor_mz"])
    pool_mz = [c["precursor_mz"] for c in pos_pool]
    print(f"  POS-matchable reference spectra: {len(pos_pool):,}")

    pos_rows = []
    for fid, spec in pos_spec.items():
        pmz, peaks = spec["precursor_mz"], spec["peaks"]
        tol = pmz * PPM_TOL / 1e6
        lo = bisect_left(pool_mz, pmz - tol)
        hi = bisect_right(pool_mz, pmz + tol)
        best = None
        for i in range(lo, hi):
            comp = pos_pool[i]
            cos, mp = cosine_similarity(peaks, comp["peaks"])
            if cos >= COSINE_MIN and mp >= 2 and (best is None or cos > best["cosine"]):
                best = {
                    "feature_id": fid, "query_mz": pmz,
                    "archlips_name": comp["name"], "archlips_mz": comp["precursor_mz"],
                    "archlips_formula": comp.get("formula", ""),
                    "archlips_adduct": comp.get("adduct", ""),
                    "archlips_smiles": comp.get("smiles", ""),
                    "archlips_confidence": comp.get("confidence", ""),
                    "mass_error_ppm": (pmz - comp["precursor_mz"]) / comp["precursor_mz"] * 1e6,
                    "cosine": cos, "matched_peaks": mp,
                    "query_n_peaks": len(peaks), "ref_n_peaks": len(comp["peaks"]),
                }
        if best:
            best["archlips_tier"] = pos_tier(best["cosine"], best["matched_peaks"])
            pos_rows.append(best)

    pos_df = pd.DataFrame(pos_rows)
    if not pos_df.empty:
        pos_df = pos_df.merge(pos_arc[["feature_id", "phylum", "kingdom"]],
                              on="feature_id", how="left")
    print(f"  matched features: {len(pos_df):,}")
    if not pos_df.empty:
        print(f"  tiers: {pos_df.archlips_tier.value_counts().to_dict()}")

    # ------------------------------------------------------------------ NEG
    print("\n--- NEG: computed adducts + diagnostic-ion validation ---")
    neg_pool = [c for c in compounds if c.get("neg_candidate")]
    print(f"  NEG-matchable reference compounds: {len(neg_pool):,}")

    det = {}
    if "detection_rate" in neg_arc.columns:
        det = dict(zip(neg_arc.feature_id, neg_arc.detection_rate.fillna(0)))

    adducts = [("[M-H]-", -PROTON), ("[M+HCOO]-", FORMATE), ("[M+CH3COO]-", ACETATE)]
    neg_rows = []
    for fid, spec in neg_spec.items():
        pmz, peaks = spec["precursor_mz"], spec["peaks"]
        n_ions, n_nls, ion_str, nl_str = scan_archaeal_ms2(pmz, peaks)
        for comp in neg_pool:
            mw = comp["mw"]
            for adduct_name, delta in adducts:
                target = mw + delta
                if target <= 0:
                    continue
                error = (pmz - target) / target * 1e6
                if abs(error) > PPM_TOL:
                    continue
                tier = neg_tier(n_ions, n_nls, error, float(det.get(fid, 0) or 0),
                                comp.get("confidence", ""))
                neg_rows.append({
                    "feature_id": fid, "consensus_mz": pmz,
                    "archlips_name": comp["name"], "archlips_mw": mw,
                    "archlips_formula": comp.get("formula", ""),
                    "archlips_confidence": comp.get("confidence", ""),
                    "archlips_smiles": comp.get("smiles", ""),
                    "adduct": adduct_name, "target_mz": target,
                    "mass_error_ppm": error,
                    "n_arch_ions": n_ions, "n_arch_nls": n_nls,
                    "arch_ions": ion_str, "arch_nls": nl_str,
                    "arch_detection_rate": float(det.get(fid, 0) or 0),
                    "archlips_tier": tier,
                })

    neg_df = pd.DataFrame(neg_rows)
    if not neg_df.empty:
        neg_df = neg_df.merge(neg_arc[["feature_id", "phylum", "kingdom"]],
                              on="feature_id", how="left")
    print(f"  mass matches: {len(neg_df):,} across "
          f"{neg_df.feature_id.nunique() if not neg_df.empty else 0} features")
    if not neg_df.empty:
        print(f"  tiers: {neg_df.archlips_tier.value_counts().to_dict()}")

    # ------------------------------------------------------------------ write
    pos_path = OUT_DIR / "archlips_pos_strict_matches.csv"
    neg_path = OUT_DIR / "archlips_neg_strict_matches.csv"
    pos_df.to_csv(pos_path, index=False)
    neg_df.to_csv(neg_path, index=False)

    def summarize(df: pd.DataFrame, total: int, mode: str) -> dict:
        validated = df[df.archlips_tier.isin(["Gold", "Silver", "Bronze"])] if not df.empty else df
        feats = validated.feature_id.nunique() if not validated.empty else 0
        gs = validated[validated.archlips_tier.isin(["Gold", "Silver"])] if not validated.empty else validated
        return {
            "mode": mode, "archaeal_strict_features": total,
            "features_with_match": int(df.feature_id.nunique()) if not df.empty else 0,
            "features_validated": int(feats),
            "validated_pct": round(100 * feats / max(total, 1), 2),
            "gold_silver_features": int(gs.feature_id.nunique()) if not gs.empty else 0,
            "gold_silver_pct": round(100 * (gs.feature_id.nunique() if not gs.empty else 0)
                                     / max(total, 1), 2),
        }

    summary = [summarize(pos_df, len(pos_arc), "POS"), summarize(neg_df, len(neg_arc), "NEG")]
    pd.DataFrame(summary).to_csv(OUT_DIR / "step8_coverage_summary.csv", index=False)

    manifest = {
        "taxonomy_release": RELEASE,
        "step": "Supplementary Method 3 Step 8 / Method 7 Round 1 (targeted)",
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "producer": str(Path(__file__).resolve()),
        "producer_sha256": sha256(Path(__file__).resolve()),
        "semantic_sources": {
            "POS": {"path": str(POS_PRODUCER), "sha256": sha256(POS_PRODUCER)},
            "NEG": {"path": str(NEG_PRODUCER), "sha256": sha256(NEG_PRODUCER)},
        },
        "libraries": {
            "high_confidence": {"path": str(ARCHLIPS_HIGH), "sha256": sha256(ARCHLIPS_HIGH)},
            "full": {"path": str(ARCHLIPS_FULL), "sha256": sha256(ARCHLIPS_FULL)},
            "polarity_note": ("ArchLips holds positive-mode reference spectra only "
                              "([M+H]+ / [M+NH4]+). NEG therefore uses computed negative "
                              "adducts plus diagnostic-ion validation, never cosine."),
        },
        "parameters": {
            "ppm_tolerance": PPM_TOL, "fragment_tolerance_da": FRAG_TOL,
            "cosine_min_pos": COSINE_MIN, "noise_threshold_fraction": NOISE_THRESHOLD,
            "neg_adducts": {"[M-H]-": -PROTON, "[M+HCOO]-": FORMATE, "[M+CH3COO]-": ACETATE},
            "pos_tiers": "Gold cos>=0.7 & >=4 peaks; Silver cos>=0.5 & >=3; Bronze cos>=0.3 & >=2",
            "neg_tiers": "ladder on |ppm|, diagnostic-ion+NL count, archaeal detection rate",
        },
        "summary": summary,
        "outputs": {},
    }
    for path in sorted(OUT_DIR.glob("*.csv")):
        manifest["outputs"][path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    (OUT_DIR / "STEP8_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    print("\n=== summary ===")
    print(pd.DataFrame(summary).to_string(index=False))
    print(f"\nwrote {OUT_DIR}")


if __name__ == "__main__":
    sys.exit(main())
