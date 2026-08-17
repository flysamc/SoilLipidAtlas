#!/usr/bin/env python3
"""Supplementary Figure 6 stage 2 - rerun the authors' NEG producers on the
LOCKED release units (ncbi-phylum-2026-08-04-v1).

Stage 1 reproduced the published panels exactly (max abs delta = 0 for panels
a, b, c), so every difference reported here is attributable to the taxonomy
correction and not to a reimplementation.

What changes, and only this:
  1. the atlas sample->phylum table is relabelled to the release's ncbi_phylum
     (Amoebozoa -> Discosea/Evosea + 1 sample reassigned to Pseudomonadota;
      Euryarchaeota+Methanobacteriota -> Methanobacteriota; Crenarchaeota+
      Thermoproteota -> Thermoproteota; the five land-plant labels ->
      Streptophyta; Mucoromycota -> Mucoromycota + Mortierellomycota);
     samples with no valid NCBI phylum ancestor (Bicosoecida) and out-of-scope
     viral samples get an empty phylum, which the framework's EXCLUDE_PHYLA
     already drops;
  2. PHYLUM_KINGDOM gains the strict units it does not know about
     (Discosea, Evosea, Heterolobosea -> Protozoa; Streptophyta -> Plantae;
      Mortierellomycota -> Fungi). Lookup extension only.

No estimator, threshold or correction logic is touched.

Two variants are produced:
  - noarch : Archaea excluded, five kingdoms - directly comparable to the
             published panel a and to the POS adopted column;
  - witharch : Archaea retained - the release gives NEG Methanobacteriota n=17
             and Thermoproteota n=4, depth the published run discarded wholesale.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[2]
BUNDLE = PROJECT_ROOT / "external" / "NEG_CLIMGRASS_FOR_P2R"
RELEASE_ROOT = PROJECT_ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1"
WS = RELEASE_ROOT / "suppfig6_neg_strict16_2026-08-12_v1"
STAGE2 = WS / "stage2_strict16_simperfix"
# Strict NEG SIMPER atlas rebuilt by suppfig6_rebuild_neg_simper.py (legacy-label
# gate passed at 1.1e-13). The delivered atlas is legacy-keyed and leaves
# Discosea/Evosea/Streptophyta with no fingerprints - see
# STAGE2_BLOCKED_SIMPER_RELABEL.md.
STRICT_SIMPER = (WS / "neg_simper_rebuild" / "stageB_strict16" / "out"
                 / "simper_fingerprint_atlas.csv")
CODE = STAGE2 / "code"
ATLAS = BUNDLE / "03_ATLAS_REFERENCE_NEG"
SOIL = BUNDLE / "01_SOIL_DATA_NEG"
LEGACY = BUNDLE / "02_SOIL_DATA_NEG_LEGACY"
INPUTS = BUNDLE / "04_CORRECTION_INPUTS"
META = BUNDLE / "07_TREATMENT_METADATA"
FBMN_NEG = (PROJECT_ROOT / "external" / "SOILMASS_PRODUCER_RECOVERY_2026-08-04"
            / "payload" / "core" / "workspace" / "analysis" / "FBMN_all_batches_NEG")
TAX = RELEASE_ROOT / "taxonomy"

STRICT_META = STAGE2 / "atlas_sample_metadata_strict16.csv"

NEW_KINGDOMS = {
    "Discosea": "Protozoa", "Evosea": "Protozoa", "Heterolobosea": "Protozoa",
    "Streptophyta": "Plantae", "Mortierellomycota": "Fungi",
}


def q(p: Path) -> str:
    return repr(str(p.resolve()))


def build_strict_metadata() -> pd.DataFrame:
    a = pd.read_csv(ATLAS / "atlas_sample_metadata.csv")
    r = pd.read_csv(TAX / "sample_metadata_NEG_ncbi_phylum.csv")
    if set(a["sample_col"]) != set(r["sample_col"]):
        sys.exit("atlas metadata and release NEG metadata do not cover the same samples")
    m = a.merge(r[["sample_col", "ncbi_phylum", "taxonomy_scope", "ecological_group"]],
                on="sample_col", how="left")
    core = m["taxonomy_scope"] == "core_candidate"
    out = m.copy()
    out["legacy_phylum"] = out["phylum"]
    out["phylum"] = out["ncbi_phylum"].where(core, "")
    out["kingdom"] = out["ecological_group"].where(core, "")
    STAGE2.mkdir(parents=True, exist_ok=True)
    cols = ["sample_col", "sample_name", "kingdom", "phylum", "batch",
            "legacy_phylum", "taxonomy_scope"]
    out[cols].to_csv(STRICT_META, index=False)
    return out


def stage_code() -> None:
    if CODE.exists():
        try:
            shutil.rmtree(CODE)
        except OSError:
            for child in sorted(CODE.rglob("*"), key=lambda p: -len(p.parts)):
                try:
                    child.unlink() if child.is_file() else child.rmdir()
                except OSError:
                    pass
    shutil.copytree(BUNDLE / "05_CODE", CODE, dirs_exist_ok=True)
    (STAGE2 / "15_neg_pos_matched_strategy").mkdir(parents=True, exist_ok=True)
    inputs_dir = STAGE2 / "00_inputs"
    inputs_dir.mkdir(parents=True, exist_ok=True)
    for f in INPUTS.iterdir():
        if f.is_file():
            shutil.copy2(f, inputs_dir / f.name)

    edits = {
        "pipeline_neg_cer.py": [
            ('CONSENSUS_NEG = A16 / "negative_mode" / "03_alignment" / "consensus_aligned_table.csv"',
             f'CONSENSUS_NEG = Path({q(ATLAS / "consensus_aligned_table.csv")})'),
            ('METADATA_NEG  = A16 / "negative_mode" / "04_biomarker_discovery" / "01_composite_scoring" / "sample_metadata.csv"',
             f'METADATA_NEG  = Path({q(STRICT_META)})'),
            ('SIMPER_NEG    = A16 / "negative_mode" / "04_biomarker_discovery" / "13_simper_fingerprint" / "simper_fingerprint_atlas.csv"',
             f'SIMPER_NEG    = Path({q(STRICT_SIMPER)})'),
            ('ANNOT_NEG     = A16 / "negative_mode" / "04_biomarker_discovery" / "03_lipidsearch_annotations" / "consensus_lipidsearch_annotations.csv"',
             f'ANNOT_NEG     = Path({q(ATLAS / "consensus_lipidsearch_annotations.csv")})'),
            ('CLIMGRASS_NEG_QUANT_REPROCESSED = HERE / "output" / "climgrass_neg_reprocess_gnps_quant.csv"',
             f'CLIMGRASS_NEG_QUANT_REPROCESSED = Path({q(SOIL / "climgrass_neg_reprocess_gnps_quant.csv")})'),
            ('CLIMGRASS_NEG_QUANT_OLD = A16 / "negative_mode" / "Mzmime3paper1-negative-with-50" / "without-background" / "oe12-3neg50-without-background_quant.csv"',
             f'CLIMGRASS_NEG_QUANT_OLD = Path({q(LEGACY / "oe12-3neg50-without-background_quant.csv")})'),
            ('CLIMGRASS_NEG_META = ROOT / "analysis-17" / "positive" / "soil_decomposition" / "climgrass-experiment" / "sample_metadata.csv"',
             f'CLIMGRASS_NEG_META = Path({q(META / "pos_sample_metadata_HAS_SAMPLE48.csv")})'),
            ('RIE_TABLE = A19 / "00_inputs" / "rie_table_s10.csv"',
             f'RIE_TABLE = Path({q(INPUTS / "rie_table_s10.csv")})'),
            ('EXPECTED_REF = A19 / "00_inputs" / "expected_kingdom_composition.csv"',
             f'EXPECTED_REF = Path({q(INPUTS / "expected_kingdom_composition.csv")})'),
        ],
        "framework/decomposition.py": [
            ("    'Amoebozoa': 'Protozoa',\n}",
             "    'Amoebozoa': 'Protozoa',\n"
             "    # strict release ncbi-phylum-2026-08-04-v1 units (lookup extension)\n"
             "    'Discosea': 'Protozoa', 'Evosea': 'Protozoa', 'Heterolobosea': 'Protozoa',\n"
             "    'Streptophyta': 'Plantae', 'Mortierellomycota': 'Fungi',\n}"),
        ],
        "neg_pos_matched_decomposition.py": [
            ('A19 = Path("/Users/rahulsamrat/Desktop/Projects/soilmass-analysis/analysis/analysis-19")',
             f'A19 = Path({q(STAGE2)})'),
            ('_meta = pd.read_csv("/Users/rahulsamrat/Desktop/Projects/soilmass-analysis/analysis/analysis-16/"\n'
             '                    "negative_mode/04_biomarker_discovery/01_composite_scoring/sample_metadata.csv")',
             f'_meta = pd.read_csv({q(STRICT_META)})\n'
             '_meta = _meta[_meta["phylum"].astype(str) != ""]'),
        ],
        "neg_ms2_confirmation.py": [
            ('A19 = Path("/Users/rahulsamrat/Desktop/Projects/soilmass-analysis/analysis/analysis-19")',
             f'A19 = Path({q(STAGE2)})'),
            ('NPHY = pd.read_csv(BASE / "analysis-16/negative_mode/04_biomarker_discovery/01_composite_scoring/sample_metadata.csv")["phylum"].value_counts().to_dict()',
             f'NPHY = pd.read_csv(Path({q(STRICT_META)}))["phylum"].value_counts().to_dict()'),
            ('SOIL_MGF = A19 / "climgrass_neg_reprocess" / "output" / "climgrass_neg_reprocess_gnps.mgf"',
             f'SOIL_MGF = Path({q(SOIL / "climgrass_neg_reprocess_gnps.mgf")})'),
            ('ATLAS_MGFS = sorted(BASE.glob("FBMN_all_batches_NEG/batch_0*/batch0*_iimn_gnps.mgf"))',
             f'ATLAS_MGFS = sorted(Path({q(FBMN_NEG)}).glob("batch_0*/batch0*_iimn_gnps*.mgf"))'),
        ],
        "robustness_checks.py": [
            ('A19 = Path("/Users/rahulsamrat/Desktop/Projects/soilmass-analysis/analysis/analysis-19")',
             f'A19 = Path({q(STAGE2)})'),
            ('NPHY = pd.read_csv(A19.parent / "analysis-16/negative_mode/04_biomarker_discovery/01_composite_scoring/sample_metadata.csv")["phylum"].value_counts().to_dict()',
             f'NPHY = pd.read_csv(Path({q(STRICT_META)}))["phylum"].value_counts().to_dict()'),
        ],
    }

    missed = []
    for fname, pairs in edits.items():
        p = CODE / fname
        s = p.read_text(encoding="utf-8")
        for old, new in pairs:
            if old not in s:
                missed.append(f"{fname}: {old.splitlines()[0][:70]}")
                continue
            s = s.replace(old, new, 1)
        p.write_text(s, encoding="utf-8")
    if missed:
        print("EDIT TARGETS NOT FOUND:")
        for m in missed:
            print("  -", m)
        sys.exit(1)


def add_archaea_variant() -> None:
    """Append an Archaea-retained run to the producer (additive only)."""
    p = CODE / "neg_pos_matched_decomposition.py"
    s = p.read_text(encoding="utf-8")
    anchor = 'neg5, np5 = neg_all, len(ph_all)'
    extra = (
        'neg_arch_all, ph_arch_all = run(mapped, drop_archaea=False, '
        'label="neg_witharch_nall", min_n=2)\n'
        'neg_arch_n5, ph_arch_n5 = run(mapped, drop_archaea=False, '
        'label="neg_witharch_n5", min_n=5)\n'
        'print("\\n=== Archaea-retained variant (strict units) ===")\n'
        'print("  n>=2 phyla:", len(ph_arch_all), sorted(ph_arch_all))\n'
        'print(pd.DataFrame({"n>=2": neg_arch_all, "n>=5": neg_arch_n5}).round(2).to_string())\n'
        + anchor
    )
    if anchor not in s:
        sys.exit("anchor for archaea variant not found")
    p.write_text(s.replace(anchor, extra, 1), encoding="utf-8")


def main() -> int:
    print("Building strict-relabelled atlas metadata...")
    m = build_strict_metadata()
    core = m[m["taxonomy_scope"] == "core_candidate"]
    vc = core["phylum"].value_counts()
    tax = json.loads((TAX / "taxonomy_summary.json").read_text())
    strict16 = sorted(tax["analysis_phyla"])
    ge2 = sorted(vc[vc >= 2].index)
    print(f"  core samples {len(core)} / {len(m)}; collection phyla {len(vc)}")
    print(f"  phyla with n>=2: {len(ge2)}")
    if ge2 != strict16:
        sys.exit(f"GATE FAIL: n>=2 set is not the strict 16\n  got {ge2}\n  want {strict16}")
    print("  GATE PASS: n>=2 set == release analysis_phyla (strict 16)")
    dropped = m[m["taxonomy_scope"] != "core_candidate"]
    print(f"  dropped (no phylum-rank ancestor / out of scope): "
          f"{dropped['legacy_phylum'].value_counts().to_dict()}")

    print("Staging code...")
    stage_code()
    add_archaea_variant()

    env_pp = str((CODE / "framework").resolve())
    for script in ("neg_pos_matched_decomposition.py", "robustness_checks.py",
                   "neg_ms2_confirmation.py"):
        print(f"\n=== running {script} ===")
        r = subprocess.run([sys.executable, script], cwd=str(CODE),
                           env={**__import__("os").environ,
                                "PYTHONPATH": env_pp, "PYTHONIOENCODING": "utf-8"},
                           capture_output=True, text=True)
        sys.stdout.write(r.stdout[-4000:])
        if r.returncode != 0:
            sys.stderr.write(r.stderr[-3000:])
            return r.returncode
    print(f"\nStage 2 outputs: {STAGE2 / '15_neg_pos_matched_strategy'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
