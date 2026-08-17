#!/usr/bin/env python3
"""
Supplementary Figure 4 — NNLS leave-one-out under the LOCKED strict release
ncbi-phylum-2026-08-04-v1 (16 analysis phyla; Euryarchaeota+Halobacteriota ->
Methanobacteriota, Crenarchaeota -> Thermoproteota, land plants -> Streptophyta).

Supersedes the 19-unit run in suppfig4_loo_2026-08-11_v1, which followed the
older 2026-08-03 coauthor-package partition. Code path is the same port of
rerun_impact.py validated there (old-label gate 97/161 = 60.2% PASS); this
script re-runs stage 1 as a code-path gate, then the strict-16 run.

Published 52.8% remains not reproducible; all numbers here are the documented
reimplementation and comparable only within that family.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import nnls

ROOT = Path(r"C:\Users\Shadow\Desktop\P2R")
POS_TABLE = ROOT / "analysis/analysis-15/03_alignment/consensus_aligned_table.csv"
POS_META_OLD = ROOT / "analysis/analysis-15/04_biomarker_discovery/01_composite_scoring/sample_metadata.csv"
POS_META_STRICT = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/sample_metadata_POS_ncbi_phylum.csv"
TAX_SUMMARY = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/taxonomy_summary.json"
GATE_OLD = ROOT / "COAUTHOR_PACKAGE_2026-08-03/08_regenerated_figures/data/loo_per_unit_old.csv"
OUT_DIR = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/suppfig4_loo_strict16_2026-08-11_v1"

OLD_UNITS = sorted([
    'Actinomycetota', 'Amoebozoa', 'Arthropoda', 'Ascomycota', 'Bacillota', 'Basidiomycota',
    'Bryophyta', 'Chlorophyta', 'Euryarchaeota', 'Marchantiophyta', 'Methanobacteriota',
    'Mollusca', 'Mucoromycota', 'Nematoda', 'Pseudomonadota', 'Thermoproteota', 'Trachaeophyta'])


def build(table_csv, meta_csv, sample_col, units, tag):
    df = pd.read_csv(table_csv, index_col=0, low_memory=False)
    meta = pd.read_csv(meta_csv)
    sample_cols = [c for c in df.columns if c.startswith('sample:')]
    s2u = {}
    for _, row in meta.iterrows():
        col = row[sample_col]
        if col in df.columns and row.get('phylum') in units:
            s2u[col] = row['phylum']
    mapped = [s for s in sample_cols if s in s2u]

    cross = df[df['n_batches'] >= 2]
    inten = cross[mapped].fillna(0)
    det = (inten > 0).sum(axis=1) / len(mapped)
    mi = inten.mean(axis=1)
    qmask = (det >= 0.05) & (mi >= 500)
    q = inten[qmask]
    print(f"  [{tag}] mapped={len(mapped)} features={int(qmask.sum()):,} units={len(units)}")
    return q, mapped, s2u


def loo_nnls_confusion(q, mapped, s2u, units):
    X = q[mapped].values.astype(float)
    labels = np.array([s2u[s] for s in mapped])
    ok, tot = 0, 0
    per = {u: [0, 0] for u in units}
    conf = pd.DataFrame(0, index=units, columns=units, dtype=int)
    unassigned = 0
    for i in range(len(mapped)):
        keep = np.ones(len(mapped), bool)
        keep[i] = False
        cents, names = [], []
        for u in units:
            m = keep & (labels == u)
            if m.sum() == 0:
                continue
            cents.append(X[:, m].mean(axis=1))
            names.append(u)
        A = np.array(cents).T
        y = X[:, i]
        na = np.linalg.norm(A, axis=0)
        na[na == 0] = 1
        coef, _ = nnls(A / na, y / (np.linalg.norm(y) or 1))
        pred = names[int(np.argmax(coef))] if coef.max() > 0 else None
        tru = labels[i]
        tot += 1
        per[tru][1] += 1
        if pred == tru:
            ok += 1
            per[tru][0] += 1
        if pred is None:
            unassigned += 1
        else:
            conf.loc[tru, pred] += 1
    return ok / tot, per, conf, unassigned


def per_unit_frame(per, units):
    return pd.DataFrame(
        [(u, per[u][0], per[u][1],
          round(100 * per[u][0] / per[u][1], 1) if per[u][1] else None)
         for u in units],
        columns=['unit', 'correct', 'n', 'accuracy_pct'])


def main():
    if OUT_DIR.exists():
        sys.exit(f"Refusing to overwrite {OUT_DIR} — delete or rename it first.")

    tax = json.loads(TAX_SUMMARY.read_text())
    units16 = sorted(tax['analysis_phyla'])
    assert len(units16) == 16 and tax['taxonomy_release'] == 'ncbi-phylum-2026-08-04-v1'

    print("STAGE 1 — code-path gate: old labels (17 units) must match coauthor package")
    q_old, m_old, s_old = build(POS_TABLE, POS_META_OLD, 'original_header', OLD_UNITS, 'POS old')
    acc_old, per_old, _, _ = loo_nnls_confusion(q_old, m_old, s_old, OLD_UNITS)
    pu_old = per_unit_frame(per_old, OLD_UNITS)
    ref = pd.read_csv(GATE_OLD)
    if not pu_old.reset_index(drop=True).equals(ref):
        sys.exit("GATE FAIL: old-label per-unit results do not match loo_per_unit_old.csv")
    print(f"  GATE PASS ({100*acc_old:.1f}% old labels)")

    print("STAGE 2 — strict-16 run (locked release metadata)")
    q16, m16, s16 = build(POS_TABLE, POS_META_STRICT, 'original_header', units16, 'POS strict16')
    acc16, per16, conf16, un16 = loo_nnls_confusion(q16, m16, s16, units16)
    pu16 = per_unit_frame(per16, units16)
    print(f"  LOO accuracy (strict 16 phyla): {100*acc16:.1f}%  ({conf16.values.trace()}/{len(m16)})")
    assert (np.diag(conf16.values) == pu16['correct'].values).all()
    assert conf16.values.sum() + un16 == len(m16)

    OUT_DIR.mkdir(parents=True)
    conf16.to_csv(OUT_DIR / "loo_confusion_matrix_strict16.csv")
    pu16.to_csv(OUT_DIR / "loo_per_unit_strict16.csv", index=False)
    json.dump({
        'status': 'reimplementation — published 52.8% remains not reproducible',
        'taxonomy_release': 'ncbi-phylum-2026-08-04-v1 (locked strict policy)',
        'supersedes': 'suppfig4_loo_2026-08-11_v1 (19-unit 2026-08-03 coauthor partition)',
        'code_path_gate_old_labels': 'PASS (60.2%)',
        'loo_accuracy_strict16_pct': round(100 * acc16, 1),
        'chance_baseline_pct': round(100 / 16, 1),
        'n_samples': len(m16),
        'unassigned_predictions': un16,
        'units': units16,
        'filters': 'n_batches>=2, detection>=5%, mean_intensity>=500',
    }, open(OUT_DIR / "RUN_SUMMARY.json", 'w'), indent=2)
    print(f"\nWrote {OUT_DIR}")


if __name__ == '__main__':
    main()
