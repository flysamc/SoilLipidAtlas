#!/usr/bin/env python3
"""
Supplementary Figure 8 rerun — cross-method fingerprint validation under the
LOCKED strict release ncbi-phylum-2026-08-04-v1 (16 analysis phyla).

The published producer (method_validation_analysis20/06_figure/
fig_cross_method_validation.py and its substrate builder) lived in a deleted
worktree. Methods were recovered by fitting the frozen per-method top-K union
sizes in figures_r/data/supp_cross_method/panel_a_mantel_curves.csv:

  SIMPER        classic pairwise Bray-Curtis contribution (vegan form), raw
                intensities                       -> EXACT (all five unions)
  SCBD          per-phylum block of the global SCBD sum-of-squares on
                Hellinger-transformed data        -> EXACT (all five unions)
  CAP           PCoA of Bray-Curtis on TSS data, first 7 axes, one-vs-rest
                least-squares discriminant, features ranked by |Pearson r|
                with the discriminant scores on Hellinger data
                                                  -> BOUNDED best fit
                (aggregate union deviation ~2.2%; exact config unrecoverable)
  L1 stability  stability selection: LARS lasso entry order over B=50
                stratified half-subsamples (seed 20260811) on standardised
                Hellinger data, ranked by selection frequency then mean entry
                step, remainder by |point-biserial r|
                                                  -> DECLARED reimplementation
                (stochastic in the original; frozen unions not exactly
                reproducible without the historical seed)

Old-label (17-unit) union sizes are computed first as gates/deviation reports;
SIMPER and SCBD must match EXACTLY or the run aborts. Then everything is
recomputed on the strict-16 substrate and the five panel CSVs are written in
the schemas the authors' supp_fig8_cross_method.R consumes unchanged.

Panel notes:
  b baseline   full-substrate NNLS LOO (ported, gated code path from
               suppfig4_loo_strict16.py; 106/164 = 64.6%)
  c            the published "769 ClimGrass spectrally verified soil features"
               is replaced by the corrected strict-16 5-ppm verified mapping
               (strict16_verified_simper_mapping_5ppm.csv, 722 rows)
  d            kingdoms from taxonomy_policy ecological groups, display-mapped
               (Viridiplantae->Plantae, Protists->Protozoa) to match
               soilmass_style.R KINGDOM_COLOURS
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.optimize import nnls
from scipy.spatial.distance import pdist, squareform
from sklearn.linear_model import lars_path

ROOT = Path(r"C:\Users\Shadow\Desktop\P2R")
POS_TABLE = ROOT / "analysis/analysis-15/03_alignment/consensus_aligned_table.csv"
POS_META_OLD = ROOT / "analysis/analysis-15/04_biomarker_discovery/01_composite_scoring/sample_metadata.csv"
POS_META_STRICT = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/sample_metadata_POS_ncbi_phylum.csv"
TAX_SUMMARY = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/taxonomy_summary.json"
POLICY = ROOT / "paper2_repro/config/taxonomy_policy.json"
VERIFIED = ROOT / ("outputs/analysis/ncbi-phylum-2026-08-04-v1/climgrass/"
                   "strict16_all_simper_sensitivity_2026-08-06/strict16_verified_simper_mapping_5ppm.csv")
FROZEN_A = ROOT / "manuscript_2_clean/06_figures/figures_r/data/supp_cross_method/panel_a_mantel_curves.csv"
OUT_DIR = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/suppfig8_cross_method_strict16_2026-08-11_v1"

OLD_UNITS = sorted([
    'Actinomycetota', 'Amoebozoa', 'Arthropoda', 'Ascomycota', 'Bacillota', 'Basidiomycota',
    'Bryophyta', 'Chlorophyta', 'Euryarchaeota', 'Marchantiophyta', 'Methanobacteriota',
    'Mollusca', 'Mucoromycota', 'Nematoda', 'Pseudomonadota', 'Thermoproteota', 'Trachaeophyta'])
KS = [100, 250, 500, 1000, 2500]
SEED = 20260811
L1_B = 50
CAP_M = 7
KINGDOM_DISPLAY = {'Viridiplantae': 'Plantae', 'Protists': 'Protozoa'}


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def build(meta_csv, units):
    df = pd.read_csv(POS_TABLE, index_col=0, low_memory=False)
    meta = pd.read_csv(meta_csv)
    s2u = {r['original_header']: r['phylum'] for _, r in meta.iterrows()
           if r['original_header'] in df.columns and r.get('phylum') in units}
    mapped = [c for c in df.columns if str(c).startswith('sample:') and c in s2u]
    labels = np.array([s2u[s] for s in mapped])
    cross = df[df['n_batches'] >= 2]
    inten = cross[mapped].fillna(0)
    det = (inten > 0).sum(axis=1) / len(mapped)
    mi = inten.mean(axis=1)
    q = inten[(det >= 0.05) & (mi >= 500)]
    X = q.values.T.astype(float)
    fids = np.asarray(q.index)
    log(f"substrate: {X.shape[0]} samples x {X.shape[1]:,} features, {len(units)} units")
    return X, labels, fids


# ---------------------------------------------------------------- methods
def rank_simper(X, labels, units):
    """Classic pairwise Bray-Curtis SIMPER, raw intensities. EXACT-recovered."""
    scores = {}
    for u in units:
        t = labels == u
        G, R = X[t], X[~t]
        acc = np.zeros(X.shape[1])
        for i in range(G.shape[0]):
            num = np.abs(G[i][None, :] - R)
            den = (G[i][None, :] + R).sum(axis=1, keepdims=True)
            den[den == 0] = 1
            acc += (num / den).sum(axis=0)
        scores[u] = acc / (G.shape[0] * R.shape[0])
    return scores


def rank_scbd(X, labels, units):
    """Per-phylum block of global SCBD SS on Hellinger data. EXACT-recovered."""
    H = np.sqrt(X / X.sum(axis=1, keepdims=True))
    gmean = H.mean(axis=0)
    return {u: ((H[labels == u] - gmean) ** 2).sum(axis=0) for u in units}


def _corr_cols(M, v):
    Mc = M - M.mean(axis=0)
    vc = v - v.mean()
    num = Mc.T @ vc
    den = np.sqrt((Mc ** 2).sum(axis=0) * (vc ** 2).sum())
    den[den == 0] = 1
    return np.abs(num / den)


def rank_cap(X, labels, units):
    """Bounded best-fit CAP: PCoA(BC on TSS), m=7, lstsq discriminant,
    |corr| with Hellinger features."""
    n = X.shape[0]
    T = X / X.sum(axis=1, keepdims=True)
    H = np.sqrt(T)
    D = squareform(pdist(T, metric='braycurtis'))
    J = np.eye(n) - np.ones((n, n)) / n
    G = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh(G)
    o = np.argsort(w)[::-1]
    w, V = w[o], V[:, o]
    pos = w > 1e-10
    C = (V[:, pos] * np.sqrt(w[pos]))[:, :CAP_M]
    A = np.column_stack([np.ones(n), C])
    scores = {}
    for u in units:
        t = (labels == u).astype(float)
        beta, *_ = np.linalg.lstsq(A, t, rcond=None)
        scores[u] = _corr_cols(H, A @ beta)
    return scores


def rank_l1(X, labels, units, seed=SEED, B=L1_B):
    """Declared stability-selection reimplementation (see module docstring)."""
    rng = np.random.default_rng(seed)
    H = np.sqrt(X / X.sum(axis=1, keepdims=True))
    mu, sd = H.mean(axis=0), H.std(axis=0)
    sd[sd == 0] = 1
    Z = (H - mu) / sd
    n, p = Z.shape
    scores = {}
    for u in units:
        t = labels == u
        freq = np.zeros(p)
        entry_sum = np.zeros(p)
        for b in range(B):
            pos_idx = np.flatnonzero(t)
            neg_idx = np.flatnonzero(~t)
            take_pos = max(2, len(pos_idx) // 2) if len(pos_idx) >= 2 else len(pos_idx)
            sub = np.concatenate([
                rng.choice(pos_idx, size=take_pos, replace=False),
                rng.choice(neg_idx, size=len(neg_idx) // 2, replace=False)])
            y = t[sub].astype(float)
            y = y - y.mean()
            try:
                _, _, coefs = lars_path(Z[sub], y, method='lasso',
                                        max_iter=min(160, len(sub) - 1))
            except Exception:
                continue
            entered = (np.abs(coefs) > 0).argmax(axis=1).astype(float)
            active = (np.abs(coefs) > 0).any(axis=1)
            freq[active] += 1
            entry_sum[active] += entered[active]
        mean_entry = np.where(freq > 0, entry_sum / np.maximum(freq, 1), np.inf)
        tie = _corr_cols(H, t.astype(float))
        # lexicographic rank: freq desc, mean entry asc, |corr| desc
        order = np.lexsort((-tie, mean_entry, -freq))
        sc = np.empty(p)
        sc[order] = -np.arange(p, dtype=float)
        scores[u] = sc
    return scores


def unions(scores, units, ks=KS):
    out = {}
    for k in ks:
        s = set()
        for u in units:
            s.update(np.argsort(-scores[u])[:k].tolist())
        out[k] = sorted(s)
    return out


# ---------------------------------------------------------------- panels
def unit_centroids(X, labels, units):
    return np.array([X[labels == u].mean(axis=0) for u in units])


def bc_condensed(C):
    return pdist(C, metric='braycurtis')


def mantel_and_coph(cent, idx, full_cond):
    sub = bc_condensed(cent[:, idx])
    r = float(np.corrcoef(sub, full_cond)[0, 1])
    Z = linkage(sub, method='average')
    c, _ = cophenet(Z, sub)
    return r, float(c)


def loo_nnls(X, labels, units, col_idx):
    Xs = X[:, col_idx]
    n = X.shape[0]
    ok = 0
    for i in range(n):
        keep = np.ones(n, bool)
        keep[i] = False
        cents, names = [], []
        for u in units:
            m = keep & (labels == u)
            if m.sum() == 0:
                continue
            cents.append(Xs[m].mean(axis=0))
            names.append(u)
        A = np.array(cents).T
        y = Xs[i]
        na = np.linalg.norm(A, axis=0)
        na[na == 0] = 1
        coef, _ = nnls(A / na, y / (np.linalg.norm(y) or 1))
        pred = names[int(np.argmax(coef))] if coef.max() > 0 else None
        ok += (pred == labels[i])
    return ok, n


def main():
    if OUT_DIR.exists():
        sys.exit(f"Refusing to overwrite {OUT_DIR}")

    frozen = pd.read_csv(FROZEN_A)
    frozen_unions = {m: dict(zip(g.K, g.n_features)) for m, g in frozen.groupby('method')}

    # ---------------- Stage 1: old-label validation --------------
    log("STAGE 1 — old-label validation (17 units)")
    Xo, lo, _ = build(POS_META_OLD, OLD_UNITS)
    validation = {}
    method_funcs = [('simper', rank_simper), ('scbd', rank_scbd),
                    ('cap', rank_cap), ('stability_l1', rank_l1)]
    for name, fn in method_funcs:
        log(f"  old-label {name}...")
        sc = fn(Xo, lo, OLD_UNITS)
        mine = {k: len(v) for k, v in unions(sc, OLD_UNITS).items()}
        ref = frozen_unions[name]
        dev = {int(k): int(mine[k] - ref[k]) for k in KS}
        rel = max(abs(mine[k] - ref[k]) / ref[k] for k in KS)
        validation[name] = {'mine': mine, 'frozen': {int(k): int(ref[k]) for k in KS},
                            'delta': dev, 'max_rel_dev': round(rel, 4)}
        log(f"    unions {mine} vs frozen {ref} (max rel dev {rel:.2%})")
        if name in ('simper', 'scbd') and any(dev[k] != 0 for k in KS):
            sys.exit(f"GATE FAIL: {name} old-label unions must match frozen exactly")
    validation['simper']['mode'] = 'exact_gate_pass'
    validation['scbd']['mode'] = 'exact_gate_pass'
    validation['cap']['mode'] = 'bounded_best_fit'
    validation['stability_l1']['mode'] = 'declared_reimplementation_stochastic'

    # ---------------- Stage 2: strict-16 production --------------
    log("STAGE 2 — strict-16 production run")
    tax = json.loads(TAX_SUMMARY.read_text())
    units16 = sorted(tax['analysis_phyla'])
    X, lab, fids = build(POS_META_STRICT, units16)

    rankings = {}
    for name, fn in method_funcs:
        log(f"  strict16 {name}...")
        rankings[name] = fn(X, lab, units16)

    cent = unit_centroids(X, lab, units16)
    full_cond = bc_condensed(cent)
    Zfull = linkage(full_cond, method='average')
    full_coph, _ = cophenet(Zfull, full_cond)
    full_coph = float(full_coph)

    # panel a: mantel curves
    log("panel a: mantel curves")
    rows_a = []
    union_sets = {}
    for name, _ in method_funcs:
        u = unions(rankings[name], units16)
        union_sets[name] = u
        for k in KS:
            r, c = mantel_and_coph(cent, u[k], full_cond)
            rows_a.append((name, k, len(u[k]), r, c, full_coph))
    df_a = pd.DataFrame(rows_a, columns=['method', 'K', 'n_features', 'mantel_r_vs_full',
                                         'cophenetic_r', 'full_cophenetic_r'])

    # panel a null: K random features per phylum, union, 1000 iterations
    log("panel a: random null (1000 iterations)")
    rng = np.random.default_rng(SEED)
    p = X.shape[1]
    rows_null = []
    for k in KS:
        rs, usz = [], []
        for _ in range(1000):
            s = set()
            for _u in units16:
                s.update(rng.choice(p, size=k, replace=False).tolist())
            idx = sorted(s)
            usz.append(len(idx))
            sub = bc_condensed(cent[:, idx])
            rs.append(float(np.corrcoef(sub, full_cond)[0, 1]))
        rs = np.array(rs)
        rows_null.append((k, int(np.mean(usz)), rs.mean(), rs.std(ddof=1),
                          np.quantile(rs, 0.05), np.quantile(rs, 0.95)))
    df_null = pd.DataFrame(rows_null, columns=['K_per_phylum', 'approx_K_union',
                                               'null_mean', 'null_sd', 'null_q05', 'null_q95'])

    # panel b: LOO accuracy
    log("panel b: LOO accuracy (baseline)")
    ok_full, n_full = loo_nnls(X, lab, units16, np.arange(p))
    baseline = ok_full / n_full
    log(f"  full-substrate baseline: {ok_full}/{n_full} = {baseline:.1%}")
    rows_b = []
    for name, _ in method_funcs:
        for k in KS:
            idx = union_sets[name][k]
            log(f"  LOO {name} K={k} ({len(idx):,} features)")
            ok, ntot = loo_nnls(X, lab, units16, np.array(idx))
            rows_b.append((name, k, len(idx), ok / ntot, ok, ntot, baseline))
    df_b = pd.DataFrame(rows_b, columns=['method', 'K', 'n_features', 'accuracy',
                                         'n_correct', 'n_total', 'full_substrate_baseline'])

    # panel c: verified soil-feature overlap (corrected 5-ppm mapping)
    log("panel c: soil overlap")
    ver = pd.read_csv(VERIFIED)
    ver_ids = set(ver['feature_id'].unique())
    fid_pos = {f: i for i, f in enumerate(fids)}
    ver_idx = {fid_pos[f] for f in ver_ids if f in fid_pos}
    n_cg = len(ver_idx)
    log(f"  verified soil features on substrate: {n_cg} (of {len(ver_ids)} unique ids)")
    rows_c = []
    for name, _ in method_funcs:
        for k in KS:
            idx = set(union_sets[name][k])
            ov = len(idx & ver_idx)
            rows_c.append((name, k, len(idx), n_cg, ov, ov / n_cg, ov / len(idx)))
    df_c = pd.DataFrame(rows_c, columns=['method', 'K', 'n_features', 'n_climgrass_features',
                                         'n_overlap', 'frac_of_climgrass_in_topK',
                                         'frac_of_topK_in_climgrass'])

    # panel d: consensus per phylum at K=500
    log("panel d: consensus per phylum")
    policy = json.loads(POLICY.read_text())
    eco = policy['ecological_group']
    rows_d = []
    for u in units16:
        sets = [set(np.argsort(-rankings[name][u])[:500].tolist())
                for name, _ in method_funcs]
        allf = set().union(*sets)
        counts = {f: sum(f in s for s in sets) for f in allf}
        vals = np.array(list(counts.values()))
        kingdom = KINGDOM_DISPLAY.get(eco[u], eco[u])
        rows_d.append((u, kingdom, len(allf), int((vals == 4).sum()),
                       int((vals >= 3).sum()), int((vals >= 2).sum()),
                       int((vals == 1).sum())))
    df_d = pd.DataFrame(rows_d, columns=['phylum', 'kingdom', 'total', 'all4',
                                         'geq3', 'geq2', 'only_1'])

    # ---------------- write ----------------
    OUT_DIR.mkdir(parents=True)
    data_dir = OUT_DIR / "r_render" / "data" / "supp_cross_method"
    data_dir.mkdir(parents=True)
    df_a.to_csv(data_dir / "panel_a_mantel_curves.csv", index=False)
    df_null.to_csv(data_dir / "panel_a_random_null.csv", index=False)
    df_b.to_csv(data_dir / "panel_b_loo_accuracy.csv", index=False)
    df_c.to_csv(data_dir / "panel_c_climgrass_overlap.csv", index=False)
    df_d.to_csv(data_dir / "panel_d_consensus_per_phylum.csv", index=False)
    for df, nm in [(df_a, 'panel_a_mantel_curves'), (df_null, 'panel_a_random_null'),
                   (df_b, 'panel_b_loo_accuracy'), (df_c, 'panel_c_climgrass_overlap'),
                   (df_d, 'panel_d_consensus_per_phylum')]:
        df.to_csv(OUT_DIR / f"{nm}.csv", index=False)
    json.dump({
        'taxonomy_release': 'ncbi-phylum-2026-08-04-v1 (locked strict policy, 16 phyla)',
        'seed': SEED,
        'method_validation_old_labels': validation,
        'full_substrate': {'n_samples': int(X.shape[0]), 'n_features': int(p),
                           'loo_baseline': f"{ok_full}/{n_full} = {baseline:.3f}",
                           'full_cophenetic_r': full_coph},
        'panel_c_verified_set': {'source': str(VERIFIED.relative_to(ROOT)),
                                 'n_on_substrate': n_cg,
                                 'replaces': 'published 769-feature list (old labels)'},
        'l1_params': {'B': L1_B, 'subsample': 'stratified half', 'path': 'LARS lasso',
                      'data': 'standardised Hellinger', 'max_steps': 160},
        'cap_params': {'distance': 'Bray-Curtis on TSS', 'axes': CAP_M,
                       'discriminant': 'one-vs-rest least squares',
                       'feature_score': '|Pearson r| vs Hellinger features'},
    }, open(OUT_DIR / "RUN_SUMMARY.json", 'w'), indent=2)
    log(f"DONE -> {OUT_DIR}")


if __name__ == '__main__':
    main()
