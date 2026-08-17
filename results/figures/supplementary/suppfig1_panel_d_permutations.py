#!/usr/bin/env python3
"""
Strict-16 Mantel permutation tests for Supp Fig 1 panel d.

Mirrors the published mantel_permutation_results.csv method column:
exact all-label permutations when the phylum family has n <= 7 (n! <= 5040,
matching every published exact row), seeded random 9,999 label permutations
otherwise. Overall-Mantel seeds reuse the published 20720247 (POS) /
21211386 (NEG); per-batch random seeds are derived deterministically.
"""
import itertools
import json
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr

ROOT = Path(r"C:\Users\Shadow\Desktop\P2R")
POS_TABLE = ROOT / "analysis/analysis-15/03_alignment/consensus_aligned_table.csv"
NEG_TABLE = ROOT / "analysis/analysis-16/negative_mode/03_alignment/consensus_aligned_table.csv"
POS_META_S = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/sample_metadata_POS_ncbi_phylum.csv"
NEG_META_S = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/sample_metadata_NEG_ncbi_phylum.csv"
TAX_SUMMARY = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/taxonomy_summary.json"
OUT = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/suppfig1_full_strict16_2026-08-11_v1/panel_d_data"
PHYLO16 = OUT / "phylo_dist_16phyla.csv"

SEEDS = {'POS': 20720247, 'NEG': 21211386}


def mantel_perm(bc_cond_matrix, ph_matrix, mode, tag):
    n = ph_matrix.shape[0]
    iu = np.triu_indices(n, k=1)
    obs, p_pearson = pearsonr(ph_matrix[iu], bc_cond_matrix[iu])
    if n <= 7:
        perms = list(itertools.permutations(range(n)))
        rs = np.array([pearsonr(ph_matrix[np.ix_(p, p)][iu], bc_cond_matrix[iu])[0]
                       for p in perms])
        method = f"exact_all_{len(perms)}_label_permutations"
        seed = None
    else:
        # crc32, not hash(): Python randomises string hashing per process, so
        # hash-derived seeds made these p values differ between runs.
        seed = SEEDS[mode] if tag == 'Overall Mantel' else (
            SEEDS[mode] + zlib.crc32(tag.encode()) % 10**6)
        rng = np.random.default_rng(seed)
        rs = np.array([pearsonr(ph_matrix[np.ix_(p, p)][iu], bc_cond_matrix[iu])[0]
                       for p in (rng.permutation(n) for _ in range(9999))])
        method = "random_9999_label_permutations"
    p_greater = (np.sum(rs >= obs - 1e-12) + (0 if n <= 7 else 1)) / (len(rs) + (0 if n <= 7 else 1))
    p_two = (np.sum(np.abs(rs) >= abs(obs) - 1e-12) + (0 if n <= 7 else 1)) / (len(rs) + (0 if n <= 7 else 1))
    return dict(observed_r=obs, pearson_p_original_method=p_pearson,
                mantel_p_greater=p_greater, mantel_p_two_sided=p_two,
                permutation_method=method, seed=seed,
                null_mean=rs.mean(), null_sd=rs.std(ddof=1),
                null_min=rs.min(), null_max=rs.max(),
                n_phyla=n, n_pairs=n * (n - 1) // 2)


def build_profiles(table_csv, meta_csv, sample_col, units):
    df = pd.read_csv(table_csv, index_col=0, low_memory=False)
    meta = pd.read_csv(meta_csv)
    sp, sb = {}, {}
    for _, row in meta.iterrows():
        col = row[sample_col]
        if col in df.columns and row.get('phylum') in units:
            sp[col] = row['phylum']
            sb[col] = row['batch']
    cols = [c for c in df.columns if c.startswith('sample:') and c in sp]
    cross = df[df['n_batches'] >= 2]
    inten = cross[cols].fillna(0)
    det = (inten > 0).sum(axis=1) / len(cols)
    mi = inten.mean(axis=1)
    q = inten[(det >= 0.05) & (mi >= 500)]
    return q, cols, sp, sb


def main():
    if not PHYLO16.exists():
        sys.exit("run suppfig1_full_strict16.py first")
    phylo = pd.read_csv(PHYLO16, index_col=0)
    units = sorted(json.loads(TAX_SUMMARY.read_text())['analysis_phyla'])
    neg_col = 'sample_col' if 'sample_col' in pd.read_csv(NEG_META_S, nrows=1).columns else 'original_header'

    rows = []
    for mode, table, meta, scol in [('POS', POS_TABLE, POS_META_S, 'original_header'),
                                    ('NEG', NEG_TABLE, NEG_META_S, neg_col)]:
        q, cols, sp, sb = build_profiles(table, meta, scol, units)
        jobs = [('Overall Mantel', cols)]
        for batch in sorted(set(sb.values())):
            bs = [s for s in cols if sb[s] == batch]
            if len(set(sp[s] for s in bs)) >= 3:
                jobs.append((f"Within-batch {batch.split('_')[1]}", bs))
        for tag, ss in jobs:
            bp = sorted(set(sp[s] for s in ss) & set(phylo.index))
            prof = np.array([q[[s for s in ss if sp[s] == p]].mean(axis=1).values for p in bp])
            bc = squareform(pdist(prof, metric='braycurtis'))
            ph = phylo.loc[bp, bp].values
            res = mantel_perm(bc, ph, mode, tag)
            res.update(figure_context='Supp Fig. 1d (strict16)', test=tag, mode=mode)
            rows.append(res)
            print(f"[{mode}] {tag}: r={res['observed_r']:.3f} "
                  f"p_mantel={res['mantel_p_greater']:.4g} ({res['permutation_method']})")

    cols_order = ['figure_context', 'test', 'mode', 'n_phyla', 'n_pairs', 'observed_r',
                  'pearson_p_original_method', 'mantel_p_greater', 'mantel_p_two_sided',
                  'permutation_method', 'seed', 'null_mean', 'null_sd', 'null_min', 'null_max']
    pd.DataFrame(rows)[cols_order].to_csv(OUT / "mantel_permutation_results.csv", index=False)
    print(f"wrote {OUT / 'mantel_permutation_results.csv'}")


if __name__ == '__main__':
    main()
