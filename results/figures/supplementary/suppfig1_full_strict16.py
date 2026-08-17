#!/usr/bin/env python3
"""
Supplementary Figure 1 (submitted 4-panel version) — strict-16 update.

Submitted panels (SLA supplementary docx):
  a  cross-batch PCoA, 168 samples, organism-group colours (silhouette 0.090)
  b  same ordination coloured by batch (silhouette 0.122)
  c  within-batch PCoA, batch 02 (silhouette 0.106)
  d  within-batch Mantel r vs phylogenetic distance per batch (POS)

Producers: manuscript_2_clean/06_figures/source_folders/analysis18_09_figures/
fig2{a,b,c,d}_* (the repo figures_r/supp_fig1_within_batch.R is an older
two-panel draft and is NOT the submitted figure).

What this script does:
  Panels a+b: sample-level audit of pcoa_coords.csv against the strict release;
    kingdom (=colour) changes abort, phylum-column relabels are applied to a
    copy. Coordinates and silhouettes untouched: the ordination itself is not
    recomputed, so the panels stay visually identical.
  Panel c: already relabelled in suppfig1_relabel_2026-08-11_v1 (copied here).
  Panel d: REPRODUCE-FIRST — a faithful port of fig2d_mantel_forest/data/
    freeze_data.py must match the frozen forest_data.csv exactly; only then is
    the strict-16 forest computed (16 phyla, corrected metadata, five-rank
    cladogram extended per the coauthor-package CLADO_NEW rules with
    Streptophyta at ('Eukaryota','Archaeplastida','Plantae','Streptophyta') —
    a judgement call for coauthor sign-off, replacing the split plant units).
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import pearsonr
from sklearn.metrics import silhouette_score

ROOT = Path(r"C:\Users\Shadow\Desktop\P2R")
SRCF = ROOT / "manuscript_2_clean/06_figures/source_folders/analysis18_09_figures"
POS_TABLE = ROOT / "analysis/analysis-15/03_alignment/consensus_aligned_table.csv"
POS_META_OLD = ROOT / "analysis/analysis-15/04_biomarker_discovery/01_composite_scoring/sample_metadata.csv"
NEG_TABLE = ROOT / "analysis/analysis-16/negative_mode/03_alignment/consensus_aligned_table.csv"
NEG_META_OLD = ROOT / "analysis/analysis-16/negative_mode/04_biomarker_discovery/01_composite_scoring/sample_metadata.csv"
# original path (analysis-15/.../phylo_dist_matrix.csv) is absent from this
# package; the identical published 17-phyla matrix survives in the fig1b freeze
PHYLO_CSV = ROOT / "manuscript_2_clean/06_figures/source_folders/analysis18_09_figures/fig1b_mantel/data/phylo_dist_17phyla.csv"
POS_META_S = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/sample_metadata_POS_ncbi_phylum.csv"
NEG_META_S = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/sample_metadata_NEG_ncbi_phylum.csv"
TAX_SUMMARY = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/taxonomy_summary.json"
POLICY = ROOT / "paper2_repro/config/taxonomy_policy.json"
PREV_C = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/suppfig1_relabel_2026-08-11_v1/r_render/data/supp_within_batch"
OUT = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/suppfig1_full_strict16_2026-08-11_v1"

TARGET_PHYLA = sorted([
    'Actinomycetota', 'Amoebozoa', 'Arthropoda', 'Ascomycota', 'Bacillota', 'Basidiomycota',
    'Bryophyta', 'Chlorophyta', 'Euryarchaeota', 'Marchantiophyta', 'Methanobacteriota',
    'Mollusca', 'Mucoromycota', 'Nematoda', 'Pseudomonadota', 'Thermoproteota', 'Trachaeophyta'])

CLADO_16 = {
    'Actinomycetota':    ('Bacteria', 'Terrabacteria', 'Bacteria', 'Actinobacteria'),
    'Bacillota':         ('Bacteria', 'Terrabacteria', 'Bacteria', 'Firmicutes'),
    'Pseudomonadota':    ('Bacteria', 'Gracilicutes', 'Bacteria', 'Proteobacteria'),
    'Methanobacteriota': ('Archaea', 'Methanobacteriota', 'Archaea', 'Methanobacteria'),
    'Thermoproteota':    ('Archaea', 'TACK', 'Archaea', 'Thermoproteota'),
    'Ascomycota':        ('Eukaryota', 'Opisthokonta', 'Fungi', 'Dikarya'),
    'Basidiomycota':     ('Eukaryota', 'Opisthokonta', 'Fungi', 'Dikarya'),
    'Mucoromycota':      ('Eukaryota', 'Opisthokonta', 'Fungi', 'Basal Fungi'),
    'Arthropoda':        ('Eukaryota', 'Opisthokonta', 'Animalia', 'Ecdysozoa'),
    'Nematoda':          ('Eukaryota', 'Opisthokonta', 'Animalia', 'Ecdysozoa'),
    'Mollusca':          ('Eukaryota', 'Opisthokonta', 'Animalia', 'Lophotrochozoa'),
    'Chlorophyta':       ('Eukaryota', 'Archaeplastida', 'Plantae', 'Green algae'),
    'Streptophyta':      ('Eukaryota', 'Archaeplastida', 'Plantae', 'Streptophyta'),
    'Discosea':          ('Eukaryota', 'Amoebozoa', 'Protozoa', 'Discosea'),
    'Evosea':            ('Eukaryota', 'Amoebozoa', 'Protozoa', 'Evosea'),
    'Heterolobosea':     ('Eukaryota', 'Excavata', 'Protozoa', 'Discoba'),
}


def phylo_matrix(units, clado):
    n = len(units)
    D = np.zeros((n, n))
    for i, a in enumerate(units):
        for j, b in enumerate(units):
            if i == j:
                continue
            shared = 0
            for x, y in zip(clado[a], clado[b]):
                if x == y:
                    shared += 1
                else:
                    break
            D[i, j] = 5 - shared
    return pd.DataFrame(D, index=units, columns=units)


def forest(table_csv, meta_csv, sample_col, polarity, units, phylo):
    """Faithful port of fig2d_mantel_forest/data/freeze_data.py::compute."""
    df = pd.read_csv(table_csv, index_col=0, low_memory=False)
    meta = pd.read_csv(meta_csv)
    sp, sk, sb = {}, {}, {}
    for _, row in meta.iterrows():
        col = row[sample_col]
        if col in df.columns and row.get('phylum') in units:
            sp[col] = row['phylum']
            sk[col] = row.get('kingdom', '')
            sb[col] = row['batch']
    sample_cols = [c for c in df.columns if c.startswith('sample:') and c in sp]
    cross = df[df['n_batches'] >= 2]
    intensity = cross[sample_cols].fillna(0)
    det = (intensity > 0).sum(axis=1) / len(sample_cols)
    mi = intensity.mean(axis=1)
    intensity_q = intensity[(det >= 0.05) & (mi >= 500)]
    phylo_set = set(phylo.index)

    rows = []
    phyla_p = sorted(set(sp.values()) & phylo_set)
    prof = {p: intensity_q[[s for s in sample_cols if sp[s] == p]].mean(axis=1).values
            for p in phyla_p}
    mat = np.array([prof[p] for p in phyla_p])
    bc = squareform(pdist(mat, metric='braycurtis'))
    ph = phylo.loc[phyla_p, phyla_p].values
    iu = np.triu_indices(len(phyla_p), k=1)
    r, p = pearsonr(ph[iu], bc[iu])
    rows.append(dict(test='Overall Mantel', mode=polarity, value=r, pvalue=p,
                     n=len(phyla_p), metric='Mantel r', group='Mantel'))

    for batch in sorted(set(sb[s] for s in sample_cols)):
        bs = [s for s in sample_cols if sb[s] == batch]
        bp = sorted(set(sp[s] for s in bs) & phylo_set)
        if len(bp) < 3:
            continue
        m = np.array([intensity_q[[s for s in bs if sp[s] == p]].mean(axis=1).values
                      for p in bp])
        b_bc = squareform(pdist(m, metric='braycurtis'))
        b_ph = phylo.loc[bp, bp].values
        iu_b = np.triu_indices(len(bp), k=1)
        r_b, p_b = pearsonr(b_ph[iu_b], b_bc[iu_b])
        bnum = batch.split('_')[1]
        rows.append(dict(test=f'Within-batch {bnum}', mode=polarity, value=r_b,
                         pvalue=p_b, n=len(bp), metric='Mantel r', group='Mantel'))

    iq_T = intensity_q.T
    s_bc = squareform(pdist(iq_T.values, metric='braycurtis'))
    sil_k = silhouette_score(s_bc, [sk[s] for s in iq_T.index], metric='precomputed')
    sil_b = silhouette_score(s_bc, [sb[s] for s in iq_T.index], metric='precomputed')
    rows.append(dict(test='Kingdom silhouette', mode=polarity, value=sil_k,
                     pvalue=np.nan, n=len(sample_cols), metric='Silhouette', group='Silhouette'))
    rows.append(dict(test='Batch silhouette', mode=polarity, value=sil_b,
                     pvalue=np.nan, n=len(sample_cols), metric='Silhouette', group='Silhouette'))
    return pd.DataFrame(rows)


def main():
    if OUT.exists():
        sys.exit(f"Refusing to overwrite {OUT}")

    policy = json.loads(POLICY.read_text())
    eco = policy['ecological_group']
    disp = {'Viridiplantae': 'Plantae', 'Protists': 'Protozoa'}
    meta_s = pd.read_csv(POS_META_S)
    ph_new = dict(zip(meta_s['original_header'], meta_s['phylum']))

    # ---------- panels a+b: audit + relabel ----------
    print("panels a+b: audit + relabel")
    src = SRCF / "fig2a_umap_kingdom/data/pcoa_coords.csv"
    df = pd.read_csv(src)
    relabels = []
    for _, r in df.iterrows():
        s = r['sample']
        np_ = ph_new.get(s)
        if np_ is None:
            print(f"  NOTE {s}: not in corrected metadata (kept as-is: {r['phylum']})")
            continue
        nk = disp.get(eco.get(np_, ''), eco.get(np_, r['kingdom']))
        if np_ != r['phylum']:
            if nk != r['kingdom'] and np_ not in ('Bicosoecida',):
                sys.exit(f"KINGDOM CHANGE for {s}: {r['kingdom']} -> {nk} — panel colours would change; stopping for review")
            relabels.append({'sample': s, 'old': r['phylum'], 'new': np_})
    out_ab = OUT / "panel_ab_data"
    out_ab.mkdir(parents=True)
    df2 = df.copy()
    df2['phylum'] = [ph_new.get(s, p) for s, p in zip(df['sample'], df['phylum'])]
    assert df2[['sample', 'pcoa1', 'pcoa2', 'kingdom', 'batch']].equals(
        df[['sample', 'pcoa1', 'pcoa2', 'kingdom', 'batch']])
    df2.to_csv(out_ab / "pcoa_coords.csv", index=False)
    print(f"  {len(relabels)} of {len(df)} samples relabelled; kingdoms unchanged")

    # ---------- panel c: copy previous relabel ----------
    out_c = OUT / "panel_c_data"
    shutil.copytree(PREV_C, out_c)

    # ---------- panel d: reproduce-first ----------
    print("panel d: reproduction check (17 phyla, original metadata)")
    phylo_old = pd.read_csv(PHYLO_CSV, index_col=0)
    rep_pos = forest(POS_TABLE, POS_META_OLD, 'original_header', 'POS', TARGET_PHYLA, phylo_old)
    rep_neg = forest(NEG_TABLE, NEG_META_OLD, 'sample_col', 'NEG', TARGET_PHYLA, phylo_old)
    mine = pd.concat([rep_pos, rep_neg], ignore_index=True)
    frozen = pd.read_csv(SRCF / "fig2d_mantel_forest/data/forest_data.csv")
    if len(mine) != len(frozen):
        sys.exit(f"GATE FAIL: row count {len(mine)} vs {len(frozen)}")
    key = ['test', 'mode']
    j = frozen.merge(mine, on=key, suffixes=('_pub', '_new'))
    dv = (j['value_pub'] - j['value_new']).abs().max()
    dn = (j['n_pub'] - j['n_new']).abs().max()
    print(f"  max |delta value| = {dv:.2e}, max |delta n| = {dn}")
    if dv > 1e-9 or dn != 0:
        print(j[(j['value_pub'] - j['value_new']).abs() > 1e-9].to_string())
        sys.exit("GATE FAIL: panel d reproduction does not match frozen forest_data.csv")
    print("  REPRODUCTION PASS")

    # ---------- panel d: strict-16 ----------
    print("panel d: strict-16 rerun")
    tax = json.loads(TAX_SUMMARY.read_text())
    units16 = sorted(tax['analysis_phyla'])
    assert set(units16) == set(CLADO_16), "cladogram does not cover the strict-16 set"
    phylo16 = phylo_matrix(units16, CLADO_16)
    neg_meta_cols = pd.read_csv(NEG_META_S, nrows=1).columns
    neg_col = 'sample_col' if 'sample_col' in neg_meta_cols else 'original_header'
    new_pos = forest(POS_TABLE, POS_META_S, 'original_header', 'POS', units16, phylo16)
    new_neg = forest(NEG_TABLE, NEG_META_S, neg_col, 'NEG', units16, phylo16)
    strict = pd.concat([new_pos, new_neg], ignore_index=True)

    out_d = OUT / "panel_d_data"
    out_d.mkdir()
    mine.to_csv(out_d / "forest_data_reproduced_17phyla.csv", index=False)
    strict.to_csv(out_d / "forest_data.csv", index=False)
    phylo16.to_csv(out_d / "phylo_dist_16phyla.csv")
    pos_rows = strict[(strict['mode'] == 'POS') & (strict['group'] == 'Mantel')]
    print(strict.to_string(index=False))

    json.dump({
        'taxonomy_release': 'ncbi-phylum-2026-08-04-v1 (16 analysis phyla)',
        'panel_ab': {'relabels': relabels, 'n_relabelled': len(relabels),
                     'visual_change': 'none (kingdom colours, coords, silhouettes untouched)'},
        'panel_c': 'copied from suppfig1_relabel_2026-08-11_v1 (25 relabels, visual unchanged)',
        'panel_d': {
            'reproduction_gate': 'PASS (exact vs frozen forest_data.csv, POS+NEG)',
            'cladogram': 'CLADO_NEW rules from coauthor package rerun_impact.py; Streptophyta placed at (Eukaryota, Archaeplastida, Plantae, Streptophyta) — coauthor sign-off needed',
            'strict16_pos_within_batch_range': [
                round(float(pos_rows[pos_rows.test.str.startswith("Within")]['value'].min()), 3),
                round(float(pos_rows[pos_rows.test.str.startswith("Within")]['value'].max()), 3)],
        },
    }, open(OUT / "RUN_SUMMARY.json", 'w'), indent=2, default=str)
    print(f"DONE -> {OUT}")


if __name__ == '__main__':
    main()
