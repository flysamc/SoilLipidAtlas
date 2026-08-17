#!/usr/bin/env python3
"""
Supplementary Figure 3 — MS2LDA motif x unit enrichment under the LOCKED strict
release ncbi-phylum-2026-08-04-v1.

POS: genuine re-run (port of the coauthor package's rerun_suppfig3.py with its
empirically recovered parameters: motif threshold 0.10, unit tested if >= 6
spectra, motif tested if >= 3 spectra, a >= 2, BH across the family).
Stage 1 gate: published motif_phylum_enrichment.csv recovered exactly
(per-unit spectrum counts, per-test a and p) before Stage 2 runs.

NEG: RE-DERIVATION, not a re-run (the published spectrum-to-unit map is
unrecoverable; documented 2026-08-03). The published NEG table's Fisher inputs
are validated by recomputing every p-value from its own counts, then units are
merged per the strict policy, viral spectra removed from the universe, and the
2x2 tables + BH recomputed on the reduced family.

Strict-16 label handling (per-FEATURE table, no sample identity):
  Euryarchaeota / Halobacteriota -> Methanobacteriota
  Crenarchaeota                  -> Thermoproteota
  Trachaeophyta / Tracheophyta / Magnoliophyta / Bryophyta / Marchantiophyta /
  Charophyta                     -> Streptophyta
  Virus / Nucleocytoviricota     -> excluded (study scope)
  Amoebozoa                      -> UNCHANGED and flagged: the Discosea/Evosea
                                    split is per-sample and cannot be applied here
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(r"C:\Users\Shadow\Desktop\P2R")
DOC_TOPIC = ROOT / ("external/SOILMASS_PRODUCER_RECOVERY_2026-08-04/payload/core/workspace/"
                    "analysis/analysis-15/04_biomarker_discovery/ms2lda/results_full/doc_topic_matrix.csv")
ATLAS = ROOT / "analysis/analysis-15/04_biomarker_discovery/04_platinum_diamond/atlas_expanded_final.csv"
PUB_POS = ROOT / "analysis/analysis-15/04_biomarker_discovery/ms2lda/motif_phylum_enrichment.csv"
PUB_NEG = ROOT / ("manuscript_2_clean/06_figures/source_folders/analysis18_09_figures/"
                  "fig3c_ms2lda_heatmap/data/motif_phylum_enrichment_neg.csv")
POS_TABLE = ROOT / "analysis/analysis-15/03_alignment/consensus_aligned_table.csv"
POS_META_S = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/sample_metadata_POS_ncbi_phylum.csv"
TAX_SUMMARY = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/taxonomy_summary.json"
OUT = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/suppfig3_ms2lda_strict16_2026-08-11_v3"

MOTIF_THRESHOLD = 0.10
MIN_SPECTRA_PER_UNIT = 6
MIN_SPECTRA_PER_MOTIF = 3
EXCLUDE_UNITS = ('Virus', 'Nucleocytoviricota')
# Labels without an NCBI phylum-rank ancestor: retained in the spectrum
# universe (the samples are in-study) but NOT tested as phylum units
# (taxonomy_policy.json: descriptive_only_labels).
DESCRIPTIVE_ONLY = ('Bicosoecida',)

RELABEL_16 = {
    'Euryarchaeota': 'Methanobacteriota',
    'Halobacteriota': 'Methanobacteriota',
    'Crenarchaeota': 'Thermoproteota',
    'Trachaeophyta': 'Streptophyta',
    'Tracheophyta': 'Streptophyta',
    'Magnoliophyta': 'Streptophyta',
    'Bryophyta': 'Streptophyta',
    'Marchantiophyta': 'Streptophyta',
    'Charophyta': 'Streptophyta',
    'Virus': 'Nucleocytoviricota',
}
KINGDOM = {
    'Bacillota': 'Bacteria', 'Actinomycetota': 'Bacteria', 'Pseudomonadota': 'Bacteria',
    'Cyanobacteriota': 'Bacteria', 'Methanobacteriota': 'Archaea', 'Thermoproteota': 'Archaea',
    'Ascomycota': 'Fungi', 'Basidiomycota': 'Fungi', 'Mucoromycota': 'Fungi',
    'Mortierellomycota': 'Fungi', 'Streptophyta': 'Plantae', 'Chlorophyta': 'Plantae',
    'Arthropoda': 'Animalia', 'Nematoda': 'Animalia', 'Mollusca': 'Animalia',
    'Amoebozoa': 'Protozoa', 'Bicosoecida': 'Protozoa', 'Cercozoa': 'Protozoa',
    'Heterolobosea': 'Protozoa', 'Discosea': 'Protozoa', 'Evosea': 'Protozoa',
}


def bh(p):
    p = np.asarray(p, float)
    n = len(p)
    o = np.argsort(p)
    q = np.empty(n)
    prev = 1.0
    for k in range(n - 1, -1, -1):
        prev = min(prev, p[o[k]] * n / (k + 1))
        q[o[k]] = prev
    return q


def enrichment(dt, unit_of, tag, skip_units=(), allowed_units=None):
    d = dt.copy()
    d['unit'] = d.feature_id.map(unit_of)
    d = d[d.unit.notna()]
    motif_cols = [c for c in d.columns if c.startswith('motif_')]
    total = len(d)
    counts = d.unit.value_counts()
    units = sorted(u for u in counts[counts >= MIN_SPECTRA_PER_UNIT].index
                   if str(u).strip() not in EXCLUDE_UNITS
                   and str(u).strip() not in skip_units
                   and (allowed_units is None or str(u).strip() in allowed_units))
    print(f"  [{tag}] spectra {total:,} | units >= {MIN_SPECTRA_PER_UNIT}: {len(units)}")
    rows = []
    for motif in motif_cols:
        mm = (d[motif] >= MOTIF_THRESHOLD).values
        nm = int(mm.sum())
        if nm < MIN_SPECTRA_PER_MOTIF:
            continue
        for u in units:
            pm = (d.unit == u).values
            nu = int(pm.sum())
            a = int((mm & pm).sum())
            b = nu - a
            c = nm - a
            dd = total - nm - b
            if a < 2 or c < 0 or dd < 0:
                continue
            odds, pval = fisher_exact([[a, b], [c, dd]], alternative='greater')
            obs, exp = a / nm, nu / total
            rows.append(dict(motif=motif, unit=u, n_motif_in_unit=a, n_motif_total=nm,
                             n_unit_total=nu, observed_frac=round(obs, 4),
                             expected_frac=round(exp, 4),
                             enrichment_ratio=round(obs / exp if exp else 0, 4),
                             odds_ratio=round(float(odds), 4), pvalue=float(pval)))
    out = pd.DataFrame(rows)
    out['padj'] = bh(out.pvalue.values)
    out['significant'] = out.padj < 0.05
    print(f"  [{tag}] tests {len(out):,} | significant: {int(out.significant.sum())} | "
          f"motifs with >=1 enriched unit: {out[out.significant].motif.nunique()}")
    return out


def renderer_frame(df):
    out = df.rename(columns={'unit': 'phylum', 'n_motif_in_unit': 'n_motif_in_phylum',
                             'n_unit_total': 'n_phylum_total'}).copy()
    out.insert(2, 'kingdom', out['phylum'].map(KINGDOM).fillna('Protozoa'))
    return out[['motif', 'phylum', 'kingdom', 'n_motif_in_phylum', 'n_motif_total',
                'n_phylum_total', 'observed_frac', 'expected_frac', 'enrichment_ratio',
                'odds_ratio', 'pvalue', 'padj', 'significant']]


def main():
    if OUT.exists():
        sys.exit(f"Refusing to overwrite {OUT}")

    dt = pd.read_csv(DOC_TOPIC, low_memory=False)
    atlas = pd.read_csv(ATLAS, usecols=['feature_id', 'phylum'],
                        low_memory=False).drop_duplicates('feature_id')
    old_map = dict(zip(atlas.feature_id, atlas.phylum))
    new_map = {k: RELABEL_16.get(v, v) for k, v in old_map.items()}

    # ---- Amoebozoa feature-level split (POS) ----
    # The Discosea/Evosea split is defined per SAMPLE. For each Amoebozoa-
    # labelled feature, compare detection fractions in the corrected Discosea
    # vs Evosea samples of the consensus table and assign the unit with the
    # higher fraction; features detected in neither stay unassigned (dropped
    # from the unit map, retained nowhere — counted and reported).
    tbl = pd.read_csv(POS_TABLE, index_col=0, low_memory=False)
    meta_s = pd.read_csv(POS_META_S)
    cols = {u: [r['original_header'] for _, r in meta_s.iterrows()
                if r.get('phylum') == u and r['original_header'] in tbl.columns]
            for u in ('Discosea', 'Evosea')}
    amoeba_feats = [f for f, v in old_map.items() if v == 'Amoebozoa' and f in tbl.index]
    split = {}
    unassigned = 0
    for f in amoeba_feats:
        fr = {u: float((tbl.loc[f, cs].fillna(0) > 0).mean()) for u, cs in cols.items()}
        if fr['Discosea'] == 0 and fr['Evosea'] == 0:
            unassigned += 1
            continue
        split[f] = 'Discosea' if fr['Discosea'] >= fr['Evosea'] else 'Evosea'
    for f, u in split.items():
        new_map[f] = u
    for f in amoeba_feats:
        if f not in split and new_map.get(f) == 'Amoebozoa':
            del new_map[f]
    n_disc = sum(1 for v in split.values() if v == 'Discosea')
    print(f"Amoebozoa split: {len(amoeba_feats)} features -> Discosea {n_disc}, "
          f"Evosea {len(split) - n_disc}, unassigned (no detection) {unassigned}")

    # ---------------- POS stage 1: reproduction gate ----------------
    print("POS STAGE 1 — reproduction (published labels)")
    old = enrichment(dt, old_map, 'POS old')
    pub = pd.read_csv(PUB_POS)
    cnt = old.groupby('unit').n_unit_total.first()
    pubcnt = pub.groupby('phylum').n_phylum_total.first()
    shared_units = sorted(set(cnt.index) & set(pubcnt.index))
    cnt_ok = (all(int(cnt[u]) == int(pubcnt[u]) for u in shared_units)
              and not (set(cnt.index) - set(pubcnt.index))
              and (set(pubcnt.index) - set(cnt.index)) <= set(EXCLUDE_UNITS))
    j = pub.rename(columns={'phylum': 'unit'}).merge(
        old, on=['motif', 'unit'], suffixes=('_pub', '_mine'))
    da = (j.n_motif_in_phylum - j.n_motif_in_unit).abs().max()
    dp = (j.pvalue_pub - j.pvalue_mine).abs().max()
    print(f"  unit counts match: {cnt_ok} | shared tests {len(j):,}/{len(pub):,} | "
          f"max |delta a| = {int(da)} | max |delta p| = {dp:.2e}")
    if not (cnt_ok and da == 0 and dp < 1e-9):
        sys.exit("GATE FAIL: published POS enrichment not recovered exactly")
    print("  POS REPRODUCTION PASS")

    # ---------------- POS stage 2: strict-16 ----------------
    analysis16 = set(json.loads(TAX_SUMMARY.read_text())['analysis_phyla'])
    print("POS STAGE 2 — strict-16 labels (Amoebozoa split; tested family = 16 analysis phyla)")
    new = enrichment(dt, new_map, 'POS strict16', skip_units=DESCRIPTIVE_ONLY,
                     allowed_units=analysis16)

    # ---------------- NEG: validate + re-derive ----------------
    print("NEG — validating published table's Fisher inputs")
    neg = pd.read_csv(PUB_NEG)
    tot_est = (neg.n_phylum_total / neg.expected_frac.replace(0, np.nan)).dropna()
    TOTAL_NEG = int(round(tot_est.median()))
    print(f"  inferred NEG spectrum universe: {TOTAL_NEG:,} "
          f"(estimates range {tot_est.min():.1f}-{tot_est.max():.1f})")
    dev = []
    for _, r in neg.iterrows():
        a = int(r.n_motif_in_phylum)
        b = int(r.n_phylum_total) - a
        c = int(r.n_motif_total) - a
        d = TOTAL_NEG - int(r.n_motif_total) - b
        _, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        dev.append(abs(pval - r.pvalue))
    dev = float(np.max(dev))
    print(f"  max |recomputed p - published p| over {len(neg):,} rows: {dev:.2e}")
    if dev > 1e-6:
        sys.exit("NEG VALIDATION FAIL: published Fisher inputs not self-consistent")

    print("NEG — strict-16 re-derivation (merge labels, drop viral spectra)")
    neg['unit'] = neg['phylum'].map(lambda x: RELABEL_16.get(x, x))
    viral = neg[neg.unit.isin(EXCLUDE_UNITS)]
    virus_total = int(viral.groupby('phylum').n_phylum_total.first().sum())
    virus_a = viral.groupby('motif').n_motif_in_phylum.sum()
    # Tested family restricted to the 16 analysis phyla (author decision
    # 2026-08-11): drops Cercozoa (below core replication) and Amoebozoa (not a
    # valid analysis unit; unsplittable here — its spectra stay in the
    # universe, untested). Bicosoecida (descriptive-only) likewise untested.
    keep = neg[~neg.unit.isin(EXCLUDE_UNITS) & neg.unit.isin(analysis16)].copy()
    unit_tot = keep.groupby(['unit', 'phylum']).n_phylum_total.first() \
                   .groupby('unit').sum().astype(int)
    total_new = TOTAL_NEG - virus_total
    rows = []
    for (motif, unit), g in keep.groupby(['motif', 'unit']):
        a = int(g.n_motif_in_phylum.sum())
        nm = int(g.n_motif_total.iloc[0]) - int(virus_a.get(motif, 0))
        nu = int(unit_tot[unit])
        b = nu - a
        c = nm - a
        d = total_new - nm - b
        if a < 2 or c < 0 or d < 0:
            continue
        odds, pval = fisher_exact([[a, b], [c, d]], alternative='greater')
        obs, exp = (a / nm if nm else 0), nu / total_new
        rows.append(dict(motif=motif, unit=unit, n_motif_in_unit=a, n_motif_total=nm,
                         n_unit_total=nu, observed_frac=round(obs, 4),
                         expected_frac=round(exp, 4),
                         enrichment_ratio=round(obs / exp if exp else 0, 4),
                         odds_ratio=round(float(odds), 4), pvalue=float(pval)))
    neg16 = pd.DataFrame(rows)
    neg16['padj'] = bh(neg16.pvalue.values)
    neg16['significant'] = neg16.padj < 0.05
    print(f"  NEG strict16: tests {len(neg16):,} | significant {int(neg16.significant.sum())} | "
          f"motifs with >=1 enriched unit: {neg16[neg16.significant].motif.nunique()} | "
          f"units {neg16.unit.nunique()}")

    # ---------------- write ----------------
    OUT.mkdir(parents=True)
    data_dir = OUT / "r_render" / "data" / "supp_ms2lda"
    data_dir.mkdir(parents=True)
    old.to_csv(OUT / "pos_enrichment_reproduced_old.csv", index=False)
    renderer_frame(new).to_csv(data_dir / "motif_phylum_enrichment_pos.csv", index=False)
    renderer_frame(neg16).to_csv(data_dir / "motif_phylum_enrichment_neg.csv", index=False)
    json.dump({
        'taxonomy_release': 'ncbi-phylum-2026-08-04-v1',
        'pos': {'mode': 'genuine re-run', 'gate': 'PASS (counts, a, p exact)',
                'n_significant': int(new.significant.sum()),
                'n_motifs_enriched': int(new[new.significant].motif.nunique()),
                'n_units_tested': int(new.unit.nunique())},
        'neg': {'mode': 'RE-DERIVATION from published table (map unrecoverable)',
                'validation': f'max |recomputed p - published p| = {dev:.2e}',
                'inferred_universe': TOTAL_NEG, 'viral_spectra_removed': virus_total,
                'n_significant': int(neg16.significant.sum()),
                'n_motifs_enriched': int(neg16[neg16.significant].motif.nunique()),
                'n_units_tested': int(neg16.unit.nunique())},
        'params': {'motif_threshold': MOTIF_THRESHOLD, 'min_spectra_per_unit': MIN_SPECTRA_PER_UNIT,
                   'min_spectra_per_motif': MIN_SPECTRA_PER_MOTIF, 'min_overlap_a': 2},
        'amoebozoa_handling': {
            'pos': 'feature-level split via corrected-sample detection fractions '
                   f'(Discosea {n_disc}, Evosea {len(split) - n_disc}, unassigned {unassigned}) '
                   '— documented reassignment, needs author sign-off',
            'neg': 'reported unmodified (re-derivation from published table; no feature or sample identity)'},
        'bicosoecida': 'not tested as a phylum unit (descriptive-only label, taxonomy_policy.json); spectra retained in universe',
    }, open(OUT / "RUN_SUMMARY.json", 'w'), indent=2)
    print(f"DONE -> {OUT}")


if __name__ == '__main__':
    main()
