#!/usr/bin/env python3
"""
Supplementary Figure 2 — clade-conserved vs clade-exclusive features under the
LOCKED strict release ncbi-phylum-2026-08-04-v1 (16 analysis phyla), POS + NEG.

Method: faithful port of the coauthor package's rerun_suppfig2.py (itself copied
from analysis-18 fig1d_clade_bars/data/freeze_data.py), extended to negative
mode. Reproduce-first: Stage 1 must recover the published clade_counts.csv
EXACTLY (all nine clades, POS and NEG, substrates 44,534 / 14,393) or the run
aborts.

Clade framework: the five-rank Adl et al. 2019 + NCBI cladogram (Table S15) —
NOT the SSU tree, which measures distance, not membership (author decision
2026-08-11). Strict-16 clade list:
  - Euryarchaeota_sg: dropped (was the same taxon twice, taxid 28890)
  - Embryophyta: dropped (its three plant units merge into Streptophyta,
    a single unit is not a clade)
  - Amoebozoa (Discosea + Evosea): enters (real clade created by the split)
  - All Archaea = Methanobacteriota + Thermoproteota (halophiles inside
    Methanobacteriota under the locked policy)
  - All Plantae = Chlorophyta + Streptophyta
  - Opisthokonta (all Fungi + all Animalia units): PROPOSED addition, marked
    in the output; answers the reviewer monophyly objection (C48). Keeps the
    clade count at nine if adopted.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(r"C:\Users\Shadow\Desktop\P2R")
POS_TABLE = ROOT / "analysis/analysis-15/03_alignment/consensus_aligned_table.csv"
NEG_TABLE = ROOT / "analysis/analysis-16/negative_mode/03_alignment/consensus_aligned_table.csv"
POS_META_OLD = ROOT / "analysis/analysis-15/04_biomarker_discovery/01_composite_scoring/sample_metadata.csv"
NEG_META_OLD = ROOT / "analysis/analysis-16/negative_mode/04_biomarker_discovery/01_composite_scoring/sample_metadata.csv"
POS_META_S = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/sample_metadata_POS_ncbi_phylum.csv"
NEG_META_S = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/sample_metadata_NEG_ncbi_phylum.csv"
TAX_SUMMARY = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/taxonomy/taxonomy_summary.json"
PUB = ROOT / "manuscript_2_clean/06_figures/figures_r/data/supp_clade/clade_counts.csv"
OUT = ROOT / "outputs/analysis/ncbi-phylum-2026-08-04-v1/suppfig2_clade_strict16_2026-08-11_v1"
DET = 0.05

UNITS_17 = sorted(['Actinomycetota', 'Amoebozoa', 'Arthropoda', 'Ascomycota', 'Bacillota',
                   'Basidiomycota', 'Bryophyta', 'Chlorophyta', 'Euryarchaeota', 'Marchantiophyta',
                   'Methanobacteriota', 'Mollusca', 'Mucoromycota', 'Nematoda', 'Pseudomonadota',
                   'Thermoproteota', 'Trachaeophyta'])

CLADES_OLD = [
    ('Dikarya', 'Fungi', ['Ascomycota', 'Basidiomycota']),
    ('Ecdysozoa', 'Animalia', ['Arthropoda', 'Nematoda']),
    ('Euryarchaeota_sg', 'Archaea', ['Euryarchaeota', 'Methanobacteriota']),
    ('All Animalia', 'Animalia', ['Arthropoda', 'Nematoda', 'Mollusca']),
    ('All Archaea', 'Archaea', ['Euryarchaeota', 'Methanobacteriota', 'Thermoproteota']),
    ('All Fungi', 'Fungi', ['Ascomycota', 'Basidiomycota', 'Mucoromycota']),
    ('All Bacteria', 'Bacteria', ['Actinomycetota', 'Bacillota', 'Pseudomonadota']),
    ('Embryophyta', 'Plantae', ['Bryophyta', 'Marchantiophyta', 'Trachaeophyta']),
    ('All Plantae', 'Plantae', ['Bryophyta', 'Marchantiophyta', 'Trachaeophyta', 'Chlorophyta']),
]

CLADES_16 = [
    ('Dikarya', 'Fungi', ['Ascomycota', 'Basidiomycota']),
    ('Ecdysozoa', 'Animalia', ['Arthropoda', 'Nematoda']),
    ('Amoebozoa', 'Protozoa', ['Discosea', 'Evosea']),
    ('All Animalia', 'Animalia', ['Arthropoda', 'Nematoda', 'Mollusca']),
    ('All Archaea', 'Archaea', ['Methanobacteriota', 'Thermoproteota']),
    ('All Fungi', 'Fungi', ['Ascomycota', 'Basidiomycota', 'Mucoromycota']),
    ('All Bacteria', 'Bacteria', ['Actinomycetota', 'Bacillota', 'Pseudomonadota']),
    ('All Plantae', 'Plantae', ['Chlorophyta', 'Streptophyta']),
    ('Opisthokonta', 'Fungi', ['Ascomycota', 'Basidiomycota', 'Mucoromycota',
                               'Arthropoda', 'Nematoda', 'Mollusca']),  # PROPOSED
]
PROPOSED = {'Opisthokonta'}


def detection_matrix(table_csv, meta_csv, sample_col, units, tag):
    df = pd.read_csv(table_csv, index_col=0, low_memory=False)
    meta = pd.read_csv(meta_csv)
    s2u = {r[sample_col]: r['phylum'] for _, r in meta.iterrows()
           if r[sample_col] in df.columns and r.get('phylum') in units}
    mapped = [c for c in df.columns if str(c).startswith('sample:') and c in s2u]
    M = df[mapped].fillna(0).values.astype(float)
    det = (M > 0).sum(axis=1) / len(mapped)
    mi = M.mean(axis=1)
    qf = (df['n_batches'].values >= 2) & (det >= 0.05) & (mi >= 500)
    Mq = M[qf]
    labels = np.array([s2u[c] for c in mapped])
    present = {}
    for u in units:
        m = labels == u
        if m.sum() == 0:
            continue
        present[u] = ((Mq[:, m] > 0).sum(axis=1) / m.sum()) >= DET
    print(f"  [{tag}] quality features: {int(qf.sum()):,} | samples: {len(mapped)} | units: {len(present)}")
    return present, int(qf.sum())


def clade_counts(present, clades, all_units, prefix):
    rows = []
    nfeat = len(next(iter(present.values())))
    for name, kingdom, members in clades:
        members = [m for m in members if m in present]
        if len(members) < 2:
            print(f"    SKIP {name}: fewer than 2 member units with data")
            continue
        shared = np.ones(nfeat, dtype=bool)
        for m in members:
            shared &= present[m]
        outside_any = np.zeros_like(shared)
        for u in all_units:
            if u in present and u not in members:
                outside_any |= present[u]
        excl = shared & ~outside_any
        rows.append(dict(clade=name, kingdom=kingdom, n_members=len(members),
                         members='; '.join(members),
                         **{f'{prefix}_shared': int(shared.sum()),
                            f'{prefix}_exclusive': int(excl.sum())}))
    return pd.DataFrame(rows)


def main():
    if OUT.exists():
        sys.exit(f"Refusing to overwrite {OUT}")

    neg_col_old = 'sample_col' if 'sample_col' in pd.read_csv(NEG_META_OLD, nrows=1).columns else 'original_header'
    neg_col_s = 'sample_col' if 'sample_col' in pd.read_csv(NEG_META_S, nrows=1).columns else 'original_header'

    print("STAGE 1 — reproduction (17 units, 9 published clades, POS + NEG)")
    p_pos, nq_pos = detection_matrix(POS_TABLE, POS_META_OLD, 'original_header', UNITS_17, 'POS old')
    p_neg, nq_neg = detection_matrix(NEG_TABLE, NEG_META_OLD, neg_col_old, UNITS_17, 'NEG old')
    old = clade_counts(p_pos, CLADES_OLD, UNITS_17, 'pos').merge(
        clade_counts(p_neg, CLADES_OLD, UNITS_17, 'neg')[['clade', 'neg_shared', 'neg_exclusive']],
        on='clade')
    pub = pd.read_csv(PUB)
    j = pub.merge(old, on='clade', suffixes=('_pub', '_mine'))
    dev = max((j[f'{m}_{k}_pub'] - j[f'{m}_{k}_mine']).abs().max()
              for m in ('pos', 'neg') for k in ('shared', 'exclusive'))
    print(f"  max |delta| across all clades/modes: {int(dev)}")
    if dev != 0 or nq_pos != 44534 or nq_neg != 14393:
        print(j.to_string())
        sys.exit("GATE FAIL: published clade counts or substrates not recovered exactly")
    print("  REPRODUCTION PASS (all 9 clades, POS + NEG, substrates 44,534 / 14,393)")

    print("STAGE 2 — strict-16 run")
    units16 = sorted(json.loads(TAX_SUMMARY.read_text())['analysis_phyla'])
    s_pos, sq_pos = detection_matrix(POS_TABLE, POS_META_S, 'original_header', units16, 'POS strict16')
    s_neg, sq_neg = detection_matrix(NEG_TABLE, NEG_META_S, neg_col_s, units16, 'NEG strict16')
    new = clade_counts(s_pos, CLADES_16, units16, 'pos').merge(
        clade_counts(s_neg, CLADES_16, units16, 'neg')[['clade', 'neg_shared', 'neg_exclusive']],
        on='clade')
    new['proposed'] = new['clade'].isin(PROPOSED)
    print(new.to_string(index=False))

    OUT.mkdir(parents=True)
    data_dir = OUT / "r_render" / "data" / "supp_clade"
    data_dir.mkdir(parents=True)
    old.to_csv(OUT / "clade_counts_reproduced_17units.csv", index=False)
    new.to_csv(OUT / "clade_counts_strict16.csv", index=False)
    # renderer schema: no 'proposed' column
    new.drop(columns=['proposed']).to_csv(data_dir / "clade_counts.csv", index=False)
    json.dump({
        'pos_n_quality_features': sq_pos, 'neg_n_quality_features': sq_neg,
        'detection_threshold_per_phylum': DET,
        'n_clades': int(len(new)),
        'method': 'Shared = detected in EVERY phylum of clade. Exclusive = shared AND not detected in any phylum outside clade. "Detected in phylum" = >=5% detection rate within that phylum.',
        'framework': 'five-rank Adl 2019 + NCBI cladogram (Table S15), strict-16 membership',
        'reproduction_gate': 'PASS — published clade_counts.csv recovered exactly (POS+NEG)',
        'dropped': ['Euryarchaeota_sg (same taxon twice)', 'Embryophyta (collapses into Streptophyta)'],
        'added': ['Amoebozoa (Discosea + Evosea)', 'Opisthokonta (PROPOSED, pending author decision)'],
        'caveat': 'clade-conserved counts scale with samples per unit (5% threshold < 1 sample for every unit); open author decision',
    }, open(data_dir / "freeze_stats.json", 'w'), indent=2)
    print(f"DONE -> {OUT}")


if __name__ == '__main__':
    main()
