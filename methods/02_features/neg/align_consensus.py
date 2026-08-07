#!/usr/bin/env python3
"""
Step 03: Cross-batch consensus alignment for NEG mode.

Algorithm:
  1. Load all 6 RT-corrected batch feature tables
  2. Keep only biological sample columns (from metadata)
  3. Use reference batch (OE23-NEG) as seed; iteratively align other batches
  4. For each unaligned feature, find best match in consensus by m/z (5 ppm) + RT (tolerance)
  5. If match found → merge intensities into existing consensus row
  6. If no match → add new consensus row
  7. Output: consensus_aligned_table.csv with one row per consensus feature

Usage:
  python align_consensus.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import time
import warnings
warnings.filterwarnings('ignore')

# ── Config ──────────────────────────────────────────────────────
BASE = Path(__file__).parent.parent
LOESS_DIR = BASE / "02_loess_correction"
META_PATH = BASE / "00_sample_mapping" / "neg_sample_metadata.csv"
OUT_DIR = Path(__file__).parent
FIG_DIR = OUT_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

PPM_TOL = 5.0          # Orbitrap mass accuracy
RT_TOL  = 0.5          # minutes — post-LOESS corrected RT tolerance
RT_MIN  = 1.5          # min RT window (chromatographic elution range)
RT_MAX  = 21.0         # max RT window
REF_BATCH = "OE23-NEG" # Reference batch (same as POS pipeline)

# Batch order: reference first, then largest to smallest for best seed coverage
BATCH_ORDER = [
    "OE23-NEG",       # reference, 46K features
    "OE21-4-NEG",     # 64K — largest
    "ALL-25-2-NEG",   # 36K
    "OE26-1-NEG",     # 31K
    "OE11-3-NEG",     # 12K
    "OE25-1-NEG",     #  5K — smallest
]

BATCH_LABELS = {
    "OE23-NEG":      "batch_03_OE23-NEG",
    "OE21-4-NEG":    "batch_02_OE21-4-NEG",
    "ALL-25-2-NEG":  "batch_06_ALL-25-2-NEG",
    "OE26-1-NEG":    "batch_04_OE26-1-NEG",
    "OE11-3-NEG":    "batch_01_OE11-3-NEG",
    "OE25-1-NEG":    "batch_05_OE25-1-NEG",
}


def load_batch(batch_name, meta):
    """Load RT-corrected batch, keep only biological sample columns."""
    path = LOESS_DIR / f"{batch_name}_rt_corrected.csv"
    df = pd.read_csv(path)

    # Identify biological sample columns for this batch
    batch_label = BATCH_LABELS[batch_name]
    bio_meta = meta[(meta['batch'] == batch_label) & (meta['include_in_analysis'] == 'Yes')]
    bio_headers = set(bio_meta['original_header'].values)

    # Find matching columns
    sample_cols = [c for c in df.columns if c in bio_headers]

    # Build clean dataframe
    out = pd.DataFrame({
        'mz': df['row m/z'].values,
        'rt': df['rt_corrected'].values,
        'row_id': df['row ID'].values,
    })
    for col in sample_cols:
        # Map to clean sample name
        name = bio_meta.loc[bio_meta['original_header'] == col, 'sample_name'].values[0]
        out[f"sample:{name}"] = df[col].values

    # Remove features with zero intensity across all bio samples
    intensity_cols = [c for c in out.columns if c.startswith('sample:')]
    mask = out[intensity_cols].sum(axis=1) > 0
    out = out[mask].reset_index(drop=True)

    # Filter to chromatographic RT window
    rt_mask = (out['rt'] >= RT_MIN) & (out['rt'] <= RT_MAX)
    n_before = len(out)
    out = out[rt_mask].reset_index(drop=True)
    if n_before - len(out) > 0:
        print(f"    RT filtered: {n_before} → {len(out)} (removed {n_before - len(out)} outside [{RT_MIN}, {RT_MAX}] min)")

    print(f"  {batch_name}: {len(df)} total → {len(out)} with bio signal, "
          f"{len(sample_cols)} bio samples")
    return out


def find_best_match(mz, rt, consensus_mz, consensus_rt, ppm_tol, rt_tol):
    """Find best matching consensus feature by combined m/z + RT distance.

    Returns index into consensus arrays, or -1 if no match.
    Uses vectorized operations for speed.
    """
    # m/z filter: ppm tolerance
    ppm_err = np.abs(consensus_mz - mz) / mz * 1e6
    mz_mask = ppm_err <= ppm_tol

    if not np.any(mz_mask):
        return -1

    # RT filter
    rt_err = np.abs(consensus_rt[mz_mask] - rt)
    rt_mask = rt_err <= rt_tol

    if not np.any(rt_mask):
        return -1

    # Combined distance: normalize ppm to [0,1] range and RT to [0,1] range
    # Weight equally: score = ppm/ppm_tol + rt_err/rt_tol
    candidates = np.where(mz_mask)[0]
    rt_errs = np.abs(consensus_rt[candidates] - rt)
    valid = rt_errs <= rt_tol
    candidates = candidates[valid]
    rt_errs = rt_errs[valid]
    ppm_errs = ppm_err[candidates]

    score = ppm_errs / ppm_tol + rt_errs / rt_tol
    best = np.argmin(score)
    return candidates[best]


def align_batches(batch_data, batch_order):
    """Align all batches into a consensus feature table.

    Strategy: Use m/z-sorted arrays for efficient windowed search.
    """
    ref_name = batch_order[0]
    ref = batch_data[ref_name]
    sample_cols_ref = [c for c in ref.columns if c.startswith('sample:')]

    # Initialize consensus from reference batch
    n_ref = len(ref)
    consensus = {
        'consensus_mz': ref['mz'].values.copy(),
        'consensus_rt': ref['rt'].values.copy(),
        'n_merged': np.ones(n_ref, dtype=np.int32),  # total merge count (can exceed 6)
        'batch_set': [{ref_name} for _ in range(n_ref)],  # set of unique batches
        'ref_batch_id': ref['row_id'].values.copy(),
    }
    # Sample intensities — start with reference
    sample_intensities = {}
    for col in sample_cols_ref:
        sample_intensities[col] = ref[col].values.copy()

    print(f"\n  Consensus seeded with {n_ref} features from {ref_name}")

    # Sort consensus by m/z for efficient windowed search
    sort_idx = np.argsort(consensus['consensus_mz'])
    consensus['consensus_mz'] = consensus['consensus_mz'][sort_idx]
    consensus['consensus_rt'] = consensus['consensus_rt'][sort_idx]
    consensus['n_merged'] = consensus['n_merged'][sort_idx]
    consensus['batch_set'] = [consensus['batch_set'][i] for i in sort_idx]
    consensus['ref_batch_id'] = consensus['ref_batch_id'][sort_idx]
    for col in sample_intensities:
        sample_intensities[col] = sample_intensities[col][sort_idx]

    stats = {'batch': [ref_name], 'features_in': [n_ref],
             'matched': [n_ref], 'new': [0], 'consensus_size': [n_ref]}

    # Align remaining batches
    for batch_name in batch_order[1:]:
        t0 = time.time()
        batch = batch_data[batch_name]
        batch_sample_cols = [c for c in batch.columns if c.startswith('sample:')]

        # Add empty columns for this batch's samples
        n_consensus = len(consensus['consensus_mz'])
        for col in batch_sample_cols:
            if col not in sample_intensities:
                sample_intensities[col] = np.zeros(n_consensus, dtype=np.float64)

        matched = 0
        new_features = []

        cmz = consensus['consensus_mz']
        crt = consensus['consensus_rt']

        for i in range(len(batch)):
            mz_i = batch['mz'].iloc[i]
            rt_i = batch['rt'].iloc[i]

            # Binary search for m/z window
            mz_lo = mz_i * (1 - PPM_TOL / 1e6)
            mz_hi = mz_i * (1 + PPM_TOL / 1e6)
            idx_lo = np.searchsorted(cmz, mz_lo, side='left')
            idx_hi = np.searchsorted(cmz, mz_hi, side='right')

            best_idx = -1
            if idx_lo < idx_hi:
                # Candidates within m/z window
                rt_errs = np.abs(crt[idx_lo:idx_hi] - rt_i)
                rt_valid = rt_errs <= RT_TOL
                if np.any(rt_valid):
                    # Pick closest combined distance
                    cand_slice = np.arange(idx_lo, idx_hi)[rt_valid]
                    ppm_errs = np.abs(cmz[cand_slice] - mz_i) / mz_i * 1e6
                    scores = ppm_errs / PPM_TOL + rt_errs[rt_valid] / RT_TOL
                    best_idx = cand_slice[np.argmin(scores)]

            if best_idx >= 0:
                # Merge: update consensus m/z as weighted average, RT as average
                n = consensus['n_merged'][best_idx]
                consensus['consensus_mz'][best_idx] = (
                    consensus['consensus_mz'][best_idx] * n + mz_i) / (n + 1)
                consensus['consensus_rt'][best_idx] = (
                    consensus['consensus_rt'][best_idx] * n + rt_i) / (n + 1)
                consensus['n_merged'][best_idx] += 1
                consensus['batch_set'][best_idx].add(batch_name)
                # Fill intensities
                for col in batch_sample_cols:
                    sample_intensities[col][best_idx] = batch[col].iloc[i]
                matched += 1
            else:
                # New feature — collect and add after loop
                new_row = {
                    'mz': mz_i, 'rt': rt_i,
                    'row_id': batch['row_id'].iloc[i],
                    'batch_name': batch_name,
                }
                for col in batch_sample_cols:
                    new_row[col] = batch[col].iloc[i]
                new_features.append(new_row)

        # Add new features to consensus
        n_new = len(new_features)
        if n_new > 0:
            new_mz = np.array([f['mz'] for f in new_features])
            new_rt = np.array([f['rt'] for f in new_features])
            new_nmerged = np.ones(n_new, dtype=np.int32)
            new_batch_sets = [{f['batch_name']} for f in new_features]
            new_refid = np.array([f['row_id'] for f in new_features])

            consensus['consensus_mz'] = np.concatenate([cmz, new_mz])
            consensus['consensus_rt'] = np.concatenate([crt, new_rt])
            consensus['n_merged'] = np.concatenate([consensus['n_merged'], new_nmerged])
            consensus['batch_set'] = consensus['batch_set'] + new_batch_sets
            consensus['ref_batch_id'] = np.concatenate([consensus['ref_batch_id'], new_refid])

            # Extend sample intensity arrays
            for col in sample_intensities:
                old = sample_intensities[col]
                ext = np.zeros(n_new, dtype=np.float64)
                sample_intensities[col] = np.concatenate([old, ext])

            # Fill new features' intensities
            offset = n_consensus
            for j, feat in enumerate(new_features):
                for col in batch_sample_cols:
                    sample_intensities[col][offset + j] = feat.get(col, 0.0)

            # Re-sort by m/z
            sort_idx = np.argsort(consensus['consensus_mz'])
            consensus['consensus_mz'] = consensus['consensus_mz'][sort_idx]
            consensus['consensus_rt'] = consensus['consensus_rt'][sort_idx]
            consensus['n_merged'] = consensus['n_merged'][sort_idx]
            consensus['batch_set'] = [consensus['batch_set'][i] for i in sort_idx]
            consensus['ref_batch_id'] = consensus['ref_batch_id'][sort_idx]
            for col in sample_intensities:
                sample_intensities[col] = sample_intensities[col][sort_idx]

        elapsed = time.time() - t0
        n_total = len(consensus['consensus_mz'])
        print(f"  {batch_name}: {len(batch)} features → "
              f"{matched} matched, {n_new} new → consensus: {n_total} "
              f"({elapsed:.1f}s)")

        stats['batch'].append(batch_name)
        stats['features_in'].append(len(batch))
        stats['matched'].append(matched)
        stats['new'].append(n_new)
        stats['consensus_size'].append(n_total)

    return consensus, sample_intensities, pd.DataFrame(stats)


def build_output_table(consensus, sample_intensities):
    """Build final consensus DataFrame."""
    n = len(consensus['consensus_mz'])

    # Count samples detected per feature
    all_sample_cols = sorted(sample_intensities.keys())
    n_detected = np.zeros(n, dtype=np.int32)
    for col in all_sample_cols:
        n_detected += (sample_intensities[col] > 0).astype(np.int32)

    # Derive n_batches and batch strings from batch_set
    n_batches = np.array([len(s) for s in consensus['batch_set']], dtype=np.int32)
    batch_strs = [','.join(sorted(s)) for s in consensus['batch_set']]

    # Build DataFrame
    out = pd.DataFrame({
        'feature_id': [f"NEG_{i}" for i in range(n)],
        'consensus_mz': consensus['consensus_mz'],
        'consensus_rt': consensus['consensus_rt'],
        'n_batches': n_batches,
        'batches': batch_strs,
        'ref_batch_id': consensus['ref_batch_id'],
        'n_samples_detected': n_detected,
    })

    # Add sample columns
    for col in all_sample_cols:
        out[col] = sample_intensities[col]

    # Sort by m/z
    out = out.sort_values('consensus_mz').reset_index(drop=True)
    out['feature_id'] = [f"NEG_{i}" for i in range(len(out))]

    return out


def make_figures(consensus_df, alignment_stats):
    """Generate QC figures."""

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # 1. Features per n_batches
    ax = axes[0, 0]
    counts = consensus_df['n_batches'].value_counts().sort_index()
    colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00']
    bars = ax.bar(counts.index, counts.values, color=colors[:len(counts)])
    ax.set_xlabel('Number of batches')
    ax.set_ylabel('Number of features')
    ax.set_title('Feature batch coverage')
    for bar, v in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width()/2, v + 500,
                f'{v:,}', ha='center', va='bottom', fontsize=9)
    total = len(consensus_df)
    multi = (consensus_df['n_batches'] >= 2).sum()
    ax.text(0.95, 0.95, f'Total: {total:,}\nMulti-batch: {multi:,} ({100*multi/total:.1f}%)',
            transform=ax.transAxes, ha='right', va='top', fontsize=9,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # 2. m/z distribution
    ax = axes[0, 1]
    ax.hist(consensus_df['consensus_mz'], bins=200, color='#0072B2', alpha=0.7, edgecolor='none')
    ax.set_xlabel('m/z')
    ax.set_ylabel('Feature count')
    ax.set_title('m/z distribution of consensus features')
    ax.axvline(consensus_df['consensus_mz'].median(), color='red', ls='--', label=f"median={consensus_df['consensus_mz'].median():.1f}")
    ax.legend()

    # 3. RT distribution
    ax = axes[1, 0]
    ax.hist(consensus_df['consensus_rt'], bins=200, color='#009E73', alpha=0.7, edgecolor='none')
    ax.set_xlabel('RT (min)')
    ax.set_ylabel('Feature count')
    ax.set_title('RT distribution of consensus features')

    # 4. Alignment cascade
    ax = axes[1, 1]
    x = range(len(alignment_stats))
    ax.bar(x, alignment_stats['matched'], label='Matched to consensus', color='#0072B2')
    ax.bar(x, alignment_stats['new'], bottom=alignment_stats['matched'],
           label='New (unmatched)', color='#E69F00')
    ax.set_xticks(x)
    ax.set_xticklabels([b.replace('-NEG', '') for b in alignment_stats['batch']],
                       rotation=45, ha='right')
    ax.set_ylabel('Features')
    ax.set_title('Alignment cascade: matched vs new per batch')
    ax.legend()

    plt.tight_layout()
    fig.savefig(FIG_DIR / 'alignment_summary.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # 5. Samples detected histogram
    fig, ax = plt.subplots(figsize=(10, 5))
    sample_cols = [c for c in consensus_df.columns if c.startswith('sample:')]
    n_detected = (consensus_df[sample_cols] > 0).sum(axis=1)
    ax.hist(n_detected, bins=min(100, n_detected.max()), color='#CC79A7', alpha=0.7, edgecolor='none')
    ax.set_xlabel('Number of samples detected')
    ax.set_ylabel('Feature count')
    ax.set_title(f'Detection breadth across {len(sample_cols)} biological samples')
    ax.axvline(n_detected.median(), color='red', ls='--',
               label=f'median={n_detected.median():.0f}')
    ax.legend()
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'detection_breadth.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    print(f"  Figures saved to {FIG_DIR}/")


def main():
    print("=" * 65)
    print("Step 03: Cross-Batch Consensus Alignment — NEG Mode")
    print("=" * 65)
    print(f"  m/z tolerance: {PPM_TOL} ppm")
    print(f"  RT tolerance:  {RT_TOL} min (post-LOESS)")
    print(f"  Reference:     {REF_BATCH}")
    print()

    # Load metadata
    meta = pd.read_csv(META_PATH)
    print("Loading RT-corrected batch files...")

    # Load all batches
    batch_data = {}
    total_features = 0
    for batch in BATCH_ORDER:
        df = load_batch(batch, meta)
        batch_data[batch] = df
        total_features += len(df)

    print(f"\n  Total input features: {total_features:,}")
    print(f"\n{'─' * 65}")
    print("Aligning batches into consensus...")

    consensus, sample_intensities, stats = align_batches(batch_data, BATCH_ORDER)

    print(f"\n{'─' * 65}")
    print("Building output table...")

    consensus_df = build_output_table(consensus, sample_intensities)

    # Save
    out_path = OUT_DIR / "consensus_aligned_table.csv"
    consensus_df.to_csv(out_path, index=False)
    stats.to_csv(OUT_DIR / "alignment_stats.csv", index=False)

    # Summary
    sample_cols = [c for c in consensus_df.columns if c.startswith('sample:')]
    n_total = len(consensus_df)
    n_multi = (consensus_df['n_batches'] >= 2).sum()
    n_detected = consensus_df['n_samples_detected']

    print(f"\n{'=' * 65}")
    print("CONSENSUS ALIGNMENT SUMMARY")
    print(f"{'=' * 65}")
    print(f"  Total consensus features:  {n_total:,}")
    print(f"  Single-batch features:     {n_total - n_multi:,} ({100*(n_total-n_multi)/n_total:.1f}%)")
    print(f"  Multi-batch features:      {n_multi:,} ({100*n_multi/n_total:.1f}%)")
    print(f"  Biological sample columns: {len(sample_cols)}")
    print(f"  Median samples detected:   {n_detected.median():.0f}")
    print(f"  Features in ≥2 samples:    {(n_detected >= 2).sum():,}")
    print()

    by_nbatch = consensus_df['n_batches'].value_counts().sort_index()
    for k, v in by_nbatch.items():
        print(f"    In {k} batch(es): {v:,}")
    print()

    # m/z and RT ranges
    print(f"  m/z range: {consensus_df['consensus_mz'].min():.4f} – {consensus_df['consensus_mz'].max():.4f}")
    print(f"  RT range:  {consensus_df['consensus_rt'].min():.2f} – {consensus_df['consensus_rt'].max():.2f} min")
    print(f"\n  Output: {out_path}")

    # Figures
    print("\nGenerating QC figures...")
    make_figures(consensus_df, stats)

    print(f"\n{'=' * 65}")
    print("Step 03 complete.")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
