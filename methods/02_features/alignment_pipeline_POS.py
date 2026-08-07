"""
Analysis-15: Cross-Batch Alignment Pipeline (6 FBMN Batches)
=============================================================
Direct cross-batch alignment without classical networking backbone.
Uses FBMN quantification tables as input.

Pipeline:
  Step 1: Load all batches, find cross-batch anchors (features in ≥3 batches)
  Step 2: Calculate RT drift per anchor using reference batch
  Step 3: Fit LOESS RT correction curves
  Step 4: Build unified consensus table with corrected RTs
  Step 5: Integrate library matches and molecular family info

Reference batch: OE23-POS (largest, best chromatography)
Total batches: 6
"""

import pandas as pd
import numpy as np
import json
import os
import re
import sys
import gc
from collections import defaultdict

try:
    from statsmodels.nonparametric.smoothers_lowess import lowess
except ImportError:
    os.system(f"{sys.executable} -m pip install statsmodels --break-system-packages -q")
    from statsmodels.nonparametric.smoothers_lowess import lowess

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError:
    os.system(f"{sys.executable} -m pip install matplotlib --break-system-packages -q")
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

# ==========================================
# CONFIGURATION
# ==========================================
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.dirname(SCRIPT_DIR)
FBMN_DIR = os.path.join(os.path.dirname(ANALYSIS_DIR), "FBMN_all_batches_POS")

# Output directories
ANCHOR_DIR = os.path.join(ANALYSIS_DIR, "01_anchor_discovery")
LOESS_DIR = os.path.join(ANALYSIS_DIR, "02_loess_correction")
ALIGN_DIR = os.path.join(ANALYSIS_DIR, "03_alignment")
FIGURE_DIR = os.path.join(ANALYSIS_DIR, "05_figures")

for d in [ANCHOR_DIR, LOESS_DIR, ALIGN_DIR, FIGURE_DIR]:
    os.makedirs(d, exist_ok=True)

# 6 FBMN batch quantification tables
BATCH_QUANT = {
    "OE11-3-POS":  os.path.join(FBMN_DIR, "batch_01_OE11-3-POS/results/OE11-3-POS_quantification_table.csv"),
    "OE21-4-POS":  os.path.join(FBMN_DIR, "batch_02_OE21-4-POS/results/OE21-4-POS_quantification_table.csv"),
    "OE23-POS":    os.path.join(FBMN_DIR, "batch_03_OE23-POS/results/OE23-POS_quantification_table.csv"),
    "OE26-1POS":   os.path.join(FBMN_DIR, "batch_04_OE26-1POS/results/OE26-1POS_quantification_table.csv"),
    "OE25-1-ALL":  os.path.join(FBMN_DIR, "batch_05_OE25-1-ALL/results/OE25-1-ALL_quantification_table.csv"),
    "ALL-25-2-POS": os.path.join(FBMN_DIR, "batch_06_ALL-25-2/results/ALL-25-2-POS_quantification_table.csv"),
}

# Library match files
BATCH_LIBRARY = {
    "OE11-3-POS":  os.path.join(FBMN_DIR, "batch_01_OE11-3-POS/results/OE11-3-POS_library_matches.tsv"),
    "OE21-4-POS":  os.path.join(FBMN_DIR, "batch_02_OE21-4-POS/results/OE21-4-POS_library_matches.tsv"),
    "OE23-POS":    os.path.join(FBMN_DIR, "batch_03_OE23-POS/results/OE23-POS_library_matches.tsv"),
    "OE26-1POS":   os.path.join(FBMN_DIR, "batch_04_OE26-1POS/results/OE26-1POS_library_matches.tsv"),
    "OE25-1-ALL":  os.path.join(FBMN_DIR, "batch_05_OE25-1-ALL/results/OE25-1-ALL_library_matches.tsv"),
    "ALL-25-2-POS": os.path.join(FBMN_DIR, "batch_06_ALL-25-2/results/ALL-25-2-POS_library_matches.tsv"),
}

# Cluster summary files (for molecular family info)
BATCH_CLUSTERS = {
    "OE11-3-POS":  os.path.join(FBMN_DIR, "batch_01_OE11-3-POS/results/OE11-3-POS_cluster_summary.tsv"),
    "OE21-4-POS":  os.path.join(FBMN_DIR, "batch_02_OE21-4-POS/results/OE21-4-POS_cluster_summary.tsv"),
    "OE23-POS":    os.path.join(FBMN_DIR, "batch_03_OE23-POS/results/OE23-POS_cluster_summary.tsv"),
    "OE26-1POS":   os.path.join(FBMN_DIR, "batch_04_OE26-1POS/results/OE26-1POS_cluster_summary.tsv"),
    "OE25-1-ALL":  os.path.join(FBMN_DIR, "batch_05_OE25-1-ALL/results/OE25-1-ALL_cluster_summary.tsv"),
    "ALL-25-2-POS": os.path.join(FBMN_DIR, "batch_06_ALL-25-2/results/ALL-25-2-POS_cluster_summary.tsv"),
}

# Metadata
METADATA_FILE = os.path.join(FBMN_DIR, "analysis15_master_metadata.tsv")

REF_BATCH = "OE23-POS"  # Reference batch for RT alignment
MZ_TOLERANCE_PPM = 5     # m/z matching tolerance (Orbitrap accuracy)
RT_TOLERANCE_ANCHOR = 2.0  # Wide RT tolerance for anchor discovery
RT_TOLERANCE_ALIGN = 0.5   # Tight RT for final alignment after LOESS
LOESS_FRAC = 0.3
DRIFT_THRESHOLD = 1.0
MIN_BATCHES_ANCHOR = 3
RT_MIN = 1.5              # Valid chromatographic window start (min)
RT_MAX = 25.0             # Valid chromatographic window end (min)


def get_sample_columns(df):
    """Extract sample columns (ending in .raw Peak area or .mzML Peak area)."""
    return [c for c in df.columns if c.endswith('.raw Peak area') or c.endswith('.mzML Peak area')]


def clean_sample_name(col_name):
    """Extract clean sample name from quantification column."""
    name = col_name.replace('.raw Peak area', '').replace('.mzML Peak area', '')
    return name.strip()


# ==========================================
# STEP 1: LOAD BATCHES & FIND ANCHORS
# ==========================================
def step1_find_anchors():
    """
    Load all 6 batches. For each batch, extract m/z and RT arrays.
    Find features that match across ≥3 batches by m/z and RT.
    Only features within the valid RT window [RT_MIN, RT_MAX] are considered.

    Strategy: Use reference batch features as candidate anchors,
    then check which ones appear in other batches.
    """
    print("\n" + "=" * 70)
    print("  STEP 1: CROSS-BATCH ANCHOR DISCOVERY")
    print("=" * 70)
    print(f"  m/z tolerance: {MZ_TOLERANCE_PPM} ppm")
    print(f"  RT tolerance: {RT_TOLERANCE_ANCHOR} min")
    print(f"  RT window: {RT_MIN} - {RT_MAX} min")
    print(f"  Min batches: {MIN_BATCHES_ANCHOR}")

    # Load all batch m/z and RT arrays (filtered to valid RT window)
    batch_data = {}
    for batch_name, path in BATCH_QUANT.items():
        if not os.path.exists(path):
            print(f"  {batch_name}: NOT FOUND at {path}")
            continue
        df = pd.read_csv(path, usecols=['row ID', 'row m/z', 'row retention time'])
        n_total = len(df)
        df = df[(df['row retention time'] >= RT_MIN) & (df['row retention time'] <= RT_MAX)]
        n_filtered = len(df)
        df = df.sort_values('row m/z').reset_index(drop=True)
        batch_data[batch_name] = {
            'ids': df['row ID'].values,
            'mzs': df['row m/z'].values,
            'rts': df['row retention time'].values,
            'n': len(df)
        }
        print(f"  {batch_name}: {n_filtered:,} features ({n_total - n_filtered:,} outside RT window removed)")
        del df
        gc.collect()

    # Use ALL features from ALL batches as anchor candidates
    # Group features across batches by m/z + RT matching
    print(f"\n  Building cross-batch feature groups...")

    # Strategy: iterate through reference batch features, find matches in all other batches
    ref = batch_data[REF_BATCH]
    other_batches = [b for b in batch_data.keys() if b != REF_BATCH]

    anchor_records = []

    for idx in range(ref['n']):
        if idx % 20000 == 0 and idx > 0:
            print(f"    Processed {idx:,} / {ref['n']:,} reference features...", flush=True)

        ref_mz = ref['mzs'][idx]
        ref_rt = ref['rts'][idx]
        ref_id = ref['ids'][idx]
        mz_tol = ref_mz * MZ_TOLERANCE_PPM / 1e6

        batch_matches = {REF_BATCH: (ref_id, ref_mz, ref_rt)}

        for batch_name in other_batches:
            bd = batch_data[batch_name]
            low = np.searchsorted(bd['mzs'], ref_mz - mz_tol, side='left')
            high = np.searchsorted(bd['mzs'], ref_mz + mz_tol, side='right')

            if low >= high:
                continue

            best_idx = None
            best_rt_diff = RT_TOLERANCE_ANCHOR + 1
            for i in range(low, high):
                rt_diff = abs(bd['rts'][i] - ref_rt)
                if rt_diff < best_rt_diff:
                    best_rt_diff = rt_diff
                    best_idx = i

            if best_idx is not None and best_rt_diff <= RT_TOLERANCE_ANCHOR:
                batch_matches[batch_name] = (
                    bd['ids'][best_idx],
                    bd['mzs'][best_idx],
                    bd['rts'][best_idx]
                )

        n_batches = len(batch_matches)
        if n_batches >= MIN_BATCHES_ANCHOR:
            record = {
                'ref_feature_id': ref_id,
                'ref_mz': ref_mz,
                'ref_rt': ref_rt,
                'n_batches': n_batches,
                'batches': ','.join(sorted(batch_matches.keys())),
            }
            # Store per-batch RTs for drift calculation
            for bn in batch_data.keys():
                if bn in batch_matches:
                    _, bmz, brt = batch_matches[bn]
                    record[f'rt_{bn}'] = brt
                    record[f'mz_{bn}'] = bmz
                    record[f'id_{bn}'] = batch_matches[bn][0]
                else:
                    record[f'rt_{bn}'] = None
                    record[f'mz_{bn}'] = None
                    record[f'id_{bn}'] = None

            anchor_records.append(record)

    # Also find anchors NOT in reference batch (features shared across ≥3 non-ref batches)
    print(f"\n  Searching for additional anchors outside reference batch...")
    ref_anchored_features = set()
    for rec in anchor_records:
        for bn in other_batches:
            fid = rec.get(f'id_{bn}')
            if fid is not None:
                ref_anchored_features.add((bn, fid))

    # For each non-ref batch, find features matching in ≥MIN_BATCHES_ANCHOR-1 other non-ref batches
    extra_anchors = 0
    for i, batch_a in enumerate(other_batches):
        bd_a = batch_data[batch_a]
        remaining = [b for b in other_batches if b != batch_a]

        for idx in range(bd_a['n']):
            if (batch_a, bd_a['ids'][idx]) in ref_anchored_features:
                continue

            a_mz = bd_a['mzs'][idx]
            a_rt = bd_a['rts'][idx]
            a_id = bd_a['ids'][idx]
            mz_tol = a_mz * MZ_TOLERANCE_PPM / 1e6

            batch_matches = {batch_a: (a_id, a_mz, a_rt)}

            for batch_b in remaining:
                bd_b = batch_data[batch_b]
                low = np.searchsorted(bd_b['mzs'], a_mz - mz_tol, side='left')
                high = np.searchsorted(bd_b['mzs'], a_mz + mz_tol, side='right')

                if low >= high:
                    continue

                best_idx = None
                best_rt_diff = RT_TOLERANCE_ANCHOR + 1
                for j in range(low, high):
                    rt_diff = abs(bd_b['rts'][j] - a_rt)
                    if rt_diff < best_rt_diff:
                        best_rt_diff = rt_diff
                        best_idx = j

                if best_idx is not None and best_rt_diff <= RT_TOLERANCE_ANCHOR:
                    batch_matches[batch_b] = (
                        bd_b['ids'][best_idx],
                        bd_b['mzs'][best_idx],
                        bd_b['rts'][best_idx]
                    )

            if len(batch_matches) >= MIN_BATCHES_ANCHOR:
                record = {
                    'ref_feature_id': a_id,
                    'ref_mz': a_mz,
                    'ref_rt': a_rt,
                    'n_batches': len(batch_matches),
                    'batches': ','.join(sorted(batch_matches.keys())),
                }
                for bn in batch_data.keys():
                    if bn in batch_matches:
                        _, bmz, brt = batch_matches[bn]
                        record[f'rt_{bn}'] = brt
                        record[f'mz_{bn}'] = bmz
                        record[f'id_{bn}'] = batch_matches[bn][0]
                    else:
                        record[f'rt_{bn}'] = None
                        record[f'mz_{bn}'] = None
                        record[f'id_{bn}'] = None
                anchor_records.append(record)
                extra_anchors += 1

        print(f"    {batch_a}: +{extra_anchors} extra anchors so far")

    anchors = pd.DataFrame(anchor_records)
    print(f"\n  === ANCHOR SUMMARY ===")
    print(f"  Total anchors: {len(anchors):,}")
    print(f"  From reference batch: {len(anchors) - extra_anchors:,}")
    print(f"  Additional (non-ref): {extra_anchors:,}")
    print(f"\n  Batch coverage:")
    for n in range(6, 0, -1):
        count = (anchors['n_batches'] >= n).sum()
        if count > 0:
            print(f"    In {n}+ batches: {count:,}")

    anchors.to_csv(os.path.join(ANCHOR_DIR, 'cross_batch_anchors.csv'), index=False)
    print(f"  Saved: cross_batch_anchors.csv")

    return batch_data, anchors


# ==========================================
# STEP 2: CALCULATE RT DRIFT
# ==========================================
def step2_rt_drift(anchors):
    """Calculate RT drift for each anchor relative to reference batch."""
    print("\n" + "=" * 70)
    print("  STEP 2: RT DRIFT CHARACTERIZATION")
    print("=" * 70)

    batches = list(BATCH_QUANT.keys())
    non_ref = [b for b in batches if b != REF_BATCH]

    # Calculate drift = batch_rt - reference_rt for each anchor
    for batch in non_ref:
        rt_col = f'rt_{batch}'
        ref_col = f'rt_{REF_BATCH}'
        drift_col = f'drift_{batch}'

        if rt_col in anchors.columns and ref_col in anchors.columns:
            anchors[drift_col] = anchors[rt_col] - anchors[ref_col]

    # Summary
    print(f"\n  RT Drift Summary (reference: {REF_BATCH}):")
    for batch in non_ref:
        drift_col = f'drift_{batch}'
        if drift_col in anchors.columns:
            valid = anchors[drift_col].dropna()
            if len(valid) > 0:
                print(f"    {batch}: n={len(valid):,}, "
                      f"mean={valid.mean():.3f}, std={valid.std():.3f}, "
                      f"median={valid.median():.3f} min")

    anchors.to_csv(os.path.join(ANCHOR_DIR, 'anchors_with_drift.csv'), index=False)
    return anchors


# ==========================================
# STEP 3: FIT LOESS CORRECTION
# ==========================================
def step3_loess_correction(anchors):
    """Fit LOESS curves for RT correction per batch. Only uses anchors within RT window."""
    print("\n" + "=" * 70)
    print("  STEP 3: LOESS RT CORRECTION")
    print("=" * 70)
    print(f"  RT window for LOESS: {RT_MIN} - {RT_MAX} min")

    non_ref = [b for b in BATCH_QUANT.keys() if b != REF_BATCH]

    # Filter to RT window + clean drift
    rt_mask = (anchors['ref_rt'] >= RT_MIN) & (anchors['ref_rt'] <= RT_MAX)
    mask = rt_mask.copy()
    for batch in non_ref:
        drift_col = f'drift_{batch}'
        if drift_col in anchors.columns:
            batch_ok = anchors[drift_col].isna() | (anchors[drift_col].abs() < DRIFT_THRESHOLD)
            mask = mask & batch_ok

    clean = anchors[mask].copy()
    n_rt_filtered = (~rt_mask).sum()
    print(f"  Anchors in RT window: {rt_mask.sum():,} ({n_rt_filtered:,} outside window removed)")
    print(f"  Clean anchors (|drift| < {DRIFT_THRESHOLD} min, in RT window): {len(clean):,}")

    loess_models = {}
    fit_stats = []

    for batch in non_ref:
        ref_rt_col = f'rt_{REF_BATCH}'
        drift_col = f'drift_{batch}'

        if drift_col not in clean.columns:
            continue

        valid = clean[ref_rt_col].notna() & clean[drift_col].notna()
        x = clean.loc[valid, ref_rt_col].values
        y = clean.loc[valid, drift_col].values

        if len(x) < 10:
            print(f"  {batch}: only {len(x)} points, skipping LOESS")
            continue

        try:
            frac = min(LOESS_FRAC, max(0.3, 30.0 / len(x)))
            smoothed = lowess(y, x, frac=frac, return_sorted=True)
            loess_x = smoothed[:, 0]
            loess_y = smoothed[:, 1]
            loess_models[batch] = (loess_x, loess_y)

            predicted = np.interp(x, loess_x, loess_y)
            residuals = y - predicted
            rmse = np.sqrt(np.mean(residuals**2))

            print(f"  {batch}: n={len(x):,}, RMSE={rmse:.4f} min, "
                  f"drift {y.mean():.3f} -> {residuals.mean():.4f}")

            fit_stats.append({
                'batch': batch,
                'n_anchors': len(x),
                'mean_drift_before': y.mean(),
                'std_drift_before': y.std(),
                'rmse_after_loess': rmse,
                'mean_residual': residuals.mean(),
                'std_residual': residuals.std(),
            })

            curve_df = pd.DataFrame({'reference_rt': loess_x, 'predicted_drift': loess_y})
            curve_df.to_csv(os.path.join(LOESS_DIR, f'loess_curve_{batch}.csv'), index=False)

        except Exception as e:
            print(f"  {batch}: LOESS failed - {e}")

    pd.DataFrame(fit_stats).to_csv(os.path.join(LOESS_DIR, 'loess_stats.csv'), index=False)

    # Create plots
    try:
        create_loess_plots(clean, loess_models, non_ref)
    except Exception as e:
        print(f"  Plot creation failed: {e}")

    return loess_models


def create_loess_plots(clean_df, loess_models, batches):
    """Create LOESS correction visualization."""
    n = len(batches)
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(6*ncols, 4*nrows))
    if nrows * ncols == 1:
        axes = np.array([axes])
    axes = axes.flatten()

    for i, batch in enumerate(batches):
        ax = axes[i]
        ref_col = f'rt_{REF_BATCH}'
        drift_col = f'drift_{batch}'

        if drift_col not in clean_df.columns:
            continue

        valid = clean_df[ref_col].notna() & clean_df[drift_col].notna()
        x = clean_df.loc[valid, ref_col].values
        y = clean_df.loc[valid, drift_col].values

        ax.scatter(x, y, alpha=0.15, s=5, c='#2563EB', label='Anchors')
        ax.axhline(y=0, color='gray', ls='--', alpha=0.5)

        if batch in loess_models:
            lx, ly = loess_models[batch]
            ax.plot(lx, ly, color='red', lw=2.5, label='LOESS fit')

        ax.set_xlabel('Reference RT (min)', fontsize=9)
        ax.set_ylabel('RT Drift (min)', fontsize=9)
        ax.set_title(f'{batch} (n={sum(valid):,})', fontsize=10, fontweight='bold')
        ax.set_ylim(-1.5, 1.5)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.2)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle(f'LOESS RT Correction (Ref: {REF_BATCH})', fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(os.path.join(FIGURE_DIR, 'loess_correction_plots.png'), dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Saved: loess_correction_plots.png")


# ==========================================
# STEP 4: BUILD UNIFIED CONSENSUS TABLE (MEMORY-OPTIMIZED)
# ==========================================
def step4_build_consensus(batch_data_dict, loess_models):
    """
    Build a unified consensus table — MEMORY-OPTIMIZED version.

    Processes one batch at a time to avoid loading all 6 DataFrames simultaneously.
    Uses numpy arrays instead of dict-per-row to minimize memory footprint.

    For each feature in the reference batch:
      1. Apply LOESS to find corresponding RT in other batches
      2. Match features across batches by m/z + corrected RT
      3. Merge sample intensities into one row

    Also include features unique to non-reference batches.
    """
    print("\n" + "=" * 70)
    print("  STEP 4: BUILD UNIFIED CONSENSUS TABLE (memory-optimized)")
    print("=" * 70)
    print(f"  m/z: {MZ_TOLERANCE_PPM} ppm | RT: {RT_TOLERANCE_ALIGN} min (post-LOESS)")

    all_batches = list(BATCH_QUANT.keys())
    non_ref = [b for b in all_batches if b != REF_BATCH]

    # ---- Phase A: Collect sample names from headers only ----
    print(f"\n  Phase A: Collecting sample names from all batches...")
    batch_sample_info = {}
    all_sample_names = []

    for batch_name, path in BATCH_QUANT.items():
        if not os.path.exists(path):
            continue
        df_head = pd.read_csv(path, nrows=0)
        sample_cols = get_sample_columns(df_head)
        sample_map = {col: clean_sample_name(col) for col in sample_cols}
        batch_sample_info[batch_name] = {
            'sample_cols': sample_cols,
            'sample_map': sample_map,
        }
        all_sample_names.extend(sample_map.values())
        print(f"    {batch_name}: {len(sample_cols)} samples")
        del df_head

    unique_samples = sorted(set(all_sample_names))
    sample_to_idx = {s: i for i, s in enumerate(unique_samples)}
    n_samples = len(unique_samples)
    print(f"  Total unique samples: {n_samples}")

    # ---- Phase B: Load reference batch, build intensity matrix ----
    print(f"\n  Phase B: Loading reference batch ({REF_BATCH}), RT window [{RT_MIN}, {RT_MAX}]...")
    ref_path = BATCH_QUANT[REF_BATCH]
    ref_df = pd.read_csv(ref_path)
    n_total_ref = len(ref_df)
    ref_df = ref_df[(ref_df['row retention time'] >= RT_MIN) & (ref_df['row retention time'] <= RT_MAX)]
    ref_df = ref_df.sort_values('row m/z').reset_index(drop=True)
    n_ref = len(ref_df)
    print(f"    {n_ref:,} features (removed {n_total_ref - n_ref:,} outside RT window)")

    ref_mzs = ref_df['row m/z'].values.copy()
    ref_rts = ref_df['row retention time'].values.copy()
    ref_ids = ref_df['row ID'].values.copy()

    # Initialize intensity matrix (float32 to save memory)
    # ~90K features × 263 samples × 4 bytes ≈ 94 MB
    intensities = np.zeros((n_ref, n_samples), dtype=np.float32)
    n_batches_arr = np.ones(n_ref, dtype=np.int8)
    batch_lists = [[REF_BATCH] for _ in range(n_ref)]
    mz_sums = ref_mzs.copy()   # For computing consensus m/z
    rt_sums = ref_rts.copy()    # For computing consensus RT

    # Fill reference batch intensities using numpy column extraction
    ref_scols = batch_sample_info[REF_BATCH]['sample_cols']
    ref_smap = batch_sample_info[REF_BATCH]['sample_map']
    for col in ref_scols:
        sidx = sample_to_idx[ref_smap[col]]
        vals = ref_df[col].values
        valid = np.isfinite(vals) & (vals > 0)
        intensities[valid, sidx] = vals[valid].astype(np.float32)

    # Free reference DataFrame (keep only arrays)
    del ref_df
    gc.collect()
    print(f"    Reference intensities loaded, DataFrame freed")

    # ---- Phase C: Match each non-ref batch one at a time ----
    matched_features = {bn: set() for bn in non_ref}

    for batch_name in non_ref:
        print(f"\n  Phase C: Matching {batch_name}...")
        path = BATCH_QUANT[batch_name]
        df = pd.read_csv(path)
        n_total_b = len(df)
        df = df[(df['row retention time'] >= RT_MIN) & (df['row retention time'] <= RT_MAX)]
        df = df.sort_values('row m/z').reset_index(drop=True)
        n_batch = len(df)
        print(f"    {n_batch:,} features (removed {n_total_b - n_batch:,} outside RT window)")

        bd_mzs = df['row m/z'].values
        bd_rts = df['row retention time'].values
        bd_ids = df['row ID'].values

        # Pre-extract sample columns as numpy arrays for fast access
        bd_scols = batch_sample_info[batch_name]['sample_cols']
        bd_smap = batch_sample_info[batch_name]['sample_map']
        bd_sample_arrays = {}
        for col in bd_scols:
            bd_sample_arrays[col] = df[col].values

        # Free the DataFrame, keep only numpy arrays
        del df
        gc.collect()

        # Get LOESS model for this batch
        has_loess = batch_name in loess_models
        if has_loess:
            lx, ly = loess_models[batch_name]

        n_matched = 0
        for idx in range(n_ref):
            if idx % 20000 == 0 and idx > 0:
                print(f"    {idx:,} / {n_ref:,} ({n_matched:,} matched)...", flush=True)

            rmz = ref_mzs[idx]
            rrt = ref_rts[idx]
            mz_tol = rmz * MZ_TOLERANCE_PPM / 1e6

            # LOESS-corrected search RT
            search_rt = rrt
            if has_loess:
                drift = np.interp(rrt, lx, ly, left=ly[0], right=ly[-1])
                search_rt = rrt + drift

            low = np.searchsorted(bd_mzs, rmz - mz_tol, side='left')
            high = np.searchsorted(bd_mzs, rmz + mz_tol, side='right')

            if low >= high:
                continue

            best_j = None
            best_rt_diff = RT_TOLERANCE_ALIGN + 1
            for j in range(low, high):
                rt_diff = abs(bd_rts[j] - search_rt)
                if rt_diff < best_rt_diff:
                    best_rt_diff = rt_diff
                    best_j = j

            if best_j is not None and best_rt_diff <= RT_TOLERANCE_ALIGN:
                n_batches_arr[idx] += 1
                batch_lists[idx].append(batch_name)
                mz_sums[idx] += bd_mzs[best_j]
                rt_sums[idx] += bd_rts[best_j]
                matched_features[batch_name].add(bd_ids[best_j])
                n_matched += 1

                # Copy intensities (take max if sample already has value)
                for col in bd_scols:
                    sidx = sample_to_idx[bd_smap[col]]
                    val = bd_sample_arrays[col][best_j]
                    if np.isfinite(val) and val > 0:
                        fval = np.float32(val)
                        if fval > intensities[idx, sidx]:
                            intensities[idx, sidx] = fval

        print(f"    {batch_name}: {n_matched:,} / {n_ref:,} ref features matched "
              f"({n_matched*100/n_ref:.1f}%)")

        del bd_mzs, bd_rts, bd_ids, bd_sample_arrays
        gc.collect()

    # ---- Phase D: Write reference-based features to CSV ----
    print(f"\n  Phase D: Writing {n_ref:,} reference-based features to CSV...")
    output_file = os.path.join(ALIGN_DIR, 'consensus_aligned_table.csv')

    # Build header
    meta_cols = ['feature_id', 'consensus_mz', 'consensus_rt', 'n_batches',
                 'batches', 'ref_batch_id', 'n_samples_detected']
    sample_col_names = [f"sample:{s}" for s in unique_samples]
    all_cols = meta_cols + sample_col_names

    # Track statistics
    total_features = 0
    batch_dist = defaultdict(int)

    with open(output_file, 'w') as fout:
        fout.write(','.join(all_cols) + '\n')

        for idx in range(n_ref):
            nb = int(n_batches_arr[idx])
            c_mz = mz_sums[idx] / nb
            c_rt = rt_sums[idx] / nb
            batches_str = ','.join(sorted(batch_lists[idx]))
            n_detected = int(np.sum(intensities[idx] > 0))

            parts = [
                f"A15_{ref_ids[idx]}",
                f"{c_mz:.6f}",
                f"{c_rt:.4f}",
                str(nb),
                f'"{batches_str}"',
                str(ref_ids[idx]),
                str(n_detected),
            ]
            # Append intensities
            for j in range(n_samples):
                v = intensities[idx, j]
                parts.append(f"{v:.0f}" if v > 0 else "0")

            fout.write(','.join(parts) + '\n')
            total_features += 1
            batch_dist[nb] += 1

            if idx % 20000 == 0 and idx > 0:
                print(f"    Written {idx:,} / {n_ref:,}...", flush=True)

    print(f"    {n_ref:,} reference features written")

    # Free large arrays
    del intensities, mz_sums, rt_sums, batch_lists, n_batches_arr
    gc.collect()

    # ---- Phase E: Add unmatched features from non-ref batches ----
    print(f"\n  Phase E: Adding unmatched features from non-ref batches...")

    with open(output_file, 'a') as fout:
        for batch_name in non_ref:
            print(f"    Loading {batch_name} for unmatched features...")
            path = BATCH_QUANT[batch_name]
            df = pd.read_csv(path)
            df = df[(df['row retention time'] >= RT_MIN) & (df['row retention time'] <= RT_MAX)]
            df = df.sort_values('row m/z').reset_index(drop=True)

            bd_ids = df['row ID'].values
            bd_mzs = df['row m/z'].values
            bd_rts = df['row retention time'].values

            bd_scols = batch_sample_info[batch_name]['sample_cols']
            bd_smap = batch_sample_info[batch_name]['sample_map']
            bd_sample_arrays = {col: df[col].values for col in bd_scols}

            has_loess = batch_name in loess_models
            if has_loess:
                lx, ly = loess_models[batch_name]

            matched_set = matched_features[batch_name]
            unmatched_count = 0

            for idx in range(len(df)):
                fid = bd_ids[idx]
                if fid in matched_set:
                    continue

                batch_mz = bd_mzs[idx]
                batch_rt = bd_rts[idx]

                # Correct RT back to reference frame
                corrected_rt = batch_rt
                if has_loess:
                    drift = np.interp(batch_rt, lx, ly, left=ly[0], right=ly[-1])
                    corrected_rt = batch_rt - drift

                # Build intensity row
                row_intensities = np.zeros(n_samples, dtype=np.float32)
                for col in bd_scols:
                    val = bd_sample_arrays[col][idx]
                    if np.isfinite(val) and val > 0:
                        row_intensities[sample_to_idx[bd_smap[col]]] = np.float32(val)

                n_detected = int(np.sum(row_intensities > 0))

                parts = [
                    f"A15_{batch_name}_{fid}",
                    f"{batch_mz:.6f}",
                    f"{corrected_rt:.4f}",
                    "1",
                    f'"{batch_name}"',
                    "",
                    str(n_detected),
                ]
                for j in range(n_samples):
                    v = row_intensities[j]
                    parts.append(f"{v:.0f}" if v > 0 else "0")

                fout.write(','.join(parts) + '\n')
                unmatched_count += 1
                total_features += 1
                batch_dist[1] += 1

            print(f"    {batch_name}: {unmatched_count:,} unique features added "
                  f"({len(df) - len(matched_set):,} unmatched)")

            del df, bd_ids, bd_mzs, bd_rts, bd_sample_arrays
            gc.collect()

    # ---- Summary ----
    print(f"\n  === CONSENSUS TABLE SUMMARY ===")
    print(f"  Total features: {total_features:,}")
    print(f"  Total samples: {n_samples}")
    print(f"\n  Cross-batch distribution:")
    for n in range(6, 0, -1):
        count = sum(v for k, v in batch_dist.items() if k >= n)
        if count > 0:
            print(f"    In {n}+ batches: {count:,}")
    print(f"    Single-batch only: {batch_dist.get(1, 0):,}")

    print(f"\n  Saved: {output_file}")
    return output_file


# ==========================================
# STEP 5: INTEGRATE LIBRARY MATCHES
# ==========================================
def step5_integrate_annotations(consensus_file):
    """Integrate library matches from all 6 FBMN runs into consensus table."""
    print("\n" + "=" * 70)
    print("  STEP 5: INTEGRATE LIBRARY MATCHES & FAMILIES")
    print("=" * 70)

    # Load all library matches
    all_lib = []
    for batch_name, path in BATCH_LIBRARY.items():
        if os.path.exists(path):
            lib = pd.read_csv(path, sep='\t')
            lib['source_batch'] = batch_name
            all_lib.append(lib)
            print(f"  {batch_name}: {len(lib)} library matches")

    if all_lib:
        combined_lib = pd.concat(all_lib, ignore_index=True)
        print(f"  Total combined: {len(combined_lib)} matches")
        combined_lib.to_csv(os.path.join(ALIGN_DIR, 'all_library_matches_combined.tsv'),
                           sep='\t', index=False)

    # Load cluster summaries for molecular family info
    all_clusters = []
    for batch_name, path in BATCH_CLUSTERS.items():
        if os.path.exists(path):
            cs = pd.read_csv(path, sep='\t')
            cs['source_batch'] = batch_name
            all_clusters.append(cs)
            n_families = cs['componentindex'].nunique() if 'componentindex' in cs.columns else 'N/A'
            print(f"  {batch_name}: {len(cs)} nodes, {n_families} families")

    if all_clusters:
        combined_clusters = pd.concat(all_clusters, ignore_index=True)
        combined_clusters.to_csv(os.path.join(ALIGN_DIR, 'all_cluster_summaries_combined.tsv'),
                                sep='\t', index=False)

    print(f"\n  Annotation integration complete.")
    print(f"  Library matches and cluster summaries saved for downstream analysis.")


# ==========================================
# MAIN
# ==========================================
def main():
    print("=" * 70)
    print("  ANALYSIS-15: CROSS-BATCH ALIGNMENT (6 FBMN BATCHES)")
    print("  Direct alignment without classical networking backbone")
    print("=" * 70)
    print(f"  Reference batch: {REF_BATCH}")
    print(f"  m/z tolerance: {MZ_TOLERANCE_PPM} ppm (Orbitrap)")
    print(f"  RT tolerance (anchors): {RT_TOLERANCE_ANCHOR} min")
    print(f"  RT tolerance (alignment): {RT_TOLERANCE_ALIGN} min")
    print(f"  RT window: {RT_MIN} - {RT_MAX} min")
    print(f"  LOESS fraction: {LOESS_FRAC}")
    print(f"  Min batches for anchor: {MIN_BATCHES_ANCHOR}")
    print("=" * 70)

    # Verify files
    print("\n  Checking input files...")
    for name, path in BATCH_QUANT.items():
        exists = os.path.exists(path)
        print(f"    {'OK' if exists else 'MISSING'}: {name}")

    # Check for resume
    anchor_file = os.path.join(ANCHOR_DIR, 'anchors_with_drift.csv')

    # Step 1 + 2: Anchors and drift
    if os.path.exists(anchor_file) and os.path.getsize(anchor_file) > 100:
        print(f"\n  [RESUME] Loading cached anchors...")
        anchors = pd.read_csv(anchor_file)
    else:
        _, anchors = step1_find_anchors()
        anchors = step2_rt_drift(anchors)

    # Step 3: LOESS
    loess_models = step3_loess_correction(anchors)

    # Free anchors from memory before Step 4
    del anchors
    gc.collect()

    # Step 4: Consensus table (memory-optimized, returns file path)
    consensus_file = step4_build_consensus(None, loess_models)

    # Step 5: Annotations
    step5_integrate_annotations(consensus_file)

    # Count consensus features for summary
    consensus_count = 0
    with open(consensus_file) as f:
        for _ in f:
            consensus_count += 1
    consensus_count -= 1  # Subtract header

    # Final summary
    print("\n" + "=" * 70)
    print("  PIPELINE COMPLETE!")
    print("=" * 70)
    print(f"  LOESS models: {len(loess_models)}")
    print(f"  Consensus features: {consensus_count:,}")
    print(f"\n  Outputs in: {ANALYSIS_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
