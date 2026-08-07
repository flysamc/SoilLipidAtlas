#!/usr/bin/env python3
"""
Step 02 v2: Iterative LOWESS RT correction for NEG mode.

Improvements over v1:
  - Uses statsmodels LOWESS (true locally-weighted regression) instead of polynomial on binned medians
  - 10 iterations with progressively tighter RT tolerance
  - Robust outlier rejection at each iteration (MAD-based)
  - Separate correction curve per batch, evaluated on ALL features

Algorithm per iteration:
  1. Match batch features to reference by m/z (5 ppm) within current RT tolerance
  2. Compute RT drift for each match
  3. Reject outliers (>3 MAD from median within RT bins)
  4. Fit LOWESS curve: drift = f(RT)
  5. Apply correction to ALL batch features (not just matched ones)
  6. Tighten RT tolerance for next iteration

Usage:
  python rt_correction_v2.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from pathlib import Path
import time
import warnings
warnings.filterwarnings('ignore')

# ── Config ──────────────────────────────────────────────────────
BASE = Path(__file__).parent.parent
QUANT_DIR = BASE.parent.parent / "FBMN_all_batches_NEG"
OUT_DIR = Path(__file__).parent
FIG_DIR = OUT_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

PPM_TOL = 5.0
RT_MIN = 1.5
RT_MAX = 21.0
N_ITERATIONS = 10

# RT tolerance schedule: starts wide, tightens each iteration
# More gradual tightening than v1
RT_TOL_SCHEDULE = [3.0, 2.5, 2.0, 1.5, 1.2, 1.0, 0.8, 0.7, 0.6, 0.5]

# LOWESS fraction — controls smoothness. Lower = more local flexibility
LOWESS_FRAC = 0.15  # 15% of data used for each local regression

REF_BATCH = "OE23-NEG"
BATCHES_TO_CORRECT = ["OE11-3-NEG", "OE21-4-NEG", "ALL-25-2-NEG", "OE26-1-NEG", "OE25-1-NEG"]

# Map batch short names to quant CSV directories
BATCH_DIRS = {
    "OE23-NEG":      "batch_03_OE23_NEG",
    "OE21-4-NEG":    "batch_02_OE21-4_NEG",
    "ALL-25-2-NEG":  "batch_06_ALL-25-2_NEG",
    "OE26-1-NEG":    "batch_04_OE26-1_NEG",
    "OE11-3-NEG":    "batch_01_OE11-3_NEG",
    "OE25-1-NEG":    "batch_05_OE25-1_NEG",
}


def load_quant_table(batch_name):
    """Load the original (blank-subtracted) quantification CSV for a batch."""
    batch_dir = BATCH_DIRS[batch_name]
    quant_dir = QUANT_DIR / batch_dir
    # Find the quant CSV (blanksub version preferred)
    candidates = list(quant_dir.glob("*_quant_blanksub.csv"))
    if not candidates:
        candidates = list(quant_dir.glob("*blanksub*quant*")) + list(quant_dir.glob("*quant*blanksub*"))
    if not candidates:
        candidates = list(quant_dir.glob("*_quant.csv"))
    if not candidates:
        raise FileNotFoundError(f"No quant CSV found in {quant_dir}")

    path = candidates[0]
    df = pd.read_csv(path)

    # Filter to RT window
    rt_mask = (df['row retention time'] >= RT_MIN) & (df['row retention time'] <= RT_MAX)
    df = df[rt_mask].reset_index(drop=True)
    return df


def match_features(batch_mz, batch_rt, ref_mz_sorted, ref_rt_sorted, ref_sort_idx,
                   ppm_tol, rt_tol):
    """Match batch features to reference features by m/z + RT.

    Returns arrays of (batch_idx, ref_idx, drift) for matches.
    """
    batch_indices = []
    ref_indices = []
    drifts = []

    for i in range(len(batch_mz)):
        mz_i = batch_mz[i]
        rt_i = batch_rt[i]

        mz_lo = mz_i * (1 - ppm_tol / 1e6)
        mz_hi = mz_i * (1 + ppm_tol / 1e6)
        idx_lo = np.searchsorted(ref_mz_sorted, mz_lo, side='left')
        idx_hi = np.searchsorted(ref_mz_sorted, mz_hi, side='right')

        if idx_lo >= idx_hi:
            continue

        # Among m/z candidates, find closest RT
        rt_diffs = ref_rt_sorted[idx_lo:idx_hi] - rt_i
        abs_rt = np.abs(rt_diffs)
        best_local = np.argmin(abs_rt)

        if abs_rt[best_local] <= rt_tol:
            batch_indices.append(i)
            ref_indices.append(ref_sort_idx[idx_lo + best_local])
            drifts.append(rt_diffs[best_local])

    return (np.array(batch_indices, dtype=np.int64),
            np.array(ref_indices, dtype=np.int64),
            np.array(drifts, dtype=np.float64))


def reject_outliers_mad(rt_values, drift_values, n_bins=50, mad_threshold=3.0):
    """Reject outlier drifts using MAD within RT bins."""
    keep = np.ones(len(rt_values), dtype=bool)

    rt_bins = np.linspace(rt_values.min(), rt_values.max(), n_bins + 1)
    for b in range(n_bins):
        mask = (rt_values >= rt_bins[b]) & (rt_values < rt_bins[b + 1])
        if mask.sum() < 5:
            continue
        bin_drift = drift_values[mask]
        med = np.median(bin_drift)
        mad = np.median(np.abs(bin_drift - med))
        if mad < 1e-6:
            mad = 1e-6
        outlier = np.abs(bin_drift - med) > mad_threshold * 1.4826 * mad
        idx = np.where(mask)[0]
        keep[idx[outlier]] = False

    return keep


def fit_lowess_correction(rt_matched, drift_matched, rt_all, frac=LOWESS_FRAC):
    """Fit LOWESS to matched drift data, predict correction for all features."""
    # Sort by RT for LOWESS
    sort_idx = np.argsort(rt_matched)
    rt_s = rt_matched[sort_idx]
    drift_s = drift_matched[sort_idx]

    # Fit LOWESS — returns (rt, smoothed_drift) pairs
    result = lowess(drift_s, rt_s, frac=frac, it=3, return_sorted=True)
    lowess_rt = result[:, 0]
    lowess_drift = result[:, 1]

    # Interpolate to predict correction for all features
    correction = np.interp(rt_all, lowess_rt, lowess_drift)
    return correction, lowess_rt, lowess_drift


def correct_batch(batch_name, ref_df, iteration_stats):
    """Run iterative LOWESS correction for one batch."""
    batch_df = load_quant_table(batch_name)
    n_features = len(batch_df)

    # Working copy of RT values (updated each iteration)
    batch_rt = batch_df['row retention time'].values.copy()
    batch_mz = batch_df['row m/z'].values

    # Reference features (sorted by m/z for binary search)
    ref_mz = ref_df['row m/z'].values
    ref_rt = ref_df['row retention time'].values
    ref_sort = np.argsort(ref_mz)
    ref_mz_s = ref_mz[ref_sort]
    ref_rt_s = ref_rt[ref_sort]

    total_correction = np.zeros(n_features)
    iter_iqrs = []

    print(f"\n  {batch_name} ({n_features} features)")

    for it in range(N_ITERATIONS):
        rt_tol = RT_TOL_SCHEDULE[it]

        # Match
        b_idx, r_idx, drifts = match_features(
            batch_mz, batch_rt, ref_mz_s, ref_rt_s, ref_sort, PPM_TOL, rt_tol)

        if len(drifts) < 50:
            print(f"    Iter {it+1}: only {len(drifts)} matches at RT_tol={rt_tol:.1f}, stopping")
            break

        # Outlier rejection
        keep = reject_outliers_mad(batch_rt[b_idx], drifts)
        b_idx_clean = b_idx[keep]
        drifts_clean = drifts[keep]
        n_rejected = (~keep).sum()

        # LOWESS fit
        correction, lw_rt, lw_drift = fit_lowess_correction(
            batch_rt[b_idx_clean], drifts_clean, batch_rt)

        # Apply correction
        batch_rt += correction
        total_correction += correction

        # Measure residual
        _, _, residuals = match_features(
            batch_mz, batch_rt, ref_mz_s, ref_rt_s, ref_sort, PPM_TOL, rt_tol)
        if len(residuals) > 0:
            q25, q75 = np.percentile(residuals, [25, 75])
            iqr = q75 - q25
            med = np.median(residuals)
        else:
            iqr = float('nan')
            med = float('nan')

        iter_iqrs.append(iqr)
        print(f"    Iter {it+1}: RT_tol={rt_tol:.1f}, matches={len(drifts):,}, "
              f"outliers={n_rejected}, residual IQR={iqr:.3f}, median={med:.4f}")

    # Store results
    batch_df['rt_corrected'] = batch_rt
    batch_df['rt_correction'] = total_correction
    iteration_stats[batch_name] = iter_iqrs

    return batch_df


def make_figures(all_results, ref_df, iteration_stats):
    """Generate comprehensive QC figures."""

    ref_mz = ref_df['row m/z'].values
    ref_rt = ref_df['row retention time'].values
    ref_sort = np.argsort(ref_mz)
    ref_mz_s = ref_mz[ref_sort]
    ref_rt_s = ref_rt[ref_sort]

    colors = {
        'OE11-3-NEG': '#0072B2',
        'OE21-4-NEG': '#E69F00',
        'OE26-1-NEG': '#D55E00',
        'OE25-1-NEG': '#CC79A7',
        'ALL-25-2-NEG': '#56B4E9',
    }

    # ── 1. Iteration convergence ──
    fig, ax = plt.subplots(figsize=(10, 6))
    for batch_name, iqrs in iteration_stats.items():
        iters = range(1, len(iqrs) + 1)
        ax.plot(iters, iqrs, 'o-', label=batch_name, color=colors.get(batch_name, 'gray'), linewidth=2)
    ax.axhline(0.5, color='green', ls='--', alpha=0.7, label='Target (0.5 min)')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Drift IQR (min)')
    ax.set_title(f'Iterative LOWESS Convergence (5 ppm, {N_ITERATIONS} iterations)')
    ax.legend()
    ax.set_xticks(range(1, N_ITERATIONS + 1))
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'iteration_convergence_v2.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ── 2. Residual drift histograms ──
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    axes_flat = axes.flatten()
    for idx, batch_name in enumerate(BATCHES_TO_CORRECT):
        ax = axes_flat[idx]
        batch_df = all_results[batch_name]
        batch_mz = batch_df['row m/z'].values
        batch_rt = batch_df['rt_corrected'].values

        _, _, residuals = match_features(
            batch_mz, batch_rt, ref_mz_s, ref_rt_s, ref_sort, PPM_TOL, 2.0)

        if len(residuals) > 0:
            q25, q50, q75 = np.percentile(residuals, [25, 50, 75])
            iqr = q75 - q25
            ax.hist(residuals, bins=200, range=(-3, 3), color=colors.get(batch_name, 'gray'),
                    alpha=0.7, edgecolor='none')
            ax.axvline(q50, color='red', ls='--', linewidth=1.5)
            ax.set_title(f'{batch_name}\nmed={q50:.3f}, IQR={iqr:.2f}')
        else:
            ax.set_title(f'{batch_name}\nno matches')
        ax.set_xlabel('ΔRT (min)')
        ax.set_ylabel('Count')

    axes_flat[-1].axis('off')
    fig.suptitle(f'Post-LOWESS v2 Residual Drift ({N_ITERATIONS} iterations, 5 ppm)', fontsize=13)
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'residual_drift_histograms_v2.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    # ── 3. Before vs After scatter (systematic drift curves) ──
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

    for batch_name in BATCHES_TO_CORRECT:
        batch_df = all_results[batch_name]
        batch_mz = batch_df['row m/z'].values
        orig_rt = batch_df['row retention time'].values
        corr_rt = batch_df['rt_corrected'].values

        c = colors.get(batch_name, 'gray')

        # Before: match uncorrected to ref
        b_idx_before, _, drifts_before = match_features(
            batch_mz, orig_rt, ref_mz_s, ref_rt_s, ref_sort, PPM_TOL, 5.0)

        if len(drifts_before) > 100:
            # Bin and plot running median
            rt_vals = orig_rt[b_idx_before]
            bins = np.linspace(RT_MIN, RT_MAX, 80)
            bin_idx = np.digitize(rt_vals, bins)
            bin_med = [np.median(drifts_before[bin_idx == b]) for b in range(1, len(bins))
                       if np.sum(bin_idx == b) >= 5]
            bin_centers = [(bins[b-1] + bins[b]) / 2 for b in range(1, len(bins))
                          if np.sum(bin_idx == b) >= 5]
            pre_iqr = np.percentile(drifts_before, 75) - np.percentile(drifts_before, 25)
            ax1.plot(bin_centers, bin_med, '-', color=c, linewidth=1.5,
                     label=f'{batch_name} (IQR={pre_iqr:.2f})')

        # After: match corrected to ref
        b_idx_after, _, drifts_after = match_features(
            batch_mz, corr_rt, ref_mz_s, ref_rt_s, ref_sort, PPM_TOL, 3.0)

        if len(drifts_after) > 100:
            rt_vals = corr_rt[b_idx_after]
            bin_idx = np.digitize(rt_vals, bins)
            bin_med = [np.median(drifts_after[bin_idx == b]) for b in range(1, len(bins))
                       if np.sum(bin_idx == b) >= 5]
            bin_centers = [(bins[b-1] + bins[b]) / 2 for b in range(1, len(bins))
                          if np.sum(bin_idx == b) >= 5]
            post_iqr = np.percentile(drifts_after, 75) - np.percentile(drifts_after, 25)
            ax2.plot(bin_centers, bin_med, '-', color=c, linewidth=1.5,
                     label=f'{batch_name} (IQR={post_iqr:.2f})')

    for ax, title in [(ax1, 'BEFORE'), (ax2, 'AFTER')]:
        ax.axhline(0, color='black', ls='-', linewidth=0.5)
        ax.set_xlabel('RT (min)')
        ax.set_ylabel('Drift (min)')
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.set_ylim(-3, 3)
        ax.grid(True, alpha=0.3)

    fig.suptitle(f'RT Drift Before vs After LOWESS v2 (5 ppm, {N_ITERATIONS} iter)', fontsize=13)
    plt.tight_layout()
    fig.savefig(FIG_DIR / 'rt_correction_before_after_v2.png', dpi=150, bbox_inches='tight')
    plt.close(fig)

    print("  Figures saved.")


def main():
    t_start = time.time()
    print("=" * 65)
    print(f"Step 02 v2: Iterative LOWESS RT Correction — NEG Mode")
    print("=" * 65)
    print(f"  m/z tolerance: {PPM_TOL} ppm")
    print(f"  RT window: [{RT_MIN}, {RT_MAX}] min")
    print(f"  Iterations: {N_ITERATIONS}")
    print(f"  RT tolerance schedule: {RT_TOL_SCHEDULE}")
    print(f"  LOWESS fraction: {LOWESS_FRAC}")
    print(f"  Reference: {REF_BATCH}")

    # Load reference batch
    ref_df = load_quant_table(REF_BATCH)
    print(f"\n  Reference: {REF_BATCH} — {len(ref_df)} features")

    # Correct each batch
    all_results = {}
    iteration_stats = {}

    for batch_name in BATCHES_TO_CORRECT:
        corrected = correct_batch(batch_name, ref_df, iteration_stats)
        all_results[batch_name] = corrected

    # Save reference (correction = 0)
    ref_out = ref_df.copy()
    ref_out['rt_corrected'] = ref_out['row retention time']
    ref_out['rt_correction'] = 0.0
    ref_out.to_csv(OUT_DIR / f"{REF_BATCH}_rt_corrected.csv", index=False)
    all_results[REF_BATCH] = ref_out

    # Save corrected batches
    for batch_name in BATCHES_TO_CORRECT:
        out_path = OUT_DIR / f"{batch_name}_rt_corrected.csv"
        all_results[batch_name].to_csv(out_path, index=False)

    # Final summary
    print(f"\n{'=' * 65}")
    print("FINAL POST-CORRECTION SUMMARY")
    print(f"{'=' * 65}")

    ref_mz = ref_df['row m/z'].values
    ref_rt = ref_df['row retention time'].values
    ref_sort = np.argsort(ref_mz)
    ref_mz_s = ref_mz[ref_sort]
    ref_rt_s = ref_rt[ref_sort]

    print(f"{'Batch':<18} {'Features':<10} {'Matches':<10} {'Median':<10} {'IQR':<10} {'±0.5min':<10} {'±1.0min':<10}")
    print("-" * 68)

    for batch_name in BATCHES_TO_CORRECT:
        batch_df = all_results[batch_name]
        batch_mz = batch_df['row m/z'].values
        batch_rt = batch_df['rt_corrected'].values

        _, _, residuals = match_features(
            batch_mz, batch_rt, ref_mz_s, ref_rt_s, ref_sort, PPM_TOL, 3.0)

        if len(residuals) > 0:
            q25, q50, q75 = np.percentile(residuals, [25, 50, 75])
            iqr = q75 - q25
            w05 = 100 * np.mean(np.abs(residuals) <= 0.5)
            w10 = 100 * np.mean(np.abs(residuals) <= 1.0)
        else:
            q50 = iqr = w05 = w10 = float('nan')

        print(f"  {batch_name:<16} {len(batch_df):<10} {len(residuals):<10} "
              f"{q50:<10.4f} {iqr:<10.3f} {w05:<10.1f}% {w10:<10.1f}%")

    elapsed = time.time() - t_start
    print(f"\n  Total time: {elapsed:.1f}s")

    # Generate figures
    print("\nGenerating QC figures...")
    make_figures(all_results, ref_df, iteration_stats)

    print(f"\n{'=' * 65}")
    print("Step 02 v2 complete.")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
