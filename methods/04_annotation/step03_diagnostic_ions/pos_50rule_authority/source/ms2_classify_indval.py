#!/usr/bin/env python3
"""
MS2 Diagnostic Ion Classification for IndVal-Only Unidentified Features
=======================================================================
Extracts MS2 spectra from batch MGF files for the 1,180 unidentified
IndVal consensus features and classifies them using the same 5% relative
intensity threshold diagnostic ion pipeline used for the composite atlas.

All IndVal consensus features are multi-batch with OE23-POS as reference,
so ref_batch_id maps directly to SCANS= in the OE23-POS MGF.

Output: indval_ms2_classification.csv (same format as platinum_ms2_reclassification_5pct.csv)
"""

import os
import sys
import re
import csv
from collections import defaultdict

import pandas as pd
import numpy as np

# ── Paths ──────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..', '..', '..'))

ATLAS_CSV = os.path.join(
    ANALYSIS_DIR, 'analysis-15', '04_biomarker_discovery',
    '04_platinum_diamond', 'atlas_expanded_final.csv'
)
CONSENSUS_CSV = os.path.join(
    ANALYSIS_DIR, 'analysis-15', '03_alignment', 'consensus_aligned_table.csv'
)

# MGF files — all 6 batch MGFs for fallback, but primary is OE23-POS
BATCH_MGFS = {
    'OE11-3-POS':    os.path.join(ANALYSIS_DIR, 'FBMN_all_batches_POS', 'batch_01_OE11-3-POS', 'OE11-3-POS_iimn_gnps.mgf'),
    'OE21-4-POS':    os.path.join(ANALYSIS_DIR, 'FBMN_all_batches_POS', 'batch_02_OE21-4-POS', 'OE21-4-POS_ALL_iimn_gnps.mgf'),
    'OE23-POS':      os.path.join(ANALYSIS_DIR, 'FBMN_all_batches_POS', 'batch_03_OE23-POS', 'OE23-POS_iimn_gnps_blanksub.mgf'),
    'OE26-1POS':     os.path.join(ANALYSIS_DIR, 'FBMN_all_batches_POS', 'batch_04_OE26-1POS', 'OE26-1POS_iimn_gnps.mgf'),
    'OE25-1-ALL':    os.path.join(ANALYSIS_DIR, 'FBMN_all_batches_POS', 'batch_05_OE25-1-ALL', 'OE25-1-ALL-POS_iimn_gnps.mgf'),
    'ALL-25-2-POS':  os.path.join(ANALYSIS_DIR, 'FBMN_all_batches_POS', 'batch_06_ALL-25-2', 'ALL-25-2-POS_iimn_gnps.mgf'),
}

OUT_DIR = os.path.join(ANALYSIS_DIR, 'analysis-17', 'positive', 'ms2_indval')
os.makedirs(OUT_DIR, exist_ok=True)
OUT_CSV = os.path.join(OUT_DIR, 'indval_ms2_classification.csv')

# ── Parameters ─────────────────────────────────────────────────────────────
MZ_TOL = 0.03          # 30 mDa tolerance for diagnostic ion matching
MIN_REL_INT = 5.0       # 5% of base peak minimum relative intensity
MIN_TERP_IONS = 3       # Minimum terpenoid ions for "possible terpenoid"
STRONG_TERP_IONS = 5    # Strong terpenoid evidence

# ── Diagnostic Ion Database (Positive Mode) ────────────────────────────────
# Format: (m/z, name, ion_type)  where ion_type = 'fragment' or 'neutral_loss'

DIAGNOSTIC_FRAGMENTS = {
    'Phospholipid': [
        (184.0733, 'Phosphocholine_C5H13NO4P', 'fragment'),
    ],
    'PE': [
        (196.0380, 'Glycerophosphoethanolamine', 'fragment'),
    ],
    'PG': [
        (171.0058, 'Glycerophosphoglycerol', 'fragment'),
        (152.9953, 'Glycerophosphoglycerol-H2O', 'fragment'),
    ],
    'PI': [
        (241.0119, 'Glycerophosphoinositol', 'fragment'),
    ],
    'Sphingolipid': [
        (264.2686, 'Sphingoid_d18:1', 'fragment'),
        (282.2791, 'Sphingoid_d18:1+H2O', 'fragment'),
        (252.2686, 'Sphingoid_d18:0', 'fragment'),
    ],
    'Betaine_Lipid': [
        (236.1492, 'DGTS_headgroup', 'fragment'),
    ],
    'Sterol': [
        (369.3516, 'Cholesterol_-H2O', 'fragment'),
        (383.3672, 'Campesterol_-H2O', 'fragment'),
        (395.3672, 'Stigmasterol_-H2O', 'fragment'),
        (397.3829, 'Sitosterol_-H2O', 'fragment'),
    ],
    'Archaeal_Lipid': [
        (653.5727, 'GDGT_core', 'fragment'),
        (743.6196, 'GDGT_fragment', 'fragment'),
        (373.3672, 'Archaeol_fragment', 'fragment'),
        (595.5310, 'Archaeol_core', 'fragment'),
        (669.5676, 'GDGT_variant', 'fragment'),
    ],
    'Quinone': [
        (197.0808, 'Ubiquinone_headgroup', 'fragment'),
    ],
    'Glycolipid': [
        (243.0863, 'Hexose_fragment', 'fragment'),
        (225.0757, 'Hexose-H2O', 'fragment'),
    ],
    'Hopanoid': [
        (191.1794, 'Hopane_skeleton', 'fragment'),
        (369.3516, 'Hopanol', 'fragment'),
        (395.3672, 'Diplopterol', 'fragment'),
    ],
}

# Neutral losses (computed from precursor - fragment)
DIAGNOSTIC_NEUTRAL_LOSSES = {
    'PE': [(141.0191, 'Ethanolamine_phosphate')],
    'PS': [(87.0320, 'Serine')],
    'Glycolipid': [(162.0528, 'Hexose_loss'), (180.0634, 'Hexose+H2O_loss')],
    'Sphingolipid': [(180.0634, 'Hexose_ceramide')],
}

# Terpenoid diagnostic ions (isoprene series + triterpenoid markers)
TERP_IONS = [
    (69.0699, 'Isoprene_unit'),
    (81.0699, 'Monoterpene_frag'),
    (93.0699, 'Toluenium'),
    (95.0855, 'Methylcyclohexenyl'),
    (105.0699, 'Xylene_cation'),
    (107.0855, 'Dimethylcyclohexenyl'),
    (109.1012, 'Sesquiterpene_frag'),
    (119.0855, 'C9H11+'),
    (121.1012, 'Sesquiterpene_frag2'),
    (123.1168, 'Sesquiterpene_frag3'),
    (133.1012, 'Diterpene_frag'),
    (135.1168, 'Diterpene_frag2'),
    (147.1168, 'C11H15+'),
    (149.1325, 'Sesterterpene_frag'),
    (161.1325, 'C12H17+'),
    (175.1481, 'C13H19+'),
    (189.1638, 'Hopane_C14'),
    (191.1794, 'Hopane_skeleton'),
    (203.1794, 'Triterpene_C15'),
    (205.1951, 'Steroid_C15'),
    (218.2035, 'RDA_triterpene'),
]


def parse_mgf_index(mgf_path):
    """Build scan→file-offset index for fast random access."""
    index = {}  # scan_id → (offset, pepmass)
    with open(mgf_path, 'r') as f:
        offset = 0
        current_scan = None
        current_pepmass = None
        for line in f:
            if line.startswith('BEGIN IONS'):
                offset_start = offset
            elif line.startswith('SCANS=') or line.startswith('FEATURE_ID='):
                val = line.strip().split('=', 1)[1]
                if line.startswith('SCANS='):
                    current_scan = int(val)
                elif current_scan is None and line.startswith('FEATURE_ID='):
                    current_scan = int(val)
            elif line.startswith('PEPMASS='):
                current_pepmass = float(line.strip().split('=', 1)[1].split()[0])
            elif line.startswith('END IONS'):
                if current_scan is not None:
                    index[current_scan] = (offset_start, current_pepmass)
                current_scan = None
                current_pepmass = None
            offset += len(line.encode('utf-8'))
    return index


def extract_spectrum(mgf_path, offset):
    """Extract peaks from MGF at given file offset. Returns list of (mz, intensity)."""
    peaks = []
    in_peaks = False
    with open(mgf_path, 'r') as f:
        f.seek(offset)
        for line in f:
            line = line.strip()
            if line == 'END IONS':
                break
            if in_peaks and line:
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        peaks.append((mz, intensity))
                    except ValueError:
                        pass
            elif line.startswith('Num peaks=') or line.startswith('NUM PEAKS='):
                in_peaks = True
            elif not line.startswith(('BEGIN', 'FEATURE', 'MSLEVEL', 'RT', 'PEPMASS',
                                      'CHARGE', 'SCANS', 'SPEC', 'MERGED', 'COLLISION',
                                      'FRAG', 'ISOL', 'ION', 'FILE', 'SOURCE')):
                # Might be peak data without explicit "Num peaks=" header
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        mz = float(parts[0])
                        intensity = float(parts[1])
                        peaks.append((mz, intensity))
                        in_peaks = True
                    except ValueError:
                        pass
    return peaks


def classify_spectrum(peaks, precursor_mz):
    """
    Classify a spectrum using diagnostic ion matching with 5% threshold.
    Returns dict with classification results.
    """
    if not peaks:
        return {
            'n_peaks_total': 0,
            'n_peaks_above_threshold': 0,
            'n_diagnostic_fragments': 0,
            'n_neutral_losses': 0,
            'n_diagnostic_total': 0,
            'assigned_class': 'Unknown',
            'assigned_subclass': 'Unknown',
            'classification_confidence': 'None',
            'classification_reasoning': 'No MS2 peaks',
            'diagnostic_ions_detail': '',
        }

    # Apply 5% relative intensity threshold
    max_int = max(p[1] for p in peaks)
    if max_int == 0:
        return {
            'n_peaks_total': len(peaks),
            'n_peaks_above_threshold': 0,
            'n_diagnostic_fragments': 0,
            'n_neutral_losses': 0,
            'n_diagnostic_total': 0,
            'assigned_class': 'Unknown',
            'assigned_subclass': 'Unknown',
            'classification_confidence': 'None',
            'classification_reasoning': 'All peaks zero intensity',
            'diagnostic_ions_detail': '',
        }

    filtered = [(mz, (intensity / max_int) * 100)
                for mz, intensity in peaks
                if (intensity / max_int) * 100 >= MIN_REL_INT]

    n_total = len(peaks)
    n_above = len(filtered)

    if not filtered:
        return {
            'n_peaks_total': n_total,
            'n_peaks_above_threshold': 0,
            'n_diagnostic_fragments': 0,
            'n_neutral_losses': 0,
            'n_diagnostic_total': 0,
            'assigned_class': 'Unknown',
            'assigned_subclass': 'Unknown',
            'classification_confidence': 'None',
            'classification_reasoning': 'No peaks above 5% threshold',
            'diagnostic_ions_detail': '',
        }

    filtered_mzs = np.array([p[0] for p in filtered])
    filtered_ints = {p[0]: p[1] for p in filtered}

    # Compute neutral losses from precursor
    neutral_losses = precursor_mz - filtered_mzs

    # ── Match diagnostic fragments ────────────────────────────────────
    matched_classes = defaultdict(list)  # class → list of (ion_name, mz, rel_int)
    detail_parts = []

    for lipid_class, ions in DIAGNOSTIC_FRAGMENTS.items():
        for ref_mz, name, ion_type in ions:
            diffs = np.abs(filtered_mzs - ref_mz)
            idx = np.argmin(diffs)
            if diffs[idx] <= MZ_TOL:
                rel_int = filtered_ints[filtered_mzs[idx]]
                matched_classes[lipid_class].append((name, filtered_mzs[idx], rel_int))
                detail_parts.append(f'fragment:{filtered_mzs[idx]:.4f}({name},{rel_int:.0f}%)')

    # ── Match diagnostic neutral losses ───────────────────────────────
    n_nl_matches = 0
    for lipid_class, losses in DIAGNOSTIC_NEUTRAL_LOSSES.items():
        for ref_loss, name in losses:
            diffs = np.abs(neutral_losses - ref_loss)
            idx = np.argmin(diffs)
            if diffs[idx] <= MZ_TOL:
                rel_int = filtered_ints[filtered_mzs[idx]]
                matched_classes[lipid_class].append((name + '_NL', filtered_mzs[idx], rel_int))
                detail_parts.append(f'loss:{ref_loss:.4f}({name},{rel_int:.0f}%)')
                n_nl_matches += 1

    # ── Count terpenoid ions ──────────────────────────────────────────
    n_terp = 0
    terp_details = []
    for ref_mz, name in TERP_IONS:
        diffs = np.abs(filtered_mzs - ref_mz)
        idx = np.argmin(diffs)
        if diffs[idx] <= MZ_TOL:
            rel_int = filtered_ints[filtered_mzs[idx]]
            n_terp += 1
            terp_details.append(f'fragment:{filtered_mzs[idx]:.4f}({name},{rel_int:.0f}%)')

    # ── Classification decision ───────────────────────────────────────
    n_frag_matches = sum(len(v) for v in matched_classes.values()) - n_nl_matches
    n_diag_total = n_frag_matches + n_nl_matches

    assigned_class = 'Unknown'
    assigned_subclass = 'Unknown'
    confidence = 'None'
    reasoning = ''

    if matched_classes:
        # Priority: strongest evidence first
        # Phospholipid (m/z 184) is strongest single diagnostic
        if 'Phospholipid' in matched_classes:
            pc_rel = matched_classes['Phospholipid'][0][2]
            assigned_class = 'Phospholipid'
            if 'PE' in matched_classes:
                assigned_subclass = 'PE'
                confidence = 'High'
                reasoning = f'PE headgroup + phosphocholine at {pc_rel:.0f}%'
            elif 'Sphingolipid' in matched_classes:
                assigned_subclass = 'SM'
                assigned_class = 'Sphingolipid'
                confidence = 'High'
                reasoning = f'Sphingoid base + phosphocholine → SM'
            else:
                assigned_subclass = 'PC/SM'
                confidence = 'High' if pc_rel >= 20 else 'Medium'
                reasoning = f'Phosphocholine at {pc_rel:.0f}% → PC or SM'
        elif 'PE' in matched_classes:
            assigned_class = 'Phospholipid'
            assigned_subclass = 'PE'
            confidence = 'Medium'
            reasoning = 'PE-specific ions without phosphocholine'
        elif 'PG' in matched_classes:
            assigned_class = 'Phospholipid'
            assigned_subclass = 'PG'
            confidence = 'Medium'
            reasoning = 'Glycerophosphoglycerol diagnostic ions'
        elif 'PI' in matched_classes:
            assigned_class = 'Phospholipid'
            assigned_subclass = 'PI'
            confidence = 'Medium'
            reasoning = 'Glycerophosphoinositol diagnostic ion'
        elif 'PS' in matched_classes:
            assigned_class = 'Phospholipid'
            assigned_subclass = 'PS'
            confidence = 'Medium'
            reasoning = 'Serine neutral loss'
        elif 'Sphingolipid' in matched_classes:
            assigned_class = 'Sphingolipid'
            assigned_subclass = 'Ceramide'
            confidence = 'Medium'
            reasoning = 'Sphingoid base fragments'
        elif 'Archaeal_Lipid' in matched_classes:
            assigned_class = 'Archaeal_Lipid'
            assigned_subclass = 'GDGT/Archaeol'
            confidence = 'High' if len(matched_classes['Archaeal_Lipid']) >= 2 else 'Medium'
            reasoning = f'{len(matched_classes["Archaeal_Lipid"])} archaeal lipid diagnostic fragments'
        elif 'Sterol' in matched_classes:
            assigned_class = 'Sterol'
            assigned_subclass = 'Sterol'
            confidence = 'Medium'
            reasoning = 'Sterol dehydration fragments'
        elif 'Betaine_Lipid' in matched_classes:
            assigned_class = 'Betaine_Lipid'
            assigned_subclass = 'DGTS'
            confidence = 'Medium'
            reasoning = 'DGTS headgroup fragment'
        elif 'Glycolipid' in matched_classes:
            assigned_class = 'Glycolipid'
            assigned_subclass = 'Glycolipid'
            confidence = 'Medium' if len(matched_classes['Glycolipid']) >= 2 else 'Low'
            reasoning = 'Hexose fragment/neutral loss'
        elif 'Quinone' in matched_classes:
            assigned_class = 'Quinone'
            assigned_subclass = 'Ubiquinone'
            confidence = 'Medium'
            reasoning = 'Ubiquinone headgroup fragment'
        elif 'Hopanoid' in matched_classes:
            assigned_class = 'Pentacyclic_Triterpenoid'
            assigned_subclass = 'Hopanoid'
            confidence = 'Medium' if len(matched_classes['Hopanoid']) >= 2 else 'Low'
            reasoning = 'Hopanoid diagnostic fragments'

    # If no class-specific match, check terpenoids
    if assigned_class == 'Unknown' and n_terp >= MIN_TERP_IONS:
        assigned_class = 'Terpenoid'
        assigned_subclass = 'Terpenoid'
        if n_terp >= STRONG_TERP_IONS:
            # Additional check: with many peaks, terp ions may be noise
            if n_above <= 50 or n_terp >= 7:
                confidence = 'Medium'
                reasoning = f'{n_terp} terpenoid fragments → terpenoid/isoprenoid'
            else:
                confidence = 'Low'
                reasoning = f'{n_terp} terpenoid fragments in complex spectrum ({n_above} peaks) → possible terpenoid'
        else:
            confidence = 'Low'
            reasoning = f'{n_terp} terpenoid fragments → possible terpenoid/isoprenoid'
        detail_parts.extend(terp_details)
        n_diag_total += n_terp
        n_frag_matches += n_terp

    detail_str = '; '.join(detail_parts)

    return {
        'n_peaks_total': n_total,
        'n_peaks_above_threshold': n_above,
        'n_diagnostic_fragments': n_frag_matches,
        'n_neutral_losses': n_nl_matches,
        'n_diagnostic_total': n_frag_matches + n_nl_matches,
        'assigned_class': assigned_class,
        'assigned_subclass': assigned_subclass,
        'classification_confidence': confidence,
        'classification_reasoning': reasoning,
        'diagnostic_ions_detail': detail_str,
    }


def main():
    # ── Load IndVal unidentified features ──────────────────────────────
    print("Loading expanded atlas...")
    atlas = pd.read_csv(ATLAS_CSV, usecols=['feature_id', 'discovery_method', 'annot_confidence', 'phylum'])
    targets = atlas[(atlas['discovery_method'] == 'indval_consensus') & (atlas['annot_confidence'] == 'Unidentified')]
    target_ids = set(targets['feature_id'].tolist())
    print(f"IndVal unidentified targets: {len(target_ids)}")

    # Also get ALL IndVal-only features for comprehensive coverage
    all_indval = atlas[atlas['discovery_method'] == 'indval_consensus']
    all_indval_ids = set(all_indval['feature_id'].tolist())
    print(f"All IndVal-only features: {len(all_indval_ids)}")

    # ── Load consensus table for ref_batch_id mapping ──────────────────
    print("Loading consensus table mapping...")
    ct = pd.read_csv(CONSENSUS_CSV, usecols=['feature_id', 'consensus_mz', 'ref_batch_id'])
    ct_map = ct.set_index('feature_id')

    # Map targets to scan IDs
    scan_to_feature = {}
    feature_to_mz = {}
    missing = 0
    for fid in all_indval_ids:
        if fid in ct_map.index:
            row = ct_map.loc[fid]
            scan_id = int(row['ref_batch_id'])
            scan_to_feature[scan_id] = fid
            feature_to_mz[fid] = float(row['consensus_mz'])
        else:
            missing += 1

    print(f"Mapped {len(scan_to_feature)} features to OE23-POS scans ({missing} not found)")

    # ── Build MGF index for OE23-POS ───────────────────────────────────
    mgf_path = BATCH_MGFS['OE23-POS']
    print(f"Building MGF index for OE23-POS ({os.path.basename(mgf_path)})...")
    mgf_index = parse_mgf_index(mgf_path)
    print(f"MGF index: {len(mgf_index)} spectra indexed")

    # Check how many targets have spectra
    found_scans = set(scan_to_feature.keys()) & set(mgf_index.keys())
    print(f"Features with MS2 spectra: {len(found_scans)} of {len(scan_to_feature)}")

    # ── Extract and classify ───────────────────────────────────────────
    print("Extracting and classifying spectra...")
    results = []
    n_classified = 0

    for scan_id in sorted(found_scans):
        fid = scan_to_feature[scan_id]
        offset, pepmass = mgf_index[scan_id]
        precursor_mz = feature_to_mz.get(fid, pepmass)

        # Extract peaks
        peaks = extract_spectrum(mgf_path, offset)

        # Classify
        result = classify_spectrum(peaks, precursor_mz)
        result['feature_id'] = fid
        result['precursor_mz'] = precursor_mz
        result['intensity_threshold_pct'] = MIN_REL_INT

        if result['assigned_class'] != 'Unknown':
            n_classified += 1

        results.append(result)

    # Also add entries for features without spectra
    for fid in all_indval_ids:
        if fid in ct_map.index:
            scan_id = int(ct_map.loc[fid, 'ref_batch_id'])
            if scan_id not in found_scans:
                results.append({
                    'feature_id': fid,
                    'precursor_mz': feature_to_mz.get(fid, 0),
                    'n_peaks_total': 0,
                    'n_peaks_above_threshold': 0,
                    'n_diagnostic_fragments': 0,
                    'n_neutral_losses': 0,
                    'n_diagnostic_total': 0,
                    'assigned_class': 'Unknown',
                    'assigned_subclass': 'Unknown',
                    'classification_confidence': 'None',
                    'classification_reasoning': 'No MS2 spectrum in OE23-POS MGF',
                    'diagnostic_ions_detail': '',
                    'intensity_threshold_pct': MIN_REL_INT,
                })

    # ── Save results ──────────────────────────────────────────────────
    df = pd.DataFrame(results)
    cols = ['feature_id', 'precursor_mz', 'n_peaks_total', 'n_peaks_above_threshold',
            'n_diagnostic_fragments', 'n_neutral_losses', 'n_diagnostic_total',
            'assigned_class', 'assigned_subclass', 'classification_confidence',
            'classification_reasoning', 'diagnostic_ions_detail', 'intensity_threshold_pct']
    df = df[cols]
    df.to_csv(OUT_CSV, index=False)
    print(f"\nSaved: {OUT_CSV}")

    # ── Summary ───────────────────────────────────────────────────────
    total = len(df)
    has_ms2 = (df['n_peaks_total'] > 0).sum()
    classified = (df['assigned_class'] != 'Unknown').sum()
    unid_classified = df[df['feature_id'].isin(target_ids) & (df['assigned_class'] != 'Unknown')]

    print(f"\n{'='*60}")
    print(f"MS2 CLASSIFICATION SUMMARY (all IndVal-only features)")
    print(f"{'='*60}")
    print(f"Total features:          {total}")
    print(f"With MS2 spectra:        {has_ms2}")
    print(f"Classified (non-Unknown): {classified} ({classified/total*100:.1f}%)")
    print(f"\nPreviously Unidentified → now classified: {len(unid_classified)}")
    print(f"\nClass distribution (classified features):")
    class_counts = df[df['assigned_class'] != 'Unknown']['assigned_class'].value_counts()
    for cls, n in class_counts.items():
        print(f"  {cls:<25s} {n:>5d}")

    print(f"\nConfidence distribution:")
    conf_counts = df[df['assigned_class'] != 'Unknown']['classification_confidence'].value_counts()
    for conf, n in conf_counts.items():
        print(f"  {conf:<10s} {n:>5d}")

    # Per-phylum breakdown for unidentified targets
    unid_merged = targets.merge(df, on='feature_id', how='left')
    print(f"\nPer-phylum upgrades (previously Unidentified → classified):")
    for phylum in sorted(unid_merged['phylum'].unique()):
        sub = unid_merged[unid_merged['phylum'] == phylum]
        classified_n = (sub['assigned_class'] != 'Unknown').sum()
        total_n = len(sub)
        print(f"  {phylum:<20s}  {classified_n:>4d} / {total_n:>4d}  ({classified_n/total_n*100:.1f}%)")

    return df


if __name__ == '__main__':
    main()
