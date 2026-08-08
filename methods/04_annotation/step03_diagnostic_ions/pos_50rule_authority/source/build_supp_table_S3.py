#!/usr/bin/env python3
"""
Build Supplementary Table S3: Diagnostic MS2 ion database.

Flattens DIAGNOSTIC_FRAGMENTS, DIAGNOSTIC_NEUTRAL_LOSSES, and TERP_IONS
from scripts/ms2_classify_indval.py into a tidy table.

Columns:
  lipid_class, type (fragment|neutral_loss|terpenoid), ion_name,
  m_z_or_neutral_loss
"""
import os, sys
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(SCRIPT_DIR, '..'))
from config.paths import SUPPLEMENTARY

# Import the data structures from the classify script
sys.path.insert(0, SCRIPT_DIR)
from ms2_classify_indval import (
    DIAGNOSTIC_FRAGMENTS, DIAGNOSTIC_NEUTRAL_LOSSES, TERP_IONS,
)

rows = []
# Fragments are stored as {class: [(mz, name, ion_type), ...]}
for cls, ion_list in DIAGNOSTIC_FRAGMENTS.items():
    for mz, name, ion_type in ion_list:
        rows.append({
            'lipid_class': cls,
            'type': ion_type,  # 'fragment'
            'ion_name': name,
            'm_z_or_neutral_loss': mz,
        })

# Neutral losses: same shape
for cls, loss_list in DIAGNOSTIC_NEUTRAL_LOSSES.items():
    for mz, name in loss_list:
        rows.append({
            'lipid_class': cls,
            'type': 'neutral_loss',
            'ion_name': name,
            'm_z_or_neutral_loss': mz,
        })

# Terpenoid diagnostic ions — single global list, no class
for mz, name in TERP_IONS:
    rows.append({
        'lipid_class': 'Terpenoid',
        'type': 'terpenoid_fragment',
        'ion_name': name,
        'm_z_or_neutral_loss': mz,
    })

df = pd.DataFrame(rows)
by_type = df['type'].value_counts().to_dict()
print(f"Flattened {len(df)} diagnostic ion entries across "
      f"{df['lipid_class'].nunique()} classes")
print(f"By type: {by_type}")
print(df.to_string(index=False))

OUT_XLSX = os.path.join(SUPPLEMENTARY, 'tables',
                        'SupplementaryTable_S3_diagnostic_ions.xlsx')
OUT_CSV  = OUT_XLSX.replace('.xlsx', '.csv')
os.makedirs(os.path.dirname(OUT_XLSX), exist_ok=True)
df.to_excel(OUT_XLSX, index=False)
df.to_csv(OUT_CSV, index=False)
print(f"\nWrote {OUT_XLSX} and .csv")
