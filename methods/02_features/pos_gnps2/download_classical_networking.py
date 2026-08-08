#!/usr/bin/env python3
"""
Download Classical Molecular Networking results from GNPS2
Task: 4b456c4abae44df0b8fea89cd434efc4

Run on server: python3 download_classical_networking.py

NOTE: Classical networking uses different output paths than FBMN:
  - GraphML lives under nf_output/networking/ (not cytoscape/)
  - Components file is component_summary.tsv (not components.tsv)
  - No written_description.tsv exists in classical workflow
"""

import os, sys, requests

TASK = "4b456c4abae44df0b8fea89cd434efc4"
BASE = "https://gnps2.org/resultfile"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "classical_networking_results")
os.makedirs(OUT, exist_ok=True)

FILES = {
    # Clustering / consensus spectra
    "nf_output/clustering/specs_ms.mgf":                         "classical_consensus_spectra.mgf",
    "nf_output/clustering/clusterinfo.tsv":                      "classical_clusterinfo.tsv",
    "nf_output/clustering/clustersummary.tsv":                   "classical_clustersummary.tsv",
    "nf_output/clustering/featuretable_reformatted_presence.csv": "classical_featuretable_presence.csv",

    # Library matching
    "nf_output/library/merged_results_with_gnps.tsv":            "classical_library_matches.tsv",

    # Networking
    "nf_output/networking/filtered_pairs.tsv":                   "classical_filtered_pairs.tsv",
    "nf_output/networking/pairs_with_components.tsv":            "classical_pairs_with_components.tsv",
    "nf_output/networking/component_summary.tsv":                "classical_components.tsv",
    "nf_output/networking/clustersummary_with_network.tsv":      "classical_clustersummary_with_network.tsv",
    "nf_output/networking/clustersummary_with_groups.tsv":       "classical_clustersummary_with_groups.tsv",
    "nf_output/networking/network.graphml":                      "classical_network.graphml",
    "nf_output/networking/network_singletons.graphml":           "classical_network_singletons.graphml",

    # Misc
    "nf_output/modifinder_input.csv":                            "classical_modifinder_input.csv",
}

print(f"Downloading classical networking results")
print(f"Task: {TASK}")
print(f"Output: {OUT}")
print(f"Files: {len(FILES)}\n")

success = 0
failed = 0

for remote_path, local_name in FILES.items():
    local_path = os.path.join(OUT, local_name)
    url = f"{BASE}?task={TASK}&file={remote_path}"
    print(f"  {local_name}...", end=" ", flush=True)
    try:
        r = requests.get(url, timeout=600, stream=True)
        if r.status_code == 200 and len(r.content) > 100:
            with open(local_path, "wb") as f:
                f.write(r.content)
            size = len(r.content)
            if size > 1048576:
                print(f"OK ({size/1048576:.1f} MB)")
            else:
                print(f"OK ({size/1024:.1f} KB)")
            success += 1
        else:
            print(f"FAILED (HTTP {r.status_code}, {len(r.content)} bytes)")
            failed += 1
    except Exception as e:
        print(f"ERROR ({e})")
        failed += 1

print(f"\nDone: {success} downloaded, {failed} failed")
print(f"Results in: {OUT}")

# Quick validation
if success > 0:
    print("\n--- Quick validation ---")
    for fname in sorted(os.listdir(OUT)):
        fpath = os.path.join(OUT, fname)
        size = os.path.getsize(fpath)
        if fname.endswith(".tsv") or fname.endswith(".csv"):
            with open(fpath) as f:
                lines = sum(1 for _ in f)
            print(f"  {fname}: {size/1024:.1f} KB, {lines} lines")
        elif fname.endswith(".mgf"):
            with open(fpath) as f:
                spectra = sum(1 for line in f if line.startswith("BEGIN IONS"))
            print(f"  {fname}: {size/1048576:.1f} MB, {spectra} spectra")
        elif fname.endswith(".graphml"):
            print(f"  {fname}: {size/1048576:.1f} MB")
        else:
            print(f"  {fname}: {size/1024:.1f} KB")
