#!/usr/bin/env python3
"""
Download all GNPS2 FBMN results for SOILMASS project.
Run this on a machine with internet access to gnps2.org.

Usage:
  python download_gnps2_results.py

  # Or specify a custom output directory:
  python download_gnps2_results.py --output /path/to/FBMN_all_batches_POS
"""

import os
import sys
import requests
from urllib.parse import quote
from pathlib import Path

# ============================================================
# Task IDs — all 6 batches completed on GNPS2
# ============================================================
BATCHES = {
    "batch_01_OE11-3-POS":  "57a91e57a53f4ab7ac5ccd69fbab2c4e",
    "batch_02_OE21-4-POS":  "8e0f0e314b814193b42bbcca4af86982",
    "batch_03_OE23-POS":    "bb4fd387ab584c32bae96e74a00b0b08",
    "batch_04_OE26-1POS":   "26e41456c73f457cbfd6d471b207f37a",
    "batch_05_OE25-1-ALL":  "5f9c17db47c94bcda93ea3cd12ac1b64",
    "batch_06_ALL-25-2":    "9f8b16574b08468f88acf4eca84d2805",
}

# ============================================================
# Result files to download from each FBMN task
# ============================================================
RESULT_FILES = {
    # Network results
    "cytoscape_network.graphml":    "nf_output/networking/network.graphml",
    "cluster_summary.tsv":          "nf_output/networking/clustersummary_with_network.tsv",

    # Library search results
    "library_matches.tsv":          "nf_output/library/merged_results_with_gnps.tsv",

    # Quantification
    "quantification_table.csv":     "nf_output/clustering/featuretable_reformated.csv",

    # Consensus spectra
    "consensus_spectra.mgf":        "nf_output/clustering/spectra_reformatted.mgf",

    # Metadata
    "merged_metadata.tsv":          "nf_output/metadata/merged_metadata.tsv",
}

GNPS2_URL = "https://gnps2.org"


def download_file(task_id, remote_path, local_path):
    """Download a single result file from GNPS2."""
    url = f"{GNPS2_URL}/resultfile?task={task_id}&file={quote(remote_path)}"

    try:
        resp = requests.get(url, timeout=600, stream=True)
        if resp.status_code == 200:
            total = int(resp.headers.get('content-length', 0))
            downloaded = 0
            with open(local_path, 'wb') as f:
                for chunk in resp.iter_content(chunk_size=65536):
                    f.write(chunk)
                    downloaded += len(chunk)

            size_mb = os.path.getsize(local_path) / (1024 * 1024)
            return True, f"{size_mb:.1f} MB"
        else:
            return False, f"HTTP {resp.status_code}"
    except Exception as e:
        return False, str(e)


def main():
    # Determine output directory
    if len(sys.argv) > 2 and sys.argv[1] == "--output":
        base_dir = Path(sys.argv[2])
    else:
        base_dir = Path(__file__).parent

    print("=" * 65)
    print("GNPS2 FBMN Results Downloader — SOILMASS Project")
    print("=" * 65)
    print(f"Output: {base_dir}")
    print(f"Batches: {len(BATCHES)}")
    print(f"Files per batch: {len(RESULT_FILES)}")
    print()

    total_ok = 0
    total_fail = 0

    for batch_name, task_id in BATCHES.items():
        print(f"\n{'─' * 65}")
        print(f"  {batch_name}")
        print(f"  Task: {task_id}")
        print(f"{'─' * 65}")

        results_dir = base_dir / batch_name / "results"
        results_dir.mkdir(parents=True, exist_ok=True)

        for local_name, remote_path in RESULT_FILES.items():
            local_path = results_dir / local_name

            # Skip if already downloaded
            if local_path.exists() and local_path.stat().st_size > 100:
                size_mb = local_path.stat().st_size / (1024 * 1024)
                print(f"    SKIP  {local_name} (already exists, {size_mb:.1f} MB)")
                total_ok += 1
                continue

            print(f"    GET   {local_name}...", end=" ", flush=True)
            ok, info = download_file(task_id, remote_path, str(local_path))

            if ok:
                print(f"OK ({info})")
                total_ok += 1
            else:
                print(f"FAIL ({info})")
                total_fail += 1
                # Remove partial file
                if local_path.exists():
                    local_path.unlink()

    print(f"\n{'=' * 65}")
    print(f"Download complete: {total_ok} succeeded, {total_fail} failed")
    print(f"Results in: {base_dir}/batch_*/results/")
    print(f"{'=' * 65}")


if __name__ == "__main__":
    main()
