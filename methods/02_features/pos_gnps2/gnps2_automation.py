#!/usr/bin/env python3
"""
GNPS2 FBMN Automation Script
=============================
Automates the end-to-end workflow for submitting FBMN jobs to GNPS2:
  1. Upload files via SFTP
  2. Submit FBMN workflow via GNPS2 API
  3. Monitor job status
  4. Download results when complete

Usage:
  # Step 1: Upload all batches and submit jobs
  python gnps2_automation.py submit --username YOUR_GNPS2_USERNAME --password YOUR_GNPS2_PASSWORD

  # Step 2: Check status of all submitted jobs
  python gnps2_automation.py status

  # Step 3: Download results for completed jobs
  python gnps2_automation.py download

  # Or do everything in one go (submit, wait, download):
  python gnps2_automation.py full --username YOUR_GNPS2_USERNAME --password YOUR_GNPS2_PASSWORD

Author: Automated for SOILMASS lipidomics project
Date: 2026-02-24
"""

import os
import sys
import json
import time
import argparse
import subprocess
from pathlib import Path
from datetime import datetime

try:
    import requests
except ImportError:
    print("Installing requests...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "requests", "--break-system-packages", "-q"])
    import requests

try:
    import paramiko
except ImportError:
    print("Installing paramiko for SFTP...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "paramiko", "--break-system-packages", "-q"])
    import paramiko

try:
    import pandas as pd
except ImportError:
    print("Installing pandas...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas", "--break-system-packages", "-q"])
    import pandas as pd


# =============================================================================
# Configuration
# =============================================================================

BASE_DIR = Path(__file__).parent
STATE_FILE = BASE_DIR / "gnps2_jobs_state.json"

# GNPS2 endpoints
GNPS2_URL = "https://gnps2.org"
GNPS2_SFTP_HOST = "sftp.gnps2.org"
GNPS2_SFTP_PORT = 443

# FBMN workflow parameters (matching your Orbitrap high-res data)
FBMN_PARAMS = {
    "workflowname": "feature_based_molecular_networking_workflow",
    "workflow_version": "release_35",
    # Tolerances for Orbitrap
    "tolerance.precursormztolerance": "0.02",
    "tolerance.mztolerance": "0.02",
    # Filtering
    "FILTER_PRECURSOR_WINDOW": "1",
    "WINDOW_FILTER": "1",
    # Networking
    "PAIRS_MIN_COSINE": "0.7",
    "MIN_MATCHED_PEAKS": "6",
    "MAX_SHIFT": "500",
    # Network topology
    "TOP_K": "10",
    "MAX_COMPONENT_SIZE": "100",
    # Library search
    "SCORE_THRESHOLD": "0.7",
    "MIN_MATCHED_PEAKS_SEARCH": "6",
    "ANALOG_SEARCH": "0",
    "MAX_SHIFT_MASS": "200",
    "TOP_K_RESULTS": "1",
    # Normalization
    "QUANT_TABLE_NORMALIZATION": "None",
}

# Batch definitions
BATCHES = [
    {
        "name": "batch_01_OE11-3-POS",
        "job_name": "SOILMASS_POS_FBMN_batch01_OE11-3-POS",
        "quant_csv": "OE11-3-POS_iimn_gnps_quant.csv",
        "mgf": "OE11-3-POS_iimn_gnps.mgf",
        "metadata": "gnps2_metadata.tsv",
    },
    {
        "name": "batch_02_OE21-4-POS",
        "job_name": "SOILMASS_POS_FBMN_batch02_OE21-4-POS",
        "quant_csv": "OE21-4-POS_ALL_iimn_gnps_quant.csv",
        "mgf": "OE21-4-POS_ALL_iimn_gnps.mgf",
        "metadata": "gnps2_metadata.tsv",
    },
    {
        "name": "batch_03_OE23-POS",
        "job_name": "SOILMASS_POS_FBMN_batch03_OE23-POS",
        "quant_csv": "OE23-POS_iimn_gnps_quant_blanksub.csv",
        "mgf": "OE23-POS_iimn_gnps_blanksub.mgf",
        "metadata": "gnps2_metadata.tsv",
    },
    {
        "name": "batch_04_OE26-1POS",
        "job_name": "SOILMASS_POS_FBMN_batch04_OE26-1POS",
        "quant_csv": "OE26-1POS_iimn_gnps_quant.csv",
        "mgf": "OE26-1POS_iimn_gnps.mgf",
        "metadata": "gnps2_metadata.tsv",
    },
    {
        "name": "batch_05_OE25-1-ALL",
        "job_name": "SOILMASS_POS_FBMN_batch05_OE25-1-ALL",
        "quant_csv": "OE25-1-ALL-POS_iimn_gnps_quant.csv",
        "mgf": "OE25-1-ALL-POS_iimn_gnps.mgf",
        "metadata": "gnps2_metadata.tsv",
    },
    {
        "name": "batch_06_ALL-25-2",
        "job_name": "SOILMASS_POS_FBMN_batch06_ALL-25-2",
        "quant_csv": "ALL-25-2-POS_iimn_gnps_quant.csv",
        "mgf": "ALL-25-2-POS_iimn_gnps.mgf",
        "metadata": "gnps2_metadata.tsv",
    },
]

# Result files to download from each completed FBMN job
RESULT_FILES = {
    # Network results
    "cytoscape_network.graphml": "nf_output/networking/network.graphml",
    "cluster_summary.tsv": "nf_output/networking/clustersummary_with_network.tsv",
    # Library search results
    "library_matches.tsv": "nf_output/library/merged_results_with_gnps.tsv",
    # Quantification
    "quantification_table.csv": "nf_output/clustering/featuretable_reformated.csv",
    # Spectra
    "consensus_spectra.mgf": "nf_output/clustering/spectra_reformatted.mgf",
    # Metadata
    "merged_metadata.tsv": "nf_output/metadata/merged_metadata.tsv",
}


# =============================================================================
# State Management
# =============================================================================

def load_state():
    """Load job submission state from JSON file."""
    if STATE_FILE.exists():
        with open(STATE_FILE) as f:
            return json.load(f)
    return {"jobs": {}, "last_updated": None}


def save_state(state):
    """Save job submission state to JSON file."""
    state["last_updated"] = datetime.now().isoformat()
    with open(STATE_FILE, "w") as f:
        json.dump(state, f, indent=2)


# =============================================================================
# SFTP Upload
# =============================================================================

def upload_files_sftp(username, password, batch):
    """Upload batch files to GNPS2 via SFTP."""
    batch_dir = BASE_DIR / batch["name"]
    remote_dir = f"/SOILMASS_FBMN/{batch['name']}"

    print(f"\n  Uploading {batch['name']} via SFTP...")

    transport = paramiko.Transport((GNPS2_SFTP_HOST, GNPS2_SFTP_PORT))
    transport.connect(username=username, password=password)
    sftp = paramiko.SFTPClient.from_transport(transport)

    # Create remote directory structure
    try:
        sftp.stat(f"/SOILMASS_FBMN")
    except FileNotFoundError:
        sftp.mkdir(f"/SOILMASS_FBMN")

    try:
        sftp.stat(remote_dir)
    except FileNotFoundError:
        sftp.mkdir(remote_dir)

    # Upload files
    files_to_upload = [batch["quant_csv"], batch["mgf"], batch["metadata"]]
    for filename in files_to_upload:
        local_path = str(batch_dir / filename)
        remote_path = f"{remote_dir}/{filename}"

        if not os.path.exists(local_path):
            print(f"    WARNING: {local_path} not found, skipping")
            continue

        file_size = os.path.getsize(local_path) / (1024 * 1024)
        print(f"    Uploading {filename} ({file_size:.1f} MB)...")
        sftp.put(local_path, remote_path)
        print(f"    Done: {filename}")

    sftp.close()
    transport.close()

    return remote_dir


# =============================================================================
# GNPS2 Workflow Submission
# =============================================================================

def login_gnps2(session, username, password):
    """Login to GNPS2 and return authenticated session."""
    print("  Logging in to GNPS2...")

    # Try the login endpoint
    login_url = f"{GNPS2_URL}/user/login"
    resp = session.post(login_url, data={
        "username": username,
        "password": password,
    }, allow_redirects=True)

    if resp.status_code == 200:
        print("  Login successful!")
        return True
    else:
        print(f"  Login returned status {resp.status_code}")
        # Try alternative login
        login_url2 = f"{GNPS2_URL}/api/user/login"
        resp2 = session.post(login_url2, json={
            "username": username,
            "password": password,
        })
        if resp2.status_code == 200:
            print("  Login successful (via API)!")
            return True

    print("  WARNING: Login may have failed. Continuing anyway...")
    return False


def submit_fbmn_job(session, username, batch, remote_dir):
    """Submit an FBMN workflow job to GNPS2."""
    print(f"\n  Submitting FBMN job: {batch['job_name']}...")

    # Build the submission payload
    params = dict(FBMN_PARAMS)
    params["desc"] = batch["job_name"]
    params["email"] = ""

    # File paths on GNPS2 server (relative to user's home)
    params["inputfeatures"] = f"SOILMASS_FBMN/{batch['name']}/{batch['quant_csv']}"
    params["inputspectra"] = f"SOILMASS_FBMN/{batch['name']}/{batch['mgf']}"
    params["inputmetadata"] = f"SOILMASS_FBMN/{batch['name']}/{batch['metadata']}"

    # Try the workflow invoke endpoint
    submit_url = f"{GNPS2_URL}/workflowinvoke"
    resp = session.post(submit_url, data=params, allow_redirects=True)

    if resp.status_code == 200:
        try:
            result = resp.json()
            task_id = result.get("task", result.get("task_id", ""))
            if task_id:
                print(f"  Job submitted! Task ID: {task_id}")
                return task_id
        except Exception:
            pass

        # Try to extract task ID from redirect URL
        if "task=" in resp.url:
            task_id = resp.url.split("task=")[-1].split("&")[0]
            print(f"  Job submitted! Task ID: {task_id}")
            return task_id

        # Check response text for task ID
        text = resp.text
        if "task" in text.lower():
            print(f"  Response received (check GNPS2 dashboard for task)")
            print(f"  Response URL: {resp.url}")
            return resp.url

    # Try alternative API endpoint
    submit_url2 = f"{GNPS2_URL}/api/workflow/invoke"
    resp2 = session.post(submit_url2, json=params)
    if resp2.status_code == 200:
        try:
            result = resp2.json()
            task_id = result.get("task", result.get("task_id", ""))
            if task_id:
                print(f"  Job submitted via API! Task ID: {task_id}")
                return task_id
        except Exception:
            pass

    print(f"  Submission response: {resp.status_code}")
    print(f"  You may need to submit manually at:")
    print(f"  {GNPS2_URL}/workflowinput?workflowname=feature_based_molecular_networking_workflow")
    return None


def check_job_status(task_id):
    """Check the status of a GNPS2 task."""
    if not task_id or task_id.startswith("http"):
        return "UNKNOWN"

    status_url = f"{GNPS2_URL}/taskstatus?task={task_id}"
    try:
        resp = requests.get(status_url, timeout=30)
        if resp.status_code == 200:
            try:
                data = resp.json()
                return data.get("status", "UNKNOWN")
            except Exception:
                if "DONE" in resp.text.upper():
                    return "DONE"
                elif "RUNNING" in resp.text.upper():
                    return "RUNNING"
                elif "FAILED" in resp.text.upper():
                    return "FAILED"
    except Exception as e:
        print(f"  Error checking status: {e}")

    return "UNKNOWN"


# =============================================================================
# Result Download
# =============================================================================

def download_results(task_id, batch_name, output_dir):
    """Download all result files from a completed GNPS2 FBMN task."""
    results_dir = output_dir / "results"
    results_dir.mkdir(exist_ok=True)

    print(f"\n  Downloading results for {batch_name} (task: {task_id})...")

    from urllib.parse import quote

    downloaded = []
    for local_name, remote_path in RESULT_FILES.items():
        url = f"{GNPS2_URL}/resultfile?task={task_id}&file={quote(remote_path)}"
        output_file = results_dir / local_name

        print(f"    Downloading {local_name}...")
        try:
            resp = requests.get(url, timeout=300, stream=True)
            if resp.status_code == 200:
                with open(output_file, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=8192):
                        f.write(chunk)
                size_mb = output_file.stat().st_size / (1024 * 1024)
                print(f"    Saved: {local_name} ({size_mb:.1f} MB)")
                downloaded.append(local_name)
            else:
                print(f"    Failed ({resp.status_code}): {local_name}")
        except Exception as e:
            print(f"    Error: {e}")

    # Also try to download the Cytoscape export with singletons (the main CSV)
    cytoscape_paths = [
        "nf_output/networking/network_nodes_with_singletons.tsv",
        "nf_output/networking/clustersummary_with_network.tsv",
    ]
    for cp in cytoscape_paths:
        url = f"{GNPS2_URL}/resultfile?task={task_id}&file={quote(cp)}"
        local_name = os.path.basename(cp)
        output_file = results_dir / local_name
        try:
            resp = requests.get(url, timeout=300, stream=True)
            if resp.status_code == 200:
                with open(output_file, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=8192):
                        f.write(chunk)
                print(f"    Saved: {local_name}")
                downloaded.append(local_name)
        except Exception:
            pass

    print(f"  Downloaded {len(downloaded)} files to {results_dir}/")
    return downloaded


# =============================================================================
# Main Commands
# =============================================================================

def cmd_submit(args):
    """Upload files and submit all 6 FBMN jobs."""
    state = load_state()
    session = requests.Session()

    print("=" * 60)
    print("GNPS2 FBMN Batch Submission")
    print("=" * 60)

    # Login
    login_gnps2(session, args.username, args.password)

    for batch in BATCHES:
        batch_dir = BASE_DIR / batch["name"]
        if not batch_dir.exists():
            print(f"\n  SKIP: {batch['name']} folder not found")
            continue

        print(f"\n{'=' * 60}")
        print(f"Processing: {batch['name']}")
        print(f"{'=' * 60}")

        # Upload via SFTP
        try:
            remote_dir = upload_files_sftp(args.username, args.password, batch)
        except Exception as e:
            print(f"  SFTP upload failed: {e}")
            print(f"  You can upload manually via: sftp -P 443 {args.username}@sftp.gnps2.org")
            remote_dir = f"SOILMASS_FBMN/{batch['name']}"

        # Submit job
        task_id = submit_fbmn_job(session, args.username, batch, remote_dir)

        # Save state
        state["jobs"][batch["name"]] = {
            "task_id": task_id,
            "job_name": batch["job_name"],
            "submitted": datetime.now().isoformat(),
            "status": "SUBMITTED" if task_id else "UPLOAD_ONLY",
            "remote_dir": remote_dir,
        }
        save_state(state)

    print(f"\n{'=' * 60}")
    print("Submission complete! State saved to: gnps2_jobs_state.json")
    print(f"{'=' * 60}")

    # Print summary
    print("\nJob Summary:")
    for name, info in state["jobs"].items():
        print(f"  {name}: {info.get('status', 'UNKNOWN')} (task: {info.get('task_id', 'N/A')})")

    print(f"\nNext steps:")
    print(f"  1. Check status:   python {__file__} status")
    print(f"  2. Download:       python {__file__} download")
    print(f"  3. Or check at:    {GNPS2_URL}/user/tasks")


def cmd_status(args):
    """Check status of all submitted jobs."""
    state = load_state()

    if not state["jobs"]:
        print("No jobs found. Run 'submit' first.")
        return

    print("=" * 60)
    print("GNPS2 Job Status Check")
    print(f"Last updated: {state.get('last_updated', 'N/A')}")
    print("=" * 60)

    all_done = True
    for name, info in state["jobs"].items():
        task_id = info.get("task_id")
        if task_id and not task_id.startswith("http"):
            status = check_job_status(task_id)
            info["status"] = status
            info["last_checked"] = datetime.now().isoformat()
        else:
            status = info.get("status", "UNKNOWN")

        emoji = {"DONE": "OK", "RUNNING": "..", "FAILED": "XX", "QUEUED": ">>"}
        marker = emoji.get(status, "??")
        print(f"  [{marker}] {name}: {status} (task: {task_id or 'N/A'})")

        if status not in ("DONE",):
            all_done = False

    save_state(state)

    if all_done:
        print(f"\nAll jobs complete! Run: python {__file__} download")
    else:
        print(f"\nSome jobs still running. Check again later.")
        print(f"  View at: {GNPS2_URL}/user/tasks")


def cmd_download(args):
    """Download results for all completed jobs."""
    state = load_state()

    if not state["jobs"]:
        print("No jobs found. Run 'submit' first.")
        return

    print("=" * 60)
    print("GNPS2 Results Download")
    print("=" * 60)

    for name, info in state["jobs"].items():
        task_id = info.get("task_id")

        # Check status first
        if task_id and not task_id.startswith("http"):
            status = check_job_status(task_id)
            info["status"] = status

            if status == "DONE":
                batch_dir = BASE_DIR / name
                downloaded = download_results(task_id, name, batch_dir)
                info["downloaded"] = downloaded
                info["download_time"] = datetime.now().isoformat()
            elif status == "RUNNING":
                print(f"\n  {name}: Still running, skipping download")
            elif status == "FAILED":
                print(f"\n  {name}: Job FAILED - check GNPS2 for details")
                print(f"    {GNPS2_URL}/status?task={task_id}")
            else:
                print(f"\n  {name}: Status {status}, skipping")
        else:
            print(f"\n  {name}: No valid task ID, skipping")

    save_state(state)
    print(f"\nDownload complete! Results saved to batch_*/results/ folders.")


def cmd_full(args):
    """Full pipeline: submit, wait, download."""
    # Submit
    cmd_submit(args)

    # Wait and poll
    print(f"\n{'=' * 60}")
    print("Waiting for jobs to complete (polling every 2 minutes)...")
    print("Press Ctrl+C to stop waiting (you can download later)")
    print(f"{'=' * 60}")

    try:
        while True:
            time.sleep(120)
            state = load_state()

            all_done = True
            for name, info in state["jobs"].items():
                task_id = info.get("task_id")
                if task_id and not task_id.startswith("http"):
                    status = check_job_status(task_id)
                    info["status"] = status
                    if status not in ("DONE", "FAILED"):
                        all_done = False

            save_state(state)

            # Print current status
            now = datetime.now().strftime("%H:%M:%S")
            statuses = [f"{n}: {i.get('status','?')}" for n, i in state["jobs"].items()]
            print(f"  [{now}] {' | '.join(statuses)}")

            if all_done:
                print("\nAll jobs finished!")
                break

    except KeyboardInterrupt:
        print(f"\nStopped waiting. Run 'download' later when jobs complete.")
        return

    # Download
    cmd_download(args)


def cmd_manual_add(args):
    """Manually add a task ID for a batch that was submitted via web UI."""
    state = load_state()

    batch_name = args.batch_name
    task_id = args.task_id

    # Find matching batch
    found = False
    for batch in BATCHES:
        if batch_name in batch["name"] or batch_name == batch["name"]:
            state["jobs"][batch["name"]] = {
                "task_id": task_id,
                "job_name": batch["job_name"],
                "submitted": datetime.now().isoformat(),
                "status": "SUBMITTED",
                "source": "manual",
            }
            save_state(state)
            print(f"Added task {task_id} for {batch['name']}")
            found = True
            break

    if not found:
        # Add as-is
        state["jobs"][batch_name] = {
            "task_id": task_id,
            "job_name": batch_name,
            "submitted": datetime.now().isoformat(),
            "status": "SUBMITTED",
            "source": "manual",
        }
        save_state(state)
        print(f"Added task {task_id} for {batch_name}")


# =============================================================================
# Entry Point
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="GNPS2 FBMN Batch Automation for SOILMASS Project",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Submit all 6 batches:
  python gnps2_automation.py submit --username myuser --password mypass

  # If you submitted manually via web, add task IDs:
  python gnps2_automation.py add-task batch_01 abc123def456
  python gnps2_automation.py add-task batch_02 xyz789ghi012

  # Check job statuses:
  python gnps2_automation.py status

  # Download results for completed jobs:
  python gnps2_automation.py download

  # Full automated pipeline (submit + wait + download):
  python gnps2_automation.py full --username myuser --password mypass
        """,
    )

    subparsers = parser.add_subparsers(dest="command", help="Command to run")

    # Submit command
    sub = subparsers.add_parser("submit", help="Upload files and submit FBMN jobs")
    sub.add_argument("--username", required=True, help="GNPS2 username")
    sub.add_argument("--password", required=True, help="GNPS2 password")

    # Status command
    subparsers.add_parser("status", help="Check status of submitted jobs")

    # Download command
    subparsers.add_parser("download", help="Download results for completed jobs")

    # Full pipeline
    sub = subparsers.add_parser("full", help="Full pipeline: submit, wait, download")
    sub.add_argument("--username", required=True, help="GNPS2 username")
    sub.add_argument("--password", required=True, help="GNPS2 password")

    # Manual task add
    sub = subparsers.add_parser("add-task", help="Manually add a task ID for web-submitted job")
    sub.add_argument("batch_name", help="Batch name (e.g., batch_01 or batch_01_OE11-3-POS)")
    sub.add_argument("task_id", help="GNPS2 task ID from the results URL")

    args = parser.parse_args()

    if args.command == "submit":
        cmd_submit(args)
    elif args.command == "status":
        cmd_status(args)
    elif args.command == "download":
        cmd_download(args)
    elif args.command == "full":
        cmd_full(args)
    elif args.command == "add-task":
        cmd_manual_add(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
