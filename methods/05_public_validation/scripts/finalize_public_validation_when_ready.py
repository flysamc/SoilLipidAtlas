#!/usr/bin/env python3
"""Wait for both FASST queues, then run strict public-validation finalization."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path


def state_counts(status_dir: Path) -> Counter[str]:
    counts: Counter[str] = Counter()
    for path in status_dir.glob("*.json"):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            counts["unreadable"] += 1
            continue
        counts[str(payload.get("state", "missing"))] += 1
    return counts


def write_status(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pos-status-dir", required=True, type=Path)
    parser.add_argument("--neg-status-dir", required=True, type=Path)
    parser.add_argument("--expected-pos", type=int, default=7408)
    parser.add_argument("--expected-neg", type=int, default=5695)
    parser.add_argument("--poll-seconds", type=float, default=20.0)
    parser.add_argument("--status-file", required=True, type=Path)
    parser.add_argument("finalizer_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    finalizer_args = list(args.finalizer_args)
    if finalizer_args and finalizer_args[0] == "--":
        finalizer_args.pop(0)
    if not finalizer_args:
        raise ValueError("public_validation arguments are required after --")
    if args.poll_seconds <= 0:
        raise ValueError("poll-seconds must be positive")

    while True:
        pos = state_counts(args.pos_status_dir)
        neg = state_counts(args.neg_status_dir)
        snapshot: dict[str, object] = {
            "updated_utc": datetime.now(timezone.utc).isoformat(),
            "state": "waiting",
            "expected_success": {"POS": args.expected_pos, "NEG": args.expected_neg},
            "observed_states": {"POS": dict(pos), "NEG": dict(neg)},
        }
        ready = pos["success"] == args.expected_pos and neg["success"] == args.expected_neg
        if ready:
            snapshot["state"] = "finalizing"
            write_status(args.status_file, snapshot)
            completed = subprocess.run(
                [sys.executable, "-m", "paper2_repro.public_validation", *finalizer_args]
            )
            snapshot["completed_utc"] = datetime.now(timezone.utc).isoformat()
            snapshot["finalizer_exit_code"] = completed.returncode
            snapshot["state"] = "complete" if completed.returncode == 0 else "failed"
            write_status(args.status_file, snapshot)
            raise SystemExit(completed.returncode)
        write_status(args.status_file, snapshot)
        time.sleep(args.poll_seconds)


if __name__ == "__main__":
    main()
