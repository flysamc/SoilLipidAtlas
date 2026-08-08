#!/usr/bin/env python3
"""Repeat a resumable fastMASST run until its completeness gate passes."""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=Path(__file__).with_name("run_fastmasst_async.py"))
    parser.add_argument("--max-rounds", type=int, default=50)
    parser.add_argument("--retry-delay", type=float, default=5.0)
    parser.add_argument("runner_args", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    runner_args = list(args.runner_args)
    if runner_args and runner_args[0] == "--":
        runner_args.pop(0)
    if not runner_args:
        raise ValueError("Runner arguments are required after --")
    if args.max_rounds < 1 or args.retry_delay < 0:
        raise ValueError("Invalid supervisor retry settings")

    for round_number in range(1, args.max_rounds + 1):
        timestamp = datetime.now(timezone.utc).isoformat()
        print(f"[{timestamp}] fastMASST supervisor round {round_number}/{args.max_rounds}", flush=True)
        completed = subprocess.run([sys.executable, str(args.runner.resolve()), *runner_args])
        if completed.returncode == 0:
            print("fastMASST completeness gate passed", flush=True)
            raise SystemExit(0)
        if completed.returncode != 2:
            print(f"fastMASST stopped with non-retryable exit code {completed.returncode}", flush=True)
            raise SystemExit(completed.returncode)
        if round_number < args.max_rounds:
            print(
                f"Completeness gate not met; retrying durable errors in {args.retry_delay:g}s",
                flush=True,
            )
            time.sleep(args.retry_delay)
    print("fastMASST retry rounds exhausted before completeness", flush=True)
    raise SystemExit(2)


if __name__ == "__main__":
    main()
