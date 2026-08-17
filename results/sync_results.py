"""Sync built figures/tables from P2R into the repo's results/ directory.

Usage:
    python sync_results.py          # dry-run: preview what would copy
    python sync_results.py --apply  # actually copy files

Maps the canonical rendered assets (the submission-package picks) from the
P2R working tree into results/. Add new assets to FIGURE_MAP.
"""

import argparse
import hashlib
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
P2R_ROOT = Path(r"C:\Users\Shadow\Desktop\P2R")
ANALYSIS_DIR = P2R_ROOT / "outputs" / "analysis" / "ncbi-phylum-2026-08-04-v1"
SUBMISSION_DIR = P2R_ROOT / "submission_package_strict16_2026-08-13"

FIGURE_MAP = {
    # figure 1 — data-driven concept skeleton (rendered straight from the release)
    "figure1_svg": {
        "source": ANALYSIS_DIR / "figure1_redesign_2026-08-11_v1" / "Figure1_concept_skeleton.svg",
        "target": REPO_ROOT / "figures" / "figure1" / "Figure1_concept_skeleton.svg",
    },
}

# main figures 2-5: submission-package picks
for _n in range(2, 6):
    FIGURE_MAP[f"figure{_n}"] = {
        "source": SUBMISSION_DIR / "main_figures" / f"Figure_{_n}.png",
        "target": REPO_ROOT / "figures" / f"figure{_n}" / f"Figure_{_n}.png",
    }

# supplementary figures 1-8
for _n in range(1, 9):
    FIGURE_MAP[f"suppfig{_n}"] = {
        "source": SUBMISSION_DIR / "supplementary_figures" / f"Supplementary_Figure_{_n}.png",
        "target": REPO_ROOT / "figures" / "supplementary" / f"Supplementary_Figure_{_n}.png",
    }

# supplementary tables S1-S15
for _n in range(1, 16):
    FIGURE_MAP[f"table_s{_n}"] = {
        "source": SUBMISSION_DIR / "supplementary_tables" / f"Table_S{_n}.xlsx",
        "target": REPO_ROOT / "tables" / f"Table_S{_n}.xlsx",
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--apply", action="store_true",
                        help="Actually copy files (default: dry-run)")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show per-file details")
    args = parser.parse_args()

    missing = []
    changed = []
    unchanged = []
    copied = 0

    for name, mapping in FIGURE_MAP.items():
        src = mapping["source"]
        dst = mapping["target"]

        if not src.exists():
            missing.append(name)
            if args.verbose:
                print(f"  [MISSING] {name}: source not found\n"
                      f"      {src}")
            continue

        dst.parent.mkdir(parents=True, exist_ok=True)

        if dst.exists() and sha256(src) == sha256(dst):
            unchanged.append(name)
            if args.verbose:
                print(f"  [OK]      {name}: unchanged")
            continue

        changed.append(name)
        src_size = src.stat().st_size

        if args.apply:
            shutil.copy2(src, dst)
            copied += 1
            status = "COPIED"
        else:
            status = "WOULD COPY"

        print(f"  [{status}]  {name} ({src_size:,} B)")
        if args.verbose:
            print(f"      src: {src}")
            print(f"      dst: {dst}")

    print()
    print(f"Summary: {len(changed)} to sync"
          f" ({copied} {'copied' if args.apply else 'would copy'}),"
          f" {len(unchanged)} already current,"
          f" {len(missing)} source files missing")

    if missing and args.verbose:
        print("\nMissing sources:")
        for name in missing:
            print(f"  - {name}: {FIGURE_MAP[name]['source']}")

    return 0 if not missing else 1


if __name__ == "__main__":
    raise SystemExit(main())