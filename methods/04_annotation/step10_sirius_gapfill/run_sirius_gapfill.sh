#!/bin/bash
#SBATCH --job-name=sirius_gapfill
#SBATCH --partition=zen2_1024
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/lisc/data/work/ter/rahul/dreams/sirius_gapfill_2026-08-12/logs/%x_%j.log
#SBATCH --error=/lisc/data/work/ter/rahul/dreams/sirius_gapfill_2026-08-12/logs/%x_%j.err
#
# Close the SIRIUS submission gap: annotate the eligible-unsubmitted strict
# biomarkers (POS 2,209 / NEG 80). Exact 7-step pipeline from sirius_atlas_rerun_v2.sh
# (P in elements, ppm 10, BIO structure DB, spectral search, CANOPUS, MSNovelist).
# Usage:  sbatch -J sirius_gapfill_POS run_sirius_gapfill.sh POS
#         sbatch -J sirius_gapfill_NEG run_sirius_gapfill.sh NEG
set -e
MODE="${1:?usage: run_sirius_gapfill.sh POS|NEG}"
BASE=/lisc/data/work/ter/rahul/dreams/sirius_gapfill_2026-08-12
SIRIUS=/lisc/home/user/samrat/sirius/sirius/bin/sirius
INPUT="$BASE/input/sirius_gapfill_${MODE}_mzle850.mgf"
OUTPUT="$BASE/results/gapfill_${MODE}"

# libbz2 fix required by the bundled CBC ILP solver
export LD_LIBRARY_PATH="/lisc/opt/sw/software/bzip2/1.0.8-GCCcore-13.3.0/lib:${LD_LIBRARY_PATH}"

echo "=== SIRIUS gapfill ${MODE}  job ${SLURM_JOB_ID}  node $(hostname)  $(date) ==="
[ -f "$INPUT" ] || { echo "ERROR: missing $INPUT"; exit 1; }
echo "Input: $INPUT  ($(grep -c '^BEGIN IONS' "$INPUT") spectra)"
echo "Output: $OUTPUT"

echo "--- 1/7 spectral-library-search $(date) ---"
$SIRIUS -i "$INPUT" -o "$OUTPUT" --cores ${SLURM_CPUS_PER_TASK} --recompute spectral-library-search
echo "--- 2/7 formula (CHNOPS+SBrCl, ppm10, mzmax850) $(date) ---"
$SIRIUS -i "$INPUT" -o "$OUTPUT" --mzmax 850 --cores ${SLURM_CPUS_PER_TASK} --recompute \
    formula -c 50 -p orbitrap --ppm-max 10 --elements-enforced CHNOPS --elements-considered SBrCl
echo "--- 3/7 fingerprint (CSI:FingerID) $(date) ---"
$SIRIUS -o "$OUTPUT" --cores ${SLURM_CPUS_PER_TASK} --recompute fingerprint
echo "--- 4/7 structure-db-search (BIO) $(date) ---"
$SIRIUS -o "$OUTPUT" --cores ${SLURM_CPUS_PER_TASK} --recompute structure-db-search --database BIO
echo "--- 5/7 canopus $(date) ---"
$SIRIUS -o "$OUTPUT" --cores ${SLURM_CPUS_PER_TASK} --recompute canopus
echo "--- 6/7 denovo-structures (MSNovelist) $(date) ---"
$SIRIUS -o "$OUTPUT" --cores ${SLURM_CPUS_PER_TASK} --recompute denovo-structures -c 128
echo "--- 7/7 write-summaries $(date) ---"
$SIRIUS -o "$OUTPUT" --recompute write-summaries --top-hit-summary --top-k-summary 5

echo "=== done ${MODE} $(date) ==="
for f in formula_identifications canopus_structure_summary canopus_formula_summary structure_identifications denovo_structure_identifications; do
    [ -f "$OUTPUT/$f.tsv" ] && echo "$f: $(($(wc -l < "$OUTPUT/$f.tsv")-1)) rows"
done
