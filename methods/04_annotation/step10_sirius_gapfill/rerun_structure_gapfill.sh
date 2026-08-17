#!/bin/bash
#SBATCH --job-name=sirius_struct_rerun
#SBATCH --partition=zen2_1024
#SBATCH --qos=normal
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/lisc/data/work/ter/rahul/dreams/sirius_gapfill_2026-08-12/logs/%x_%j.log
#SBATCH --error=/lisc/data/work/ter/rahul/dreams/sirius_gapfill_2026-08-12/logs/%x_%j.err
#
# CSI:FingerID structure recovery v2. First attempts returned 0 for two reasons:
#   (1) -o gapfill_POS resolved to the empty TSV *directory* (both gapfill_POS/ and
#       gapfill_POS.sirius exist) -> 0 compounds. FIX: pass the .sirius project file.
#   (2) --database BIO is too narrow for these microbial/archaeal lipids.
# Interactive test with `-o gapfill_POS.sirius --database PUBCHEM` genuinely ran
# (48 instances in 2m15s, ~21/min). This job recomputes ONLY structure-db-search
# over a broad DB union on the real projects (fingerprints cached), then rewrites
# summaries. POS then NEG SEQUENTIALLY (one SIRIUS process, no token race).
BASE=/lisc/data/work/ter/rahul/dreams/sirius_gapfill_2026-08-12
SIRIUS=/lisc/home/user/samrat/sirius/sirius/bin/sirius
DB="BIO,PUBCHEM,HMDB,GNPS,YMDB,PLANTCYC,KNAPSACK"
export LD_LIBRARY_PATH="/lisc/opt/sw/software/bzip2/1.0.8-GCCcore-13.3.0/lib:${LD_LIBRARY_PATH}"

echo "=== structure rerun v2 (FIX: .sirius path + PubChem)  job ${SLURM_JOB_ID}  node $(hostname)  $(date) ==="
echo "DB union: $DB"
$SIRIUS login --show 2>&1 | grep -iE "Logged in|expires" || true

for MODE in POS NEG; do
  PROJ="$BASE/results/gapfill_${MODE}.sirius"
  SUM="$BASE/results/gapfill_${MODE}"
  if [ ! -f "$PROJ" ]; then echo "ERROR: missing project $PROJ"; continue; fi
  echo "--- [$MODE] structure-db-search  $(date) ---"
  $SIRIUS -o "$PROJ" --cores ${SLURM_CPUS_PER_TASK} --recompute structure-db-search --database "$DB" || echo "[$MODE] structure-db-search returned non-zero"
  echo "--- [$MODE] write-summaries  $(date) ---"
  $SIRIUS -o "$PROJ" --recompute write-summaries --top-hit-summary --top-k-summary 5 || echo "[$MODE] write-summaries returned non-zero"
  echo "--- [$MODE] result rows (from $SUM) ---"
  for f in structure_identifications canopus_structure_summary formula_identifications canopus_formula_summary denovo_structure_identifications; do
    [ -f "$SUM/$f.tsv" ] && echo "$MODE $f: $(($(wc -l < "$SUM/$f.tsv")-1)) rows"
  done
  echo "--- [$MODE] actual structure TSV location ---"
  find "$BASE/results" -name "structure_identifications.tsv" -path "*gapfill_${MODE}*" -printf "%p  %s bytes\n" 2>/dev/null
done
echo "=== done  $(date) ==="
