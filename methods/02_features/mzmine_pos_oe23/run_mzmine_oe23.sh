#!/bin/bash
#SBATCH --job-name=mzmine_OE23POS
#SBATCH --output=/lisc/home/user/samrat/snakeautomatic/results/OE23-POS/mzmine_%j.log
#SBATCH --error=/lisc/home/user/samrat/snakeautomatic/results/OE23-POS/mzmine_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00

echo "=========================================="
echo "MZmine OE23-POS Batch Processing"
echo "=========================================="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Start time: $(date)"
echo "=========================================="

# Load MZmine
echo "[$(date '+%H:%M:%S')] Loading MZmine 4.9.0 module..."
module load mzmine/4.9.0
echo "[$(date '+%H:%M:%S')] MZmine loaded: $(which mzmine)"

# Fix ThermoRawFileParser permission issue - use FORK instead of posix_spawn
export JAVA_TOOL_OPTIONS="-Djdk.lang.Process.launchMechanism=FORK"
echo "[$(date '+%H:%M:%S')] Set JAVA_TOOL_OPTIONS: $JAVA_TOOL_OPTIONS"

# Verify external_tools are accessible
echo "[$(date '+%H:%M:%S')] Checking external_tools..."
if [ -x /lisc/home/user/samrat/external_tools/thermo_raw_file_parser/ThermoRawFileParser ]; then
    echo "[$(date '+%H:%M:%S')] ✓ ThermoRawFileParser is executable"
else
    echo "[$(date '+%H:%M:%S')] ✗ ThermoRawFileParser NOT executable - attempting fix"
    chmod +x /lisc/home/user/samrat/external_tools/thermo_raw_file_parser/ThermoRawFileParser
fi

# Verify raw files exist
RAW_COUNT=$(ls /lisc/home/user/samrat/snakeautomatic/data/raw/OE23-POS/*.raw 2>/dev/null | wc -l)
echo "[$(date '+%H:%M:%S')] Found $RAW_COUNT raw files in OE23-POS"

# Verify output directory
mkdir -p /lisc/home/user/samrat/snakeautomatic/results/OE23-POS
echo "[$(date '+%H:%M:%S')] Output dir: /lisc/home/user/samrat/snakeautomatic/results/OE23-POS"

# Verify batch file
BATCH="/lisc/home/user/samrat/snakeautomatic/soilmass-pos-oe23.mzbatch"
if [ -f "$BATCH" ]; then
    echo "[$(date '+%H:%M:%S')] ✓ Batch file found: $BATCH"
else
    echo "[$(date '+%H:%M:%S')] ✗ Batch file NOT found!"
    exit 1
fi

echo "=========================================="
echo "[$(date '+%H:%M:%S')] Starting MZmine batch processing..."
echo "=========================================="

# Change to home dir where external_tools exists
cd /lisc/home/user/samrat

# Run MZmine in batch (headless) mode with verbose output
mzmine -batch "$BATCH" 2>&1 | while IFS= read -r line; do
    echo "[$(date '+%H:%M:%S')] $line"
done

EXIT_CODE=${PIPESTATUS[0]}

echo "=========================================="
echo "[$(date '+%H:%M:%S')] MZmine batch finished with exit code: $EXIT_CODE"
echo "=========================================="

# Check what output files were created
echo "[$(date '+%H:%M:%S')] Output files:"
ls -lh /lisc/home/user/samrat/snakeautomatic/results/OE23-POS/ 2>/dev/null

echo "=========================================="
echo "End time: $(date)"
echo "=========================================="
