#!/bin/bash
#SBATCH --job-name=pbmc3k_prof
#SBATCH --output=profiling_results/slurm_pbmc3k_%j.out
#SBATCH --error=profiling_results/slurm_pbmc3k_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=batch

# SLURM batch script for pbmc3k profiling (coarse + instrumented + cProfile)
# Usage on cluster: sbatch slurm/run_pbmc3k.sbatch

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p profiling_results data

# Activate environment (adjust path on cluster as needed)
if [ -f "$HOME/bmi500/.venv/bin/activate" ]; then
  source "$HOME/bmi500/.venv/bin/activate"
elif [ -f ".venv/bin/activate" ]; then
  source .venv/bin/activate
elif command -v conda >/dev/null 2>&1; then
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate bmi500
fi

echo "Host: $(hostname)"
echo "Date: $(date)"
echo "Dataset: pbmc3k"
echo "Python: $(which python)"

# Ensure data is decompressed
if [ ! -d data/pbmc3k/filtered_gene_bc_matrices ]; then
  tar xvfz data/pbmc3k.tgz -C data
fi

echo "===== COARSE-GRAIN: /usr/bin/time -v ====="
/usr/bin/time -v python scanpy_pbmc.py \
  --data-dir data \
  --data-set pbmc3k \
  --out-dir data \
  --num-threads 1 \
  --profile-dir profiling_results \
  2> profiling_results/pbmc3k_time_v.txt

# Copy instrumented timings note
echo "===== DONE pbmc3k ====="
cat profiling_results/pbmc3k_section_timings.json || true
echo "Coarse-grain /usr/bin/time -v output:"
cat profiling_results/pbmc3k_time_v.txt || true
