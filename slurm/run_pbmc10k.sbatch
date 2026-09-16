#!/bin/bash
#SBATCH --job-name=pbmc10k_prof
#SBATCH --output=profiling_results/slurm_pbmc10k_%j.out
#SBATCH --error=profiling_results/slurm_pbmc10k_%j.err
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --partition=batch

# SLURM batch script for pbmc10k profiling (coarse + instrumented + cProfile)
# Usage on cluster: sbatch slurm/run_pbmc10k.sbatch

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p profiling_results data

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
echo "Dataset: pbmc10k"
echo "Python: $(which python)"

if [ ! -d data/pbmc10k/filtered_gene_bc_matrices ]; then
  tar xvfz data/pbmc10k.tgz -C data
fi

echo "===== COARSE-GRAIN: /usr/bin/time -v ====="
/usr/bin/time -v python scanpy_pbmc.py \
  --data-dir data \
  --data-set pbmc10k \
  --out-dir data \
  --num-threads 1 \
  --profile-dir profiling_results \
  2> profiling_results/pbmc10k_time_v.txt

echo "===== DONE pbmc10k ====="
cat profiling_results/pbmc10k_section_timings.json || true
echo "Coarse-grain /usr/bin/time -v output:"
cat profiling_results/pbmc10k_time_v.txt || true
