#!/bin/bash
# Local profiling runner (Mac/Linux) — mirrors SLURM scripts without sbatch.
# Usage: ./scripts/run_local_profiling.sh pbmc3k

set -euo pipefail
DATASET="${1:?usage: $0 pbmc3k|pbmc6k|pbmc10k}"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
mkdir -p profiling_results data

source .venv/bin/activate

if [ ! -d "data/${DATASET}/filtered_gene_bc_matrices" ]; then
  tar xvfz "data/${DATASET}.tgz" -C data
fi

# Prefer GNU time (-v); on macOS use brew's gtime if present, else /usr/bin/time -l
if command -v gtime >/dev/null 2>&1; then
  TIME_CMD=(gtime -v)
elif /usr/bin/time -v true >/dev/null 2>&1; then
  TIME_CMD=(/usr/bin/time -v)
else
  TIME_CMD=(/usr/bin/time -l)
fi

OUT_TIME="profiling_results/${DATASET}_time_v.txt"
OUT_LOG="profiling_results/${DATASET}_console.log"

echo "===== Profiling ${DATASET} on $(hostname) at $(date) =====" | tee "$OUT_LOG"
echo "TIME_CMD: ${TIME_CMD[*]}" | tee -a "$OUT_LOG"

# Coarse-grain wraps the full instrumented+cProfile run
{
  "${TIME_CMD[@]}" python scanpy_pbmc.py \
    --data-dir data \
    --data-set "$DATASET" \
    --out-dir data \
    --num-threads 1 \
    --profile-dir profiling_results
} > >(tee -a "$OUT_LOG") 2> >(tee "$OUT_TIME" | tee -a "$OUT_LOG" >&2)

echo "===== DONE ${DATASET} =====" | tee -a "$OUT_LOG"
