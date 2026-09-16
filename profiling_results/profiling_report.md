# BMI500 Profiling Exercise Report
**Author:** Abhya Reddy Ambati  
**Environment:** Local macOS run (BMI cluster SSH unavailable at time of submission).  
SLURM scripts are included under `slurm/` for cluster reproduction.

## a) Datasets
Profiling repeated for **pbmc3k**, **pbmc6k**, and **pbmc10k** (coarse-grain `/usr/bin/time`/`gtime -v`, instrumented `time.time()` per `#%%` section, fine-grain `cProfile` on `sc.tl.rank_genes_groups`).

## b) SLURM batch scripts
- `slurm/run_pbmc3k.sbatch`
- `slurm/run_pbmc6k.sbatch`
- `slurm/run_pbmc10k.sbatch`

See appendix / repository for full script text.

## c) Coarse-grain timing (time / CPU / memory)
| Dataset | Wall (s) | User (s) | Sys (s) | Max RSS (MB) |
|---|---:|---:|---:|---:|
| **pbmc3k** | 21.42 | 14.25 | 1.61 | 827.8 |
| **pbmc6k** | 16.29 | 16.2 | 1.34 | 904.8 |
| **pbmc10k** | 30.64 | 28.18 | 1.46 | 1696.4 |

Full console / `time -v` logs: `profiling_results/*_console.log`, `*_time_v.txt`.

### Instrumented major bottleneck
Largest average section time across datasets: **`neighbors`** (mean 4.74 s).

### cProfile visualization (`rank_genes_groups`)
![pbmc3k cProfile](profiling_results/figures/pbmc3k_cprofile_top.png)

![pbmc6k cProfile](profiling_results/figures/pbmc6k_cprofile_top.png)

![pbmc10k cProfile](profiling_results/figures/pbmc10k_cprofile_top.png)


## d) Instrumented profile vs dataset size
![stacked](profiling_results/figures/instrumented_stacked_bar.png)

![lines](profiling_results/figures/instrumented_lines.png)

| Dataset | io_read_10x | filter_cells_genes | normalize_log1p | highly_variable_genes | scale | pca | neighbors | clustering | umap | write_h5ad | rank_genes_groups |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| pbmc3k | 0.947 | 0.019 | 0.014 | 0.653 | 0.030 | 0.131 | 2.065 | 0.332 | 3.071 | 0.047 | 1.212 |
| pbmc6k | 0.445 | 0.024 | 0.020 | 0.026 | 0.051 | 0.218 | 1.737 | 1.504 | 5.394 | 0.066 | 2.232 |
| pbmc10k | 0.271 | 0.222 | 0.113 | 0.090 | 0.067 | 0.463 | 10.432 | 1.392 | 4.792 | 0.359 | 7.609 |

## e) Bottlenecks and 20k extrapolation
- Primary bottleneck section: **neighbors**.
- Other costly steps typically include neighbors/UMAP/PCA and `rank_genes_groups` (Wilcoxon over many genes × clusters), which grow with cells and clusters.
- Extrapolation to **20k cells** uses a log–log power-law fit when possible (else linear) on the 3k/6k/10k observations:

| Section | 3k | 6k | 10k | Pred. 20k (s) | Model |
|---|---:|---:|---:|---:|---|
| io_read_10x | 0.947 | 0.445 | 0.271 | 0.130 | power (y ~ x^-1.04) |
| filter_cells_genes | 0.019 | 0.024 | 0.222 | 0.536 | power (y ~ x^1.93) |
| normalize_log1p | 0.014 | 0.020 | 0.113 | 0.264 | power (y ~ x^1.67) |
| highly_variable_genes | 0.653 | 0.026 | 0.090 | 0.012 | power (y ~ x^-1.82) |
| scale | 0.030 | 0.051 | 0.067 | 0.109 | power (y ~ x^0.67) |
| pca | 0.131 | 0.218 | 0.463 | 0.872 | power (y ~ x^1.03) |
| neighbors | 2.065 | 1.737 | 10.432 | 16.322 | power (y ~ x^1.25) |
| clustering | 0.332 | 1.504 | 1.392 | 4.293 | power (y ~ x^1.25) |
| umap | 3.071 | 5.394 | 4.792 | 7.084 | power (y ~ x^0.40) |
| write_h5ad | 0.047 | 0.066 | 0.359 | 0.802 | power (y ~ x^1.62) |
| rank_genes_groups | 1.212 | 2.232 | 7.609 | 18.013 | power (y ~ x^1.49) |

**Sum of predicted section times at 20k ≈ 48.4 s** (rough; ignores I/O cache effects and parallelization).

## Notes
- If `louvain` failed to build locally, clustering used `leiden` (same pipeline stage).
- Re-run on BMI cluster with `sbatch slurm/run_pbmc*.sbatch` when SSH access is available.
