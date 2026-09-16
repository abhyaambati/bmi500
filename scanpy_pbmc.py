# %%
"""PBMC Scanpy pipeline with instrumented and fine-grain profiling."""

from typing_extensions import ParamSpecArgs  # noqa: F401  # kept from upstream stub
import numpy as np  # noqa: F401
import pandas as pd  # noqa: F401
import scanpy as sc

import time
import sys
import argparse
import cProfile
import pstats
import io
import json
import os

# Section timings collected for instrumented profiling report
SECTION_TIMINGS = {}


def time_section(name):
    """Context-manager style helper: record wall time for a named section."""

    class _Timer:
        def __enter__(self):
            self.t0 = time.time()
            print(f"[TIMING] START {name}")
            return self

        def __exit__(self, exc_type, exc, tb):
            elapsed = time.time() - self.t0
            SECTION_TIMINGS[name] = elapsed
            print(f"[TIMING] END   {name}: {elapsed:.4f} s")

    return _Timer()


# %%
with time_section("settings"):
    sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
    sc.logging.print_header()
    sc.settings.set_figure_params(dpi=80, facecolor="white")
    sc.settings.n_jobs = 1
    print(f"using {sc.settings.n_jobs} threads")

# %%
with time_section("argparse"):
    parser = argparse.ArgumentParser(description="Process arguments.")
    parser.add_argument(
        "--data-dir",
        type=str,
        help="Directory containing the dataset subdirectories",
        default="data",
    )
    parser.add_argument(
        "--data-set",
        type=str,
        help="Dataset name, which is the subdirectory name",
        default="pbmc3k",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        help="Output directory",
        required=False,
        default="data",
    )
    parser.add_argument(
        "--num-threads",
        type=int,
        help="Number of threads",
        default=1,
        required=False,
    )
    parser.add_argument(
        "--profile-dir",
        type=str,
        help="Directory for timing JSON and cProfile outputs",
        default="profiling_results",
    )

    args = parser.parse_args()

    datadir = args.data_dir if args.data_dir.endswith("/") else args.data_dir + "/"
    dataset = args.data_set
    outdir = args.out_dir if args.out_dir.endswith("/") else args.out_dir + "/"
    nthreads = args.num_threads
    profile_dir = args.profile_dir
    os.makedirs(profile_dir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)

# %%
# I/O
with time_section("io_read_10x"):
    results_file = "/".join([outdir, dataset + ".scanpy.h5ad"])

    adata = sc.read_10x_mtx(
        "/".join([datadir, dataset, "filtered_gene_bc_matrices"]),
        var_names="gene_symbols",
        cache=True,
    )

    adata.var_names_make_unique()

# %%
# preprocessing — basic filtering
with time_section("filter_cells_genes"):
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

# %%
# normalize + log1p
with time_section("normalize_log1p"):
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

# %%
# highly variable genes
with time_section("highly_variable_genes"):
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]

# %%
# scale
with time_section("scale"):
    sc.pp.scale(adata)

# %%
# pca
with time_section("pca"):
    sc.tl.pca(adata, svd_solver="arpack", n_comps=30)

# %%
# neighborhood graph
with time_section("neighbors"):
    sc.pp.neighbors(adata, n_pcs=30)

# %%
# clustering
# Prefer louvain (assignment default). Fall back to leiden if louvain is not installed
# (common on macOS when the louvain C core fails to build).
with time_section("clustering"):
    try:
        import louvain  # noqa: F401

        sc.tl.louvain(adata, resolution=0.5)
        cluster_key = "louvain"
    except ImportError:
        print("[WARN] louvain not available; using leiden clustering instead")
        sc.tl.leiden(adata, resolution=0.5)
        cluster_key = "leiden"
    print(f"[INFO] cluster key: {cluster_key}")

# %%
# umap
with time_section("umap"):
    sc.tl.umap(adata, n_components=30)

# %%
with time_section("write_h5ad"):
    adata.write(results_file)
    print(adata)

# %%
# find marker genes — fine-grain profile with cProfile
with time_section("rank_genes_groups"):
    pr = cProfile.Profile()
    pr.enable()
    sc.tl.rank_genes_groups(adata, cluster_key, method="wilcoxon", use_raw=True)
    pr.disable()

    # Human-readable stats
    s = io.StringIO()
    ps = pstats.Stats(pr, stream=s).sort_stats("cumulative")
    ps.print_stats(40)
    stats_text = s.getvalue()
    print(stats_text)

    stats_path = os.path.join(profile_dir, f"{dataset}_cprofile_rank_genes.txt")
    with open(stats_path, "w") as f:
        f.write(stats_text)
    print(f"[PROFILE] wrote {stats_path}")

    # Binary profile for snakeviz / pstats visualization
    pstats_path = os.path.join(profile_dir, f"{dataset}_cprofile_rank_genes.pstats")
    pr.dump_stats(pstats_path)
    print(f"[PROFILE] wrote {pstats_path}")

# Persist instrumented section timings
SECTION_TIMINGS["dataset"] = dataset
SECTION_TIMINGS["n_obs"] = int(adata.n_obs)
SECTION_TIMINGS["n_vars"] = int(adata.n_vars)
timing_path = os.path.join(profile_dir, f"{dataset}_section_timings.json")
with open(timing_path, "w") as f:
    json.dump(SECTION_TIMINGS, f, indent=2)
print(f"[TIMING] wrote {timing_path}")
print("[TIMING] SUMMARY")
for k, v in SECTION_TIMINGS.items():
    if isinstance(v, float):
        print(f"  {k:24s} {v:10.4f} s")
    else:
        print(f"  {k:24s} {v}")
