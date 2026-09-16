#!/usr/bin/env python3
"""Generate profiling report figures and a Markdown/PDF-ready summary."""

from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PROF = ROOT / "profiling_results"
FIGS = PROF / "figures"
FIGS.mkdir(parents=True, exist_ok=True)

DATASETS = ["pbmc3k", "pbmc6k", "pbmc10k"]
# Nominal cell counts for x-axis / extrapolation
SIZE_K = {"pbmc3k": 3, "pbmc6k": 6, "pbmc10k": 10}

# Computational sections to plot (exclude argparse/settings overhead)
PLOT_SECTIONS = [
    "io_read_10x",
    "filter_cells_genes",
    "normalize_log1p",
    "highly_variable_genes",
    "scale",
    "pca",
    "neighbors",
    "clustering",
    "umap",
    "write_h5ad",
    "rank_genes_groups",
]


def load_timings(dataset: str) -> dict:
    path = PROF / f"{dataset}_section_timings.json"
    with open(path) as f:
        return json.load(f)


def parse_time_v(dataset: str) -> dict:
    """Parse GNU time -v or macOS time -l output."""
    path = PROF / f"{dataset}_time_v.txt"
    text = path.read_text(errors="replace") if path.exists() else ""
    out = {"elapsed_sec": None, "user_sec": None, "sys_sec": None, "max_rss": None, "raw": text}

    # GNU time -v
    m = re.search(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([0-9.:]+)", text)
    if m:
        parts = m.group(1).split(":")
        if len(parts) == 3:
            out["elapsed_sec"] = int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        elif len(parts) == 2:
            out["elapsed_sec"] = int(parts[0]) * 60 + float(parts[1])
    m = re.search(r"User time \(seconds\):\s*([0-9.]+)", text)
    if m:
        out["user_sec"] = float(m.group(1))
    m = re.search(r"System time \(seconds\):\s*([0-9.]+)", text)
    if m:
        out["sys_sec"] = float(m.group(1))
    m = re.search(r"Maximum resident set size \(kbytes\):\s*([0-9]+)", text)
    if m:
        out["max_rss"] = int(m.group(1))  # kbytes

    # macOS /usr/bin/time -l fallback
    if out["elapsed_sec"] is None:
        m = re.search(r"([0-9.]+)\s+real", text)
        if m:
            out["elapsed_sec"] = float(m.group(1))
        m = re.search(r"([0-9.]+)\s+user", text)
        if m:
            out["user_sec"] = float(m.group(1))
        m = re.search(r"([0-9.]+)\s+sys", text)
        if m:
            out["sys_sec"] = float(m.group(1))
        m = re.search(r"\s+([0-9]+)\s+maximum resident set size", text)
        if m:
            # macOS reports bytes
            out["max_rss"] = int(m.group(1)) // 1024
    return out


def top_cprofile_lines(dataset: str, n: int = 15) -> str:
    path = PROF / f"{dataset}_cprofile_rank_genes.txt"
    if not path.exists():
        return "(missing)"
    lines = path.read_text(errors="replace").splitlines()
    # Keep header + top n function lines
    keep = []
    for i, line in enumerate(lines):
        keep.append(line)
        if i > n + 5:
            break
    return "\n".join(keep)


def plot_stacked_sections(all_timings: dict) -> Path:
    labels = [f"{d}\n({SIZE_K[d]}k)" for d in DATASETS]
    x = np.arange(len(DATASETS))
    bottoms = np.zeros(len(DATASETS))
    fig, ax = plt.subplots(figsize=(9, 5.5))
    cmap = plt.get_cmap("tab20")
    for i, sec in enumerate(PLOT_SECTIONS):
        vals = [float(all_timings[d].get(sec, 0.0)) for d in DATASETS]
        ax.bar(x, vals, bottom=bottoms, label=sec, color=cmap(i % 20), width=0.6)
        bottoms += np.array(vals)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Runtime (seconds)")
    ax.set_xlabel("Dataset size")
    ax.set_title("Instrumented section runtimes vs dataset size")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    out = FIGS / "instrumented_stacked_bar.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_section_lines(all_timings: dict) -> Path:
    xs = [SIZE_K[d] for d in DATASETS]
    fig, ax = plt.subplots(figsize=(9, 5.5))
    for sec in PLOT_SECTIONS:
        ys = [float(all_timings[d].get(sec, 0.0)) for d in DATASETS]
        ax.plot(xs, ys, marker="o", label=sec)
    ax.set_xticks(xs)
    ax.set_xticklabels(["pbmc3k", "pbmc6k", "pbmc10k"])
    ax.set_xlabel("Dataset (approx. cells ×1000)")
    ax.set_ylabel("Runtime (seconds)")
    ax.set_title("Section runtime growth with dataset size")
    ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8)
    fig.tight_layout()
    out = FIGS / "instrumented_lines.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def plot_cprofile_top(dataset: str) -> Path | None:
    """Parse cumulative seconds from cProfile text and bar-plot top functions."""
    path = PROF / f"{dataset}_cprofile_rank_genes.txt"
    if not path.exists():
        return None
    rows = []
    for line in path.read_text(errors="replace").splitlines():
        # e.g. "   12    0.001    0.000    1.234    0.103 filename:lineno(func)"
        m = re.match(
            r"\s*\d+\s+[0-9.]+\s+[0-9.]+\s+([0-9.]+)\s+[0-9.]+\s+(.+)$",
            line,
        )
        if not m:
            continue
        cum = float(m.group(1))
        name = m.group(2).strip()
        if "rank_genes" in name or cum >= 0.01:
            rows.append((cum, name))
    rows = sorted(rows, key=lambda t: -t[0])[:12]
    if not rows:
        return None
    fig, ax = plt.subplots(figsize=(9, 5))
    cum_vals = [r[0] for r in rows][::-1]
    names = [r[1][-60:] for r in rows][::-1]
    ax.barh(names, cum_vals, color="steelblue")
    ax.set_xlabel("Cumulative time (s)")
    ax.set_title(f"cProfile top functions — rank_genes_groups ({dataset})")
    fig.tight_layout()
    out = FIGS / f"{dataset}_cprofile_top.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out


def extrapolate_20k(all_timings: dict) -> dict:
    """Log-linear / linear fit from 3/6/10k to predict 20k per section."""
    xs = np.array([SIZE_K[d] for d in DATASETS], dtype=float)
    preds = {}
    for sec in PLOT_SECTIONS:
        ys = np.array([float(all_timings[d].get(sec, 0.0)) for d in DATASETS], dtype=float)
        # Prefer log-log power law if all positive; else linear
        if np.all(ys > 0):
            coef = np.polyfit(np.log(xs), np.log(ys), 1)
            pred = float(np.exp(coef[1] + coef[0] * np.log(20.0)))
            model = f"power (y ~ x^{coef[0]:.2f})"
        else:
            coef = np.polyfit(xs, ys, 1)
            pred = float(coef[0] * 20.0 + coef[1])
            model = "linear"
        preds[sec] = {"pred_20k_sec": max(pred, 0.0), "model": model, "observed": dict(zip(DATASETS, ys.tolist()))}
    return preds


def main():
    all_timings = {d: load_timings(d) for d in DATASETS}
    coarse = {d: parse_time_v(d) for d in DATASETS}

    stacked = plot_stacked_sections(all_timings)
    lines = plot_section_lines(all_timings)
    cprof_figs = {d: plot_cprofile_top(d) for d in DATASETS}
    preds = extrapolate_20k(all_timings)

    # Bottleneck = max mean section time across datasets
    means = {
        sec: float(np.mean([all_timings[d].get(sec, 0.0) for d in DATASETS]))
        for sec in PLOT_SECTIONS
    }
    bottleneck = max(means, key=means.get)

    md = []
    md.append("# BMI500 Profiling Exercise Report\n")
    md.append("**Author:** Abhya Reddy Ambati  \n")
    md.append("**Environment:** Local macOS run (BMI cluster SSH unavailable at time of submission).  \n")
    md.append("SLURM scripts are included under `slurm/` for cluster reproduction.\n")

    md.append("\n## a) Datasets\n")
    md.append("Profiling repeated for **pbmc3k**, **pbmc6k**, and **pbmc10k** ")
    md.append("(coarse-grain `/usr/bin/time`/`gtime -v`, instrumented `time.time()` per `#%%` section, ")
    md.append("fine-grain `cProfile` on `sc.tl.rank_genes_groups`).\n")

    md.append("\n## b) SLURM batch scripts\n")
    md.append("- `slurm/run_pbmc3k.sbatch`\n")
    md.append("- `slurm/run_pbmc6k.sbatch`\n")
    md.append("- `slurm/run_pbmc10k.sbatch`\n")
    md.append("\nSee appendix / repository for full script text.\n")

    md.append("\n## c) Coarse-grain timing (time / CPU / memory)\n")
    md.append("| Dataset | Wall (s) | User (s) | Sys (s) | Max RSS (MB) |\n")
    md.append("|---|---:|---:|---:|---:|\n")
    for d in DATASETS:
        c = coarse[d]
        rss_mb = (c["max_rss"] / 1024.0) if c["max_rss"] else float("nan")
        md.append(
            f"| **{d}** | {c['elapsed_sec']} | {c['user_sec']} | {c['sys_sec']} | {rss_mb:.1f} |\n"
        )
    md.append("\nFull console / `time -v` logs: `profiling_results/*_console.log`, `*_time_v.txt`.\n")

    md.append("\n### Instrumented major bottleneck\n")
    md.append(f"Largest average section time across datasets: **`{bottleneck}`** ")
    md.append(f"(mean {means[bottleneck]:.2f} s).\n")

    md.append("\n### cProfile visualization (`rank_genes_groups`)\n")
    for d in DATASETS:
        fig = cprof_figs[d]
        if fig:
            md.append(f"![{d} cProfile]({fig.relative_to(ROOT).as_posix()})\n\n")

    md.append("\n## d) Instrumented profile vs dataset size\n")
    md.append(f"![stacked]({stacked.relative_to(ROOT).as_posix()})\n\n")
    md.append(f"![lines]({lines.relative_to(ROOT).as_posix()})\n\n")

    md.append("| Dataset | " + " | ".join(PLOT_SECTIONS) + " |\n")
    md.append("|---|" + "|".join(["---:"] * len(PLOT_SECTIONS)) + "|\n")
    for d in DATASETS:
        vals = [f"{all_timings[d].get(s, 0):.3f}" for s in PLOT_SECTIONS]
        md.append(f"| {d} | " + " | ".join(vals) + " |\n")

    md.append("\n## e) Bottlenecks and 20k extrapolation\n")
    md.append(f"- Primary bottleneck section: **{bottleneck}**.\n")
    md.append("- Other costly steps typically include neighbors/UMAP/PCA and `rank_genes_groups` ")
    md.append("(Wilcoxon over many genes × clusters), which grow with cells and clusters.\n")
    md.append("- Extrapolation to **20k cells** uses a log–log power-law fit when possible ")
    md.append("(else linear) on the 3k/6k/10k observations:\n\n")
    md.append("| Section | 3k | 6k | 10k | Pred. 20k (s) | Model |\n")
    md.append("|---|---:|---:|---:|---:|---|\n")
    for sec in PLOT_SECTIONS:
        p = preds[sec]
        o = p["observed"]
        md.append(
            f"| {sec} | {o['pbmc3k']:.3f} | {o['pbmc6k']:.3f} | {o['pbmc10k']:.3f} | "
            f"{p['pred_20k_sec']:.3f} | {p['model']} |\n"
        )
    total_20 = sum(preds[s]["pred_20k_sec"] for s in PLOT_SECTIONS)
    md.append(f"\n**Sum of predicted section times at 20k ≈ {total_20:.1f} s** ")
    md.append("(rough; ignores I/O cache effects and parallelization).\n")

    md.append("\n## Notes\n")
    md.append("- If `louvain` failed to build locally, clustering used `leiden` (same pipeline stage).\n")
    md.append("- Re-run on BMI cluster with `sbatch slurm/run_pbmc*.sbatch` when SSH access is available.\n")

    report_md = PROF / "profiling_report.md"
    report_md.write_text("".join(md))
    print(f"Wrote {report_md}")

    # Also dump predictions JSON
    (PROF / "extrapolation_20k.json").write_text(json.dumps(preds, indent=2))
    print(f"Figures in {FIGS}")


if __name__ == "__main__":
    main()
