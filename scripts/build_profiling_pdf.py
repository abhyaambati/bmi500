#!/usr/bin/env python3
"""Build profiling_report.pdf from profiling_results/ (matplotlib PdfPages)."""

from pathlib import Path
import json
import re

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROOT = Path(__file__).resolve().parents[1]
PROF = ROOT / "profiling_results"
FIGS = PROF / "figures"
OUT = PROF / "profiling_report.pdf"


def parse_time(ds: str) -> dict:
    text = (PROF / f"{ds}_time_v.txt").read_text(errors="replace")

    def grab(pat):
        m = re.search(pat, text)
        return m.group(1) if m else "n/a"

    return {
        "wall": grab(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([0-9.:]+)"),
        "user": grab(r"User time \(seconds\):\s*([0-9.]+)"),
        "sys": grab(r"System time \(seconds\):\s*([0-9.]+)"),
        "rss": grab(r"Maximum resident set size \(kbytes\):\s*([0-9]+)"),
    }


def text_page(pdf, title, body):
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.08, 0.95, title, fontsize=14, fontweight="bold", va="top")
    fig.text(0.08, 0.90, body, fontsize=9, va="top", family="monospace")
    pdf.savefig(fig)
    plt.close(fig)


def image_page(pdf, title, path):
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.08, 0.96, title, fontsize=12, fontweight="bold", va="top")
    ax = fig.add_axes([0.08, 0.15, 0.84, 0.75])
    ax.imshow(plt.imread(path))
    ax.axis("off")
    pdf.savefig(fig)
    plt.close(fig)


def main():
    timings = {
        d: json.loads((PROF / f"{d}_section_timings.json").read_text())
        for d in ["pbmc3k", "pbmc6k", "pbmc10k"]
    }
    preds = json.loads((PROF / "extrapolation_20k.json").read_text())
    sections = [
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
    means = {
        s: sum(float(timings[d].get(s, 0)) for d in timings) / 3 for s in sections
    }
    bottleneck = max(means, key=means.get)

    with PdfPages(OUT) as pdf:
        lines = [
            "BMI500 Profiling Exercise Report",
            "Author: Abhya Reddy Ambati",
            "Environment: Local macOS (cluster SSH unavailable).",
            "SLURM scripts: slurm/run_pbmc{3,6,10}k.sbatch",
            "",
            "=== c) Coarse-grain (gtime -v) ===",
        ]
        for ds in ["pbmc3k", "pbmc6k", "pbmc10k"]:
            t = parse_time(ds)
            rss = f"{int(t['rss'])/1024:.1f} MB" if t["rss"] != "n/a" else "n/a"
            lines.append(
                f"{ds}: wall={t['wall']} user={t['user']}s sys={t['sys']}s maxRSS={rss}"
            )
        lines += [
            "",
            "=== Instrumented section timings (s) ===",
            f"{'section':28s}{'3k':>8}{'6k':>8}{'10k':>8}",
        ]
        for s in sections:
            lines.append(
                f"{s:28s}{timings['pbmc3k'].get(s, 0):8.2f}"
                f"{timings['pbmc6k'].get(s, 0):8.2f}{timings['pbmc10k'].get(s, 0):8.2f}"
            )
        lines += [
            "",
            f"Major bottleneck (avg): {bottleneck} ({means[bottleneck]:.2f}s)",
            "cProfile hot path: pandas.rank / Wilcoxon in rank_genes_groups",
        ]
        text_page(pdf, "Summary", "\n".join(lines))

        image_page(
            pdf,
            "d) Stacked section runtimes vs dataset size",
            FIGS / "instrumented_stacked_bar.png",
        )
        image_page(
            pdf,
            "Section runtime lines vs dataset size",
            FIGS / "instrumented_lines.png",
        )
        for ds in ["pbmc3k", "pbmc6k", "pbmc10k"]:
            image_page(
                pdf,
                f"cProfile top - rank_genes_groups ({ds})",
                FIGS / f"{ds}_cprofile_top.png",
            )

        lines = [
            "e) Bottlenecks and 20k extrapolation",
            "",
            f"Primary bottleneck: {bottleneck}.",
            "neighbors and rank_genes_groups grow with cells/clusters;",
            "UMAP is also expensive.",
            "",
            f"{'section':28s}{'pred20k':>10}  model",
        ]
        total = 0.0
        for s in sections:
            p = float(preds[s]["pred_20k_sec"])
            total += p
            lines.append(f"{s:28s}{p:10.2f}  {preds[s]['model']}")
        lines += [
            "",
            f"Sum predicted section times at 20k ~ {total:.1f} s",
            "",
            "b) SLURM: run_pbmc3k/6k/10k.sbatch in slurm/",
            "Note: leiden used locally (louvain build failed on macOS).",
        ]
        text_page(pdf, "Analysis", "\n".join(lines))

    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
