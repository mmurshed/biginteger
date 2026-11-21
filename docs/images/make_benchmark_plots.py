#!/usr/bin/env python3
"""Plot the single-thread and default-build ratios from the checked-in benchmark CSVs."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

HERE = Path(__file__).resolve().parent
BENCHMARKS = HERE.parent / "benchmarks"

# Session progress: BigMath/GMP ratio at the 2026-05-30 baseline vs after the
# 2026-06 runs (PRs #82-#99, #107-#114, then #116-#119).
BEFORE_AFTER = [
    ("mul 1M×1M", 1.15, 0.85),
    ("mul 10M×10M", 1.05, 0.31),
    ("mul 100M×100M", 1.80, 0.91),
    ("div 1M×200k", 3.36, 1.44),
    ("div 10M×2M", 2.79, 1.18),
    ("div 100M×20M", 2.85, 0.49),
    ("tostr 1M", 4.48, 1.29),
    ("parse 20M", 2.33, 0.38),
]


def ratio_plot(path, results, figure_title, threads, legend_location="upper right"):
    with results.open(newline="") as src:
        rows = list(csv.DictReader(line for line in src if not line.startswith("#")))

    def series(op, phase, shape=None):
        points = []
        for row in rows:
            if row["op"] != op:
                continue
            digits = [int(n) for n in row["label"].split()[0].split("x")]
            if shape == "equal" and digits[0] != digits[1]:
                continue
            if shape == "skew" and digits[0] == digits[1]:
                continue
            ratio = float(row[f"bigmath_{phase}_ms"]) / float(row[f"gmp_{phase}_ms"])
            points.append((digits[0], ratio))
        return zip(*sorted(points))

    panels = [
        ("Multiplication", [("Equal lengths", "mul", "warm", "equal"),
                            ("10:1 lengths", "mul", "warm", "skew")]),
        ("Skewed division", [("Quotient", "div_q", "warm", "skew"),
                             ("Quotient + remainder", "div_qr", "warm", "skew")]),
        ("Decimal parsing", [("First call", "parse", "first", None),
                             ("Warm", "parse", "warm", None)]),
        ("Decimal formatting", [("First call", "tostr", "first", None),
                                ("Warm", "tostr", "warm", None)]),
    ]

    fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.6), dpi=180, sharex=True, sharey=True)
    for ax, (title, lines) in zip(axes.flat, panels):
        ax.set_title(title, loc="left", fontsize=13, fontweight="bold", pad=6)
        ax.set_xscale("log")
        ax.set_yscale("log", base=2)
        ax.axhline(1, color="#424b57", linewidth=1, linestyle=(0, (4, 3)))
        ax.grid(axis="y", color="#dce1e6", linewidth=0.6)
        for (label, op, phase, shape), color, marker, style in zip(
                lines, ["#2165a1", "#a54718"], ["o", "s"], ["-", "--"]):
            ax.plot(*series(op, phase, shape), label=label, color=color, marker=marker,
                    linestyle=style, markersize=4.4, linewidth=1.7,
                    markerfacecolor="white" if style == "--" else color)
        ax.set_xlim(700, 3e8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["bottom", "left"]].set_color("#a9b1b9")
        ax.tick_params(labelsize=9)
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.legend(loc=legend_location, frameon=False, fontsize=9)
    ticks = [1 / 16, 1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8]
    lo, hi = axes[0, 0].get_ylim()
    if lo > 0.5:
        lo = 0.5
        axes[0, 0].set_ylim(bottom=lo)
    visible = [t for t in ticks if lo <= t <= hi]
    axes[0, 0].set_yticks(visible, [f"{t:g}×" for t in visible])
    axes[0, 0].set_xticks([1e3, 1e4, 1e5, 1e6, 1e7, 1e8],
                          ["1k", "10k", "100k", "1M", "10M", "100M"])
    fig.suptitle(figure_title, fontsize=17, fontweight="bold", y=0.98)
    fig.text(0.5, 0.918, "Time relative to GMP  ·  Below 1× means BigMath is faster", ha="center", fontsize=11)
    fig.supxlabel("Decimal digits in the first operand", fontsize=11, y=0.066)
    fig.text(0.5, 0.021, f"November 21, 2026  ·  {threads}  ·  "
             "Arithmetic panels show warm calls", ha="center", fontsize=9, color="#384655")
    fig.subplots_adjust(left=0.075, right=0.978, top=0.845, bottom=0.15, hspace=0.3, wspace=0.13)
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)


def small_ops_plot(path, results):
    """BigMath/GMP ratio on 1-16-limb operands, before and after inline storage."""
    with results.open(newline="") as src:
        rows = list(csv.DictReader(line for line in src if not line.startswith("#")))

    def series(version, op):
        points = sorted((int(r["limbs"]), float(r["ratio"])) for r in rows
                        if r["version"] == version and r["op"] == op)
        return zip(*points)

    panels = [("Multiplication", "mul"), ("Addition", "add"), ("Subtraction", "sub"),
              ("Division (2n ÷ n limbs)", "div"), ("Remainder (2n ÷ n limbs)", "mod"), ("Copy", "copy")]
    fig, axes = plt.subplots(2, 3, figsize=(11.4, 6.2), dpi=180, sharex=True, sharey=True)
    for ax, (title, op) in zip(axes.flat, panels):
        ax.set_title(title, loc="left", fontsize=12, fontweight="bold", pad=6)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
        ax.axhline(1, color="#424b57", linewidth=1, linestyle=(0, (4, 3)))
        ax.grid(axis="y", color="#dce1e6", linewidth=0.6)
        for version, label, color, marker, style in [
                ("before", "Before (heap-only storage)", "#a54718", "s", "--"),
                ("after", "After (inline storage)", "#2165a1", "o", "-")]:
            ax.plot(*series(version, op), label=label, color=color, marker=marker, linestyle=style,
                    markersize=4.4, linewidth=1.7, markerfacecolor="white" if style == "--" else color)
        ax.spines[["top", "right"]].set_visible(False)
        ax.spines[["bottom", "left"]].set_color("#a9b1b9")
        ax.tick_params(labelsize=9)
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.xaxis.set_minor_formatter(NullFormatter())
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.895), ncol=2,
               frameon=False, fontsize=10)
    ratios = [float(r["ratio"]) for r in rows]
    lo, hi = min(ratios) / 1.3, max(ratios) * 1.3
    axes[0, 0].set_ylim(lo, hi)
    ticks = [t for t in [1 / 16, 1 / 8, 1 / 4, 1 / 2, 1, 2, 4, 8, 16, 32, 64] if lo <= t <= hi]
    axes[0, 0].set_yticks(ticks, [f"{t:g}×" for t in ticks])
    axes[0, 0].set_xticks([1, 2, 4, 8, 16], ["1", "2", "4", "8", "16"])
    fig.suptitle("Small operands: BigMath vs GMP 6.3.0 on Apple M1 Max", fontsize=16, fontweight="bold", y=0.98)
    fig.text(0.5, 0.915, "Time relative to GMP  ·  Below 1× means BigMath is faster  ·  "
             "Values of up to 4 limbs are now stored inline", ha="center", fontsize=10.5)
    fig.supxlabel("Limbs per operand (64-bit limbs)", fontsize=11, y=0.07)
    fig.text(0.5, 0.02, "November 21, 2026  ·  1 thread  ·  BigMath operators return a new value; "
             "GMP reuses its output", ha="center", fontsize=9, color="#384655")
    fig.subplots_adjust(left=0.07, right=0.98, top=0.8, bottom=0.155, hspace=0.32, wspace=0.12)
    fig.savefig(path)
    plt.close(fig)
    print("wrote", path)


def progress_plot(path):
    fig, ax = plt.subplots(figsize=(8.4, 4.6), dpi=150)
    labels = [r[0] for r in BEFORE_AFTER]
    before = [r[1] for r in BEFORE_AFTER]
    after = [r[2] for r in BEFORE_AFTER]
    x = range(len(labels))
    w = 0.38
    ax.bar([i - w / 2 for i in x], before, w, label="2026-05-30 baseline",
           color="#bbbbbb")
    ax.bar([i + w / 2 for i in x], after, w, label="after PRs #82–#119",
           color="#1f77b4")
    ax.axhline(1.0, color="black", linewidth=1.0, linestyle="--", alpha=0.7)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("BigMath / GMP ratio  (lower is better)")
    ax.set_title("The 2026-06 optimization runs: PRs #82–#99 + #107–#114 + #116–#119",
                 fontsize=11)
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    print("wrote", path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ratio-only", action="store_true", help="Generate both GMP comparison figures and skip the history figure")
    args = parser.parse_args()
    ratio_plot(HERE / "bigmath_vs_gmp_single_thread.png",
               BENCHMARKS / "2026-11-21-factored-mfa-m1-max-single-thread.csv",
               "Single-thread BigMath vs GMP 6.3.0 on Apple M1 Max", "1 thread per library",
               legend_location="best")
    ratio_plot(HERE / "bigmath_vs_gmp.png",
               BENCHMARKS / "2026-11-21-factored-mfa-m1-max.csv",
               "BigMath vs GMP 6.3.0 on Apple M1 Max", "Up to 8 BigMath threads, 1 GMP thread")
    small = BENCHMARKS / "2026-11-21-small-operands.csv"
    if small.exists():
        small_ops_plot(HERE / "bigmath_vs_gmp_small.png", small)
    if not args.ratio_only:
        progress_plot(HERE / "session_2026_06_11.png")
