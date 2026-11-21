#!/usr/bin/env python3
"""Generate the BENCHMARK.md / README.md performance figures.

Ratio figure: canonical BENCHMARK.md tables, multiplication 2026-11-20;
division and decimal I/O 2026-11-21. M1 Max, paired same-run BigMath vs
GMP 6.3; BigMath uses up to 8 threads, GMP is single-threaded. Multiplication
at 50M-200M digits uses warm steady-state ratios. ToString cold and the three
measured warm ratios are shown separately. Ratios are copied from the tables
rather than recomputed from rounded timings; see BENCHMARK.md for methodology.
The session-progress figure retains its historical June data.
Regenerate after updating the benchmark tables and ratios below:
python3 docs/images/make_benchmark_plots.py
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# digits -> BigMath/GMP wall-clock ratio (lower is better; < 1 = BigMath faster)
MUL_BAL = {
    1e3: 1.39, 5e3: 1.48, 1e4: 0.90, 5e4: 0.63, 1e5: 0.42, 5e5: 0.60,
    1e6: 0.50, 2e6: 0.55, 5e6: 0.31, 1e7: 0.21, 2e7: 0.34,
    5e7: 0.63, 1e8: 0.68, 2e8: 0.72,  # 50M+ warm steady-state
}
DIV_SKEW = {  # 4:1, 5:1 and 10:1 digit ratios, keyed by dividend digits
    4e4: 1.36, 1e5: 2.01, 2e5: 0.67, 5e5: 0.56, 1e6: 0.49, 2e6: 0.35,
    5e6: 0.36, 1e7: 0.35, 2e7: 0.30, 5e7: 0.19, 1e8: 0.27, 2e8: 0.32,
}
PARSE = {
    1e3: 1.49, 1e4: 1.71, 5e4: 1.37, 1e5: 0.54, 5e5: 0.37, 1e6: 0.34,
    2e6: 0.32, 5e6: 0.25, 1e7: 0.24, 2e7: 0.23, 5e7: 0.28,
}
TOSTR_COLD = {
    1e3: 1.80, 1e4: 2.33, 5e4: 1.80, 1e5: 1.37, 2e5: 0.74, 5e5: 0.52,
    1e6: 0.42, 2e6: 0.38, 5e6: 0.44, 1e7: 0.40, 2e7: 0.39,
}
TOSTR_WARM = {1e5: 0.57, 1e6: 0.30, 5e6: 0.32}

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


def ratio_plot(path):
    fig, ax = plt.subplots(figsize=(8.8, 5.6), dpi=180)
    series = [
        (MUL_BAL, "multiply (balanced; ≥50M warm)", "#1f77b4", "o"),
        (DIV_SKEW, "divide (4:1–10:1 digit ratios)", "#d62728", "s"),
        (PARSE, "parse (str → int)", "#2ca02c", "^"),
        (TOSTR_COLD, "ToString (cold)", "#9467bd", "d"),
    ]
    for data, label, color, marker in series:
        xs = sorted(data)
        ax.plot(xs, [data[x] for x in xs], label=label, color=color,
                marker=marker, markersize=4.5, linewidth=1.6)
    xs = sorted(TOSTR_WARM)
    ax.plot(xs, [TOSTR_WARM[x] for x in xs], label="ToString (warm; measured sizes)",
            color="#9467bd", marker="D", markersize=6, markerfacecolor="white",
            markeredgewidth=1.5, linestyle="none", zorder=3)
    ax.axhline(1.0, color="black", linewidth=1.0, linestyle="--", alpha=0.7)
    ax.text(1.25e3, 1.04, "GMP parity", fontsize=8.5, alpha=0.8)
    ax.text(1.25e3, 0.12, "BigMath faster", fontsize=8.5, color="#35603a")
    ax.axhspan(0, 1.0, color="green", alpha=0.06)
    ax.set_xscale("log")
    ax.set_xlim(7e2, 3e8)
    ax.set_ylim(0, 2.6)
    ax.set_xlabel("operand size (decimal digits; division: dividend size)")
    ax.set_ylabel("BigMath / GMP wall-clock ratio   (lower is better)")
    ax.set_title("BigMath vs GMP 6.3 — Apple M1 Max\n"
                 "2026-11-20: multiplication · 2026-11-21: division and decimal I/O",
                 fontsize=11)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper right", fontsize=8.5)
    fig.text(0.5, 0.02,
             "BigMath: up to 8 threads · GMP: single-threaded. "
             "Source and timing methodology: BENCHMARK.md",
             ha="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.045, 1, 1))
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
    import os
    here = os.path.dirname(os.path.abspath(__file__))
    ratio_plot(os.path.join(here, "bigmath_vs_gmp.png"))
    progress_plot(os.path.join(here, "session_2026_06_11.png"))
