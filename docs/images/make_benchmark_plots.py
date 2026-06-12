#!/usr/bin/env python3
"""Generate the BENCHMARK.md / README.md performance figures.

Data: 2026-06-12 evening canonical run, post PR #107-#114 (M1 Max, paired same-run
BigMath vs GMP 6.3, quiet machine; 50M-200M mul rows use warm steady-state ratios,
see BENCHMARK.md methodology notes). Regenerate by rerunning bench_vs_gmp and
updating the tables below, then: python3 docs/images/make_benchmark_plots.py
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# digits -> BigMath/GMP wall-clock ratio (lower is better; < 1 = BigMath faster)
MUL_BAL = {
    1e3: 2.17, 5e3: 2.55, 1e4: 3.19, 5e4: 1.20, 1e5: 1.08, 5e5: 0.93,
    1e6: 0.81, 2e6: 0.79, 5e6: 0.51, 1e7: 0.32, 2e7: 0.47,
    5e7: 0.90, 1e8: 0.94, 2e8: 0.97,  # 50M+ warm steady-state
}
DIV_SKEW = {  # a = 5b shapes, keyed by dividend digits
    4e4: 3.35, 1e5: 4.84, 2e5: 2.26, 5e5: 1.65, 1e6: 1.46, 2e6: 1.15,
    5e6: 1.01, 1e7: 1.14, 2e7: 0.82, 5e7: 0.45, 1e8: 0.49, 2e8: 0.67,
}
PARSE = {
    1e3: 1.60, 1e4: 2.26, 5e4: 2.75, 1e5: 2.31, 5e5: 1.69, 1e6: 1.58,
    2e6: 1.46, 5e6: 1.19, 1e7: 1.10, 2e7: 1.01, 5e7: 1.04,
}
TOSTR = {
    1e3: 1.84, 1e4: 3.38, 5e4: 2.80, 1e5: 4.01, 2e5: 3.27, 5e5: 2.57,
    1e6: 2.16, 2e6: 1.85, 5e6: 1.28, 1e7: 1.17, 2e7: 1.47,  # 5M/10M warm
}

# Session progress: BigMath/GMP ratio at the 2026-05-30 baseline vs after the
# two 2026-06 runs (PRs #82-#99, then #107-#114).
BEFORE_AFTER = [
    ("mul 1M×1M", 1.15, 0.81),
    ("mul 10M×10M", 1.05, 0.32),
    ("mul 100M×100M", 1.80, 0.94),
    ("div 1M×200k", 3.36, 1.46),
    ("div 10M×2M", 2.79, 1.14),
    ("div 100M×20M", 2.85, 0.49),
    ("tostr 1M", 4.48, 2.16),
    ("parse 20M", 2.33, 1.01),
]


def ratio_plot(path):
    fig, ax = plt.subplots(figsize=(8.4, 5.0), dpi=150)
    series = [
        (MUL_BAL, "multiply (balanced)", "#1f77b4", "o"),
        (DIV_SKEW, "divide (5:1 skew)", "#d62728", "s"),
        (PARSE, "parse (str → int)", "#2ca02c", "^"),
        (TOSTR, "to-string (int → str)", "#9467bd", "d"),
    ]
    for data, label, color, marker in series:
        xs = sorted(data)
        ax.plot(xs, [data[x] for x in xs], label=label, color=color,
                marker=marker, markersize=4.5, linewidth=1.6)
    ax.axhline(1.0, color="black", linewidth=1.0, linestyle="--", alpha=0.7)
    ax.text(1.25e3, 0.93, "GMP parity", fontsize=8.5, alpha=0.8)
    ax.fill_between([7e2, 3e8], 0.3, 1.0, color="green", alpha=0.06)
    ax.set_xscale("log")
    ax.set_xlim(7e2, 3e8)
    ax.set_ylim(0.3, 5.2)
    ax.set_xlabel("operand size (decimal digits)")
    ax.set_ylabel("BigMath / GMP wall-clock ratio   (lower is better)")
    ax.set_title("BigMath vs GMP 6.3 — Apple M1 Max, 2026-06-12\n"
                 "below the dashed line, BigMath is faster than GMP",
                 fontsize=11)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
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
    ax.bar([i + w / 2 for i in x], after, w, label="after PRs #82–#114",
           color="#1f77b4")
    ax.axhline(1.0, color="black", linewidth=1.0, linestyle="--", alpha=0.7)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("BigMath / GMP ratio  (lower is better)")
    ax.set_title("The 2026-06 optimization runs: PRs #82–#99 + #107–#114",
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
