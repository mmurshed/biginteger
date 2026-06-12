#!/usr/bin/env python3
"""Generate the BENCHMARK.md / README.md performance figures.

Data: 2026-06-12 full refresh (M1 Max, paired same-run BigMath vs GMP 6.3,
quiet machine). Regenerate by rerunning the refresh harness and updating the
tables below, then: python3 docs/images/make_benchmark_plots.py
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# digits -> BigMath/GMP wall-clock ratio (lower is better; < 1 = BigMath faster)
MUL_BAL = {
    1e3: 2.18, 5e3: 2.68, 1e4: 3.31, 5e4: 1.16, 1e5: 1.09, 5e5: 0.59,
    1e6: 0.65, 2e6: 0.58, 5e6: 0.41, 1e7: 0.49, 2e7: 1.05,
}
DIV_SKEW = {  # a = 5b shapes, keyed by dividend digits
    4e4: 3.38, 1e5: 4.92, 2e5: 2.23, 5e5: 1.64, 1e6: 1.43, 2e6: 1.18,
    5e6: 1.03, 1e7: 0.94, 2e7: 1.02,
}
PARSE = {
    1e3: 1.58, 1e4: 2.26, 5e4: 2.64, 1e5: 2.37, 5e5: 1.64, 1e6: 1.55,
    2e6: 1.45, 5e6: 1.20, 1e7: 1.11, 2e7: 1.09,
}
TOSTR = {
    1e3: 1.83, 1e4: 3.35, 5e4: 3.06, 1e5: 4.31, 2e5: 3.29, 5e5: 2.39,
    1e6: 2.05, 2e6: 1.77, 5e6: 1.33, 1e7: 1.19, 2e7: 1.16,
}

# Session progress: BigMath/GMP ratio at the start of the 2026-06-11 session
# (2026-05-30 baseline tables) vs after PRs #82-#99.
BEFORE_AFTER = [
    ("mul 1M×1M", 1.15, 0.65),
    ("mul 5M×5M", 0.75, 0.41),
    ("div 1M×200k", 3.36, 1.43),
    ("div 5M×1M", 2.87, 1.03),
    ("div 10M×2M", 2.79, 0.94),
    ("tostr 1M", 4.48, 2.05),
    ("tostr 20M", 2.52, 1.16),
    ("parse 1M", 2.33, 1.55),
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
    ax.fill_between([7e2, 3e7], 0.3, 1.0, color="green", alpha=0.06)
    ax.set_xscale("log")
    ax.set_xlim(7e2, 3e7)
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
    ax.bar([i + w / 2 for i in x], after, w, label="after PRs #82–#99",
           color="#1f77b4")
    ax.axhline(1.0, color="black", linewidth=1.0, linestyle="--", alpha=0.7)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel("BigMath / GMP ratio  (lower is better)")
    ax.set_title("One optimization day: 2026-06-11 → 06-12 (PRs #82–#99)",
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
