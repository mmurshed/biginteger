#!/usr/bin/env python3
"""Generate the BENCHMARK.md / README.md performance figures.

Data: 2026-06-12 evening canonical run at v13.0 (post #120-#126: audit fixes,
CRT squaring, dedups, header hygiene). M1 Max, paired same-run BigMath vs
GMP 6.3, quiet machine; 50M-200M mul rows and 5M-20M ToString rows use warm
steady-state ratios (first call discarded, best-of-3) — see BENCHMARK.md
methodology notes. Regenerate by rerunning bench_vs_gmp and updating the
tables below, then:
python3 docs/images/make_benchmark_plots.py
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# digits -> BigMath/GMP wall-clock ratio (lower is better; < 1 = BigMath faster)
MUL_BAL = {
    1e3: 2.33, 5e3: 2.56, 1e4: 3.19, 5e4: 1.24, 1e5: 1.29, 5e5: 0.89,
    1e6: 0.85, 2e6: 0.79, 5e6: 0.52, 1e7: 0.31, 2e7: 0.47,
    5e7: 0.89, 1e8: 0.91, 2e8: 0.97,  # 50M+ warm steady-state
}
DIV_SKEW = {  # a = 5b shapes, keyed by dividend digits
    4e4: 2.64, 1e5: 3.81, 2e5: 2.31, 5e5: 1.59, 1e6: 1.44, 2e6: 1.18,
    5e6: 1.00, 1e7: 1.18, 2e7: 0.82, 5e7: 0.45, 1e8: 0.49, 2e8: 0.69,
}
PARSE = {
    1e3: 1.55, 1e4: 2.26, 5e4: 2.77, 1e5: 1.62, 5e5: 0.69, 1e6: 0.60,
    2e6: 0.56, 5e6: 0.44, 1e7: 0.40, 2e7: 0.38, 5e7: 0.53,
}
TOSTR = {
    1e3: 1.83, 1e4: 3.37, 5e4: 2.81, 1e5: 3.38, 2e5: 2.27, 5e5: 1.62,
    1e6: 1.29, 2e6: 1.07, 5e6: 0.58, 1e7: 0.54, 2e7: 0.49,  # 5M+ warm
}

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
