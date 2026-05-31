#!/usr/bin/env python3
"""Generate the near-balanced division dispatch speedup plot for BENCHMARK.md.

Reads the BEFORE/AFTER dispatch_ms numbers produced by
tests/performance/division_balanced_bench.cpp (run against main and against the
perf/newton-balanced-division branch) and renders a grouped bar chart of
dispatch wall-clock plus the per-shape speedup. Hard-coded with the measured
M1 Max numbers so the figure is reproducible without re-running the ~6-min
benchmark.

Usage: python3 docs/images/make_division_balanced_plot.py
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# shape label, divisor pow2?, BEFORE dispatch ms (main), AFTER dispatch ms (branch)
rows = [
    ("200k/100k\n(non-pow2)", False, 496.6, 249.9),
    ("500k/250k\n(non-pow2)", False, 5177.7, 588.3),
    ("1M/500k\n(non-pow2)", False, 10554.7, 1485.7),
    ("2M/1M\n(non-pow2)", False, 21006.9, 3747.6),
    ("512k/256k\n(pow2)", True, 656.4, 684.0),
    ("524290/262145\n(2^18+1)", False, 92060.5, 941.9),
]

labels = [r[0] for r in rows]
before = np.array([r[2] for r in rows])
after = np.array([r[3] for r in rows])
speedup = before / after

x = np.arange(len(labels))
w = 0.38

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), height_ratios=[3, 1.4],
                                sharex=True)

ax1.bar(x - w / 2, before, w, label="before (main → BZ)", color="#c0504d")
ax1.bar(x + w / 2, after, w, label="after (Newton balanced band)", color="#4f81bd")
ax1.set_ylabel("dispatch wall-clock (ms), log scale")
ax1.set_yscale("log")
ax1.set_title("Near-balanced (ratio 2) division dispatch — BigMath M1 Max, Base2_64\n"
              "Newton balanced band (b ≥ 98304) vs Burnikel-Ziegler")
ax1.legend(loc="upper left")
ax1.grid(True, axis="y", which="both", alpha=0.3)
for i, (b, a) in enumerate(zip(before, after)):
    ax1.text(i - w / 2, b, f"{b:.0f}", ha="center", va="bottom", fontsize=8)
    ax1.text(i + w / 2, a, f"{a:.0f}", ha="center", va="bottom", fontsize=8)

colors = ["#9bbb59" if s >= 1 else "#c0504d" for s in speedup]
ax2.bar(x, speedup, 0.5, color=colors)
ax2.axhline(1.0, color="black", lw=0.8, ls="--")
ax2.set_ylabel("speedup\n(before / after)")
ax2.set_xticks(x)
ax2.set_xticklabels(labels, fontsize=8)
ax2.grid(True, axis="y", alpha=0.3)
for i, s in enumerate(speedup):
    ax2.text(i, s, f"{s:.1f}×", ha="center", va="bottom", fontsize=8)

fig.tight_layout()
fig.savefig("docs/images/division_balanced_speedup.png", dpi=130)
fig.savefig("docs/images/division_balanced_speedup.svg")
print("wrote docs/images/division_balanced_speedup.{png,svg}")
