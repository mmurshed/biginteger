"""Plot representative paired comparisons from the memory optimization runs."""
import csv
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
CASES = [
    ("skew-default", "mul", 10_000_000, 1_000_000, "Skewed product\n10M × 1M"),
    ("weights-single", "mul", 200_000_000, 200_000_000, "Single-thread\n200M × 200M"),
    ("planes-default", "mul", 200_000_000, 200_000_000, "Threaded\n200M × 200M"),
    ("conversion-default", "tostr", 10_000_000, 0, "Format\n10M digits"),
]
fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0), layout="constrained")
for ax, field, title in zip(axes, ("warm_peak_bytes", "warm_ms"),
                            ("Peak resident memory", "Warm operation time")):
    for i, (file, op, digits, rhs, label) in enumerate(CASES):
        with (HERE.parent / "benchmarks" / f"memory-fix-{file}.csv").open() as f:
            rows = [r for r in csv.DictReader(f) if r["op"] == op and
                    int(r["a_digits"]) == digits and int(r["b_digits"]) == rhs]
        def value(build):
            return median(float(r[field]) for r in rows if r["build"] == build)
        for name, offset, color in (("before", -.18, "#94adc4"), ("after", .18, "#2165a1")):
            ratio = value(name) / value("gmp")
            ax.bar(i + offset, ratio, .34, color=color, label=name.capitalize() if i == 0 else None)
            ax.annotate(f"{ratio:.2g}×", (i + offset, ratio), xytext=(0, 3),
                        textcoords="offset points", ha="center", fontsize=8)
    ax.axhline(1, color="#55616b", ls="--", lw=1)
    ax.set_xticks(range(len(CASES)), [case[-1] for case in CASES], fontsize=8)
    ax.set_title(title, loc="left")
    ax.set_ylabel("Ratio to GMP")
    ax.margins(y=.19)
    ax.spines[["right", "top"]].set_visible(False)
axes[0].legend(frameon=False)
fig.suptitle("Memory optimizations on Apple M1 Max", weight="bold")
fig.savefig(HERE / "bigmath_memory_fixes.png", dpi=180)
