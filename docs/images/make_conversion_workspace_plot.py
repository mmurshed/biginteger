"""Plot matched decimal conversions before and after the scratch-buffer change."""
import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "benchmarks" / "2026-11-21-decimal-conversion-scratch.csv"


def main():
    with RESULTS.open(newline="") as source:
        rows = [r for r in csv.DictReader(source)
                if r["op"] in ("parse", "tostr") and int(r["a_digits"]) >= 1_000_000]
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.8), layout="constrained")
    series = [
        ("before_threaded", "Up to 8 threads · before", "#2868a0", "--", "o"),
        ("after_threaded", "Up to 8 threads · after", "#2868a0", "-", "o"),
        ("before_single", "1 thread · before", "#ae541e", "--", "s"),
        ("after_single", "1 thread · after", "#ae541e", "-", "s"),
    ]
    for column, (field, unit) in enumerate([("warm_peak_bytes", "Peak RSS"), ("warm_ms", "Warm time")]):
        samples = defaultdict(list)
        for row in rows:
            samples[row["op"], int(row["a_digits"]), row["build"]].append(float(row[field]))
        for ax, op, title in zip(axes[:, column], ("parse", "tostr"), ("Parsing", "Formatting")):
            digits = sorted({n for operation, n, build in samples if operation == op})
            for build, label, color, style, marker in series:
                ratios = [median(samples[op, n, build]) / median(samples[op, n, "gmp"]) for n in digits]
                ax.plot(digits, ratios, ls=style, color=color, marker=marker, ms=4,
                        mfc="white" if style == "--" else color, label=label)
            ax.axhline(1, color="#525b62", ls=":", lw=1)
            ax.set_xscale("log")
            ax.set_xticks(digits, [f"{n // 1_000_000}M" for n in digits])
            ax.minorticks_off()
            ax.set_ylim(bottom=0)
            ax.set_ylabel(f"{unit} / GMP")
            ax.set_title(f"{title} · {unit}", loc="left", weight="bold")
            ax.grid(axis="y", alpha=0.25)
            ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Decimal conversion on Apple M1 Max", weight="bold", fontsize=16)
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc="outside lower center",
               ncol=2, frameon=False)
    for ax in axes[1]:
        ax.set_xlabel("Decimal digits")
    fig.savefig(HERE / "bigmath_conversion_scratch.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
