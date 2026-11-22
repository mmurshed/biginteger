"""Plot the paired, warm BigMath/GMP product timings from the M1 Max sweep."""
import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "benchmarks" / "2026-11-22-submillion-products-m1-max.csv"


def main():
    ratios = defaultdict(list)
    with DATA.open(newline="") as source:
        for row in csv.DictReader(source):
            digits = int(row["digits_each"])
            if digits <= 1_000_000:
                key = (int(row["bigmath_threads"]), row["build"], digits)
                ratios[key].append(float(row["warm_us"]) / float(row["gmp_warm_us"]))

    sizes = sorted({digits for _, _, digits in ratios})
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True, layout="constrained")
    for ax, threads, title in zip(axes, (1, 8), ("One BigMath thread", "Up to eight BigMath threads")):
        ax.axhline(1, color="#595f67", linewidth=1, linestyle=":")
        for build, label, color, marker in (("before", "Before", "#777777", "o"),
                                            ("after", "Optimized", "#0072B2", "s")):
            ax.plot(sizes, [median(ratios[threads, build, d]) for d in sizes],
                    color=color, marker=marker, markersize=4, label=label)
        ax.set_xscale("log")
        ax.set_xticks([1000, 5000, 20000, 100000, 500000, 1000000],
                      ["1k", "5k", "20k", "100k", "500k", "1M"])
        ax.set_xlabel("Decimal digits per operand")
        ax.set_title(title, loc="left", weight="bold")
        ax.grid(axis="y", alpha=0.2)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Warm time / GMP time")
    fig.suptitle("Balanced multiplication on Apple M1 Max", weight="bold")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center",
               ncol=2, frameon=False)
    fig.savefig(HERE / "bigmath_submillion_products.png", dpi=160)
    plt.close(fig)


if __name__ == "__main__":
    main()
