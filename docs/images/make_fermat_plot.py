"""Plot the local Phase 2 comparison against forced wide-CRT blocks."""

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "benchmarks" / "2026-11-22-fermat-outer-m1-max.csv"
MODES = (
    ("default", "Up to 8 threads", "#0072B2", "o"),
    ("single", "Single thread", "#B34815", "s"),
)


def main():
    with RESULTS.open(newline="") as source:
        rows = {(int(row["a_digits"]), row["library"]): row for row in csv.DictReader(source)}
    sizes = sorted({size for size, _ in rows})

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.7), layout="constrained")
    for ax, field, title in zip(axes, ("warm_ms", "warm_peak_bytes"),
                                ("Warm multiplication time", "Peak resident memory")):
        ax.axhline(1, color="#596777", linestyle="--", linewidth=1)
        for mode, label, color, marker in MODES:
            measured = [size for size in sizes if (size, f"fermat_{mode}_auto") in rows]
            ratios = [float(rows[size, f"fermat_{mode}_auto"][field]) /
                      float(rows[size, f"wide_blocks3_{mode}"][field]) for size in measured]
            ax.plot(measured, ratios, color=color, marker=marker, label=label)
        ax.set_xscale("log")
        ax.set_xticks(sizes, [f"{n / 1e9:g}B" if n >= 1e9 else f"{n / 1e6:g}M" for n in sizes])
        ax.minorticks_off()
        ax.set_xlabel("Decimal digits per operand")
        ax.set_ylabel("Fermat / forced blocks")
        ax.set_title(title, loc="left", weight="bold")
        ax.grid(axis="y", alpha=0.22)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Fermat outer FFT on the 32 GiB M1 Max", fontsize=14, weight="bold")
    fig.legend(*axes[0].get_legend_handles_labels(), loc="outside lower center", ncol=2, frameon=False)
    fig.savefig(HERE / "bigmath_fermat_outer.png", dpi=170)
    plt.close(fig)


if __name__ == "__main__":
    main()
