"""Plot the release-profile comparison from GitHub Actions run 45151452463."""
import csv
from collections import defaultdict
from pathlib import Path
from random import Random
from statistics import geometric_mean

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import NullLocator

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / "benchmarks" / "2026-11-23-platform-tuning-45151452463.csv"
PLATFORMS = {
    "linux-x86_64-gcc": "Linux x86-64 · GCC",
    "linux-x86_64-clang": "Linux x86-64 · Clang",
    "macos-arm64-clang": "macOS · Apple silicon",
    "macos-x86_64-clang": "macOS · Intel",
    "windows-x86_64-clang": "Windows x86-64 · Clang",
}
BLUE = "#0072B2"
ORANGE = "#B34815"
INK = "#19283B"


def load_ratios():
    ratios = defaultdict(list)
    with RESULTS.open(newline="") as source:
        for row in csv.DictReader(source):
            key = (row["platform"], row["operation"], row["state"])
            ratios[key].append(float(row["paired_speedup"]))
    return ratios


def format_axes(ax, ticks, limits):
    ax.axvline(1, color="#576675", ls="--", lw=1, zorder=2)
    ax.set_xticks(ticks, [f"{tick:g}×" for tick in ticks])
    ax.set_xlim(limits)
    ax.set_ylim(len(PLATFORMS) - 0.45, -0.65)
    ax.set_yticks(range(len(PLATFORMS)), PLATFORMS.values())
    ax.tick_params(length=0, pad=7)
    ax.xaxis.set_minor_locator(NullLocator())
    ax.grid(axis="x", color="#e3e8ee", linewidth=0.7)
    ax.set_axisbelow(True)
    ax.spines[["top", "right", "left", "bottom"]].set_visible(False)


def save_figure(fig, title, label, handles, filename):
    fig.subplots_adjust(left=0.185, right=0.985, bottom=0.235, top=0.79, wspace=0.16)
    fig.suptitle(title, y=0.96, fontsize=15, weight="bold")
    fig.supxlabel(label, x=0.585, y=0.035, fontsize=10)
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.585, 0.085),
               ncol=2, frameon=False, borderaxespad=0)
    fig.savefig(HERE / filename, dpi=170)
    plt.close(fig)


def plot_summary(ratios):
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 3.65), sharey=True)
    for ax, operation, title in zip(axes, ("mul-balanced", "square"),
                                    ("Balanced multiplication", "Squaring")):
        format_axes(ax, (0.9, 1, 1.1, 1.2, 1.3, 1.4), (0.86, 1.51))
        ax.set_title(title, loc="left", weight="bold")
        for state, color, marker, offset in (
            ("burst", BLUE, "o", -0.16), ("parked", ORANGE, "s", 0.16),
        ):
            values = [geometric_mean(ratios[platform, operation, state]) for platform in PLATFORMS]
            y = [index + offset for index in range(len(PLATFORMS))]
            ax.scatter(values, y, color=color, marker=marker, s=30, zorder=3,
                       label=state.capitalize())
            for value, position in zip(values, y):
                ax.annotate(f"{value:.2f}×", (value, position), xytext=(7, 0),
                            textcoords="offset points", va="center", color=color, fontsize=9)
    handles, _ = axes[0].get_legend_handles_labels()
    save_figure(fig, "Speedup from BigMath release profiles",
                "Geometric mean of baseline time ÷ tuned time", handles,
                "bigmath_platform_tuning.png")


def plot_cases(ratios):
    fig, axes = plt.subplots(1, 4, figsize=(13.2, 3.9), sharey=True)
    jitter = Random(44998103595)
    families = (
        (("mul-",), "Multiplication"),
        (("square",), "Squaring"),
        (("div-",), "Quotient and remainder"),
        (("add-", "sub-"), "Addition and subtraction"),
    )
    for ax, (family, title) in zip(axes, families):
        ax.set_xscale("log", base=2)
        format_axes(ax, (0.5, 1, 2, 3), (0.45, 3.85))
        ax.set_title(title, loc="left", weight="bold")
        for index, platform in enumerate(PLATFORMS):
            values = [ratio for (key, operation, _), sample in ratios.items()
                      if key == platform and operation.startswith(family) for ratio in sample]
            ax.scatter(values, [index + jitter.uniform(-0.27, 0.27) for _ in values],
                       s=11, color=BLUE, alpha=0.48, linewidths=0, zorder=3)
            ax.scatter(geometric_mean(values), index, color=ORANGE, marker="D", s=29,
                       edgecolor="white", linewidth=0.5, zorder=4)
    handles = [
        Line2D([], [], color=BLUE, linestyle="none", marker="o", markersize=4, label="One case"),
        Line2D([], [], color=ORANGE, linestyle="none", marker="D", markersize=5, label="Geometric mean"),
    ]
    save_figure(fig, "Speedup across individual validation cases",
                "Baseline time ÷ tuned time (log scale)", handles,
                "bigmath_platform_tuning_cases.png")


def main():
    plt.rcParams.update({"font.size": 10, "text.color": INK, "axes.labelcolor": INK,
                         "xtick.color": INK, "ytick.color": INK, "savefig.facecolor": "white"})
    ratios = load_ratios()
    plot_summary(ratios)
    plot_cases(ratios)


if __name__ == "__main__":
    main()
