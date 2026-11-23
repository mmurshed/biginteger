#!/usr/bin/env python3
"""Plot peak resident memory and large-product timings from the isolated benchmark."""
import argparse
import csv
import math
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, NullFormatter

HERE = Path(__file__).resolve().parent
DEFAULT_CSV = HERE.parent / "benchmarks" / "2026-11-23-full-ed9a01c-m1-max-memory.csv"
CEILING_CSV = HERE.parent / "benchmarks" / "2026-11-21-crt-ceiling-m1-max.csv"
# Seeded balanced products remeasured after convolution residues replaced
# block products above 646M digits; these replace the seeded rows above.
LARGE_CSV = HERE.parent / "benchmarks" / "2026-11-21-wide-crt-m1-max-memory.csv"
COLORS = ("#2165a1", "#a54718", "#438456")
LIBRARIES = {"bigmath_default": "BigMath, up to 8 threads", "bigmath_single": "BigMath, 1 thread",
             "gmp": "GMP, 1 thread"}


def read_memory_rows(path=DEFAULT_CSV):
    with path.open(newline="") as src:
        return list(csv.DictReader(src))


def read_published_rows():
    """Manifest rows from DEFAULT_CSV and seeded balanced rows from LARGE_CSV."""
    rows = [r for r in read_memory_rows(DEFAULT_CSV) if r["source"] != "seeded_binary"]
    for r in read_memory_rows(LARGE_CSV):
        rows.append(dict(r, id="large-" + r["id"]))  # case ids restart in each sweep
    return rows


def axes_style(ax):
    ax.set_xscale("log")
    ax.set_yscale("log", base=2)
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:g}"))
    ax.grid(axis="y", color="#dce1e6", linewidth=0.7)
    ax.spines[["right", "top"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color("#a9b1b9")


def plot_ceiling(path):
    with CEILING_CSV.open(newline="") as src:
        times = {}
        for row in csv.DictReader(src):
            key = row["library"], row["variant"], int(row["a_digits"]) / 1e6
            times.setdefault(key, []).append(float(row["warm_ms"]) / 1000)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True, layout="constrained")
    for ax, lib in zip(axes, ["bigmath_single", "bigmath_default"]):
        for key, name, color, marker in [
                ((lib, "previous"), "Previous BigMath", "#747e88", "s"),
                ((lib, "current"), "Revised BigMath", COLORS[0], "o"),
                (("gmp", "reference"), "GMP, 1 thread", COLORS[2], "^")]:
            series = sorted((digits, median(v)) for (library, variant, digits), v in times.items()
                            if (library, variant) == key)
            ax.plot(*zip(*series), color=color, marker=marker, ms=4, label=name)
        for boundary in (323.23, 484.84):
            ax.axvline(boundary, ls=":", color="#838c96", lw=1)
        ax.set_title(LIBRARIES[lib], loc="left", fontsize=12, weight="bold")
        ax.set_xticks([320, 360, 400, 440, 480, 520])
        ax.set_ylim(bottom=0)
        ax.grid(axis="y", color="#dce1e6", lw=.6)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("Warm product time (seconds)")
    axes[0].legend(frameon=False, loc="upper left", bbox_to_anchor=(.08, 1))
    fig.supxlabel("Millions of decimal digits per operand", fontsize=10)
    fig.suptitle("Balanced multiplication across the 323M and 485M boundaries · Apple M1 Max",
                 fontsize=13, weight="bold")
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_ratios(rows, library, path, host, gmp_name, default_name):
    ref = {(r["id"], r["op"]): r for r in rows if r["library"] == "gmp"}
    data = [r for r in rows if r["library"] == library and r["status"] == "ok"]
    warm, first = "warm_peak_bytes", "first_peak_bytes"
    panels = [
        ("Multiplication", [("Equal lengths", "mul", warm, True), ("10:1 lengths", "mul", warm, False)]),
        ("Skewed division", [("Quotient", "div_q", warm, False), ("Quotient + remainder", "div_qr", warm, False)]),
        ("Decimal parsing", [("Through first call", "parse", first, None), ("Through warm calls", "parse", warm, None)]),
        ("Decimal formatting", [("Through first call", "tostr", first, None), ("Through warm calls", "tostr", warm, None)]),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(10.2, 6.6), sharex=True, sharey=True)
    for ax, (title, series) in zip(axes.flat, panels):
        axes_style(ax)
        ax.axhline(1, color="#424b57", ls="--", lw=1)
        for i, (label, op, field, equal) in enumerate(series):
            points = sorted((int(r["a_digits"]), int(r[field]) / int(ref[r["id"], op][field]))
                            for r in data if r["op"] == op and
                            (equal is None or (r["a_digits"] == r["b_digits"]) == equal))
            ax.plot(*zip(*points), label=label, color=COLORS[i], ls=["-", "--"][i],
                    marker=["o", "s"][i], ms=4, mfc=["white", COLORS[i]][i], lw=1.6)
        ax.set_title(title, loc="left", weight="bold", fontsize=12)
        ax.legend(frameon=False, fontsize=8, loc="best")
    ticks = [1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
    labels = ["1k", "10k", "100k", "1M", "10M", "100M"]
    largest = max(int(r["a_digits"]) for r in data)
    for tick, label in [(1e9, "1B"), (5e9, "5B")]:
        if largest >= tick:
            ticks.append(tick)
            labels.append(label)
    axes[0, 0].set_xticks(ticks, labels)
    axes[0, 0].set_xlim(700, max(3e8, largest * 1.35))
    fig.suptitle("Peak resident memory relative to GMP", fontsize=17, weight="bold", y=.98)
    name = default_name if library == "bigmath_default" else LIBRARIES[library]
    fig.text(.5, .908, "  ·  ".join(filter(None, (name, f"{gmp_name}, 1 thread", host))), ha="center")
    fig.supylabel("BigMath peak RSS / GMP peak RSS", fontsize=10, x=.025)
    fig.supxlabel("Decimal digits in the first operand", fontsize=11, y=.063)
    fig.text(.5, .022, "Fresh process per case  ·  Includes inputs and runtime  ·  Dashed horizontal line = equal memory",
             ha="center", fontsize=9)
    fig.subplots_adjust(top=.82, bottom=.16, left=.09, right=.98, hspace=.34, wspace=.14)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def plot_large(rows, path, host, default_name, *, show_wide_boundary=False):
    samples = [r for r in rows if r["op"] == "mul" and r["a_digits"] == r["b_digits"]
               and int(r["a_digits"]) >= 10_000_000]
    largest = max(int(r["a_digits"]) for r in samples) / 1e6
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), sharex=True)
    ceiling = 2**24 * 64 * math.log10(2) / 1e6
    # Shorter operand of 185,794,560 coefficients: the most one set of
    # signed residues resolves (CrtMaxBiasedTerms).
    wide_boundary = 185_794_560 * 32 * math.log10(2) / 1e6
    limited = any(r["status"] == "crt_length_limit" for r in samples)
    for ax, field, divisor, title in zip(axes, ["warm_ms", "warm_peak_bytes"], [1000, 2**30],
                                         ["Warm product time (seconds)", "Peak resident memory (GiB)"]):
        axes_style(ax)
        for i, (lib, label) in enumerate(LIBRARIES.items()):
            points = sorted((int(r["a_digits"])/1e6, float(r[field])/divisor) for r in samples
                            if r["library"] == lib and r["status"] == "ok")
            ax.plot(*zip(*points), marker=["o", "s", "^"][i], ms=4, color=COLORS[i],
                    lw=1.6, ls=["-", "--", "-."][i], label=default_name if lib == "bigmath_default" else label)
        ax.axvline(ceiling, color="#555", ls=":", lw=1.2)
        if show_wide_boundary:
            ax.axvline(wide_boundary, color="#555", ls="--", lw=1.0)
        if limited:
            ax.axvspan(ceiling, largest * 1.1, color="#f0f1f4", zorder=-2)
        if field == "warm_ms":
            lo, hi = ax.get_ylim()
            ax.set_yticks([t for t in (.01, .03, .1, .3, 1, 3, 10, 30, 100, 300) if lo <= t <= hi])
        ax.set_xlim(8.5, largest * 1.1)
        ax.set_title(title, loc="left", weight="bold", fontsize=12)
        ticks = [(10, "10M"), (20, "20M"), (50, "50M"), (100, "100M"), (200, "200M"), (500, "500M"),
                 (1000, "1B"), (2500, "2.5B"), (5000, "5B")]
        ax.set_xticks(*zip(*[(t, l) for t, l in ticks if t <= largest * 1.05]))
        ax.get_xticklabels()[-1].set_ha("left")
        ax.tick_params(labelsize=9)
    reach = "five billion" if largest >= 5000 else "one billion"
    fig.suptitle(f"Balanced multiplication through {reach} digits", fontsize=16, weight="bold", y=.97)
    fig.legend(*axes[0].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(.5,.90),
               ncol=3, frameon=False, fontsize=9)
    fig.supxlabel("Decimal digits per operand", fontsize=11, y=.09)
    note = "Shaded sizes exceed BigMath’s CRT limit (323.2M)" if limited else "Dotted line = single-transform CRT limit (323.2M)"
    if show_wide_boundary:
        note += f"  ·  Dashed line = wide CRT crossover ({wide_boundary / 1000:.2f}B)"
    fig.text(.5,.035, "  ·  ".join(filter(None, (host, note, "Points are measured products"))),
             ha="center", fontsize=9)
    fig.subplots_adjust(top=.74, bottom=.23, left=.077, right=.98, wspace=.18)
    fig.savefig(path, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, help="plot just this CSV in place of the published measurements")
    parser.add_argument("--outdir", type=Path, default=HERE)
    parser.add_argument("--host", help="machine name for a new measurement")
    parser.add_argument("--gmp-version", help="GMP version used for a new measurement")
    args = parser.parse_args()
    published = args.csv is None
    published_snapshot = published or args.csv.resolve() == DEFAULT_CSV.resolve()
    host = args.host or ("Apple M1 Max" if published_snapshot else "")
    version = args.gmp_version or ("6.3.0" if published_snapshot else "")
    gmp_name = f"GMP {version}".rstrip()
    default_name = LIBRARIES["bigmath_default"] if published_snapshot else "BigMath, default build"
    args.outdir.mkdir(parents=True, exist_ok=True)
    results = read_memory_rows(args.csv) if args.csv else read_published_rows()
    if published:
        plot_ceiling(args.outdir / "bigmath_crt_ceiling.png")
    for name, suffix in [("bigmath_default", ""), ("bigmath_single", "_single_thread")]:
        plot_ratios(results, name, args.outdir / f"bigmath_vs_gmp_memory{suffix}.png", host, gmp_name, default_name)
    if any(r["source"] == "seeded_binary" and r["a_digits"] == "1000000000" for r in results):
        plot_large(results, args.outdir / "bigmath_vs_gmp_large_products.png", host, default_name,
                   show_wide_boundary=published_snapshot)
