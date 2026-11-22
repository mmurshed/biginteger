"""Regenerate the historical PR #78 figure from the published M1 Max sweep."""
from pathlib import Path
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
# BENCHMARK.md at d6752bc, "MFA transpose fusion". The original sweep reports
# milliseconds and percentage gains separately, each rounded for publication.
# Columns: limbs per operand, baseline ms, fusion ms, reported speedup percent.
SAMPLES = [
    (3_000, 0.51, 0.54, -5.0),
    (10_000, 1.89, 1.92, -1.5),
    (50_000, 8.69, 9.08, -4.2),
    (200_000, 40.02, 40.06, -0.1),
    (700_000, 250.91, 241.81, 3.8),
    (1_300_000, 648.92, 649.04, 0.0),
    (2_000_000, 689.23, 678.93, 1.5),
    (2_600_000, 1121.25, 1011.62, 10.8),
    (4_000_000, 1159.97, 1046.13, 10.9),
    (5_200_000, 2379.93, 2163.77, 10.0),
    (6_000_000, 2393.66, 2181.19, 9.7),
    (8_000_000, 2473.24, 2231.99, 10.8),
    (10_000_000, 4978.10, 4636.81, 7.4),
]


def main():
    digits_per_limb = 64 * math.log10(2)
    digits = [limbs * digits_per_limb / 1e6 for limbs, *_ in SAMPLES]
    fig, (throughput, speedup) = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    for column, label, style in [(1, "Before fusion", "--"), (2, "With fusion", "-")]:
        rate = [d * 1000 / row[column] for d, row in zip(digits, SAMPLES)]
        throughput.plot(digits, rate, marker="o", ms=4, ls=style, label=label)
    throughput.set_ylabel("Millions of digits per second")
    throughput.set_title("Multiplication throughput", loc="left")
    throughput.legend(frameon=False)
    speedup.plot(digits, [row[3] for row in SAMPLES], marker="o", ms=4, color="#438456")
    speedup.axhline(0, linewidth=0.8, color="gray")
    speedup.set_ylabel("Speedup from fusion (%)")
    speedup.set_title("Change in throughput", loc="left")
    for ax in (throughput, speedup):
        ax.set_xscale("log")
        ax.axvline(2_000_000 * digits_per_limb / 1e6, ls=":", color="gray", label="First MFA sample")
        ax.set_xlabel("Millions of decimal digits per operand")
        ax.grid(axis="y", alpha=0.25)
        ax.spines[["top", "right"]].set_visible(False)
    speedup.legend(frameon=False)
    fig.suptitle("Transpose fusion on Apple M1 Max · PR #78", weight="bold")
    fig.savefig(HERE / "mfa_fusion_speedup.png", dpi=150)
    svg = HERE / "mfa_fusion_speedup.svg"
    with plt.rc_context({"svg.hashsalt": "bigmath-mfa-fusion", "svg.fonttype": "none"}):
        fig.savefig(svg, metadata={"Date": None})
    # Matplotlib leaves trailing spaces inside multiline SVG paths.
    svg.write_text("\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n")
    plt.close(fig)


if __name__ == "__main__":
    main()
