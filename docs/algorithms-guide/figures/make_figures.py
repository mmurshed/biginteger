#!/usr/bin/env python3
"""Generate the data figures for the BigMath algorithms guide.

Every chart here is either computed from the library's dispatch rules
(src/algorithms/Multiplication.cpp, src/algorithms/Division.cpp,
include/biginteger/build/DispatchThresholds.h) or plotted from measurements
already recorded in docs/MULTIPLICATION.md and docs/DIVISION.md.
If a threshold changes, update the constants below and re-run:

    python3 docs/algorithms-guide/figures/make_figures.py
"""

import math
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import numpy as np

OUT = os.path.dirname(os.path.abspath(__file__))

# Reference palette (light mode), fixed slot order.
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df"

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 9,
    "axes.edgecolor": MUTED,
    "axes.labelcolor": INK,
    "axes.titlesize": 10,
    "axes.titlecolor": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.frameon": False,
    "lines.linewidth": 2.0,
    "lines.markersize": 5,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})

# --- Multiplication dispatch constants (64-bit limbs) -----------------------
CLASSIC_SUM = 56
CLASSIC_MIN = 24
NTT_SUM = 640
NTT_MIN = 64
NTT_AREA = 96000
NTT_COLD_AREA = 262144


def mul_algo(la, lb, warm):
    """0 = schoolbook, 1 = Karatsuba family, 2 = NTT. Mirrors Multiply()."""
    if la == 1 or lb == 1:
        return 0
    s, mn = la + lb, min(la, lb)
    if s <= CLASSIC_SUM or mn <= CLASSIC_MIN:
        return 0
    area = la * lb
    if s >= NTT_SUM and mn >= NTT_MIN and area >= NTT_AREA:
        if area >= NTT_COLD_AREA or warm:
            return 2
        return 1
    return 1


def fig_mul_map():
    xs = np.unique(np.round(np.logspace(0, math.log10(4096), 400)).astype(int))
    ys = np.unique(np.round(np.logspace(0, math.log10(65536), 500)).astype(int))
    cmap = ListedColormap([AQUA, BLUE, ORANGE])
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6), sharey=True)
    for ax, warm, title in ((axes[0], True, "Thread pool awake (warm)"),
                            (axes[1], False, "Thread pool parked (cold) or 1 thread")):
        grid = np.full((len(ys), len(xs)), np.nan)
        for i, y in enumerate(ys):
            for j, x in enumerate(xs):
                if y >= x:
                    grid[i, j] = mul_algo(x, y, warm)
        ax.pcolormesh(xs, ys, grid, cmap=cmap, vmin=0, vmax=2, shading="nearest", rasterized=True)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log", base=2)
        ax.set_xlim(1, 4096)
        ax.set_ylim(1, 65536)
        ax.set_title(title)
        ax.set_xlabel("shorter operand (limbs)")
        ax.grid(False)
        ax.text(1.9, 6000, "schoolbook", color="white", fontsize=9, weight="bold", rotation=90, va="center")
        ax.text(60, 300, "Karatsuba\n(+ Toom-3\nlevels)", color="white", fontsize=9, weight="bold", ha="center")
        ax.text(700, 9000, "NTT", color="white", fontsize=10, weight="bold", ha="center")
        ax.plot([1, 4096], [1, 4096], color=INK2, lw=0.8, ls="--")
        # Example shapes discussed in the text.
        for (x, y, lab) in ((100, 100, "A"), (400, 400, "B"), (64, 2000, "C"), (10, 5000, "D")):
            ax.plot([x], [y], marker="o", ms=6, mfc="white", mec=INK, mew=1.2)
            ax.annotate(lab, (x, y), xytext=(5, 3), textcoords="offset points", fontsize=8, color=INK)
    axes[0].set_ylabel("longer operand (limbs)")
    fig.savefig(os.path.join(OUT, "mul_dispatch_map.pdf"))
    plt.close(fig)


# --- Division dispatch constants (64-bit limbs, quotient + remainder) -------
DIV_BZ_MIN_B = 96
DIV_NEWTON_MIN_B = 1536
DIV_NEWTON_SKEW_MIN_B = 768
DIV_NEWTON_SKEW_RATIO = 3
DIV_NEWTON_NEAR2_MIN_B = 3072
DIV_QSIZED_MIN_QUOTIENT = 3072
QS_GUARD = 4


def div_algo(na, nb):
    """0 schoolbook, 1 Burnikel-Ziegler, 2 Newton, 3 quotient-sized."""
    delta = na - nb
    if nb < DIV_BZ_MIN_B:
        return 0
    if delta >= nb:
        if nb >= DIV_NEWTON_MIN_B or (nb >= DIV_NEWTON_SKEW_MIN_B and na >= DIV_NEWTON_SKEW_RATIO * nb):
            return 2
        return 1
    if delta >= DIV_QSIZED_MIN_QUOTIENT and delta + 2 * QS_GUARD < nb:
        return 3
    if nb >= DIV_NEWTON_NEAR2_MIN_B and 4 * delta >= 3 * nb:
        return 2
    return 1


def fig_div_map():
    nbs = np.unique(np.round(np.logspace(math.log10(16), math.log10(262144), 420)).astype(int))
    ratios = np.linspace(1.0, 6.0, 400)
    grid = np.zeros((len(ratios), len(nbs)))
    for i, r in enumerate(ratios):
        for j, nb in enumerate(nbs):
            na = max(nb + 1, int(round(nb * r)))
            grid[i, j] = div_algo(na, nb)
    cmap = ListedColormap([AQUA, BLUE, ORANGE, VIOLET])
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ax.pcolormesh(nbs, ratios, grid, cmap=cmap, vmin=0, vmax=3, shading="nearest", rasterized=True)
    ax.set_xscale("log", base=2)
    ax.set_xlim(16, 262144)
    ax.set_ylim(1, 6)
    ax.grid(False)
    ax.set_xlabel("divisor length  $n_b$  (limbs)")
    ax.set_ylabel("dividend / divisor length")
    ax.text(40, 3.5, "school-\nbook\n(Knuth D)", color="white", fontsize=9, weight="bold", ha="center")
    ax.text(400, 3.5, "Burnikel–\nZiegler", color="white", fontsize=9, weight="bold", ha="center")
    ax.text(20000, 4.0, "Newton", color="white", fontsize=11, weight="bold", ha="center")
    ax.text(40000, 1.25, "quotient-sized", color="white", fontsize=9, weight="bold", ha="center")
    ax.text(700, 1.25, "Burnikel–Ziegler", color="white", fontsize=8, weight="bold", ha="center")
    ax.axhline(2.0, color="white", lw=0.6, ls=":")
    ax.text(17, 2.05, "ratio 2", color="white", fontsize=7)
    fig.savefig(os.path.join(OUT, "div_dispatch_map.pdf"))
    plt.close(fig)


# --- Limb-multiply counts ----------------------------------------------------
def count_karatsuba(n, leaf=24, toom=None, memo=None):
    if memo is None:
        memo = {}
    if n in memo:
        return memo[n]
    # Balanced n x n product; mirrors the split sizes of MultiplyRecursive.
    if n <= leaf:
        v = n * n
    elif toom is not None and n >= toom:
        k = (n + 2) // 3
        # v0 and vinf on the pieces, v(1), v(-1), v(2) on (k+1)-limb sums.
        v = (count_karatsuba(k, leaf, toom, memo) + count_karatsuba(n - 2 * k, leaf, toom, memo)
             + 3 * count_karatsuba(k + 1, leaf, toom, memo))
    else:
        m = (n + 1) // 2
        # low halves, high halves, and the (m+1)-limb sums.
        v = (count_karatsuba(m, leaf, toom, memo) + count_karatsuba(n - m, leaf, toom, memo)
             + count_karatsuba(m + 1, leaf, toom, memo))
    memo[n] = v
    return v


def fig_mul_counts():
    ns = np.unique(np.round(np.logspace(math.log10(16), math.log10(8192), 80)).astype(int))
    school = [1.0 for n in ns]
    kara = [count_karatsuba(int(n)) / (n * n) for n in ns]
    toom = [count_karatsuba(int(n), toom=168) / (n * n) for n in ns]
    fig, ax = plt.subplots(figsize=(5.6, 2.9))
    ax.plot(ns, school, color=AQUA, label="schoolbook ($n^2$)")
    ax.plot(ns, kara, color=BLUE, label="Karatsuba, leaf 24")
    ax.plot(ns, toom, color=VIOLET, label="Karatsuba + Toom-3 from 168")
    ax.set_xscale("log", base=2)
    ax.set_ylim(0, 1.1)
    ax.set_xlabel("operand length $n$ (limbs, both operands)")
    ax.set_ylabel("limb products / $n^2$")
    ax.legend(loc="center right")
    ax.text(5200, kara[-1] + 0.05, f"{kara[-1]:.0%}", color=INK2, fontsize=8, ha="center")
    ax.text(5200, toom[-1] - 0.09, f"{toom[-1]:.0%}", color=INK2, fontsize=8, ha="center")
    fig.savefig(os.path.join(OUT, "mul_counts.pdf"))
    plt.close(fig)


# --- Pool temperature (docs/MULTIPLICATION.md, "Pool temperature") ----------
def fig_pool():
    limbs = [320, 384, 448, 512, 640, 1024]
    kara = [23.9, 32.7, 42.3, 51.2, 72.2, 154.0]
    warm = [24.7, 23.6, 23.8, 25.0, 27.1, 44.9]
    cold = [52.1, 50.9, 52.6, 52.8, 59.4, 80.4]
    fig, ax = plt.subplots(figsize=(5.6, 2.9))
    ax.plot(limbs, kara, color=BLUE, marker="o", label="Karatsuba (+Toom-3)")
    ax.plot(limbs, warm, color=ORANGE, marker="o", label="NTT, pool awake")
    ax.plot(limbs, cold, color=AQUA, marker="o", label="NTT, pool parked")
    ax.axvspan(310, 512, color=GRID, alpha=0.6, lw=0)
    ax.text(415, 100, "decided by\npool state", ha="center", fontsize=8, color=INK2)
    ax.set_xlabel("limbs per operand (balanced product)")
    ax.set_ylabel("µs per multiply")
    ax.set_ylim(0, 160)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3)
    fig.savefig(os.path.join(OUT, "pool_temperature.pdf"))
    plt.close(fig)


# --- Truncated transform staircase ------------------------------------------
def fig_trunc():
    cs = np.arange(1024, 8193, 8)
    full, trunc = [], []
    for c in cs:
        n = 1 << int(math.ceil(math.log2(c)))
        full.append(n)
        H = n // 2
        if c > H and c <= H + H // 2 and H >= 2048:
            L = 1 << int(math.ceil(math.log2(c - H)))
            trunc.append(H + L)
        else:
            trunc.append(n)
    fig, ax = plt.subplots(figsize=(5.6, 2.8))
    ax.plot(cs, cs, color=MUTED, lw=1.0, ls="--", label="ideal (exactly $c$)")
    ax.plot(cs, full, color=AQUA, label="power-of-two length")
    ax.plot(cs, trunc, color=ORANGE, label="with $H+L$ truncation")
    ax.set_xlabel("coefficients in the product, $c$")
    ax.set_ylabel("transform work (points)")
    ax.legend(loc="lower right")
    fig.savefig(os.path.join(OUT, "trunc_staircase.pdf"))
    plt.close(fig)


# --- Newton convergence for 1/7 ----------------------------------------------
def fig_newton():
    from fractions import Fraction
    d = Fraction(7)
    x = Fraction(1, 10)
    errs = []
    for _ in range(6):
        errs.append(abs(float(Fraction(1, 7) - x)))
        x = x * (2 - d * x)
    digits = [-math.log10(e) for e in errs]
    fig, ax = plt.subplots(figsize=(5.2, 2.6))
    ax.bar(range(len(digits)), digits, color=BLUE, width=0.55)
    for i, v in enumerate(digits):
        ax.text(i, v + 0.4, f"{v:.1f}", ha="center", fontsize=8, color=INK)
    ax.set_xlabel("Newton step $k$")
    ax.set_ylabel("correct decimal digits")
    ax.set_ylim(0, max(digits) + 3)
    fig.savefig(os.path.join(OUT, "newton_convergence.pdf"))
    plt.close(fig)


if __name__ == "__main__":
    fig_mul_map()
    fig_div_map()
    fig_mul_counts()
    fig_pool()
    fig_trunc()
    fig_newton()
    print("figures written to", OUT)
