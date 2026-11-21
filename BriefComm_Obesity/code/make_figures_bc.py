# -*- coding: utf-8 -*-
"""
make_figures_bc.py -- Brief Communication (Obesity) Figure 1 / Figure 2
Data are hard-coded from the project's result CSVs; each value is annotated
with its source file so the figure can be re-derived if a CSV changes.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "figures")
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 7.5,
    "axes.linewidth": 0.7,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

INK = "#1a1a1a"
GREY = "#8a8a8a"
RED = "#b2182b"
BLUE = "#2166ac"
TEAL = "#1b7837"
SAND = "#c9b79c"


def panel_tag(ax, letter, x=-0.30, y=1.05):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="bottom", ha="left", color=INK)


# ----------------------------------------------------------------------------
# Figure 1 -- three ways to manufacture the alarm
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(6.5, 2.35))
gs = fig.add_gridspec(1, 3, wspace=0.55, left=0.075, right=0.985,
                      top=0.86, bottom=0.17)

# --- panel a: confounding -------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
labels = ["Crude", "Causal forest", "Parametric"]
vals = [1.87, -0.33, -0.52]          # S1_subject_level / E2-32 / table2_main
ypos = [2, 1, 0]
ax.axvline(0, color=GREY, lw=0.8, zorder=1)
ax.hlines(ypos[0], 0, vals[0], color=RED, lw=1.4, zorder=2)
ax.hlines(ypos[1:], 0, vals[1:], color=BLUE, lw=1.4, zorder=2)
ax.plot([vals[0]], [ypos[0]], "o", ms=6, color=RED, zorder=3)
ax.plot(vals[1:], ypos[1:], "o", ms=6, color=BLUE, zorder=3)
for y, v in zip(ypos, vals):
    ax.text(v + (0.24 if v > 0 else -0.20), y, f"{v:+.2f}",
            va="center", ha="left" if v > 0 else "right", fontsize=7.5, color=INK)
ax.annotate("", xy=(0.0, 2.55), xytext=(1.87, 2.55),
            arrowprops=dict(arrowstyle="<->", color=INK, lw=0.7))
ax.text(0.94, 2.64, "shift 2.20 years, opposite direction", ha="center",
        va="bottom", fontsize=6.6, color=INK, style="italic")
ax.set_yticks(ypos)
ax.set_yticklabels(labels)
ax.set_ylim(-0.55, 3.05)
ax.set_xlim(-1.5, 2.8)
ax.set_xlabel("Difference in PhenoAge (years)")
ax.set_title("Confounding", fontsize=8.5, pad=4)
panel_tag(ax, "a")

# --- panel b: prescribing -------------------------------------------------
ax = fig.add_subplot(gs[0, 1])
groups = ["All men", "Drugs excluded"]
lo = [-4.6, 0.1]
hi = [-3.4, 0.8]
mid = [np.mean([-4.6, -3.4]), np.mean([0.1, 0.8])]
ypos = [1, 0]
ax.axvline(0, color=GREY, lw=0.8, zorder=1)
for y, l, h, m, c in zip(ypos, lo, hi, mid, [RED, BLUE]):
    ax.hlines(y, l, h, color=c, lw=2.6, alpha=0.85, zorder=2)
    ax.plot([m], [y], "|", ms=11, color=INK, mew=1.2, zorder=3)
ax.text(-3.9, 1.24, "P=0.012\u20130.091", ha="center", fontsize=6.8, color=RED)
ax.text(0.45, 0.24, "P=0.383\u20130.670", ha="center", fontsize=6.8, color=BLUE)
ax.set_yticks(ypos)
ax.set_yticklabels(groups)
ax.set_ylim(-0.55, 1.65)
ax.set_xlim(-5.4, 2.2)
ax.set_xlabel("\u0394CATE, highest vs lowest PFAS decile (pp)")
ax.set_title("Prescribing (167/2,244 = 7.4%)", fontsize=8.5, pad=4)
panel_tag(ax, "b")

# --- panel c: BMI-bearing formulas ---------------------------------------
ax = fig.add_subplot(gs[0, 2])
pairs = [("HSI", 1.2e-6, 0.68, 2), ("FLI", 0.017, 0.14, 1)]
ax.axvline(0.05, color=INK, ls="--", lw=0.7)
for nm, before, after, y in pairs:
    ax.plot([before, after], [y, y], "-", color=GREY, lw=0.9, zorder=1)
    ax.plot([before], [y], "o", ms=6, color=RED, zorder=3)
    ax.plot([after], [y], "o", ms=6, color=BLUE, zorder=3)
ax.plot([0.61], [0], "o", ms=6, color=TEAL, zorder=3)
ax.set_yticks([2, 1, 0])
ax.set_yticklabels(["HSI", "FLI", "ALT"])
ax.set_ylim(-0.6, 2.6)
ax.set_xscale("log")
ax.set_xlim(2e-7, 30)
ax.set_xticks([1e-6, 1e-4, 1e-2, 1])
ax.set_xlabel("Interaction P (log scale)")
ax.set_title("BMI-bearing formulas", fontsize=8.5, pad=4)
ax.text(0.055, 2.42, "P=0.05", fontsize=6.4, ha="left", va="bottom")
ax.text(1.35, 0.0, "no BMI\nterm", fontsize=6.0, va="center", ha="left", color=TEAL)
panel_tag(ax, "c")

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"fig1_artefacts.{ext}"), dpi=300,
                bbox_inches="tight", facecolor="white")
plt.close(fig)

# ----------------------------------------------------------------------------
# Figure 2 -- the signal that survives
# ----------------------------------------------------------------------------
fig = plt.figure(figsize=(6.5, 2.5))
gs = fig.add_gridspec(1, 2, width_ratios=[1.45, 1.0], wspace=0.42,
                      left=0.115, right=0.985, top=0.86, bottom=0.19)

# --- panel a: forest of the androgen panel (table3_androgen_axis.csv) -----
ax = fig.add_subplot(gs[0, 0])
rows = [
    ("Albumin (control)", 0.20, -0.16, 0.56, "control"),
    ("Free androgen index", -0.72, -2.81, 1.41, "bio"),
    ("Free testosterone", 1.62, -0.65, 3.94, "bio"),
    ("%Free testosterone", -2.31, -3.21, -1.40, "bio"),
    ("SHBG", 4.78, 3.15, 6.43, "binding"),
    ("Total testosterone", 4.02, 1.70, 6.39, "total"),
]
colmap = {"control": GREY, "bio": BLUE, "binding": SAND, "total": RED}
ax.axvline(0, color=GREY, lw=0.8, zorder=1)
for i, (nm, est, lo, hi, kind) in enumerate(rows):
    c = colmap[kind]
    ax.plot([lo, hi], [i, i], "-", color=c, lw=1.1, zorder=2)
    ax.plot([est], [i], "o", ms=5.5, color=c, zorder=3)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([r[0] for r in rows])
ax.set_xlim(-4.5, 7.5)
ax.set_xlabel("Change per 5% weight loss (%)")
ax.set_title("Dose-dependent androgen response", fontsize=8.5, pad=4)
ax.text(-0.16, 1.06, "a", transform=ax.transAxes, fontsize=10,
        fontweight="bold", va="bottom")
ax.set_ylim(-0.6, 5.6)

# --- panel b: the arithmetic ---------------------------------------------
ax = fig.add_subplot(gs[0, 1])
bars = [("Total T", 4.02, RED), ("SHBG", 4.78, SAND), ("Albumin", 0.20, GREY)]
x = np.arange(3)
ax.bar(x, [b[1] for b in bars], width=0.55,
       color=[b[2] for b in bars], edgecolor="none")
for xi, b in zip(x, bars):
    ax.text(xi, b[1] + 0.22, f"{b[1]:+.1f}", ha="center", fontsize=7)
ax.set_xticks(x)
ax.set_xticklabels([b[0] for b in bars])
ax.set_ylim(0, 6.6)
ax.set_ylabel("Change per 5% (%)")
ax.set_title("Why the free fraction falls", fontsize=8.5, pad=4)
ax.text(-0.20, 1.06, "b", transform=ax.transAxes, fontsize=10,
        fontweight="bold", va="bottom")
ax.set_xlim(-0.65, 2.6)
ax.set_ylim(0, 6.4)

for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"fig2_androgen.{ext}"), dpi=300,
                bbox_inches="tight", facecolor="white")
plt.close(fig)

print("written:", sorted(os.listdir(OUT)))
