"""Figure 1 and 2 for the reworked Brief Communication."""
import os, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

BASE = os.environ.get("AGE22_ROOT", os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
OUT = os.path.join(BASE, "figures"); os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.linewidth": 0.7, "axes.spines.top": False, "axes.spines.right": False,
    "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "savefig.dpi": 400, "savefig.bbox": "tight",
})
C_CRUDE, C_ADJ, C_BAND, C_REF = "#9aa5b1", "#1f4e79", "#dbe6f0", "#b03a2e"

# ---------------------------------------------------------------- Figure 1
panels = [
    dict(key="PhenoAge", unit="years",
         est=[1.877, -0.180], se=[0.230, 0.206], margin=1.0,
         xlim=(-2.2, 2.6), xticks=[-2, -1, 0, 1, 2],
         note="harm >1 y excluded\nP = 3.5$\\times$10$^{-5}$", note_xy=(1.35, 0.62)),
    dict(key="Serum $\\Sigma$PFAS", unit="%",
         est=[-10.90, -1.243], se=[2.17, 1.50], margin=5.0,
         xlim=(-16, 8), xticks=[-15, -10, -5, 0, 5],
         note="rise >5% excluded\nP = 0.0060", note_xy=(3.0, 0.62)),
]

fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.55))
for ax, p in zip(axes, panels):
    ax.axvspan(-p["margin"], p["margin"], color=C_BAND, zorder=0, lw=0)
    ax.axvline(0, color="0.45", lw=0.7, ls=(0, (3, 2)), zorder=1)
    ys = [1, 0]
    for y, e, s, c, lab in [(1, p["est"][0], p["se"][0], C_CRUDE, "Crude"),
                            (0, p["est"][1], p["se"][1], C_ADJ, "Adjusted")]:
        ax.plot([e - 1.96*s, e + 1.96*s], [y, y], color=c, lw=1.6, solid_capstyle="butt", zorder=3)
        ax.plot(e, y, "o", ms=5.2, mfc=c, mec="white", mew=0.9, zorder=4)
        ax.text(e, y + (0.15 if y == 1 else -0.17),
                f"{e:+.2f}" if p["key"].startswith("Pheno") else f"{e:+.1f}",
                va="bottom" if y == 1 else "top", ha="center", fontsize=6.8, color=c)
    ax.annotate("", xy=(p["est"][1], 0.34), xytext=(p["est"][0], 0.34),
                arrowprops=dict(arrowstyle="->", lw=0.8, color="0.45",
                                connectionstyle="arc3,rad=0.12"))
    ax.set_yticks(ys); ax.set_yticklabels(["Crude", "Adjusted"], fontsize=7.5)
    ax.set_ylim(-0.45, 1.45); ax.set_xlim(*p["xlim"]); ax.set_xticks(p["xticks"])
    ax.set_xlabel(f"Difference per $\\geq$5% weight loss ({p['unit']})", fontsize=7.5)
    ax.set_title(p["key"], fontsize=8.5, loc="left", pad=6)
    ax.text(*p["note_xy"], p["note"], fontsize=6.4, color=C_REF, ha="left", va="center")
    ax.tick_params(labelsize=7)
axes[0].text(0.02, 0.03, "shaded band: equivalence margin ($\\pm$1 year)",
             transform=axes[0].transAxes, fontsize=6.2, color="0.4")
axes[1].text(0.02, 0.03, "shaded band: equivalence margin ($\\pm$5%)",
             transform=axes[1].transAxes, fontsize=6.2, color="0.4")
fig.tight_layout(w_pad=2.2)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"fig1_two_predictions.{ext}"))
plt.close(fig)

# ---------------------------------------------------------------- Figure 2
fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.5))

# (a) forest heterogeneity
ax = axes[0]
labels = ["All men\n(n=2,243)", "Excluding opioid\nusers (n=2,113)", "Excluding all 167\n(n=2,077)"]
pvals = [0.055, 0.548, 0.673]
cols = [C_CRUDE, C_ADJ, C_ADJ]
y = np.arange(3)[::-1]
ax.barh(y, pvals, height=0.52, color=cols, edgecolor="white", lw=0.6)
ax.axvline(0.05, color=C_REF, lw=0.8, ls=(0, (3, 2)))
ax.text(0.068, 2.56, "P = 0.05", fontsize=6.2, color=C_REF, va="bottom")
for yy, v in zip(y, pvals):
    ax.text(max(v + 0.025, 0.105), yy, f"{v:.2f}", va="center", fontsize=7, color="0.25")
ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=6.8)
ax.set_ylim(-0.55, 2.8); ax.set_xlim(0, 0.8); ax.set_xlabel("One-sided $P$ for forest heterogeneity", fontsize=7.5)
ax.set_title("(a) Prescription screening", fontsize=8.5, loc="left", pad=6)

# (b) formula decomposition
ax = axes[1]
names = ["Hepatic steatosis\nindex", "Fatty liver\nindex", "log ALT"]
before = [1.2e-6, 0.017, 0.61]
after  = [0.68, 0.14, 0.61]
x = np.arange(3); w = 0.34
ax.bar(x - w/2, before, w, color=C_CRUDE, edgecolor="white", lw=0.6, label="With BMI terms")
ax.bar(x + w/2, after,  w, color=C_ADJ,   edgecolor="white", lw=0.6, label="BMI terms removed")
ax.set_yscale("log"); ax.set_ylim(1e-7, 3)
ax.axhline(0.05, color=C_REF, lw=0.8, ls=(0, (3, 2)))
ax.text(-0.55, 0.115, "P = 0.05", fontsize=6.2, color=C_REF, ha="left", va="bottom")
ax.set_xticks(x); ax.set_xticklabels(names, fontsize=6.8)
ax.set_ylabel("Interaction $P$ (log scale)", fontsize=7.5)
ax.legend(frameon=False, fontsize=6.4, loc="upper center", ncol=2,
          handlelength=1.1, columnspacing=1.2)
ax.set_title("(b) Formula decomposition", fontsize=8.5, loc="left", pad=6)
fig.tight_layout(w_pad=2.4)
for ext in ("png", "pdf"):
    fig.savefig(os.path.join(OUT, f"fig2_three_routes.{ext}"))
plt.close(fig)
print("wrote", OUT)
