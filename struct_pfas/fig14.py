import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from geom import read_sdf, read_dht, shape, VDW
import csv, io

plt.rcParams.update({"font.family":"DejaVu Sans","font.size":8,
                     "axes.linewidth":0.8,"xtick.major.width":0.7,
                     "ytick.major.width":0.7,"axes.spines.top":False,
                     "axes.spines.right":False})

C_DHT, C_PFOA, C_PFOS = "#4C6A92", "#E08A2E", "#B33A3A"
C_HSA, RED = "#2E7D5B", "#B00020"

fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 2, hspace=0.50, wspace=0.34,
                      left=0.10, right=0.97, top=0.93, bottom=0.10)

# ================= A: 截面宽度 =================
axA = fig.add_subplot(gs[0, 0])
ligs = [("DHT", read_dht("1D2S.pdb"), C_DHT, "native steroid"),
        ("PFOA", read_sdf("pc_PFOA.sdf"), C_PFOA, "C8 perfluoro carboxylate"),
        ("PFOS", read_sdf("pc_PFOS.sdf"), C_PFOS, "C8 perfluoro sulfonate")]
SP = 6.0
for k, (nm, at, col, desc) in enumerate(ligs):
    y0 = k * SP
    s = shape(at)
    X = np.array([c for el, c in at if el != "H"]) - s["centroid"]
    Y = X @ s["axes"]
    axA.scatter(Y[:, 1], Y[:, 2] + y0, s=17, c=col, alpha=.85,
                edgecolors="white", linewidths=.35, zorder=3)
    e = s["ext"]
    axA.plot([-2.6, -2.6 + e[1]], [y0 - 2.15] * 2, color=col, lw=1.8,
             solid_capstyle="butt", zorder=2)
    axA.text(3.4, y0, f"{e[1]:.2f} Å", fontsize=7.4, color=col,
             va="center", ha="left", fontweight="bold")
    axA.text(-2.6, y0 + 2.05, desc, fontsize=6.0, color=col, style="italic", va="bottom")
axA.set_xlim(-3.0, 4.9); axA.set_ylim(-2.9, 3 * SP - 1.2)
axA.set_yticks([k * SP for k in range(3)])
axA.set_yticklabels(["DHT", "PFOA", "PFOS"], fontsize=8.5, fontweight="bold")
axA.set_xticks([]); axA.spines["bottom"].set_visible(False)
axA.set_title("a  Width at the thickest cross-section", loc="left",
              fontsize=8.5, fontweight="bold")

# ================= B: 口袋半径剖面 =================
axB = fig.add_subplot(gs[0, 1])
d = np.load("profile.npz", allow_pickle=True)
pro, dpr = d["prof"], d["dht_prof"]
z, wall = pro[:, 0], pro[:, 1]
axB.fill_between(z, wall, 6.25, color="#DDE3EA", alpha=.95, lw=0, zorder=1)
axB.plot(z, wall, color="#3D5A80", lw=1.5, zorder=3, label="pocket wall")
ok = ~np.isnan(dpr)
axB.plot(z[ok], dpr[ok], "o-", color=C_DHT, ms=3.2, lw=1.0, zorder=4, label="DHT (native)")
for nm, col, ls in [("PFOA", C_PFOA, "--"), ("PFOS", C_PFOS, ":")]:
    at = read_sdf(f"pc_{nm}.sdf"); s = shape(at)
    X = np.array([c for el, c in at if el != "H"]) - s["centroid"]
    Y = X @ s["axes"]
    r = np.array([VDW.get(el, 1.7) for el, _ in at if el != "H"])
    rr = (np.hypot(Y[:, 1], Y[:, 2]) + r).max()
    axB.axhline(rr, color=col, ls=ls, lw=1.2, zorder=5, label=f"{nm} {rr:.2f} Å")
axB.axhline(np.nanmin(wall), color=RED, lw=.9, ls="-.", alpha=.65, zorder=2)
axB.annotate(f"narrowest {np.nanmin(wall):.2f} Å", xy=(-4.6, np.nanmin(wall)),
             xytext=(-5.9, 3.05), fontsize=6.3, color=RED,
             arrowprops=dict(arrowstyle="-", lw=.6, color=RED))
axB.set_xlim(-6.3, 5.9); axB.set_ylim(1.5, 6.25)
axB.set_xlabel("position along steroid long axis (Å)", fontsize=7)
axB.set_ylabel("radius (Å)", fontsize=7)
axB.legend(fontsize=5.9, frameon=False, loc="lower right", handlelength=1.8,
           borderaxespad=.25, labelspacing=.26, ncol=2)
axB.set_title("b  Every PFAS fits — steric exclusion fails", loc="left",
              fontsize=8.5, fontweight="bold")

# ================= C: 电荷环境（横向） =================
axC = fig.add_subplot(gs[1, 0])
axC.barh([1, 0], [14, -1], height=.44, color=[C_HSA, C_DHT], alpha=.88, zorder=3)
axC.axvline(0, color="#333", lw=.9, zorder=4)
axC.set_yticks([1, 0])
axC.set_yticklabels(["HSA\nfatty-acid sites\n(1E7I)", "SHBG\nsteroid pocket\n(1D2S)"],
                    fontsize=7)
axC.text(14.6, 1, "+14", va="center", fontsize=11, fontweight="bold", color=C_HSA)
axC.text(-1.6, 0, "−1", va="center", ha="right", fontsize=11, fontweight="bold", color=C_DHT)
axC.text(6.6, 1.34, "6 Arg + 1 Lys  →  binds the −1 head group",
         fontsize=6.1, color=C_HSA, ha="center")
axC.text(-2.6, -0.34, "one Asp65  →  repels the −1 head group",
         fontsize=6.1, color=C_DHT, ha="left")
axC.set_xlim(-11, 20); axC.set_ylim(-0.62, 1.68)
axC.set_xlabel("net side-chain charge within 4.5 Å of the ligand", fontsize=7)
axC.set_title("c  Opposite electrostatic design", loc="left",
              fontsize=8.5, fontweight="bold")

# ================= D: 占据率 =================
axD = fig.add_subplot(gs[1, 1])
Kd = np.logspace(0, 4, 200)
for lab, L, f, col, ls in [("PFOA median, 1% free", 5.00, 0.01, C_PFOA, "-"),
                           ("PFOA median, 5% free", 5.00, 0.05, C_PFOA, "--"),
                           ("PFOS P95, 5% free", 68.18, 0.05, C_PFOS, "-")]:
    axD.plot(Kd, 100 * (L * f) / (Kd + L * f), color=col, ls=ls, lw=1.35, label=lab)
axD.axhline(4.0, color=RED, lw=.9, ls="-.")
axD.text(2.6, 5.6, "4% occupancy", fontsize=5.9, color=RED,
         bbox=dict(facecolor="white", edgecolor="none", pad=0.8))
axD.axvline(1.2, color=RED, lw=.8, ls=":")
axD.text(1.5, 0.012, "Kd ≈ 1 nM would be needed\neven for 4%: tighter than\ntestosterone, not attainable\nfor a same-charge anion",
         fontsize=5.5, color=RED, ha="left", va="bottom")
axD.set_xscale("log"); axD.set_yscale("log")
axD.set_xlabel("Kd for PFAS–SHBG (nM)", fontsize=7)
axD.set_ylabel("SHBG sites occupied (%)", fontsize=7)
axD.set_xlim(1, 1e4); axD.set_ylim(1e-4, 1e2)
axD.legend(fontsize=5.9, frameon=False, loc="lower left", handlelength=1.8,
           borderaxespad=.25, labelspacing=.28)
axD.set_title("d  Why the population assay sees nothing", loc="left",
              fontsize=8.5, fontweight="bold")

fig.savefig("../figures/fig14_pfas_shbg_structure.png", dpi=300, bbox_inches="tight")
fig.savefig("../figures/fig14_pfas_shbg_structure.pdf", bbox_inches="tight")
print("saved")
