import numpy as np, csv
from geom import read_sdf, read_dht, shape, VDW

MW = {"PFOA":414.07, "PFOS":500.13, "PFHxS":400.11, "PFBS":300.10, "DHT":290.44}
CONC_NGML = {"PFOA":2.07, "PFOS":6.80, "PFHxS":1.30, "PFBS":np.nan}

# ---------- 表1：形状 ----------
rows = []
dht = read_dht("1D2S.pdb"); sD = shape(dht)
ref_vol = sD["volume"]
for nm, fn in [("DHT","1D2S.pdb")] + [(n, f"pc_{n}.sdf") for n in ["PFOA","PFOS","PFHxS","PFBS"]]:
    at = read_dht(fn) if fn.endswith(".pdb") else read_sdf(fn)
    s = shape(at); e = s["ext"]
    X = np.array([c for el, c in at if el != "H"]) - s["centroid"]
    Y = X @ s["axes"]
    r = np.array([VDW.get(el, 1.7) for el, _ in at if el != "H"])
    rho = np.hypot(Y[:, 1], Y[:, 2])
    rows.append(dict(
        molecule=nm, n_heavy=s["n_heavy"],
        max_span_A=round(s["dmax"], 2),
        axis1_A=round(e[0], 2), axis2_A=round(e[1], 2), axis3_A=round(e[2], 2),
        thinnest_A=round(e[2], 2),
        max_radius_A=round((rho + r).max(), 2),
        vdw_volume_A3=round(s["volume"], 0),
        vol_ratio_vs_DHT=round(s["volume"] / ref_vol, 3),
        thinnest_ratio_vs_DHT=round(e[2] / sD["ext"][2], 2),
        comment=("天然配体（参考）" if nm == "DHT" else
                 "C8 全氟羧酸" if nm == "PFOA" else
                 "C8 全氟磺酸" if nm == "PFOS" else
                 "C6 全氟磺酸" if nm == "PFHxS" else "C4 全氟磺酸（短链对照）")))
with open("t1_ligand_shape.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
print("t1_ligand_shape.csv  →", len(rows), "行")

# ---------- 表2：口袋剖面 ----------
prof = np.load("profile.npz", allow_pickle=True)
pr, dpr = prof["prof"], prof["dht_prof"]
with open("t2_pocket_profile.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["z_along_axis_A", "pocket_available_radius_A", "DHT_radius_A", "clearance_A"])
    for (z, wall), rl in zip(pr, dpr):
        w.writerow([round(z, 2), round(wall, 2) if not np.isnan(wall) else "",
                    round(rl, 2) if not np.isnan(rl) else "",
                    round(wall - rl, 2) if not (np.isnan(wall) or np.isnan(rl)) else ""])
print("t2_pocket_profile.csv →", len(pr), "行；最窄可用半径",
      f"{np.nanmin(pr[:,1]):.2f} A")

# ---------- 表3：电荷环境 ----------
charge = [
    ("SHBG", "steroid pocket (1D2S, DHT-bound)", 0, 1, -1,
     "ASP65 (bidentate, OD1/OD2 at 2.73/2.72 A from DHT 17-OH)",
     "为中性 17beta-OH 设计：负电羧基当氢键受体；对阴离子头基同号相斥"),
    ("HSA", "fatty-acid sites FA1-FA7 (1E7I, stearate)", 16, 2, 14,
     "ARG117/209/257/348/410/485, LYS525",
     "为阴离子脂肪酸羧酸根设计：正电胍基/氨基当锚点"),
]
with open("t3_pocket_charge.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["protein", "pocket", "n_positive", "n_negative", "n_charged",
                "net_charge", "key_residues", "design_logic"])
    for prot, pk, npos, nneg, net, key, logic in charge:
        w.writerow([prot, pk, npos, nneg, npos + nneg, net, key, logic])
print("t3_pocket_charge.csv → 2 行（SHBG 净 -1 / HSA 净 +14）")

# ---------- 表4：占据率敏感性 ----------
SHBG_NM, ALB_UM = 46.0, 632.0
def occupancy(L_tot_nM, f_free, Kd_nM):
    Lf = L_tot_nM * f_free
    return 100.0 * Lf / (Kd_nM + Lf)
rows = []
for an, c_ng in [("PFOA", 2.07), ("PFOS", 6.80)]:
    for tag, pct in [("median", 50), ("P95", 95)]:
        c = c_ng if pct == 50 else {"PFOA":6.82, "PFOS":34.10}[an]
        L = c / MW[an] * 1000
        for f_free in (0.01, 0.05):
            for Kd in (3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0):
                rows.append(dict(
                    analyte=an, level=tag, conc_ng_ml=c, conc_nM=round(L, 2),
                    free_fraction=f_free, free_conc_nM=round(L*f_free, 3),
                    Kd_nM=Kd, occupancy_pct=round(occupancy(L, f_free, Kd), 3),
                    enough_for_4pct=("YES" if occupancy(L, f_free, Kd) >= 4 else "no")))
with open("t4_occupancy.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]))
    w.writeheader(); w.writerows(rows)
ok = [r for r in rows if r["enough_for_4pct"] == "YES"]
print(f"t4_occupancy.csv → {len(rows)} 行；达到 4% 占据的组合 {len(ok)} 个")
print("\n=== 关键摩尔比 ===")
print(f"PFOA 5.00 nM : SHBG {SHBG_NM} nM = {5.00/SHBG_NM:.3f}  (即使 100% 占据也只有 {5.0/SHBG_NM*100:.0f}%)")
print(f"白蛋白 {ALB_UM*1000:.0f} nM : PFOA 5.00 nM = {ALB_UM*1000/5.0:,.0f} 倍")
