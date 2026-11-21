import numpy as np
from geom import read_dht

prot = []
for L in open("1D2S.pdb"):
    if L.startswith("ATOM"):
        prot.append(dict(name=L[12:16].strip(), resi=int(L[22:26]),
                         resn=L[17:20].strip(), el=L[76:78].strip(),
                         xyz=np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])])))

dht = read_dht("1D2S.pdb")

# DHT 的极性原子：3-酮 O、17-OH O
print("=== DHT 的极性原子 ===")
for el, c in dht:
    if el in ("O", "N"):
        print(f"  {el}  坐标 {np.round(c,2)}")

# 口袋内所有可能形成氢键/带电的侧链原子
POLAR = {
    "ASP": ["OD1","OD2"], "GLU": ["OE1","OE2"],
    "ASN": ["OD1","ND2"], "GLN": ["OE1","NE2"],
    "SER": ["OG"], "THR": ["OG1"], "TYR": ["OH"],
    "ARG": ["NE","NH1","NH2"], "LYS": ["NZ"], "HIS": ["ND1","NE2"],
    "TRP": ["NE1"],
}
print("\n=== 口袋内极性/带电侧链原子（距任一 DHT 极性原子 < 5 Å）===")
print(f"{'残基':14s}{'原子':6s}{'到最近 DHT-O':>13s}{'电荷':>6s}")
rows = []
for a in prot:
    if a["resn"] in POLAR and a["name"] in POLAR[a["resn"]]:
        d = min(np.linalg.norm(a["xyz"] - c) for el, c in dht if el == "O")
        if d < 5.0:
            q = {"ASP": -1, "GLU": -1, "LYS": +1, "ARG": +1,
                 "HIS": +0.1}.get(a["resn"], 0.0)
            rows.append((a["resn"], a["resi"], a["name"], d, q))
            print(f"{a['resn']+str(a['resi']):14s}{a['name']:6s}{d:13.2f}{q:6.1f}")

# 关键：Asp65 的羧基氧到 DHT 17-OH 氧的距离
print("\n=== 关键接触（≤3.5 Å）===")
for r, i, n, d, q in sorted(rows, key=lambda t: t[3]):
    if d <= 3.5:
        print(f"  {r}{i} {n}  —  {d:.2f} Å   电荷 {q:+.0f}")

# 库仑粗算：若把 DHT 17-OH 换成 PFAS 的羧酸根/磺酸根（-1）
print("\n=== 若阴离子头基占据 DHT 17-OH 位置：与口袋侧链的静电作用 ===")
print("  两个 -1 电荷相距 d 时，ΔG ≈ 332·q1·q2/(ε·d) kcal/mol")
for d in (2.7, 3.5, 4.5):
    for eps in (10, 20, 40):
        dG = 332.0 * ((-1) * (-1)) / (eps * d)
        print(f"    d={d:.1f} Å, ε={eps:2d} → ΔG = +{dG:.1f} kcal/mol")
