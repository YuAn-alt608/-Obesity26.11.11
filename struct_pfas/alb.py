import numpy as np, collections
P, L = [], []
for l in open("1E7I.pdb"):
    if l.startswith("ATOM"):
        P.append((int(l[22:26]), l[17:20].strip(), l[12:16].strip(),
                  np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), l[76:78].strip()))
    elif l.startswith("HETATM") and l[17:20].strip() == "STE":
        L.append((int(l[22:26]), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]), l[76:78].strip()))
print(f"HSA 蛋白原子 {len(P)}，硬脂酸(STE) 重原子 {len(L)}")

# 每个 STE 分子单独分析
mols = sorted({m for m, _, _ in L})
print("STE 分子编号:", mols)

CHG = {"ASP":-1, "GLU":-1, "LYS":+1, "ARG":+1}
for mid in mols:
    C = np.array([c for m, c, e in L if m == mid])
    if len(C) < 15: continue
    Pp = np.array([a[3] for a in P])
    dmin = np.linalg.norm(Pp[:, None] - C[None], axis=-1).min(1)
    near = [(P[i][0], P[i][1], P[i][2], dmin[i]) for i in range(len(P)) if dmin[i] < 4.5]
    # 侧链带电原子（排除主链 N/O）
    charged = [(r, n, nm, d) for r, nm, n, d in near
               if nm in CHG and n in ("OD1","OD2","OE1","OE2","NZ","NE","NH1","NH2")]
    pos = sum(1 for _, _, nm, _ in charged if CHG[nm] > 0)
    neg = sum(1 for _, _, nm, _ in charged if CHG[nm] < 0)
    if charged:
        print(f"\nSTE {mid}: 4.5 Å 内带电侧链 {len(charged)} 个（+{pos} / -{neg}）")
        for r, n, nm, d in sorted(charged, key=lambda t: t[3])[:8]:
            print(f"   {nm}{r:<4d} {n:4s} {d:.2f} Å  {CHG[nm]:+.0f}")

# 汇总：所有 STE 位点合计
allc = []
for mid in mols:
    C = np.array([c for m, c, e in L if m == mid])
    if len(C) < 15: continue
    Pp = np.array([a[3] for a in P])
    dmin = np.linalg.norm(Pp[:, None] - C[None], axis=-1).min(1)
    for i in range(len(P)):
        if dmin[i] < 4.5 and P[i][1] in CHG and P[i][2] in ("OD1","OD2","OE1","OE2","NZ","NE","NH1","NH2"):
            allc.append(CHG[P[i][1]])
print(f"\n=== 全部脂肪酸位点合计：正电 {sum(1 for x in allc if x>0)} 个 / 负电 {sum(1 for x in allc if x<0)} 个 ===")
print(f"    净电荷 = {sum(allc):+.0f}")
