import numpy as np
from geom import read_sdf, read_dht, shape, VDW

prot = []
for L in open("1D2S.pdb"):
    if L.startswith("ATOM"):
        prot.append((L[12:16].strip(), int(L[22:26]), L[17:20].strip(),
                     np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])]),
                     L[76:78].strip()))
P = np.array([a[3] for a in prot])

dht = read_dht("1D2S.pdb")
sD = shape(dht)
cen, A = sD["centroid"], sD["axes"]          # A[:,0]=长轴
Lig = np.array([c for _, c in dht]) - cen
Loc = Lig @ A                                 # DHT 自身主轴坐标

def pocket_profile(center, axes, span, nb=20, rmax=9.0):
    """沿 axes[:,0] 切片，返回 (z_center, 该切片内蛋白原子到轴的最小垂直距离)"""
    Pp = (P - center) @ axes
    zmin, zmax = span
    edges = np.linspace(zmin, zmax, nb + 1)
    out = []
    for i in range(nb):
        m = (Pp[:, 0] >= edges[i]) & (Pp[:, 0] < edges[i + 1])
        if not m.any():
            out.append(((edges[i] + edges[i + 1]) / 2, np.nan)); continue
        rho = np.hypot(Pp[m, 1], Pp[m, 2])
        wall = rho[rho < rmax]
        out.append(((edges[i] + edges[i + 1]) / 2, wall.min() if len(wall) else np.nan))
    return np.array(out)

prof = pocket_profile(cen, A, (Loc[:, 0].min() - 0.5, Loc[:, 0].max() + 0.5), nb=20)

# DHT 在各切片的实际半径（含 vdW）
R_dht = np.array([VDW.get(e, 1.7) for e, _ in dht])
rho_dht = np.hypot(Loc[:, 1], Loc[:, 2])
edges = np.linspace(Loc[:, 0].min() - 0.5, Loc[:, 0].max() + 0.5, 21)
dht_prof = []
for i in range(20):
    m = (Loc[:, 0] >= edges[i]) & (Loc[:, 0] < edges[i + 1])
    dht_prof.append((rho_dht[m] + R_dht[m]).max() if m.any() else np.nan)
dht_prof = np.array(dht_prof)

print("沿 DHT 长轴的口袋半径剖面（Å）")
print(f"{'z (Å)':>8s}{'口袋壁 r':>10s}{'DHT r':>9s}{'余量':>8s}")
for (z, w), rl in zip(prof, dht_prof):
    print(f"{z:8.2f}{w:10.2f}{rl:9.2f}{w - rl:8.2f}")

narrow = np.nanmin(prof[:, 1])
print(f"\n最窄处 口袋可用半径 = {narrow:.2f} Å（对应可容纳的配体最大半径）")

# 各 PFAS 的最粗处半径（沿自身长轴的等效半径）
print("\n=== 各配体的'最粗半径'（决定能否通过窄口）===")
res = {}
for nm in ["PFOA", "PFOS", "PFHxS", "PFBS"]:
    at = read_sdf(f"pc_{nm}.sdf")
    s = shape(at)
    X = np.array([c for e, c in at if e != "H"]) - s["centroid"]
    Y = X @ s["axes"]
    r = np.array([VDW.get(e, 1.7) for e, _ in at if e != "H"])
    rho = np.hypot(Y[:, 1], Y[:, 2])
    res[nm] = (rho + r).max()
    print(f"{nm:6s} 最粗半径 {res[nm]:.2f} Å   三轴 {s['ext'][0]:.2f}×{s['ext'][1]:.2f}×{s['ext'][2]:.2f}")

print(f"\nDHT 最粗半径 {np.nanmax(dht_prof):.2f} Å")
np.savez("profile.npz", prof=prof, dht_prof=dht_prof,
         lig=np.array([[k, v] for k, v in res.items()], dtype=object))
