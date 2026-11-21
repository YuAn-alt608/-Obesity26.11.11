import numpy as np

VDW = {'C':1.70,'F':1.47,'O':1.52,'S':1.80,'H':1.20,'N':1.55}

def read_sdf(fn):
    lines = open(fn).read().splitlines()
    na = int(lines[3][:3])
    at = []
    for L in lines[4:4+na]:
        at.append((L[31:34].strip(),
                   np.array([float(L[0:10]), float(L[10:20]), float(L[20:30])])))
    return at

def read_dht(fn):
    out = []
    for L in open(fn):
        if L.startswith("HETATM") and L[17:20].strip() == "DHT":
            out.append((L[76:78].strip(),
                        np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])])))
    return out

def vdw_volume(atoms, spacing=0.20):
    """网格法：计数被任一原子 vdW 球覆盖的格点"""
    R = np.array([VDW.get(e, 1.7) for e, _ in atoms])
    X = np.array([c for _, c in atoms])
    lo = (X - R[:, None]).min(0) - 1.0
    hi = (X + R[:, None]).max(0) + 1.0
    g = [np.arange(lo[i], hi[i], spacing) for i in range(3)]
    G = np.stack(np.meshgrid(*g, indexing="ij"), -1).reshape(-1, 3)
    inside = np.zeros(len(G), bool)
    for k in range(len(X)):
        d2 = ((G - X[k])**2).sum(1)
        inside |= d2 <= R[k]**2
    return inside.sum() * spacing**3

def shape(atoms):
    """重原子（排除 H）的惯性主轴与尺寸"""
    heavy = [(e, c) for e, c in atoms if e != "H"]
    X = np.array([c for _, c in heavy])
    c0 = X.mean(0)
    Xc = X - c0
    w, V = np.linalg.eigh(Xc.T @ Xc)
    order = np.argsort(w)[::-1]
    V = V[:, order]
    Y = Xc @ V
    ext = Y.max(0) - Y.min(0)          # 沿三主轴的跨度
    rmax = np.linalg.norm(Xc, axis=1).max()
    dmax = np.linalg.norm(X[:, None] - X[None], axis=-1).max()
    r = np.array([VDW.get(e, 1.7) for e, _ in heavy])
    return dict(n_heavy=len(heavy), centroid=c0, axes=V, ext=ext,
                rmax=rmax, dmax=dmax, volume=vdw_volume(heavy),
                vdw_vol_sum=(4/3*np.pi*r**3).sum(),
                mean_r=r.mean())

def axial_profile(atoms, nbins=22):
    """沿最长主轴的截面半径剖面（重原子 + vdW）"""
    heavy = [(e, c) for e, c in atoms if e != "H"]
    s = shape(atoms)
    X = np.array([c for _, c in heavy]) - s["centroid"]
    Y = X @ s["axes"]
    z = Y[:, 0]; rho = np.hypot(Y[:, 1], Y[:, 2])
    R = np.array([VDW.get(e, 1.7) for e, _ in heavy])
    edges = np.linspace(z.min() - 0.5, z.max() + 0.5, nbins + 1)
    prof = []
    for i in range(nbins):
        m = (z >= edges[i]) & (z < edges[i+1])
        if m.any():
            prof.append(((edges[i]+edges[i+1])/2, (rho[m] + R[m]).max()))
        else:
            prof.append(((edges[i]+edges[i+1])/2, np.nan))
    return np.array(prof)

print(f"{'分子':8s}{'重原子':>5s}{'最长跨度':>9s}{'三轴跨度 (Å)':>22s}{'vdW体积':>10s}")
print("-"*62)
res = {}
for nm, fn in [("DHT","1D2S.pdb")] + [(n, f"pc_{n}.sdf") for n in ["PFOA","PFOS","PFHxS","PFBS"]]:
    at = read_dht(fn) if fn.endswith(".pdb") else read_sdf(fn)
    s = shape(at); res[nm] = s
    e = s["ext"]
    print(f"{nm:8s}{s['n_heavy']:>5d}{s['dmax']:>9.2f}"
          f"{e[0]:>8.2f}{e[1]:>7.2f}{e[2]:>7.2f}{s['volume']:>10.0f}")

print("\n=== 与 DHT 的相对形状（体积比 / 截面积比）===")
d = res["DHT"]
for nm in ["PFOA","PFOS","PFHxS","PFBS"]:
    s = res[nm]
    # 中轴截面：第二、三轴跨度的乘积 作为"扁/圆"的代理
    a_d = d["ext"][1]*d["ext"][2]; a_p = s["ext"][1]*s["ext"][2]
    print(f"{nm:8s} 体积 {s['volume']/d['volume']:.2f}×   "
          f"长轴 {s['ext'][0]/d['ext'][0]:.2f}×   "
          f"截面(轴2×轴3) {np.sqrt(a_p/a_d):.2f}×")
np.save("shapecache.npy", np.array([1]))
