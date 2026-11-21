import numpy as np, collections

def read_pdb(fn):
    prot, lig, het = [], [], []
    for L in open(fn):
        if L.startswith("ATOM"):
            prot.append((L[21], int(L[22:26]), L[12:16].strip(),
                         L[17:20].strip(),
                         np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])]),
                         L[76:78].strip()))
        elif L.startswith("HETATM"):
            res = L[17:20].strip()
            rec = (res, int(L[22:26]), L[12:16].strip(),
                   np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])]),
                   L[76:78].strip())
            (lig if res == "DHT" else het).append(rec)
    return prot, lig, het

prot, lig, het = read_pdb("1D2S.pdb")
P = np.array([a[4] for a in prot])
L = np.array([a[3] for a in lig])
print(f"蛋白原子 {len(P)}，DHT 原子 {len(L)}（重原子 {sum(1 for a in lig if a[4]!='H')}）")

# DHT 质心与尺寸
c = L.mean(0)
d = np.linalg.norm(L - c, axis=1)
print(f"DHT 质心 {np.round(c,2)}，最大半径 {d.max():.2f} Å")
print(f"DHT 最大原子间距 {np.linalg.norm(L[:,None]-L[None],axis=-1).max():.2f} Å")

# 口袋残基：任一原子距 DHT 任一原子 < 4.5 Å
dmin = np.linalg.norm(P[:,None]-L[None], axis=-1).min(1)
sel = dmin < 4.5
res = collections.defaultdict(lambda: 9e9)
for i, a in enumerate(prot):
    if sel[i]:
        k = (a[1], a[3])
        res[k] = min(res[k], dmin[i])
print(f"\n4.5 Å 内残基 {len(res)} 个：")
for (n, nm), dd in sorted(res.items(), key=lambda x: x[0][0]):
    print(f"  {nm}{n:<4d} 最近距离 {dd:.2f} Å")
