import numpy as np
prot = []
for L in open("1D2S.pdb"):
    if L.startswith("ATOM"):
        prot.append((int(L[22:26]), L[17:20].strip(), L[12:16].strip(),
                     np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])]), L[76:78].strip()))
P = np.array([a[3] for a in prot])
Os = []
for L in open("1D2S.pdb"):
    if L.startswith("HETATM") and L[17:20].strip() == "DHT" and L[76:78].strip() == "O":
        Os.append(np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])]))

for k, O in enumerate(Os, 1):
    d = np.linalg.norm(P - O, axis=1)
    idx = np.argsort(d)[:4]
    print(f"\nDHT 氧原子 {k}  ({np.round(O,2)})  4 Å 内最近残基：")
    seen = set()
    for i in idx:
        r, nm, at, _, el = prot[i]
        if (r, nm) in seen: continue
        seen.add((r, nm))
        print(f"   {nm}{r:<4d} {at:4s} {d[i]:.2f} Å")

# DHT 碳骨架端点：找最远的两个碳（判断 A 环 vs D 环方位）
Cs = [np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])])
      for L in open("1D2S.pdb") if L.startswith("HETATM") and L[17:20].strip()=="DHT" and L[76:78].strip()=="C"]
from itertools import combinations
best = max(combinations(range(len(Cs)), 2), key=lambda ij: np.linalg.norm(Cs[ij[0]]-Cs[ij[1]]))
print(f"\nDHT 最长碳-碳距离 {np.linalg.norm(Cs[best[0]]-Cs[best[1]]):.2f} Å")
