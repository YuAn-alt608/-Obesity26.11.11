# -*- coding: utf-8 -*-
"""
V01b 第一层（加强）：与 CDC 官方服务器逐字节比对
--------------------------------------------------
不满足于"本地清单自证"：把 CDC 现行发布的同一文件下载下来，与本地副本
逐字节比较 SHA256。相同 => 排除传输截断、镜像替换、人为改动；
不同 => 立即锁定该文件（并给出官方与本地大小差）。

URL 规则（2024 年后 CDC 改版）：
  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{YEAR}/DataFiles/{FILE}.XPT
周期后缀 -> 年份：C=2003 D=2005 F=2009 G=2011 H=2013 I=2015 J=2017 L=2021
"""
import hashlib, os, csv, sys, time, urllib.request

BASE = r"D:\NewBeginning\BIOinformation\Age2.2"
RAW = os.path.join(BASE, "data", "raw")
OUT = os.path.join(BASE, "verification")
TMP = os.path.join(OUT, "_tmp")
os.makedirs(TMP, exist_ok=True)

YEAR = dict(C="2003", D="2005", F="2009", G="2011", H="2013", I="2015", J="2017", L="2021")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(8 << 20), b""):
            h.update(b)
    return h.hexdigest()


def stream_sha(url):
    """不落盘：边下载边哈希"""
    h = hashlib.sha256()
    n = 0
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=120) as r:
        for b in iter(lambda: r.read(1 << 20), b""):
            h.update(b)
            n += len(b)
    return h.hexdigest(), n


def main():
    files = sorted(f for f in os.listdir(RAW) if f.lower().endswith(".xpt"))
    rows = []
    t0 = time.time()
    for i, f in enumerate(files, 1):
        sfx = f.rsplit("_", 1)[-1].replace(".xpt", "").upper()
        yr = YEAR.get(sfx)
        if not yr:
            rows.append(dict(file=f, status="NO_YEAR_MAP", official_url="", match=""))
            continue
        url = f"https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/{yr}/DataFiles/{f}"
        try:
            h_off, n_off = stream_sha(url)
        except Exception as e:
            rows.append(dict(file=f, status=f"HTTP_FAIL:{type(e).__name__}", official_url=url, match=""))
            continue
        h_local = sha256(os.path.join(RAW, f))
        n_local = os.path.getsize(os.path.join(RAW, f))
        rows.append(dict(file=f,
                         status="OK",
                         official_url=url,
                         local_bytes=n_local,
                         official_bytes=n_off,
                         local_sha256=h_local, official_sha256=h_off,
                         match="IDENTICAL" if (h_local == h_off and n_local == n_off) else "DIFFER"))
        if i % 25 == 0:
            print(f"  {i}/{len(files)} ...", flush=True)

    keys = ["file", "status", "match", "local_bytes", "official_bytes", "local_sha256", "official_sha256", "official_url"]
    for r in rows:
        for k in keys:
            r.setdefault(k, "")
    with open(os.path.join(OUT, "official_bytewise_comparison.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    ident = [r for r in rows if r["match"] == "IDENTICAL"]
    diff = [r for r in rows if r["match"] == "DIFFER"]
    fail = [r for r in rows if r["status"] != "OK"]
    print("=" * 66)
    print(f"逐字节比对 {len(rows)} 个文件：一致 {len(ident)}，不一致 {len(diff)}，无法取回 {len(fail)}")
    for r in diff:
        print("  DIFFER:", r["file"], r["local_bytes"], "vs", r["official_bytes"])
    for r in fail:
        print("  FAIL  :", r["file"], r["status"])
    print(f"耗时 {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
