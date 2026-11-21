# -*- coding: utf-8 -*-
"""
V01  第一层：原始文件指纹清单
--------------------------------
目的：排除网络传输截断、镜像替换、文件被误编辑。
产出：verification/sha256_manifest.csv  （可直接作为 Data Availability 附件）

覆盖
  data/raw/*.xpt                  NHANES 传输文件
  data/raw/mort/*.dat             NCHS 公共链接死亡随访
  data/raw/geo/GSE295708_RAW.tar  GEO 原始归档
  data/raw/geo/extracted/*        （仅登记大小 + gzip 完整性，不重复哈希 2.5 GB）
"""
import hashlib, os, sys, csv, gzip, time

BASE = r"D:\NewBeginning\BIOinformation\Age2.2"
RAW = os.path.join(BASE, "data", "raw")
OUT = os.path.join(BASE, "verification")
os.makedirs(OUT, exist_ok=True)


def digest(path, algo="sha256", buf=8 << 20):
    h = hashlib.new(algo)
    n = 0
    with open(path, "rb") as f:
        while True:
            b = f.read(buf)
            if not b:
                break
            h.update(b)
            n += len(b)
    return h.hexdigest(), n


def target_files():
    for f in sorted(os.listdir(RAW)):
        p = os.path.join(RAW, f)
        if os.path.isfile(p) and f.lower().endswith(".xpt"):
            yield p, "NHANES_XPT"
    mort = os.path.join(RAW, "mort")
    if os.path.isdir(mort):
        for f in sorted(os.listdir(mort)):
            p = os.path.join(mort, f)
            if os.path.isfile(p):
                yield p, "NCHS_MORT"
    tar = os.path.join(RAW, "geo", "GSE295708_RAW.tar")
    if os.path.isfile(tar):
        yield tar, "GEO_ARCHIVE"


def main():
    rows = []
    t0 = time.time()
    for p, kind in target_files():
        sha, size = digest(p)
        md5, _ = digest(p, "md5")
        st = os.stat(p)
        rows.append(dict(kind=kind,
                         file=os.path.relpath(p, BASE).replace("\\", "/"),
                         bytes=size,
                         sha256=sha, md5=md5,
                         mtime=time.strftime("%Y-%m-%dT%H:%M:%S", time.localtime(st.st_mtime))))
    with open(os.path.join(OUT, "sha256_manifest.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["kind", "file", "bytes", "sha256", "md5", "mtime"])
        w.writeheader()
        w.writerows(rows)

    # GEO 解压后的 54 个文件：只做 gzip 完整性 + 大小登记
    ex = os.path.join(RAW, "geo", "extracted")
    ok = bad = 0
    ex_rows = []
    if os.path.isdir(ex):
        for f in sorted(os.listdir(ex)):
            p = os.path.join(ex, f)
            if not os.path.isfile(p):
                continue
            good = True
            if f.endswith(".gz"):
                try:
                    with gzip.open(p, "rb") as g:
                        while g.read(4 << 20):
                            pass
                except Exception:
                    good = False
            ok += good
            bad += (not good)
            ex_rows.append((f, os.path.getsize(p), "ok" if good else "CORRUPT"))
    with open(os.path.join(OUT, "geo_extracted_integrity.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["file", "bytes", "gzip_status"])
        w.writerows(ex_rows)

    xpt = [r for r in rows if r["kind"] == "NHANES_XPT"]
    mort = [r for r in rows if r["kind"] == "NCHS_MORT"]
    print("=" * 70)
    print("V01 原始文件指纹清单")
    print("=" * 70)
    print(f"NHANES XPT        : {len(xpt):>4d} 个，共 {sum(r['bytes'] for r in xpt)/1e6:,.1f} MB")
    print(f"NCHS 死亡文件     : {len(mort):>4d} 个，共 {sum(r['bytes'] for r in mort)/1e6:,.1f} MB")
    print(f"GEO 归档          : {len(rows)-len(xpt)-len(mort)} 个")
    print(f"GEO 解压文件      : {ok} 个 gzip 完整性通过，{bad} 个损坏")
    print(f"清单 -> verification/sha256_manifest.csv（{len(rows)} 行）")
    print(f"耗时 {time.time()-t0:.1f}s")
    if bad:
        print("!! 存在损坏文件，见 geo_extracted_integrity.csv")
        sys.exit(1)


if __name__ == "__main__":
    main()
