# -*- coding: utf-8 -*-
"""
V06  第二层：Vermeulen 游离睾酮方程的独立复算
--------------------------------------------------
R 端 (07_build_horm.R) 用的是 Vermeulen 1999 的**闭式二次方程解**。
这里用完全不同的数值方法——**二分法求根**——在质量作用平衡方程上求解：
    T_f * (1 + Kt*S_free + Ka*C_alb) = T_total,  S_free = S_total/(1+Kt*T_f)
不依赖任何代数变形，因此可以检验 R 端是否把系数/单位/公式抄错。

常数（Vermeulen 1999, J Clin Endocrinol Metab）：
    Kt(T-SHBG) = 1.0e9  L/mol
    Ka(T-alb)  = 3.6e4  L/mol
    albumin MW = 66430 g/mol
产出：verification/vermeulen_check.csv、vermeulen_sample.csv
"""
import csv, os, math

BASE = r"D:\NewBeginning\BIOinformation\Age2.2"
V = os.path.join(BASE, "verification")
KT, KA, MW_ALB = 1.0e9, 3.6e4, 66430.0


def free_testosterone_bisect(tt_ngdl, shbg_nmol, alb_gL, tol=1e-14, itmax=400):
    """二分法：在 [0, Ttot] 上解质量作用方程，返回游离睾酮 nmol/L"""
    ## 全部换算到 SI 基本单位（mol/L），否则 Kt=1e9 L/mol 的量纲不成立
    Ttot = tt_ngdl * 0.034671e-9       # ng/dL -> nmol/L -> mol/L
    Stot = shbg_nmol * 1e-9            # nmol/L -> mol/L
    Calb = alb_gL / MW_ALB             # g/L -> mol/L

    def f(Tf):
        Sfree = Stot / (1.0 + KT * Tf)
        bound = KT * Tf * Sfree + KA * Tf * Calb
        return Tf + bound - Ttot       # 守恒：游离 + 结合 = 总

    lo, hi = 1e-18, Ttot
    if f(hi) < 0:
        return hi
    for _ in range(itmax):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
        if hi - lo < tol * hi:
            break
    return 0.5 * (lo + hi) * 1e9        # mol/L -> nmol/L


def main():
    rows = list(csv.DictReader(open(os.path.join(V, "_horm_full.csv"), encoding="utf-8")))
    out = []
    for r in rows:
        tt, sh, al = float(r["tt_ngdl"]), float(r["shbg"]), float(r["alb_gL"])
        ft_py = free_testosterone_bisect(tt, sh, al)
        ft_r = float(r["ft_nmol"])
        out.append(dict(SEQN=r["SEQN"], cycle=r["cycle"], sex=r["sex"], age=r["age"],
                        tt_ngdl=tt, shbg=sh, alb_gL=al,
                        ft_r_nmol=ft_r, ft_py_nmol=ft_py,
                        abs_diff=abs(ft_py - ft_r),
                        rel_diff_pct=100 * abs(ft_py - ft_r) / ft_r if ft_r else 0,
                        pctfree_r=100 * ft_r / (tt * 0.034671),
                        pctfree_py=100 * ft_py / (tt * 0.034671)))

    keys = list(out[0].keys())
    with open(os.path.join(V, "vermeulen_check.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(out)

    import statistics as st
    absmax = max(o["abs_diff"] for o in out)
    relmax = max(o["rel_diff_pct"] for o in out)
    rel_med = st.median(o["rel_diff_pct"] for o in out)
    pf = sorted(o["pctfree_py"] for o in out)

    print("=" * 66)
    print(f"V06 Vermeulen 独立复算：{len(out):,} 例（覆盖全部有 TT+SHBG+Alb 的样本）")
    print(f"  最大绝对偏差 : {absmax:.3e} nmol/L")
    print(f"  最大相对偏差 : {relmax:.3e} %")
    print(f"  中位相对偏差 : {rel_med:.3e} %")
    print(f"  %游离睾酮中位 : {st.median(pf):.2f}%  （文献人群值 1.5–2.0%）")
    print(f"  10-90 百分位  : {pf[int(.1*len(pf))]:.2f}% – {pf[int(.9*len(pf))]:.2f}%")

    # 教科书校验点：TT 20 nmol/L, SHBG 30, Alb 43 -> 约 1.5–1.6% 游离
    ft = free_testosterone_bisect(20 / 0.034671, 30, 43)
    print(f"  校验点 TT=20 nmol/L, SHBG=30, Alb=43 -> FT={ft:.3f} nmol/L, "
          f"游离比例={100*ft/20:.2f}%")

    # 抽样 12 例另存
    samp = [o for o in out if o["SEQN"] in
            {r["SEQN"] for r in csv.DictReader(open(os.path.join(V, "_horm_sample.csv"), encoding="utf-8"))}]
    with open(os.path.join(V, "vermeulen_sample.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(samp)
    print(f"  抽样 {len(samp)} 例明细 -> verification/vermeulen_sample.csv")
    print("=" * 66)


if __name__ == "__main__":
    main()
