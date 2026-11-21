# -*- coding: utf-8 -*-
"""
V10b  第三层：Python 独立复现（statsmodels CR1 聚类稳健标准误）
---------------------------------------------------------------
输入与 R 端完全相同（verification/_py_*.csv），独立实现：
  OLS + cov_type="cluster", use_correction=True  ->  CR1 = G/(G-1) * (n-1)/(n-k)
比较 beta、聚类稳健 SE、P 值是否在 1e-10 量级内一致。
产出：verification/python_crosscheck.csv
"""
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm

BASE = os.environ.get("AGE22_ROOT", os.getcwd())
V = os.path.join(BASE, "verification")

COV = ["age", "female", "pir", "bmi_p", "log_cot", "drink2",
       "pa_any", "comorbid", "tc", "hdl", "egfr", "yr"]
SPECS = [("androgen", "y_tt", "Total testosterone (log), per 5%"),
         ("liver", "y_fli", "FLI (log), per 5%"),
         ("inflam", "y_wbc", "WBC (log), per 5%")]

rows = []
for tag, y, label in SPECS:
    d = pd.read_csv(os.path.join(V, f"_py_{tag}.csv"))
    cols = [y, "dw5"] + COV + ["strata", "psu"]
    d = d[cols].dropna()
    X = d[["dw5"] + COV].astype(float)
    X = sm.add_constant(X)
    # 剔除零方差列（与 R 一致）；但必须保留截距 const
    # （V16 审计：初版此处把 const 一并删掉，R 端犯同一错，故"一致"不可作为正确性证据）
    keep = [c for c in X.columns if X[c].std(ddof=1) > 1e-12 or c == "const"]
    X = X[keep]
    grp = (d["strata"].astype(float).astype(int).astype(str) + "_" +
           d["psu"].astype(float).astype(int).astype(str))
    m = sm.OLS(d[y].astype(float), X).fit(cov_type="cluster",
                                          cov_kwds={"groups": grp, "use_correction": True})
    rows.append(dict(dataset=tag, term=label, N=int(m.nobs), n_cluster=int(grp.nunique()),
                     beta=m.params["dw5"], se_cr1=m.bse["dw5"], p_cr1=m.pvalues["dw5"],
                     pct_per5=100 * (np.exp(m.params["dw5"]) - 1)))

out = pd.DataFrame(rows)
out.to_csv(os.path.join(V, "python_crosscheck.csv"), index=False)

r = pd.read_csv(os.path.join(V, "python_crosscheck_r.csv"))
mrg = out.merge(r, on="dataset", suffixes=("_py", "_r"))
mrg["d_beta"] = (mrg.beta_py - mrg.beta_r).abs()
mrg["d_se"] = (mrg.se_cr1_py - mrg.se_cr1_r).abs()
mrg["d_pct"] = (mrg.pct_per5_py - mrg.pct_per5_r).abs()
mrg[["dataset", "N_py", "N_r", "beta_py", "beta_r", "d_beta",
     "se_cr1_py", "se_cr1_r", "d_se", "p_cr1_py", "p_cr1_r", "d_pct"]].to_csv(
    os.path.join(V, "python_vs_r.csv"), index=False)

print("=" * 78)
print("V10  Python vs R 交叉复算（同一份 CSV，两套独立代码）")
print("=" * 78)
print(mrg[["dataset", "N_py", "N_r", "beta_py", "beta_r", "d_beta",
           "se_cr1_py", "se_cr1_r", "d_se"]].to_string(index=False, float_format=lambda x: f"{x:.3e}"))
print("\n最大 |Δbeta| =", mrg.d_beta.max(), " 最大 |ΔSE| =", mrg.d_se.max(),
      " 最大 |Δ百分点| =", mrg.d_pct.max())
print("两套代码一致性判定：", "完全一致（<1e-9）" if max(mrg.d_beta.max(), mrg.d_se.max()) < 1e-9
      else "存在差异，需排查")
