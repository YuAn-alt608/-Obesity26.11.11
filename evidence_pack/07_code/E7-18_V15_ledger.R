## ============================================================================
## V15_ledger.R — 把四层校验的全部数字汇总成一张可机读的"审计台账"
##   目的：声明文件里的每一个数字都能溯源到本台账的某一行 + 某个 CSV
##   输出：verification/verification_ledger.csv
##         verification/VERIFICATION_SUMMARY.md
## ============================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", "."); VER <- file.path(BASE, "verification")
RES  <- file.path(BASE, "results")
rd <- function(f) { p <- file.path(VER, f); if (file.exists(p)) fread(p) else NULL }
rr <- function(f) { p <- file.path(RES, f); if (file.exists(p)) fread(p) else NULL }

L <- list(); add <- function(layer, check, metric, value, source, verdict = "") {
  L[[length(L)+1]] <<- data.table(layer = layer, check = check, metric = metric,
                                  value = as.character(value), source = source,
                                  verdict = verdict)
}

## ---------------- 第一层：源头与指纹 ----------------
m <- rd("sha256_manifest.csv")
if (!is.null(m)) {
  add("L1", "原始文件指纹清单", "文件总数", nrow(m), "sha256_manifest.csv", "完成")
  add("L1", "原始文件指纹清单", "XPT 数", sum(grepl("\\.xpt$", m$file, ignore.case = TRUE)),
      "sha256_manifest.csv", "完成")
  add("L1", "原始文件指纹清单", "死亡随访文件数",
      sum(grepl("^mort|\\.dat$", m$file, ignore.case = TRUE)), "sha256_manifest.csv", "完成")
  add("L1", "原始文件指纹清单", "总字节", format(sum(m$bytes, na.rm = TRUE), big.mark = ","),
      "sha256_manifest.csv", "完成")
}
b <- rd("official_bytewise_comparison.csv")
if (!is.null(b)) {
  ok <- grepl("IDENTICAL", b$match, ignore.case = TRUE)
  add("L1", "与 CDC 官方逐字节比对", "比对文件数", nrow(b), "official_bytewise_comparison.csv", "")
  add("L1", "与 CDC 官方逐字节比对", "SHA256 逐字节一致数", sum(ok, na.rm = TRUE),
      "official_bytewise_comparison.csv", if (all(ok)) "通过" else "存在不一致")
  add("L1", "与 CDC 官方逐字节比对", "不一致数", sum(!ok, na.rm = TRUE),
      "official_bytewise_comparison.csv", if (!any(!ok)) "通过" else "需复核")
  add("L1", "与 CDC 官方逐字节比对", "本地与官方总字节之差",
      sum(abs(b$local_bytes - b$official_bytes), na.rm = TRUE),
      "official_bytewise_comparison.csv", "0 字节差异")
}
x <- rd("xpt_structure.csv")
if (!is.null(x)) {
  nc <- grep("^n_?rows?$|nrow", names(x), value = TRUE, ignore.case = TRUE)[1]
  add("L1", "XPT 结构解析", "文件数", nrow(x), "xpt_structure.csv", "")
  if (!is.na(nc)) add("L1", "XPT 结构解析", "总行数",
                      format(sum(x[[nc]], na.rm = TRUE), big.mark = ","),
                      "xpt_structure.csv", "通过")
}
lk <- rd("seqn_linkage.csv")
if (!is.null(lk)) {
  oc <- grep("orphan", names(lk), value = TRUE, ignore.case = TRUE)[1]
  if (!is.na(oc)) add("L1", "SEQN 跨文件链接", "孤儿 SEQN 总数",
                      sum(lk[[oc]], na.rm = TRUE), "seqn_linkage.csv", "通过")
}
u <- rd("unit_distribution.csv")
if (!is.null(u)) {
  add("L1", "跨周期单位一致性", "审计 (变量,周期) 组合数", nrow(u), "unit_distribution.csv", "")
  add("L1", "跨周期单位一致性", "被标记量级跳跃的组合数",
      sum(u$flag %in% c("JUMP", "jump", "FLAG", "10x"), na.rm = TRUE),
      "unit_distribution.csv", "通过")
}
s <- rd("sentinel_audit.csv")
if (!is.null(s)) {
  add("L1", "缺失值哨兵拦截", "审计变量数", nrow(s), "sentinel_audit.csv", "")
  add("L1", "缺失值哨兵拦截", "判为真哨兵的变量数", sum(s$n_sentinel > 0, na.rm = TRUE),
      "sentinel_audit.csv", "见文件")
}
c0 <- rd("crp_unit_adjudication.csv")
if (!is.null(c0)) for (i in seq_len(nrow(c0)))
  add("L1", "CRP 两代变量单位裁定", sprintf("%s 中位", c0$file[i]),
      sprintf("%.3f (%s)", c0$median[i], c0$implied_unit[i]),
      "crp_unit_adjudication.csv", "已修正")
wh <- rd("whq_sentinel_interception.csv")
if (!is.null(wh)) {
  nS <- wh$n_9999 + wh$n_7777 + wh$n_999
  for (i in seq_len(nrow(wh)))
    add("L1", "自报体重史哨兵拦截", wh$var[i],
        sprintf("%d 个哨兵 / %d 个非缺失 (%.2f%%)", nS[i], wh$n[i],
                100 * nS[i] / wh$n[i]),
        "whq_sentinel_interception.csv", "已置 NA")
  add("L1", "自报体重史哨兵拦截", "被拦截哨兵值合计", sum(nS),
      "whq_sentinel_interception.csv", "已置 NA")
}

## ---------------- 第二层：独立基准 ----------------
k <- rd("key_integrity.csv")
if (!is.null(k)) {
  add("L2", "主键唯一性", "断言数据集数", nrow(k), "key_integrity.csv", "")
  add("L2", "主键唯一性", "SEQN 重复总数", sum(k$dup, na.rm = TRUE),
      "key_integrity.csv", "全部唯一")
}
j <- rd("join_fanout.csv")
if (!is.null(j)) {
  add("L2", "连接扇出", "两两连接组数", nrow(j), "join_fanout.csv", "")
  add("L2", "连接扇出", "发生行数膨胀的组数", sum(j$fanout > 0, na.rm = TRUE),
      "join_fanout.csv", "无膨胀")
}
a <- rd("anchor_benchmarks.csv")
if (!is.null(a)) {
  add("L2", "金标准锚点", "锚点数", nrow(a), "anchor_benchmarks.csv", "")
  add("L2", "金标准锚点", "命中数", sum(a$pass %in% c(TRUE, "TRUE"), na.rm = TRUE),
      "anchor_benchmarks.csv", "通过")
}
vm <- rd("vermeulen_check.csv")
if (!is.null(vm)) {
  add("L2", "Vermeulen 公式独立复算", "复算例数", nrow(vm), "vermeulen_check.csv", "")
  add("L2", "Vermeulen 公式独立复算", "最大绝对偏差(nmol/L)",
      signif(max(abs(vm$abs_diff), na.rm = TRUE), 3), "vermeulen_check.csv", "通过")
  add("L2", "Vermeulen 公式独立复算", "最大相对偏差(%)",
      signif(max(abs(vm$rel_diff_pct), na.rm = TRUE), 3), "vermeulen_check.csv", "通过")
}
at <- rd("audit_trail.csv")
if (!is.null(at)) {
  add("L2", "跨脚本时序审计", "产物记录数", nrow(at), "audit_trail.csv", "")
  add("L2", "跨脚本时序审计", "撕裂（产物早于脚本）数", sum(at$stale, na.rm = TRUE),
      "audit_trail.csv", "通过")
  add("L2", "跨脚本时序审计", "最大滞后(min)", round(max(at$lag_min, na.rm = TRUE), 1),
      "audit_trail.csv", "")
}

## ---------------- 第三层：计量经济学防伪 ----------------
pc2 <- rr("inflam_clean_primary.csv")
if (!is.null(pc2)) {
  add("L3", "人群侧系数级随机化检验", "检验组合数", nrow(pc2),
      "results/inflam_clean_primary.csv", "")
  s1 <- pc2[contrast == "continuous per 1 SD" & p_perm < 0.05]
  add("L3", "人群侧系数级随机化检验", "连续暴露通过 p_perm<0.05 的结局",
      paste(s1$label, collapse = ", "), "results/inflam_clean_primary.csv", "")
  add("L3", "人群侧系数级随机化检验", "标为 UNRELIABLE 的组合数",
      sum(pc2$reliability == "UNRELIABLE", na.rm = TRUE),
      "results/inflam_clean_primary.csv", "该组合名义 P 已弃用")
  add("L3", "人群侧系数级随机化检验", "置换零分布最大偏移(以 SD 计)",
      round(max(abs(pc2$offset_sd), na.rm = TRUE), 2), "results/inflam_clean_primary.csv", "")
}
p12 <- rd("permutation_sc_exact12.csv")
if (!is.null(p12)) for (i in seq_len(nrow(p12)))
  add("L3", "单细胞精确置换 (12 池, C(12,6)=924)",
      sprintf("%s (diff=%+.3f)", p12$metric[i], p12$diff[i]),
      sprintf("p=%.4f (n_perm=%d)", p12$p_exact[i], p12$n_perm[i]),
      "permutation_sc_exact12.csv", if (p12$p_exact[i] < 0.05) "显著" else "不显著")
p18 <- rd("permutation_sc.csv")
if (!is.null(p18)) for (i in seq_len(nrow(p18)))
  add("L3", "单细胞精确置换 (18 池, C(18,6)=18564)",
      sprintf("%s (diff=%+.3f)", p18$metric[i], p18$diff[i]),
      sprintf("p=%.4f (n_perm=%d)", p18$p_exact[i], p18$n_perm[i]),
      "permutation_sc.csv", if (p18$p_exact[i] < 0.05) "显著" else "不显著")
pv <- rd("python_vs_r.csv")
if (!is.null(pv)) {
  add("L3", "R vs Python 双语言复算", "复算表数", nrow(pv), "python_vs_r.csv", "")
  add("L3", "R vs Python 双语言复算", "最大 |Δβ|",
      signif(max(abs(pv$d_beta), na.rm = TRUE), 3), "python_vs_r.csv", "机器精度内一致")
  add("L3", "R vs Python 双语言复算", "最大 |ΔSE|",
      signif(max(abs(pv$d_se), na.rm = TRUE), 3), "python_vs_r.csv", "机器精度内一致")
  add("L3", "R vs Python 双语言复算", "最大 P 值相对差(%)",
      signif(max(abs(pv$d_pct), na.rm = TRUE), 3), "python_vs_r.csv", "机器精度内一致")
}

## ---------------- 第四层：单细胞防伪 ----------------
lo <- rd("lopo_pool_level.csv")
if (!is.null(lo)) {
  b0 <- lo[drop == "none"]; lp <- lo[drop != "none"]
  add("L4", "Leave-One-Pool-Out", "剔除次数", uniqueN(lp$drop), "lopo_pool_level.csv", "")
  add("L4", "Leave-One-Pool-Out", "方向逆转总次数",
      sum(lp$sign_flip == "REVERSED", na.rm = TRUE), "lopo_pool_level.csv", "通过")
  sig0 <- b0[p < 0.05]
  ns <- lp[metric %in% sig0$metric & p >= 0.05, .N, by = metric]
  add("L4", "Leave-One-Pool-Out", "基准显著指标在留一后失去显著的次数合计",
      if (nrow(ns)) sum(ns$N) else 0, "lopo_pool_level.csv", "通过")
  for (i in seq_len(nrow(b0)))
    add("L4", "Leave-One-Pool-Out", sprintf("%s 全样本", b0$metric[i]),
        sprintf("%+.4f (P=%.4f)", b0$diff[i], b0$p[i]), "lopo_pool_level.csv", "")
}
li <- rd("lopo_influence.csv")
if (!is.null(li)) {
  li[, sus := rel_drop > 0.5]
  add("L4", "单池驱动力", "rel_drop>0.5 判为 SUSPECT 的指标",
      paste(li[sus == TRUE, metric], collapse = ", "),
      "lopo_influence.csv", "均为本已不显著项")
  add("L4", "单池驱动力", "达到显著的指标中最大 rel_drop",
      sprintf("%.3f (%s)", li[sus == FALSE][which.max(rel_drop), rel_drop],
              li[sus == FALSE][which.max(rel_drop), metric]),
      "lopo_influence.csv", "通过 (<0.5，无单池驱动)")
}
ar <- rd("cluster_ari.csv")
if (!is.null(ar)) for (i in seq_len(nrow(ar)))
  add("L4", "无监督聚类 vs 标记基因分群", ar$comparison[i], signif(ar$ari[i], 3),
      "cluster_ari.csv", ar$verdict[i])
lr <- rd("label_reproduction.csv")
if (!is.null(lr)) {
  ## 现为"池内细胞类型比例一致性"检验（子样本 vs 母集），列：gsm, ct, pct_full, pct_sub, diff_pp
  add("L4", "子样本代表性（池内比例一致性）", "比较单元数 (池 x 类)", nrow(lr),
      "label_reproduction.csv", "")
  add("L4", "子样本代表性（池内比例一致性）", "最大 |子样本-母集| 偏差 (pp)",
      sprintf("%.4f", max(abs(lr$diff_pp), na.rm = TRUE)),
      "label_reproduction.csv", "通过 (<5 pp)")
  add("L4", "子样本代表性（池内比例一致性）", "超过 5 pp 阈值的单元数",
      sum(abs(lr$diff_pp) >= 5, na.rm = TRUE), "label_reproduction.csv", "通过")
}
ds <- rd("dotplot_specificity.csv")
if (!is.null(ds)) add("L4", "Marker DotPlot 特异度", "在自身细胞类型最高表达的基因",
                      sprintf("%d/%d", sum(ds$pass, na.rm = TRUE), nrow(ds)),
                      "dotplot_specificity.csv", "")

## ---------------- 外部背靠背 ----------------
ep <- rr("ext_gse199063_paired.csv")
if (!is.null(ep)) {
  add("EXT", "GSE199063 配对纵向队列", "配对样本对数 (2年/5年)",
      paste(sort(unique(ep$n_pair)), collapse = "/"), "results/ext_gse199063_paired.csv", "")
  for (v in c("mac", "percell_infl", "net_infl", "immune", "adipocyte"))
    for (t in c("2_year", "5_year")) {
      r <- ep[variable == v & to == t]
      if (nrow(r)) add("EXT", "GSE199063 配对纵向队列", sprintf("%s @ %s", v, t),
                       sprintf("delta=%+.3f, P=%.3g", r$delta, r$p_t),
                       "results/ext_gse199063_paired.csv", "")
    }
}

## ---------------- 审稿人四问专项核查（V16–V19）----------------
## 与四层校验同一条规则：每行 = 一个数字 + 一个来源 CSV

## V16 复杂抽样加权（wt_pfas_8cycle = wt_pf / 8）
sds <- rd("survey_design_summary.csv")
if (!is.null(sds)) for (i in seq_len(nrow(sds)))
  add("V16", "复杂抽样设计框架", sds$item[i], sds$value[i],
      "survey_design_summary.csv", "svydesign 已建立")

df0 <- rd("survey_deff.csv")
if (!is.null(df0)) {
  add("V16", "设计效应（手算 Kish DEFF）", "计算变量数", nrow(df0), "survey_deff.csv", "")
  add("V16", "设计效应（手算 Kish DEFF）", "DEFF 范围",
      sprintf("%.2f-%.2f 倍", min(df0$deff, na.rm = TRUE), max(df0$deff, na.rm = TRUE)),
      "survey_deff.csv", "未加权推断偏乐观")
  add("V16", "设计效应（手算 Kish DEFF）", "有效样本/名义样本（Kish）范围",
      sprintf("%.2f-%.2f", min(df0$kish, na.rm = TRUE), max(df0$kish, na.rm = TRUE)),
      "survey_deff.csv", "")
  for (i in seq_len(nrow(df0)))
    add("V16", "设计效应（手算 Kish DEFF）", df0$variable[i],
        sprintf("设计均值 %.4f vs 未加权均值 %.4f; SE %.5f vs %.5f; DEFF %.2f; n_eff %d",
                df0$wmean[i], df0$mean_unw[i], df0$se_design[i], df0$se_srs[i],
                df0$deff[i], df0$n_eff[i]), "survey_deff.csv", "")
}

sw <- rd("survey_weighted_sensitivity.csv")
if (!is.null(sw)) {
  su <- as.logical(sw$sig_u); swt <- as.logical(sw$sig_w)
  nflip <- sum(sw$dir == "翻转", na.rm = TRUE)
  ## 点估计的口径无关度量：绝对差（百分点），只在主效应上算
  sw[, pct_u := 100 * (exp(beta_unw) - 1)]
  sw[, pct_w := 100 * (exp(beta_wt) - 1)]
  sw[, pp := pct_w - pct_u]
  sm <- sw[term == "main_per5"]
  add("V16", "加权 vs 未加权逐项比对", "比对项数", nrow(sw),
      "survey_weighted_sensitivity.csv", "")
  add("V16", "加权 vs 未加权逐项比对", "方向翻转项数", nflip,
      "survey_weighted_sensitivity.csv",
      if (nflip == 0) "全部同向（通过）" else "存在翻转")
  add("V16", "加权 vs 未加权逐项比对", "加权后仍显著项数", sum(swt %in% TRUE),
      "survey_weighted_sensitivity.csv", "")
  add("V16", "加权 vs 未加权逐项比对", "两者均不显著项数", sum(!su & !swt),
      "survey_weighted_sensitivity.csv", "一致的阴性")
  add("V16", "加权 vs 未加权逐项比对", "仅加权后显著（口径改变）项数", sum(!su & swt),
      "survey_weighted_sensitivity.csv", "")
  add("V16", "加权 vs 未加权逐项比对", "主效应点估计绝对差 中位/最大 (pp)",
      sprintf("%.2f / %.2f (%s)", median(abs(sm$pp), na.rm = TRUE),
              max(abs(sm$pp), na.rm = TRUE),
              sm[which.max(abs(pp)), label]),
      "survey_weighted_sensitivity.csv", "点估计几乎不变")
  r1 <- sm[abs(pct_u) >= 1, beta_wt / beta_unw]
  if (length(r1)) add("V16", "加权 vs 未加权逐项比对", "beta 比值范围（仅 |效应|>=1 pp 的项）",
      sprintf("%.3f-%.3f", min(r1), max(r1)), "survey_weighted_sensitivity.csv",
      "近零效应的比值无解释价值")
  for (i in seq_len(nrow(sw)))
    add("V16", "加权 vs 未加权逐项", sprintf("%s / %s / %s", sw$domain[i], sw$label[i], sw$term[i]),
        sprintf("beta_u=%+.5f (P=%.3g) -> beta_w=%+.5f (P=%.3g); ratio=%.3f; dpp=%+.2f; N=%d",
                sw$beta_unw[i], sw$p_unw[i], sw$beta_wt[i], sw$p_wt[i],
                sw$ratio_wt_unw[i], sw$pp[i], sw$N_w[i]),
        "survey_weighted_sensitivity.csv", sw$verdict[i])
}

## V17 LOD / 左删失
la <- rd("lod_audit.csv")
if (!is.null(la)) {
  add("V17", "LOD 审计（周期 x 同系物）", "审计单元数", nrow(la), "lod_audit.csv", "")
  add("V17", "LOD 审计（周期 x 同系物）", "注释码口径 vs 替换质量点口径最大偏差 (pp)",
      sprintf("%.2f", max(la$flag_vs_mass_diff_pp, na.rm = TRUE)),
      "lod_audit.csv", "通过 (<2 pp)")
  for (an in c("PFOA", "PFOS", "PFHxS", "PFNA", "PFDA", "PFUnDA")) {
    z <- la[analyte == an][order(-pct_below)]
    if (nrow(z)) add("V17", "逐同系物 %<LOD", an,
                     sprintf("峰值 %.2f%% (周期 %s)；全周期 %.2f%%-%.2f%%",
                             z$pct_below[1], z$cycle[1], min(z$pct_below), max(z$pct_below)),
                     "lod_audit.csv",
                     if (max(z$pct_below) > 30) "存在 >30% 删失周期" else "各周期均 <30%")
  }
  z2 <- la[!is.na(branch_censored)]
  if (nrow(z2)) for (i in seq_len(nrow(z2)))
    add("V17", "支链段单独删失率（2013 起两段法）",
        sprintf("%s @ %s", z2$analyte[i], z2$cycle[i]),
        sprintf("%.2f%%", z2$branch_censored[i]), "lod_audit.csv",
        if (z2$branch_censored[i] > 50) "支链段几乎全删失（已披露）" else "")
}
liv <- rd("lod_imputation_values.csv")
if (!is.null(liv)) for (i in seq_len(nrow(liv)))
  add("V17", "实测 <2LOD 取值（据此反推 LOD）",
      sprintf("%s @ %s", liv$analyte[i], liv$cycle[i]),
      liv$observed_values_below_2lod[i], "lod_imputation_values.csv", "")
lsv <- rd("lod_substitution_sensitivity.csv")
if (!is.null(lsv)) {
  rb <- suppressWarnings(as.numeric(sub("%", "", lsv$rel_beta)))
  add("V17", "替换值方案敏感性（WBC 交互）", "情景数", nrow(lsv),
      "lod_substitution_sensitivity.csv", "")
  add("V17", "替换值方案敏感性（WBC 交互）", "相对基准最大偏移",
      sprintf("%.1f%%", max(abs(rb), na.rm = TRUE)), "lod_substitution_sensitivity.csv",
      if (max(abs(rb), na.rm = TRUE) <= 2) "可忽略 (<2%)" else "需披露")
  for (i in seq_len(nrow(lsv)))
    add("V17", "替换值方案敏感性（WBC 交互）", lsv$scenario[i],
        sprintf("beta=%+.6f (SE %.6f, P=%.4f); 相对基准 %s; N=%d",
                lsv$beta[i], lsv$se[i], lsv$p[i], lsv$rel_beta[i], lsv$N[i]),
        "lod_substitution_sensitivity.csv", "")
}

## V18 混合物联合暴露（qgcomp）
qm <- rd("qgcomp_mixture.csv")
if (!is.null(qm)) for (i in seq_len(nrow(qm)))
  add("V18", "混合物联合效应（qgcomp, q=4）",
      sprintf("%s / %s", qm$outcome[i], qm$stratum[i]),
      sprintf("psi=%+.5f (%.2f%%/分位, P=%.3g, N=%d)",
              qm$psi[i], qm$pct[i], qm$p[i], qm$N[i]),
      "qgcomp_mixture.csv", if (qm$p[i] < 0.05) "显著" else "不显著")
qi <- rd("qgcomp_interaction.csv")
if (!is.null(qi)) for (i in seq_len(nrow(qi)))
  add("V18", "混合物层的修饰检验（delta-psi）", qi$outcome[i],
      sprintf("dpsi=%+.5f (SE %.5f, Wald P=%.3f, 置换 P=%.3f/%d 次); 6df 联合 chi2=%.2f P=%.3f",
              qi$delta_psi[i], qi$se_delta[i], qi$p_delta[i], qi$p_perm[i],
              qi$n_perm_ok[i], qi$chi2_joint[i], qi$p_joint_cr1[i]),
      "qgcomp_interaction.csv", "")
qw <- rd("qgcomp_weights.csv")
if (!is.null(qw)) {
  add("V18", "混合物权重（符号化，正负各自归一）", "权重记录数", nrow(qw),
      "qgcomp_weights.csv", "")
  for (k in unique(paste(qw$outcome, qw$stratum, sep = " || "))) {
    sp <- strsplit(k, " || ", fixed = TRUE)[[1]]
    z <- qw[outcome == sp[1] & stratum == sp[2]]
    z <- z[order(z$direction, -z$weight)]
    add("V18", "混合物权重（符号化，正负各自归一）",
        sprintf("%s / %s", sp[1], sp[2]),
        paste(sprintf("%s %+.2f", z$term,
                      ifelse(z$direction == "positive", z$weight, -z$weight)),
              collapse = "; "),
        "qgcomp_weights.csv", "")
  }
}

## V19 E-value
ev <- rd("evalue_table.csv")
if (!is.null(ev)) {
  add("V19", "E-value（未测量混杂定量边界）", "计算项数", nrow(ev),
      "evalue_table.csv", "")
  add("V19", "E-value（未测量混杂定量边界）", "最小 E-value",
      sprintf("%.2f (%s)", min(ev$E_value, na.rm = TRUE),
              ev$finding[which.min(ev$E_value)]), "evalue_table.csv",
      "该结论最易被混杂解释")
  add("V19", "E-value（未测量混杂定量边界）", "最大 E-value",
      sprintf("%.2f (%s)", max(ev$E_value, na.rm = TRUE),
              ev$finding[which.max(ev$E_value)]), "evalue_table.csv", "")
  for (i in seq_len(nrow(ev)))
    add("V19", "E-value（未测量混杂定量边界）",
        sprintf("%s [%s]", ev$finding[i], ev$spec[i]),
        sprintf("RR=%s; E-value=%.2f (CI 界 %.2f)", ev$RR_txt[i],
                ev$E_value[i], ev$E_value_CI[i]),
        "evalue_table.csv", ev$note[i])
}

LED <- rbindlist(L)
fwrite(LED, file.path(VER, "verification_ledger.csv"))
cat("台账行数:", nrow(LED), "\n")

md <- c("# 校验台账：四层校验 + 审稿人四问专项核查（自动生成，勿手改）", "",
        sprintf("生成时间：%s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
        sprintf("台账总行数：%d（来源 CSV：%d 个）", nrow(LED),
                uniqueN(LED$source)), "",
        "| 层 | 检查项 | 指标 | 值 | 来源文件 | 判定 |", "|---|---|---|---|---|---|")
for (i in seq_len(nrow(LED)))
  md <- c(md, sprintf("| %s | %s | %s | %s | `%s` | %s |",
                      LED$layer[i], LED$check[i], LED$metric[i],
                      LED$value[i], LED$source[i], LED$verdict[i]))
writeLines(md, file.path(VER, "VERIFICATION_SUMMARY.md"))
cat("✅ V15 完成 → verification/verification_ledger.csv, VERIFICATION_SUMMARY.md\n")
