## ============================================================================
## V11a_lopo_pool.R — 第四层①：Leave-One-Pool-Out（留一池交叉验证）
##   目的：证明单细胞结论不是被某个极端汇池（outlier pool）拽出来的
##   做法：18 个汇池逐一剔除，每次重算
##     (a) 巨噬细胞占比 WL - Obese
##     (b) 留存巨噬细胞炎症 z 评分 WL - Obese
##     (c) 促炎输出（配体口径）WL - Obese
##     (d) CSF1 -> CSF1R 配体-受体活性 WL - Obese
##     (e) 净炎症负荷（比例 × 每细胞 z）
##   判据：18 次剔除中，效应方向逆转 0 次；显著性丧失需列明具体池号与原因
## ============================================================================
suppressMessages({library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
RES  <- file.path(BASE, "results"); VER <- file.path(BASE, "verification")
dir.create(VER, showWarnings = FALSE, recursive = TRUE)

## ---------- 读入池水平数据 --------------------------------------------------
cmp <- fread(file.path(RES, "geo_composition.csv"))          # gsm,group,celltype,N,pct
mac <- fread(file.path(RES, "geo_macrophage_modules.csv"))   # gsm,group,lipolysis..antioxid
m12 <- fread(file.path(RES, "sc_m1m2.csv"))                  # gsm,group,M1,M2,inflam_output
lr  <- fread(file.path(RES, "sc_lr_pairs.csv"))              # 汇总，不够用，需重算

## CSF1 -> CSF1R 需从假体 bulk 重算（LR 汇总表只有组均值）
pb <- fread(file.path(RES, "geo_pseudobulk_counts.csv.gz"))
pb[, lib := sum(expr), by = .(gsm, celltype)]
pb[, lcpm := log10(expr / lib * 1e6 + 1)]
agg <- pb[, .(lcpm = mean(lcpm)), by = .(gsm, celltype, gene)]
M <- dcast(agg, gene ~ gsm + celltype, value.var = "lcpm")
G <- M$gene; M <- as.matrix(M[, -1]); rownames(M) <- G
GRP <- unique(pb[, .(gsm, group)])[order(gsm)]
E <- function(gsm, ct, genes) {
  cols <- paste0(gsm, "_", ct); genes <- intersect(genes, rownames(M))
  if (!length(genes)) return(rep(NA_real_, length(gsm)))
  colMeans(M[genes, cols, drop = FALSE], na.rm = TRUE)
}
csf <- data.table(gsm = GRP$gsm, group = GRP$group,
                  score = (E(GRP$gsm, "adipocyte", "CSF1") +
                           E(GRP$gsm, "macrophage", "CSF1R")) / 2)

## ---------- 组装池水平长表 --------------------------------------------------
P <- Reduce(function(a, b) merge(a, b, by = c("gsm", "group")), list(
  cmp[celltype == "macrophage", .(gsm, group, mac_pct = pct)],
  mac[, .(gsm, group, inflam_z = inflamma)],
  m12[, .(gsm, group, m1m2 = M1 - M2, inflam_out = inflam_output)],
  csf[, .(gsm, group, csf1_csf1r = score)]
))
P[, burden := mac_pct / 100 * inflam_z]
cat("池水平数据：", nrow(P), "池\n"); print(P[, .N, by = group])

VARS <- c("mac_pct", "inflam_z", "m1m2", "inflam_out", "csf1_csf1r", "burden")

## ---------- 全样本基准 ------------------------------------------------------
stat <- function(d, v) {
  a <- d[[v]][d$group == "WL"]; b <- d[[v]][d$group == "Obese"]
  if (length(a) < 2 || length(b) < 2) return(c(diff = NA, p = NA, t = NA))
  tt <- t.test(a, b)
  c(diff = mean(a) - mean(b), p = tt$p.value, t = unname(tt$statistic))
}
full <- rbindlist(lapply(VARS, function(v) {
  s <- stat(P, v)
  data.table(metric = v, drop = "none", diff = s[["diff"]], p = s[["p"]], t = s[["t"]],
             n_WL = sum(P$group == "WL"), n_Ob = sum(P$group == "Obese"))
}))

## ---------- 逐一剔除 --------------------------------------------------------
loo <- rbindlist(lapply(seq_len(nrow(P)), function(i) {
  d <- P[-i, ]
  rbindlist(lapply(VARS, function(v) {
    s <- stat(d, v)
    data.table(metric = v, drop = P$gsm[i], drop_grp = P$group[i],
               diff = s[["diff"]], p = s[["p"]], t = s[["t"]],
               n_WL = sum(d$group == "WL"), n_Ob = sum(d$group == "Obese"),
               drop_mac_pct = P$mac_pct[i])
  }))
}))

## ---------- 判定 ------------------------------------------------------------
OUT <- rbind(full, loo, fill = TRUE)
setcolorder(OUT, c("metric", "drop", "drop_grp", "n_WL", "n_Ob", "diff", "t", "p",
                   "drop_mac_pct"))
OUT[, `:=`(sign_flip = NA_character_, still_sig = NA_character_)]
OUT[drop != "none", sign_flip := {
  base <- full$diff[match(metric, full$metric)]
  fifelse(sign(diff) != sign(base), "REVERSED", "same")
}, by = metric]
OUT[drop != "none", still_sig := fifelse(p < 0.05, "P<0.05", "P>=0.05")]

fwrite(OUT, file.path(VER, "lopo_pool_level.csv"))

## ---------- 控制台报告 ------------------------------------------------------
cat("\n================ Leave-One-Pool-Out（18 池逐一剔除）================\n")
for (v in VARS) {
  sub <- OUT[metric == v]
  base <- sub[drop == "none"]
  lo  <- sub[drop != "none"]
  nrev <- sum(lo$sign_flip == "REVERSED", na.rm = TRUE)
  sig0 <- base$p < 0.05
  nsig <- sum(lo$p >= 0.05, na.rm = TRUE)
  worst <- lo[which.max(p)]
  cat(sprintf("\n[%s]\n", v))
  cat(sprintf("  全样本 : diff=%+.4f  t=%+.2f  P=%.4f  (WL n=%d, Obese n=%d)\n",
              base$diff, base$t, base$p, base$n_WL, base$n_Ob))
  cat(sprintf("  LOPO   : diff 范围 [%+.4f, %+.4f]  P 范围 [%.4f, %.4f]\n",
              min(lo$diff), max(lo$diff), min(lo$p), max(lo$p)))
  cat(sprintf("           方向逆转 %d/18 次；显著性丧失 %d/18 次（基准 P%s0.05）\n",
              nrev, if (sig0) nsig else 0, if (sig0) "<" else ">="))
  cat(sprintf("           最弱一次：剔除 %s (%s)，P=%.4f, diff=%+.4f\n",
              worst$drop, worst$drop_grp, worst$p, worst$diff))
}

## ---------- 极端池识别（影响最大的池） --------------------------------------
cat("\n================ 单池影响力排名（按剔除后 |t| 变化）================\n")
infl <- OUT[drop != "none", .(metric, drop, drop_grp, p,
                              dt = abs(t - full$t[match(metric, full$metric)])),
            by = .(metric, drop)]
infl <- OUT[drop != "none", .(metric, drop, drop_grp, p, t,
                              t_full = full$t[match(metric, full$metric)])]
infl[, delta_abs_t := abs(t) - abs(t_full)]
top <- infl[, .SD[which.max(delta_abs_t)], by = metric]
print(top[, .(metric, drop, drop_grp, p = signif(p, 3),
              t = round(t, 2), t_full = round(t_full, 2),
              d_abs_t = round(delta_abs_t, 2))])

## ---------- 敏感性：是否只有一个池在"撑"结果 -------------------------------
## 若结论由单池驱动，则剔除该池后 |t| 应大幅下降（|Δt|/|t_full| > 0.5）
top[, rel_drop := delta_abs_t / abs(t_full)]
cat("\n结论由单池驱动的判定（|Δt|/|t| > 0.5 视为可疑）：\n")
print(top[, .(metric, drop, rel_drop = round(rel_drop, 3),
              verdict = fifelse(rel_drop > 0.5, "SUSPECT", "robust"))])

fwrite(top, file.path(VER, "lopo_influence.csv"))
cat("\n✅ V11a 完成 → verification/lopo_pool_level.csv, lopo_influence.csv\n")
