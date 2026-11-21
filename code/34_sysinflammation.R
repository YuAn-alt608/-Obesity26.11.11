## ============================================================================
## 34_sysinflammation.R — 突破口二的人群层面验证
##   单细胞发现：WL 后巨噬细胞比例 25.3%→13.5%（撤退），但留存细胞促炎输出翻倍；
##              净炎症负荷不变；CCL2→CCR2 受体端上调 112%。
##   可检验推论：
##     H2a 减重降低循环单核细胞（巨噬细胞招募的入口）→ 单核细胞计数/百分比下降
##     H2b 减重降低系统炎症（CRP、NLR）→ 但幅度小于"巨噬细胞撤退"的预期
##     H2c 高 PFAS 负荷者减重后炎症下降更少（撤退受阻 + 留存细胞高炎症）
##     H2d "净负荷不变"在人群层面复现：单核细胞↓ 与 CRP 变化脱钩
##   结局：CRP / WBC / 中性粒 / 淋巴 / 单核 / NLR / MLR / PLR / SII
##   周期：C,D,F,G,H,I,J,L（CRP 仅 D,F,I,J,L 五周期）
## 输出：results/inflam_per5.csv / inflam_tertile.csv / inflam_decouple.csv
##       data/derived/inflam.rds
## ============================================================================
suppressMessages({library(foreign); library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", "."); RES <- file.path(BASE, "results")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))

crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X)
  if (n <= k + 2) return(matrix(NA_real_, k, k))
  bread <- tryCatch(solve(crossprod(X)), error = function(e) MASS::ginv(crossprod(X)))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
wald <- function(b, V) { if (any(is.na(V))) return(NA_real_); as.numeric(t(b) %*% MASS::ginv(V) %*% b) }
## 设计矩阵：剔除零方差列，但**必须保留截距**
## （V16 审计：初版 sd(v)>0 判据会把常数列 (Intercept) 删掉，使模型在无截距下
##   拟合，全部系数系统性偏移。此处已修正；修正后结果见 36_inflam_clean.R）
mkX <- function(f, dd) {
  X <- model.matrix(f, data = dd)
  ok <- apply(X, 2, function(v) all(is.finite(v)))
  vary <- apply(X, 2, function(v) sd(v) > 0)
  X[, ok & (vary | colnames(X) == "(Intercept)"), drop = FALSE]
}

## ---------- ① 从原始 XPT 抽取炎症指标 ---------------------------------------
STEM <- data.table(
  cycle = c("C","D","F","G","H","I","J","L"),
  cbc   = c("L25_C","CBC_D","CBC_F","CBC_G","CBC_H","CBC_I","CBC_J","CBC_L"),
  crp   = c(NA, "CRP_D", "CRP_F", NA, NA, "HSCRP_I", "HSCRP_J", "HSCRP_L")
)
CV <- c(wbc="LBXWBCSI", lym_p="LBXLYPCT", mono_p="LBXMOPCT", neut_p="LBXNEPCT",
        lym_n="LBDLYMNO", mono_n="LBDMONO", neut_n="LBDNENO", plt="LBXPLTSI")

pick <- function(f, map) {
  p <- file.path(BASE, "data/raw", paste0(f, ".xpt"))
  if (!file.exists(p)) return(NULL)
  x <- as.data.table(read.xport(p))
  keep <- intersect(unname(map), names(x))
  if (!length(keep)) return(NULL)
  out <- x[, c("SEQN", keep), with = FALSE]
  setnames(out, keep, names(map)[match(keep, map)])
  for (v in names(map)) if (!v %in% names(out)) out[, (v) := NA_real_]
  out[, c("SEQN", names(map)), with = FALSE]
}

IN <- rbindlist(lapply(seq_len(nrow(STEM)), function(i) {
  s <- STEM$cycle[i]
  a <- pick(STEM$cbc[i], CV)
  if (is.null(a)) return(NULL)
  ## CRP：LBXCRP (mg/dL? 实为 mg/L in NHANES 2003-06) 与 LBXHSCRP (mg/L)
  if (!is.na(STEM$crp[i])) {
    b <- pick(STEM$crp[i], c(crp = "LBXCRP", crp_hs = "LBXHSCRP"))
    if (!is.null(b)) a <- merge(a, b, by = "SEQN", all.x = TRUE)
  }
  if (!"crp" %in% names(a)) a[, crp := NA_real_]
  if (!"crp_hs" %in% names(a)) a[, crp_hs := NA_real_]
  a[, crp_all := fifelse(is.na(crp_hs), crp, crp_hs)]
  a[, cycle := s]
  a
}))
cat("炎症原始记录 N =", nrow(IN), " 周期 =", paste(unique(IN$cycle), collapse=","), "\n")
cat("  WBC 非缺失:", sum(!is.na(IN$wbc)), "  CRP 非缺失:", sum(!is.na(IN$crp_all)), "\n\n")

## ---------- ② 派生指标 ------------------------------------------------------
IN[, `:=`(
  nlr = neut_n/pmax(lym_n, 1e-3),
  mlr = mono_n/pmax(lym_n, 1e-3),
  plr = plt/pmax(lym_n, 1e-3),
  sii = plt*neut_n/pmax(lym_n, 1e-3)
)]
IN[, `:=`(
  y_crp   = log(pmax(crp_all, 0.05)),
  y_wbc   = log(pmax(wbc, 0.5)),
  y_neut  = log(pmax(neut_n, 0.05)),
  y_lym   = log(pmax(lym_n, 0.05)),
  y_mono  = log(pmax(mono_n, 0.01)),
  y_monop = log(pmax(mono_p, 0.1)/(100 - pmin(mono_p, 99.9))),
  y_nlr   = log(pmax(nlr, 0.05)),
  y_mlr   = log(pmax(mlr, 0.01)),
  y_plr   = log(pmax(plr, 1)),
  y_sii   = log(pmax(sii, 1))
)]
IN[, ok_infl := !is.na(wbc)]

D <- IN[A$D, on = "SEQN"]
## ---- 单位归一（V03 审计发现并修正）----------------------------------------
## LBXCRP 单位为 mg/dL（2003-2010），LBXHSCRP 为 mg/L（2015-2020），
## 二者直接合并会造成 10 倍量级断裂（D 周期中位 0.20 vs I 周期 1.90）。
## A$D$crp_mgl 已在 05_build 按 1 mg/dL = 10 mg/L 统一到 mg/L，此处优先采用。
D[, crp_all := fifelse(!is.na(crp_mgl), crp_mgl, crp_all)]
D[, y_crp   := log(pmax(crp_all, 0.05))]
cat("CRP 单位归一后 各周期中位(mg/L):",
    paste(D[, .(m = round(median(crp_all, na.rm = TRUE), 2)), by = cycle][order(cycle)][[2]], collapse = " / "), "\n")
saveRDS(D, file.path(BASE, "data/derived/inflam.rds"))

base <- D[ok_infl == TRUE & !is.na(w5) & !is.na(pfas_z)]
base[, `:=`(mv = as.numeric(pfas_z), dw5 = dw1/5)]
base[, mv2 := mv^2]
base[, pfas_t := cut(mv, unique(quantile(mv, c(0,1/3,2/3,1), na.rm=TRUE)),
                     labels = c("T1","T2","T3"), include.lowest = TRUE)]
CONF <- A$CONF
cf <- CONF[sapply(CONF, function(v) length(unique(base[[v]])) > 1)]
cat("分析样本 N =", nrow(base), " 处理 =", sum(base$w5), " 聚类 =",
    uniqueN(paste(base$strata, base$psu)), " 协变量 =", length(cf), "\n\n")

YK <- c("y_crp","y_wbc","y_neut","y_lym","y_mono","y_monop","y_nlr","y_mlr","y_plr","y_sii")
LAB <- c(y_crp="CRP", y_wbc="WBC", y_neut="Neutrophil", y_lym="Lymphocyte",
         y_mono="Monocyte count", y_monop="Monocyte %", y_nlr="NLR",
         y_mlr="MLR", y_plr="PLR", y_sii="SII")

fml <- function(y, rhs) as.formula(paste(y, "~", rhs, "+", paste(cf, collapse = " + ")))

## ---------- ③ 减重主效应（每 5%） -------------------------------------------
cat("========== A. 每 5% 减重的主效应（全部结局） ==========\n")
RP <- rbindlist(lapply(YK, function(k) {
  dd <- base[!is.na(get(k))]
  X <- mkX(fml(k, "dw5 + w5"), dd); cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl)
  j <- grep("^dw5$", colnames(X)); b <- m$coefficients[j]; se <- sqrt(V[j,j])
  data.table(outcome = k, label = LAB[k], N = nrow(dd),
             per5 = 100*(exp(b)-1), lo = 100*(exp(b-1.96*se)-1), hi = 100*(exp(b+1.96*se)-1),
             p = 2*pnorm(-abs(b/se)))
}))
RP[, p_bh := p.adjust(p, "BH")]
print(RP[, .(label, N, `每5%`=sprintf("%+6.2f%%", per5),
             `95%CI`=sprintf("%+.2f, %+.2f", lo, hi), p=format.pval(p, digits=3),
             p_bh=format.pval(p_bh, digits=3))], nrows = 20)

## ---------- ④ PFAS 三分位交互 -------------------------------------------------
cat("\n========== B. PFAS 三档的减重效应（不假设单调） ==========\n")
RT <- rbindlist(lapply(YK, function(k) {
  dd <- base[!is.na(get(k)) & !is.na(pfas_t)]
  X <- mkX(fml(k, "w5*pfas_t + dw5*pfas_t"), dd); cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  nm <- colnames(X)
  jD  <- which(nm == "dw5")
  jIn <- grep("pfas_t", grep("dw5", nm, value = TRUE), value = TRUE)   # 交互列名与顺序无关
  jIn <- match(jIn, nm); jIn <- jIn[!is.na(jIn)]
  eff <- sapply(c("T1","T2","T3"), function(t) {
    j <- which(nm == paste0("pfas_t", t, ":dw5") | nm == paste0("dw5:pfas_t", t))
    if (t == "T1") return(100*(exp(b[jD])-1))
    if (!length(j)) return(NA_real_)
    100*(exp(b[jD] + b[j]) - 1)
  })
  data.table(outcome=k, label=LAB[k], T1=eff[1], T2=eff[2], T3=eff[3],
             p_int = if (length(jIn)) pchisq(wald(b[jIn], V[jIn,jIn,drop=FALSE]), length(jIn),
                                             lower.tail=FALSE) else NA_real_)
}))
RT[, p_bh := p.adjust(p_int, "BH")]
print(RT[, .(label, T1=sprintf("%+6.2f%%",T1), T2=sprintf("%+6.2f%%",T2), T3=sprintf("%+6.2f%%",T3),
             `交互P`=format.pval(p_int, digits=3))], nrows=20)

## ---------- ⑤ 连续剂量交互（线性 + 二次） -------------------------------------
cat("\n========== C. 连续 PFAS 剂量交互 ==========\n")
RQ <- rbindlist(lapply(YK, function(k) {
  dd <- base[!is.na(get(k))]
  X <- mkX(fml(k, "dw5*(mv+mv2) + mv + mv2"), dd); cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  j1 <- grep("^dw5:mv$", colnames(X)); j2 <- grep("^dw5:mv2$", colnames(X))
  data.table(outcome=k, label=LAB[k],
             p_lin = if (length(j1)) 2*pnorm(-abs(b[j1]/sqrt(V[j1,j1]))) else NA_real_,
             p_quad = if (length(c(j1,j2))==2) pchisq(wald(b[c(j1,j2)], V[c(j1,j2),c(j1,j2)]), 2, lower.tail=FALSE) else NA_real_)
}))
print(RQ[, .(label, `线性P`=format.pval(p_lin, digits=3), `二次P`=format.pval(p_quad, digits=3))], nrows=20)

## ---------- ⑥ 脱钩检验：单核细胞变化 vs CRP 变化 -------------------------------
cat("\n========== D. 脱钩检验（标准化后交互系数之差） ==========\n")
## 若"巨噬细胞撤退"与"系统炎症消退"同步，两者的 PFAS 交互应同向等量
dsub <- base[!is.na(y_crp) & !is.na(y_mono)]
z <- function(v) as.numeric(scale(v))
dsub[, `:=`(z_crp = z(y_crp), z_mono = z(y_mono), z_nlr = z(y_nlr))]
dec <- function(k) {
  dd <- dsub[!is.na(get(k))]
  X <- mkX(fml(k, "dw5*(mv+mv2) + mv + mv2"), dd); cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  j1 <- grep("^dw5:mv$", colnames(X))
  c(est = unname(b[j1]), se = unname(sqrt(V[j1,j1])))
}
pairs <- list(c("z_mono","z_crp"), c("z_nlr","z_crp"), c("z_mono","z_nlr"))
RD <- rbindlist(lapply(pairs, function(pp) {
  a <- dec(pp[1]); b <- dec(pp[2])
  data.table(pair = paste(pp[1], "-", pp[2]), diff = a["est"]-b["est"],
             se = sqrt(a["se"]^2 + b["se"]^2),
             p = 2*pnorm(-abs((a["est"]-b["est"])/sqrt(a["se"]^2+b["se"]^2))))
}))
print(RD[, .(pair, diff=round(diff,3), se=round(se,3), p=format.pval(p, digits=3))])

## ---------- ⑦ 关键诊断：交互是否只是"剂量构成差异" ---------------------------
cat("\n========== E. 关键诊断：PFAS 是否修饰减重幅度本身 ==========\n")
base[, dbmi := bmi - bmi_p]      # 实测 BMI 变化（负 = 下降）
diag <- rbindlist(lapply(c("dw1","dbmi","waist"), function(k) {
  dd <- base[!is.na(get(k))]
  X <- mkX(fml(k, "dw5*(mv+mv2) + mv + mv2"), dd); cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  jD <- grep("^dw5$", colnames(X)); j1 <- grep("^dw5:mv$", colnames(X)); j2 <- grep("^dw5:mv2$", colnames(X))
  data.table(term=k, per5=b[jD], p_lin=2*pnorm(-abs(b[j1]/sqrt(V[j1,j1]))),
             p_joint = if (length(c(j1,j2))==2)
               pchisq(wald(b[c(j1,j2)], V[c(j1,j2),c(j1,j2)]), 2, lower.tail=FALSE) else NA_real_)
}))
print(diag[, .(term, `每5%`=round(per5,3), `线性P`=format.pval(p_lin, digits=3),
               `联合P`=format.pval(p_joint, digits=3))])

cat("\n-- 处理组内，按 PFAS 三分位的实际减重幅度 --\n")
dose <- base[w5 == 1 & !is.na(pfas_t),
             .(N = .N, `ΔW1y 中位(%)` = round(median(dw1, na.rm=TRUE),1),
               `≥10% 比例` = sprintf("%.1f%%", 100*mean(dw1 >= 10, na.rm=TRUE)),
               `ΔBMI 均值` = round(mean(dbmi, na.rm=TRUE),2),
               `腰围 均值(cm)` = round(mean(waist, na.rm=TRUE),1)), by = pfas_t]
print(dose)

cat("\n-- 改用实测 ΔBMI 作为剂量，交互是否仍存在 --\n")
dose2 <- rbindlist(lapply(c("y_wbc","y_neut","y_monop","y_sii","y_crp"), function(k) {
  dd <- base[!is.na(get(k)) & !is.na(pfas_t) & !is.na(dbmi)]
  X <- mkX(fml(k, "dbmi*pfas_t"), dd); cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  nm <- colnames(X); jD <- which(nm == "dbmi")
  jIn <- match(grep("pfas_t", grep("dbmi", nm, value = TRUE), value = TRUE), nm)
  jIn <- jIn[!is.na(jIn)]
  eff <- sapply(c("T1","T2","T3"), function(t) {
    j <- which(nm == paste0("pfas_t", t, ":dbmi") | nm == paste0("dbmi:pfas_t", t))
    if (t == "T1") return(100*(exp(-b[jD])-1))     # 每 -1 kg/m2（即减重）
    if (!length(j)) return(NA_real_)
    100*(exp(-(b[jD] + b[j])) - 1)
  })
  data.table(label=LAB[k], T1=eff[1], T2=eff[2], T3=eff[3],
             p_int = if (length(jIn)) pchisq(wald(b[jIn], V[jIn,jIn,drop=FALSE]),
                                             length(jIn), lower.tail=FALSE) else NA_real_)
}))
print(dose2[, .(label, `T1/1单位`=sprintf("%+6.2f%%",T1), `T2/1单位`=sprintf("%+6.2f%%",T2),
                `T3/1单位`=sprintf("%+6.2f%%",T3), `交互P`=format.pval(p_int, digits=3))])

fwrite(diag,  file.path(RES, "inflam_weight_diag.csv"))
fwrite(dose,  file.path(RES, "inflam_weight_dose.csv"))
fwrite(dose2, file.path(RES, "inflam_dbmi_interaction.csv"))

fwrite(RP, file.path(RES, "inflam_per5.csv"))
fwrite(RT, file.path(RES, "inflam_tertile.csv"))
fwrite(RQ, file.path(RES, "inflam_quad.csv"))
fwrite(RD, file.path(RES, "inflam_decouple.csv"))
cat("\n✅ 完成\n")
