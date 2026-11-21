## ============================================================================
## 36_inflam_clean.R — 突破口二：修正共线规格后的权威结果（supersedes 34 的 B 段）
## ----------------------------------------------------------------------------
## 背景（V08 置换检验暴露的问题）：
##   34_sysinflammation.R 的 B 段模型为
##       y ~ w5*pfas_t + dw5*pfas_t + 混杂
##   其中 dw5 = dw1/5 = w5 * (1 年体重变化% / 5)，对对照组恒为 0，
##   与 w5 构成"零膨胀"共线结构；该规格的 2 df 联合 Wald 在 CR1 下给出
##   SII P = 6.0e-11、Mono% P = 7.8e-8，而把 PFAS 随机打乱后仍给出 P~1e-6
##   （verification/specification_sensitivity.csv 的 S4）。→ 统计量被放大，是伪影。
##
## 本脚本给出的三套规格，全部移除零膨胀共线项：
##   S1  主报（DiD 型）  : y ~ w5*pfas_t + 混杂      交互 = "减重幅度随 PFAS 档位的变化"
##   S1b 主报（连续）     : y ~ w5*mv     + 混杂      交互 = 每 +1 SD PFAS
##   S2  处理组内剂量     : 仅 w5==1，y ~ dw1*pfas_t   交互 = 每减 1% 体重随 PFAS 的变化
## 每个系数同时报告：
##   p_cr1  —— CR1 聚类稳健 Wald（与 34 同口径，便于对照）
##   p_perm —— 系数级随机化检验（暴露在 周期×性别 层内打乱 B 次，不经过方差估计）
## 判据：以 p_perm 为金标准；p_cr1 与 p_perm 相差 > 1 个数量级时标记 UNRELIABLE
## 输出：results/inflam_clean_primary.csv / inflam_clean_tertile.csv
##       results/inflam_clean_dose.csv / verification/coef_stability_diag.csv
## ============================================================================
suppressMessages({library(data.table); library(MASS)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
RES  <- file.path(BASE, "results"); VER <- file.path(BASE, "verification")
set.seed(20261117)

A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
D <- readRDS(file.path(BASE, "data/derived/inflam.rds"))
CONF <- A$CONF

crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X)
  if (n <= k + 2) return(matrix(NA_real_, k, k))
  bread <- tryCatch(solve(crossprod(X)), error = function(e) ginv(crossprod(X)))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx)
  list(V = bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k)),
       meat = meat, G = G)
}
mkX <- function(f, dd) {
  X <- model.matrix(f, data = dd)
  ## 只剔除"含非有限值"或"零方差"的列，但**必须保留截距**
  ## （V16 审计发现：初版 apply(... sd(v)>0) 会把常数列 (Intercept) 一并删掉，
  ##   导致全部模型在无截距下拟合，系数系统性偏移约 20%）
  ok <- apply(X, 2, function(v) all(is.finite(v)))
  vary <- apply(X, 2, function(v) sd(v) > 0)
  keep <- ok & (vary | colnames(X) == "(Intercept)")
  X[, keep, drop = FALSE]
}
w2p <- function(b, V, j) if (!length(j) || any(is.na(V[j, j]))) NA_real_ else
  as.numeric(t(b[j]) %*% ginv(V[j, j, drop = FALSE]) %*% b[j])

base <- D[ok_infl == TRUE & !is.na(w5) & !is.na(pfas_z)]
base[, `:=`(mv = as.numeric(pfas_z), dw5 = dw1/5)]
base[, pfas_t := cut(mv, unique(quantile(mv, c(0, 1/3, 2/3, 1), na.rm = TRUE)),
                     labels = c("T1","T2","T3"), include.lowest = TRUE)]
base[, stratum := paste(cycle, sex)]
cf <- CONF[sapply(CONF, function(v) length(unique(base[[v]])) > 1)]
cat("分析样本 N =", nrow(base), " 处理 =", sum(base$w5),
    " 聚类 =", uniqueN(paste(base$strata, base$psu)), " 协变量 =", length(cf), "\n\n")

YK <- c(y_crp="CRP", y_wbc="WBC", y_neut="Neutrophil", y_lym="Lymphocyte",
        y_mono="Monocyte count", y_monop="Monocyte %", y_nlr="NLR",
        y_mlr="MLR", y_plr="PLR", y_sii="SII")
fml <- function(y, rhs) as.formula(paste(y, "~", rhs, "+", paste(cf, collapse=" + ")))

## ============================================================================
## 0. 共线诊断：S4（旧，坏）vs S1（新，好）
## ============================================================================
cat("========== 0. 规格稳定性诊断（为什么旧规格不可用）==========\n")
diag <- rbindlist(lapply(names(YK), function(k) {
  dd <- base[!is.na(get(k))]
  cl <- factor(paste(dd$strata, dd$psu))
  out <- list()
  for (sp in c("S4_old", "S1_new")) {
    rhs <- if (sp == "S4_old") "w5*pfas_t + dw5*pfas_t" else "w5*pfas_t"
    X <- mkX(fml(k, rhs), dd)
    m <- lm.fit(X, dd[[k]]); cv <- crsV(X, m$residuals, cl)
    nm <- colnames(X)
    jI <- grep("^(w5:pfas_t|pfas_t.*:w5|dw5:pfas_t)", nm)
    jI <- jI[!grepl("^pfas_tT[23]$", nm[jI])]
    Vs <- cv$V[jI, jI, drop = FALSE]
    ev <- eigen(Vs, symmetric = TRUE, only.values = TRUE)$values
    ## meat 的 PSU 集中度：前 5 个 PSU 贡献占比
    mt <- cv$meat
    contrib <- diag(cv$meat)  # 仅近似，用于排序
    out[[sp]] <- data.table(
      spec = sp, outcome = k, label = YK[[k]], N = nrow(dd),
      kappa_X = kappa(X),
      kappa_Vsub = if (length(ev) > 1) max(ev)/min(ev) else 1,
      min_eig_Vsub = min(ev),
      p_joint = pchisq(w2p(m$coefficients, cv$V, jI), length(jI), lower.tail = FALSE),
      n_int_terms = length(jI))
  }
  rbindlist(out)
}))
DIAG <- dcast(diag, outcome + label + N + n_int_terms ~ spec,
              value.var = c("kappa_X", "kappa_Vsub", "min_eig_Vsub", "p_joint"))
DIAG[, p_ratio := p_joint_S4_old / p_joint_S1_new]
fwrite(DIAG, file.path(VER, "coef_stability_diag.csv"))
print(DIAG[, .(label, N, kappa_X_S4_old = round(kappa_X_S4_old,1),
               kappa_X_S1_new = round(kappa_X_S1_new,1),
               kappa_V_S4 = round(kappa_Vsub_S4_old,1),
               p_S4_old = signif(p_joint_S4_old,3),
               p_S1_new = signif(p_joint_S1_new,3),
               ratio = signif(p_ratio,3))], row.names = FALSE)

## ============================================================================
## 1. 系数级随机化检验（金标准）
## ============================================================================
B <- 1000
coefs_of <- function(dat, y) {
  dat <- dat[is.finite(dat[[y]]), ]
  dat$pfas_t <- cut(dat$mv, unique(quantile(dat$mv, c(0,1/3,2/3,1), na.rm = TRUE)),
                    labels = c("T1","T2","T3"), include.lowest = TRUE)
  X1 <- mkX(fml(y, "w5*pfas_t"), dat); m1 <- lm.fit(X1, dat[[y]])
  X2 <- mkX(fml(y, "w5*mv"),     dat); m2 <- lm.fit(X2, dat[[y]])
  X3 <- mkX(fml(y, "dw5*mv"),    dat); m3 <- lm.fit(X3, dat[[y]])
  n1 <- colnames(X1); n2 <- colnames(X2); n3 <- colnames(X3)
  jt <- grep("^w5:pfas_tT3$|^pfas_tT3:w5$", n1)
  jc <- grep("^w5:mv$|^mv:w5$", n2)
  jd <- grep("^dw5:mv$|^mv:dw5$", n3)
  c(tert = if (length(jt)) unname(m1$coefficients[jt]) else NA_real_,
    cont = if (length(jc)) unname(m2$coefficients[jc]) else NA_real_,
    dose = if (length(jd)) unname(m3$coefficients[jd]) else NA_real_)
}
KEY <- names(YK)
obs <- t(sapply(KEY, function(y) coefs_of(base, y)))
colnames(obs) <- c("tert", "cont", "dose")
cat("\n========== 1. 系数级随机化检验（B =", B, "）==========\n")
cat("  置换：PFAS 在 周期×性别 层内打乱，保持层结构与边缘分布\n")
cat("  三个估计量：\n")
cat("    tert = w5:pfas_tT3  （DiD 型，三档）\n")
cat("    cont = w5:mv        （DiD 型，连续，每 +1 SD）\n")
cat("    dose = dw5:mv       （剂量-反应交互，每 5% 减重 x 每 +1 SD PFAS）\n")

KIND <- c("tert", "cont", "dose")
nulls <- array(NA_real_, dim = c(B, length(KEY), length(KIND)),
               dimnames = list(NULL, KEY, KIND))
t0 <- Sys.time()
for (b in seq_len(B)) {
  dd <- copy(base)
  dd[, mv := ave(mv, stratum, FUN = function(x) sample(x))]
  for (k in KEY) nulls[b, k, ] <- coefs_of(dd, k)
  if (b %% 200 == 0) cat("    置换", b, "/", B,
                         sprintf(" [%.0fs]\n", as.numeric(difftime(Sys.time(), t0, units="secs"))))
}
saveRDS(nulls, file.path(VER, "permutation_pop_clean_null.rds"))

## CR1 名义 P（同一批规格）
LABK <- c(tert = "T3 vs T1 (DiD)", cont = "continuous per 1 SD",
          dose = "dose-response (per 5% x per SD)")
RHS  <- c(tert = "w5*pfas_t", cont = "w5*mv", dose = "dw5*mv")
JPAT <- c(tert = "^w5:pfas_tT3$|^pfas_tT3:w5$",
          cont = "^w5:mv$|^mv:w5$",
          dose = "^dw5:mv$|^mv:dw5$")
cr1p <- rbindlist(lapply(names(YK), function(k) {
  dd <- base[!is.na(get(k))]; cl <- factor(paste(dd$strata, dd$psu))
  r <- list()
  for (kind in KIND) {
    X <- mkX(fml(k, RHS[[kind]]), dd); m <- lm.fit(X, dd[[k]]); cv <- crsV(X, m$residuals, cl)
    j <- grep(JPAT[[kind]], colnames(X))
    r[[kind]] <- data.table(outcome = k, label = YK[[k]], contrast = LABK[[kind]],
                            N = nrow(dd), coef = m$coefficients[j],
                            se = sqrt(cv$V[j, j]),
                            p_cr1 = 2*pnorm(-abs(m$coefficients[j]/sqrt(cv$V[j, j]))))
  }
  rbindlist(r)
}))

PERM <- rbindlist(lapply(KEY, function(k) rbindlist(lapply(KIND, function(kind) {
  np <- nulls[, k, kind]; np <- np[is.finite(np)]
  mu <- mean(np)
  data.table(outcome = k, label = YK[[k]], contrast = LABK[[kind]],
             coef_obs = obs[k, kind],
             null_mean = mu, null_sd = sd(np),
             ## 原始（未校正）置换 P：零分布若偏移，此值不可用，仅作诊断
             p_perm_raw = (1 + sum(abs(np) >= abs(obs[k, kind])))/(length(np) + 1),
             ## 位置校正置换 P：以零分布自身均值为参照
             p_perm = (1 + sum(abs(np - mu) >= abs(obs[k, kind] - mu)))/(length(np) + 1),
             offset_sd = mu/sd(np))
}))))
PRIM <- merge(cr1p, PERM[, .(outcome, contrast, null_mean, null_sd, offset_sd,
                             p_perm, p_perm_raw)],
              by = c("outcome", "contrast"), all.x = TRUE)
PRIM[, `:=`(z_perm = (coef - null_mean)/null_sd,
            log10_ratio = log10(p_cr1/p_perm),
            reliability = fifelse(abs(log10(p_cr1/p_perm)) > 1, "UNRELIABLE",
                            fifelse(abs(log10(p_cr1/p_perm)) > 0.5, "marginal", "ok")))]
setorder(PRIM, contrast, outcome)
fwrite(PRIM, file.path(RES, "inflam_clean_primary.csv"))
fwrite(PERM, file.path(VER, "permutation_pop_clean.csv"))

cat("\n-- 零分布形态诊断（null_mean 应为 0；|offset_sd| > 0.5 说明该规格有偏）--\n")
print(PRIM[, .(label, contrast, null_mean = signif(null_mean,3),
               null_sd = signif(null_sd,3), offset_sd = round(offset_sd,2))],
      row.names = FALSE)

cat("\n-- 主报结果（S1 规格，CR1 名义 P vs 随机化 P）--\n")
print(PRIM[, .(label, contrast, N, coef = signif(coef,3), se = signif(se,3),
               p_cr1 = signif(p_cr1,3), p_perm = signif(p_perm,3),
               z_perm = round(z_perm,2), reliability)], row.names = FALSE)

## ============================================================================
## 2. 三档分档效应（S1 规格，无共线项）
## ============================================================================
cat("\n========== 2. 各 PFAS 档位的减重效应（S1：y ~ w5*pfas_t）==========\n")
RT <- rbindlist(lapply(names(YK), function(k) {
  dd <- base[!is.na(get(k)) & !is.na(pfas_t)]
  X <- mkX(fml(k, "w5*pfas_t"), dd); cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); cv <- crsV(X, m$residuals, cl); b <- m$coefficients
  nm <- colnames(X)
  jT <- grep("^w5:pfas_tT[23]$|^pfas_tT[23]:w5$", nm)
  eff <- sapply(c("T1","T2","T3"), function(t) {
    j <- which(nm == paste0("w5:pfas_t", t) | nm == paste0("pfas_t", t, ":w5"))
    if (t == "T1") return(100*(exp(b[["w5"]]) - 1))
    if (!length(j)) return(NA_real_)
    100*(exp(b[["w5"]] + b[j]) - 1)
  })
  data.table(outcome = k, label = YK[[k]], N = nrow(dd),
             T1 = eff[1], T2 = eff[2], T3 = eff[3],
             spread = eff[3] - eff[1],
             p_int_joint = pchisq(w2p(b, cv$V, jT), length(jT), lower.tail = FALSE))
}))
RT[, p_bh := p.adjust(p_int_joint, "BH")]
print(RT[, .(label, N, T1 = sprintf("%+6.2f%%", T1), T2 = sprintf("%+6.2f%%", T2),
             T3 = sprintf("%+6.2f%%", T3), spread = sprintf("%+5.2f", spread),
             `联合P` = format.pval(p_int_joint, digits=3))], row.names = FALSE)
fwrite(RT, file.path(RES, "inflam_clean_tertile.csv"))

## ============================================================================
## 3. 处理组内剂量规格（S2，独立复核）
## ============================================================================
cat("\n========== 3. 处理组内剂量规格（S2：w5==1，y ~ dw1*pfas_t）==========\n")
tr <- base[w5 == 1]
RD <- rbindlist(lapply(names(YK), function(k) {
  dd <- tr[!is.na(get(k)) & !is.na(pfas_t)]
  if (nrow(dd) < 200) return(NULL)
  X <- mkX(fml(k, "dw1*pfas_t"), dd); cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); cv <- crsV(X, m$residuals, cl); b <- m$coefficients
  nm <- colnames(X)
  jT <- grep("^dw1:pfas_tT[23]$|^pfas_tT[23]:dw1$", nm)
  j1 <- grep("^dw1:mv$", nm)
  data.table(outcome = k, label = YK[[k]], N = nrow(dd),
             p_int_joint = pchisq(w2p(b, cv$V, jT), length(jT), lower.tail = FALSE))
}))
RD[, p_bh := p.adjust(p_int_joint, "BH")]
print(RD[, .(label, N, p = format.pval(p_int_joint, digits=3))], row.names = FALSE)
fwrite(RD, file.path(RES, "inflam_clean_dose.csv"))

## ============================================================================
## 4. 结论摘要
## ============================================================================
cat("\n================ 结论摘要 ================\n")
sig <- PRIM[p_perm < 0.05 & contrast == "continuous per 1 SD"]
cat("在连续暴露尺度上通过随机化检验（p_perm < 0.05）的结局：\n")
if (nrow(sig)) print(sig[, .(label, coef = signif(coef,3), z_perm = round(z_perm,2),
                             p_perm = signif(p_perm,3))], row.names = FALSE) else cat("  （无）\n")
bad <- PRIM[reliability == "UNRELIABLE"]
cat("\nCR1 名义 P 与随机化 P 相差 >10 倍的组合（名义 P 不可信）：",
    nrow(bad), "/", nrow(PRIM), "\n")
if (nrow(bad)) print(bad[, .(label, contrast, p_cr1 = signif(p_cr1,3),
                             p_perm = signif(p_perm,3))], row.names = FALSE)
cat("\n✅ 36 完成\n")
