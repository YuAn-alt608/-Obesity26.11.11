## ============================================================================
## V19_evalue.R — 未测量混杂的定量边界（E-value 敏感性分析）
##                 审稿人第 4 号必问点
## ----------------------------------------------------------------------------
## 问题：§6 承认"横断面回顾性体重，无法完全排除残留混杂"。
##       审稿人会问：需要多强的未测量混杂，才能把唯一存活的 WBC 交互推翻？
##
## 关键便利：本项目全部结局都做了**对数变换**（log CRP / log WBC / log FLI /
##   log SHBG / log ALT），因此 exp(β) 本身就是**比值比**，
##   可直接用 EValue::evalues.RR，无需 Chinn 近似或标准化。
##
## 对交互项的处理（如实说明）：
##   β_int 是"减重者 vs 未减重者"两个对数-RR 之差，exp(β_int) 是**比值比之比（RRR）**。
##   E-value 公式对 RRR 只在"把 RRR 当作 RR 处理"的意义下成立，属近似；
##   因此同时报告两个**分层简单斜率**（各自是严格的 RR）的 E-value，
##   交互的 E-value 不可能低于其组成效应中最弱的一支。
##
## 输出：verification/evalue_table.csv
## ============================================================================
suppressMessages({library(data.table); library(EValue); library(MASS)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
VER  <- file.path(BASE, "verification")
set.seed(20261117)

A <- readRDS(file.path(BASE, "data/derived/analysis.rds")); CONF <- A$CONF
D <- readRDS(file.path(BASE, "data/derived/inflam.rds"))
L <- readRDS(file.path(BASE, "data/derived/liver.rds"))
M <- merge(D, L[, .(SEQN, fli, alt, ok_liver)], by = "SEQN", all.x = TRUE)
M[, `:=`(y_fli = log(pmax(fli, 0.1)), y_alt = log(pmax(alt, .5)),
         mv = as.numeric(pfas_z), dw5 = dw1/5)]
M[, tt_nmol := tt_ngdl / 28.84][, y_shbg := log(shbg)]

crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X)
  bread <- tryCatch(solve(crossprod(X)), error = function(e) ginv(crossprod(X)))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
mkX <- function(f, dd) {
  X <- model.matrix(f, data = dd)
  ok <- apply(X, 2, function(v) all(is.finite(v)))
  vary <- apply(X, 2, function(v) sd(v) > 0)
  X[, ok & (vary | colnames(X) == "(Intercept)"), drop = FALSE]
}

EV <- list()
addEV <- function(finding, spec, est, se, kind = c("RR", "RRR"), note = "") {
  kind <- match.arg(kind)
  rr <- exp(est); lo <- exp(est - 1.96*se); hi <- exp(est + 1.96*se)
  ## evalues.RR 要求 est 为"偏离无效值"的方向；<1 时取倒数，E-value 对称
  e <- tryCatch(evalues.RR(rr, lo, hi), error = function(x) NULL)
  if (is.null(e)) return(NULL)
  EV[[length(EV)+1]] <<- data.table(
    finding = finding, spec = spec, scale = kind,
    est_log = est, se = se,
    RR = rr, RR_lo = lo, RR_hi = hi,
    RR_txt = sprintf("%.3f (%.3f, %.3f)", rr, lo, hi),
    E_value = as.numeric(e["E-values", "point"]),
    ## CI 侧 E-value 取"离无效值最近"的那一端：RR>1 取 lower，RR<1 取 upper
    E_value_CI = as.numeric(if (rr > 1) e["E-values", "lower"] else e["E-values", "upper"]),
    note = note)
}

## ---------------- 1. 主效应（每 5% 减重），未加权与加权 ----------------
SW <- fread(file.path(VER, "survey_weighted_sensitivity.csv"))
KEY <- c(y_crp = "CRP", y_wbc = "WBC", y_fli = "Fatty liver index",
         y_shbg = "SHBG", y_alt = "ALT", y_neut = "Neutrophil")
for (i in seq_len(nrow(SW[term == "main_per5" & outcome %in% names(KEY)]))) {
  r <- SW[term == "main_per5" & outcome %in% names(KEY)][i]
  addEV(sprintf("%s main effect, per 5%% weight loss", KEY[[r$outcome]]),
        "unweighted CR1", r$beta_unw, r$se_unw, "RR",
        sprintf("N=%d, P=%.3g", r$N_unw, r$p_unw))
  addEV(sprintf("%s main effect, per 5%% weight loss", KEY[[r$outcome]]),
        "survey-weighted svyglm", r$beta_wt, r$se_wt, "RR",
        sprintf("P=%.3g", r$p_wt))
}

## ---------------- 2. WBC / 中性粒 交互（重拟合以取简单斜率） ----------------
base <- M[ok_infl == TRUE & !is.na(w5) & !is.na(mv)]
cf <- CONF[sapply(CONF, function(v) length(unique(base[[v]])) > 1)]
for (k in c("y_wbc", "y_neut")) {
  lab <- c(y_wbc = "WBC", y_neut = "Neutrophil")[[k]]
  dd <- base[!is.na(get(k))]
  f <- as.formula(paste(k, "~ w5*mv +", paste(cf, collapse = " + ")))
  X <- mkX(f, dd); m <- lm.fit(X, dd[[k]])
  cl <- factor(paste(dd$strata, dd$psu)); V <- crsV(X, m$residuals, cl)
  nm <- colnames(X); b <- m$coefficients
  jm <- which(nm == "mv"); ji <- grep("^w5:mv$|^mv:w5$", nm)
  ## 交互
  addEV(sprintf("%s x PFAS interaction (continuous DiD)", lab),
        "unweighted CR1", b[ji], sqrt(V[ji, ji]), "RRR",
        sprintf("N=%d；RRR = exp(beta_int)，比值比之比", nrow(dd)))
  ## 简单斜率：未减重者（w5=0）与减重者（w5=1）
  s0 <- b[jm]; s0se <- sqrt(V[jm, jm])
  s1 <- b[jm] + b[ji]
  s1se <- sqrt(V[jm, jm] + V[ji, ji] + 2*V[jm, ji])
  addEV(sprintf("PFAS -> %s, among non-losers (w5 = 0)", lab),
        "simple slope", s0, s0se, "RR", "严格的 RR（结局为对数尺度）")
  addEV(sprintf("PFAS -> %s, among weight losers (w5 = 1)", lab),
        "simple slope", s1, s1se, "RR", "严格的 RR")
}

## ---------------- 3. 书写辅助：E-value 的含义 ------------------------------
E <- rbindlist(EV)
fwrite(E, file.path(VER, "evalue_table.csv"))

cat("================ E-value 表 ================\n")
print(E[, .(finding, spec, RR = RR_txt,
            `E-value` = round(E_value, 2), `E-value(CI)` = round(E_value_CI, 2))],
      row.names = FALSE)

cat("\n================ 一句话结论（可直接写进正文）================\n")
pick <- function(p) E[grepl(p, finding) & grepl("unweighted", spec)][1]
wbc_i <- pick("WBC x PFAS")
wbc_m <- E[grepl("^WBC main", finding) & grepl("unweighted", spec)][1]
crp   <- E[grepl("^CRP main", finding) & grepl("unweighted", spec)][1]
fli   <- E[grepl("^Fatty liver index main", finding) & grepl("unweighted", spec)][1]
shbg  <- E[grepl("^SHBG main", finding) & grepl("unweighted", spec)][1]
cat(sprintf("WBC 交互（RRR = %.3f）：E-value = %.2f（CI 界 %.2f）\n",
            wbc_i$RR, wbc_i$E_value, wbc_i$E_value_CI))
cat(sprintf("WBC 主效应（RR = %.3f）：E-value = %.2f\n", wbc_m$RR, wbc_m$E_value))
cat(sprintf("CRP 主效应（RR = %.3f）：E-value = %.2f（CI 界 %.2f）\n",
            crp$RR, crp$E_value, crp$E_value_CI))
cat(sprintf("FLI 主效应（RR = %.3f）：E-value = %.2f\n", fli$RR, fli$E_value))
cat(sprintf("SHBG 主效应（RR = %.3f）：E-value = %.2f\n", shbg$RR, shbg$E_value))
cat("\n解释模板：要把 WBC 交互完全解释掉，未测量混杂必须与 PFAS 暴露和 WBC 结局\n")
cat("同时达到 RR ≥ ", round(wbc_i$E_value, 2), "（且已调整 12 个协变量之外的额外混杂）。\n", sep = "")
cat("\n✅ V19 完成 → verification/evalue_table.csv\n")
