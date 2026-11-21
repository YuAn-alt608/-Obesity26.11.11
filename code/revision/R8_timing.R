## ============================================================================
## R8_timing.R — 回应"抽血时点"：减重发生得早或晚，能否解释 R6 的零结果
##   ① 轨迹分型（dw10 x dw1）作为时间代理，参照 = 十年与一年都未减重
##   ② 意图分层（WHQ060）
##   ③ 参照组污染：现行参照（当前未减重）vs 干净参照（十年未减重）
##   ④ 药代上界：一次释放脉冲在 1–36 个月后的残留分数
## 输出: results/revision/R8_trajectory.csv, R8_pk_decay.csv
## ============================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", ".")
A  <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
I  <- readRDS(file.path(BASE, "data/derived/inflam.rds"))
OUT <- file.path(BASE, "results/revision"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CONF <- A$CONF

crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X); bread <- solve(crossprod(X))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}

## term: 单一因子变量名（每组一个哑变量）
fitY <- function(dd, y, term, covs, label) {
  vars <- unique(c(y, term, covs, "race_f", "edu_f", "smk", "strata", "psu"))
  dd <- dd[complete.cases(dd[, ..vars])]
  f <- as.formula(paste(y, "~", paste(c(term, covs, "race_f", "edu_f", "smk"), collapse = " + ")))
  X <- model.matrix(f, data = dd)
  keep <- apply(X, 2, function(z) sd(z) > 0) | colnames(X) == "(Intercept)"
  X <- X[, keep, drop = FALSE]
  m <- lm.fit(X, dd[[y]]); V <- crsV(X, m$residuals, factor(paste(dd$strata, dd$psu)))
  b <- m$coefficients; nm <- colnames(X)
  rr <- rbindlist(lapply(grep(paste0("^", term), nm, value = TRUE), function(j) {
    i <- match(j, nm); se <- sqrt(V[i, i])
    data.table(outcome = label, group = sub(paste0("^", term), "", j),
               N_model = nrow(dd), n_group = sum(dd[[term]] == sub(paste0("^", term), "", j)),
               pct = 100 * (exp(b[i]) - 1),
               lo = 100 * (exp(b[i] - 1.96 * se) - 1),
               hi = 100 * (exp(b[i] + 1.96 * se) - 1),
               p = 2 * pnorm(-abs(b[i] / se)))
  }))
  rr
}

d <- I[!is.na(log_pfas) & !is.na(w5) & !is.na(dw1)]
cat("分析域 N =", nrow(d), "; dw10 非缺 =", sum(!is.na(d$dw10)), "\n")
cat("符号核对：dw1>0 = 体重下降；w5==1 时 dw1 中位 =", median(d[w5 == 1, dw1]), "\n\n")

## ---------- ① 轨迹分型 ------------------------------------------------------
d[, traj := fcase(
  !is.na(dw10) & abs(dw1) < 5 & abs(dw10) < 5, "never_lost",
  !is.na(dw10) & abs(dw1) < 5 & dw10 >= 5,     "long_past_loss_now_stable",
  !is.na(dw10) & dw1 >= 5 & abs(dw10) < 5,     "recent_loss",
  !is.na(dw10) & dw1 >= 5 & dw10 >= 5,         "sustained_loss",
  !is.na(dw10) & dw1 >= 5 & dw10 <= -5,        "loss_after_gain",
  default = NA_character_)]
lev <- c("never_lost", "long_past_loss_now_stable", "recent_loss",
         "sustained_loss", "loss_after_gain")
d[, traj := factor(traj, levels = lev)]
cat("=== ① 轨迹分型 ===\n")
print(d[!is.na(traj), .(N = .N, geo = round(exp(mean(log_pfas)), 1),
                        med = round(median(pfas_molar), 1)), by = traj][order(traj)])
cat("\n调整后（参照 = never_lost，n =", d[traj == "never_lost", .N], "）\n")
R1 <- fitY(d[!is.na(traj)], "log_pfas", "traj", CONF, "trajectory")
print(R1[, .(group, n_group, pct = sprintf("%+.2f%%", pct),
             ci = sprintf("%+.2f, %+.2f", lo, hi), p = signif(p, 3))])
fwrite(R1, file.path(OUT, "R8_trajectory.csv"))

cat("\n=== ①b 与 R6 同口径（w5 二分，参照 = 全样本未减重）===\n")
d[, bin := factor(fifelse(w5 == 1, "lost5", "not_lost"), levels = c("not_lost", "lost5"))]
R1b <- fitY(d, "log_pfas", "bin", CONF, "w5 binary")
print(R1b[, .(group, n_group, pct = sprintf("%+.2f%%", pct),
              ci = sprintf("%+.2f, %+.2f", lo, hi), p = signif(p, 3))])

## ---------- ② 意图分层 ------------------------------------------------------
cat("\n=== ② 意图（WHQ060；参照 = 全样本未减重）===\n")
for (v in c("intentional", "unintentional")) {
  code <- if (v == "intentional") 1 else 2
  dd <- d[w5 == 0 | (w5 == 1 & !is.na(wt_intent) & wt_intent == code)]
  dd[, trt := factor(fifelse(w5 == 1, "lost5", "not_lost"), levels = c("not_lost", "lost5"))]
  r <- fitY(dd, "log_pfas", "trt", CONF, v)
  cat(sprintf("  %-14s", v))
  print(r[group == "lost5", .(n_lost = n_group, pct = round(pct, 2),
                              ci = sprintf("%+.2f, %+.2f", lo, hi), p = signif(p, 3))])
}

## ---------- ③ 参照组污染 ----------------------------------------------------
cat("\n=== ③ 参照组：当前未减重（现行）vs 十年与一年都未减重 ===\n")
d[, refset := fcase(!is.na(dw10) & abs(dw10) < 5 & abs(dw1) < 5, "never_lost_ref",
                    !is.na(dw10) & dw1 >= 5, "lost5_clean",
                    default = NA_character_)]
d[, refset := factor(refset, levels = c("never_lost_ref", "lost5_clean"))]
R3 <- fitY(d[!is.na(refset)], "log_pfas", "refset", CONF, "clean reference")
print(R3[, .(group, n_group, pct = sprintf("%+.2f%%", pct),
             ci = sprintf("%+.2f, %+.2f", lo, hi), p = signif(p, 3))])

## ---------- ④ 药代上界 ------------------------------------------------------
cat("\n=== ④ 一次释放脉冲的残留分数（按报告半衰期，年）===\n")
hl <- data.table(congener = c("PFOA", "PFOS", "PFHxS", "PFNA", "PFDA", "PFUnDA"),
                 t_half = c(3.5, 5.4, 8.5, 3.5, 6.0, 7.0))
mon <- c(1, 3, 6, 12, 18, 24, 36)
PK <- copy(hl)
for (m in mon) PK[, (paste0("m", m)) := round(100 * 0.5^(m/12/t_half), 1)]
print(PK)
fwrite(PK, file.path(OUT, "R8_pk_decay.csv"))
cat("\n敏感性（PFOA，12 个月残留）：")
for (h in c(1.75, 3.5, 7.0)) cat(sprintf("  t1/2=%.2fy -> %.0f%%;", h, 100*0.5^(1/h)))
cat("\n")
cat("✅ R8 完成\n")
