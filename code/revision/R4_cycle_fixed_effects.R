## ============================================================================
## R4_cycle_fixed_effects.R — 回应审稿意见 #3/#5（周期只有 2–3 个、yr 只是线性项）
##   现行模型把周期当作单一连续变量 yr 进入（线性趋势）。本脚本把 yr 换成
##   周期固定效应（各周期哑变量），检验：
##     ① PhenoAge 主效应（粗 +1.87 → 调整后 −0.33）是否还站得住
##     ② PhenoAge 森林 ATE(overlap) 是否还站得住
##     ③ 雄激素每 5% 主效应（TT / SHBG / fT / %free）是否还站得住
## 输出: results/revision/R4_cycle_fe.csv
## ============================================================================
suppressMessages({library(data.table); library(grf)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
A  <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
FL <- readRDS(file.path(BASE, "data/derived/rx_flags.rds"))
OUT <- file.path(BASE, "results/revision"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CONF <- A$CONF

crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X); bread <- solve(crossprod(X))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
rows <- list()
add <- function(...) rows[[length(rows)+1]] <<- data.table(...)

## ---------- ① PhenoAge 参数模型 ---------------------------------------------
S1 <- copy(A$S1)
S1 <- S1[!is.na(phenoage_accel) & !is.na(w5)]
S1[, cyc := factor(cycle)]
CONF_N <- setdiff(CONF, "yr")
CONF_S1 <- CONF[sapply(CONF, function(v) length(unique(S1[[v]])) > 1)]
CONF_S1N <- setdiff(CONF_S1, "yr")

## 粗
m <- lm(phenoage_accel ~ w5, data = S1)
V <- crsV(model.matrix(m), residuals(m), factor(paste(S1$strata, S1$psu)))
add(model = "PhenoAge", spec = "crude (no covariates)", N = nrow(S1),
    est = coef(m)["w5"], se = sqrt(V["w5","w5"]),
    p = 2*pnorm(-abs(coef(m)["w5"]/sqrt(V["w5","w5"]))))
## 现行：yr 线性
f1 <- as.formula(paste("phenoage_accel ~ w5 +", paste(CONF_S1, collapse=" + "), "+ race_f + edu_f + smk"))
m1 <- lm(f1, data = S1); X1 <- model.matrix(m1)
V1 <- crsV(X1, residuals(m1), factor(paste(S1$strata, S1$psu)))
add(model = "PhenoAge", spec = "adjusted, year entered linearly (current)", N = nrow(S1),
    est = coef(m1)["w5"], se = sqrt(V1["w5","w5"]),
    p = 2*pnorm(-abs(coef(m1)["w5"]/sqrt(V1["w5","w5"]))))
## 周期固定效应
f2 <- as.formula(paste("phenoage_accel ~ w5 +", paste(CONF_S1N, collapse=" + "), "+ race_f + edu_f + smk + cyc"))
m2 <- lm(f2, data = S1); X2 <- model.matrix(m2)
V2 <- crsV(X2, residuals(m2), factor(paste(S1$strata, S1$psu)))
add(model = "PhenoAge", spec = "adjusted, cycle fixed effects", N = nrow(S1),
    est = coef(m2)["w5"], se = sqrt(V2["w5","w5"]),
    p = 2*pnorm(-abs(coef(m2)["w5"]/sqrt(V2["w5","w5"]))))

## ---------- ② PhenoAge 森林：yr 线性 vs 周期固定效应 -------------------------
mkX <- function(dt, cycle_fe = FALSE) {
  cc <- if (cycle_fe) CONF_N else CONF
  X <- cbind(as.matrix(dt[, ..cc]))
  if (cycle_fe) {
    X <- cbind(X, model.matrix(~ race_f + edu_f + smk + cycle - 1, data = dt))
  } else {
    X <- cbind(X, model.matrix(~ race_f + edu_f + smk - 1, data = dt))
  }
  X <- cbind(X, as.matrix(dt[, .(pfas_z)]))
  mc <- grep("^miss_", names(dt), value = TRUE)
  if (length(mc)) X <- cbind(X, as.matrix(dt[, ..mc]))
  storage.mode(X) <- "double"; X
}
run_cf_ate <- function(dt, cycle_fe, tag) {
  dt <- copy(dt)[!is.na(phenoage_accel) & !is.na(pfas_z) & !is.na(w5)]
  dt[, mv := as.numeric(pfas_z)]
  X <- mkX(dt, cycle_fe); Y <- dt$phenoage_accel; W <- as.numeric(dt$w5)
  cl <- as.integer(factor(paste(dt$strata, dt$psu)))
  e0 <- as.numeric(regression_forest(X, W, num.trees = 2000, tune.parameters = "all",
                                     seed = 20261117)$predictions)
  k <- e0 >= 0.05 & e0 <= 0.95
  cf <- causal_forest(X[k, ], Y[k], W[k], W.hat = e0[k], num.trees = 4000, honesty = TRUE,
                      tune.parameters = "all", clusters = cl[k], seed = 20261117)
  ate <- average_treatment_effect(cf, target.sample = "overlap")
  cat(sprintf("  %-34s N=%5d  ATE = %+.3f (SE %.3f)  p = %.3g\n", tag, sum(k), ate[1], ate[2],
              2*pnorm(-abs(ate[1]/ate[2]))))
  add(model = "PhenoAge forest", spec = tag, N = sum(k), est = ate[1], se = ate[2],
      p = 2*pnorm(-abs(ate[1]/ate[2])))
}
cat("=== PhenoAge 森林 ===\n")
run_cf_ate(A$S1, FALSE, "AIPW forest, year linear (current)")
run_cf_ate(A$S1, TRUE,  "AIPW forest, cycle fixed effects")

## ---------- ③ 雄激素：周期固定效应 ------------------------------------------
d <- FL[copy(A$S4), on = "SEQN"]
d[is.na(castration), castration := 0L]; d[is.na(suppress), suppress := 0L]
d <- d[castration == 0 & suppress == 0]
d[, `:=`(tt_nmol = tt_ngdl/28.84, e2_pmol = e2_pgml*3.671)]
d[, `:=`(y_tt = log(tt_nmol), y_shbg = log(shbg), y_ft = log(ft_nmol), y_pctfree = log(pct_free))]
YK <- c("y_tt","y_shbg","y_ft","y_pctfree")
cc <- complete.cases(d[, c("w5","pfas_z","dw1",YK,CONF,"race_f","edu_f","smk"), with = FALSE])
d <- d[cc]; d[, `:=`(mv = as.numeric(pfas_z), dw5 = dw1/5)]; d[, mv2 := mv^2]; d[, cyc := factor(cycle)]
CONF_A  <- CONF[sapply(CONF, function(v) length(unique(d[[v]])) > 1)]
CONF_AN <- setdiff(CONF_A, "yr")
cat("\n=== 雄激素（清洁样本 N =", nrow(d), ", 周期 =", paste(levels(d$cyc), collapse="/"), "）===\n")
cat("周期分布:\n"); print(d[, .N, by=cycle])
for (k in YK) {
  ## 现行
  fa <- as.formula(paste(k, "~ dw5*(mv+mv2) + mv + mv2 +", paste(CONF_A, collapse=" + "), "+ race_f + edu_f + smk"))
  ma <- lm(fa, data = d); Xa <- model.matrix(ma)
  Va <- crsV(Xa, residuals(ma), factor(paste(d$strata, d$psu)))
  ## 周期固定效应
  fb <- as.formula(paste(k, "~ dw5*(mv+mv2) + mv + mv2 +", paste(CONF_AN, collapse=" + "), "+ race_f + edu_f + smk + cyc"))
  mb <- lm(fb, data = d); Xb <- model.matrix(mb)
  Vb <- crsV(Xb, residuals(mb), factor(paste(d$strata, d$psu)))
  ia <- grep("^dw5$", colnames(Xa)); ib <- grep("^dw5$", colnames(Xb))
  j1 <- grep("^dw5:mv$", colnames(Xa)); j2 <- grep("^dw5:mv2$", colnames(Xa))
  k1 <- grep("^dw5:mv$", colnames(Xb)); k2 <- grep("^dw5:mv2$", colnames(Xb))
  wald <- function(b, V) as.numeric(t(b) %*% solve(V) %*% b)
  cat(sprintf("  %-9s 现行 %+.2f%% (p=%.2g) | 周期FE %+.2f%% (p=%.2g) | 交互P 现行 %.3f → FE %.3f\n",
              k, 100*(exp(coef(ma)[ia])-1), 2*pnorm(-abs(coef(ma)[ia]/sqrt(Va[ia,ia]))),
              100*(exp(coef(mb)[ib])-1), 2*pnorm(-abs(coef(mb)[ib]/sqrt(Vb[ib,ib]))),
              pchisq(wald(coef(ma)[c(j1,j2)], Va[c(j1,j2),c(j1,j2)]), 2, lower.tail=FALSE),
              pchisq(wald(coef(mb)[c(k1,k2)], Vb[c(k1,k2),c(k1,k2)]), 2, lower.tail=FALSE)))
  add(model = k, spec = "per 5%, year linear (current)", N = nrow(d),
      est = 100*(exp(coef(ma)[ia])-1), se = 100*sqrt(Va[ia,ia]),
      p = 2*pnorm(-abs(coef(ma)[ia]/sqrt(Va[ia,ia]))))
  add(model = k, spec = "per 5%, cycle fixed effects", N = nrow(d),
      est = 100*(exp(coef(mb)[ib])-1), se = 100*sqrt(Vb[ib,ib]),
      p = 2*pnorm(-abs(coef(mb)[ib]/sqrt(Vb[ib,ib]))))
}

RES <- rbindlist(rows)
fwrite(RES, file.path(OUT, "R4_cycle_fe.csv"))
cat("\n"); print(RES[, .(model, spec, N, est = round(est, 3), se = round(se, 3),
                        p = format.pval(p, digits = 3))])
cat("✅ R4 完成\n")
