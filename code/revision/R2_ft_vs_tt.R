## ============================================================================
## R2_ft_vs_tt.R — 回应审稿意见 #4（临床靶子应为 fT 而非 %free；缺 fT−TT 差值检验）
##   把 TT 与 fT 堆叠进同一个模型（受试者内两行），聚类稳健协方差在 PSU 层面，
##   直接检验 H0: β(fT) = β(TT)，而不是拿"一个显著、一个不显著"当作差异存在的证据。
##   同时检验 β(SHBG) vs β(TT)（机制：SHBG 涨得比 TT 快）。
##   并做 PSU 聚类自助法（1000 次）交叉验证。
## 输出: results/revision/R2_ft_tt_difference.csv, logs/revision_R2.log
## ============================================================================
suppressMessages({library(data.table)})
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

d <- FL[copy(A$S4), on = "SEQN"]
d[is.na(castration), castration := 0L]; d[is.na(suppress), suppress := 0L]
d <- d[castration == 0 & suppress == 0]                       # 主口径：清洁样本
d[, `:=`(tt_nmol = tt_ngdl / 28.84, e2_pmol = e2_pgml * 3.671)]
d[, `:=`(y_tt = log(tt_nmol), y_shbg = log(shbg), y_ft = log(ft_nmol),
         y_pctfree = log(pct_free))]
YK <- c("y_tt", "y_ft", "y_shbg", "y_pctfree")
cf <- CONF[sapply(CONF, function(v) length(unique(d[[v]])) > 1)]
cc <- complete.cases(d[, c("w5", "pfas_z", "dw1", YK, cf, "race_f", "edu_f", "smk"), with = FALSE])
d  <- d[cc]
d[, `:=`(mv = as.numeric(pfas_z), dw5 = dw1/5)]
d[, mv2 := mv^2]
cat("清洁样本 N =", nrow(d), " 处理组 =", sum(d$w5), " 聚类 =", uniqueN(paste(d$strata, d$psu)), "\n\n")

## ---------- 堆叠设计 --------------------------------------------------------
Z <- model.matrix(as.formula(paste("~ dw5*(mv+mv2) + mv + mv2 +", paste(cf, collapse = " + "),
                                   "+ race_f + edu_f + smk")), data = d)
n <- nrow(d); p <- ncol(Z)
Xstack <- rbind(cbind(Z, matrix(0, n, p)), cbind(matrix(0, n, p), Z))
Ystack <- c(d$y_tt, d$y_ft)
cl <- factor(rep(paste(d$strata, d$psu), 2))
m <- lm.fit(Xstack, Ystack); V <- crsV(Xstack, m$residuals, cl)
b <- m$coefficients
iw <- grep("^dw5$", colnames(Z))
bt <- b[iw]; bf <- b[p + iw]; bs <- NULL

## 差值检验（delta 法）
L <- rep(0, 2 * p); L[p + iw] <- 1; L[iw] <- -1
diff_est <- as.numeric(L %*% b); diff_se <- sqrt(as.numeric(t(L) %*% V %*% L))
zt <- diff_est / diff_se
cat("=== fT vs TT：每 5% 减重的对数差 ===\n")
cat(sprintf("  β(TT) = %+.4f  → %+.2f%%\n", bt, 100 * (exp(bt) - 1)))
cat(sprintf("  β(fT) = %+.4f  → %+.2f%%\n", bf, 100 * (exp(bf) - 1)))
cat(sprintf("  差值  = %+.4f  → %+.2f%%  (SE %.4f, z = %.2f, 双侧 P = %.3g)\n",
            diff_est, 100 * (exp(diff_est) - 1), diff_se, zt, pnorm(-abs(zt)) * 2))
cat(sprintf("  同名恒等式检查：β(%%free) = %+.4f  vs  β(fT)-β(TT) = %+.4f\n\n",
            b[p + iw] - bt, diff_est))

## TT vs SHBG
Ystack2 <- c(d$y_tt, d$y_shbg)
m2 <- lm.fit(Xstack, Ystack2); V2 <- crsV(Xstack, m2$residuals, cl); b2 <- m2$coefficients
ds <- b2[p + iw] - b2[iw]; ds_se <- sqrt(as.numeric(t(L) %*% V2 %*% L))
cat("=== SHBG vs TT ===\n")
cat(sprintf("  β(SHBG) = %+.4f → %+.2f%% ; 差值 %+.4f (%+.2f%%), SE %.4f, P = %.3g\n\n",
            b2[p + iw], 100 * (exp(b2[p + iw]) - 1), ds, 100 * (exp(ds) - 1), ds_se,
            2 * pnorm(-abs(ds / ds_se))))

## ---------- PSU 聚类自助法 --------------------------------------------------
set.seed(20261117)
ps <- unique(cl); B <- 1000
boot <- matrix(NA_real_, B, 3)
for (k in seq_len(B)) {
  bs <- sample(ps, length(ps), replace = TRUE)
  idx <- unlist(lapply(bs, function(p_) which(cl == p_)), use.names = FALSE)
  ## 单人两行必须同时进入
  idx1 <- idx[idx <= n]; idx2 <- idx[idx > n]
  Xi <- rbind(Xstack[idx1, , drop = FALSE], Xstack[idx2, , drop = FALSE])
  Yi <- c(d$y_tt[idx1], d$y_ft[idx2 - n])
  mm <- tryCatch(lm.fit(Xi, Yi), error = function(e) NULL)
  if (is.null(mm)) next
  boot[k, ] <- c(mm$coefficients[iw], mm$coefficients[p + iw],
                 mm$coefficients[p + iw] - mm$coefficients[iw])
}
boot <- boot[complete.cases(boot), , drop = FALSE]
cat(sprintf("=== PSU 自助法（%d 次有效）===\n", nrow(boot)))
cat(sprintf("  β(fT)-β(TT) 自助 95%%CI: %+.4f – %+.4f  (点估计 %+.4f)\n",
            quantile(boot[,3], .025), quantile(boot[,3], .975), diff_est))
cat(sprintf("  自助 P(差值≥0) = %.4f\n", mean(boot[,3] >= 0)))

## ---------- 绝对尺度 & 交互差分 ---------------------------------------------
## %free 本身的模型（与差值同源，但作为独立报告项）
Z2 <- Z
m3 <- lm.fit(Z, d$y_pctfree); V3 <- crsV(Z, m3$residuals, factor(paste(d$strata, d$psu)))
b3 <- m3$coefficients
cat(sprintf("\nβ(%%free) = %+.4f → %+.2f%% (P = %.3g)\n",
            b3[iw], 100 * (exp(b3[iw]) - 1), 2 * pnorm(-abs(b3[iw] / sqrt(V3[iw, iw])))))

RES <- data.table(
  contrast = c("TT", "fT", "%free", "SHBG", "fT - TT", "SHBG - TT"),
  beta = c(bt, bf, b3[iw], b2[p + iw], diff_est, ds),
  pct = 100 * (exp(c(bt, bf, b3[iw], b2[p + iw], diff_est, ds)) - 1),
  se = c(sqrt(V[iw, iw]), sqrt(V[p + iw, p + iw]), sqrt(V3[iw, iw]),
         sqrt(V2[p + iw, p + iw]), diff_se, ds_se),
  p = c(2*pnorm(-abs(bt/sqrt(V[iw,iw]))), 2*pnorm(-abs(bf/sqrt(V[p+iw,p+iw]))),
        2*pnorm(-abs(b3[iw]/sqrt(V3[iw,iw]))), 2*pnorm(-abs(b2[p+iw]/sqrt(V2[p+iw,p+iw]))),
        2*pnorm(-abs(zt)), 2*pnorm(-abs(ds/ds_se))),
  N = nrow(d), treated = sum(d$w5))
print(RES)
fwrite(RES, file.path(OUT, "R2_ft_tt_difference.csv"))
fwrite(data.table(diff = boot[,3]), file.path(OUT, "R2_bootstrap_diff.csv"))
cat("✅ R2 完成\n")
