## ============================================================================
## 10_causal_forest.R — 广义随机森林 / 因果森林估计 CATE
##   处理: 1 年内体重下降 >= 5%   (W)
##   结局: PhenoAge 加速(年) / log 总睾酮 / log 游离睾酮
##   效应修饰: ΣPFAS (8周期) / ΣPCB / ΣOCP (2003-04)
## ============================================================================
suppressMessages({library(data.table); library(grf); library(mgcv); library(splines)})
set.seed(20261117)
BASE <- Sys.getenv("AGE22_ROOT", ".")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
CONF <- A$CONF
mkX <- function(dt, extra) {
  X <- cbind(as.matrix(dt[, ..CONF]), model.matrix(~ race_f + edu_f + smk - 1, data = dt))
  if (length(extra)) X <- cbind(X, as.matrix(dt[, ..extra]))
  mc <- grep("^miss_", names(dt), value = TRUE)
  if (length(mc)) X <- cbind(X, as.matrix(dt[, ..mc]))
  storage.mode(X) <- "double"; X
}

run_cf <- function(dt, yvar, modvar, label, ntree = 4000, cluster_boot = 200) {
  dt <- copy(dt)[!is.na(get(yvar)) & !is.na(get(modvar)) & !is.na(w5)]
  dt[, mv := as.numeric(get(modvar))]
  X <- mkX(dt, modvar); Y <- dt[[yvar]]; W <- as.numeric(dt$w5)
  cl <- as.integer(factor(paste(dt$strata, dt$psu)))  # 层-半PSU 为聚类单位
  cat("\n################", label, "################\n")
  cat("N =", nrow(dt), " 处理组 =", sum(W), " 修饰因子 SD =", round(sd(dt$mv), 3), "\n")

  ## 先估计倾向得分, 修剪至重叠区 (避免极端倾向得分导致的 AIPW 不稳)
  e0 <- as.numeric(regression_forest(X, W, num.trees = 2000, tune.parameters = "all", seed = 20261117)$predictions)
  keep <- e0 >= 0.05 & e0 <= 0.95
  cat("重叠区保留 ", sum(keep), "/", length(keep), " (", round(mean(keep) * 100, 1), "% )\n", sep = "")
  X <- X[keep, , drop = FALSE]; Y <- Y[keep]; W <- W[keep]; e0 <- e0[keep]
  cl <- cl[keep]; dt <- dt[keep]
  cat("修剪后 N =", nrow(dt), " 处理组 =", sum(W), "\n")
  cf <- causal_forest(X, Y, W, Y.hat = NULL, W.hat = e0, num.trees = ntree, honesty = TRUE,
                      tune.parameters = "all", clusters = cl, seed = 20261117)
  pr <- predict(cf, estimate.variance = TRUE)
  dt[, `:=`(tau = pr$predictions, tau_se = sqrt(pr$variance.estimates))]
  cat("倾向得分范围:", round(range(cf$W.hat), 3), " 1/99分位:",
      round(quantile(cf$W.hat, c(.01, .99)), 3), "\n")

  ate <- average_treatment_effect(cf, target.sample = "overlap")
  ate_all <- average_treatment_effect(cf, target.sample = "all")
  cat(sprintf("ATE(overlap) = %.3f (SE %.3f), p = %s ;  ATE(all) = %.3f\n", ate[1], ate[2],
               format.pval(2 * pnorm(-abs(ate[1] / ate[2])), digits = 3), ate_all[1]))
  cal <- test_calibration(cf)
  cat("\n-- 校准检验 --\n"); print(round(cal, 4))

  ## BLP: 二次型 (直接检验 U/拐点)
  blp <- best_linear_projection(cf, A = cbind(mv = dt$mv, mv2 = dt$mv^2))
  cat("\n-- BLP (线性 + 二次) --\n"); print(round(blp, 4))

  ## ---- AIPW 伪结局 + 平滑曲线 ----
  rf1 <- regression_forest(X[W == 1, , drop = FALSE], Y[W == 1], num.trees = 2000, tune.parameters = "all")
  rf0 <- regression_forest(X[W == 0, , drop = FALSE], Y[W == 0], num.trees = 2000, tune.parameters = "all")
  mu1 <- as.numeric(predict(rf1, X)$predictions); mu0 <- as.numeric(predict(rf0, X)$predictions)
  eh  <- pmin(pmax(e0, 0.02), 0.98)
  dt[, aipw := mu1 - mu0 + W * (Y - mu1) / eh - (1 - W) * (Y - mu0) / (1 - eh)]

  g <- gam(aipw ~ s(mv, k = 5, bs = "tp"), data = dt, method = "REML")
  grid <- data.table(mv = seq(quantile(dt$mv, .05), quantile(dt$mv, .95), length.out = 80))
  lpm <- predict(g, newdata = grid, type = "lpmatrix")
  V <- vcov(g); fit <- as.numeric(lpm %*% coef(g)); se <- sqrt(rowSums((lpm %*% V) * lpm))
  curve <- data.table(mod = grid$mv, cate = fit, se = se,
                      lo = fit - 1.96 * se, hi = fit + 1.96 * se)

  cat("\n-- 曲线形状 --\n")
  cat(sprintf("有效自由度 = %.2f, 非线性 F = %.2f, p = %s\n", sum(g$edf),
              g$s.table[1, 3], format.pval(g$s.table[1, 4], digits = 3)))
  ## 高暴露 vs 低暴露 CATE 差异 (GAM 线性组合)
  ii <- c(which.min(abs(grid$mv - quantile(dt$mv, .1))), which.min(abs(grid$mv - quantile(dt$mv, .9))))
  L <- lpm[ii[2], ] - lpm[ii[1], ]
  d_est <- sum(L * coef(g)); d_se <- sqrt(as.numeric(t(L) %*% V %*% L))
  cat(sprintf("CATE(p90) - CATE(p10) = %.3f (SE %.3f), p = %s\n", d_est, d_se,
              format.pval(2 * pnorm(-abs(d_est / d_se)), digits = 3)))

  ## 二次项 (AIPW)
  q <- lm(aipw ~ mv + I(mv^2), data = dt)
  b1 <- coef(q)[2]; b2 <- coef(q)[3]
  cat(sprintf("二次项 b2 = %.4f (SE %.4f), p = %s, 拐点 = %.2f SD\n", b2,
              sqrt(vcov(q)[3, 3]), format.pval(summary(q)$coef[3, 4], digits = 3),
              ifelse(b2 != 0, -b1 / (2 * b2), NA)))

  ## 十分位
  dt[, dec := cut(mv, unique(quantile(mv, 0:10/10)), include.lowest = TRUE, labels = 1:10)]
  bydec <- dt[, .(N = .N, treated = sum(w5), mod_med = round(median(mv), 2),
                  tau_cf = mean(tau), aipw = mean(aipw),
                  aipw_se = sd(aipw) / sqrt(.N)), by = dec][order(dec)]
  print(round(bydec[, 3:6], 3))

  ## 聚类自助法 CI (重抽 PSU, 固定 AIPW 伪结局)
  ps <- unique(cl)
  bidx <- lapply(seq_len(cluster_boot), function(b) sample(ps, length(ps), replace = TRUE))
  bm <- sapply(bidx, function(bs) {
    id <- unlist(lapply(bs, function(p) which(cl == p)))
    gg <- tryCatch(gam(aipw ~ s(mv, k = 5, bs = "tp"), data = dt[id], method = "REML"),
                   error = function(e) NULL)
    if (is.null(gg)) return(rep(NA_real_, nrow(grid)))
    predict(gg, newdata = grid)
  })
  curve[, `:=`(boot_lo = apply(bm, 1, quantile, .025, na.rm = TRUE),
               boot_hi = apply(bm, 1, quantile, .975, na.rm = TRUE))]

  list(cf = cf, data = dt, bydec = bydec, curve = curve, ate = ate, ate_all = ate_all, cal = cal, blp = blp,
       edf = sum(g$edf), nf = g$s.table[1, 3], nf_p = g$s.table[1, 4],
       d_est = d_est, d_se = d_se, d_p = 2 * pnorm(-abs(d_est / d_se)),
       b2 = b2, b2_p = summary(q)$coef[3, 4], turn = ifelse(b2 != 0, -b1 / (2 * b2), NA),
       modvar = modvar, yvar = yvar, label = label)
}

res <- list()
res$S1  <- run_cf(A$S1,  "phenoage_accel", "pfas_z", "S1  ΣPFAS × PhenoAge加速 (8周期)")
res$S2  <- run_cf(A$S2,  "phenoage_accel", "pcb_z",  "S2  ΣPCB × PhenoAge加速 (2003-04)")
res$S2b <- run_cf(A$S2b, "phenoage_accel", "ocp_z",  "S2b ΣOCP × PhenoAge加速 (2003-04)")
A$S3[, log_tt := log(tt_ngdl)]; A$S4[, log_ft := log(ft_nmol)]
res$S3  <- run_cf(A$S3, "log_tt", "pfas_z", "S3  ΣPFAS × log总睾酮 (男性)")
res$S4  <- run_cf(A$S4, "log_ft", "pfas_z", "S4  ΣPFAS × log游离睾酮 (男性)")

saveRDS(res, file.path(BASE, "results/cf_results.rds"))
sm <- rbindlist(lapply(names(res), function(k) {
  r <- res[[k]]
  data.table(model = k, label = r$label, N = nrow(r$data), treated = sum(r$data$w5),
             ATE = round(r$ate[1], 3), ATE_se = round(r$ate[2], 3),
             edf = round(r$edf, 2), nonlin_p = signif(r$nf_p, 3),
             d_hi_lo = round(r$d_est, 3), d_p = signif(r$d_p, 3),
             quad_b2 = round(r$b2, 4), quad_p = signif(r$b2_p, 3), turn = round(r$turn, 2))
}))
print(sm); fwrite(sm, file.path(BASE, "results/cf_summary.csv"))
cat("\nALL DONE\n")
