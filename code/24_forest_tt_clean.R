## ============================================================================
## 24_forest_tt_clean.R — 睾酮因果森林：药物混杂排除前后对比
##   检验 "PFAS 修饰减重对睾酮的效应" 是否只是 3 名药物去势/抑制者的产物
## ============================================================================
suppressMessages({library(data.table); library(grf)})
set.seed(20261117)
BASE <- Sys.getenv("AGE22_ROOT", ".")
A  <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
FL <- readRDS(file.path(BASE, "data/derived/rx_flags.rds"))
CONF <- A$CONF

mkX <- function(dt, extra) {
  X <- cbind(as.matrix(dt[, ..CONF]), model.matrix(~ race_f + edu_f + smk - 1, data = dt))
  if (length(extra)) X <- cbind(X, as.matrix(dt[, ..extra]))
  mc <- grep("^miss_", names(dt), value = TRUE)
  if (length(mc)) X <- cbind(X, as.matrix(dt[, ..mc]))
  storage.mode(X) <- "double"; X
}

run <- function(dt, label) {
  dt <- copy(dt)[!is.na(log_tt) & !is.na(pfas_z) & !is.na(w5)]
  dt[, mv := as.numeric(pfas_z)]
  X <- mkX(dt, "pfas_z"); Y <- dt$log_tt; W <- as.numeric(dt$w5)
  cl <- as.integer(factor(paste(dt$strata, dt$psu)))
  e0 <- as.numeric(regression_forest(X, W, num.trees = 2000, tune.parameters = "all", seed = 20261117)$predictions)
  keep <- e0 >= 0.05 & e0 <= 0.95
  X <- X[keep, , drop = FALSE]; Y <- Y[keep]; W <- W[keep]; e0 <- e0[keep]
  cl <- cl[keep]; dt <- dt[keep]
  cf <- causal_forest(X, Y, W, W.hat = e0, num.trees = 4000, honesty = TRUE,
                      tune.parameters = "all", clusters = cl, seed = 20261117)
  pr <- predict(cf, estimate.variance = TRUE)
  dt[, `:=`(tau = pr$predictions, se = sqrt(pr$variance.estimates))]
  ate <- average_treatment_effect(cf, target.sample = "overlap")
  cal <- test_calibration(cf)
  blp <- best_linear_projection(cf, A = cbind(mv = dt$mv, mv2 = dt$mv^2))
  cat("\n########", label, "########\n")
  cat("N =", nrow(dt), " 处理 =", sum(W), " 聚类 =", length(unique(cl)), "\n")
  cat(sprintf("ATE(overlap) = %+.4f (SE %.4f) = %+.1f%%\n", ate[1], ate[2], 100*(exp(ate[1])-1)))
  cat("-- test_calibration --\n"); print(round(cal, 4))
  cat("BLP  PFAS 线性  系数 =", round(blp["mv","Estimate"], 4), " p =", format.pval(blp["mv","Pr(>|t|)"], digits=3), "\n")
  cat("BLP  PFAS 二次  系数 =", round(blp["mv2","Estimate"], 4), " p =", format.pval(blp["mv2","Pr(>|t|)"], digits=3), "\n")
  ## 十分位 CATE（加权）
  br <- quantile(dt$mv, seq(0, 1, .1)); br[1] <- -Inf; br[11] <- Inf
  dt[, `:=`(dec = cut(mv, br, labels = 1:10), Wc = W, ec = e0)]
  dq <- dt[, .(n = .N, mv_med = median(mv),
               cate = mean(tau),
               se = mean(se)), by = dec][order(dec)]
  dq[, cate_pct := 100*(exp(cate)-1)]
  print(dq[, .(dec, n, mv_med = round(mv_med,2), pct = round(cate_pct,1),
               lo = round(100*(exp(cate-1.96*se)-1),1),
               hi = round(100*(exp(cate+1.96*se)-1),1))])
  ## 高/低暴露两端的加权平均 CATE
  lo <- dt[mv <= quantile(mv,.2)]; hi <- dt[mv >= quantile(mv,.8)]
  cat(sprintf("最低 20%% 暴露: CATE = %+.1f%%   最高 20%%: CATE = %+.1f%%   差 = %+.1f 个百分点\n",
              100*(exp(mean(lo$tau))-1), 100*(exp(mean(hi$tau))-1),
              100*(exp(mean(hi$tau))-1) - 100*(exp(mean(lo$tau))-1)))
  invisible(list(cal = cal, blp = blp, dec = dq))
}

d3 <- FL[copy(A$S3), on = "SEQN"]
d3[is.na(castration), castration := 0L]; d3[is.na(suppress), suppress := 0L]
d3[, log_tt := log(tt_ngdl)]

cat("=== 样本构成 ===\n")
cat("全样本 N =", nrow(d3), "\n")
cat("药物去势/抗雄/雌激素 =", sum(d3$castration), "\n")
cat("糖皮质激素/阿片/酮康唑 =", sum(d3$suppress), "\n")
cat("两者任一 =", sum(d3$castration==1 | d3$suppress==1), "\n")

r_all <- run(d3, "A. 全样本（原报告口径）")
r_cln <- run(d3[castration==0 & suppress==0], "B. 排除药物影响（清洁样本）")

saveRDS(list(all = r_all, clean = r_cln), file.path(BASE, "results/forest_tt_clean.rds"))
cat("\n✅ 完成\n")
