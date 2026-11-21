## ============================================================================
## 25_seed_stability.R — 森林异质性检验的随机种子稳定性
##   同一份数据、同一套超参数，只换随机种子，看 P 值波动范围
##   若 P 值在种子间大幅漂移 => 不能作为"显著异质性"的证据
## ============================================================================
suppressMessages({library(data.table); library(grf)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
A  <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
FL <- readRDS(file.path(BASE, "data/derived/rx_flags.rds"))
CONF <- A$CONF

mkX <- function(dt) {
  X <- cbind(as.matrix(dt[, ..CONF]), model.matrix(~ race_f + edu_f + smk - 1, data = dt))
  X <- cbind(X, as.matrix(dt[, .(pfas_z)]))
  mc <- grep("^miss_", names(dt), value = TRUE)
  if (length(mc)) X <- cbind(X, as.matrix(dt[, ..mc]))
  storage.mode(X) <- "double"; X
}

one <- function(d3, seed, ntree = 3000) {
  set.seed(seed)
  dt <- copy(d3)[!is.na(log_tt) & !is.na(pfas_z) & !is.na(w5)]
  dt[, mv := as.numeric(pfas_z)]
  X <- mkX(dt); Y <- dt$log_tt; W <- as.numeric(dt$w5)
  cl <- as.integer(factor(paste(dt$strata, dt$psu)))
  e0 <- as.numeric(regression_forest(X, W, num.trees = 1500, tune.parameters = "all", seed = seed)$predictions)
  keep <- e0 >= 0.05 & e0 <= 0.95
  X <- X[keep, ]; Y <- Y[keep]; W <- W[keep]; e0 <- e0[keep]; cl <- cl[keep]; dt <- dt[keep]
  cf <- causal_forest(X, Y, W, W.hat = e0, num.trees = ntree, honesty = TRUE,
                      tune.parameters = "all", clusters = cl, seed = seed)
  cal <- test_calibration(cf)
  blp <- best_linear_projection(cf, A = cbind(mv = dt$mv))
  tau <- predict(cf)$predictions
  lo <- mean(tau[dt$mv <= quantile(dt$mv, .2)]); hi <- mean(tau[dt$mv >= quantile(dt$mv, .8)])
  data.table(seed = seed, N = nrow(dt),
             p_hetero = cal["differential.forest.prediction", "Pr(>t)"],
             p_hetero_2s = 2 * cal["differential.forest.prediction", "Pr(>t)"],
             blp_mv = blp["mv", "Estimate"], p_blp = blp["mv", "Pr(>|t|)"],
             cate_lo20 = lo, cate_hi20 = hi, gap_pct = 100 * (exp(hi) - exp(lo)))
}

d3 <- FL[copy(A$S3), on = "SEQN"]
d3[is.na(castration), castration := 0L]; d3[is.na(suppress), suppress := 0L]
d3[, log_tt := log(tt_ngdl)]

SEEDS <- c(20261117, 1, 2, 3, 4)
cat("=== 全样本 ===\n")
RA <- rbindlist(lapply(SEEDS, function(s) one(d3, s)))
print(RA[, .(seed, N, p_hetero = round(p_hetero, 3), p_blp = round(p_blp, 3),
             blp_mv = round(blp_mv, 3), gap_pct = round(gap_pct, 1))])

cat("\n=== 清洁样本（排除药物去势/抗雄/雌激素/糖皮质激素/阿片/酮康唑）===\n")
RC <- rbindlist(lapply(SEEDS, function(s) one(d3[castration == 0 & suppress == 0], s)))
print(RC[, .(seed, N, p_hetero = round(p_hetero, 3), p_blp = round(p_blp, 3),
             blp_mv = round(blp_mv, 3), gap_pct = round(gap_pct, 1))])

cat("\n=== 汇总 ===\n")
cat(sprintf("全样本   异质性 p: %.3f – %.3f (中位 %.3f)；显著(<0.05)种子数 %d/%d\n",
            min(RA$p_hetero), max(RA$p_hetero), median(RA$p_hetero), sum(RA$p_hetero < .05), nrow(RA)))
cat(sprintf("清洁样本 异质性 p: %.3f – %.3f (中位 %.3f)；显著(<0.05)种子数 %d/%d\n",
            min(RC$p_hetero), max(RC$p_hetero), median(RC$p_hetero), sum(RC$p_hetero < .05), nrow(RC)))
cat(sprintf("全样本   BLP p: %.3f – %.3f；清洁 BLP p: %.3f – %.3f\n",
            min(RA$p_blp), max(RA$p_blp), min(RC$p_blp), max(RC$p_blp)))
cat(sprintf("低-高暴露 CATE 差: 全 %.1f – %.1f pp；清洁 %.1f – %.1f pp\n",
            min(RA$gap_pct), max(RA$gap_pct), min(RC$gap_pct), max(RC$gap_pct)))

fwrite(rbind(RA[, sample := "all"], RC[, sample := "clean"]),
       file.path(BASE, "results/forest_seed_stability.csv"))
cat("\n✅ 完成\n")
