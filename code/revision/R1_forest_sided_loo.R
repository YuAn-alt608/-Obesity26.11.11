## ============================================================================
## R1_forest_sided_loo.R — 回应审稿意见 #3（森林 P 单侧 / 不可复现 / 无留一法）
##   ① 复现归档规格（4000 树, e0 2000 树, seed 20261117）并给出双侧 P
##   ② 树数 × 倾向得分树数的 2×2 归因：P 到底是种子不稳还是蒙特卡洛不稳
##   ③ 树数网格 500–8000，看 P 的分布
##   ④ 留一周期的留一法（G / H / I 各剔除一次）
##   ⑤ 剔除 BLP 回归 Cook 距离最大的 2 人
## 输出: results/revision/R1_*.csv, logs/revision_R1.log
## ============================================================================
suppressMessages({library(data.table); library(grf)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
A  <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
FL <- readRDS(file.path(BASE, "data/derived/rx_flags.rds"))
CONF <- A$CONF
OUT <- file.path(BASE, "results/revision"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
TWOS <- function(p) pmin(1, 2 * p)

mkX <- function(dt) {
  X <- cbind(as.matrix(dt[, ..CONF]), model.matrix(~ race_f + edu_f + smk - 1, data = dt))
  X <- cbind(X, as.matrix(dt[, .(pfas_z)]))
  mc <- grep("^miss_", names(dt), value = TRUE)
  if (length(mc)) X <- cbind(X, as.matrix(dt[, ..mc]))
  storage.mode(X) <- "double"; X
}

## d 必须已含 log_tt / pfas_z / w5
fit_cf <- function(d, seed = 20261117, ntree = 4000, ntree_e = 2000, tag = "") {
  d <- copy(d)[!is.na(log_tt) & !is.na(pfas_z) & !is.na(w5)]
  d[, mv := as.numeric(pfas_z)]
  X <- mkX(d); Y <- d$log_tt; W <- as.numeric(d$w5)
  cl <- as.integer(factor(paste(d$strata, d$psu)))
  set.seed(seed)
  e0 <- as.numeric(regression_forest(X, W, num.trees = ntree_e, tune.parameters = "all",
                                     seed = seed)$predictions)
  keep <- e0 >= 0.05 & e0 <= 0.95
  d <- d[keep]; X <- X[keep, , drop = FALSE]; Y <- Y[keep]; W <- W[keep]
  e0 <- e0[keep]; cl <- cl[keep]
  cf <- causal_forest(X, Y, W, W.hat = e0, num.trees = ntree, honesty = TRUE,
                      tune.parameters = "all", clusters = cl, seed = seed)
  cal <- test_calibration(cf)
  blp <- best_linear_projection(cf, A = cbind(mv = d$mv))
  tau <- predict(cf)$predictions
  lo <- mean(tau[d$mv <= quantile(d$mv, .2)]); hi <- mean(tau[d$mv >= quantile(d$mv, .8)])
  p1 <- cal["differential.forest.prediction", "Pr(>t)"]
  list(cf = cf, d = d, tau = tau, p1 = p1, cal = cal, blp = blp,
       row = data.table(tag = tag, seed = seed, ntree = ntree, ntree_e = ntree_e,
                        N = nrow(d), treated = sum(W), clusters = length(unique(cl)),
                        p_hetero_1s = p1, p_hetero_2s = TWOS(p1),
                        blp_slope = blp["mv", "Estimate"], p_blp = blp["mv", "Pr(>|t|)"],
                        cate_lo20 = lo, cate_hi20 = hi, gap_pct = 100 * (exp(hi) - exp(lo))))
}

## ---------- 数据 ------------------------------------------------------------
d3 <- FL[copy(A$S3), on = "SEQN"]
d3[is.na(castration), castration := 0L]; d3[is.na(suppress), suppress := 0L]
d3[, log_tt := log(tt_ngdl)]
cat("S3 全样本 N =", nrow(d3), " 处理组 =", sum(d3$w5, na.rm = TRUE),
    " 周期 =", paste(sort(unique(d3$cycle)), collapse = "/"), "\n")

rows <- list()

## ---------- ① 归档规格复现 --------------------------------------------------
cat("\n=== ① 归档规格（4000/2000/20261117）复现 ===\n")
a <- fit_cf(d3, 20261117, 4000, 2000, "archival_4000_2000")
print(a$cal); rows[[length(rows)+1]] <- a$row
cat(sprintf("  单侧 P = %.5f ; 双侧 P = %.5f ; BLP p = %.4f\n",
            a$p1, TWOS(a$p1), a$row$p_blp))

## ---------- ② 2×2 归因 ------------------------------------------------------
cat("\n=== ② 树数 2×2（蒙特卡洛 vs 倾向得分树数）===\n")
grid22 <- CJ(ntree = c(3000, 4000), ntree_e = c(1500, 2000))
for (i in seq_len(nrow(grid22))) {
  g <- grid22[i]
  r <- fit_cf(d3, 20261117, g$ntree, g$ntree_e,
              sprintf("mc_%d_%d", g$ntree, g$ntree_e))
  rows[[length(rows)+1]] <- r$row
  cat(sprintf("  ntree=%5d e0=%5d : 单侧 %.4f  双侧 %.4f\n", g$ntree, g$ntree_e, r$p1, TWOS(r$p1)))
}

## ---------- ③ 树数网格 ------------------------------------------------------
cat("\n=== ③ 树数网格（seed 20261117, e0 固定 2000）===\n")
for (nt in c(500, 1000, 2000, 6000, 8000)) {
  r <- fit_cf(d3, 20261117, nt, 2000, sprintf("grid_%d", nt))
  rows[[length(rows)+1]] <- r$row
  cat(sprintf("  ntree=%5d : 单侧 %.4f  双侧 %.4f  BLP p %.3f\n", nt, r$p1, TWOS(r$p1), r$row$p_blp))
}

## ---------- ④ 留一周期 ------------------------------------------------------
cat("\n=== ④ 留一周期（每周期剔除一次，4000 树）===\n")
for (cy in sort(unique(d3$cycle))) {
  r <- fit_cf(d3[cycle != cy], 20261117, 4000, 2000, paste0("leaveout_", cy))
  rows[[length(rows)+1]] <- r$row
  cat(sprintf("  剔除 %s : N=%4d 单侧 %.4f  双侧 %.4f\n", cy, r$row$N, r$p1, TWOS(r$p1)))
}

## ---------- ⑤ 剔除 Cook 距离最大的两人 --------------------------------------
cat("\n=== ⑤ 剔除 BLP 回归 Cook 距离最大的 2 人 ===\n")
dd <- a$d
bl <- lm(tau ~ mv, data = data.table(tau = a$tau, mv = dd$mv))
cd <- cooks.distance(bl)
top2 <- dd$SEQN[order(-cd)][1:2]
cat("  受影响最大 SEQN:", paste(top2, collapse = ", "),
    " tau =", paste(round(a$tau[dd$SEQN %in% top2], 3), collapse = ", "),
    " mv =", paste(round(dd$mv[dd$SEQN %in% top2], 2), collapse = ", "), "\n")
r <- fit_cf(d3[!SEQN %in% top2], 20261117, 4000, 2000, "drop_cook2")
rows[[length(rows)+1]] <- r$row
cat(sprintf("  剔除 2 人 : N=%4d 单侧 %.4f  双侧 %.4f\n", r$row$N, r$p1, TWOS(r$p1)))

## ---------- 汇总 ------------------------------------------------------------
RES <- rbindlist(rows)
print(RES[, .(tag, N, p_hetero_1s = round(p_hetero_1s, 4), p_hetero_2s = round(p_hetero_2s, 4),
              p_blp = round(p_blp, 3), gap_pct = round(gap_pct, 1))])
fwrite(RES, file.path(OUT, "R1_forest_heterogeneity.csv"))

## 树数网格的不稳定性摘要
gr <- RES[grepl("^grid_|^mc_", tag)]
cat(sprintf("\n蒙特卡洛网格：单侧 P %.4f–%.4f（中位 %.4f）；双侧 P %.4f–%.4f\n",
            min(gr$p_hetero_1s), max(gr$p_hetero_1s), median(gr$p_hetero_1s),
            min(gr$p_hetero_2s), max(gr$p_hetero_2s)))
cat("✅ R1 完成\n")
