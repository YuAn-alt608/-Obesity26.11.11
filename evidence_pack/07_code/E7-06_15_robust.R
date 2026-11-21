## ============================================================================
## 15_robust.R — 阴性对照结局、安慰剂检验、阈值敏感性、分周期复现
## ============================================================================
suppressMessages({library(data.table); library(grf)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
CONF <- A$CONF
set.seed(20261117)

crs <- function(m, cl) {
  X <- model.matrix(m); e <- residuals(m); n <- nrow(X); k <- ncol(X)
  bread <- solve(crossprod(X)); meat <- matrix(0, k, k)
  for (g in unique(cl)) { i <- which(cl == g)
    s <- crossprod(X[i, , drop = FALSE], e[i]); meat <- meat + s %*% t(s) }
  G <- length(unique(cl)); V <- bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
  list(coef = coef(m), V = V, se = sqrt(diag(V)), n = n, G = G)
}
quad <- function(d, y, mv, treat = "w5", ctrl_only_stable = FALSE) {
  d <- copy(d)[!is.na(get(y)) & !is.na(get(mv)) & !is.na(get(treat))]
  if (ctrl_only_stable) d <- d[get(treat) == 1 | abs(dw1) < 2]
  d[, mz := get(mv)]; d[, mz2 := mz^2]
  cf <- CONF[sapply(CONF, function(v) length(unique(d[[v]])) > 1)]
  f <- as.formula(paste0(y, " ~ ", treat, " * (mz + mz2) + mz + mz2 + ",
                         paste(cf, collapse = " + "), " + race_f + edu_f + smk"))
  m <- lm(f, data = d); r <- crs(m, factor(paste(d$strata, d$psu)))
  nm <- names(r$coef)
  i1 <- grep(paste0("^", treat, ":mz$"), nm); i2 <- grep(paste0("^", treat, ":mz2$"), nm)
  b <- r$coef[c(i1, i2)]; V <- r$V[c(i1, i2), c(i1, i2)]
  W <- as.numeric(t(b) %*% solve(V) %*% b)
  list(n = r$n, b1 = r$coef[i1], p1 = 2*pnorm(-abs(r$coef[i1]/r$se[i1])),
       b2 = r$coef[i2], p2 = 2*pnorm(-abs(r$coef[i2]/r$se[i2])),
       W = W, pjoint = pchisq(W, 2, lower.tail = FALSE))
}

cat("=========== 1) 阴性对照结局: 身高 ===========\n")
h <- A$S1[!is.na(ht_meas)]
r_h <- quad(h, "ht_meas", "pfas_z")
cat(sprintf("身高 ~ 处理×PFAS: b1=%.3f (p=%.3f), b2=%.3f (p=%.3f), 联合 p=%.3f, N=%d\n",
            r_h$b1, r_h$p1, r_h$b2, r_h$p2, r_h$pjoint, r_h$n))
## 森林版
X <- cbind(as.matrix(h[, ..CONF]), model.matrix(~ race_f + edu_f + smk - 1, data = h),
           pfas_z = h$pfas_z)
storage.mode(X) <- "double"
cfh <- causal_forest(X, h$ht_meas, as.numeric(h$w5), num.trees = 2000,
                     tune.parameters = "all", seed = 1)
cat("森林阴性对照 差异预测 t =", round(test_calibration(cfh)[2, 3], 3),
    " (|t|<2 表示无虚假异质性)\n")

cat("\n=========== 2) 安慰剂检验 (随机置换处理) ===========\n")
pl <- sapply(1:100, function(b) {
  d <- copy(A$S1); d[, w5 := sample(w5)]
  rr <- quad(d, "phenoage_accel", "pfas_z")
  c(rr$pjoint, rr$b2)
})
cat(sprintf("置换后联合检验 p<0.05 的比例 = %.2f (期望 0.05)\n", mean(pl[1, ] < 0.05)))
cat(sprintf("真实数据联合 p = %.4f ；置换分布 p 的中位数 = %.3f\n",
            quad(A$S1, "phenoage_accel", "pfas_z")$pjoint, median(pl[1, ])))

cat("\n=========== 3) 阈值 / 对照定义敏感性 ===========\n")
A$S1[, w10 := as.integer(dw1 >= 10)]
tab <- rbind(
  data.table(spec = "主分析: >=5% vs 其他", quad(A$S1, "phenoage_accel", "pfas_z", "w5")[c("n","b1","p1","b2","p2","pjoint")]),
  data.table(spec = ">=10% vs 其他",        quad(A$S1, "phenoage_accel", "pfas_z", "w10")[c("n","b1","p1","b2","p2","pjoint")]),
  data.table(spec = ">=5% vs 稳定体重(|Δ|<2%)", quad(A$S1, "phenoage_accel", "pfas_z", "w5", TRUE)[c("n","b1","p1","b2","p2","pjoint")])
)
print(tab, digits = 3)

cat("\n=========== 4) 分周期复现 (交互项) ===========\n")
byc <- rbindlist(lapply(sort(unique(A$S1$cycle)), function(cc) {
  d <- A$S1[cycle == cc]
  r <- quad(d, "phenoage_accel", "pfas_z")
  data.table(cycle = cc, n = r$n, b1 = r$b1, p1 = r$p1, b2 = r$b2, p2 = r$p2, pjoint = r$pjoint)
}))
print(byc, digits = 3)

cat("\n=========== 5) 加权 AIPW (调查权重) ===========\n")
d <- copy(A$S1)
d[, wgt := wt_mec / 6]     # 6 个周期的 MEC 权重
X <- cbind(as.matrix(d[, ..CONF]), model.matrix(~ race_f + edu_f + smk - 1, data = d))
storage.mode(X) <- "double"
e0 <- as.numeric(regression_forest(X, d$w5, num.trees = 2000, tune.parameters = "all")$predictions)
rf1 <- regression_forest(X[d$w5 == 1, ], d$phenoage_accel[d$w5 == 1], num.trees = 2000, tune.parameters = "all")
rf0 <- regression_forest(X[d$w5 == 0, ], d$phenoage_accel[d$w5 == 0], num.trees = 2000, tune.parameters = "all")
mu1 <- as.numeric(predict(rf1, X)$predictions); mu0 <- as.numeric(predict(rf0, X)$predictions)
d[, aipw := mu1 - mu0 + w5 * (phenoage_accel - mu1) / e0 - (1 - w5) * (phenoage_accel - mu0) / (1 - e0)]
wsum <- function(x) sum(x * d$wgt) / sum(d$wgt)
cat("加权 ATE =", round(wsum(d$aipw), 3), " 未加权 ATE =", round(mean(d$aipw), 3), "\n")
q <- quantile(d$pfas_z, c(.1, .5, .9))
for (k in seq_along(q)) {
  s <- abs(d$pfas_z - q[k]) < 0.35
  cat(sprintf("  PFAS z=%.2f 附近加权 CATE = %.3f (n=%d)\n", q[k], wsum(d$aipw[s]), sum(s)))
}

RES <- list(neg_control = r_h, placebo = pl, thresh = tab, bycycle = byc)
saveRDS(RES, file.path(BASE, "results/robust.rds"))
fwrite(tab, file.path(BASE, "results/robust_threshold.csv"))
fwrite(byc, file.path(BASE, "results/robust_bycycle.csv"))
cat("\ndone\n")
