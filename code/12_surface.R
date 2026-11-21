## ============================================================================
## 12_surface.R — 减重幅度 × POP 负荷 的二维剂量-反应曲面 (g-computation)
## ============================================================================
suppressMessages({library(data.table); library(mgcv)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
CONF <- A$CONF

crs <- function(m, cl) {
  X <- model.matrix(m); e <- residuals(m); n <- nrow(X); k <- ncol(X)
  bread <- solve(crossprod(X)); meat <- matrix(0, k, k)
  for (g in unique(cl)) { i <- which(cl == g)
    s <- crossprod(X[i, , drop = FALSE], e[i]); meat <- meat + s %*% t(s) }
  G <- length(unique(cl)); V <- bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
  list(coef = coef(m), V = V, se = sqrt(diag(V)), n = n, G = G)
}

dose_cat <- function(x) {
  factor(fifelse(x <= -2, "1_gain>2%",
         fifelse(x < 2, "2_stable",
         fifelse(x < 5, "3_loss2-5%",
         fifelse(x < 10, "4_loss5-10%", "5_loss>=10%")))),
         levels = c("1_gain>2%","2_stable","3_loss2-5%","4_loss5-10%","5_loss>=10%"))
}

cell_means <- function(dt, yvar, modvar, label) {
  d <- copy(dt)[!is.na(get(yvar)) & !is.na(get(modvar)) & !is.na(dw1)]
  d[, dc := dose_cat(dw1)]
  d[, mt := cut(get(modvar), quantile(get(modvar), c(0, 1/3, 2/3, 1)),
                include.lowest = TRUE, labels = c("T1_low","T2_mid","T3_high"))]
  cf <- CONF[sapply(CONF, function(v) length(unique(d[[v]])) > 1)]
  f <- as.formula(paste0(yvar, " ~ dc * mt + ", paste(cf, collapse = " + "),
                         " + race_f + edu_f + smk"))
  m <- lm(f, data = d); r <- crs(m, factor(paste(d$strata, d$psu)))
  ## g-computation: 每个 (dc, mt) 组合对全样本平均
  grid <- CJ(dc = levels(d$dc), mt = levels(d$mt))
  X0 <- model.matrix(m)
  out <- rbindlist(lapply(seq_len(nrow(grid)), function(k) {
    dd <- copy(d); dd[, dc := grid$dc[k]]; dd[, mt := grid$mt[k]]
    Xk <- model.matrix(delete.response(terms(m)), data = dd)[, colnames(X0), drop = FALSE]
    w <- colMeans(Xk)                       # 平均设计向量
    est <- sum(w * r$coef); se <- sqrt(as.numeric(t(w) %*% r$V %*% w))
    data.table(dc = grid$dc[k], mt = grid$mt[k],
               n = nrow(d[dc == grid$dc[k] & mt == grid$mt[k]]),
               est = est, se = se, lo = est - 1.96*se, hi = est + 1.96*se)
  }))
  cat("\n#####", label, " (N =", r$n, ") #####\n")
  print(dcast(out, dc ~ mt, value.var = "est")[, lapply(.SD, function(z) round(z, 2)), by = dc])
  ## 相对稳定体重的差值
  ref <- out[dc == "2_stable", .(mt, ref = est)]
  out <- merge(out, ref, by = "mt")
  out[, `:=`(diff = est - ref)]
  out[, model := label]
  out
}

C1 <- cell_means(A$S1, "phenoage_accel", "pfas_z", "S1 PhenoAge × PFAS三分位")
C2 <- cell_means(A$S2, "phenoage_accel", "pcb_z",  "S2 PhenoAge × PCB三分位")
A$S3[, log_tt := log(tt_ngdl)]
C3 <- cell_means(A$S3, "log_tt", "pfas_z", "S3 log睾酮 × PFAS三分位")
CELL <- rbindlist(list(C1, C2, C3), fill = TRUE)
fwrite(CELL, file.path(BASE, "results/cell_means.csv"))

## ---------- 连续二维曲面 (GAM 张量积) ----------
surf <- function(dt, yvar, modvar, label) {
  d <- copy(dt)[!is.na(get(yvar)) & !is.na(get(modvar)) & !is.na(dw1)]
  d[, mv := get(modvar)]; d[, dw := pmin(pmax(dw1, -15), 20)]
  cf <- CONF[sapply(CONF, function(v) length(unique(d[[v]])) > 1)]
  m <- gam(as.formula(paste0(yvar, " ~ te(dw, mv, k = c(5,4)) + ",
                             paste(cf, collapse = " + "), " + race_f + edu_f + smk")),
           data = d, method = "REML")
  gdw <- seq(-10, 20, length.out = 40); gmv <- seq(-1.8, 1.8, length.out = 37)
  gr <- CJ(dw = gdw, mv = gmv)
  ## 构造预测矩阵: 其余协变量取均值
  base_d <- d[rep(1, nrow(gr))]
  for (v in cf) base_d[[v]] <- mean(d[[v]])
  base_d[, race_f := factor("NHWhite", levels = levels(d$race_f))]
  base_d[, edu_f  := factor("HS", levels = levels(d$edu_f))]
  base_d[, smk    := 0]
  base_d[, `:=`(dw = gr$dw, mv = gr$mv)]
  pr <- predict(m, newdata = base_d, se.fit = TRUE)
  gr[, `:=`(fit = as.numeric(pr$fit), se = as.numeric(pr$se.fit))]
  ## 以 dw = 0 为参照的差值
  ref <- gr[dw == 0][order(mv)]
  gr <- merge(gr, ref[, .(mv, ref_fit = fit)], by = "mv")
  gr[, delta := fit - ref_fit]
  gr[, model := label]
  cat("\n", label, ": 曲面拟合完成, te() p =",
      format.pval(m$s.table["te(dw,mv)", "p-value"], digits = 3), "\n")
  gr
}
S1s <- surf(A$S1, "phenoage_accel", "pfas_z", "S1")
S2s <- surf(A$S2, "phenoage_accel", "pcb_z",  "S2")
S3s <- surf(A$S3, "log_tt", "pfas_z", "S3")
SURF <- rbindlist(list(S1s, S2s, S3s))
fwrite(SURF, file.path(BASE, "results/surface_grid.csv"))
cat("\ndone\n")
