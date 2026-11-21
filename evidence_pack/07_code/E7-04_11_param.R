## ============================================================================
## 11_param.R — 参数化剂量-反应交互模型（二次型，聚类稳健 SE）
##   CATE(mv) = b_w + b1*mv + b2*mv^2 ；b2>0 表示 U 形，b1>0 表示随暴露单调恶化
## ============================================================================
suppressMessages({library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))

crs <- function(m, cl) {
  X <- model.matrix(m); e <- residuals(m); n <- nrow(X); k <- ncol(X)
  bread <- solve(crossprod(X)); meat <- matrix(0, k, k)
  for (g in unique(cl)) {
    i <- which(cl == g); s <- crossprod(X[i, , drop = FALSE], e[i]); meat <- meat + s %*% t(s)
  }
  G <- length(unique(cl))
  V <- bread %*% meat %*% bread * (G / (G - 1)) * ((n - 1) / (n - k))
  list(coef = coef(m), V = V, se = sqrt(diag(V)), n = n, G = G, X = X)
}
fmt <- function(est, se) sprintf("%.3f (%.3f), p=%s", est, se,
                                 format.pval(2 * pnorm(-abs(est / se)), digits = 3))
evalue <- function(beta, sd_y) { d <- abs(beta) / sd_y; rr <- exp(0.91 * d)
  if (rr < 1) 1 else rr + sqrt(rr * (rr - 1)) }

run_quad <- function(dt, yvar, modvar, label) {
  d <- copy(dt)[!is.na(get(yvar)) & !is.na(get(modvar)) & !is.na(w5)]
  d[, mv := as.numeric(get(modvar))]; d[, mv2 := mv^2]
  ## 剔除在本样本内恒定的协变量（如单周期数据中的周期年份）
  cf <- A$CONF[sapply(A$CONF, function(v) length(unique(d[[v]])) > 1)]
  f <- as.formula(paste0(yvar, " ~ w5 * (mv + mv2) + mv + mv2 + ",
                         paste(cf, collapse = " + "), " + race_f + edu_f + smk"))
  m <- lm(f, data = d); r <- crs(m, factor(paste(d$strata, d$psu)))
  nm <- names(r$coef); V <- r$V
  iw <- grep("^w5$", nm); i1 <- grep("^w5:mv$", nm); i2 <- grep("^w5:mv2$", nm)
  cat("\n######", label, " N =", r$n, " 聚类 =", r$G, " ######\n")
  q <- quantile(d$mv, c(.1, .5, .9))
  th <- r$coef[iw] + r$coef[i1] * q + r$coef[i2] * q^2
  Lm <- cbind(1, q, q^2)                             # 3 x 3
  B <- r$coef[c(iw, i1, i2)]; Vs <- V[c(iw, i1, i2), c(iw, i1, i2)]
  se_th <- sqrt(diag(Lm %*% Vs %*% t(Lm)))
  ## 联合异质性检验 (2 df)
  b_int <- r$coef[c(i1, i2)]; V_int <- V[c(i1, i2), c(i1, i2)]
  W <- as.numeric(t(b_int) %*% solve(V_int) %*% b_int)
  p_joint <- pchisq(W, 2, lower.tail = FALSE)
  cat("  减重主效应 (暴露=0)      :", fmt(r$coef[iw], r$se[iw]), "\n")
  cat("  交互×暴露 (b1)           :", fmt(r$coef[i1], r$se[i1]), "\n")
  cat("  交互×暴露² (b2)          :", fmt(r$coef[i2], r$se[i2]), "\n")
  cat("  联合异质性检验 (2 df)    : Wald =", round(W, 2), ", p =", format.pval(p_joint, digits = 3), "\n")
  cat("  暴露主效应               :", fmt(r$coef[grep("^mv$", nm)], r$se[grep("^mv$", nm)]), "\n")
  cat("  减重效应 @ P10/P50/P90:\n")
  print(data.table(exposure_z = round(q, 2), effect = round(th, 3), se = round(se_th, 3),
                   p = signif(2 * pnorm(-abs(th / se_th)), 3)))
  ## 零交叉
  disc <- r$coef[i1]^2 - 4 * r$coef[i2] * r$coef[iw]
  xzero <- NA_real_
  if (r$coef[i2] != 0 && disc >= 0) {
    rr <- sort((-r$coef[i1] + c(-1, 1) * sqrt(disc)) / (2 * r$coef[i2]))
    inr <- rr[rr >= quantile(d$mv, .02) & rr <= quantile(d$mv, .98)]
    if (length(inr)) xzero <- inr[1]
    cat("  曲线过零点 z =", paste(round(rr, 2), collapse = " / "), " (区间内:", round(xzero, 2), ")\n")
  }
  ## 转折点
  turn <- if (r$coef[i2] != 0) -r$coef[i1] / (2 * r$coef[i2]) else NA
  cat("  转折点 z =", round(turn, 2), "( b2>0 为极小值 = U 形 )\n")
  sd_y <- sd(d[[yvar]])
  cat("  P90-P10 效应差 E-value =", round(evalue((th[3] - th[1]) / 2, sd_y), 2), "\n")
  data.table(model = label, N = r$n, G = r$G,
             eff0 = r$coef[iw], se0 = r$se[iw], b1 = r$coef[i1], se_b1 = r$se[i1],
             b2 = r$coef[i2], se_b2 = r$se[i2], wald = W, p_joint = p_joint,
             e10 = th[1], se10 = se_th[1], e50 = th[2], se50 = se_th[2],
             e90 = th[3], se90 = se_th[3], turn = turn, xzero = xzero,
             ev = evalue((th[3] - th[1]) / 2, sd_y))
}

R <- list()
R[[1]] <- run_quad(A$S1,  "phenoage_accel", "pfas_z", "S1  PhenoAge × ΣPFAS")
R[[2]] <- run_quad(A$S2,  "phenoage_accel", "pcb_z",  "S2  PhenoAge × ΣPCB")
R[[3]] <- run_quad(A$S2b, "phenoage_accel", "ocp_z",  "S2b PhenoAge × ΣOCP")
A$S3[, log_tt := log(tt_ngdl)]; A$S4[, log_ft := log(ft_nmol)]
R[[4]] <- run_quad(A$S3, "log_tt", "pfas_z", "S3  log总睾酮 × ΣPFAS")
R[[5]] <- run_quad(A$S4, "log_ft", "pfas_z", "S4  log游离睾酮 × ΣPFAS")

## 连续处理剂量反应（每 5%）
cat("\n\n=== 连续剂量反应 ===\n")
for (lab in c("S1","S2","S3")) {
  d <- copy(A[[lab]])
  yy <- c(S1="phenoage_accel", S2="phenoage_accel", S3="log_tt")[[lab]]
  if (lab == "S3") d[, log_tt := log(tt_ngdl)]
  mm <- c(S1="pfas_z", S2="pcb_z", S3="pfas_z")[[lab]]
  d <- d[!is.na(get(yy)) & !is.na(get(mm)) & !is.na(dw1)]
  d[, mv := get(mm)]; d[, dw5 := dw1 / 5]
  cf <- A$CONF[sapply(A$CONF, function(v) length(unique(d[[v]])) > 1)]
  m <- lm(as.formula(paste0(yy, " ~ dw5 * mv + mv + ", paste(cf, collapse=" + "),
                            " + race_f + edu_f + smk")), data = d)
  r <- crs(m, factor(paste(d$strata, d$psu))); nm <- names(r$coef)
  cat(sprintf("%s: 每 5%% 减重 %s ; 交互 %s\n", lab,
              fmt(r$coef[grep("^dw5$", nm)], r$se[grep("^dw5$", nm)]),
              fmt(r$coef[grep("^dw5:mv$", nm)], r$se[grep("^dw5:mv$", nm)])))
}

RES <- rbindlist(R)
print(RES[, .(model, N, G, b1 = round(b1,3), b2 = round(b2,3), p_joint = signif(p_joint,3),
              e10 = round(e10,2), e90 = round(e90,2), turn = round(turn,2),
              xzero = round(xzero,2), ev = round(ev,2))])
fwrite(RES, file.path(BASE, "results/param_interaction.csv"))
cat("done\n")
