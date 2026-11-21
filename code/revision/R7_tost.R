## ============================================================================
## R7_tost.R — 等价性检验：把"未检出"升级为"可以排除"
##   两个主结局：① PhenoAge 加速（年）② log(ΣPFAS)
##   对每个 margin 报 TOST 双单侧 p 与 90% CI（等价于双单侧 5% 的区间）
## 输出: results/revision/R7_tost.csv
## ============================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", ".")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
I <- readRDS(file.path(BASE, "data/derived/inflam.rds"))
CONF <- A$CONF
OUT <- file.path(BASE, "results/revision"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X); bread <- solve(crossprod(X))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}

tost <- function(est, se, margins, label, unit) {
  rbindlist(lapply(margins, function(m) {
    z1 <- (est - m) / se; z2 <- (est + m) / se
    p_up <- pnorm(z1)         # H0: theta >= +m ; 拒绝需 theta_hat 显著低于 +m
    p_lo <- 1 - pnorm(z2)     # H0: theta <= -m ; 拒绝需 theta_hat 显著高于 -m
    pt <- max(p_up, p_lo)
    data.table(outcome = label, unit = unit, margin = m, est = est, se = se,
               ci90_lo = est - 1.645 * se, ci90_hi = est + 1.645 * se,
               p_lower = p_lo, p_upper = p_up, p_tost = pt,
               equivalent = pt < 0.05)
  }))
}

## ① PhenoAge 调整后主效应（与 R4 同规格）
S1 <- copy(A$S1)[!is.na(phenoage_accel) & !is.na(w5)]
S1[, cyc := factor(cycle)]
CONF_S1 <- CONF[sapply(CONF, function(v) length(unique(S1[[v]])) > 1)]
X <- model.matrix(as.formula(paste("phenoage_accel ~ w5 +", paste(CONF_S1, collapse=" + "), "+ race_f + edu_f + smk")), data = S1)
V <- crsV(X, residuals(lm.fit(X, S1$phenoage_accel)), factor(paste(S1$strata, S1$psu)))
i <- grep("^w5$", colnames(X))
est1 <- lm.fit(X, S1$phenoage_accel)$coefficients[i]; se1 <- sqrt(V[i, i])
cat(sprintf("PhenoAge: %.3f 年 (SE %.3f), N=%d\n", est1, se1, nrow(S1)))

## ② PFAS（与 R6 全协变量同规格）
d <- I[!is.na(log_pfas) & !is.na(w5)]
d <- d[complete.cases(d[, c("log_pfas","w5",CONF,"race_f","edu_f","smk","strata","psu"), with=FALSE])]
X2 <- model.matrix(as.formula(paste("log_pfas ~ w5 +", paste(CONF, collapse=" + "), "+ race_f + edu_f + smk")), data = d)
V2 <- crsV(X2, residuals(lm.fit(X2, d$log_pfas)), factor(paste(d$strata, d$psu)))
j <- grep("^w5$", colnames(X2))
b2 <- lm.fit(X2, d$log_pfas)$coefficients[j]
est2 <- 100 * (exp(b2) - 1); se2_log <- sqrt(V2[j, j]); se2 <- 100 * se2_log * exp(b2)  # delta
cat(sprintf("PFAS: %+.2f%% (SE %.2f), N=%d\n", est2, se2, nrow(d)))

R <- rbind(
  tost(est1, se1, c(0.5, 1.0), "PhenoAge acceleration", "years"),
  tost(est2, se2, c(5, 10), "Serum sum PFAS", "percent")
)
R[, `:=`(ci95_lo = est - 1.96*se, ci95_hi = est + 1.96*se)]
print(R[, .(outcome, margin, est = round(est,3), ci95 = sprintf("%.3f to %.3f", ci95_lo, ci95_hi),
            ci90 = sprintf("%.3f to %.3f", ci90_lo, ci90_hi),
            p_tost = signif(p_tost, 3), equivalent)])
fwrite(R, file.path(OUT, "R7_tost.csv"))
cat("\n判读：p_tost<0.05 且 90%CI 完全落在 ±margin 内，方可在该 margin 上主张等价。\n")
cat("✅ R7 完成\n")
