## ============================================================================
## R6_pfas_vs_loss.R — 检验"减重释放储存污染物"这一叙事在人群层面的前提
##   若前提成立，减重者血清 PFAS 应更高（释放效应 + 体重下降的浓缩效应同向）。
##   实算：是否更高、剂量关系、各同系物、对协变量集的敏感性。
## 输出: results/revision/R6_pfas_by_loss.csv, R6_congeners.csv, R6_dose.csv
## ============================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", ".")
A  <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
I  <- readRDS(file.path(BASE, "data/derived/inflam.rds"))
OUT <- file.path(BASE, "results/revision"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CONF <- A$CONF

crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X); bread <- solve(crossprod(X))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}

d <- copy(I)[!is.na(log_pfas) & !is.na(w5)]
cat("分析域 N =", nrow(d), " 减重组 =", sum(d$w5), "\n")
cat("dw1 符号检查：w5==1 时 dw1 中位 =", median(d[w5==1, dw1], na.rm=TRUE),
    "；w5==0 时 =", median(d[w5==0, dw1], na.rm=TRUE), "\n\n")

fit <- function(dd, y, covs, tag, label) {
  vars <- unique(c("w5", y, covs, "race_f", "edu_f", "smk", "strata", "psu"))
  dd <- dd[complete.cases(dd[, ..vars])]
  f <- as.formula(paste(y, "~ w5 +", paste(covs, collapse = " + "), "+ race_f + edu_f + smk"))
  X <- model.matrix(f, data = dd)
  keep <- apply(X, 2, function(z) sd(z, na.rm = TRUE) > 0) | colnames(X) == "(Intercept)"
  X <- X[, keep, drop = FALSE]
  m <- lm.fit(X, dd[[y]])
  V <- crsV(X, m$residuals, factor(paste(dd$strata, dd$psu)))
  i <- grep("^w5$", colnames(X))
  b <- m$coefficients[i]; se <- sqrt(V[i, i])
  data.table(tag = tag, spec = label, N = nrow(dd), treated = sum(dd$w5),
             pct = 100 * (exp(b) - 1), se_pct = 100 * se,
             lo = 100 * (exp(b - 1.96 * se) - 1), hi = 100 * (exp(b + 1.96 * se) - 1),
             p = 2 * pnorm(-abs(b / se)))
}

## ---------- ⓿ 逐步诊断（确认符号翻转来自哪里）------------------------------
cat("=== ⓿ 逐步加协变量 ===\n")
step1 <- function(f, tag) {
  X <- model.matrix(f, data = d)
  keep <- apply(X, 2, function(z) sd(z) > 0) | colnames(X) == "(Intercept)"
  X <- X[, keep, drop = FALSE]
  m <- lm.fit(X, d$log_pfas); V <- crsV(X, m$residuals, factor(paste(d$strata, d$psu)))
  i <- grep("^w5$", colnames(X)); b <- m$coefficients[i]; se <- sqrt(V[i,i])
  cat(sprintf("  %-46s %+7.2f%% (SE %.2f) p=%.2g\n", tag, 100*(exp(b)-1), 100*se, 2*pnorm(-abs(b/se))))
}
step1(log_pfas ~ w5, "raw")
for (v in c("race_f","edu_f","smk","age","female","pir","bmi_p")) {
  step1(as.formula(paste("log_pfas ~ w5 +", v)), paste("+", v))
}
step1(as.formula(paste("log_pfas ~ w5 +", paste(CONF, collapse=" + "), "+ race_f + edu_f + smk")), "full CONF")
cat("\n")

## ---------- ① 总 PFAS：四种协变量集 ----------------------------------------
cat("=== ① log(ΣPFAS molar) ~ 减重≥5% ===\n")
R <- list()
R[[1]] <- fit(d, "log_pfas", character(0), "crude", "crude")
R[[2]] <- fit(d, "log_pfas", CONF,  "full",  "adjusted: age, sex, PIR, prior BMI, cotinine, alcohol, activity, comorbidity, lipids, eGFR, year")
CONF_NB <- setdiff(CONF, "bmi_p")
R[[3]] <- fit(d, "log_pfas", CONF_NB, "noprior", "adjusted WITHOUT prior BMI")
R[[4]] <- fit(d, "log_pfas", c(setdiff(CONF, "bmi_p"), "bmi"), "nowbmi", "adjusted WITH current BMI (collider on the w5 path; shown for transparency)")
RES <- rbindlist(R)
for (i in seq_len(nrow(RES))) cat(sprintf("  %-8s %-72s %+6.2f%% (%+6.2f, %+6.2f) p=%.3g\n",
  RES$tag[i], substr(RES$spec[i],1,72), RES$pct[i], RES$lo[i], RES$hi[i], RES$p[i]))
fwrite(RES, file.path(OUT, "R6_pfas_by_loss.csv"))

## ---------- ② 各同系物 ------------------------------------------------------
cat("\n=== ② 各同系物 ~ 减重≥5%（全协变量）===\n")
cg <- c("PFOA","PFOS","PFHxS","PFNA","PFDA","PFUnDA")
RC <- rbindlist(lapply(cg, function(v) {
  y <- paste0("log_", tolower(v))
  if (!y %in% names(d)) { d[, (y) := log(get(v))] }
  r <- fit(d[!is.na(get(y))], y, CONF, v, v); r[, congener := v]; r[, spec := NULL]; r
}))
for (i in seq_len(nrow(RC))) cat(sprintf("  %-7s %+6.2f%% (%+6.2f, %+6.2f)  p=%.3g\n",
  RC$congener[i], RC$pct[i], RC$lo[i], RC$hi[i], RC$p[i]))
fwrite(RC, file.path(OUT, "R6_congeners.csv"))

## ---------- ③ 剂量关系：按减重幅度分档 --------------------------------------
cat("\n=== ③ 剂量关系（参照组：体重稳定 −5%~+5%）===\n")
d[, band := fcase(dw1 >= 10, ">=10% loss",
                  dw1 >= 5 & dw1 < 10, "5-10% loss",
                  dw1 > -5 & dw1 < 5, "stable (+-5%)",
                  dw1 <= -5, ">=5% gain",
                  default = NA_character_)]
d[, band := factor(band, levels = c("stable (+-5%)", "5-10% loss", ">=10% loss", ">=5% gain"))]
cat("分档计数:\n"); print(d[!is.na(band), .N, by = band])
db <- d[!is.na(band)]
db <- db[complete.cases(db[, c("band","log_pfas",CONF,"race_f","edu_f","smk","strata","psu"), with=FALSE])]
X <- model.matrix(as.formula(paste("~ band +", paste(CONF, collapse=" + "), "+ race_f + edu_f + smk")), data = db)
keep <- apply(X, 2, function(z) sd(z, na.rm = TRUE) > 0) | colnames(X) == "(Intercept)"
X <- X[, keep, drop = FALSE]
m <- lm.fit(X, db$log_pfas)
V <- crsV(X, m$residuals, factor(paste(db$strata, db$psu)))
RD <- data.table()
for (lv in levels(d$band)[-1]) {
  j <- grep(paste0("^band", gsub("([%()+-])", "\\\\\\1", lv), "$"), colnames(X))
  if (!length(j)) { cat("  未匹配:", lv, "\n"); next }
  b <- m$coefficients[j]; se <- sqrt(V[j, j])
  RD <- rbind(RD, data.table(band = lv, N = sum(db$band == lv),
                             pct = 100*(exp(b)-1), lo = 100*(exp(b-1.96*se)-1), hi = 100*(exp(b+1.96*se)-1),
                             p = 2*pnorm(-abs(b/se))))
}
print(RD)
fwrite(RD, file.path(OUT, "R6_dose.csv"))

## ---------- ④ PFAS 与体重的横断面关系（解释稀释效应）------------------------
cat("\n=== ④ log(ΣPFAS) 与当前 BMI、一年前 BMI 的偏相关（同协变量集）===\n")
for (v in c("bmi", "bmi_p")) {
  cv <- setdiff(CONF, v)
  dc <- d[complete.cases(d[, c("log_pfas", cv, v, "race_f","edu_f","smk","strata","psu"), with=FALSE])]
  X2 <- model.matrix(as.formula(paste("~", paste(c(cv, v), collapse=" + "), "+ race_f + edu_f + smk")), data = dc)
  keep2 <- apply(X2, 2, function(z) sd(z, na.rm=TRUE) > 0) | colnames(X2) == "(Intercept)"
  X2 <- X2[, keep2, drop = FALSE]
  m2 <- lm.fit(X2, dc$log_pfas); V2 <- crsV(X2, m2$residuals, factor(paste(dc$strata, dc$psu)))
  j <- grep(paste0("^", v, "$"), colnames(X2)); b <- m2$coefficients[j]; se <- sqrt(V2[j,j])
  cat(sprintf("  %-6s: 每 +1 kg/m2 → %+0.3f%% (SE %.3f, p=%.3g)\n", v, 100*(exp(b)-1), 100*se, 2*pnorm(-abs(b/se))))
}
cat("✅ R6 完成\n")
