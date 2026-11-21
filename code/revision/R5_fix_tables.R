## ============================================================================
## R5_fix_tables.R — 两处自查发现的表内错误
##  ① results/liver_decomp.csv 的 dbmi 行与 bmi 行系数完全一致
##     原因不是复制 bug，而是恒等式：dbmi = bmi - bmi_p 且 bmi_p 在协变量集中，
##     故 w5 系数必然等于 bmi 行。真正的分解必须在"不含 bmi_p"的协变量集下估计。
##     → 保留原行（标注 identical_by_construction），另加一行 dbmi_nobmip。
##  ② p_hetero_2s 由 2×单侧 P 得来，未截断，出现 1.149 / 1.341 的"概率"
##     → 截断至 1。
## 输出: results/liver_decomp.csv（增列）、results/forest_seed_stability.csv（截断）
## ============================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", ".")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds")); CONF <- A$CONF
d <- readRDS(file.path(BASE, "data/derived/liver.rds"))
d <- d[ok_liver == TRUE]
d[, `:=`(mv = as.numeric(pfas_z))]
d[, mv2 := mv^2]
d[, `:=`(hsi_liver = 8*(alt/ast),
         fli_liver = 0.953*log(tg2) + 0.718*log(ggt) + 0.053*waist,
         dbmi = bmi - bmi_p)]

crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X); bread <- solve(crossprod(X))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
wald <- function(b, V) as.numeric(t(b) %*% solve(V) %*% b)
fit1 <- function(k, covs, tag) {
  dd <- d[!is.na(get(k))]
  X <- model.matrix(as.formula(paste(k, "~ w5*(mv+mv2) + mv + mv2 +", paste(covs, collapse=" + "),
                                     "+ race_f + edu_f + smk")), data = dd)
  cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  iw <- grep("^w5$", colnames(X)); i1 <- grep("^w5:mv$", colnames(X)); i2 <- grep("^w5:mv2$", colnames(X))
  q <- quantile(dd$mv, c(.1,.5,.9)); L <- cbind(1,q,q^2); th <- as.numeric(L %*% b[c(iw,i1,i2)])
  sdy <- sd(dd[[k]])
  data.table(term = tag, label = ifelse(tag == "dbmi_nobmip",
                                        "BMI change (covariate set WITHOUT prior BMI)",
                                        "BMI change (covariate set WITH prior BMI)"),
             N = nrow(dd), effect0 = b[iw], sd_y = sdy, effect0_sd = b[iw]/sdy,
             p10 = th[1], p50 = th[2], p90 = th[3],
             p_joint = pchisq(wald(b[c(i1,i2)], V[c(i1,i2),c(i1,i2)]), 2, lower.tail = FALSE))
}
CONF_NB <- setdiff(CONF, "bmi_p")
newrows <- rbind(fit1("dbmi", CONF, "dbmi"), fit1("dbmi", CONF_NB, "dbmi_nobmip"))
cat("=== dbmi 两种协变量集 ===\n"); print(newrows[, .(term, N, effect0 = round(effect0,4),
                                                          sd_y = round(sd_y,3), p_joint = signif(p_joint,3))])
cat(sprintf("  两者差 = %.3g（同一模型下应为 0，因为 bmi_p 被部分化）\n",
            newrows[1, effect0] - newrows[2, effect0]))

old <- fread(file.path(BASE, "results/liver_decomp.csv"))
old[, note := fifelse(term == "bmi", "current BMI",
              fifelse(term == "dbmi", "identical_by_construction: dbmi = bmi - bmi_p and bmi_p is a covariate, so the w5 coefficient must equal the bmi row", ""))]
add <- newrows[, .(term, label, N, effect0, sd_y, effect0_sd, p10, p50, p90, p_joint,
                   note = fifelse(term == "dbmi_nobmip",
                                  "re-estimated with bmi_p removed from the covariate set", ""))]
out <- rbind(old[term != "dbmi_nobmip"], add[term == "dbmi_nobmip"])
fwrite(out, file.path(BASE, "results/liver_decomp.csv"))
cat("\n✅ liver_decomp.csv 已增列 dbmi_nobmip 与 note\n")

## ---------- ② p_hetero_2s 截断 -------------------------------------------------
for (f in c(file.path(BASE, "results/forest_seed_stability.csv"),
            file.path(BASE, "evidence_pack/03_modification/E3-03_forest_seed_stability.csv"))) {
  if (!file.exists(f)) { cat("skip (missing):", f, "\n"); next }
  x <- fread(f)
  if ("p_hetero_2s" %in% names(x)) {
    n_bad <- sum(x$p_hetero_2s > 1, na.rm = TRUE)
    x[, p_hetero_2s := pmin(1, p_hetero_2s)]
    fwrite(x, f)
    cat(sprintf("  %s：截断 %d 个 >1 的值\n", basename(f), n_bad))
  }
}
cat("✅ R5 完成\n")
