## ============================================================================
## 32_subgroup_highsens.R — 突破口一：高敏亚群 + 单同系物
##   预设亚组：男性 20–50 岁、ΔW1y ≥ 10%、PFOS/PFOA 单测
##   结局：总睾酮 / SHBG / 游离睾酮 / %游离 / FAI / T:E2
##   同时给出：① 每组实际 N 与可检测效应量（功效） ② 全年龄/全减重幅度对照
## ============================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", "."); RES <- file.path(BASE, "results")
A  <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
FL <- readRDS(file.path(BASE, "data/derived/rx_flags.rds"))

crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X)
  if (n <= k + 2) return(matrix(NA_real_, k, k))
  bread <- tryCatch(solve(crossprod(X)), error = function(e) NULL)
  if (is.null(bread)) return(matrix(NA_real_, k, k))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
wald <- function(b, V) { if (any(is.na(V))) return(NA_real_); as.numeric(t(b) %*% solve(V) %*% b) }

prep <- function(x) {
  x <- FL[copy(x), on = "SEQN"]
  x[is.na(castration), castration := 0L]; x[is.na(suppress), suppress := 0L]
  x[, `:=`(tt_nmol = tt_ngdl/28.84, e2_pmol = e2_pgml*3.671)]
  x[, `:=`(fai = tt_nmol/shbg*100, te2 = tt_nmol/e2_pmol)]
  x[, `:=`(y_tt=log(tt_nmol), y_shbg=log(shbg), y_ft=log(ft_nmol),
           y_pctfree=log(pct_free), y_fai=log(fai), y_te2=log(te2))]
  for (v in c("PFOA","PFOS")) x[, (paste0("z_",v)) := as.numeric(scale(log(get(v)))), by = cycle]
  x
}
dall <- prep(A$S4)
d    <- dall[castration == 0 & suppress == 0]

YK <- c("y_tt","y_shbg","y_ft","y_pctfree","y_fai","y_te2")
LAB <- c(y_tt="Total testosterone", y_shbg="SHBG", y_ft="Free testosterone",
         y_pctfree="%Free testosterone", y_fai="FAI", y_te2="T:E2 ratio")

run <- function(dd, modvar, k, tag) {
  dd <- dd[!is.na(get(k)) & !is.na(get(modvar)) & w5 == 1]   # 仅处理组：直接看反应
  if (nrow(dd) < 25) return(data.table(tag=tag, outcome=k, N=nrow(dd), est=NA_real_, se=NA_real_,
                                       p=NA_real_, lo=NA_real_, hi=NA_real_))
  dd[, mv := as.numeric(get(modvar))]
  ## 小样本用精简协变量（避免过拟合与奇异矩阵）：age, 基线BMI, 吸烟, 周期, 共病
  CV <- if (nrow(dd) < 200) c("age","bmi_p","log_cot","yr","comorbid") else A$CONF
  CV <- CV[sapply(CV, function(v) length(unique(dd[[v]])) > 1)]
  f <- as.formula(paste(k, "~ mv +", paste(CV, collapse=" + "),
                        if (nrow(dd) >= 200) "+ race_f + edu_f + smk" else ""))
  X <- model.matrix(f, data = dd); cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl)
  j <- grep("^mv$", colnames(X)); b <- m$coefficients[j]; se <- sqrt(V[j,j])
  data.table(tag=tag, outcome=k, N=nrow(dd), est=100*(exp(b)-1), se=100*(exp(b)-1)*se,
             p=2*pnorm(-abs(b/se)), lo=100*(exp(b-1.96*se)-1), hi=100*(exp(b+1.96*se)-1))
}

## ---------- 亚组定义 ---------------------------------------------------------
SUB <- list(
  "全部（对照）"                  = d,
  "20-50岁"                       = d[age >= 20 & age <= 50],
  "20-50岁 且 减重>=10%"          = d[age >= 20 & age <= 50 & dw1 >= 10],
  "20-50岁 且 减重>=10%（含用药）" = dall[age >= 20 & age <= 50 & dw1 >= 10]
)

cat("========== 各亚组样本量 ==========\n")
for (nm in names(SUB)) {
  s <- SUB[[nm]]
  cat(sprintf("%-30s 总 N=%4d  处理组 N=%3d\n", nm, nrow(s), sum(s$w5 == 1, na.rm=TRUE)))
}

cat("\n========== 处理组内：PFOS 每 +1 SD 对各结局的效应（%） ==========\n")
R <- rbindlist(lapply(names(SUB), function(nm) rbindlist(lapply(YK, function(k) run(SUB[[nm]], "z_PFOS", k, nm)))))
R[, lab := LAB[outcome]]
print(dcast(R, outcome ~ tag, value.var = "est", fun.aggregate = mean)[order(outcome)], digits=1)

cat("\nN 与 95%CI：\n")
print(R[outcome %in% c("y_tt","y_shbg","y_ft"),
        .(tag, lab, N, est=round(est,1), lo=round(lo,1), hi=round(hi,1),
          p=format.pval(p, digits=2))][order(lab, tag)], nrows = 30)

cat("\n========== PFOA 同样分析 ==========\n")
R2 <- rbindlist(lapply(names(SUB), function(nm) rbindlist(lapply(YK, function(k) run(SUB[[nm]], "z_PFOA", k, nm)))))
R2[, lab := LAB[outcome]]
print(R2[outcome %in% c("y_tt","y_shbg","y_ft"),
         .(tag, lab, N, est=round(est,1), lo=round(lo,1), hi=round(hi,1),
           p=format.pval(p, digits=2))][order(lab, tag)], nrows = 30)

## ---------- 最小可检测效应（功效说明） --------------------------------------
cat("\n========== 功效：这些亚组能检出多大的效应？ ==========\n")
PW <- rbindlist(lapply(names(SUB), function(nm) {
  s <- SUB[[nm]][!is.na(y_tt) & w5 == 1 & !is.na(z_PFOS)]
  if (nrow(s) < 8) return(data.table(tag = nm, N = nrow(s), se = NA_real_, mde = NA_real_))
  se <- summary(lm(y_tt ~ z_PFOS + age + bmi_p, data = s))$coefficients["z_PFOS","Std. Error"]
  mde <- (qnorm(.975) + qnorm(.80)) * se
  cat(sprintf("%-30s N=%3d  最小可检出效应(80%%功效) = %+5.1f%%\n", nm, nrow(s), 100*(exp(mde)-1)))
  data.table(tag = nm, N = nrow(s), se = se, mde = 100*(exp(mde)-1))
}))

RR <- rbind(R[, marker := "PFOS"], R2[, marker := "PFOA"])
RR <- merge(RR, PW[, .(tag, mde)], by = "tag", all.x = TRUE)
fwrite(RR, file.path(RES, "subgroup_highsens.csv"))
fwrite(PW, file.path(RES, "subgroup_power.csv"))
cat("\n✅ 完成\n")
