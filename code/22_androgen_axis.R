## ============================================================================
## 22_androgen_axis.R — 雄激素轴：减重效应及其 PFAS 效应修饰
##   主样本：排除药物性去势/抑制（GnRH、抗雄、雌激素、糖皮质激素、阿片、酮康唑）
##   敏感性：全样本纳入
##   结局：总睾酮 / SHBG / 游离睾酮 / %游离 / FAI / 雌二醇 / T:E2 比值 / 白蛋白(对照)
##   检验：① 二次型交互 ② 三分位交互(2 df Wald, 不假设单调) ③ 连续剂量 ④ 单味同系物
##         ⑤ 配对聚类自助法（阻尼是否成比例 = "脱钩"）
## ============================================================================
suppressMessages({library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
A    <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
FL   <- readRDS(file.path(BASE, "data/derived/rx_flags.rds"))
OUT  <- file.path(BASE, "results")

crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X); bread <- solve(crossprod(X))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
wald <- function(b, V) as.numeric(t(b) %*% solve(V) %*% b)

## ---------- 派生 ------------------------------------------------------------
d0 <- FL[copy(A$S4), on = "SEQN"]
d0[is.na(castration), castration := 0L]; d0[is.na(suppress), suppress := 0L]
d0[, `:=`(
  tt_nmol = tt_ngdl / 28.84,
  e2_pmol = e2_pgml * 3.671
)]
d0[, `:=`(
  fai = tt_nmol / shbg * 100,
  te2 = tt_nmol / e2_pmol
)]
d0[, `:=`(
  y_tt = log(tt_nmol), y_shbg = log(shbg), y_ft = log(ft_nmol),
  y_pctfree = log(pct_free), y_fai = log(fai), y_e2 = log(e2_pmol),
  y_te2 = log(te2), y_alb = log(alb_gL)
)]
YK  <- c("y_tt","y_shbg","y_ft","y_pctfree","y_fai","y_e2","y_te2","y_alb")
LAB <- c(y_tt="Total testosterone", y_shbg="SHBG", y_ft="Free testosterone",
         y_pctfree="%Free testosterone", y_fai="FAI (T/SHBG)", y_e2="Estradiol",
         y_te2="T:E2 ratio", y_alb="Albumin (control)")

CONF <- A$CONF
prep <- function(dat, tag) {
  cf <- CONF[sapply(CONF, function(v) length(unique(dat[[v]])) > 1)]
  cc <- complete.cases(dat[, c("w5","pfas_z", YK, cf, "race_f","edu_f","smk"), with = FALSE])
  d  <- dat[cc]
  d[, `:=`(mv = as.numeric(pfas_z))]
  d[, mv2 := mv^2]
  d[, dw5 := dw1 / 5]
  d[, pfas_t := cut(mv, unique(quantile(mv, c(0,1/3,2/3,1))), labels = c("T1","T2","T3"), include.lowest = TRUE)]
  attr(d, "cf") <- cf; attr(d, "tag") <- tag
  d
}

## ---------- 拟合函数 --------------------------------------------------------
fit_quad <- function(d, k) {
  cf <- attr(d, "cf")
  f <- as.formula(paste(k, "~ w5*(mv+mv2) + mv + mv2 +", paste(cf, collapse=" + "), "+ race_f + edu_f + smk"))
  X <- model.matrix(f, data = d); cl <- factor(paste(d$strata, d$psu))
  m <- lm.fit(X, d[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  iw <- grep("^w5$", colnames(X)); i1 <- grep("^w5:mv$", colnames(X)); i2 <- grep("^w5:mv2$", colnames(X))
  q <- quantile(d$mv, c(.1,.5,.9)); L <- cbind(1,q,q^2)
  bt <- b[c(iw,i1,i2)]; Vs <- V[c(iw,i1,i2), c(iw,i1,i2)]
  th <- as.numeric(L %*% bt); se <- sqrt(diag(L %*% Vs %*% t(L)))
  data.table(outcome=k, label=LAB[k], N=nrow(d), n_trt=sum(d$w5),
             b1=b[i1], p_b1=2*pnorm(-abs(b[i1]/sqrt(V[i1,i1]))),
             b2=b[i2], p_b2=2*pnorm(-abs(b[i2]/sqrt(V[i2,i2]))),
             b1_std=b[i1]*sd(d$mv)/sd(d[[k]]),
             p_joint=pchisq(wald(b[c(i1,i2)], V[c(i1,i2),c(i1,i2)]), 2, lower.tail=FALSE),
             p10=th[1], se10=se[1], p50=th[2], se50=se[2], p90=th[3], se90=se[3],
             drop=th[3]-th[1])
}

fit_tert <- function(d, k) {
  cf <- attr(d, "cf")
  f <- as.formula(paste(k, "~ w5*pfas_t + pfas_t +", paste(cf, collapse=" + "), "+ race_f + edu_f + smk"))
  X <- model.matrix(f, data = d); cl <- factor(paste(d$strata, d$psu))
  m <- lm.fit(X, d[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  j <- grep("^w5:pfas_t", colnames(X))
  W <- wald(b[j], V[j,j])
  ## 各层减重效应
  est <- sapply(c("T1","T2","T3"), function(g) {
    b[grep("^w5$", colnames(X))] + if (g=="T1") 0 else b[grep(paste0("^w5:pfas_t", g, "$"), colnames(X))]
  })
  data.table(outcome=k, label=LAB[k], N=nrow(d), wald=W, df=length(j),
             p_tert=pchisq(W, length(j), lower.tail=FALSE),
             cate_T1=est[1], cate_T2=est[2], cate_T3=est[3])
}

## ---------- 主分析：清洁样本 -------------------------------------------------
d  <- prep(d0[castration == 0 & suppress == 0], "clean")
da <- prep(d0, "all")
cat("清洁样本 N =", nrow(d), " 处理 =", sum(d$w5), " 聚类 =", uniqueN(paste(d$strata, d$psu)), "\n")
cat("全样本   N =", nrow(da), " 处理 =", sum(da$w5), "\n\n")

RA  <- rbindlist(lapply(YK, function(k) fit_quad(d,  k)))
RAa <- rbindlist(lapply(YK, function(k) fit_quad(da, k)))

cat("========== A1. 清洁样本：二次型交互 ==========\n")
print(RA[, .(label, N, b1_std=sprintf("%+.3f", b1_std),
             `P(交互)`=format.pval(p_joint, digits=3),
             P10=sprintf("%+.3f", p10), P50=sprintf("%+.3f", p50),
             P90=sprintf("%+.3f", p90), `P90-P10`=sprintf("%+.3f", drop))], nrows=20)

cat("\n========== A2. 全样本（含药物影响）对比 ==========\n")
cmp <- merge(RA[, .(outcome, label, clean_b1=b1_std, clean_p=p_joint, clean_drop=drop)],
             RAa[, .(outcome, all_b1=b1_std, all_p=p_joint, all_drop=drop)], by="outcome")
print(cmp[, .(label, all_b1=sprintf("%+.3f", all_b1), all_p=format.pval(all_p, digits=3),
              clean_b1=sprintf("%+.3f", clean_b1), clean_p=format.pval(clean_p, digits=3),
              all_drop=sprintf("%+.3f", all_drop), clean_drop=sprintf("%+.3f", clean_drop))], nrows=20)

## 临床换算
g <- function(k, col) RA[outcome==k][[col]]
cat("\n临床换算（P10 -> P90，清洁样本）：\n")
for (k in c("y_tt","y_shbg","y_ft","y_pctfree","y_fai","y_e2","y_te2","y_alb")) {
  cat(sprintf("  %-22s %+6.1f%% -> %+6.1f%%\n", LAB[k], 100*(exp(g(k,"p10"))-1), 100*(exp(g(k,"p90"))-1)))
}
cat(sprintf("  总睾酮绝对值: %.1f -> %.1f nmol/L (%.0f -> %.0f ng/dL)\n",
            median(d$tt_nmol)*exp(g("y_tt","p10")), median(d$tt_nmol)*exp(g("y_tt","p90")),
            median(d$tt_nmol)*28.84*exp(g("y_tt","p10")), median(d$tt_nmol)*28.84*exp(g("y_tt","p90"))))

## ---------- B. 三分位交互（不假设单调） --------------------------------------
RT  <- rbindlist(lapply(YK, function(k) fit_tert(d,  k)))
RTa <- rbindlist(lapply(YK, function(k) fit_tert(da, k)))
cat("\n========== B. 三分位交互检验（2 df，不假设单调） ==========\n")
print(RT[, .(label, N, `P(交互)`=format.pval(p_tert, digits=3),
             T1=sprintf("%+.3f", cate_T1), T2=sprintf("%+.3f", cate_T2),
             T3=sprintf("%+.3f", cate_T3))], nrows=20)
cat("\n全样本同一检验：\n")
print(RTa[, .(label, `P(交互)`=format.pval(p_tert, digits=3),
              T1=sprintf("%+.3f", cate_T1), T2=sprintf("%+.3f", cate_T2),
              T3=sprintf("%+.3f", cate_T3))], nrows=20)

## ---------- C. 连续剂量 ------------------------------------------------------
cat("\n========== C. 每 5% 减重（清洁样本） ==========\n")
cfd <- attr(d, "cf")
Xc <- model.matrix(as.formula(paste("~ dw5*(mv+mv2) + mv + mv2 +",
                                    paste(cfd, collapse=" + "), "+ race_f + edu_f + smk")), data = d)
jD <- grep("^dw5$", colnames(Xc)); jI <- grep("^dw5:mv$", colnames(Xc)); jI2 <- grep("^dw5:mv2$", colnames(Xc))
clc <- factor(paste(d$strata, d$psu))
RC2 <- rbindlist(lapply(YK, function(k) {
  m <- lm.fit(Xc, d[[k]]); V <- crsV(Xc, m$residuals, clc)
  data.table(outcome=k, label=LAB[k], per5=m$coefficients[jD],
             p=2*pnorm(-abs(m$coefficients[jD]/sqrt(V[jD,jD]))),
             int=m$coefficients[jI], p_int=2*pnorm(-abs(m$coefficients[jI]/sqrt(V[jI,jI]))),
             pint_joint=pchisq(wald(m$coefficients[c(jI,jI2)], V[c(jI,jI2),c(jI,jI2)]), 2, lower.tail=FALSE))
}))
print(RC2[, .(label, `每5%`=sprintf("%+.4f", per5), p=format.pval(p, digits=3),
              `交互`=sprintf("%+.4f", int), `P(交互2df)`=format.pval(pint_joint, digits=3))], nrows=20)

## ---------- D. 配对自助法：阻尼是否成比例 ------------------------------------
set.seed(20261117)
B <- 800
Y <- as.matrix(d[, ..YK]); cid <- as.integer(clc); ncl <- nlevels(clc)
Xq <- model.matrix(as.formula(paste("~ w5*(mv+mv2) + mv + mv2 +", paste(cfd, collapse=" + "),
                                    "+ race_f + edu_f + smk")), data = d)
iq1 <- grep("^w5:mv$", colnames(Xq)); spl <- split(seq_along(cid), cid)
bt <- matrix(NA_real_, B, length(YK)); ok <- 0L
for (b in seq_len(B)) {
  ii <- unlist(spl[as.character(sample.int(ncl, ncl, replace=TRUE))], use.names=FALSE)
  res <- tryCatch(sapply(YK, function(k) {
    m <- lm.fit(Xq[ii, , drop=FALSE], Y[ii, k])
    if (any(is.na(m$coefficients))) NA_real_ else m$coefficients[iq1]
  }), error=function(e) rep(NA_real_, length(YK)))
  if (any(is.na(res))) next
  ok <- ok + 1L; bt[ok, ] <- res
}
bt <- bt[seq_len(ok), , drop=FALSE]
Z <- sapply(seq_along(YK), function(j) bt[, j] * sd(d$mv) / sd(d[[YK[j]]]))
Z <- as.data.table(matrix(Z, ncol = length(YK)))
setnames(Z, YK)
contr <- list(c("y_tt","y_shbg","TT - SHBG"), c("y_ft","y_tt","FT - TT"),
              c("y_fai","y_shbg","FAI - SHBG"), c("y_te2","y_e2","T:E2 - E2"),
              c("y_shbg","y_alb","SHBG - albumin"))
RC <- rbindlist(lapply(contr, function(z) {
  dd <- Z[[z[1]]] - Z[[z[2]]]
  data.table(contrast=z[3], est=mean(dd), lo=as.numeric(quantile(dd,.025)),
             hi=as.numeric(quantile(dd,.975)), p=2*min(mean(dd<=0), mean(dd>=0)))
}))
cat("\n========== D. 配对自助法（标准化 b1 之差，", ok, " 次） ==========\n")
print(RC[, .(contrast, diff=sprintf("%+.3f", est), `95%CI`=sprintf("%+.3f, %+.3f", lo, hi),
             p=format.pval(p, digits=3))])
cat(sprintf("\nCATE 比 (TT/SHBG): P10 = %.2f  P90 = %.2f\n",
            g("y_tt","p10")/g("y_shbg","p10"), g("y_tt","p90")/g("y_shbg","p90")))

## ---------- E. 单味同系物 ----------------------------------------------------
CONG <- c("PFOA","PFOS","PFHxS","PFNA")
for (v in CONG) d[, (paste0("z_",v)) := as.numeric(scale(log(get(v)))), by = cycle]
run_cong <- function(zvar, k) {
  dd <- copy(d)[!is.na(get(zvar))]; dd[, `:=`(mv=as.numeric(get(zvar)))]; dd[, mv2 := mv^2]
  f <- as.formula(paste(k, "~ w5*(mv+mv2) + mv + mv2 +", paste(cfd, collapse=" + "), "+ race_f + edu_f + smk"))
  X <- model.matrix(f, data = dd); cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  a <- grep("^w5:mv$", colnames(X)); bb <- grep("^w5:mv2$", colnames(X)); w <- grep("^w5$", colnames(X))
  qq <- quantile(dd$mv, c(.1,.9)); L <- cbind(1,qq,qq^2)
  Vs <- V[c(w,a,bb), c(w,a,bb)]; th <- as.numeric(L %*% b[c(w,a,bb)])
  se <- sqrt(diag(L %*% Vs %*% t(L)))
  data.table(congener=zvar, outcome=k, N=nrow(dd), b1_std=b[a]*sd(dd$mv)/sd(dd[[k]]),
             se_b1_std=sqrt(V[a,a])*sd(dd$mv)/sd(dd[[k]]),
             p_joint=pchisq(wald(b[c(a,bb)], V[c(a,bb),c(a,bb)]), 2, lower.tail=FALSE),
             drop=th[2]-th[1], se_drop=sqrt(se[1]^2+se[2]^2),
             p_drop=2*pnorm(-abs((th[2]-th[1])/sqrt(se[1]^2+se[2]^2))))
}
RD <- rbindlist(lapply(CONG, function(v) rbindlist(lapply(YK[1:7], function(k) run_cong(paste0("z_",v), k)))))
RD[, congener := sub("^z_","",congener)]
cat("\n========== E1. 单味同系物：标准化交互 b1 ==========\n")
print(dcast(RD, congener ~ outcome, value.var="b1_std")[order(congener)], digits=2)
cat("\n-- P90-P10 效应差 --\n"); print(dcast(RD, congener ~ outcome, value.var="drop")[order(congener)], digits=2)
cat("\n-- 联合检验 P --\n"); print(dcast(RD, congener ~ outcome, value.var="p_joint")[order(congener)], digits=3)

cat("\n========== E2. 互调模型 ==========\n")
for (v in CONG) d[, (paste0("mv_",v)) := as.numeric(scale(log(get(v)))), by = cycle]
fM <- as.formula(paste0("~ w5*(", paste0("mv_",CONG,collapse=" + "), ") + ",
                        paste0("mv_",CONG,collapse=" + "), " + ", paste(cfd, collapse=" + "),
                        " + race_f + edu_f + smk"))
XM <- model.matrix(fM, data=d); okm <- complete.cases(XM); XM <- XM[okm,]
clM <- factor(paste(d$strata, d$psu))[okm]
RM <- rbindlist(lapply(YK[1:7], function(k) {
  m <- lm.fit(XM, d[[k]][okm]); V <- crsV(XM, m$residuals, clM); j <- grep("^w5:mv_", colnames(XM))
  data.table(outcome=k, term=sub("^w5:","",colnames(XM)[j]),
             est=m$coefficients[j], p=2*pnorm(-abs(m$coefficients[j]/sqrt(diag(V)[j]))))
}))
print(RM[, .(outcome, term, est=sprintf("%+.4f", est), p=format.pval(p, digits=3))])

## ---------- F. 影响点：T3 估计对单个处理个体的敏感性 -------------------------
cat("\n========== F. 影响点分析（TT 的三分位效应） ==========\n")
for (nm in c("all","clean")) {
  dd <- if (nm=="all") da else d
  Xt <- model.matrix(as.formula(paste("y_tt ~ w5*pfas_t + pfas_t +", paste(attr(dd,"cf"), collapse=" + "),
                                      "+ race_f + edu_f + smk")), data=dd)
  m0 <- lm.fit(Xt, dd$y_tt)
  jj <- grep("^w5:pfas_tT3$", colnames(Xt)); jw <- grep("^w5$", colnames(Xt))
  base <- m0$coefficients[jw] + m0$coefficients[jj]
  t3 <- which(dd$pfas_t=="T3" & dd$w5==1)
  jk <- sapply(t3, function(i) {
    mm <- lm.fit(Xt[-i,], dd$y_tt[-i]); mm$coefficients[jw] + mm$coefficients[jj]
  })
  cat(sprintf("  %-6s T3 效应 = %+.3f；逐一剔除 T3 处理者后范围 %+.3f ~ %+.3f（%d 人）\n",
              nm, base, min(jk), max(jk), length(t3)))
}

fwrite(RA,  file.path(OUT, "androgen_axis.csv"))
fwrite(RAa, file.path(OUT, "androgen_axis_allsample.csv"))
fwrite(RT,  file.path(OUT, "androgen_tertile.csv"))
fwrite(RC2, file.path(OUT, "androgen_per5.csv"))
fwrite(RC,  file.path(OUT, "androgen_decouple.csv"))
fwrite(RD,  file.path(OUT, "androgen_congener.csv"))
fwrite(RM,  file.path(OUT, "androgen_congener_multi.csv"))
cat("\n✅ 完成\n")
