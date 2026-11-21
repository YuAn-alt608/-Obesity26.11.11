## ============================================================================
## 31_liver_decomp.R — HSI/FLI 交互的来源分解
##   疑问：HSI = 8*(ALT/AST) + BMI + 2*female + 2*diab，本身含 BMI；
##         若 PFAS 修饰的是"减重幅度"而非"肝脏脂肪"，则 HSI 的交互是假象
##   做法：① 直接检验 当前BMI / 体重变化 × PFAS
##         ② 把 HSI、FLI 拆成"肝指标部分"与"人体测量部分"
##         ③ 对 HSI 做不含 BMI 的改造版 HSI_liver = 8*(ALT/AST)
## ============================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", "."); RES <- file.path(BASE, "results")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
d <- readRDS(file.path(BASE, "data/derived/liver.rds"))
CONF <- A$CONF

crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X); bread <- solve(crossprod(X))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
wald <- function(b, V) as.numeric(t(b) %*% solve(V) %*% b)

base <- d[ok_liver == TRUE]
base[, `:=`(mv = as.numeric(pfas_z), dw5 = dw1/5)]
base[, mv2 := mv^2]
base[, `:=`(
  hsi_liver = 8*(alt/ast),                    # 去掉 BMI / 性别 / 糖尿病
  fli_liver = 0.953*log(tg2) + 0.718*log(ggt) + 0.053*waist,  # 去掉 BMI 项
  dbmi = bmi - bmi_p,                         # 体重变化（kg/m2）
  dwt  = dw1
)]
YK <- c("bmi","dbmi","hsi","hsi_liver","fli","fli_liver","waist","y_alt2")
base[, y_alt2 := log(pmax(alt,.5))]
LAB <- c(bmi="Current BMI", dbmi="BMI change (current - 1y prior)", hsi="HSI (with BMI)",
         hsi_liver="HSI liver part 8*(ALT/AST)", fli="FLI (with BMI)",
         fli_liver="FLI liver part (TG+GGT+waist)", waist="Waist circumference",
         y_alt2="log ALT")

fml <- function(k, rhs) as.formula(paste(k, "~", rhs, "+", paste(CONF, collapse=" + "), "+ race_f + edu_f + smk"))
out <- rbindlist(lapply(YK, function(k) {
  dd <- base[!is.na(get(k))]
  X <- model.matrix(fml(k, "w5*(mv+mv2) + mv + mv2"), data = dd)
  cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  iw <- grep("^w5$", colnames(X)); i1 <- grep("^w5:mv$", colnames(X)); i2 <- grep("^w5:mv2$", colnames(X))
  q <- quantile(dd$mv, c(.1,.5,.9)); L <- cbind(1,q,q^2)
  th <- as.numeric(L %*% b[c(iw,i1,i2)])
  sdy <- sd(dd[[k]])
  data.table(term=k, label=LAB[k], N=nrow(dd),
             effect0=b[iw], sd_y=sdy, effect0_sd=b[iw]/sdy,
             p10=th[1], p50=th[2], p90=th[3],
             p_joint=pchisq(wald(b[c(i1,i2)], V[c(i1,i2),c(i1,i2)]), 2, lower.tail=FALSE))
}))
cat("========== 效应修饰的来源分解（每 5% 减重？否，此处为 w5 二分类） ==========\n")
print(out[, .(label, N, `效应@均值`=round(effect0,3), `标准化`=round(effect0_sd,3),
              P10=round(p10,3), P90=round(p90,3), `交互P`=format.pval(p_joint, digits=3))])

## ---------- 关键诊断：PFAS 是否修饰"体重/腰围"本身 ---------------------------
cat("\n========== 关键诊断：PFAS 是否修饰减重幅度 ==========\n")
diag <- rbindlist(lapply(c("dbmi","dwt","waist"), function(k) {
  dd <- base[!is.na(get(k))]
  X <- model.matrix(fml(k, "dw5*(mv+mv2) + mv + mv2"), data = dd)
  cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  jD <- grep("^dw5$", colnames(X)); j1 <- grep("^dw5:mv$", colnames(X)); j2 <- grep("^dw5:mv2$", colnames(X))
  data.table(term=k, per5=b[jD], p_int=pchisq(wald(b[c(j1,j2)], V[c(j1,j2),c(j1,j2)]), 2, lower.tail=FALSE),
             int_lin=b[j1], p_lin=2*pnorm(-abs(b[j1]/sqrt(V[j1,j1]))))
}))
print(diag[, .(term, `每5%`=round(per5,3), `联合交互P`=format.pval(p_int, digits=3),
               `线性交互`=round(int_lin,3), `P`=format.pval(p_lin, digits=3))])

## ---------- 结论性判断：把 BMI 从结局里"扣掉"后交互是否还在 ------------------
cat("\n========== 判定 ==========\n")
a <- out[term=="hsi"]; b2 <- out[term=="hsi_liver"]; c1 <- out[term=="bmi"]; c2 <- out[term=="dbmi"]
cat(sprintf("HSI（含BMI）      交互P = %.2g\n", a$p_joint))
cat(sprintf("HSI 肝脏部分      交互P = %.2g\n", b2$p_joint))
cat(sprintf("当前 BMI          交互P = %.2g\n", c1$p_joint))
cat(sprintf("BMI 变化          交互P = %.2g\n", c2$p_joint))

fwrite(out, file.path(RES, "liver_decomp.csv"))
fwrite(diag, file.path(RES, "liver_weight_diag.csv"))
cat("\n✅ 完成\n")
