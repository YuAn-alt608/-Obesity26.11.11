## ============================================================================
## 30_liver_analysis.R — 肝-脂轴：减重效应及其 PFAS/PCB/OCP 效应修饰
##   结局：ALT / AST / GGT / ALP / 胆红素 / FLI / HSI / FIB-4 / AST:ALT
##   模型：① 二次型交互 ② 三分位交互(2 df) ③ 连续剂量(每 5%)
##   样本：8 周期，N≈12,737（处理 2,513）
## ============================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", "."); RES <- file.path(BASE, "results")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
d <- readRDS(file.path(BASE, "data/derived/liver.rds"))
CONF <- A$CONF

crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X); bread <- MASS::ginv(crossprod(X))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
## 子集内自动剔除常量/共线列，避免 solve() 奇异；但必须保留截距
## （V16 审计：初版 sd(v)>0 判据会删掉常数列 (Intercept)，使模型无截距拟合）
mkX <- function(f, dd) {
  X <- model.matrix(f, data = dd)
  ok <- apply(X, 2, function(v) all(is.finite(v)))
  vary <- apply(X, 2, function(v) sd(v) > 0)
  X[, ok & (vary | colnames(X) == "(Intercept)"), drop = FALSE]
}
wald <- function(b, V) as.numeric(t(b) %*% MASS::ginv(V) %*% b)

## ---------- 结局 -------------------------------------------------------------
d[, `:=`(
  y_alt = log(pmax(alt, .5)), y_ast = log(pmax(ast, .5)), y_ggt = log(pmax(ggt, .5)),
  y_alp = log(pmax(alp2, .5)), y_tb = log(pmax(tb, .05)),
  y_fli = log(pmax(fli, 0.1)), y_hsi = log(pmax(hsi, 0.1)),
  y_fib4 = log(pmax(fib4, 1e-3)), y_astalt = log(pmax(ast_alt, 1e-3))
)]
YK <- c("y_alt","y_ast","y_ggt","y_alp","y_tb","y_fli","y_hsi","y_fib4","y_astalt")
LAB <- c(y_alt="ALT", y_ast="AST", y_ggt="GGT", y_alp="ALP", y_tb="Bilirubin",
         y_fli="Fatty liver index", y_hsi="Hepatic steatosis index",
         y_fib4="FIB-4", y_astalt="AST:ALT ratio")

base <- d[ok_liver == TRUE]
base[, mv := as.numeric(pfas_z)]; base[, mv2 := mv^2]
base[, dw5 := dw1/5]
base[, pfas_t := cut(mv, unique(quantile(mv, c(0,1/3,2/3,1))), labels = c("T1","T2","T3"), include.lowest = TRUE)]
cf <- CONF[sapply(CONF, function(v) length(unique(base[[v]])) > 1)]
cat("N =", nrow(base), " 处理 =", sum(base$w5), " 聚类 =",
    uniqueN(paste(base$strata, base$psu)), " 协变量:", length(cf), "\n\n")

fml <- function(k, rhs) as.formula(paste(k, "~", rhs, "+", paste(cf, collapse=" + "), "+ race_f + edu_f + smk"))
fitq <- function(k) {
  dd <- base[!is.na(get(k))]
  X <- model.matrix(fml(k, "w5*(mv+mv2) + mv + mv2"), data = dd)
  cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  iw <- grep("^w5$", colnames(X)); i1 <- grep("^w5:mv$", colnames(X)); i2 <- grep("^w5:mv2$", colnames(X))
  q <- quantile(dd$mv, c(.1,.5,.9)); L <- cbind(1,q,q^2)
  bt <- b[c(iw,i1,i2)]; Vs <- V[c(iw,i1,i2), c(iw,i1,i2)]
  th <- as.numeric(L %*% bt); se <- sqrt(diag(L %*% Vs %*% t(L)))
  data.table(outcome=k, label=LAB[k], b1=b[i1], se_b1=sqrt(V[i1,i1]),
             p_joint=pchisq(wald(b[c(i1,i2)], V[c(i1,i2),c(i1,i2)]), 2, lower.tail=FALSE),
             eff0=b[iw], se_eff0=sqrt(V[iw,iw]),
             p10=th[1], p50=th[2], p90=th[3], se10=se[1], se50=se[2], se90=se[3])
}
fitp <- function(k) {  # 连续剂量
  dd <- base[!is.na(get(k))]
  X <- model.matrix(fml(k, "dw5*(mv+mv2) + mv + mv2"), data = dd)
  cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  jD <- grep("^dw5$", colnames(X)); j1 <- grep("^dw5:mv$", colnames(X)); j2 <- grep("^dw5:mv2$", colnames(X))
  data.table(outcome=k, label=LAB[k], per5=b[jD], p=2*pnorm(-abs(b[jD]/sqrt(V[jD,jD]))),
             p_int=pchisq(wald(b[c(j1,j2)], V[c(j1,j2),c(j1,j2)]), 2, lower.tail=FALSE))
}
fitt <- function(k) {  # 三分位
  dd <- base[!is.na(get(k))]
  X <- model.matrix(fml(k, "w5*pfas_t + pfas_t"), data = dd)
  cl <- factor(paste(dd$strata, dd$psu))
  m <- lm.fit(X, dd[[k]]); V <- crsV(X, m$residuals, cl); b <- m$coefficients
  j <- grep("^w5:pfas_t", colnames(X)); iw <- grep("^w5$", colnames(X))
  est <- c(b[iw], b[iw] + b[j[1]], b[iw] + b[j[2]])
  data.table(outcome=k, label=LAB[k], p_tert=pchisq(wald(b[j], V[j,j]), 2, lower.tail=FALSE),
             cate_T1=est[1], cate_T2=est[2], cate_T3=est[3])
}

RQ <- rbindlist(lapply(YK, fitq))
RP <- rbindlist(lapply(YK, fitp))
RT <- rbindlist(lapply(YK, fitt))

cat("========== A. 每 5% 减重与 PFAS 交互（N=", nrow(base), "） ==========\n")
m <- merge(RP, RQ[, .(outcome, p_joint)], by="outcome")
m[, `:=`(pct = 100*(exp(per5)-1))]
print(m[, .(label, `每5%`=sprintf("%+6.2f%%", pct), p=format.pval(p, digits=2),
            `交互P`=format.pval(p_int, digits=3))], nrows=12)

cat("\n========== B. 二次型交互：PFAS 三档的条件效应（%） ==========\n")
RQ[, `:=`(p10p=100*(exp(p10)-1), p50p=100*(exp(p50)-1), p90p=100*(exp(p90)-1))]
print(RQ[, .(label, `P10`=sprintf("%+6.2f", p10p), `P50`=sprintf("%+6.2f", p50p),
             `P90`=sprintf("%+6.2f", p90p), `交互P`=format.pval(p_joint, digits=3))], nrows=12)

cat("\n========== C. 三分位（不假设单调） ==========\n")
RT[, `:=`(T1p=100*(exp(cate_T1)-1), T2p=100*(exp(cate_T2)-1), T3p=100*(exp(cate_T3)-1))]
print(RT[, .(label, `T1`=sprintf("%+6.2f", T1p), `T2`=sprintf("%+6.2f", T2p),
             `T3`=sprintf("%+6.2f", T3p), `P`=format.pval(p_tert, digits=3))], nrows=12)

## ---------- D. 二分类临床终点（logistic） ------------------------------------
cat("\n========== D. 临床二分类终点 ==========\n")
base[, `:=`(nafld = as.integer(fli >= 60), fib_hi = as.integer(fib4 > 1.3),
            alt_hi = as.integer(alt > 40), ggt_hi = as.integer(ggt > 60))]
BIN <- c(nafld="FLI>=60 (steatosis)", fib_hi="FIB-4>1.3", alt_hi="ALT>40 U/L", ggt_hi="GGT>60 U/L")
RB <- rbindlist(lapply(names(BIN), function(k) {
  dd <- base[!is.na(get(k))]
  X <- mkX(fml(k, "w5*(mv+mv2) + mv + mv2"), dd)
  cl <- factor(paste(dd$strata, dd$psu))
  fit <- glm.fit(X, dd[[k]], family = binomial())
  ## logistic 的聚类稳健 sandwich：meat 用得分残差 X*u, u = (y-p)*p*(1-p)... 标准 GLM 形式为 X*(y-p)
  phat <- fit$fitted.values
  U <- X * (dd[[k]] - phat)
  V <- crsV(X, U, cl)          # crsV 内 bread = (X'X)^-1，对 logistic 需换成 (X'WX)^-1
  Wm <- diag(phat*(1-phat), nrow = length(phat))
  bread <- MASS::ginv(crossprod(X, Wm %*% X))
  idx <- split(seq_len(nrow(X)), cl); meat <- matrix(0, ncol(X), ncol(X))
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], (dd[[k]] - phat)[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); n <- nrow(X); kk <- ncol(X)
  V <- bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-kk))
  b <- fit$coefficients
  i1 <- grep("^w5:mv$", colnames(X)); i2 <- grep("^w5:mv2$", colnames(X)); iw <- grep("^w5$", colnames(X))
  data.table(outcome=k, label=BIN[k], prev=mean(dd[[k]], na.rm=TRUE),
             OR0=exp(b[iw]), p_int=pchisq(wald(b[c(i1,i2)], V[c(i1,i2),c(i1,i2)]), 2, lower.tail=FALSE))
}))
print(RB[, .(label, prevalence=round(100*prev,1), `OR@mean exposure`=round(OR0,3),
             `交互P`=format.pval(p_int, digits=3))])

## ---------- E. 分层：性别与基线 BMI ------------------------------------------
cat("\n========== E. 分层（每 5% × PFAS 交互） ==========\n")
SE <- rbindlist(lapply(c("female","bmi_cat"), function(sv) {
  if (sv == "bmi_cat") base[, bmi_cat := cut(bmi_p, c(0,25,30,100), labels=c("<25","25-30",">=30"))]
  rbindlist(lapply(if (sv=="female") c(0,1) else levels(base$bmi_cat), function(g) {
    dd <- if (sv=="female") base[female==g] else base[bmi_cat==g & !is.na(bmi_cat)]
    if (nrow(dd) < 400) return(NULL)
    rbindlist(lapply(c("y_alt","y_ggt","y_fli"), function(k) {
      dk <- dd[!is.na(get(k))]
      X <- mkX(fml(k, "dw5*(mv+mv2) + mv + mv2"), dk)
      jD <- grep("^dw5$", colnames(X)); j1 <- grep("^dw5:mv$", colnames(X))
      if (!length(jD) || !length(j1)) return(NULL)
      cl <- factor(paste(dk$strata, dk$psu))
      mm <- lm.fit(X, dk[[k]]); V <- crsV(X, mm$residuals, cl); b <- mm$coefficients
      data.table(strat=sv, group=as.character(g), outcome=k, N=nrow(dk), trt=sum(dk$w5),
                 per5=100*(exp(b[jD])-1), p_int=2*pnorm(-abs(b[j1]/sqrt(V[j1,j1]))))
    }))
  }))
}))
print(SE[, .(strat, group, outcome, N, trt, `每5%`=sprintf("%+6.2f%%", per5),
             `交互P`=format.pval(p_int, digits=3))])

fwrite(RQ, file.path(RES, "liver_quad.csv"))
fwrite(RP, file.path(RES, "liver_per5.csv"))
fwrite(RT, file.path(RES, "liver_tertile.csv"))
fwrite(RB, file.path(RES, "liver_binary.csv"))
fwrite(SE, file.path(RES, "liver_strat.csv"))
cat("\n✅ 完成\n")
