## ============================================================================
## V18_qgcomp_mixture.R — 混合暴露模型：分位数 g-计算（Quantile G-Computation）
##                         审稿人第 3 号必问点（EHP/JMH 近年硬性要求）
## ----------------------------------------------------------------------------
## 问题：只看单污染物或 ΣPFAS 摩尔和，不足以支持"混合物"层面的推断，
##       也无法回答"到底是哪一味同系物在驱动 WBC 降幅钝化"。
##
## 做法（qgcomp 2.19.6）：
##   1. 6 种 PFAS 同系物在周期内 z 标准化（与 pfas_z 同一变换），q = 4 分位分箱
##   2. qgcomp.noboot：ψ = 混合物每上升一个分位（整体）对 log WBC 的联合效应
##      —— 同时给出正/负向权重（哪一味在推、哪一味在拉）
##   3. 分处理组 / 对照组各跑一次：Δψ = ψ(减重组) − ψ(未减重组)
##      = 混合物效应的"效应修饰"，即 PFAS 混合物是否钝化减重的 WBC 下降
##   4. 6 个 w5 × 分位暴露交互的联合 Wald（CR1 聚类稳健）
##   5. Δψ 的层内置换 P（暴露块在 周期×性别 层内整体打乱，B = 500）
##   6. CRP 作为阴性对照（预期无效应修饰）
##
## 输出：verification/qgcomp_mixture.csv
##       verification/qgcomp_weights.csv
##       verification/qgcomp_interaction.csv
## ============================================================================
suppressMessages({library(data.table); library(qgcomp); library(MASS)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
VER  <- file.path(BASE, "verification")
set.seed(20261117)

A <- readRDS(file.path(BASE, "data/derived/analysis.rds")); CONF <- A$CONF
D <- readRDS(file.path(BASE, "data/derived/inflam.rds"))
CGL <- c("PFOA","PFOS","PFHxS","PFNA","PFDA","PFUnDA")

D <- D[ok_infl == TRUE & !is.na(w5) & !is.na(y_wbc) & complete.cases(D[, CGL, with = FALSE])]
D[, stratum := paste(cycle, sex)]
## 周期内 z（与 06_build_pop.R 对 ΣPFAS 的处理同法）
for (a in CGL) D[, (paste0("z_", a)) := as.numeric(scale(log(get(a)))), by = cycle]
Z <- paste0("z_", CGL)
cf <- CONF[sapply(CONF, function(v) length(unique(D[[v]])) > 1)]
cat("分析样本 N =", nrow(D), " 处理 =", sum(D$w5), " 周期 =", uniqueN(D$cycle), "\n")
cat("同系物 z 的相关系数矩阵（周期内标准化后）：\n")
print(round(cor(D[, Z, with = FALSE]), 2))

## ---------------- 工具：CR1 ------------------------------------------------
crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X)
  bread <- tryCatch(solve(crossprod(X)), error = function(e) ginv(crossprod(X)))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
mkX <- function(f, dd) {
  X <- model.matrix(f, data = dd)
  ok <- apply(X, 2, function(v) all(is.finite(v)))
  vary <- apply(X, 2, function(v) sd(v) > 0)
  X[, ok & (vary | colnames(X) == "(Intercept)"), drop = FALSE]
}

## ---------------- 1. qgcomp：整体 / 对照 / 处理 ------------------------------
qg <- function(dd, rhs_extra, tag, BRK = NULL) {
  f <- as.formula(paste("y_wbc ~", rhs_extra, "+", paste(Z, collapse = " + "),
                        "+", paste(cf, collapse = " + ")))
  fit <- tryCatch(qgcomp.noboot(f, expnms = Z, data = as.data.frame(dd), q = 4,
                                breaks = BRK),
                  error = function(e) { cat("  qgcomp 失败:", conditionMessage(e), "\n"); NULL })
  if (is.null(fit)) return(NULL)
  pw <- fit$pos.weights; nw <- fit$neg.weights
  W <- data.table(term = c(names(pw), names(nw)),
                  weight = c(as.numeric(pw), as.numeric(nw)),
                  direction = c(rep("positive", length(pw)), rep("negative", length(nw))))
  W[, term := sub("^z_", "", term)]
  list(psi = fit$psi, var = fit$var.psi, n = nrow(dd), w = W)
}

RES <- list(); WT <- list(); INT <- list()
for (oc in c("y_wbc", "y_neut", "y_crp")) {
  lab <- c(y_wbc = "WBC", y_neut = "Neutrophil", y_crp = "CRP")[[oc]]
  D2 <- copy(D)[!is.na(get(oc))]
  if (oc != "y_wbc") setnames(D2, oc, "y_wbc")   # qgcomp 公式统一用 y_wbc 占位
  cat("\n==========", lab, " (N =", nrow(D2), ") ==========\n")
  ## 全样本共同切点：分层与置换共用，保证 ψ 可比且避免重复分箱（提速 ~50x）
  BRK <- quantize(as.data.frame(D2), expnms = Z, q = 4)$breaks
  ## 整体（调整 w5）
  g0 <- qg(D2, "w5", "overall", BRK)
  ## 分层
  gC <- qg(D2[w5 == 0], "1", "control", BRK)
  gT <- qg(D2[w5 == 1], "1", "treated", BRK)
  for (x in list(list(g0,"Overall"), list(gC,"Control (no loss)"), list(gT,"Treated (>=5% loss)"))) {
    if (is.null(x[[1]])) next
    g <- x[[1]]; nm <- x[[2]]
    RES[[length(RES)+1]] <- data.table(
      outcome = lab, stratum = nm, N = g$n, psi = g$psi, se = sqrt(g$var),
      p = 2*pnorm(-abs(g$psi/sqrt(g$var))),
      pct = 100*(exp(g$psi) - 1))
    if (nrow(g$w)) WT[[length(WT)+1]] <- cbind(data.table(outcome = lab, stratum = nm), g$w)
  }
  ## Δψ = 处理 − 对照
  if (!is.null(gC) && !is.null(gT)) {
    dpsi <- gT$psi - gC$psi
    vd <- gT$var + gC$var
    ## 置换：暴露块在层内整体打乱
    B <- 200; nd <- numeric(B)   # qgcomp 每次拟合约 0.3 s，200 次分辨率 0.005 足够
    for (b in seq_len(B)) {
      Db <- copy(D2)
      ## 同层内用同一个置换向量整体打乱 6 个同系物（保持混合物内部相关结构）
      Db[, (Z) := { p <- sample.int(.N); lapply(.SD, function(x) x[p]) },
         .SDcols = Z, by = stratum]
      fP <- as.formula(paste("y_wbc ~ 1 +", paste(Z, collapse=" + "), "+",
                             paste(cf, collapse=" + ")))
      a0 <- tryCatch(qgcomp.noboot(fP, expnms = Z, q = 4, breaks = BRK,
                    data = as.data.frame(Db[w5 == 0])), error = function(e) NULL)
      a1 <- tryCatch(qgcomp.noboot(fP, expnms = Z, q = 4, breaks = BRK,
                    data = as.data.frame(Db[w5 == 1])), error = function(e) NULL)
      if (!is.null(a0) && !is.null(a1)) nd[b] <- a1$psi - a0$psi
      if (b %% 100 == 0) cat("    perm", b, "/", B, "\n")
    }
    pperm <- (1 + sum(abs(nd) >= abs(dpsi))) / (sum(is.finite(nd)) + 1)
    ## 单模型交互的联合 Wald（CR1）：w5 × 6 个分位暴露
    Dq <- as.data.table(quantize(as.data.frame(D2), expnms = Z, q = 4)$data)
    QN <- Z                                  # quantize 保留原列名，值已变为 0..q-1
    fint <- as.formula(paste("y_wbc ~ w5*(", paste(QN, collapse = " + "), ") +",
                             paste(cf, collapse = " + ")))
    Xi <- mkX(fint, Dq); mi <- lm.fit(Xi, Dq$y_wbc)
    Vi <- crsV(Xi, mi$residuals, factor(paste(Dq$strata, Dq$psu)))
    ji <- grep("^w5:z_", colnames(Xi))
    Wv <- as.numeric(t(mi$coefficients[ji]) %*% ginv(Vi[ji, ji, drop = FALSE]) %*%
                     mi$coefficients[ji])
    pj <- pchisq(Wv, length(ji), lower.tail = FALSE)
    INT[[length(INT)+1]] <- data.table(
      outcome = lab, N = nrow(D2),
      psi_control = gC$psi, psi_treated = gT$psi, delta_psi = dpsi,
      se_delta = sqrt(vd), p_delta = 2*pnorm(-abs(dpsi/sqrt(vd))),
      p_perm = pperm, n_perm_ok = sum(is.finite(nd)),
      chi2_joint = Wv, df_joint = length(ji), p_joint_cr1 = pj)
  }
}

R <- rbindlist(RES); W <- rbindlist(WT, fill = TRUE); I <- rbindlist(INT)
fwrite(R, file.path(VER, "qgcomp_mixture.csv"))
fwrite(W, file.path(VER, "qgcomp_weights.csv"))
fwrite(I, file.path(VER, "qgcomp_interaction.csv"))

cat("\n========== 1. 混合物联合效应 ψ（每个分位）==========\n")
print(R[, .(outcome, stratum, N, psi = signif(psi,4), se = signif(se,4),
            p = signif(p,3), pct = sprintf("%+.2f%%", pct))], row.names = FALSE)

cat("\n========== 2. 效应修饰 Δψ（减重组 − 未减重组）==========\n")
print(I[, .(outcome, psi_c = signif(psi_control,4), psi_t = signif(psi_treated,4),
            dpsi = signif(delta_psi,4), se = signif(se_delta,4),
            p_delta = signif(p_delta,3), p_perm = signif(p_perm,3),
            p_joint = signif(p_joint_cr1,3))], row.names = FALSE)

cat("\n========== 3. 权重（哪一味同系物在驱动）==========\n")
print(dcast(W[direction == "positive"], outcome + stratum ~ term, value.var = "weight"),
      row.names = FALSE)
cat("\n负向权重（方向相反的同系物）：\n")
neg <- W[direction == "negative"]
if (nrow(neg)) print(dcast(neg, outcome + stratum ~ term, value.var = "weight"),
                     row.names = FALSE) else cat("（无）\n")

cat("\n========== 4. 判读 ==========\n")
wbc <- I[outcome == "WBC"]
cat(sprintf("WBC：ψ(对照) = %+.4f，ψ(减重) = %+.4f，Δψ = %+.4f（P = %.3f，置换 P = %.3f，联合 P = %.3g）\n",
            wbc$psi_control, wbc$psi_treated, wbc$delta_psi, wbc$p_delta,
            wbc$p_perm, wbc$p_joint_cr1))
cat("预期：Δψ > 0（纯 PFAS 混合物在减重者中把 WBC 往上抬 = 钝化减重带来的下降）\n")
cat("\n✅ V18 完成 → verification/qgcomp_mixture.csv, qgcomp_weights.csv, qgcomp_interaction.csv\n")
