## ============================================================================
## V16_survey_weighted.R — 复杂抽样加权敏感性核查（审稿人第 1 号必问点）
## ----------------------------------------------------------------------------
## 问题：主分析为兼容因果森林采用"未加权 + 层×PSU 聚类稳健 SE"。
##       统计审稿人会问：不用 CDC 规定的 PFAS 子样本权重（WTSA2YR / WTSB2YR /
##       WTSPF2YR）与"多周期合并除以周期数"，是否造成全美代表性偏倚？
##
## 做法：
##   0. 构建 8 周期合并权重 wt_pfas_8cycle = wt_pf / 8
##      （wt_pf 由 06_build_pop.R 按周期取 CDC 的 PFAS 专项子样本权重）
##   1. svydesign(ids=~psu, strata=~strata, weights=~wt_pfas_8cycle, nest=TRUE) + svyglm
##      重跑三块核心结果：
##        A. 炎症域：10 结局主效应（每 5% 减重）+ WBC 连续/三档交互
##        B. 雄激素轴：8 结局主效应（22 的 S1 规格）
##        C. 肝-脂轴：9 结局主效应（30 的规格）
##   2. 逐项对比"未加权 CR1" vs "加权 design-based"
##
## 输出：verification/survey_weighted_sensitivity.csv
##       verification/survey_design_summary.csv
##       verification/survey_deff.csv
## ============================================================================
suppressMessages({library(data.table); library(survey); library(MASS)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
VER  <- file.path(BASE, "verification")
options(survey.lonely.psu = "adjust")   # 域内单 PSU 层：中心化调整（保守）
set.seed(20261117)

## ---------------- 0. 主框架与权重 -------------------------------------------
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
D <- readRDS(file.path(BASE, "data/derived/inflam.rds"))
L <- readRDS(file.path(BASE, "data/derived/liver.rds"))
F <- readRDS(file.path(BASE, "data/derived/rx_flags.rds"))
CONF <- A$CONF

M <- merge(D, L[, .(SEQN, alt, ast, ggt, alp2, tb, fli, hsi, fib4, ast_alt, ok_liver)],
           by = "SEQN", all.x = TRUE)
M <- merge(M, F[, .(SEQN, castration, suppress)], by = "SEQN", all.x = TRUE)
M[is.na(castration), castration := 0L]; M[is.na(suppress), suppress := 0L]

NCYC <- uniqueN(D[!is.na(wt_pf) & wt_pf > 0, cycle])
M[, wt_pfas_8cycle := wt_pf / NCYC]
cat("合并周期数 =", NCYC, " → wt_pfas_8cycle = wt_pf /", NCYC, "\n")

M[, `:=`(tt_nmol = tt_ngdl / 28.84, e2_pmol = e2_pgml * 3.671)]
M[, `:=`(fai = tt_nmol / shbg * 100, te2 = tt_nmol / e2_pmol)]
M[, `:=`(y_tt = log(tt_nmol), y_shbg = log(shbg), y_ft = log(ft_nmol),
         y_pctfree = log(pct_free), y_fai = log(fai), y_e2 = log(e2_pmol),
         y_te2 = log(te2), y_alb = log(alb_gL))]
M[, `:=`(y_alt = log(pmax(alt, .5)), y_ast = log(pmax(ast, .5)),
         y_ggt = log(pmax(ggt, .5)), y_alp = log(pmax(alp2, .5)),
         y_tb = log(pmax(tb, .05)), y_fli = log(pmax(fli, 0.1)),
         y_hsi = log(pmax(hsi, 0.1)), y_fib4 = log(pmax(fib4, 1e-3)),
         y_astalt = log(pmax(ast_alt, 1e-3)))]
M[, `:=`(mv = as.numeric(pfas_z), dw5 = dw1 / 5)]
M[, mv2 := mv^2]

## 设计框架 = 所有有 PFAS 子样本权重的个体（子集化时保留完整层/PSU）
FR <- M[!is.na(wt_pfas_8cycle) & wt_pfas_8cycle > 0]
## 三分位（未加权切点，与 36 一致）与加权切点（补充口径）
FR[, pfas_t := cut(mv, unique(quantile(mv, c(0, 1/3, 2/3, 1), na.rm = TRUE)),
                   labels = c("T1","T2","T3"), include.lowest = TRUE)]
wq <- function(x, w, p) { o <- order(x); x <- x[o]; w <- w[o]
  cs <- cumsum(w) / sum(w); sapply(p, function(q) x[which(cs >= q)[1]]) }
okq <- !is.na(FR$mv) & !is.na(FR$wt_pfas_8cycle)
bw  <- unique(wq(FR$mv[okq], FR$wt_pfas_8cycle[okq], c(0, 1/3, 2/3, 1)))
stopifnot(length(bw) == 4)
FR[, pfas_tw := cut(mv, bw, labels = c("T1","T2","T3"), include.lowest = TRUE)]
cat("加权三分位切点：", paste(signif(bw, 4), collapse = " | "), "\n")
cat("未加权切点(36 口径)：", paste(signif(unique(quantile(FR$mv, c(0,1/3,2/3,1), na.rm=TRUE)), 4),
    collapse = " | "), "\n")
FR <- as.data.frame(FR)
des <- svydesign(ids = ~psu, strata = ~strata, weights = ~wt_pfas_8cycle,
                 nest = TRUE, data = FR, deff = TRUE)
nF <- nrow(FR)

cat("设计框架 N =", nF, " 层 =", length(unique(FR$strata)),
    " PSU 组合 =", length(unique(paste(FR$strata, FR$psu))), "\n")
cat("加权人口估计 =", format(round(sum(FR$wt_pfas_8cycle)/1e6, 1), big.mark = ","),
    "百万（8 周期平均）；加权与未加权样本量比 =",
    round(sum(FR$wt_pfas_8cycle)/nF, 1), "\n\n")

DS <- data.table(
  item = c("合并周期数", "设计框架 N", "层数 (SDMVSTRA)", "PSU 组合数",
           "加权人口估计（百万）", "权重均值", "权重中位数", "权重 CV",
           "设计自由度 (degf)", "域内单 PSU 层处理"),
  value = c(NCYC, nF, length(unique(FR$strata)),
            length(unique(paste(FR$strata, FR$psu))),
            round(sum(FR$wt_pfas_8cycle)/1e6, 1),
            round(mean(FR$wt_pfas_8cycle), 0),
            round(median(FR$wt_pfas_8cycle), 0),
            round(sd(FR$wt_pfas_8cycle)/mean(FR$wt_pfas_8cycle), 3),
            degf(des), "survey.lonely.psu = adjust"))
fwrite(DS, file.path(VER, "survey_design_summary.csv"))

## ---------------- 1. 工具函数 -----------------------------------------------
crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X)
  bread <- tryCatch(solve(crossprod(X)), error = function(e) ginv(crossprod(X)))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx)
  bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
mkX <- function(f, dd) {
  X <- model.matrix(f, data = dd)
  ok <- apply(X, 2, function(v) all(is.finite(v)))
  vary <- apply(X, 2, function(v) sd(v) > 0)
  X[, ok & (vary | colnames(X) == "(Intercept)"), drop = FALSE]
}
## 所有模型共用同一组协变量；子集必须同时滤掉协变量缺失，
## 否则 model.matrix 的 na.action 会静默删行，与 lm.fit 的 y 长度不一致
SELCOV <- unique(c(CONF, "race_f", "edu_f", "smk"))
mksel <- function(vars, extra = rep(TRUE, nF)) {
  ok <- rep(TRUE, nF)
  for (v in unique(c(vars, SELCOV))) ok <- ok & !is.na(FR[[v]])
  ok & extra
}

OUT <- list()
add <- function(...) OUT[[length(OUT) + 1]] <<- data.table(...)

cmp <- function(domain, outcome, label, term, rhs, yk, sel, pat) {
  f <- as.formula(paste(yk, "~", rhs))
  ## 未加权 CR1
  dd <- FR[sel, ]
  X <- mkX(f, dd); m <- lm.fit(X, dd[[yk]])
  cl <- factor(paste(dd$strata, dd$psu)); V <- crsV(X, m$residuals, cl)
  j <- grep(pat, colnames(X))
  ub <- if (length(j) == 1) unname(m$coefficients[j]) else NA_real_
  us <- if (length(j) == 1) sqrt(V[j, j]) else NA_real_
  up <- if (length(j) == 1) 2 * pnorm(-abs(ub/us)) else NA_real_
  ## 加权 design-based
  sub <- des[sel, ]
  mw <- tryCatch(svyglm(f, design = sub), error = function(e) NULL)
  if (is.null(mw)) { wb <- ws <- wp <- NA_real_ } else {
    cf <- summary(mw)$coefficients; jj <- grep(pat, rownames(cf))
    if (length(jj) == 1) { wb <- cf[jj,1]; ws <- cf[jj,2]; wp <- cf[jj,4] }
    else { wb <- ws <- wp <- NA_real_ }
  }
  add(domain = domain, outcome = outcome, label = label, term = term,
      N_unw = nrow(dd), N_w = sum(sel),
      beta_unw = ub, se_unw = us, p_unw = up,
      beta_wt = wb, se_wt = ws, p_wt = wp,
      ratio_wt_unw = ifelse(is.na(ub) | ub == 0, NA, wb/ub))
}

## ---------------- 2. A 域：炎症 ---------------------------------------------
cat("========== A. 炎症域 ==========\n")
selA <- mksel(c("w5", "mv", "dw5"), FR$ok_infl %in% TRUE)
cf <- CONF[sapply(CONF, function(v) length(unique(FR[selA, v])) > 1)]
CFS  <- paste(cf, collapse = " + ")
## 22 / 30 的规格额外纳入 race_f / edu_f / smk；34 / 36 没有 —— 分别照搬以严格对齐
CFS2 <- paste(c(cf, "race_f", "edu_f", "smk"), collapse = " + ")
cat("分析样本 N =", sum(selA), " 处理 =", sum(FR$w5[selA]), " 协变量 =", length(cf), "\n")

YK <- c(y_crp="CRP", y_wbc="WBC", y_neut="Neutrophil", y_lym="Lymphocyte",
        y_mono="Monocyte count", y_monop="Monocyte %", y_nlr="NLR",
        y_mlr="MLR", y_plr="PLR", y_sii="SII")
for (k in names(YK))
  cmp("inflammation", k, YK[[k]], "main_per5",
      paste("dw5 + w5 +", CFS), k,
      mksel(c(k, "w5", "mv", "dw5"), FR$ok_infl %in% TRUE), "^dw5$")

## WBC / 中性粒 交互（连续 DiD，36 的 S1b）
for (k in c("y_wbc", "y_neut")) {
  lab <- YK[[k]]
  cmp("inflammation", k, lab, "int_cont",
      paste("w5*mv +", CFS), k, mksel(c(k, "w5", "mv"), FR$ok_infl %in% TRUE),
      "^w5:mv$|^mv:w5$")
}
## WBC 三档交互（未加权切点 / 加权切点）
cmp("inflammation", "y_wbc", "WBC [unweighted cut]", "int_tert_T3vT1",
    paste("w5*pfas_t +", CFS), "y_wbc",
    mksel(c("y_wbc", "w5", "mv", "pfas_t"), FR$ok_infl %in% TRUE),
    "^w5:pfas_tT3$|^pfas_tT3:w5$")
cmp("inflammation", "y_wbc", "WBC [weighted cut]", "int_tert_T3vT1",
    paste("w5*pfas_tw +", CFS), "y_wbc",
    mksel(c("y_wbc", "w5", "mv", "pfas_tw"), FR$ok_infl %in% TRUE),
    "^w5:pfas_twT3$|^pfas_twT3:w5$")

## ---------------- 3. B 域：雄激素轴 -----------------------------------------
cat("\n========== B. 雄激素轴 ==========\n")
## 与 22 完全同一样本：A$S4 是 2013-2016 睾酮子样本（H/I），N=1,337
S4ID <- A$S4$SEQN
extraB <- FR$castration == 0 & FR$suppress == 0 & FR$SEQN %in% S4ID
selB <- mksel(c("w5", "mv", "tt_ngdl"), extraB)
cat("清洁样本 N =", sum(selB), " 处理 =", sum(FR$w5[selB]), "\n")
YB <- c(y_tt="Total testosterone", y_shbg="SHBG", y_ft="Free testosterone",
        y_pctfree="%Free testosterone", y_fai="FAI (T/SHBG)", y_e2="Estradiol",
        y_te2="T:E2 ratio", y_alb="Albumin (control)")
for (k in names(YB))
  cmp("androgen", k, YB[[k]], "main_per5",
      paste("dw5*(mv+mv2) + mv + mv2 +", CFS2),
      k, mksel(c(k, "w5", "mv", "tt_ngdl"), extraB), "^dw5$")

## ---------------- 4. C 域：肝-脂轴 ------------------------------------------
cat("\n========== C. 肝-脂轴 ==========\n")
selC <- mksel(c("w5", "mv"), FR$ok_liver %in% TRUE)
cat("样本 N =", sum(selC), " 处理 =", sum(FR$w5[selC]), "\n")
YL <- c(y_alt="ALT", y_ast="AST", y_ggt="GGT", y_alp="ALP", y_tb="Bilirubin",
        y_fli="Fatty liver index", y_hsi="Hepatic steatosis index",
        y_fib4="FIB-4", y_astalt="AST:ALT ratio")
for (k in names(YL))
  cmp("liver", k, YL[[k]], "main_per5",
      paste("dw5*(mv+mv2) + mv + mv2 +", CFS2),
      k, mksel(c(k, "w5", "mv"), FR$ok_liver %in% TRUE), "^dw5$")

## ---------------- 5. 设计效应（DEFF）---------------------------------------
## 手算 Kish 设计效应：DEFF = Var_design(均值) / (s^2 / n)
DEFF <- rbindlist(lapply(c("y_wbc", "y_crp", "y_fli", "y_shbg", "w5", "mv"),
  function(k) {
    sv <- tryCatch(svymean(as.formula(paste0("~", k)), des, na.rm = TRUE),
                   error = function(e) NULL)
    if (is.null(sv)) return(NULL)
    x <- FR[[k]]; ok <- !is.na(x); w <- FR$wt_pfas_8cycle[ok]; x <- x[ok]
    n <- length(x)
    v_srs <- var(x) / n
    data.table(variable = k, n = n,
               mean_unw = mean(x),
               wmean    = as.numeric(coef(sv)),
               se_design = as.numeric(SE(sv)),
               se_srs = sqrt(v_srs),
               deff = as.numeric(SE(sv))^2 / v_srs,
               n_eff = round(n * v_srs / as.numeric(SE(sv))^2),
               kish = sum(w)^2 / sum(w^2) / n)
  }))
cat("\n========== 5. 设计效应 ==========\n"); print(DEFF, row.names = FALSE)
fwrite(DEFF, file.path(VER, "survey_deff.csv"))

## ---------------- 6. 汇总 --------------------------------------------------
S <- rbindlist(OUT)
S[, dir := fifelse(is.na(beta_unw) | is.na(beta_wt), "NA",
            fifelse(sign(beta_unw) == sign(beta_wt), "一致", "翻转"))]
S[, `:=`(sig_u = !is.na(p_unw) & p_unw < 0.05,
         sig_w = !is.na(p_wt)  & p_wt  < 0.05)]
S[, verdict := fifelse(dir == "NA", "无法估计",
               fifelse(dir == "翻转", "方向翻转（需复核）",
                fifelse(sig_u & sig_w, "加权后仍显著（通过）",
                 fifelse(!sig_u & !sig_w, "两者均不显著（一致）",
                  "显著性口径改变（需披露）"))))]
fwrite(S, file.path(VER, "survey_weighted_sensitivity.csv"))

cat("\n========== 主效应：未加权 vs 加权（每 5% 减重）==========\n")
Pr <- S[term == "main_per5"]
print(Pr[, .(domain, label, N = N_unw,
             unw_pct = sprintf("%+7.2f", 100*(exp(beta_unw)-1)),
             p_unw = format.pval(p_unw, digits = 2),
             wt_pct  = sprintf("%+7.2f", 100*(exp(beta_wt)-1)),
             p_wt = format.pval(p_wt, digits = 2),
             ratio = round(ratio_wt_unw, 2), verdict = verdict)], row.names = FALSE)

cat("\n========== WBC / 中性粒 交互：未加权 vs 加权 ==========\n")
Pi <- S[grepl("^int_", term)]
print(Pi[, .(label, term, N = N_unw,
             b_unw = sprintf("%+.5f", beta_unw),
             p_unw = format.pval(p_unw, digits = 2),
             b_wt  = sprintf("%+.5f", beta_wt),
             p_wt = format.pval(p_wt, digits = 2),
             ratio = round(ratio_wt_unw, 2), verdict = verdict)], row.names = FALSE)

cat("\n========== 判定计数 ==========\n")
print(S[, .N, by = .(domain, verdict)][order(domain, verdict)], row.names = FALSE)
cat("\n方向翻转：", S[dir == "翻转", .N], "/", nrow(S), "\n")
cat("✅ V16 完成 → verification/survey_weighted_sensitivity.csv\n")
