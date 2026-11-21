## ============================================================================
## V17_lod_audit.R — 检出限（LOD / 左删失）核查与替换值敏感性（审稿人第 2 号必问点）
## ----------------------------------------------------------------------------
## 问题：NHANES 2003-2020 跨 20 年，PFAS 各同系物的 LOD 随周期变化。
##       若 WBC × PFAS 交互由不同周期的 LOD/√2 替换值人为制造，结论不成立。
##
## 关键事实（本脚本从原始 XPT 逐一核实）：
##   · LBD***L / SD***L 注释码：0 = ≥LOD，1 = <LOD
##   · NHANES 对 <LOD 记录填报 LOD/√2（如 0.1 = 0.1414/√2）
##   · 同一周期内可能并存两个 LOD（如 PFDA 在 2003-04 同时有 0.1 与 0.2）
##   · 2013 起 PFOA/PFOS 拆为线性 + 支链两段，支链段绝大多数 <LOD
##     → "合计 <LOD" 当且仅当两段均 <LOD（不是任一段）
##
## 输出：verification/lod_audit.csv
##       verification/lod_substitution_sensitivity.csv
## ============================================================================
suppressMessages({library(data.table); library(haven); library(MASS)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
RAW  <- file.path(BASE, "data/raw"); VER <- file.path(BASE, "verification")
set.seed(20261117)
MW <- c(PFOA = 414.07, PFOS = 500.13, PFHxS = 400.11, PFNA = 464.08,
        PFDA = 514.11, PFUnDA = 564.14)

rd <- function(s) {
  p <- file.path(RAW, paste0(s, ".xpt"))
  if (!file.exists(p)) return(NULL)
  as.data.table(suppressWarnings(read_xpt(p)))
}

## ---------------- A. 逐周期抽取值 + 注释码 -----------------------------------
## 直接型（2003-2012）：LBX*** + LBD***L
ext_direct <- function(file, cyc) {
  d <- rd(file); if (is.null(d)) return(NULL)
  T <- data.table(SEQN = as.numeric(d$SEQN), cycle = cyc)
  M <- c(PFOA="LBXPFOA", PFOS="LBXPFOS", PFHxS="LBXPFHS", PFNA="LBXPFNA",
         PFDA="LBXPFDE", PFUnDA="LBXPFUA")
  for (a in names(M)) {
    v <- M[[a]]; l <- sub("^LBX", "LBD", v); l <- paste0(l, "L")
    T[, (a) := if (v %in% names(d)) as.numeric(d[[v]]) else NA_real_]
    T[, (paste0(a, "_l")) := if (l %in% names(d)) as.numeric(d[[l]]) else NA_real_]
    T[, (paste0(a, "_l2")) := NA_real_]
  }
  T
}
## 拆分型（2013-2020）：PFOA = 线性 + 支链，PFOS = 线性 + 支链
ext_split <- function(file, cyc, surplus = NULL) {
  d <- rd(file); if (is.null(d)) return(NULL)
  T <- data.table(SEQN = as.numeric(d$SEQN), cycle = cyc)
  for (a in c("PFHxS","PFNA","PFDA","PFUnDA")) {
    v <- c(PFHxS="LBXPFHS", PFNA="LBXPFNA", PFDA="LBXPFDE", PFUnDA="LBXPFUA")[[a]]
    l <- paste0(sub("^LBX", "LBD", v), "L")
    T[, (a) := as.numeric(d[[v]])]
    T[, (paste0(a, "_l")) := as.numeric(d[[l]])]
    T[, (paste0(a, "_l2")) := NA_real_]
  }
  src <- d
  if (!is.null(surplus)) {
    sp <- rd(surplus)
    if (!is.null(sp)) {
      T <- merge(T, data.table(SEQN = as.numeric(sp$SEQN),
                               PFOA_n = as.numeric(sp$SSNPFOA), PFOA_b = as.numeric(sp$SSBPFOA),
                               PFOA_ln = as.numeric(sp$SDNPFOAL), PFOA_lb = as.numeric(sp$SDBPFOAL),
                               PFOS_n = as.numeric(sp$SSNPFOS), PFOS_b = as.numeric(sp$SSMPFOS),
                               PFOS_ln = as.numeric(sp$SDNPFOSL), PFOS_lb = as.numeric(sp$SDMPFOSL)),
                  by = "SEQN", all.x = TRUE)
    }
  } else {
    T <- merge(T, data.table(SEQN = as.numeric(src$SEQN),
                             PFOA_n = as.numeric(src$LBXNFOA), PFOA_b = as.numeric(src$LBXBFOA),
                             PFOA_ln = as.numeric(src$LBDNFOAL), PFOA_lb = as.numeric(src$LBDBFOAL),
                             PFOS_n = as.numeric(src$LBXNFOS), PFOS_b = as.numeric(src$LBXMFOS),
                             PFOS_ln = as.numeric(src$LBDNFOSL), PFOS_lb = as.numeric(src$LBDMFOSL)),
                by = "SEQN", all.x = TRUE)
  }
  for (a in c("PFOA","PFOS")) {
    vn <- T[[paste0(a, "_n")]]; vb <- T[[paste0(a, "_b")]]
    ln <- T[[paste0(a, "_ln")]]; lb <- T[[paste0(a, "_lb")]]
    tot <- rowSums(cbind(vn, vb), na.rm = TRUE); tot[is.na(vn) & is.na(vb)] <- NA
    ## 合计 <LOD 当且仅当两段均 <LOD；支链段单独删失率另记
    fl <- ifelse(is.na(ln) & is.na(lb), NA_real_,
          ifelse(!is.na(ln) & !is.na(lb), as.numeric(ln == 1 & lb == 1),
            ifelse(!is.na(ln), ln, lb)))
    T[, (a) := tot]
    T[, (paste0(a, "_l")) := fl]
    T[, (paste0(a, "_l2")) := lb]      # 支链段注释码（用于报告分段删失）
  }
  T[, c(grep("_(n|b|ln|lb)$", names(T), value = TRUE)) := NULL]
  T
}
CYC <- list(C = list(f = "L24PFC_C", k = "d"),
            D = list(f = "PFC_D",   k = "d"),
            F = list(f = "PFC_F",   k = "d"),
            G = list(f = "PFC_G",   k = "d"),
            H = list(f = "PFAS_H",  k = "s", x = "SSPFAS_H"),
            I = list(f = "PFAS_I",  k = "s"),
            J = list(f = "PFAS_J",  k = "s"),
            L = list(f = "PFAS_L",  k = "s"))
B6 <- rbindlist(lapply(names(CYC), function(s) {
  z <- CYC[[s]]
  if (z$k == "d") ext_direct(z$f, s) else ext_split(z$f, s, z$x)
}), fill = TRUE)
B6[, n_ok := rowSums(!is.na(.SD)), .SDcols = names(MW)]

INF <- readRDS(file.path(BASE, "data/derived/inflam.rds"))
SAMPLESEQ <- INF[ok_infl == TRUE & !is.na(w5) & !is.na(pfas_z), SEQN]

## 质量点识别：某取值若"只出现在 <LOD 记录里"，即为替换值
mass_info <- function(v, fl) {
  ok <- !is.na(v) & !is.na(fl)
  if (!any(ok)) return(list(n = 0, pct = NA_real_, k = 0L, lo = NA_real_, hi = NA_real_,
                            vals = character(0)))
  u <- sort(unique(round(v[ok], 8)))
  issub <- vapply(u, function(uu) all(fl[ok][abs(v[ok] - uu) < 1e-9] == 1), logical(1))
  if (!any(issub)) return(list(n = 0, pct = 0, k = 0L, lo = NA_real_, hi = NA_real_,
                               vals = character(0)))
  sub <- u[issub]
  n <- sum(v[ok] %in% sub)
  lst <- sub("0+\\.0+$", "0", formatC(sub, format = "fg", digits = 6))
  list(n = n, pct = 100 * n / sum(ok), k = length(sub),
       lo = min(sub) * sqrt(2), hi = max(sub) * sqrt(2), vals = lst)
}

AUD <- list()
for (s in names(CYC)) for (a in names(MW)) {
  T <- B6[cycle == s]
  v <- T[[a]]; fl <- T[[paste0(a, "_l")]]
  if (is.null(v) || all(is.na(v))) next
  ok <- !is.na(v) & v > 0
  nb <- sum(fl[ok] == 1, na.rm = TRUE)
  mi <- mass_info(v, fl)
  inS <- T$SEQN %in% SAMPLESEQ
  AUD[[length(AUD)+1]] <- data.table(
    cycle = s, analyte = a, n_measured = sum(ok),
    n_below = nb, pct_below = round(100 * nb / sum(ok), 2),
    n_imputed = mi$n, pct_imputed = round(mi$pct, 2),
    n_lod_levels = mi$k,
    lod_min = signif(mi$lo, 4), lod_max = signif(mi$hi, 4),
    subst_values = paste(mi$vals, collapse = "/"),
    n_in_analytic = sum(inS & ok),
    pct_below_analytic = round(100 * sum(fl[inS & ok] == 1, na.rm = TRUE) / sum(inS & ok), 2),
    branch_censored = if (a %in% c("PFOA","PFOS") && s %in% c("H","I","J","L"))
      round(100 * sum(T[[paste0(a, "_l2")]] == 1, na.rm = TRUE) /
            sum(!is.na(T[[paste0(a, "_l2")]])), 2) else NA_real_)
}
AUDT <- rbindlist(AUD); setorder(AUDT, analyte, cycle)
AUDT[, flag_vs_mass_diff_pp := round(abs(pct_below - pct_imputed), 2)]
fwrite(AUDT, file.path(VER, "lod_audit.csv"))

cat("========== A. % < LOD（注释码口径）==========\n")
print(dcast(AUDT, analyte ~ cycle, value.var = "pct_below"), row.names = FALSE)
cat("\n========== A2. % < LOD（分析样本内口径，N=12,775）==========\n")
print(dcast(AUDT, analyte ~ cycle, value.var = "pct_below_analytic"), row.names = FALSE)
cat("\n========== A3. 反推 LOD 区间（ng/mL）==========\n")
print(AUDT[, .(cycle, analyte, n_lod_levels, lod_min, lod_max, subst_values)][
        n_lod_levels > 0][order(analyte, cycle)], row.names = FALSE)

cat("\n========== B. 双口径交叉验证（注释码 vs 质量点）==========\n")
cat("最大偏差 =", max(AUDT$flag_vs_mass_diff_pp, na.rm = TRUE), "pp\n")
cat("偏差 > 1 pp 的组合：", AUDT[flag_vs_mass_diff_pp > 1, .N], "个\n")
if (AUDT[flag_vs_mass_diff_pp > 1, .N]) print(
  AUDT[flag_vs_mass_diff_pp > 1, .(cycle, analyte, pct_below, pct_imputed,
                                   flag_vs_mass_diff_pp)], row.names = FALSE)

cat("\n========== 删失率 > 30% 的同系物×周期 ==========\n")
hi <- AUDT[pct_below > 30]
print(hi[, .(cycle, analyte, pct_below, lod_min = lod_min)], row.names = FALSE)
cat("\n涉及周期：", paste(sort(unique(hi$cycle)), collapse = ", "), "\n")
cat("涉及同系物：", paste(sort(unique(hi$analyte)), collapse = ", "), "\n")
cat("\n支链段单独删失率（2013 起 PFOA/PFOS 拆分为线性+支链）：\n")
print(AUDT[!is.na(branch_censored), .(cycle, analyte, branch_censored)][order(analyte, cycle)],
      row.names = FALSE)

## ---------------- C. 替换值敏感性 -------------------------------------------
## <LOD 记录分别用 LOD/√2（现行）、LOD（上限）、LOD/2（下限）替换，
## 重建 ΣPFAS 摩尔和 → log → 周期内 z（与 06_build_pop.R 完全同法）→ 重跑交互
mkz <- function(DT, scheme = c("lod_sqrt2", "lod", "lod_half")) {
  scheme <- match.arg(scheme)
  X <- copy(DT)
  mult <- switch(scheme, lod_sqrt2 = 1, lod = sqrt(2), lod_half = 1/sqrt(2))
  if (mult != 1) for (a in names(MW)) {
    l <- X[[paste0(a, "_l")]]; v <- X[[a]]
    idx <- !is.na(l) & l == 1 & !is.na(v)
    X[idx, (a) := v[idx] * mult]
  }
  X[, molar := {
    m <- matrix(NA_real_, .N, length(MW))
    for (j in seq_along(MW)) m[, j] <- get(names(MW)[j]) / MW[j] * 1000
    rowSums(m, na.rm = TRUE)
  }]
  X[n_ok < 6, molar := NA_real_]
  X[, logp := log(molar)]
  X[, z := as.numeric(scale(logp)), by = cycle]
  X[, .(SEQN, z)]
}

A <- readRDS(file.path(BASE, "data/derived/analysis.rds")); CONF <- A$CONF
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
BASE0 <- INF[ok_infl == TRUE & !is.na(w5) & !is.na(y_wbc),
             c("SEQN", "cycle", "y_wbc", "w5", "strata", "psu", CONF), with = FALSE]

runInt <- function(DZ, tag, drop_cycle = NULL) {
  dd <- merge(BASE0, DZ, by = "SEQN")[!is.na(z)]
  if (!nrow(dd)) return(NULL)
  if (!is.null(drop_cycle)) dd <- dd[!cycle %in% drop_cycle]
  dd[, mv := z]
  cf <- CONF[sapply(CONF, function(v) length(unique(dd[[v]])) > 1)]
  X <- mkX(as.formula(paste("y_wbc ~ w5*mv +", paste(cf, collapse = " + "))), dd)
  m <- lm.fit(X, dd$y_wbc); cl <- factor(paste(dd$strata, dd$psu))
  V <- crsV(X, m$residuals, cl)
  j <- grep("^w5:mv$|^mv:w5$", colnames(X))
  b <- unname(m$coefficients[j]); se <- sqrt(V[j, j])
  data.table(scenario = tag, N = nrow(dd), n_cycle = uniqueN(dd$cycle),
             beta = b, se = se, p = 2*pnorm(-abs(b/se)))
}

SENS <- list()
for (sc in c("lod_sqrt2", "lod", "lod_half"))
  SENS[[length(SENS)+1]] <- runInt(mkz(B6, sc),
    paste0("subst = ", switch(sc, lod_sqrt2 = "LOD/sqrt2 (current)",
                              lod = "LOD (upper)", lod_half = "LOD/2 (lower)")))
## 注：删失 >30% 只发生在 PFDA / PFUnDA 上（且 PFUnDA 在全部 8 周期都 >30%），
## 因此"剔除高删失周期"不可行 —— 正确做法是"剔除高删失同系物"（见下）与"逐同系物单独建模"
## 另加：只保留 PFOA/PFOS/PFHxS/PFNA 四种低删失同系物重建 ΣPFAS
B4 <- copy(B6)
B4[, n_ok := rowSums(!is.na(.SD)), .SDcols = c("PFOA","PFOS","PFHxS","PFNA")]
mkz4 <- function(DT) {
  X <- copy(DT)
  X[, molar := {
    m <- matrix(NA_real_, .N, 4)
    for (j in 1:4) m[, j] <- get(c("PFOA","PFOS","PFHxS","PFNA")[j]) / MW[j] * 1000
    rowSums(m, na.rm = TRUE)
  }]
  X[n_ok < 4, molar := NA_real_]
  X[, z := as.numeric(scale(log(molar))), by = cycle]
  X[, .(SEQN, z)]
}
SENS[[length(SENS)+1]] <- runInt(mkz4(B4), "only 4 low-censoring congeners (drop PFDA/PFUnDA)")

## 逐同系物单独建模：哪个同系物在驱动交互？（同时回答"高删失同系物是否制造结论"）
for (a in names(MW)) {
  Xa <- copy(B6)
  Xa <- Xa[!is.na(get(a)) & get(a) > 0]
  Xa[, z := as.numeric(scale(log(get(a)))), by = cycle]
  Xa <- Xa[!is.na(z), .(SEQN, z)]
  pb <- AUDT[analyte == a, max(pct_below)]
  SENS[[length(SENS)+1]] <- runInt(Xa, sprintf("single congener: %s (max %%=<LOD %.0f%%)", a, pb))
}

SENS <- rbindlist(SENS[!sapply(SENS, is.null)])
SENS[, `:=`(rel_beta = sprintf("%+.1f%%", 100*(beta/beta[1] - 1)),
            pct_per5_perSD = sprintf("%+.2f%%", 100*(exp(beta)-1)))]
fwrite(SENS, file.path(VER, "lod_substitution_sensitivity.csv"))

cat("\n========== C/D. 替换值与周期敏感性（WBC 连续 DiD 交互）==========\n")
print(SENS[, .(scenario, N, n_cycle, beta = signif(beta,4), se = signif(se,4),
               p = signif(p,3), rel = rel_beta)], row.names = FALSE)
cat("\n点估计最大相对偏差 =", sprintf("%.1f%%", max(abs(100*(SENS$beta/SENS$beta[1]-1)))),
    "；P 值范围 =", paste(signif(range(SENS$p), 3), collapse = " ~ "), "\n")
cat("✅ V17 完成 → verification/lod_audit.csv, lod_substitution_sensitivity.csv\n")
