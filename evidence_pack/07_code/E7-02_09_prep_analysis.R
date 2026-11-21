## ============================================================================
## 09_prep_analysis.R — 构建分析样本 + 混杂变量矩阵 + 缺失填补
## ============================================================================
suppressMessages({library(data.table); library(splines)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
CORE <- readRDS(file.path(BASE, "data/derived/core.rds"))
POP  <- readRDS(file.path(BASE, "data/derived/pop.rds"))
HORM <- readRDS(file.path(BASE, "data/derived/horm.rds"))

D <- merge(CORE, POP[, .(SEQN, PFOA, PFOS, PFHxS, PFNA, PFDA, PFUnDA,
                         pfas_molar, log_pfas, pfas_z, pfos_z, pfoa_z, pfas_wet,
                         pcb_lipid, log_pcb, pcb_z, ocp_lipid, log_ocp, ocp_z,
                         wt_pf, wt_pop)], by = "SEQN", all.x = TRUE)
D <- merge(D, HORM[, .(SEQN, tt_ngdl, shbg, e2_pgml, ft_nmol, ft_ngdl, pct_free)],
           by = "SEQN", all.x = TRUE)
D <- D[adult == TRUE]

## ---------- 生物学年龄加速 ----------
D[, phenoage_accel := {
  r <- rep(NA_real_, .N)
  ok <- !is.na(phenoage) & !is.na(age)
  if (sum(ok) > 100) {
    fit <- lm(phenoage ~ ns(age, df = 4), data = .SD[ok])
    r[ok] <- resid(fit)
  }
  r
}]
D[, age_accel := phenoage - age]

## ---------- 处理变量 ----------
D[, `:=`(
  w5   = as.integer(dw1 >= 5),
  w10  = as.integer(dw1 >= 10),
  w5s  = as.integer(dw1 >= 5 & !is.na(dw1)),       # 严格版(非缺失)
  stable = as.integer(abs(dw1) < 2)
)]
D[dw1 >= 5, w5 := 1L]; D[dw1 < 5, w5 := 0L]

## ---------- 混杂变量 ----------
D[, `:=`(
  log_cot = log(pmax(cotinine, 0.01)),
  yr      = srvyr,
  bmi_p   = bmi_prior,
  drink2  = pmin(drinks, 10)
)]
## 共病负担
D[, comorbid := as.integer(diab == 1) + as.integer(htn == 1)]
D[is.na(comorbid), comorbid := NA_integer_]

CONF <- c("age", "female", "pir", "bmi_p", "log_cot", "drink2", "pa_any",
          "comorbid", "tc", "hdl", "egfr", "yr")
CATS <- c("race_f", "edu_f", "smk")

## 中位数/众数填补 + 缺失指示
imp <- function(dt, vars) {
  for (v in vars) {
    x <- dt[[v]]
    if (is.na(mean(is.na(x))) ) next
    dt[, (paste0("miss_", v)) := as.integer(is.na(get(v)))]
    if (is.numeric(x)) {
      m <- median(x, na.rm = TRUE); if (!is.finite(m)) m <- 0
      dt[is.na(get(v)), (v) := m]
    } else {
      tb <- sort(table(x), decreasing = TRUE)
      lv <- if (length(tb)) names(tb)[1] else levels(x)[1]
      dt[is.na(get(v)), (v) := lv]
    }
  }
  dt
}
D <- imp(D, c(CONF, CATS))

## 设计矩阵（数值）
mkX <- function(dt, extra = character()) {
  X <- cbind(as.matrix(dt[, ..CONF]),
             model.matrix(~ race_f + edu_f + smk - 1, data = dt))
  if (length(extra)) X <- cbind(X, as.matrix(dt[, ..extra]))
  miss_cols <- grep("^miss_", names(dt), value = TRUE)
  if (length(miss_cols)) X <- cbind(X, as.matrix(dt[, ..miss_cols]))
  storage.mode(X) <- "double"
  X
}

## ---------- 分析样本 ----------
## S1: PFAS × 生物学年龄 (8 周期)
S1 <- D[!is.na(phenoage_accel) & !is.na(w5) & !is.na(pfas_z)]
## S2: PCB/OCP × 生物学年龄 (2003-2004)
S2 <- D[!is.na(phenoage_accel) & !is.na(w5) & !is.na(pcb_z)]
## S2b: OCP 子样本
S2b <- D[!is.na(phenoage_accel) & !is.na(w5) & !is.na(ocp_z)]
## S3: 男性总睾酮 × PFAS
S3 <- D[female == 0 & !is.na(tt_ngdl) & !is.na(w5) & !is.na(pfas_z) & tt_ngdl > 0]
## S4: 男性游离睾酮 × PFAS
S4 <- D[female == 0 & !is.na(ft_nmol) & !is.na(w5) & !is.na(pfas_z) & ft_nmol > 0]
## S0: 全样本 体重变化 × 生物学年龄（大样本背景分析）
S0 <- D[!is.na(phenoage_accel) & !is.na(w5)]

cat("=== 分析样本量 ===\n")
cat("S0 全样本(体重变化×PhenoAge)  N =", nrow(S0), " 处理组 =", sum(S0$w5), "\n")
cat("S1 PFAS×PhenoAge               N =", nrow(S1), " 处理组 =", sum(S1$w5), "\n")
cat("S2 PCB×PhenoAge (2003-04)      N =", nrow(S2), " 处理组 =", sum(S2$w5), "\n")
cat("S2b OCP×PhenoAge (2003-04)     N =", nrow(S2b), " 处理组 =", sum(S2b$w5), "\n")
cat("S3 男性总睾酮×PFAS             N =", nrow(S3), " 处理组 =", sum(S3$w5), "\n")
cat("S4 男性游离睾酮×PFAS           N =", nrow(S4), " 处理组 =", sum(S4$w5), "\n")
cat("\n各周期 S1 分布:\n"); print(S1[, .N, by = cycle][order(cycle)])

saveRDS(list(D = D, S0 = S0, S1 = S1, S2 = S2, S2b = S2b, S3 = S3, S4 = S4,
             CONF = CONF, CATS = CATS), file.path(BASE, "data/derived/analysis.rds"))

## 描述统计
desc <- S1[, .(N = .N,
  age = round(mean(age), 1), female = round(mean(female) * 100, 1),
  bmi = round(mean(bmi, na.rm = TRUE), 1), bmi_prior = round(mean(bmi_p), 1),
  pfas_med = round(median(pfas_molar, na.rm = TRUE), 1),
  accel = round(mean(phenoage_accel), 2)), by = w5]
print(desc)
fwrite(S1[, .(SEQN, cycle, age, female, w5, dw1, phenoage, phenoage_accel,
              pfas_molar, pfas_z, pfas_z_q = cut(pfas_z, quantile(pfas_z, 0:10/10), include.lowest = TRUE))],
       file.path(BASE, "results/S1_subject_level.csv"))
cat("done\n")
