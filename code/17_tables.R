## ============================================================================
## 17_tables.R — 基线特征表 + 主结果汇总
## ============================================================================
suppressMessages({library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", "."); RES <- file.path(BASE, "results")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
S1 <- A$S1

num <- function(x) sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
med <- function(x) { q <- quantile(x, c(.25, .5, .75), na.rm = TRUE)
  sprintf("%.2f (%.2f-%.2f)", q[2], q[1], q[3]) }
pct <- function(x, lv = 1) sprintf("%.1f", mean(x == lv, na.rm = TRUE) * 100)

t1 <- data.table(
  Variable = c("N", "Age, y", "Female, %", "NH White, %", "NH Black, %", "Hispanic, %",
               "College or above, %", "Family PIR", "BMI, kg/m2", "BMI 1 y ago, kg/m2",
               "Waist, cm", "Current smoker, %", "Serum cotinine, ng/mL",
               "Physically active, %", "Diabetes, %", "Hypertension, %",
               "Total cholesterol, mg/dL", "HDL, mg/dL", "eGFR, mL/min/1.73m2",
               "1-y weight change, %", "Lost >=10% in 1 y, %",
               "PhenoAge, y", "PhenoAge acceleration, y",
               "Sum PFAS, nmol/L", "PFOS, ng/mL", "PFOA, ng/mL", "PFHxS, ng/mL",
               "PFNA, ng/mL", "PFDA, ng/mL", "PFUnDA, ng/mL"),
  `No rapid loss` = c(
    format(nrow(S1[w5 == 0]), big.mark = ","), num(S1[w5 == 0]$age), pct(S1[w5 == 0]$female),
    pct(S1[w5 == 0]$race, 3), pct(S1[w5 == 0]$race, 4),
    pct(S1[w5 == 0]$race, 1) , pct(S1[w5 == 0]$edu, 5),
    num(S1[w5 == 0]$pir), num(S1[w5 == 0]$bmi), num(S1[w5 == 0]$bmi_p),
    num(S1[w5 == 0]$waist), pct(S1[w5 == 0]$smk, 2), med(S1[w5 == 0]$cotinine),
    pct(S1[w5 == 0]$pa_any, 1), pct(S1[w5 == 0]$diab, 1), pct(S1[w5 == 0]$htn, 1),
    num(S1[w5 == 0]$tc), num(S1[w5 == 0]$hdl), num(S1[w5 == 0]$egfr),
    num(S1[w5 == 0]$dw1), pct(as.integer(S1[w5 == 0]$dw1 >= 10), 1),
    num(S1[w5 == 0]$phenoage), num(S1[w5 == 0]$phenoage_accel),
    med(S1[w5 == 0]$pfas_molar), med(S1[w5 == 0]$PFOS), med(S1[w5 == 0]$PFOA),
    med(S1[w5 == 0]$PFHxS), med(S1[w5 == 0]$PFNA), med(S1[w5 == 0]$PFDA),
    med(S1[w5 == 0]$PFUnDA)),
  `Rapid loss >=5%` = c(
    format(nrow(S1[w5 == 1]), big.mark = ","), num(S1[w5 == 1]$age), pct(S1[w5 == 1]$female),
    pct(S1[w5 == 1]$race, 3), pct(S1[w5 == 1]$race, 4),
    pct(S1[w5 == 1]$race, 1), pct(S1[w5 == 1]$edu, 5),
    num(S1[w5 == 1]$pir), num(S1[w5 == 1]$bmi), num(S1[w5 == 1]$bmi_p),
    num(S1[w5 == 1]$waist), pct(S1[w5 == 1]$smk, 2), med(S1[w5 == 1]$cotinine),
    pct(S1[w5 == 1]$pa_any, 1), pct(S1[w5 == 1]$diab, 1), pct(S1[w5 == 1]$htn, 1),
    num(S1[w5 == 1]$tc), num(S1[w5 == 1]$hdl), num(S1[w5 == 1]$egfr),
    num(S1[w5 == 1]$dw1), pct(as.integer(S1[w5 == 1]$dw1 >= 10), 1),
    num(S1[w5 == 1]$phenoage), num(S1[w5 == 1]$phenoage_accel),
    med(S1[w5 == 1]$pfas_molar), med(S1[w5 == 1]$PFOS), med(S1[w5 == 1]$PFOA),
    med(S1[w5 == 1]$PFHxS), med(S1[w5 == 1]$PFNA), med(S1[w5 == 1]$PFDA),
    med(S1[w5 == 1]$PFUnDA))
)
fwrite(t1, file.path(RES, "table1_baseline.csv"))

## 主结果汇总 (参数模型为主, 森林异质性检验为辅)
cf <- readRDS(file.path(RES, "cf_results.rds"))
par <- fread(file.path(RES, "param_interaction.csv"))
lab <- c("S1  PhenoAge \u00d7 \u03a3PFAS", "S2  PhenoAge \u00d7 \u03a3PCB",
         "S2b PhenoAge \u00d7 \u03a3OCP", "S3  log\u603b\u777e\u916e \u00d7 \u03a3PFAS",
         "S4  log\u6e38\u79bb\u777e\u916e \u00d7 \u03a3PFAS")
par[, key := sub(" .*", "", model)]
cal_p <- sapply(names(cf), function(k) cf[[k]]$cal[2, 4])
names(cal_p) <- names(cf)
tab2 <- data.table(
  Model = c("PhenoAge x PFAS (8 cycles)", "PhenoAge x PCB (2003-04)",
            "PhenoAge x OCP (2003-04)", "Total testosterone x PFAS (men)",
            "Free testosterone x PFAS (men)"),
  N = par$N, Treated = c(1877, 221, 271, 439, 270),
  Effect_at_mean = fifelse(grepl("^S[34]", par$key),
                    sprintf("%+.1f%%", (exp(par$eff0) - 1) * 100),
                    sprintf("%.2f y", par$eff0)),
  Interaction_b2 = sprintf("%.2f (%.2f)", par$b2, par$se_b2),
  P_interaction = signif(par$p_joint, 3),
  CATE_P10 = fifelse(grepl("^S[34]", par$key),
                     sprintf("%+.1f%%", (exp(par$e10) - 1) * 100), sprintf("%.2f y", par$e10)),
  CATE_P90 = fifelse(grepl("^S[34]", par$key),
                     sprintf("%+.1f%%", (exp(par$e90) - 1) * 100), sprintf("%.2f y", par$e90)),
  Turning_point = round(par$turn, 2),
  Forest_heterogeneity_P = signif(cal_p[par$key], 3)
)
fwrite(tab2, file.path(RES, "table2_main.csv"))
cat("\nTable 2:\n"); print(tab2)
