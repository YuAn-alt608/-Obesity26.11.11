## ============================================================================
## 13_mortality.R — 全因死亡随访
##   公共链接死亡文件仅提供死亡状态/死因(无具体日期)，故:
##   · 汇总分析: 泊松率模型 + log(随访人年) 偏移 (等价分段常数风险)
##   · 2003-2004 单周期: 随访长度近似相同, 用 logistic (近似 Cox)
## ============================================================================
suppressMessages({library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
MDIR <- file.path(BASE, "data/raw/mort")
cyc <- data.table(cycle = c("C","D","F","G","H","I","J"),
                  f = c("2003_2004","2005_2006","2009_2010","2011_2012","2013_2014","2015_2016","2017_2018"),
                  yr = c(2003,2005,2009,2011,2013,2015,2017))

M <- rbindlist(lapply(seq_len(nrow(cyc)), function(i) {
  L <- gsub("\r", "", readLines(file.path(MDIR, paste0("NHANES_", cyc$f[i],
             "_MORT_2019_PUBLIC.dat")), warn = FALSE))
  L <- L[nchar(L) >= 26]
  data.table(SEQN = suppressWarnings(as.numeric(trimws(substr(L, 1, 6)))),
             mortstat = suppressWarnings(as.integer(substr(L, 16, 16))),
             ucod = trimws(substr(L, 17, 19)),
             cycle = cyc$cycle[i], yr = cyc$yr[i])
}))
M <- M[!is.na(SEQN)]
END <- as.Date("2019-12-31")
M[, entry := as.Date(paste0(yr + 1, "-01-01"))]
M[, fu_years := as.numeric(END - entry) / 365.25]
M[, event := mortstat]
cat("死亡记录:", nrow(M), " 死亡:", sum(M$event, na.rm = TRUE), "\n")
print(M[, .(N = .N, deaths = sum(event, na.rm = TRUE), fu = round(mean(fu_years), 1)), by = cycle][order(cycle)])

A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
D <- merge(A$D, M[, .(SEQN, fu_years, event, ucod)], by = "SEQN", all.x = TRUE)
D <- D[!is.na(fu_years) & !is.na(event)]
cat("\n可分析成人:", nrow(D), " 死亡:", sum(D$event), "\n")

fmt <- function(b, se) sprintf("%.3f (95%%CI %.3f-%.3f) p=%s", exp(b),
                               exp(b - 1.96*se), exp(b + 1.96*se),
                               format.pval(2 * pnorm(-abs(b / se)), digits = 3))

## ---------- 1) PhenoAge 验证 ----------
cat("\n=== 1) PhenoAge 加速 -> 全因死亡 ===\n")
v <- D[!is.na(phenoage_accel)]
m1 <- glm(event ~ phenoage_accel + age + female + race_f + edu_f + bmi + smk +
            factor(cycle), family = poisson, offset = log(fu_years), data = v)
s1 <- summary(m1)$coefficients
cat("每 1 年加速的率比 RR =", fmt(s1["phenoage_accel", 1], s1["phenoage_accel", 2]), "\n")
v[, dec := cut(phenoage_accel, quantile(phenoage_accel, 0:5/5), include.lowest = TRUE, labels = 1:5)]
m1b <- glm(event ~ dec + age + female + race_f + edu_f + bmi + smk + factor(cycle),
           family = poisson, offset = log(fu_years), data = v)
cat("按加速五分位的率比(以Q1为参照):\n")
print(round(exp(c(0, coef(m1b)[grep("^dec", names(coef(m1b)))] )), 3))
## 区分度: 按周期内的 AUC
auc <- function(y, x) { r <- rank(x); n1 <- sum(y == 1); n0 <- sum(y == 0)
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
cat("AUC(实足年龄) =", round(auc(v$event, v$age), 3),
    " AUC(PhenoAge) =", round(auc(v$event, v$phenoage), 3), "\n")

## ---------- 2) 减重 × PFAS -> 死亡 ----------
cat("\n=== 2) 减重 × ΣPFAS -> 全因死亡 ===\n")
m <- D[!is.na(pfas_z) & !is.na(w5)]
cat("N =", nrow(m), " 死亡 =", sum(m$event), " 处理组死亡 =", sum(m$event[m$w5 == 1]), "\n")
m2 <- glm(event ~ w5 * pfas_z + age + female + race_f + edu_f + bmi_p + smk +
            comorbid + factor(cycle), family = poisson, offset = log(fu_years), data = m)
s2 <- summary(m2)$coefficients
b <- coef(m2); V <- vcov(m2); q <- quantile(m$pfas_z, c(.1, .5, .9))
rr <- exp(b["w5"] + b["w5:pfas_z"] * q)
se <- sqrt(V["w5","w5"] + q^2 * V["w5:pfas_z","w5:pfas_z"] + 2*q*V["w5","w5:pfas_z"])
tab2 <- data.table(pfas_z = round(q, 2), RR = round(rr, 3),
                   lo = round(exp(log(rr) - 1.96*se), 3), hi = round(exp(log(rr) + 1.96*se), 3))
print(tab2)
cat("交互项 p =", format.pval(s2["w5:pfas_z", 4], digits = 3), "\n")

## ---------- 3) 2003-2004 单周期 logistic ----------
cat("\n=== 3) 2003-2004: 减重 × PCB -> 16 年死亡 (logistic) ===\n")
c0 <- D[cycle == "C" & !is.na(pcb_z) & !is.na(w5)]
cat("N =", nrow(c0), " 死亡 =", sum(c0$event), "\n")
m3 <- glm(event ~ w5 * pcb_z + age + female + race_f + edu_f + bmi_p + smk + comorbid,
          family = binomial, data = c0)
s3 <- summary(m3)$coefficients
b3 <- coef(m3); V3 <- vcov(m3); q3 <- quantile(c0$pcb_z, c(.1, .5, .9))
or <- exp(b3["w5"] + b3["w5:pcb_z"] * q3)
se3 <- sqrt(V3["w5","w5"] + q3^2 * V3["w5:pcb_z","w5:pcb_z"] + 2*q3*V3["w5","w5:pcb_z"])
print(data.table(pcb_z = round(q3, 2), OR = round(or, 3),
                 lo = round(exp(log(or) - 1.96*se3), 3), hi = round(exp(log(or) + 1.96*se3), 3)))
cat("交互项 p =", format.pval(s3["w5:pcb_z", 4], digits = 3), "\n")

## ---------- 4) 死亡原因的探索 (2003-2004) ----------
cat("\n=== 4) 减重组死因构成 (2003-2004, 前三位) ===\n")
print(head(c0[w5 == 1 & event == 1, .N, by = ucod][order(-N)], 6))

fwrite(data.table(var = rownames(s1), RR = exp(s1[,1]), lo = exp(s1[,1]-1.96*s1[,2]),
                  hi = exp(s1[,1]+1.96*s1[,2]), p = s1[,4]),
       file.path(BASE, "results/mort_phenoage.csv"))
fwrite(tab2, file.path(BASE, "results/mort_w5_pfas.csv"))
saveRDS(list(m1 = m1, m2 = m2, m3 = m3, tab2 = tab2,
             auc = c(age = auc(v$event, v$age), pheno = auc(v$event, v$phenoage))),
        file.path(BASE, "results/mort_models.rds"))
cat("done\n")
