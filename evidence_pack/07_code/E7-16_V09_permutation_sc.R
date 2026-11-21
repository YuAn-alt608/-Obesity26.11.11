## =====================================================================
## V09  第三层：单细胞/汇池侧精确置换检验（Exact Permutation）
## ---------------------------------------------------------------------
## 每组仅 6 个汇池，Welch t 检验的正态近似不可靠。这里做精确枚举：
##   (a) WL vs Obese：12 池分成 6+6，全部 C(12,6) = 924 种分组
##   (b) 18 池中任取 6 池标为 "WL"，与其余 12 池比较，C(18,6) = 18,564 种
## 统计量：巨噬细胞占比、每细胞炎症评分 z、每细胞促炎输出（假体 bulk）
## 产出：verification/permutation_sc.csv、permutation_sc_exact12.csv
## =====================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", ".")
OUT  <- file.path(BASE, "verification")

comp <- fread(file.path(BASE, "results/geo_composition.csv"))
scores <- fread(file.path(BASE, "results/geo_cell_scores.csv.gz"))
m12 <- fread(file.path(BASE, "results/sc_m1m2.csv"))

## 每池一个观测：巨噬细胞占比
mac <- comp[celltype == "macrophage", .(gsm, group, pct)]
## 每池一个观测：巨噬细胞炎症 z（逐细胞均值）
zi <- scores[celltype == "macrophage", .(z = mean(inflamma, na.rm = TRUE)), by = .(gsm, group)]

D <- merge(mac, zi, by = c("gsm", "group"))
D <- merge(D, m12[, .(gsm, inflam_output)], by = "gsm", all.x = TRUE)
setorder(D, group, gsm)
cat("汇池数:", nrow(D), " 分组:", paste(D[, .N, by = group][[1]], D[, .N, by = group][[2]], collapse = "/"), "\n")

## ---------- (a) 12 池精确枚举 6+6 --------------------------------------------
exact12 <- function(x, grp) {
  idx <- which(grp != "Lean")
  v <- x[idx]
  n <- length(v); k <- sum(grp[idx] == "WL")
  cmb <- combn(n, k)
  obs <- mean(v[grp[idx] == "WL"]) - mean(v[grp[idx] == "Obese"])
  dif <- apply(cmb, 2, function(ii) mean(v[ii]) - mean(v[-ii]))
  c(diff = obs, p = mean(abs(dif) >= abs(obs) - 1e-12), n_perm = ncol(cmb))
}
R12 <- rbindlist(lapply(c("pct", "z", "inflam_output"), function(k) {
  r <- exact12(D[[k]], D$group)
  data.table(metric = c(pct = "Macrophage share (%)", z = "Inflammation z per macrophage",
                        inflam_output = "Pro-inflammatory output")[[k]],
             diff = r["diff"], p_exact = r["p"], n_perm = r["n_perm"])
}))

## ---------- (b) 18 池 C(18,6) 全枚举 -----------------------------------------
exact18 <- function(x, grp) {
  n <- length(x); k <- sum(grp == "WL")
  cmb <- combn(n, k)
  obs <- mean(x[grp == "WL"]) - mean(x[grp != "WL"])
  dif <- apply(cmb, 2, function(ii) mean(x[ii]) - mean(x[-ii]))
  c(diff = obs, p = mean(abs(dif) >= abs(obs) - 1e-12), n_perm = ncol(cmb))
}
R18 <- rbindlist(lapply(c("pct", "z"), function(k) {
  r <- exact18(D[[k]], D$group)
  data.table(metric = c(pct = "Macrophage share (%)", z = "Inflammation z per macrophage")[[k]],
             diff = r["diff"], p_exact = r["p"], n_perm = r["n_perm"])
}))

fwrite(R12, file.path(OUT, "permutation_sc_exact12.csv"))
fwrite(R18, file.path(OUT, "permutation_sc.csv"))

cat("\n========== (a) WL vs Obese 精确枚举 C(12,6)=924 ==========\n")
print(R12[, .(metric, diff = round(diff, 4), p_exact = signif(p_exact, 3), n_perm)], row.names = FALSE)
cat("\n========== (b) 18 池任取 6 池 C(18,6)=18564 ==========\n")
print(R18[, .(metric, diff = round(diff, 4), p_exact = signif(p_exact, 3), n_perm)], row.names = FALSE)
cat("\n（Welch t 检验对照：巨噬占比 P=0.010，炎症 z P=0.0003；见报告 §4）\n")
