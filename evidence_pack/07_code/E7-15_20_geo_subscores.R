## ============================================================================
## 20_geo_subscores.R — 脂肪细胞通路子评分 + 脂解-外排耦联检验
##   动机: 14_geo.R 的 lipolysis 模块混入了脂滴结构蛋白 (PLIN1/2/4, CIDEC),
##         这些在减重后随脂滴缩小而下调, 会把真正的脂解酶信号抵消掉。
##   这里按功能拆成纯子评分, 并直接检验"脂解激活 ⇄ 核受体外排"是否共变。
## ============================================================================
suppressMessages({library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
OUT  <- file.path(BASE, "results")

PB <- fread(file.path(OUT, "geo_pseudobulk_counts.csv.gz"))
pa <- PB[celltype == "adipocyte"]
pa[, cpm := expr / sum(expr) * 1e6, by = gsm]

SET <- list(
  lipase_core = c("PNPLA2","LIPE","MGLL","ABHD5"),          # 真正的脂解酶
  droplet     = c("PLIN1","PLIN2","PLIN4","CIDEC"),          # 脂滴结构/包被
  ahr_axis    = c("AHR","ARNT","CYP1A1","CYP1B1","AHRR"),    # 芳香烃受体通路
  car_pxr     = c("NR1I2","NR1I3","CYP2B6","CYP3A4","ABCB1","ABCG2"),  # 核受体+外排转运
  nrf2        = c("NQO1","HMOX1","GCLC","GSTA1"),            # 抗氧化反应
  inflam      = c("NLRP3","IL1B","CASP1","CCL2","TNF","IL6","PYCARD","TLR4"),
  sasp        = c("CDKN1A","CDKN2A","SERPINE1","GLB1","TP53"),
  thermogenic = c("UCP1","PPARGC1A","DIO2","PRDM16")
)
SET <- lapply(SET, function(g) intersect(g, unique(pa$gene)))
cat("各子评分可用基因数:\n"); print(sapply(SET, length))

## 几何均值 CPM (log 后取均值再还原) → 对低表达基因更稳健
sc <- rbindlist(lapply(names(SET), function(nm) {
  pa[gene %in% SET[[nm]], .(score = mean(log10(cpm + 1))), by = .(gsm, group)][, set := nm]
}))
w <- dcast(sc, gsm + group ~ set, value.var = "score")

cat("\n=== 脂肪细胞通路子评分 (log10 CPM 均值) ===\n")
print(w[, lapply(.SD, mean), by = group, .SDcols = names(SET)], digits = 3)

cat("\n=== WL vs Obese (pool 为单位, Welch t) ===\n")
res <- rbindlist(lapply(names(SET), function(v) {
  a <- w[[v]][w$group == "Obese"]; b <- w[[v]][w$group == "WL"]; c3 <- w[[v]][w$group == "Lean"]
  data.table(set = v, Obese = mean(a), WL = mean(b), Lean = mean(c3),
             delta = mean(b) - mean(a),
             p_WL_Obese = tryCatch(t.test(b, a)$p.value, error = function(e) NA_real_),
             p_WL_Lean  = tryCatch(t.test(b, c3)$p.value, error = function(e) NA_real_))
}))
print(res, digits = 3)

## 逐 pool 展示 WL 内部的异质性 (是否存在"高 AHR 亚群")
cat("\n=== 逐 pool: 脂解酶 与 AHR 轴 (检查 WL 内部异质性) ===\n")
print(w[order(group, gsm), .(gsm, group, lipase_core = round(lipase_core, 2),
                             ahr_axis = round(ahr_axis, 2), inflam = round(inflam, 2),
                             sasp = round(sasp, 2))], digits = 2)

## 核心机制检验: 脂解激活是否伴随核受体/外排程序上调
cat("\n=== 机制耦合检验 (18 个 pool 跨样本相关) ===\n")
for (v in c("ahr_axis","car_pxr","nrf2","inflam","sasp","droplet")) {
  ct <- tryCatch(cor.test(w$lipase_core, w[[v]]), error = function(e) NULL)
  if (!is.null(ct))
    cat(sprintf("lipase_core ~ %-11s  r = %+.3f   P = %.4f\n", v, ct$estimate, ct$p.value))
}
cat("\n--- 仅在 WL 组内 (n=6) ---\n")
ww <- w[group == "WL"]
for (v in c("ahr_axis","car_pxr","inflam","sasp")) {
  ct <- tryCatch(cor.test(ww$lipase_core, ww[[v]]), error = function(e) NULL)
  if (!is.null(ct))
    cat(sprintf("lipase_core ~ %-11s  r = %+.3f   P = %.4f\n", v, ct$estimate, ct$p.value))
}

fwrite(w,   file.path(OUT, "geo_adipocyte_subscores.csv"))
fwrite(res, file.path(OUT, "geo_adipocyte_subscore_tests.csv"))
cat("\nALL DONE\n")
