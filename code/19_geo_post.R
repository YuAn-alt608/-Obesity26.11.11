## ============================================================================
## 19_geo_post.R — 从已保存的细胞级结果补算汇总表
##   (14_geo.R 在 print(round(dcast(...))) 处因字符列崩溃, 数据本身已完整写出)
## ============================================================================
suppressMessages({library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
OUT  <- file.path(BASE, "results"); LOG <- file.path(BASE, "logs")

CELL  <- fread(file.path(OUT, "geo_cell_scores.csv.gz"))
PBULK <- fread(file.path(OUT, "geo_pseudobulk_counts.csv.gz"))

cat("=== 细胞总数 ===\n")
print(CELL[, .(ncells = .N), by = .(group)][order(-ncells)])
cat("总细胞:", nrow(CELL), "  样本数:", uniqueN(CELL$gsm), "\n\n")

## ---------- 1. 细胞组成 ----------
comp <- CELL[, .N, by = .(gsm, group, celltype)]
comp[, pct := N / sum(N) * 100, by = gsm]
cat("=== 细胞组成 (%); 各 pool 均值 ===\n")
tab <- dcast(comp[, .(pct = mean(pct)), by = .(group, celltype)],
             celltype ~ group, value.var = "pct")
print(tab[order(-Obese)], digits = 1)

cat("\n=== 组成组间比较 (WL vs Obese, 以 pool 为单位) ===\n")
for (t in unique(comp$celltype)) {
  x <- comp[celltype == t]
  a <- x$pct[x$group == "Obese"]; b <- x$pct[x$group == "WL"]; c3 <- x$pct[x$group == "Lean"]
  p <- tryCatch(t.test(b, a)$p.value, error = function(e) NA_real_)
  cat(sprintf("%-11s Obese %5.1f%%  WL %5.1f%%  Lean %5.1f%%   WL-vs-Obese P=%.4f\n",
              t, mean(a), mean(b), mean(c3), p))
}
fwrite(comp, file.path(OUT, "geo_composition.csv"))

## ---------- 2. 脂肪细胞模块评分 ----------
SC <- c("lipolysis","xenobiotic","inflamma","senescence","adipogenic","antioxid")
ad <- CELL[celltype == "adipocyte", lapply(.SD, mean), by = .(gsm, group), .SDcols = SC]
cat("\n=== 脂肪细胞模块评分 (pool 均值) ===\n")
print(dcast(melt(ad, id.vars = c("gsm","group")), variable ~ group,
            value.var = "value", fun.aggregate = mean), digits = 3)

cat("\n=== 模块评分组间检验 ===\n")
res <- rbindlist(lapply(SC, function(v) {
  x <- ad[[v]]
  a <- x[ad$group == "Obese"]; b <- x[ad$group == "WL"]; c3 <- x[ad$group == "Lean"]
  data.table(module = v,
             Obese = mean(a), WL = mean(b), Lean = mean(c3),
             delta_WL_vs_Obese = mean(b) - mean(a),
             p_WL_vs_Obese = tryCatch(t.test(b, a)$p.value, error = function(e) NA_real_),
             p_WL_vs_Lean  = tryCatch(t.test(b, c3)$p.value, error = function(e) NA_real_))
}))
print(res, digits = 3)
fwrite(ad,  file.path(OUT, "geo_adipocyte_modules.csv"))
fwrite(res, file.path(OUT, "geo_adipocyte_module_tests.csv"))

## 其他细胞类型的模块评分（用于判断"炎症风暴"是否来自巨噬细胞而非脂肪细胞）
cat("\n=== 巨噬细胞模块评分 ===\n")
mac <- CELL[celltype == "macrophage", lapply(.SD, mean), by = .(gsm, group), .SDcols = SC]
print(dcast(melt(mac, id.vars = c("gsm","group")), variable ~ group,
            value.var = "value", fun.aggregate = mean), digits = 3)
fwrite(mac, file.path(OUT, "geo_macrophage_modules.csv"))

## ---------- 3. 关键基因 (脂肪细胞假体 bulk, CPM) ----------
pa <- PBULK[celltype == "adipocyte"]
pa[, cpm := expr / sum(expr) * 1e6, by = gsm]
key <- c("PNPLA2","LIPE","ABHD5","PLIN1","PLIN2","MGLL","CIDEC",
         "AHR","AHRR","ARNT","CYP1A1","CYP1B1","NR1I2","NR1I3","CYP2B6","CYP3A4",
         "ABCB1","ABCG2","NQO1","HMOX1","GSTA1","GCLC","SOD2",
         "NLRP3","IL1B","CASP1","CCL2","TNF","IL6","PYCARD","TLR4",
         "CDKN1A","CDKN2A","TP53","GLB1","SERPINE1",
         "ADIPOQ","LEP","PPARG","FABP4","SLC2A4","UCP1")
k <- pa[gene %in% key, .(mean_cpm = mean(cpm), sd_cpm = sd(cpm)), by = .(gene, group)]
cat("\n=== 脂肪细胞关键基因 CPM ===\n")
w <- dcast(k, gene ~ group, value.var = "mean_cpm")
print(w[order(-Obese)], digits = 1)

cat("\n=== 关键基因 WL vs Obese (基因水平 t 检验, 名义 P) ===\n")
kt <- rbindlist(lapply(intersect(key, unique(pa$gene)), function(g) {
  x <- pa[gene == g]
  a <- x$cpm[x$group == "Obese"]; b <- x$cpm[x$group == "WL"]; c3 <- x$cpm[x$group == "Lean"]
  data.table(gene = g, Obese = mean(a), WL = mean(b), Lean = mean(c3),
             log2FC_WL_Obese = log2((mean(b) + 1) / (mean(a) + 1)),
             p = tryCatch(t.test(log10(b + 1), log10(a + 1))$p.value,
                          error = function(e) NA_real_))
}))
kt[, p_adj := p.adjust(p, "BH")]
print(kt[order(p)], digits = 3)
fwrite(k,  file.path(OUT, "geo_key_genes.csv"))
fwrite(kt, file.path(OUT, "geo_key_gene_tests.csv"))

## ---------- 4. 脂解-外排耦合 (跨 pool 相关) ----------
cat("\n=== 脂肪细胞: 脂解评分 与 外源物代谢/炎症评分的相关 (pool 级) ===\n")
for (v in c("xenobiotic","inflamma","senescence")) {
  ct <- tryCatch(cor.test(ad$lipolysis, ad[[v]]), error = function(e) NULL)
  if (!is.null(ct)) cat(sprintf("lipolysis ~ %-11s r=%+.3f  P=%.4f\n", v, ct$estimate, ct$p.value))
}
cat("\nALL DONE\n")
