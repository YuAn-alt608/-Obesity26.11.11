## ============================================================================
## 14_geo.R — GSE295708 人腹部皮下脂肪 snRNA-seq (Obese / WeightLoss / Lean)
##   逐样本流式处理以控制内存; 标记基因打分做粗分类; 模块评分 + 假体bulk
## ============================================================================
suppressMessages({library(Matrix); library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
GDIR <- file.path(BASE, "data/raw/geo/extracted")
OUT  <- file.path(BASE, "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

fs <- list.files(GDIR, pattern = "_features", full.names = TRUE)
feat <- fread(fs[1], header = FALSE)
genes <- feat[[2]]
cat("基因数:", length(genes), "\n")

samp <- data.table(
  gsm = sub("_.*", "", basename(list.files(GDIR, pattern = "_matrix"))))
samp[, grp := fifelse(grepl("Weight", gsm), NA_character_, NA_character_)]
## 从样本标题映射分组
tt <- readLines(file.path(BASE, "data/meta/gse295708_samples.txt"), warn = FALSE)
cur <- NA
map <- list()
for (l in tt) {
  if (grepl("^!Sample_title", l)) cur <- trimws(sub("^!Sample_title = ", "", l))
  if (grepl("^!Sample_geo_accession", l)) map[[trimws(sub("^!Sample_geo_accession = ", "", l))]] <- cur
}
samp[, title := unlist(map[gsm])]
samp[, group := fifelse(grepl("Weight Loss", title), "WL",
                fifelse(grepl("Lean", title), "Lean", "Obese"))]
print(samp[, .N, by = group])

MOD <- list(
  lipolysis   = c("PNPLA2","LIPE","ABHD5","PLIN1","PLIN4","MGLL","CIDEC","PLIN2"),
  xenobiotic  = c("AHR","AHRR","ARNT","CYP1A1","CYP1B1","NR1I2","NR1I3","CYP2B6","CYP3A4","ABCB1","ABCG2"),
  inflamma    = c("NLRP3","IL1B","CASP1","CCL2","TNF","IL6","PYCARD","TLR4"),
  senescence  = c("CDKN1A","CDKN2A","TP53","GLB1","SERPINE1","IGFBP3","MKI67"),
  adipogenic  = c("PPARG","CEBPA","FABP4","ADIPOQ","LEP","PLIN1","SLC2A4"),
  antioxid    = c("NQO1","HMOX1","GSTA1","GCLC","GCLM","SOD2","CAT"),
  macrophage  = c("CD68","CD163","LYZ","MRC1","C1QA","ITGAM"),
  adipocyte   = c("ADIPOQ","PLIN1","PLIN4","LEP","FABP4","CIDEC","CAV1"),
  ascell      = c("PDGFRA","DCN","LUM","DLK1","GSN","THY1"),
  endo        = c("PECAM1","VWF","EMCN","CDH5","FLT1"),
  tcell       = c("CD3D","CD3E","CD2","IL7R","NKG7"),
  pericyte    = c("RGS5","ACTA2","PDGFRB","NOTCH3","MYH11"),
  mast        = c("TPSAB1","CPA3","KIT","MS4A2")
)
MOD <- lapply(MOD, function(g) intersect(g, genes))
print(sapply(MOD, length))

CELLS <- list(); PB <- list()
for (i in seq_len(nrow(samp))) {
  g <- samp$gsm[i]
  mf <- list.files(GDIR, pattern = paste0(g, ".*_matrix"), full.names = TRUE)
  bf <- list.files(GDIR, pattern = paste0(g, ".*_barcodes"), full.names = TRUE)
  X <- readMM(mf); bc <- fread(bf, header = FALSE)[[1]]
  rownames(X) <- genes; colnames(X) <- bc
  ## QC
  ncount <- colSums(X); nfeat <- colSums(X > 0)
  mt <- colSums(X[grepl("^MT-", genes), , drop = FALSE]) / pmax(ncount, 1)
  keep <- nfeat >= 300 & nfeat <= 8000 & ncount >= 500 & mt < 20
  X <- X[, keep, drop = FALSE]; ncount <- ncount[keep]
  ## 归一化 CP10k + log1p
  X <- X %*% Diagonal(x = 1e4 / pmax(ncount, 1))
  X@x <- log1p(X@x)
  ## 模块与标记评分 (每基因按本样本细胞 z 标准化后取均值)
  sc <- sapply(MOD, function(gs) {
    S <- as.matrix(X[gs, , drop = FALSE])
    Z <- (S - rowMeans(S)) / pmax(apply(S, 1, sd), 1e-6)
    colMeans(Z)
  })
  ct <- colnames(sc)[max.col(sc[, c("adipocyte","ascell","endo","tcell","macrophage",
                                    "pericyte","mast"), drop = FALSE], ties.method = "first")]
  ct <- c("adipocyte","ascell","endo","tcell","macrophage","pericyte","mast")[
        max.col(sc[, c("adipocyte","ascell","endo","tcell","macrophage","pericyte","mast"),
                   drop = FALSE], ties.method = "first")]
  CELLS[[g]] <- data.table(gsm = g, group = samp$group[i], celltype = ct,
                           nFeature = nfeat[keep],
                           sc[, c("lipolysis","xenobiotic","inflamma","senescence",
                                  "adipogenic","antioxid")])
  ## 假体bulk: 细胞类型内基因计数求和
  for (t in unique(ct)) {
    idx <- which(ct == t)
    if (length(idx) < 30) next
    s <- rowSums(X[, idx, drop = FALSE])
    PB[[paste(g, t)]] <- data.table(gsm = g, group = samp$group[i], celltype = t,
                                    ncell = length(idx),
                                    gene = genes, expr = as.numeric(s))
  }
  cat(sprintf("%s %-6s cells=%d (keep %.0f%%)\n", g, samp$group[i], length(ct),
              100 * length(ct) / length(keep)))
  rm(X); gc()
}

CELL <- rbindlist(CELLS)
PBULK <- rbindlist(PB)
fwrite(CELL, file.path(OUT, "geo_cell_scores.csv.gz"))
fwrite(PBULK, file.path(OUT, "geo_pseudobulk_counts.csv.gz"))

## ---------- 组间比较 ----------
cat("\n=== 细胞组成 (%); Obese vs WL vs Lean ===\n")
comp <- CELL[, .N, by = .(gsm, group, celltype)]
comp[, pct := N / sum(N) * 100, by = gsm]
print(round(dcast(comp[, .(pct = mean(pct)), by = .(group, celltype)],
                  celltype ~ group, value.var = "pct"), 1))
fwrite(comp, file.path(OUT, "geo_composition.csv"))

cat("\n=== 脂肪细胞模块评分 (样本均值) ===\n")
ad <- CELL[celltype == "adipocyte", lapply(.SD, mean), by = .(gsm, group),
           .SDcols = c("lipolysis","xenobiotic","inflamma","senescence","adipogenic","antioxid")]
print(round(dcast(melt(ad, id.vars = c("gsm","group")), variable ~ group,
                  value.var = "value", fun.aggregate = mean), 3))
for (v in c("lipolysis","xenobiotic","inflamma","senescence")) {
  x <- ad[[v]]
  p1 <- t.test(x[ad$group == "WL"], x[ad$group == "Obese"])$p.value
  cat(sprintf("%-11s Obese %.3f -> WL %.3f  (P=%.4f)\n", v,
              mean(x[ad$group == "Obese"]), mean(x[ad$group == "WL"]), p1))
}
fwrite(ad, file.path(OUT, "geo_adipocyte_modules.csv"))

cat("\n=== 关键基因 (脂肪细胞假体bulk, CPM) ===\n")
key <- c("PNPLA2","LIPE","ABHD5","PLIN1","PLIN2","AHR","AHRR","CYP1A1","CYP1B1",
         "NR1I2","NR1I3","CYP2B6","ABCB1","ABCG2","NQO1","HMOX1","NLRP3","IL1B",
         "CCL2","CDKN1A","CDKN2A","ADIPOQ","LEP","PPARG","FABP4","UCP1")
pa <- PBULK[celltype == "adipocyte"]
pa[, cpm := expr / sum(expr) * 1e6, by = .(gsm)]
k <- pa[gene %in% key, .(mean_cpm = mean(cpm)), by = .(gene, group)]
print(round(dcast(k, gene ~ group, value.var = "mean_cpm"), 1))
fwrite(k, file.path(OUT, "geo_key_genes.csv"))
cat("\nALL DONE\n")
