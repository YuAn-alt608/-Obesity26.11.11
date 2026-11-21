## ============================================================================
## 33_singlecell_lr.R — GSE295708：巨噬细胞"数量减半、残余亚群信号改变"的解析
##   ① M1/M2 极化（受体/活化标志）
##   ② 促炎输出（配体：IL6/CCL2/TNF/IL1B）——与 z 打分口径交叉核对
##   ③ 配体-受体对（人工校正方向，仅保留经典配对；非 CellChat）
##   ④ 净炎症负荷
##   分析单位：池（每组 6 个）
## ============================================================================
suppressMessages({library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", "."); RES <- file.path(BASE, "results")

pb <- fread(file.path(RES, "geo_pseudobulk_counts.csv.gz"))
pb[, lib := sum(expr), by = .(gsm, celltype)]
pb[, lcpm := log10(expr / lib * 1e6 + 1)]
agg <- pb[, .(lcpm = mean(lcpm)), by = .(gsm, celltype, gene)]
M <- dcast(agg, gene ~ gsm + celltype, value.var = "lcpm")
G <- M$gene; M <- as.matrix(M[, -1]); rownames(M) <- G
cat("矩阵:", nrow(M), "基因 x", ncol(M), "样本×细胞类型\n")
GRP <- unique(pb[, .(gsm, group)])[order(gsm)]

E <- function(gsm, ct, genes) {
  cols <- paste0(gsm, "_", ct); genes <- intersect(genes, rownames(M))
  if (!length(genes)) return(rep(NA_real_, length(gsm)))
  colMeans(M[genes, cols, drop = FALSE], na.rm = TRUE)
}
cmp_grp <- function(v, lab) {
  a <- v[GRP$group=="WL"]; b <- v[GRP$group=="Obese"]; c0 <- v[GRP$group=="Lean"]
  cat(sprintf("%-22s Obese=%.3f  WL=%.3f  Lean=%.3f | WL-Obese P=%.4f  WL-Lean P=%.4f  WL-Obese%%=%+5.1f\n",
              lab, mean(b), mean(a), mean(c0), t.test(a,b)$p.value, t.test(a,c0)$p.value,
              100*(10^(mean(a)-mean(b))-1)))
}

## ---------- ① M1/M2 ----------------------------------------------------------
M1 <- c("CD80","CD86","NOS2","CXCL9","CXCL10","SOCS3","TLR2","TLR4")
M2 <- c("MRC1","CD163","IL10","STAB1","F13A1","PDGFB","CCL18","MAF")
vM1 <- E(GRP$gsm, "macrophage", M1); vM2 <- E(GRP$gsm, "macrophage", M2)
cat("\n=== ① 巨噬细胞活化表型（log10 CPM 基因均值） ===\n")
cmp_grp(vM1, "M1 markers")
cmp_grp(vM2, "M2 markers")
cmp_grp(vM1 - vM2, "M1 - M2")

## ---------- ② 促炎输出（配体口径） -------------------------------------------
INFL <- c("IL6","CCL2","TNF","IL1B","CCL3","CCL4","CXCL8","NLRP3","CASP1","PYCARD")
vINF <- E(GRP$gsm, "macrophage", INFL)
cat("\n=== ② 巨噬细胞促炎输出（配体/炎症小体，log10 CPM） ===\n")
cmp_grp(vINF, "Inflammatory output")
cat("\n单基因（WL vs Obese，%）：\n")
for (g in INFL) {
  a <- E(GRP$gsm, "macrophage", g)
  cat(sprintf("  %-8s Obese=%.3f WL=%.3f  %+6.1f%%  P=%.3f\n", g,
              mean(a[GRP$group=="Obese"]), mean(a[GRP$group=="WL"]),
              100*(10^(mean(a[GRP$group=="WL"]) - mean(a[GRP$group=="Obese"]))-1),
              tryCatch(t.test(a[GRP$group=="WL"], a[GRP$group=="Obese"])$p.value,
                       error=function(e) NA)))
}

## ---------- ③ 配体-受体（方向人工校正） --------------------------------------
LR <- data.table(
  ligand   = c("CCL2","CSF1","LEP","ADIPOQ","RETN","C3","CXCL12","LIF",
               "TNF","IL1B","IL6","SPP1","TGFB1"),
  receptor = c("CCR2","CSF1R","LEPR","ADIPOR1","CAP1","C3AR1","CXCR4","LIFR",
               "TNFRSF1A","IL1R1","IL6R","CD44","TGFBR1"),
  dir      = c(rep("adipo->mac", 8), rep("mac->adipo", 5))
)
lrp <- rbindlist(lapply(seq_len(nrow(LR)), function(i) {
  s <- if (LR$dir[i] == "adipo->mac") "adipocyte" else "macrophage"
  r <- if (LR$dir[i] == "adipo->mac") "macrophage" else "adipocyte"
  lig <- E(GRP$gsm, s, LR$ligand[i]); rec <- E(GRP$gsm, r, LR$receptor[i])
  data.table(pair = paste0(LR$ligand[i], " -> ", LR$receptor[i]),
             dir = LR$dir[i], gsm = GRP$gsm, group = GRP$group,
             score = (lig + rec)/2, lig = lig, rec = rec)
}))
lrp <- lrp[!is.na(score)]
summ <- lrp[, {
  a <- score[group=="WL"]; b <- score[group=="Obese"]; c0 <- score[group=="Lean"]
  data.table(obese=mean(b), wl=mean(a), lean=mean(c0), pct=100*(10^(mean(a)-mean(b))-1),
             p=tryCatch(t.test(a,b)$p.value, error=function(e) NA_real_))
}, by = .(pair, dir)]
summ[, p_bh := p.adjust(p, "BH")]; setorder(summ, p)
cat("\n=== ③ 配体-受体活性（WL vs Obese） ===\n")
print(summ[, .(pair, dir, obese=round(obese,2), wl=round(wl,2), lean=round(lean,2),
               pct=round(pct,1), p=signif(p,2), p_bh=signif(p_bh,2))], nrows=20)

cat("\n=== 配体与受体分别的变化（%） ===\n")
lg <- lrp[, .(wl = mean(lig[group=="WL"]), ob = mean(lig[group=="Obese"])), by = pair]
rc <- lrp[, .(wl = mean(rec[group=="WL"]), ob = mean(rec[group=="Obese"])), by = pair]
x <- merge(lg[, .(pair, lig_pct = 100*(10^(wl-ob)-1))], rc[, .(pair, rec_pct = 100*(10^(wl-ob)-1))], by="pair")
setorder(x, -lig_pct)
print(x[, .(pair, ligand=round(lig_pct,1), receptor=round(rec_pct,1))])

## ---------- ④ 净炎症负荷 -----------------------------------------------------
cmp <- fread(file.path(RES, "geo_composition.csv"))
mac <- fread(file.path(RES, "geo_macrophage_modules.csv"))
b <- merge(cmp[celltype=="macrophage", .(gsm, group, frac = pct/100)],
           mac[, .(gsm, percell_z = inflamma)], by = "gsm")
b[, burden := frac * percell_z]
cat("\n=== ④ 净炎症负荷（比例 × 每细胞 z 评分） ===\n")
print(b[, .(mean=round(mean(burden),4), se=round(sd(burden)/sqrt(.N),4)), by=group])
cat("P(WL vs Obese) =", signif(t.test(b[group=="WL"]$burden, b[group=="Obese"]$burden)$p.value,3), "\n")

cat("\n=== 口径对照：同一批基因，z 打分 vs log10CPM ===\n")
cat(sprintf("  逐细胞 z 打分（14_geo.R 口径）: Obese %.3f -> WL %.3f （%+0.1f%%）\n",
            mean(mac[group=="Obese"]$inflamma), mean(mac[group=="WL"]$inflamma),
            100*(mean(mac[group=="WL"]$inflamma)/mean(mac[group=="Obese"]$inflamma)-1)))
cat(sprintf("  假体 bulk log10CPM（本脚本口径）: Obese %.3f -> WL %.3f （%+0.1f%%）\n",
            mean(vINF[GRP$group=="Obese"]), mean(vINF[GRP$group=="WL"]),
            100*(10^(mean(vINF[GRP$group=="WL"]) - mean(vINF[GRP$group=="Obese"]))-1)))

pol <- data.table(gsm=GRP$gsm, group=GRP$group, M1=vM1, M2=vM2, inflam_output=vINF)
fwrite(pol, file.path(RES, "sc_m1m2.csv"))
fwrite(summ, file.path(RES, "sc_lr_pairs.csv"))
fwrite(x, file.path(RES, "sc_lig_rec_pct.csv"))
fwrite(b, file.path(RES, "sc_burden.csv"))
cat("\n✅ 完成\n")
