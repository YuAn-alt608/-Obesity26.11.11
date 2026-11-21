## ============================================================================
## V11b_cluster_ari.R — 第四层②③：无监督聚类一致性 + Marker DotPlot
##   问题：报告 §4.1 用"标记基因 z 打分取最大"把细胞粗分为 7 类，
##         审稿人会问：这个分群靠不靠得住？
##   做法：
##     ① 重跑标记基因 z 打分（复刻 14_geo.R 算法），与已存 CSV 标签逐细胞核对
##     ② 在分层抽样的细胞子集上跑标准无监督流程：
##        HVG(2000) -> z-score -> PCA(30) -> ① KMeans(k=7) ② Louvain(kNN k=20)
##     ③ 计算无监督标签与 z 打分标签的调整兰德指数 ARI
##     ④ 输出 Marker DotPlot 数据（平均表达 + 表达细胞比例）
##   判据：ARI > 0.5 视为两种分群高度一致；> 0.7 视为几乎等价
## ============================================================================
suppressMessages({library(Matrix); library(data.table); library(irlba)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
GDIR <- file.path(BASE, "data/raw/geo/extracted")
RES  <- file.path(BASE, "results"); VER <- file.path(BASE, "verification")

set.seed(20261117)
NPER <- 800L        # 每池抽样细胞数
K    <- 7L
NHVG <- 2000L
NPC  <- 30L

## ---------- 样本表 ----------------------------------------------------------
fs <- list.files(GDIR, pattern = "_features", full.names = TRUE)
genes <- fread(fs[1], header = FALSE)[[2]]
cat("基因数:", length(genes), "\n")

tt <- readLines(file.path(BASE, "data/meta/gse295708_samples.txt"), warn = FALSE)
cur <- NA; map <- list()
for (l in tt) {
  if (grepl("^!Sample_title", l)) cur <- trimws(sub("^!Sample_title = ", "", l))
  if (grepl("^!Sample_geo_accession", l)) map[[trimws(sub("^!Sample_geo_accession = ", "", l))]] <- cur
}
gsms <- sort(sub("_.*", "", basename(list.files(GDIR, pattern = "_matrix"))))
samp <- data.table(gsm = gsms)
samp[, title := unlist(map[gsm])]
samp[, group := fifelse(grepl("Weight Loss", title), "WL",
                fifelse(grepl("Lean", title), "Lean", "Obese"))]
print(samp[, .N, by = group])

MOD <- list(
  macrophage = c("CD68","CD163","LYZ","MRC1","C1QA","ITGAM"),
  adipocyte  = c("ADIPOQ","PLIN1","PLIN4","LEP","FABP4","CIDEC","CAV1"),
  ascell     = c("PDGFRA","DCN","LUM","DLK1","GSN","THY1"),
  endo       = c("PECAM1","VWF","EMCN","CDH5","FLT1"),
  tcell      = c("CD3D","CD3E","CD2","IL7R","NKG7"),
  pericyte   = c("RGS5","ACTA2","PDGFRB","NOTCH3","MYH11"),
  mast       = c("TPSAB1","CPA3","KIT","MS4A2")
)
MOD <- lapply(MOD, function(g) intersect(g, genes))
CT7 <- names(MOD)

## Marker DotPlot 基因（经典特异标记，覆盖 7 类）
DOT <- list(
  adipocyte = c("ADIPOQ","PLIN1","PLIN4","LEP","FABP4","CIDEC","CAV1","LIPE"),
  ascell    = c("PDGFRA","DCN","LUM","DLK1","GSN"),
  endo      = c("PECAM1","VWF","EMCN","CDH5"),
  pericyte  = c("RGS5","ACTA2","PDGFRB","NOTCH3"),
  tcell     = c("CD3D","CD3E","CD2","IL7R","NKG7","PTPRC"),
  macrophage= c("CD68","CD163","MRC1","LYZ","C1QA","ITGAM","CSF1R","PTPRC"),
  mast      = c("TPSAB1","CPA3","KIT","MS4A2")
)
DOTG <- unique(unlist(DOT)); DOTG <- intersect(DOTG, genes)
cat("DotPlot 基因数:", length(DOTG), "\n")

## ---------- 通行 1：QC + 重算 z 打分标签 + DotPlot 统计 ----------------------
CACHE <- file.path(VER, "_tmp", "v11b"); dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)
SUB <- list(); DOTSTAT <- list(); LABCHK <- list()
for (i in seq_len(nrow(samp))) {
  g <- samp$gsm[i]
  cf <- file.path(CACHE, paste0(g, ".rds"))
  if (file.exists(cf)) {
    z <- readRDS(cf)
    SUB[[g]] <- z$sub; DOTSTAT[[g]] <- z$dot
    cat(sprintf("%s %-5s [cached] sub=%4d  mac=%4.1f%%\n", g, z$grp,
                ncol(z$sub$X), 100 * mean(z$sub$ct == "macrophage")))
    next
  }
  mf <- list.files(GDIR, pattern = paste0(g, ".*_matrix"), full.names = TRUE)
  bf <- list.files(GDIR, pattern = paste0(g, ".*_barcodes"), full.names = TRUE)
  ## 注意：本数据集部分汇池的 barcode 数达 6.8e6（绝大多数为空液滴）。
  ## readMM 返回三元组 dgTMatrix（无 @p 槽，按列子集极慢），必须先转 CSC。
  ## 同时：**先做 QC 与列子集、再挂 dimnames**——若在 6.8e6 列上先挂 barcode，
  ## 单是 680 万条字符串就吃掉约 1 GB，叠加矩阵本身会触发换页。
  X <- as(readMM(mf), "CsparseMatrix")
  bc <- fread(bf, header = FALSE)[[1]]
  ncount <- Matrix::colSums(X)
  nfeat  <- diff(X@p)          # 计数矩阵非负 ⇒ 每列非零数 = diff(@p)，O(列数)
  mt <- Matrix::colSums(X[which(grepl("^MT-", genes)), , drop = FALSE]) / pmax(ncount, 1)
  keep <- which(nfeat >= 300 & nfeat <= 8000 & ncount >= 500 & mt < 20)
  cat(sprintf("    barcode=%d -> kept=%d\n", ncol(X), length(keep)))
  X <- X[, keep, drop = FALSE]; ncount <- ncount[keep]; nfeat <- nfeat[keep]
  bc <- bc[keep]
  rownames(X) <- genes; colnames(X) <- bc
  X <- X %*% Diagonal(x = 1e4 / pmax(ncount, 1))
  X@x <- log1p(X@x)

  ## 复刻 14_geo.R 的 z 打分分类
  sc <- sapply(MOD, function(gs) {
    S <- as.matrix(X[gs, , drop = FALSE])
    Z <- (S - rowMeans(S)) / pmax(apply(S, 1, sd), 1e-6)
    colMeans(Z)
  })
  ct <- CT7[max.col(sc[, CT7, drop = FALSE], ties.method = "first")]

  ## DotPlot：平均表达 + 表达比例（全 QC 后细胞）
  dg <- intersect(DOTG, rownames(X))
  Dm <- as.matrix(X[dg, , drop = FALSE]); Dp <- as.matrix(X[dg, , drop = FALSE] > 0)
  DOTSTAT[[g]] <- rbindlist(lapply(CT7, function(t) {
    idx <- which(ct == t); if (length(idx) < 1) return(NULL)
    data.table(gsm = g, group = samp$group[i], celltype = t, ncell = length(idx),
               gene = dg,
               mean_expr = rowMeans(Dm[, idx, drop = FALSE]),
               pct_expr  = rowMeans(Dp[, idx, drop = FALSE]) * 100)
  }))

  ## 抽样（分层于细胞类型，保证稀有类不被抽没）
  n <- min(NPER, length(ct))
  idx <- unlist(lapply(CT7, function(t) {
    ii <- which(ct == t)
    if (!length(ii)) return(integer(0))
    take <- max(1L, round(n * length(ii) / length(ct)))
    sample(ii, min(take, length(ii)))
  }))
  idx <- idx[seq_len(min(length(idx), n))]
  SUB[[g]] <- list(X = X[, idx, drop = FALSE], ct = ct[idx],
                   cellid = paste0(g, "_", bc[idx]), group = samp$group[i])
  LABCHK[[g]] <- data.table(gsm = g, cellid = paste0(g, "_", bc),
                            ct_recomputed = ct)
  saveRDS(list(sub = SUB[[g]], dot = DOTSTAT[[g]], lab = LABCHK[[g]],
               grp = samp$group[i]), cf)
  cat(sprintf("%s %-5s kept=%5d  sub=%4d  mac=%4.1f%%\n", g, samp$group[i],
              length(ct), length(idx), 100 * mean(ct[idx] == "macrophage")))
  rm(X, Dm, Dp); gc()
}

## ---------- ① 标签复刻核对 --------------------------------------------------
## 已存 geo_cell_scores.csv 不含 barcode，无法逐细胞对齐；改用"池内细胞类型比例
## 一致性"检验：把本次复刻的（抽样）比例与已存全量比例逐池逐类比对。
## 抽样 800 细胞/池时，比例的抽样标准误约 1–2 个百分点，故 |差| < 5 pp 视为一致。
cat("\n================ ① 标记基因 z 打分复刻核对 ================\n")
cmp0 <- fread(file.path(RES, "geo_composition.csv"))
sub_pct <- rbindlist(lapply(SUB, function(s) {
  data.table(gsm = sub("_.*", "", s$cellid[1]), ct = s$ct)
}))[, .N, by = .(gsm, ct)][, pct_sub := 100 * N / sum(N), by = gsm]
LK <- merge(cmp0[, .(gsm, ct = celltype, pct_full = pct)],
            sub_pct[, .(gsm, ct, pct_sub)], by = c("gsm", "ct"), all = TRUE)
LK[is.na(pct_sub), pct_sub := 0][is.na(pct_full), pct_full := 0]
LK[, diff_pp := pct_sub - pct_full]
cat(sprintf("比对 (池 x 细胞类型) 组合数: %d\n", nrow(LK)))
cat(sprintf("比例差 |diff|: 中位 %.2f pp, P90 %.2f pp, 最大 %.2f pp\n",
            median(abs(LK$diff_pp)), quantile(abs(LK$diff_pp), 0.9), max(abs(LK$diff_pp))))
cat(sprintf("相关系数 r = %.4f\n", cor(LK$pct_full, LK$pct_sub)))
cat(sprintf("判定：%s（阈值：|diff| 最大 < 5 pp）\n",
            if (max(abs(LK$diff_pp)) < 5) "通过" else "需复核"))
fwrite(LK, file.path(VER, "label_reproduction.csv"))

## ---------- 合并稀疏矩阵 ----------------------------------------------------
X <- do.call(cbind, lapply(SUB, `[[`, "X"))
CT <- unlist(lapply(SUB, `[[`, "ct"))
CID <- unlist(lapply(SUB, `[[`, "cellid"))
GRP <- unlist(lapply(SUB, function(s) rep(s$group, ncol(s$X))))
names(CT) <- CID
cat("\n合并矩阵:", nrow(X), "基因 x", ncol(X), "细胞; nnz =",
    format(length(X@x), big.mark = ","), "\n")
print(table(CT, GRP))

## ---------- ② HVG + PCA -----------------------------------------------------
mu  <- Matrix::rowMeans(X)
mu2 <- Matrix::rowMeans(X^2)
vr  <- pmax(mu2 - mu^2, 0)
## 表达基因 & 去线粒体/核糖体
ok <- mu > 0.01 & !grepl("^MT-|^RP[SL]|^HB[AB]", rownames(X))
bin <- cut(mu, quantile(mu[ok], seq(0, 1, length.out = 21)), include.lowest = TRUE,
           labels = FALSE)
zdisp <- ave(sqrt(vr), bin, FUN = function(v) (v - mean(v)) / pmax(sd(v), 1e-9))
hvg <- rownames(X)[ok][order(-zdisp[ok])[seq_len(NHVG)]]
cat("HVG:", length(hvg), "\n")

Xs <- as.matrix(X[hvg, ])
Xs <- t(scale(t(Xs))); Xs[!is.finite(Xs)] <- 0
Xs[Xs > 10] <- 10
pc <- prcomp_irlba(t(Xs), n = NPC, center = FALSE, scale. = FALSE)
PCS <- pc$x; rownames(PCS) <- CID
ve <- pc$sdev^2 / sum(pc$sdev^2) * 100
cat(sprintf("PC1-30 累计方差: %.1f%%\n", sum(ve)))
rm(Xs); gc()

## ---------- ③ KMeans(k=7) ---------------------------------------------------
set.seed(20261117)
km <- kmeans(PCS, centers = K, nstart = 25, iter.max = 100)
## ---------- ④ Louvain -------------------------------------------------------
suppressMessages(library(igraph)); suppressMessages(library(FNN))
nn <- get.knn(PCS, k = 20)
ed <- data.table(from = rep(seq_len(nrow(PCS)), 20), to = as.vector(nn$nn.index))
ed[, w := 1 / (1 + as.vector(nn$nn.dist))]
g2 <- graph_from_data_frame(ed, directed = FALSE, vertices = data.frame(name = seq_len(nrow(PCS))))
E(g2)$weight <- ed$w
lv <- cluster_louvain(g2, resolution = 1)
lvlab <- as.integer(membership(lv))
cat("Louvain 得到群数:", length(unique(lvlab)), "\n")

## ---------- ⑤ ARI -----------------------------------------------------------
ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  c2 <- function(z) sum(choose(z, 2))
  sij <- sum(choose(tab, 2))
  si <- c2(rowSums(tab)); sj <- c2(colSums(tab))
  expct <- si * sj / choose(n, 2)
  mx <- (si + sj) / 2
  (sij - expct) / (mx - expct)
}
## 匈牙利式最优标签映射（用贪婪 + 列联表最大化；K 小，暴力全排列）
best_map <- function(pred, truth) {
  lev <- sort(unique(truth))
  gs <- expand.grid(rep(list(unique(pred)), length(lev)))
  gs <- gs[apply(gs, 1, function(r) length(unique(r)) == length(lev)), , drop = FALSE]
  sc <- apply(gs, 1, function(r) mean(truth == lev[as.integer(r)]))
  gs[which.max(sc), ]
}
REP <- data.table(
  comparison = c("KMeans(k=7) vs z-score", "Louvain vs z-score", "KMeans vs Louvain"),
  ari = c(ari(km$cluster, CT), ari(lvlab, CT), ari(km$cluster, lvlab))
)
REP[, verdict := fifelse(ari > 0.7, "几乎等价", fifelse(ari > 0.5, "高度一致",
                    fifelse(ari > 0.3, "中等一致", "弱一致")))]
cat("\n================ ② 无监督聚类 vs 标记基因 z 打分 ================\n")
print(REP)
fwrite(REP, file.path(VER, "cluster_ari.csv"))

## 列联表（KMeans 标签 vs 标记标签），用于判断是否只是标签重排
tab <- as.data.table(as.data.frame.matrix(table(km = km$cluster, marker = CT)),
                     keep.rownames = "km")
fwrite(tab, file.path(VER, "cluster_contingency_kmeans.csv"))
cat("\nKMeans x 标记基因 列联表：\n"); print(tab)

## 每群的主导标记类型占比（purity）
pur <- data.table(km = km$cluster, ct = CT)[, .N, by = .(km, ct)]
pur[, frac := N / sum(N), by = km]
best <- pur[, .SD[which.max(frac)], by = km]
cat("\n每群主导细胞类型（纯度）：\n")
print(best[order(km), .(km, dominant = ct, frac = round(frac, 3))])
cat(sprintf("总体纯度（按群大小加权）: %.3f\n", sum(best$frac * best$N) / sum(best$N)))
fwrite(best, file.path(VER, "cluster_purity.csv"))

## ---------- ⑥ DotPlot 数据 --------------------------------------------------
DS <- rbindlist(DOTSTAT)
dot <- DS[, .(mean_expr = mean(mean_expr), pct_expr = mean(pct_expr)),
          by = .(celltype, gene)]
dot[, z := (mean_expr - mean(mean_expr)) / pmax(sd(mean_expr), 1e-9), by = gene]
dot[, family := sapply(gene, function(g) {
  hit <- names(DOT)[sapply(DOT, function(v) g %in% v)]
  paste(hit, collapse = "/")
})]
fwrite(dot, file.path(VER, "dotplot_markers.csv"))
cat("\n================ ③ Marker DotPlot（平均表达，按基因 z 标准化）================\n")
DP <- dcast(dot, gene + family ~ celltype, value.var = "z")
ord <- c("adipocyte","ascell","endo","pericyte","tcell","macrophage","mast")
show <- intersect(c("gene","family", ord), names(DP))
print(DP[order(-abs(adipocyte))][seq_len(min(30, .N)), ..show], digits = 2)

## 特异度检验：每个家族标记在"自己的"细胞类型上是否最高
spec <- dot[, {
  fam <- strsplit(family, "/")[[1]]
  own <- mean_expr[match(intersect(fam, celltype), celltype)]
  data.table(own_type = paste(intersect(fam, celltype), collapse = "|"),
             own_expr = mean(own), max_other = max(mean_expr[!celltype %in% fam]),
             rank = rank(-mean_expr)[celltype %in% fam][1])
}, by = gene]
spec[, pass := own_expr > max_other]
cat(sprintf("\n标记基因特异度: %d/%d 在自身细胞类型上表达最高\n",
            sum(spec$pass), nrow(spec)))
print(spec[pass == FALSE])
fwrite(spec, file.path(VER, "dotplot_specificity.csv"))

saveRDS(list(ari = REP, purity = best, pcs = PCS, ct = CT, grp = GRP,
             km = km$cluster, lv = lvlab),
        file.path(VER, "cluster_objects.rds"))
cat("\n✅ V11b 完成 → cluster_ari.csv, cluster_purity.csv, dotplot_markers.csv,\n")
cat("             dotplot_specificity.csv, label_reproduction.csv\n")
