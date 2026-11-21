## ============================================================================
## V12_external_GSE199063.R — 外部队列"背靠背"校验（纵向减重队列）
## ----------------------------------------------------------------------------
## 数据集：GSE199063（人腹部皮下脂肪，女性，Roux-en-Y 胃旁路术）
##   同一批受试者：术前 n=50、术后 2 年 n=49、术后 5 年 n=38
##   另有 28 例年龄匹配、未手术对照（单时点）
##   平台 GPL23126（Clariom D Human），series matrix 为 log2 归一信号
##
## 目的：GSE295708 单细胞得到的三条结论，能否在独立的**个体配对**组织中复现？
##   单细胞结论（WL vs Obese）：
##     S-1 巨噬细胞比例 25.3% -> 13.5%（大幅下降）
##     S-2 留存巨噬细胞每细胞炎症评分翻倍
##     S-3 净炎症负荷（比例 x 每细胞）不变
##   在 bulk 组织中的对应量：
##     E1 mac_score        = 巨噬细胞标志基因模块分（代理巨噬细胞含量）
##     E2 percell_infl     = infl_score - mac_score（对数尺度上剔除含量后的每细胞炎症输出）
##     E3 net_infl         = infl_score（全组织炎症信号 = 含量 x 每细胞）
##   预测：E1 显著下降；E2 显著上升；E3 变化远小于 E1/E2 且不显著
##
## 统计：个体配对 t 检验 + Wilcoxon 符号秩（术前 vs 2 年、术前 vs 5 年）
##       并以未手术对照的单时点水平作参照
## 输出：results/ext_gse199063_paired.csv / ext_gse199063_scores.csv
##       verification/ext_gse199063_map.csv
## ============================================================================
suppressMessages({library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
EXT  <- file.path(BASE, "data/external"); RES <- file.path(BASE, "results")
VER  <- file.path(BASE, "verification")

## ---------- ① 探针 -> 基因符号 ----------------------------------------------
## GPL23126 的 soft 文件无 "Gene Symbol" 列，基因名在 gene_assignment 字段，
## 形如 "NM_015000 // STK38L // serine/threonine kinase 38 like"（" /// " 分隔多条）。
## 已用 awk 流式抽取为 ID<TAB>gene_assignment 两列（GPL23126_symbol.tsv）。
sf <- file.path(EXT, "GPL23126_symbol.tsv")
stopifnot(file.exists(sf))
MAP <- fread(sf, sep = "\t", header = FALSE, quote = "", col.names = c("ID", "assign"))
first_sym <- function(a) {
  a <- sub(" ?///.*$", "", a)              # 只取第一条 assignment
  p <- strsplit(a, " // ", fixed = TRUE)[[1]]
  if (length(p) < 2) return(NA_character_)
  s <- trimws(p[2])
  if (s %in% c("---", "")) return(NA_character_)
  toupper(s)
}
MAP[, symbol := vapply(assign, first_sym, character(1))]
MAP <- MAP[!is.na(symbol)]
MAP <- unique(MAP, by = "ID")
cat("平台行数:", nrow(MAP), " 独立基因:", uniqueN(MAP$symbol), "\n")
fwrite(head(MAP, 50), file.path(VER, "ext_gse199063_map.csv"))

## ---------- ② 表达矩阵 ------------------------------------------------------
fm <- file.path(EXT, "GSE199063_series_matrix.txt.gz")
hdrl <- readLines(gzfile(fm), n = 200, warn = FALSE)
## 样本元数据
getrow <- function(key) {
  l <- grep(paste0("^", key), hdrl, value = TRUE)[1]
  if (is.na(l)) return(NULL)
  v <- strsplit(l, "\t")[[1]][-1]
  gsub('^"|"$', "", v)
}
gsm  <- getrow("!Sample_geo_accession")
sid  <- getrow("!Sample_characteristics_ch1.*subject_id")
## subject_id 不是第一行 characteristics，逐行找
chr_rows <- grep("^!Sample_characteristics_ch1", hdrl, value = TRUE)
pick_chr <- function(tag) {
  for (r in chr_rows) {
    v <- gsub('^"|"$', "", strsplit(r, "\t")[[1]][-1])
    if (any(grepl(paste0("^", tag, ":"), v))) return(sub(paste0("^", tag, ": *"), "", v))
  }
  NULL
}
sid  <- pick_chr("subject_id")
tp   <- pick_chr("time_point")
sexv <- pick_chr("gender")
META <- data.table(gsm = gsm, subject = sid, time_point = tp)
META[, grp := fifelse(time_point == "Control", "Control", "RYGB")]
cat("\n样本:", nrow(META), "\n"); print(META[, .N, by = .(grp, time_point)])

X <- fread(fm, sep = "\t", header = TRUE, skip = "ID_REF", quote = "")
setnames(X, 1, "ID")
## fread(quote="") 会保留字段两端的双引号，需显式剥离
X[, ID := gsub('"', "", ID)]
setnames(X, names(X)[-1], gsub('"', "", names(X)[-1]))
cat("矩阵列示例:", paste(head(names(X)[-1], 3), collapse = ", "), "\n")
X <- X[!grepl("^AFFX", ID)]
cat(sprintf("矩阵探针 %d 个；平台可映射覆盖 %d 个（缺失 %d，%.3f%%）\n",
            nrow(X), sum(X$ID %in% MAP$ID), sum(!X$ID %in% MAP$ID),
            100 * mean(!X$ID %in% MAP$ID)))
X <- merge(MAP, X, by = "ID")
cat("合并注释后:", nrow(X), "探针\n")
## 同基因多探针：取平均表达最高的探针（简洁且常用）
X[, avg := rowMeans(.SD), .SDcols = META$gsm]
setorder(X, symbol, -avg)
X <- unique(X, by = "symbol")
cat("去重后基因:", nrow(X), "\n")

E <- as.matrix(X[, META$gsm, with = FALSE])
stopifnot(ncol(E) == nrow(META), nrow(E) == length(X$symbol))
rownames(E) <- X$symbol
storage.mode(E) <- "double"
cat("表达矩阵:", nrow(E), "基因 x", ncol(E), "样本\n")
cat("NA 比例:", sprintf("%.4f%%", 100*mean(is.na(E))), "\n")
E <- E[!is.na(rownames(E)) & rownames(E) != "", , drop = FALSE]
E <- E[apply(E, 1, function(v) all(is.finite(v))), , drop = FALSE]
cat("可分析基因:", nrow(E), "\n")

## ---------- ③ 模块评分 ------------------------------------------------------
MOD <- list(
  mac = c("CD68","CD163","MRC1","C1QA","C1QB","MS4A4A","F13A1","MSR1","MARCO",
          "ITGAM","LYZ","CSF1R","MRC2","TMEM37","PLA2G7"),
  immune = c("PTPRC","CD53","LAPTM5","CORO1A","ARHGDIB","TYROBP","FCER1G","CD52"),
  inflam = c("IL6","IL1B","TNF","CCL2","CCL3","CCL4","CXCL8","NLRP3","CASP1",
             "PYCARD","TLR4","SOCS3","NOS2","CXCL9","CXCL10","IL18","TNFRSF1B"),
  adipocyte = c("ADIPOQ","PLIN1","PLIN4","FABP4","LEP","CIDEC","CAV1","PPARG",
                "LIPE","PNPLA2","SLC2A4","FABP5"),
  fibro = c("PDGFRA","DCN","LUM","COL1A1","COL3A1","THY1","FBN1")
)
MOD <- lapply(MOD, function(g) intersect(g, rownames(E)))
cat("\n各模块命中基因数:\n"); print(sapply(MOD, length))

## 每个基因在全部样本内 z 标准化，再对模块取均值（模块分单位 = SD）
modscore <- function(genes) {
  S <- E[genes, , drop = FALSE]
  Z <- (S - rowMeans(S)) / pmax(apply(S, 1, sd), 1e-9)
  colMeans(Z)
}
SC <- data.table(gsm = META$gsm, subject = META$subject, time_point = META$time_point,
                 grp = META$grp)
for (m in names(MOD)) SC[, (m) := modscore(MOD[[m]])]
SC[, percell_infl := inflam - mac]          # 对数尺度剔除巨噬含量的"每细胞"炎症
SC[, net_infl := inflam]                   # 全组织炎症信号
SC[, percell_raw := inflam - immune]       # 稳健性：改用泛免疫标志剔除
fwrite(SC, file.path(RES, "ext_gse199063_scores.csv"))

## ---------- ④ 配对检验 ------------------------------------------------------
paired <- function(var, t1, t2) {
  w <- dcast(SC[grp == "RYGB" & time_point %in% c(t1, t2)],
             subject ~ time_point, value.var = var)
  if (!all(c(t1, t2) %in% names(w))) return(NULL)
  d <- w[[t2]] - w[[t1]]; d <- d[is.finite(d)]
  if (length(d) < 5) return(NULL)
  tt <- t.test(d); wt <- suppressWarnings(wilcox.test(d))
  data.table(variable = var, from = t1, to = t2, n_pair = length(d),
             mean_from = mean(w[[t1]], na.rm = TRUE),
             mean_to   = mean(w[[t2]], na.rm = TRUE),
             delta = mean(d), se = sd(d)/sqrt(length(d)),
             pct = 100*(exp(mean(d)) - 1),      # 模块分是 z 尺度，按 log 量级解释
             p_t = tt$p.value, p_wilcox = wt$p.value,
             dz = mean(d)/sd(d))
}
VARS <- c("mac","immune","inflam","percell_infl","net_infl","percell_raw","adipocyte","fibro")
OUT <- rbindlist(lapply(VARS, function(v) rbindlist(list(
  paired(v, "Baseline", "2_year"), paired(v, "Baseline", "5_year")))))

## 多重比较（每个时点内 8 个变量）
OUT[, p_bh := p.adjust(p_t, "BH"), by = to]
setorder(OUT, to, variable)
fwrite(OUT, file.path(RES, "ext_gse199063_paired.csv"))

## ---------- ⑤ 未手术对照的单时点参照 ----------------------------------------
CTRL <- SC[grp == "Control", lapply(.SD, mean), .SDcols = VARS]
BASE0 <- SC[grp == "RYGB" & time_point == "Baseline", lapply(.SD, mean), .SDcols = VARS]
Y5 <- SC[grp == "RYGB" & time_point == "5_year", lapply(.SD, mean), .SDcols = VARS]
CMP <- rbindlist(lapply(VARS, function(v) data.table(
  variable = v, control = CTRL[[v]], baseline = BASE0[[v]], year5 = Y5[[v]],
  d_ctrl_base = BASE0[[v]] - CTRL[[v]], d_ctrl_y5 = Y5[[v]] - CTRL[[v]])))
fwrite(CMP, file.path(RES, "ext_gse199063_vs_control.csv"))

## ---------- ⑥ 报告 ----------------------------------------------------------
cat("\n================ 配对分析（术前 vs 术后）================\n")
print(OUT[, .(variable, to, n_pair = n_pair,
              base = round(mean_from, 3), post = round(mean_to, 3),
              delta = round(delta, 3), dz = round(dz, 2),
              p_t = signif(p_t, 3), p_bh = signif(p_bh, 3))], row.names = FALSE)

cat("\n================ 与未手术对照的比较 ================\n")
print(CMP[, .(variable, control = round(control,3), baseline = round(baseline,3),
              year5 = round(year5,3), d_base = round(d_ctrl_base,2),
              d_y5 = round(d_ctrl_y5,2))], row.names = FALSE)

cat("\n================ 三条预测的判定 ================\n")
chk <- function(v, tpv, dir) {
  r <- OUT[variable == v & to == tpv]
  if (!nrow(r)) return("NO DATA")
  ok <- if (dir == "up") (r$delta > 0 && r$p_t < 0.05) else
        if (dir == "down") (r$delta < 0 && r$p_t < 0.05) else (r$p_t >= 0.05)
  sprintf("%s (delta=%+.3f, P=%.4g)", if (ok) "支持" else "不支持", r$delta, r$p_t)
}
for (tpv in c("2_year", "5_year")) {
  cat(sprintf("\n[术前 -> %s]\n", tpv))
  cat("  E1 巨噬含量下降    :", chk("mac", tpv, "down"), "\n")
  cat("  E2 每细胞炎症上升  :", chk("percell_infl", tpv, "up"), "\n")
  cat("  E3 净炎症负荷下降  :", chk("net_infl", tpv, "down"), "\n")
  cat("  (稳健性) 泛免疫剔除:", chk("percell_raw", tpv, "up"), "\n")
  d1 <- OUT[variable == "mac" & to == tpv]$delta
  d2 <- OUT[variable == "percell_infl" & to == tpv]$delta
  d3 <- OUT[variable == "net_infl" & to == tpv]$delta
  cat(sprintf("  分解：net = mac + percell = %+.3f + %+.3f = %+.3f（自洽）\n", d1, d2, d3))
  cat(sprintf("  抵消比例：若每细胞炎症不变，净负荷应降 %+.3f；实测 %+.3f\n", d1, d3))
  cat(sprintf("            → 每细胞炎症上升抵掉了巨噬撤退效应的 %.1f%%\n",
              100 * (1 - d3/d1)))
}

## 效应量对照
cmp <- dcast(OUT, variable ~ to, value.var = "delta")
cat("\n效应量对照（|delta|，SD 单位）:\n"); print(cmp)
fwrite(cmp, file.path(VER, "ext_gse199063_effect_size.csv"))
cat("\n✅ V12 完成\n")
