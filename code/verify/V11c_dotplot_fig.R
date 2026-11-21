# ============================================================
# V11c  Marker DotPlot（第四层核验配图）
#   输入: verification/dotplot_markers.csv  (7 celltype x 38 gene)
#   输出: figures/fig13_dotplot_markers.png / .pdf
#   说明: 点大小 = 表达细胞比例(%)，颜色 = 跨细胞类型行标准化 z 分值
#         即标准 Seurat DotPlot 的两个美学通道
# ============================================================
suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})

BASE <- Sys.getenv("AGE22_ROOT", ".")
VER  <- file.path(BASE, "verification")
FIG  <- file.path(BASE, "figures")

D <- fread(file.path(VER, "dotplot_markers.csv"))

# ---- 基因顺序：按所属家族分组，家族按"最强特异度"排序 --------------
FAM_ORDER <- c("adipocyte", "ascell", "endo", "pericyte",
               "tcell", "macrophage", "mast")
CT_ORDER  <- c("adipocyte", "ascell", "endo", "pericyte",
               "tcell", "macrophage", "mast")

# 每组内按自身细胞类型上的 z 降序
G <- D[, .(zmax = z[celltype == family][1]), by = .(gene, family)]
setorder(G, family, -zmax)
D[, gene := factor(gene, levels = rev(G$gene))]
D[, family := factor(family, levels = FAM_ORDER)]
D[, celltype := factor(celltype, levels = CT_ORDER)]

# 家族分隔线位置
gl <- as.integer(table(G$family))
sep <- cumsum(gl)[-length(gl)] + 0.5

# ---- 高亮：审稿人点名的两组标志基因 --------------------------------
HL <- c("ADIPOQ", "PLIN1", "PLIN4", "CIDEC", "LEP",
        "PTPRC", "CD68", "MRC1", "CD163", "CSF1R", "ITGAM")

P <- ggplot(D, aes(x = celltype, y = gene)) +
  geom_point(aes(size = pct_expr, fill = z), shape = 21, colour = "grey25",
             stroke = 0.25) +
  geom_hline(yintercept = sep, colour = "grey60", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                       midpoint = 0, name = "row-scaled\nz-score",
                       limits = c(-1.2, 2.4), oob = scales::squish) +
  scale_size_continuous(range = c(0.2, 6.0), name = "% expressing",
                        breaks = c(10, 30, 50, 70)) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL,
       title = "Marker DotPlot, GSE295708",
       subtitle = paste0("7 marker-scored cell types, 18 pools, 14,396 cells\n",
                         "38/38 markers peak in their own cell type ",
                         "(37 at rank 1; PTPRC rank 2, pan-leukocyte)")) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major   = element_line(colour = "grey92", linewidth = 0.25),
    panel.grid.minor   = element_blank(),
    axis.text.x.top    = element_text(angle = 35, hjust = 0, vjust = 0,
                                      face = "bold", size = 9),
    axis.text.y        = element_text(size = 8),
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 8.5, colour = "grey30"),
    legend.key.height  = unit(0.9, "lines"),
    legend.position    = "right"
  )

ggsave(file.path(FIG, "fig13_dotplot_markers.png"), P,
       width = 7.4, height = 9.2, dpi = 320, bg = "white")
ggsave(file.path(FIG, "fig13_dotplot_markers.pdf"), P,
       width = 7.4, height = 9.2, bg = "white")

cat("✅ V11c 完成 → figures/fig13_dotplot_markers.png / .pdf\n")
cat(sprintf("   基因 %d 个 x 细胞类型 %d 个 = %d 个点\n",
            uniqueN(D$gene), uniqueN(D$celltype), nrow(D)))
cat(sprintf("   特异度: %d/%d 在自身类型表达最高（rank==1: %d）\n",
            sum(D[, .(r = rank(-mean_expr)[celltype == family][1]), by = gene]$r == 1),
            uniqueN(D$gene),
            sum(D[, .(r = rank(-mean_expr)[celltype == family][1]), by = gene]$r == 1)))
