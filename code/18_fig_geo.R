## ============================================================================
## 18_fig_geo.R — Fig 7: GSE295708 single-cell mechanism figure  (all-English labels)
##   a) adipocyte pathway sub-scores: Obese / WeightLoss / Lean
##   b) cell composition
##   c) key-gene heatmap (adipocyte pseudobulk)
##   d) coupling test: lipase core vs AHR axis
## ============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
OUT  <- file.path(BASE, "results"); FIG <- file.path(BASE, "figures")

sub  <- fread(file.path(OUT, "geo_adipocyte_subscores.csv"))
comp <- fread(file.path(OUT, "geo_composition.csv"))
kg   <- fread(file.path(OUT, "geo_key_genes.csv"))
st   <- fread(file.path(OUT, "geo_adipocyte_subscore_tests.csv"))

GRP <- c("Obese", "WL", "Lean")
COL <- c(Obese = "#C0504D", WL = "#4F81BD", Lean = "#9BBB59")

## ---------------- a) pathway sub-scores ----------------
LAB <- c(lipase_core = "Lipase core",
         droplet     = "Lipid droplet",
         ahr_axis    = "AHR axis",
         car_pxr     = "CAR / PXR efflux",
         nrf2        = "NRF2 antioxidant",
         inflam      = "Inflammation",
         sasp        = "SASP")
VS <- names(LAB)

ml <- melt(sub[, c("gsm", "group", VS), with = FALSE], id.vars = c("gsm", "group"))
ml[, group := factor(group, GRP)]
ml[, variable := factor(variable, VS, LAB[VS])]

s   <- st[match(VS, set)]
tag <- vapply(seq_along(VS), function(i) {
  p <- s$p_WL_Obese[i]
  mk <- if (VS[i] == "sasp") "*" else ""          # 与逐细胞 z 打分方向相反
  sprintf("%s  (%s)%s", LAB[VS[i]],
          ifelse(is.na(p), "NA", ifelse(p < 0.001, "P<0.001", sprintf("P=%.3f", p))), mk)
}, character(1))
names(tag) <- LAB[VS]

pa <- ggplot(ml, aes(group, value, fill = group)) +
  geom_boxplot(width = .6, outlier.shape = NA, alpha = .8, colour = "grey25") +
  geom_jitter(width = .12, size = 1.35, alpha = .8, colour = "grey10") +
  facet_wrap(~variable, nrow = 1, scales = "free_y",
             labeller = as_labeller(tag, default = label_value)) +
  scale_fill_manual(values = COL, guide = "none") +
  labs(x = NULL, y = expression(log[10]~"CPM (pathway mean)"), tag = "a") +
  theme_bw(base_size = 10) +
  theme(strip.text = element_text(size = 7.1, margin = margin(2, 1, 2, 1)),
        strip.background = element_rect(fill = "grey93", colour = NA),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(.28, "cm"),
        axis.text.x = element_text(size = 8.2),
        plot.tag = element_text(face = "bold", size = 13))

## ---------------- b) cell composition ----------------
cp <- comp[, .(pct = mean(pct)), by = .(group, celltype)]
cp[, group := factor(group, GRP)]
ord <- c("adipocyte", "ascell", "endo", "pericyte", "macrophage", "tcell", "mast")
cp[, celltype := factor(celltype, ord)]
CC <- setNames(c("#7A5195", "#BC5090", "#EF5675", "#FF764A", "#FFA600", "#4F81BD", "#9BBB59"), ord)

pb <- ggplot(cp, aes(group, pct, fill = celltype)) +
  geom_col(width = .64, colour = "white", linewidth = .25) +
  scale_fill_manual(values = CC, name = NULL, guide = guide_legend(ncol = 2, keyheight = .74)) +
  labs(x = NULL, y = "Cell fraction (%)", tag = "b") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "right", legend.text = element_text(size = 7.8),
        axis.text.x = element_text(size = 8.2),
        plot.tag = element_text(face = "bold", size = 13))

## ---------------- c) key-gene heatmap ----------------
GENE <- c("PNPLA2", "LIPE", "MGLL", "ABHD5", "PLIN1", "PLIN2", "CIDEC",
          "AHR", "ARNT", "CYP1B1", "NR1I2", "NR1I3", "ABCB1", "ABCG2",
          "NQO1", "HMOX1", "GCLC", "NLRP3", "IL1B", "CCL2", "IL6",
          "CDKN1A", "TP53", "SERPINE1", "ADIPOQ", "LEP", "UCP1")
h  <- kg[gene %in% GENE]
h[, group := factor(group, GRP)]
hg <- dcast(h, gene ~ group, value.var = "mean_cpm")
hg <- hg[order(match(gene, rev(GENE[GENE %in% hg$gene])))]
M  <- as.matrix(hg[, ..GRP]); rownames(M) <- hg$gene
Z  <- log10(M + 1); Z <- Z - rowMeans(Z)
Zl <- melt(as.data.table(Z, keep.rownames = "gene"), id.vars = "gene")
Zl[, gene := factor(gene, hg$gene)]

pc <- ggplot(Zl, aes(variable, gene, fill = value)) +
  geom_tile(colour = "white", linewidth = .55) +
  geom_text(aes(label = sprintf("%+.1f", value)), size = 2.35, colour = "grey12") +
  scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0,
                       name = "log10 CPM\n(centred)", limits = c(-.45, .45),
                       oob = scales::squish) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL, tag = "c") +
  theme_minimal(base_size = 8.6) +
  theme(panel.grid = element_blank(),
        axis.text.x.top = element_text(size = 8.8, face = "bold"),
        axis.text.y = element_text(size = 7.4, face = "italic"),
        legend.position = "right", legend.key.height = unit(.85, "cm"),
        legend.title = element_text(size = 7.4), legend.text = element_text(size = 7),
        plot.tag = element_text(face = "bold", size = 13))

## ---------------- d) coupling scatter ----------------
sm <- sub[, .(gsm, group, lipase_core, ahr_axis)]
sm[, group := factor(group, GRP)]
ct <- cor.test(sm$lipase_core, sm$ahr_axis)
rlab <- sprintf("r = %+.2f\nP = %.3f", ct$estimate, ct$p.value)

pd <- ggplot(sm, aes(lipase_core, ahr_axis)) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
              colour = "grey35", fill = "grey82", linewidth = .7) +
  geom_point(aes(fill = group), shape = 21, size = 3.1, colour = "grey15", stroke = .5) +
  scale_fill_manual(values = COL, name = NULL) +
  labs(x = expression("Lipase core score ("*log[10]*" CPM)"),
       y = expression("AHR axis score ("*log[10]*" CPM)"),
       subtitle = rlab, tag = "d") +
  theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.6, colour = "grey25", hjust = 0),
        legend.position = c(.17, .14), legend.background = element_blank(),
        legend.key.size = unit(.42, "cm"), legend.text = element_text(size = 8),
        plot.tag = element_text(face = "bold", size = 13))

CAP <- paste(
  "GSE295708 | human abdominal subcutaneous adipose snRNA-seq | Obese n=6 / post-weight-loss n=6 / Lean n=6 pools (4-5 participants per pool) |",
  "186,221 nuclei",
  "\nGene sets - lipase core: PNPLA2, LIPE, MGLL, ABHD5 | droplet: PLIN1/2/4, CIDEC | AHR: AHR, ARNT, CYP1A1, CYP1B1 |",
  "CAR/PXR: NR1I2, NR1I3, CYP2B6, CYP3A4, ABCB1, ABCG2",
  "\nNRF2: NQO1, HMOX1, GCLC, GSTA1 | inflammation: NLRP3, IL1B, CASP1, CCL2, TNF, IL6, PYCARD, TLR4 |",
  "SASP: CDKN1A, CDKN2A, SERPINE1, GLB1, TP53",
  "\n* direction discordant with the cell-level z-score metric; not interpreted (see report section 4.5)")

p <- pa / (pb | pc | pd) +
  plot_layout(heights = c(1, 1.08), widths = c(1, 1.06, 1)) +
  plot_annotation(caption = CAP,
    theme = theme(plot.caption = element_text(size = 6.2, colour = "grey35", hjust = 0, lineheight = 1.25)))

ggsave(file.path(FIG, "fig7_singlecell.png"), p, width = 13.6, height = 9.2, dpi = 300)
ggsave(file.path(FIG, "fig7_singlecell.pdf"), p, width = 13.6, height = 9.2)
cat("saved fig7\n")
