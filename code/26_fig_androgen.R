## ============================================================================
## 26_fig_androgen.R — 图 8：雄激素轴与 PFAS 效应修饰（全英文标注）
##   a 每 5% 减重对各雄激素结局的效应
##   b 四个核心结局在 PFAS P10/P50/P90 的条件效应（交互全阴性）
##   c 单味同系物的标准化交互（无亲和力梯度）
##   d 森林异质性检验的随机种子稳定性（对数刻度）
## ============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
RES  <- file.path(BASE, "results"); FIG <- file.path(BASE, "figures")

per5 <- fread(file.path(RES, "androgen_per5.csv"))
ax   <- fread(file.path(RES, "androgen_axis.csv"))
cg   <- fread(file.path(RES, "androgen_congener.csv"))
sdd  <- fread(file.path(RES, "forest_seed_stability.csv"))

COL <- c("#2C6E8F", "#C4623C", "#6B8E5A")
TH  <- theme_bw(base_size = 9.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = .3),
        plot.tag = element_text(face = "bold", size = 12),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8.1, colour = "grey35"),
        legend.key.size = unit(.36, "cm"))

ORD <- c("Total testosterone","SHBG","Free testosterone","%Free testosterone",
         "FAI (T/SHBG)","Estradiol","T:E2 ratio","Albumin (control)")

## ---------- a ---------------------------------------------------------------
pd4 <- per5[match(ORD, label)]
pd4[, se := abs(per5) / qnorm(1 - p/2)]
pd4[, `:=`(est = 100*(exp(per5)-1),
           lo  = 100*(exp(per5 - 1.96*se)-1),
           hi  = 100*(exp(per5 + 1.96*se)-1))]
pd4[, label := factor(label, rev(ORD))]
pd4[, ctrl := label == "Albumin (control)"]
pd4[, sig := ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))]

pa <- ggplot(pd4, aes(est, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .17, colour = "grey40", linewidth = .55) +
  geom_point(aes(fill = ctrl), shape = 21, size = 2.6, colour = "grey20", stroke = .5) +
  geom_text(aes(x = max(pd4$hi) + 0.9, label = sig), hjust = 0, size = 3.0, colour = "grey25") +
  scale_fill_manual(values = c(`FALSE` = COL[1], `TRUE` = "white"), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(.10, .16))) +
  labs(x = "Change per 5% body-weight loss (%)", y = NULL,
       title = "a  Androgen-axis response to weight loss",
       subtitle = "N = 1,232 men (NHANES 2013-2016); drug-induced androgen suppression excluded") +
  TH

## ---------- b ---------------------------------------------------------------
B4  <- c("Total testosterone","SHBG","Free testosterone","%Free testosterone")
a4  <- ax[match(B4, label), .(label, p10, se10, p50, se50, p90, se90, p_joint)]
cb  <- rbind(a4[, .(label, lev = "P10", est = p10, se = se10)],
             a4[, .(label, lev = "P50", est = p50, se = se50)],
             a4[, .(label, lev = "P90", est = p90, se = se90)])
cb[, `:=`(lev = factor(lev, c("P10","P50","P90")),
          pct = 100*(exp(est)-1),
          lo  = 100*(exp(est - 1.96*se)-1),
          hi  = 100*(exp(est + 1.96*se)-1))]
cb[, label := factor(label, B4)]

pb <- ggplot(cb, aes(pct, lev, colour = lev)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .14, linewidth = .5) +
  geom_point(size = 2.2) +
  geom_line(aes(group = label), colour = "grey62", linewidth = .35) +
  facet_wrap(~label, nrow = 2, scales = "free_x") +
  scale_colour_manual(values = c(COL[1], "grey45", COL[2]), guide = "none") +
  scale_y_discrete(limits = rev(c("P10","P50","P90"))) +
  scale_x_continuous(expand = expansion(mult = c(.14, .14))) +
  labs(x = "Weight-loss effect (%)", y = "PFAS percentile",
       title = "b  Effect by PFAS level",
       subtitle = sprintf("no attenuation: interaction P = %.2f-%.2f",
                          min(a4$p_joint), max(a4$p_joint))) +
  TH + theme(strip.text = element_text(size = 8),
             strip.background = element_rect(fill = "grey94", colour = NA),
             panel.spacing = unit(.35, "cm"))

## ---------- c ---------------------------------------------------------------
cc <- cg[outcome %in% c("y_tt","y_shbg","y_ft")]
cc[, OUT := factor(outcome, c("y_tt","y_shbg","y_ft"),
                   c("Total testosterone","SHBG","Free testosterone"))]
cc[, congener := factor(congener, c("PFOS","PFOA","PFHxS","PFNA"))]
setorder(cc, OUT, congener)
cc[, `:=`(lo = b1_std - 1.96*se_b1_std, hi = b1_std + 1.96*se_b1_std)]

pc <- ggplot(cc, aes(b1_std, congener, colour = OUT)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .17, linewidth = .5,
                 position = position_dodge(width = .6)) +
  geom_point(size = 2.1, position = position_dodge(width = .6)) +
  scale_colour_manual(values = c(COL[1], COL[3], COL[2]), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(.10, .10))) +
  labs(x = "Standardised interaction coefficient (per SD)", y = NULL,
       title = "c  Single-congener analysis",
       subtitle = "no attenuation by PFOS or PFOA (highest albumin affinity)") +
  TH + theme(legend.position = "top", legend.direction = "horizontal",
             legend.margin = margin(0, 0, -4, 0), legend.text = element_text(size = 8))

## ---------- d ---------------------------------------------------------------
sdd[, sample := factor(sample, c("all","clean"),
                       c("All men (n = 2,244)", "Drug-affected excluded (n = 2,077)"))]
sdd[, seed_f := factor(seed)]
setorder(sdd, seed_f)

pd_ <- ggplot(sdd, aes(seed_f, p_hetero, colour = sample)) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = "grey35", linewidth = .45) +
  geom_line(aes(group = sample), linewidth = .4, alpha = .55,
            position = position_dodge(width = .34)) +
  geom_point(size = 2.6, position = position_dodge(width = .34)) +
  scale_y_log10(limits = c(0.008, 1.6),
                breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1),
                labels = c("0.01","0.02","0.05","0.1","0.2","0.5","1")) +
  scale_colour_manual(values = c("grey45", COL[2]), name = NULL) +
  labs(x = "Random seed", y = "Heterogeneity P (one-sided, log)",
       title = "d  Fragility of the heterogeneity test",
       subtitle = "P spans 0.012-0.091 across seeds (dashed = 0.05); null after excluding 167 drug users") +
  TH + theme(legend.position = "top", legend.direction = "horizontal",
             legend.margin = margin(0, 0, -2, 0), legend.text = element_text(size = 7.8),
             plot.margin = margin(6, 10, 6, 14))

## ---------- 组装 ------------------------------------------------------------
p <- pa / (pb | pc) / pd_ +
  plot_layout(heights = c(1.22, 1.05, 0.78)) +
  plot_annotation(
    caption = paste(
      "NHANES 2013-2016 (cycles H, I); men aged >= 20 with complete testosterone, SHBG and estradiol.",
      "Models: outcome ~ weight loss x (PFAS + PFAS^2) + covariates; CR1 cluster-robust SE (strata x PSU, 60 clusters).",
      "\nPFAS = molar sum of 6 congeners (within-cycle z-score); single-congener models use within-cycle z-scores of each analyte.",
      "* P<0.05, ** P<0.01, *** P<0.001. Albumin is a negative control for hepatic protein synthesis."),
    theme = theme(plot.caption = element_text(size = 6.9, colour = "grey35", hjust = 0)))

ggsave(file.path(FIG, "fig8_androgen_axis.png"), p, width = 12.6, height = 10.2, dpi = 300)
ggsave(file.path(FIG, "fig8_androgen_axis.pdf"), p, width = 12.6, height = 10.2)
cat("saved fig8\n")
