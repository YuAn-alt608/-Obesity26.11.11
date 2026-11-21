## ============================================================================
## 35_fig_breakthroughs.R — 三个突破口的图（全英文标注）
##   fig9  肝-脂轴：强主效应 + 被拆穿的交互
##   fig10 高敏亚组：功效不足与方向相反
##   fig11 单细胞脂肪免疫重塑
##   fig12 人群层面炎症验证：髓系被钝化、CRP 不受影响
## ============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork); library(scales)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
RES <- file.path(BASE, "results"); FIG <- file.path(BASE, "figures")

COL <- c("#2C6E8F", "#C4623C", "#6B8E5A", "#8A6FA8", "#B8912F")
TH <- theme_bw(base_size = 9.5) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(colour = "grey92", linewidth = .3),
        plot.tag = element_text(face = "bold", size = 12),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8.1, colour = "grey35"),
        legend.key.size = unit(.36, "cm"),
        strip.background = element_rect(fill = "grey94", colour = NA))
star <- function(p) ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))
fP <- function(p) ifelse(is.na(p), "NA", ifelse(p < .001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p)))
pctCI <- function(est, p) {                      # est 为 log 尺度
  se <- abs(est)/qnorm(1 - pmax(p, 1e-300)/2)
  data.table(est = 100*(exp(est)-1), lo = 100*(exp(est-1.96*se)-1), hi = 100*(exp(est+1.96*se)-1))
}

## ============================================================================
## fig9 — 肝-脂轴
## ============================================================================
per5 <- fread(file.path(RES, "liver_per5.csv"))
dec  <- fread(file.path(RES, "liver_decomp.csv"))
ter  <- fread(file.path(RES, "liver_tertile.csv"))

ORD9 <- c("ALT","AST","GGT","ALP","Bilirubin","Fatty liver index",
          "Hepatic steatosis index","FIB-4","AST:ALT ratio")
pa <- local({
  L <- readRDS(file.path(BASE, "data/derived/liver.rds"))
  L <- L[ok_liver == TRUE & !is.na(w5)]
  x <- per5[match(ORD9, label)]
  ci <- pctCI(x$per5, x$p)
  x <- cbind(x, ci)
  x[, label := factor(label, rev(ORD9))]
  x[, sig := star(p)]
  ggplot(x, aes(est, label)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = .17, colour = "grey40", linewidth = .55) +
    geom_point(shape = 21, size = 2.5, fill = COL[1], colour = "grey20", stroke = .45) +
    geom_text(aes(x = max(hi) + 1.4, label = sig), hjust = 0, size = 3, colour = "grey25") +
    scale_x_continuous(expand = expansion(mult = c(.10, .14))) +
    labs(x = "Effect per 5% weight loss (%)", y = NULL,
         title = "a  Weight loss improves every liver outcome",
         subtitle = sprintf("8 cycles, N = %s (treated %s)", comma(nrow(L)), comma(sum(L$w5)))) +
    TH
})

pb <- local({
  x <- dec[term %in% c("hsi","hsi_liver","fli","fli_liver","y_alt2","bmi","dbmi")]
  x[, grp := fifelse(term %in% c("hsi","fli"), "Score incl. BMI",
              fifelse(term %in% c("hsi_liver","fli_liver","y_alt2"), "Liver-only part", "Body size"))]
  x[, lab := c(hsi = "HSI", hsi_liver = "HSI, liver part", fli = "FLI",
               fli_liver = "FLI, liver part", y_alt2 = "log ALT",
               bmi = "Current BMI", dbmi = "BMI change")[term]]
  setorder(x, -p_joint)                       # 最显著者排最后 -> 显示在最上方
  x[, lab := factor(lab, lab)]
  ggplot(x, aes(-log10(p_joint), lab, fill = grp)) +
    geom_vline(xintercept = -log10(.05), linetype = 2, colour = "grey45", linewidth = .35) +
    geom_col(width = .6, colour = "grey25", linewidth = .25) +
    geom_text(aes(label = paste0("P = ", fP(p_joint))), hjust = -0.10, size = 2.6, colour = "grey25") +
    scale_fill_manual(values = c(`Score incl. BMI` = COL[2], `Liver-only part` = COL[1],
                                 `Body size` = "grey55"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0, .46))) +
    labs(x = expression(-log[10]~"interaction P (weight loss x PFAS)"), y = NULL,
         title = "b  The interaction is a BMI artefact",
         subtitle = "HSI / FLI lose all PFAS interaction once their BMI term is removed") +
    TH + theme(legend.position = "bottom", legend.margin = margin(0,0,0,0))
})

pc <- local({
  x <- ter[outcome %in% c("y_alt","y_hsi","y_fli")]
  x[, lab := c(y_alt = "ALT", y_hsi = "HSI", y_fli = "FLI")[outcome]]
  m <- melt(x[, .(lab, p_tert, T1 = T1p, T2 = T2p, T3 = T3p)], id.vars = c("lab","p_tert"),
            variable.name = "tert", value.name = "pct")
  m[, tert := factor(tert, c("T1","T2","T3"), c("Low","Mid","High"))]
  m[, lab := factor(lab, c("ALT","HSI","FLI"))]
  m[, ytxt := pct - sign(pct)*max(abs(m$pct))*.10]
  ORDL <- c("ALT","HSI","FLI")
  m[, labP := factor(paste0(lab, "\nP = ", fP(p_tert)),
                     paste0(ORDL, "\nP = ", fP(p_tert[match(ORDL, lab)])))]
  ggplot(m, aes(tert, pct, fill = tert)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
    geom_col(width = .7, colour = "grey25", linewidth = .25) +
    geom_text(aes(y = ytxt, label = sprintf("%+.1f", pct)), size = 2.5, colour = "white") +
    facet_wrap(~ labP, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = COL[1:3], guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(.10, .10))) +
    labs(x = "PFAS tertile (L/M/H)", y = "Effect of 5% weight loss (%)",
         title = "c  No monotonic exposure-response",
         subtitle = "Tertile-specific effects are flat or non-monotonic") +
    TH
})

pd <- local({
  L <- readRDS(file.path(BASE, "data/derived/liver.rds"))
  L <- L[ok_liver == TRUE & w5 == 1 & !is.na(pfas_z)]
  L[, mv := as.numeric(pfas_z)]
  L[, pfas_t := cut(mv, unique(quantile(mv, c(0,1/3,2/3,1), na.rm = TRUE)),
                    labels = c("Low","Mid","High"), include.lowest = TRUE)]
  L[, dbmi := bmi - bmi_p]
  s <- L[, .(med = median(dw1, na.rm = TRUE), ge10 = 100*mean(dw1 >= 10, na.rm = TRUE),
             dbmi = mean(dbmi, na.rm = TRUE), N = .N), by = pfas_t]
  s[, pfas_t := factor(pfas_t, c("Low","Mid","High"))]
  print(s)
  yl <- max(s$med)*1.26
  s[, xt := factor(sprintf("%s\n(%.0f%% lost >=10%%)", pfas_t, ge10),
                   sprintf("%s\n(%.0f%% lost >=10%%)", c("Low","Mid","High"),
                           ge10[match(c("Low","Mid","High"), as.character(pfas_t))]))]
  ggplot(s, aes(xt, med)) +
    geom_col(width = .58, fill = COL[4], colour = "grey25", linewidth = .25) +
    geom_text(aes(label = sprintf("%.1f%%", med)), vjust = -.75, size = 3.0) +
    scale_y_continuous(limits = c(0, yl), expand = expansion(mult = c(0, 0))) +
    labs(x = "PFAS tertile (treated group only)", y = "Median 1-year loss (%)",
         title = "d  Why the interaction looked real",
         subtitle = "High-PFAS treated subjects lost less weight to begin with") +
    TH
})

f9 <- (pa | pb) / (pc | pd) + plot_layout(heights = c(1.05, 1))
ggsave(file.path(FIG, "fig9_liver_axis.png"), f9, width = 9.6, height = 7.2, dpi = 330)
ggsave(file.path(FIG, "fig9_liver_axis.pdf"), f9, width = 9.6, height = 7.2)

## ============================================================================
## fig10 — 高敏亚组：功效不足
## ============================================================================
S  <- fread(file.path(RES, "subgroup_highsens.csv"))
PW <- fread(file.path(RES, "subgroup_power.csv"))
TAGMAP <- c("全部（对照）" = "All men\n(N=247)", "20-50岁" = "Age 20-50\n(N=140)",
            "20-50岁 且 减重>=10%" = "Age 20-50 &\n>=10% loss (N=56)",
            "20-50岁 且 减重>=10%（含用药）" = "Same, incl.\nmedicated (N=57)")
LV <- unname(TAGMAP)
S[,  tag := factor(unname(TAGMAP[tag]), LV)]
PW[, tag := factor(unname(TAGMAP[tag]), LV)]

pa10 <- local({
  x <- S[marker == "PFOS" & outcome %in% c("y_tt","y_shbg","y_ft","y_pctfree")]
  x[, out := factor(c(y_tt = "Total T", y_shbg = "SHBG", y_ft = "Free T",
                      y_pctfree = "%Free T")[outcome], c("Total T","SHBG","Free T","%Free T"))]
  x[, sig := star(p)]
  ggplot(x, aes(est, tag, colour = out)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = .18, linewidth = .5, colour = "grey55") +
    geom_point(size = 1.9) +
    facet_wrap(~ out, nrow = 1) +
    scale_colour_manual(values = COL[1:4], guide = "none") +
    labs(x = "Change in outcome per +1 SD PFOS, treated only (%)", y = NULL,
         title = "a  The high-sensitivity subgroup points the wrong way",
         subtitle = "Total T rises, not falls, in men 20-50 losing >=10%") +
    TH
})

pb10 <- local({
  x <- merge(PW, S[marker == "PFOS" & outcome == "y_tt", .(tag, est, p)], by = "tag")
  x[, sig := fifelse(p < .05, "P < 0.05", "n.s.")]
  ggplot(x, aes(est, tag)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
    geom_errorbarh(aes(xmin = -mde, xmax = mde), height = .22,
                   colour = "grey78", linewidth = 1.6) +
    geom_point(shape = 21, size = 2.6, fill = COL[1], colour = "grey25", stroke = .5) +
    geom_text(aes(x = mde + 1.5, label = sprintf("MDE %+.0f%%", mde)), hjust = 0,
              size = 2.5, colour = "grey45") +
    scale_x_continuous(expand = expansion(mult = c(.10, .24))) +
    labs(x = "Observed (dot) vs minimum detectable (band) effect on total T (%)", y = NULL,
         title = "b  Power, not biology",
         subtitle = "No estimate reaches the 80%-power threshold of its own subgroup") +
    TH + theme(axis.text.y = element_text(size = 7.6))
})

pc10 <- local({
  x <- S[outcome %in% c("y_tt","y_ft") & grepl(">=10%$", as.character(tag)) |
         outcome %in% c("y_tt","y_ft") & grepl("medicated", as.character(tag))]
  x <- S[outcome %in% c("y_tt","y_ft")][grepl(">=10%", tag)]
  x[, cong := marker]
  x[, out := c(y_tt = "Total T", y_ft = "Free T")[outcome]]
  ggplot(x, aes(marker, est, fill = marker)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
    geom_col(width = .55, colour = "grey25", linewidth = .25, alpha = .9) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = .16, linewidth = .45) +
    geom_text(aes(label = paste0(est = round(est, 0), "%\n", star(p))),
              vjust = -.35, size = 2.6, colour = "grey20") +
    facet_wrap(~ out) +
    scale_fill_manual(values = c(PFOS = COL[1], PFOA = COL[2]), guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(.12, .30))) +
    labs(x = "Congener", y = "Effect per +1 SD (%)",
         title = "c  Both congeners, same direction",
         subtitle = "PFOS and PFOA agree: no androgen suppression") +
    TH
})

f10 <- (pa10 / pb10) | pc10
f10 <- f10 + plot_layout(widths = c(1.55, 1))
ggsave(file.path(FIG, "fig10_subgroup_power.png"), f10, width = 9.4, height = 6.0, dpi = 330)
ggsave(file.path(FIG, "fig10_subgroup_power.pdf"), f10, width = 9.4, height = 6.0)

## ============================================================================
## fig11 — 单细胞脂肪免疫重塑
## ============================================================================
cmp <- fread(file.path(RES, "geo_composition.csv"))
mac <- fread(file.path(RES, "geo_macrophage_modules.csv"))
pol <- fread(file.path(RES, "sc_m1m2.csv"))
lrp <- fread(file.path(RES, "sc_lr_pairs.csv"))
kg  <- fread(file.path(RES, "geo_key_genes.csv"))
GRP <- c("Lean","Obese","WL"); GC <- c(Lean = "grey55", Obese = COL[2], WL = COL[1])

pa11 <- local({
  x <- cmp[celltype == "macrophage"]
  x[, group := factor(group, GRP)]
  pv <- t.test(x[group=="WL"]$pct, x[group=="Obese"]$pct)$p.value
  s <- x[, .(m = mean(pct), se = sd(pct)/sqrt(.N)), by = group]
  ggplot(x, aes(group, pct, fill = group)) +
    geom_col(data = s, aes(group, m), width = .48, alpha = .35, colour = NA) +
    geom_errorbar(data = s, aes(group, y = m, ymin = m - se, ymax = m + se), width = .16, linewidth = .45) +
    geom_jitter(width = .09, size = 1.7, shape = 21, colour = "grey20", stroke = .4) +
    geom_text(data = s, aes(group, y = 38.5, label = sprintf("%.1f%%", m)), size = 2.7, colour = "grey20") +
    scale_fill_manual(values = GC, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, .08))) +
    labs(x = NULL, y = "Macrophage share of cells (%)",
         title = "a  Macrophages retreat after weight loss",
         subtitle = sprintf("25.3%% -> 13.5%% (P = %.3f); Lean = 14.4%%", pv)) +
    TH
})

pb11 <- local({
  z <- merge(cmp[celltype == "macrophage", .(gsm, group, frac = pct)],
             mac[, .(gsm, percell_z = inflamma)], by = "gsm")
  z <- rbind(z[, .(group, metric = "Per-cell\ninflammation (z)", value = percell_z)],
             pol[, .(group, metric = "Pro-inflammatory\noutput (log10CPM)", value = inflam_output)])
  z[, group := factor(group, GRP)]
  s <- z[, .(m = mean(value), se = sd(value)/sqrt(.N)), by = .(metric, group)]
  s[, metric := factor(metric, c("Per-cell\ninflammation (z)", "Pro-inflammatory\noutput (log10CPM)"))]
  ggplot(s, aes(group, m, fill = group)) +
    geom_col(width = .6, colour = "grey25", linewidth = .25) +
    geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .17, linewidth = .45) +
    facet_wrap(~ metric, scales = "free_y") +
    scale_fill_manual(values = GC, guide = "none") +
    geom_text(aes(label = sprintf("%.2f", m), y = m + se), vjust = -.5, size = 2.5) +
    labs(x = NULL, y = NULL,
         title = "b  Fewer cells, but each one angrier",
         subtitle = "Two metrics disagree: per-cell score up, total output down") +
    TH
})

pc11 <- local({
  x <- lrp[order(-abs(pct))]
  x[, pair := factor(pair, rev(pair))]
  x[, sig := fifelse(p_bh < .05, "BH < 0.05", "n.s.")]
  ggplot(x, aes(pct, pair, fill = sig)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
    geom_col(width = .62, colour = "grey25", linewidth = .22) +
    geom_text(aes(label = sprintf("%+.0f%%", pct), x = pct + sign(pct)*2.5,
                  hjust = ifelse(pct < 0, 1.08, -0.08)), size = 2.4, colour = "grey25") +
    scale_fill_manual(values = c(`BH < 0.05` = COL[1], `n.s.` = "grey72"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(.16, .16))) +
    labs(x = "Change in ligand-receptor activity, WL vs Obese (%)", y = NULL,
         title = "c  Which communication routes close",
         subtitle = "IL6 / IL1B / SPP1 routes fall; CSF1-CSF1R is the one that rises") +
    TH + theme(legend.position = c(.80, .16))
})

pd11 <- local({
  G <- c("IL6","IL1B","CXCL8","NLRP3","CCL2","CASP1","PYCARD")
  x <- kg[gene %in% G]
  x[, grp := factor(group, GRP)]
  s <- dcast(x, gene ~ group, value.var = "mean_cpm")
  s[, pct := 100*(WL/Obese - 1)]
  s[, gene := factor(gene, s[order(pct)]$gene)]
  ggplot(s, aes(pct, gene, fill = pct > 0)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
    geom_col(width = .62, colour = "grey25", linewidth = .22) +
    geom_text(aes(label = sprintf("%+.0f%%", pct), x = pct + sign(pct)*5,
                  hjust = ifelse(pct < 0, 1.08, -0.08)), size = 2.5, colour = "grey25") +
    scale_fill_manual(values = c(`TRUE` = COL[2], `FALSE` = COL[1]), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(.18, .18))) +
    labs(x = "Change in adipocyte-free (stromal) expression, WL vs Obese (%)", y = NULL,
         title = "d  Effector genes split in two directions",
         subtitle = "Inflammasome adaptors (CASP1, PYCARD) and CCL2 up, effectors down") +
    TH
})

f11 <- (pa11 | pb11) / (pc11 | pd11) + plot_layout(heights = c(1, 1.15))
ggsave(file.path(FIG, "fig11_scrna_remodel.png"), f11, width = 9.4, height = 7.2, dpi = 330)
ggsave(file.path(FIG, "fig11_scrna_remodel.pdf"), f11, width = 9.4, height = 7.2)

## ============================================================================
## fig12 — 人群层面炎症验证
## ============================================================================
IP <- fread(file.path(RES, "inflam_per5.csv"))
IT <- fread(file.path(RES, "inflam_tertile.csv"))
ID <- fread(file.path(RES, "inflam_dbmi_interaction.csv"))

ORDI <- c("CRP","WBC","Neutrophil","Lymphocyte","Monocyte count","Monocyte %",
          "NLR","MLR","PLR","SII")
pa12 <- local({
  x <- IP[match(ORDI, label)]
  x[, label := factor(label, rev(ORDI))]
  x[, sig := star(p)]
  x[, myeloid := label %in% c("WBC","Neutrophil","Monocyte count","Monocyte %","SII","NLR","MLR")]
  ggplot(x, aes(per5, label)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = .17, colour = "grey40", linewidth = .55) +
    geom_point(aes(fill = myeloid), shape = 21, size = 2.5, colour = "grey20", stroke = .45) +
    geom_text(aes(x = hi + .35, label = sig), hjust = 0, size = 2.9, colour = "grey25") +
    scale_fill_manual(values = c(`TRUE` = COL[2], `FALSE` = COL[1]),
                      labels = c(`TRUE` = "Myeloid / cell count", `FALSE` = "Acute-phase / ratio"),
                      name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(.10, .16))) +
    labs(x = "Effect per 5% weight loss (%)", y = NULL,
         title = "a  Weight loss lowers systemic inflammation",
         subtitle = sprintf("N = %s (CRP n = %s); every outcome except monocyte count moves",
                            comma(max(IP$N)), comma(IP[label=="CRP"]$N))) +
    TH + theme(legend.position = "top")
})

pb12 <- local({
  ORDB <- c("WBC","Neutrophil","SII","CRP")
  x <- IT[match(ORDB, label)]
  m <- melt(x[, .(label, p_int, T1, T2, T3)], id.vars = c("label","p_int"))
  m[, tert := factor(variable, c("T1","T2","T3"), c("L","M","H"))]
  m[, label := factor(label, ORDB)]
  an <- m[, .(p_int = p_int[1], y = max(abs(value))*1.05), by = label]
  ggplot(m, aes(tert, value, fill = tert)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
    geom_col(width = .68, colour = "grey25", linewidth = .25) +
    geom_text(data = an, aes(x = 2, y = y, label = paste0("P = ", fP(p_int))),
              inherit.aes = FALSE, size = 2.6, colour = "grey25") +
    facet_wrap(~ label, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = COL[c(1,4,2)], guide = "none") +
    labs(x = "PFAS tertile (L/M/H)", y = "Effect of 5% weight loss (%)",
         title = "b  High PFAS blunts the cell-count response, not CRP",
         subtitle = "WBC / neutrophils / SII keep the gradient; CRP does not") +
    TH
})

pc12 <- local({
  ORDC <- c("WBC","Neutrophil","SII","Monocyte %","CRP")
  m <- melt(ID[, .(label, p_int, T1, T2, T3)], id.vars = c("label","p_int"))
  m[, tert := factor(variable, c("T1","T2","T3"), c("L","M","H"))]
  m[, label := factor(label, ORDC)]
  an <- m[, .(p_int = p_int[1], y = max(abs(value))*1.05), by = label]
  ggplot(m, aes(tert, value, fill = tert)) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey45", linewidth = .35) +
    geom_col(width = .68, colour = "grey25", linewidth = .25) +
    geom_text(data = an, aes(x = 2, y = y, label = paste0("P = ", fP(p_int))),
              inherit.aes = FALSE, size = 2.6, colour = "grey25") +
    facet_wrap(~ label, nrow = 1, scales = "free_y") +
    scale_fill_manual(values = COL[c(1,4,2)], guide = "none") +
    labs(x = "PFAS tertile (L/M/H)", y = NULL,
         title = "c  Same test with measured BMI change",
         subtitle = "Per -1 kg/m2; the myeloid gradient survives") +
    TH
})

f12 <- pa12 / (pb12 | pc12) + plot_layout(heights = c(1.15, 1))
f12 <- f12 + plot_annotation(caption = "Myeloid outcomes are the ones PFAS blunts; the acute-phase protein CRP is not.")
ggsave(file.path(FIG, "fig12_inflammation.png"), f12, width = 9.6, height = 6.6, dpi = 330)
ggsave(file.path(FIG, "fig12_inflammation.pdf"), f12, width = 9.6, height = 6.6)

cat("✅ 图 9-12 已输出\n")
