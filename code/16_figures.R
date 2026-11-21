## ============================================================================
## 16_figures.R — NHANES 主图 (出版级, 300 dpi)
## ============================================================================
suppressMessages({library(data.table); library(ggplot2); library(patchwork);
                  library(scales); library(viridis)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
FIG <- file.path(BASE, "figures"); RES <- file.path(BASE, "results")
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

th <- theme_classic(base_size = 10) +
  theme(axis.text = element_text(colour = "grey20"),
        axis.title = element_text(colour = "grey10"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, colour = "grey35"),
        legend.position = "top", legend.title = element_text(size = 9),
        strip.background = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92", linewidth = 0.25))
COL <- c(benefit = "#2166AC", harm = "#B2182B", neutral = "grey45")
COLG <- c(Obese = "#D6604D", WL = "#4393C3", Lean = "#92C5DE")

D <- readRDS(file.path(BASE, "data/derived/analysis.rds"))$D
cf <- readRDS(file.path(RES, "cf_results.rds"))
par <- fread(file.path(RES, "param_interaction.csv"))

## ---------------- Fig 1: 暴露分布 + 研究人群 ----------------
p1a <- ggplot(D[adult == TRUE & !is.na(pfas_molar)],
              aes(factor(cycle, levels = c("C","D","F","G","H","I","J","L")), pfas_molar)) +
  geom_boxplot(outlier.size = .2, fill = "grey85", colour = "grey30", linewidth = .3) +
  scale_y_log10(labels = label_number()) +
  labs(x = NULL, y = expression(Sigma*"PFAS (nmol/L, log scale)"), title = "a  Serum PFAS by cycle") +
  th
p1b <- ggplot(melt(D[adult == TRUE & !is.na(pcb_lipid) | !is.na(ocp_lipid),
                     .(PCB = pcb_lipid, OCP = ocp_lipid)], id.vars = character()),
              aes(variable, value)) +
  geom_boxplot(outlier.size = .2, fill = "grey85", colour = "grey30", linewidth = .3) +
  scale_y_log10() + labs(x = NULL, y = "ng/g lipid (log)", title = "b  Lipophilic POPs, 2003-04") + th
p1c <- ggplot(D[adult == TRUE & !is.na(dw1)], aes(dw1)) +
  geom_histogram(bins = 70, fill = "grey70", colour = NA) +
  geom_vline(xintercept = c(-2, 2, 5, 10), linetype = c(3,3,2,2),
             colour = rep(c("grey45", COL[2]), each = 2), linewidth = .4) +
  coord_cartesian(xlim = c(-20, 25)) +
  labs(x = "1-year weight change (%)", y = "Participants",
       title = "c  Weight-change distribution", subtitle = "positive = weight loss") + th
p1d <- ggplot(D[adult == TRUE & !is.na(phenoage)], aes(age, phenoage)) +
  geom_hex(bins = 40) + scale_fill_viridis_c(option = "magma", trans = "log10") +
  geom_abline(slope = 1, intercept = 0, colour = COL[1], linewidth = .5, linetype = 2) +
  labs(x = "Chronological age (y)", y = "PhenoAge (y)", title = "d  PhenoAge vs age",
       subtitle = expression(italic(r) == 0.95)) + th + theme(legend.position = "none")
f1 <- (p1a | p1b) / (p1c | p1d) + plot_annotation(tag_levels = NULL)
try(ggsave(file.path(FIG, "fig1_exposure.png"), f1, width = 8, height = 6.2, dpi = 300)); cat("f1 ok\n")

## ---------------- Fig 2: 平均效应 (naive vs adjusted) ----------------
naive <- D[adult == TRUE & !is.na(phenoage_accel) & !is.na(w5),
           .(m = mean(phenoage_accel), se = sd(phenoage_accel)/sqrt(.N)), by = w5]
naive[, lab := factor(w5, levels = 0:1, labels = c("No rapid loss", "Rapid loss >=5%"))]
ate <- rbindlist(lapply(names(cf), function(k) {
  data.table(model = k, est = cf[[k]]$ate[1], se = cf[[k]]$ate[2]) }))
lab_map <- c(S1 = "PhenoAge | PFAS (n=9,656)", S2 = "PhenoAge | PCB (1,202)",
             S2b = "PhenoAge | OCP (1,294)", S3 = "Testosterone | PFAS (2,244)",
             S4 = "Free T | PFAS (1,337)")
ate[, model := factor(model, levels = names(lab_map), labels = lab_map[names(lab_map)])]
p2a <- ggplot(naive, aes(lab, m)) +
  geom_col(fill = c("grey65", COL[2]), width = .6) +
  geom_errorbar(aes(ymin = m - 1.96*se, ymax = m + 1.96*se), width = .15, linewidth = .4) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey30", linewidth = .3) +
  labs(x = NULL, y = "PhenoAge acceleration (years)", title = "a  Crude comparison",
       subtitle = "rapid losers appear 1.87 y biologically older") + th
p2b <- ggplot(ate[model != "Free T | PFAS (1,337)"], aes(model, est)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey30", linewidth = .3) +
  geom_pointrange(aes(ymin = est - 1.96*se, ymax = est + 1.96*se),
                  colour = COL[1], size = .4, linewidth = .5) +
  coord_flip() + labs(x = NULL,
                      y = "Adjusted ATE (AIPW, 95% CI)\nPhenoAge: years;  testosterone: log units",
                      title = "b  Doubly-robust ATE",
                      subtitle = "confounding explains the crude difference") + th
f2 <- p2a | p2b
try(ggsave(file.path(FIG, "fig2_ate.png"), f2, width = 8, height = 3.6, dpi = 300)); cat("f2 ok\n")

## ---------------- Fig 3: CATE 曲线 ----------------
curve_plot <- function(k, xlab, ttl, ylim = NULL) {
  cv <- cf[[k]]$curve; bd <- cf[[k]]$bydec
  p <- ggplot() +
    geom_ribbon(data = cv, aes(mod, ymin = boot_lo, ymax = boot_hi), alpha = .16, fill = COL[1]) +
    geom_line(data = cv, aes(mod, cate), colour = COL[1], linewidth = .8) +
    geom_hline(yintercept = 0, linetype = 2, colour = "grey25", linewidth = .35) +
    geom_point(data = bd, aes(mod_med, tau_cf), shape = 21, fill = "white", size = 1.8, stroke = .6) +
    labs(x = xlab, y = "CATE (years)", title = ttl,
         subtitle = sprintf("N=%s; heterogeneity P=%s", format(nrow(cf[[k]]$data), big.mark = ","),
                            format.pval(cf[[k]]$nf_p, digits = 2))) + th
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  p
}
f3 <- curve_plot("S1", expression(Sigma*"PFAS (z-score)"), "a  Biological ageing") |
      curve_plot("S2", expression(Sigma*"PCB (z-score)"), "b  Lipophilic PCB", c(-4, 4)) |
      curve_plot("S2b", expression(Sigma*"OCP (z-score)"), "c  Organochlorine pesticides", c(-4, 4))
try(ggsave(file.path(FIG, "fig3_cate.png"), f3, width = 9, height = 3.2, dpi = 300)); cat("f3 ok\n")

## ---------------- Fig 4: 睾酮 + 剂量反应曲面 ----------------
cv3 <- cf$S3$curve; bd3 <- cf$S3$bydec
p4a <- ggplot() +
  geom_ribbon(data = cv3, aes(mod, ymin = boot_lo, ymax = boot_hi), alpha = .16, fill = COL[1]) +
  geom_line(data = cv3, aes(mod, cate), colour = COL[1], linewidth = .8) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey25", linewidth = .35) +
  geom_point(data = bd3, aes(mod_med, tau_cf), shape = 21, fill = "white", size = 1.8, stroke = .6) +
  labs(x = expression(Sigma*"PFAS (z-score)"), y = expression(Delta*" log testosterone"),
       title = "a  Testosterone response",
       subtitle = {
         e <- cf$S3$ate["estimate"]; s <- cf$S3$ate["std.err"]
         pv <- 2 * pnorm(-abs(e / s))
         sprintf("Forest ATE %+.1f%% (%s); benefit attenuates with PFAS",
                 100 * (exp(e) - 1),
                 if (pv < 0.001) "P<0.001" else sprintf("P=%.3f", pv))
       }) + th
## 分层剂量反应: 按 PFAS 三分位的调整后曲线
suppressMessages(library(mgcv))
A1 <- readRDS(file.path(BASE, "data/derived/analysis.rds"))$S1
CONF <- readRDS(file.path(BASE, "data/derived/analysis.rds"))$CONF
A1[, mt := cut(pfas_z, quantile(pfas_z, c(0, 1/3, 2/3, 1)), include.lowest = TRUE,
               labels = c("PFAS T1 (low)", "PFAS T2", "PFAS T3 (high)"))]
cf1 <- CONF[sapply(CONF, function(v) length(unique(A1[[v]])) > 1)]
mod <- gam(as.formula(paste0("phenoage_accel ~ te(dw1, pfas_z, k = c(5,4)) + ",
                             paste(cf1, collapse = " + "), " + race_f + edu_f + smk")),
           data = A1, method = "REML")
lv <- A1[, .(mv = median(pfas_z)), by = mt]
gr <- CJ(dw = seq(-8, 20, by = .5), mt = levels(A1$mt))
gr <- merge(gr, lv, by = "mt")
base <- A1[rep(1, nrow(gr))]
for (v in cf1) base[[v]] <- mean(A1[[v]])
base[, `:=`(race_f = factor("NHWhite", levels = levels(A1$race_f)),
            edu_f = factor("HS", levels = levels(A1$edu_f)), smk = 0,
            dw1 = gr$dw, pfas_z = gr$mv)]
gr[, fit := as.numeric(predict(mod, newdata = base))]
gr[, delta := fit - gr[dw == 0, .(mt, ref = fit)][gr$mt, on = "mt"]$ref]
p4b <- ggplot(gr, aes(dw, delta, colour = mt)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey35", linewidth = .35) +
  geom_vline(xintercept = 5, linetype = 3, colour = "grey55", linewidth = .3) +
  geom_line(linewidth = .8) +
  scale_colour_manual(values = c("#4393C3", "#F4A582", "#B2182B"), name = NULL) +
  labs(x = "1-year weight change (%)", y = "Adjusted contrast in PhenoAge accel (y)",
       title = "b  Dose-response by PFAS stratum",
       subtitle = "reference: stable weight (0%)") + th
f4 <- p4a | p4b
try(ggsave(file.path(FIG, "fig4_testosterone_surface.png"), f4, width = 8.4, height = 3.4, dpi = 300)); cat("f4 ok\n")

## ---------------- Fig 5: 死亡验证 ----------------
mp <- fread(file.path(RES, "mort_phenoage.csv"))[var != "(Intercept)"]
labmap <- c(phenoage_accel = "PhenoAge acceleration (per y)", age = "Chronological age (per y)",
            female = "Female", race_fNHBlack = "Black vs White", bmi = "BMI (per kg/m2)",
            smk = "Current smoker")
mp <- mp[var %in% names(labmap)]
mp[, term := factor(labmap[var], levels = rev(labmap))]
p5a <- ggplot(mp, aes(RR, term)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey30", linewidth = .3) +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = COL[2], size = .35, linewidth = .5) +
  scale_x_log10(breaks = c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3),
                labels = function(x) sprintf("%.1f", x)) +
  labs(x = "Rate ratio for all-cause death (95% CI)", y = NULL,
       title = "a  PhenoAge predicts mortality",
       subtitle = sprintf("RR %.3f per year; AUC %.3f vs %.3f for age",
                          mp$RR[mp$var == "phenoage_accel"], 0.857, 0.856)) + th
mw <- fread(file.path(RES, "mort_w5_pfas.csv"))
p5b <- ggplot(mw, aes(pfas_z, RR)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey30", linewidth = .3) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .15, fill = COL[2]) +
  geom_line(colour = COL[2], linewidth = .8) + geom_point(size = 1.6, colour = COL[2]) +
  labs(x = expression(Sigma*"PFAS (z-score)"), y = "Rate ratio (95% CI)",
       title = "b  Rapid weight loss and mortality",
       subtitle = "no significant effect modification (P=0.46)") + th
f5 <- p5a | p5b
try(ggsave(file.path(FIG, "fig5_mortality.png"), f5, width = 8, height = 3.3, dpi = 300)); cat("f5 ok\n")

## ---------------- Fig 6: 稳健性 ----------------
bc <- fread(file.path(RES, "robust_bycycle.csv"))
p6a <- ggplot(bc, aes(b2, factor(cycle, levels = rev(bc$cycle)))) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey30", linewidth = .3) +
  geom_point(aes(colour = p2 < 0.05), size = 2) +
  scale_colour_manual(values = c(`TRUE` = COL[2], `FALSE` = "grey45"), guide = "none") +
  labs(x = expression("Quadratic interaction"~beta[2]), y = "Survey cycle",
       title = "a  Cycle-specific interaction") + th
pl <- readRDS(file.path(RES, "robust.rds"))
p6b <- ggplot(data.table(x = pl$placebo[1, ]), aes(x)) +
  geom_histogram(bins = 25, fill = "grey75", colour = NA) +
  geom_vline(xintercept = 0.05, colour = COL[1], linetype = 2, linewidth = .4) +
  labs(x = "Placebo interaction P value", y = "Count",
       title = "b  100 placebo permutations",
       subtitle = sprintf("%.0f%% of permutations reach P<0.05 (expected 5%%)",
                          100 * mean(pl$placebo[1, ] < 0.05))) + th
f6 <- p6a | p6b
try(ggsave(file.path(FIG, "fig6_robustness.png"), f6, width = 8, height = 3.2, dpi = 300)); cat("f6 ok\n")

cat("figures written:\n"); print(list.files(FIG, pattern = "\\.png$"))
