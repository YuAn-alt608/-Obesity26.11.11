## =====================================================================
## V08  第三层：人群侧置换检验（Placebo Test）
## ---------------------------------------------------------------------
## 统计量只取**回归系数本身**（不经过任何方差估计），因此不受 CR1 聚类协方差
## 数值不稳定的影响：随机化检验下 p = P(|coef_perm| >= |coef_obs|)。
## 两个预设定对比：
##   T3 vs T1  DiD  : y ~ w5*et + 混杂，取 w5:etT3 系数
##   连续暴露       : y ~ w5*mv + 混杂，取 w5:mv 系数
## 暴露在 "周期 x 性别" 层内随机打乱 B 次（保留层结构与边缘分布）。
## 产出：permutation_pop.csv、permutation_pop_null.csv
## =====================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", ".")
OUT  <- file.path(BASE, "verification")
set.seed(20261117)

A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
D <- as.data.frame(as.data.table(readRDS(file.path(BASE, "data/derived/inflam.rds"))))
CONF <- A$CONF

base <- D[D$ok_infl & !is.na(D$w5) & !is.na(D$pfas_z), ]
base$mv <- as.numeric(base$pfas_z)
base$grp <- paste(base$cycle, base$sex)

coef_of <- function(dat, y, kind) {
  dat <- dat[is.finite(dat[[y]]), ]
  dat$et <- cut(dat$mv, unique(quantile(dat$mv, c(0, 1/3, 2/3, 1), na.rm = TRUE)),
                labels = c("T1", "T2", "T3"), include.lowest = TRUE)
  rhs <- if (kind == "tert") "w5*et" else "w5*mv"
  f <- as.formula(paste0(y, " ~ ", rhs, " + ", paste(CONF, collapse = " + ")))
  cf <- coef(lm(f, data = dat))
  if (kind == "tert") cf[["w5:etT3"]] else cf[["w5:mv"]]
}

YK <- c(wbc = "y_wbc", neut = "y_neut", sii = "y_sii", crp = "y_crp")
B <- 1000
obs <- sapply(YK, function(y) c(tert = coef_of(base, y, "tert"), cont = coef_of(base, y, "cont")))

nulls <- array(NA_real_, dim = c(B, length(YK), 2),
               dimnames = list(NULL, names(YK), c("tert", "cont")))
for (b in seq_len(B)) {
  dd <- base
  dd$mv <- ave(base$mv, base$grp, FUN = function(x) sample(x))
  for (k in names(YK)) {
    nulls[b, k, "tert"] <- coef_of(dd, YK[[k]], "tert")
    nulls[b, k, "cont"] <- coef_of(dd, YK[[k]], "cont")
  }
  if (b %% 200 == 0) cat("  置换", b, "/", B, "\n")
}

RES <- rbindlist(lapply(names(YK), function(k) rbindlist(lapply(c("tert", "cont"), function(kind) {
  np <- nulls[, k, kind]
  data.table(outcome = c(wbc = "WBC", neut = "Neutrophil", sii = "SII", crp = "CRP")[[k]],
             contrast = ifelse(kind == "tert", "T3 vs T1 (DiD)", "continuous per 1 SD"),
             coef_obs = obs[kind, k], null_sd = sd(np, na.rm = TRUE),
             p_perm = (1 + sum(abs(np) >= abs(obs[kind, k]), na.rm = TRUE)) / (B + 1),
             B = B)
}))))
fwrite(RES, file.path(OUT, "permutation_pop.csv"))

ND <- as.data.table(as.data.frame.table(nulls, responseName = "coef"))
setnames(ND, c("perm", "outcome", "contrast", "coef"))
fwrite(ND, file.path(OUT, "permutation_pop_null.csv"))

cat("\n========== V08 随机化检验（系数级，B =", B, "）==========\n")
print(RES[, .(outcome, contrast, coef_obs = signif(coef_obs, 3), null_sd = signif(null_sd, 3),
              z = round(coef_obs / null_sd, 2), p_perm = signif(p_perm, 3))], row.names = FALSE)
