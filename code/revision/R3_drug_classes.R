## ============================================================================
## R3_drug_classes.R — 回应审稿意见 #3（167 人构成未列、逐类贡献未量化）
##   逐类剔除 / 逐类保留 / 全剔除，对 ① 参数交互（2df 联合）② 每 5% 主效应
##   ③ 森林异质性 P（仅最大三类，4000 树） 分别量化
## 输出: results/revision/R3_drug_classes.csv, R3_drug_composition.csv
## ============================================================================
suppressMessages({library(data.table); library(grf)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
A  <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
FL <- readRDS(file.path(BASE, "data/derived/rx_flags.rds"))
OUT <- file.path(BASE, "results/revision"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
CONF <- A$CONF

dropconst <- function(X) { s <- apply(X, 2, sd)
  keep <- (s > 0 & !is.na(s)) | colnames(X) == "(Intercept)"
  X[, keep, drop = FALSE] }
crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X); bread <- solve(crossprod(X))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
wald <- function(b, V) as.numeric(t(b) %*% solve(V) %*% b)

d <- FL[copy(A$S3), on = "SEQN"]
for (v in c("adt","antiand","estrogen","glucoc","opioid","ketoc","fiveari","castration","suppress"))
  if (!v %in% names(d)) d[, (v) := 0L]
for (v in c("adt","antiand","estrogen","glucoc","opioid","ketoc","fiveari","castration","suppress"))
  d[is.na(get(v)), (v) := 0L]
d[, `:=`(log_tt = log(tt_ngdl), mv = as.numeric(pfas_z), dw5 = dw1/5)]
d[, mv2 := mv^2]
d[, anydrug := as.integer(castration == 1 | suppress == 1)]

## ---------- ① 构成表 --------------------------------------------------------
comp <- data.table(
  class = c("Opioid analgesic", "5-alpha reductase inhibitor", "Glucocorticoid",
            "GnRH agonist (ADT)", "Antiandrogen", "Ketoconazole", "Estrogen",
            "Any of the 167 (castration OR suppress)", "Five-ARI (not in the 167)"),
  flag  = c("opioid", "fiveari", "glucoc", "adt", "antiand", "ketoc", "estrogen", "anydrug", "fiveari"),
  n_all = c(d[, sum(opioid)], d[, sum(fiveari)], d[, sum(glucoc)], d[, sum(adt)],
            d[, sum(antiand)], d[, sum(ketoc)], d[, sum(estrogen)], d[, sum(anydrug)],
            d[, sum(fiveari)]),
  n_treated = c(d[w5 == 1, sum(opioid)], d[w5 == 1, sum(fiveari)], d[w5 == 1, sum(glucoc)],
                d[w5 == 1, sum(adt)], d[w5 == 1, sum(antiand)], d[w5 == 1, sum(ketoc)],
                d[w5 == 1, sum(estrogen)], d[w5 == 1, sum(anydrug)], d[w5 == 1, sum(fiveari)])
)
comp[, pct_all := round(100 * n_all / nrow(d), 2)]
print(comp)
fwrite(comp, file.path(OUT, "R3_drug_composition.csv"))

## ---------- ② 逐类剔除的参数模型 --------------------------------------------
fit_par <- function(dd, tag) {
  dd <- dd[!is.na(log_tt) & !is.na(pfas_z) & !is.na(dw1)]
  X <- model.matrix(as.formula(paste("~ dw5*(mv+mv2) + mv + mv2 +", paste(CONF, collapse=" + "),
                                     "+ race_f + edu_f + smk")), data = dd)
  X <- dropconst(X)
  m <- lm.fit(X, dd$log_tt); V <- crsV(X, m$residuals, factor(paste(dd$strata, dd$psu)))
  b <- m$coefficients; iw <- grep("^dw5$", colnames(X))
  i1 <- grep("^dw5:mv$", colnames(X)); i2 <- grep("^dw5:mv2$", colnames(X))
  data.table(tag = tag, N = nrow(dd), treated = sum(dd$w5),
             per5 = 100 * (exp(b[iw]) - 1),
             per5_se = 100 * sqrt(V[iw, iw]),
             p_per5 = 2 * pnorm(-abs(b[iw] / sqrt(V[iw, iw]))),
             p_int = pchisq(wald(b[c(i1, i2)], V[c(i1, i2), c(i1, i2)]), 2, lower.tail = FALSE))
}
R <- list(fit_par(d, "none excluded (N=2,244)"))
for (cl in c("opioid", "fiveari", "glucoc", "adt", "antiand", "ketoc")) {
  R[[length(R)+1]] <- fit_par(d[get(cl) == 0], paste0("exclude ", cl))
}
R[[length(R)+1]] <- fit_par(d[anydrug == 0], "exclude the 167")
R[[length(R)+1]] <- fit_par(d[anydrug == 0 & fiveari == 0], "exclude any of 7 classes")
RES <- rbindlist(R)
print(RES[, .(tag, N, per5 = sprintf("%+.2f%%", per5), p_per5 = format.pval(p_per5, digits=3),
              p_int = format.pval(p_int, digits=3))])
fwrite(RES, file.path(OUT, "R3_drug_classes.csv"))

## ---------- ③ 森林异质性：只对三类大药各剔一次 ------------------------------
mkX <- function(dt) {
  X <- cbind(as.matrix(dt[, ..CONF]), model.matrix(~ race_f + edu_f + smk - 1, data = dt))
  X <- cbind(X, as.matrix(dt[, .(pfas_z)]))
  mc <- grep("^miss_", names(dt), value = TRUE)
  if (length(mc)) X <- cbind(X, as.matrix(dt[, ..mc]))
  storage.mode(X) <- "double"
  s <- apply(X, 2, sd); X[, (s > 0 & !is.na(s)) | colnames(X) == "(Intercept)", drop = FALSE]
}
forest_p <- function(dd, tag) {
  dd <- copy(dd)[!is.na(log_tt) & !is.na(pfas_z) & !is.na(w5)]
  X <- mkX(dd); Y <- dd$log_tt; W <- as.numeric(dd$w5)
  cl <- as.integer(factor(paste(dd$strata, dd$psu)))
  e0 <- as.numeric(regression_forest(X, W, num.trees = 2000, tune.parameters = "all",
                                     seed = 20261117)$predictions)
  k <- e0 >= 0.05 & e0 <= 0.95
  cf <- causal_forest(X[k, ], Y[k], W[k], W.hat = e0[k], num.trees = 4000, honesty = TRUE,
                      tune.parameters = "all", clusters = cl[k], seed = 20261117)
  cal <- test_calibration(cf); p1 <- cal["differential.forest.prediction", "Pr(>t)"]
  blp <- best_linear_projection(cf, A = cbind(mv = dd$mv[k]))
  data.table(tag = tag, N = sum(k), p_1s = p1, p_2s = min(1, 2 * p1),
             p_blp = blp["mv", "Pr(>|t|)"])
}
F <- list(forest_p(d, "none excluded"))
for (cl in c("opioid", "glucoc", "fiveari")) F[[length(F)+1]] <- forest_p(d[get(cl) == 0], paste0("exclude ", cl))
F[[length(F)+1]] <- forest_p(d[anydrug == 0], "exclude the 167")
FR <- rbindlist(F); print(FR)
fwrite(FR, file.path(OUT, "R3_drug_classes_forest.csv"))
cat("✅ R3 完成\n")
