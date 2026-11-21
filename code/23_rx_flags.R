## ============================================================================
## 23_rx_flags.R — 从 RXQ_RX 识别影响雄激素轴的药物
##   目的：药物性去势/抑制是睾酮结局的强混杂，必须先flag再决定主分析口径
##   输出：data/derived/rx_flags.rds（SEQN + 6 类 flag）
## ============================================================================
suppressMessages({library(foreign); library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
CY <- c("C","D","F","G","H","I","J","L")
rx <- rbindlist(lapply(CY, function(s) {
  f <- file.path(BASE, "data/raw", paste0("RXQ_RX_", s, ".xpt"))
  if (!file.exists(f)) { cat("missing", f, "\n"); return(NULL) }
  x <- as.data.table(read.xport(f))
  if (!"RXDDRUG" %in% names(x)) { cat("skip", s, "：无 RXDDRUG 列（仅问卷变量）\n"); return(NULL) }
  x <- x[, .(SEQN, drug = toupper(trimws(RXDDRUG)))]
  x <- x[!is.na(drug) & drug != ""]
  x[, cycle := s]
  x
}))
cat("处方记录:", nrow(rx), " 人:", uniqueN(rx$SEQN), "\n")

PAT <- list(
  adt      = "LEUPROLIDE|GOSERELIN|DEGARELIX|TRIPTORELIN|BUSERELIN|HISTRELIN|RELUGOLIX|ABARELIX|GNRH",
  antiand  = "BICALUTAMIDE|ENZALUTAMIDE|APALUTAMIDE|DAROLUTAMIDE|NILUTAMIDE|FLUTAMIDE|CYPROTERONE",
  estrogen = "ESTRADIOL|PREMARIN|CONJUGATED ESTROGEN|ESTROGEN|ETHINYL",
  glucoc   = "PREDNISONE|PREDNISOLONE|DEXAMETHASONE|METHYLPREDNISOLONE|CORTISONE|HYDROCORTISONE|BUDESONIDE|TRIAMCINOLONE",
  opioid   = "OXYCODONE|HYDROCODONE|MORPHINE|METHADONE|TRAMADOL|FENTANYL|BUPRENORPHINE|HYDROMORPHONE|CODEINE|TAPENTADOL",
  ketoc    = "KETOCONAZOLE",
  fiveari  = "FINASTERIDE|DUTASTERIDE"
)
for (k in names(PAT)) rx[, (k) := as.integer(grepl(PAT[[k]], drug))]

FL <- rx[, lapply(.SD, max), by = SEQN, .SDcols = names(PAT)]
FL[, castration := as.integer(adt | antiand | estrogen)]
FL[, suppress   := as.integer(glucoc | opioid | ketoc)]
cat("\n各类药物人数（全程 1999-2018 任意周期）:\n")
print(FL[, lapply(.SD, sum), .SDcols = c(names(PAT), "castration", "suppress")])

saveRDS(FL, file.path(BASE, "data/derived/rx_flags.rds"))

## ---------- 在 S3 / S4 中看影响 -------------------------------------------------
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
for (nm in c("S3","S4")) {
  d <- FL[copy(A[[nm]]), on = "SEQN"]
  d[is.na(castration), castration := 0L]; d[is.na(suppress), suppress := 0L]
  cat("\n===", nm, " N =", nrow(d), "======\n")
  cat("  药物去势/抗雄/雌激素:", sum(d$castration), " (处理组中", sum(d$castration[d$w5==1],na.rm=TRUE), ")\n")
  cat("  糖皮质激素/阿片/酮康唑:", sum(d$suppress), " (处理组中", sum(d$suppress[d$w5==1],na.rm=TRUE), ")\n")
  cat("  睾酮分布（去势与否）:\n")
  print(d[, .(N=.N, tt_med = round(median(tt_ngdl,na.rm=TRUE),0),
              tt_min = round(suppressWarnings(min(tt_ngdl,na.rm=TRUE)),1), p1 = round(quantile(tt_ngdl,.01,na.rm=TRUE),0)),
          by = .(castration)])
  cat("  TT < 50 ng/dL 且非去势:", sum(d$tt_ngdl < 50 & d$castration==0, na.rm=TRUE), "\n")
  cat("  TT < 50 ng/dL 且去势  :", sum(d$tt_ngdl < 50 & d$castration==1, na.rm=TRUE), "\n")
  cat("  TT < 50 明细（非去势）:\n")
  print(d[tt_ngdl < 50 & castration==0, .(SEQN, tt_ngdl, shbg, ft_nmol, w5, pfas_z, glucoc, opioid, fiveari, adt, antiand)])
}

## ---------- 排除前后对主效应的影响 ---------------------------------------------
cat("\n\n========== 排除药物影响个体前后：每 5% 减重效应 ==========\n")
crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X); bread <- solve(crossprod(X))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
CONF <- A$CONF
for (nm in c("S3","S4")) {
  d <- FL[copy(A[[nm]]), on = "SEQN"]
  d[is.na(castration), castration := 0L]; d[is.na(suppress), suppress := 0L]
  for (cut in c("included","excluded_any")) {
    dd <- if (cut == "included") copy(d) else d[castration==0 & suppress==0]
    dd <- dd[!is.na(tt_ngdl) & !is.na(dw1) & !is.na(pfas_z)]
    dd[, `:=`(log_tt = log(tt_ngdl), dw5 = dw1/5, mv = pfas_z)]
    cf <- CONF[sapply(CONF, function(v) length(unique(dd[[v]])) > 1)]
    X <- model.matrix(as.formula(paste("~ dw5*mv + mv +", paste(cf, collapse=" + "), "+ race_f + edu_f + smk")), data = dd)
    m <- lm.fit(X, dd$log_tt); V <- crsV(X, m$residuals, factor(paste(dd$strata, dd$psu)))
    j <- grep("^dw5$", colnames(X)); ji <- grep("^dw5:mv$", colnames(X))
    cat(sprintf("  %s %-14s N=%4d  每5%%: %+.4f (p=%.2g)   交互: %+.4f (p=%.3f)\n", nm, cut,
                nrow(dd), m$coefficients[j],
                2*pnorm(-abs(m$coefficients[j]/sqrt(V[j,j]))),
                m$coefficients[ji], 2*pnorm(-abs(m$coefficients[ji]/sqrt(V[ji,ji])))))
  }
}
cat("\n✅ 完成\n")
