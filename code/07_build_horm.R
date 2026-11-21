## ============================================================================
## 07_build_horm.R — 性激素（总睾酮 / SHBG / 雌二醇 / 游离睾酮）
##   Vermeulen 平衡解离模型计算游离睾酮
## ============================================================================
suppressMessages({library(haven); library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", "."); RAW <- file.path(BASE, "data/raw")
rd <- function(s) { f <- file.path(RAW, paste0(s, ".xpt")); if (!file.exists(f)) return(NULL)
                    as.data.table(suppressWarnings(read_xpt(f))) }
gv <- function(d, cands) { for (v in cands) if (!is.null(d) && v %in% names(d)) return(as.numeric(d[[v]]))
                           rep(NA_real_, nrow(d)) }

## Vermeulen (1999): FT 解二次方程
free_T <- function(tt_ngdl, shbg_nmol, alb_gL) {
  ok <- !is.na(tt_ngdl) & !is.na(shbg_nmol) & !is.na(alb_gL)
  out <- rep(NA_real_, length(tt_ngdl))
  if (!any(ok)) return(out)
  Kt <- 1e9; Ka <- 3.6e4
  Tt <- tt_ngdl[ok] * 0.03467 * 1e-9        # ng/dL -> mol/L  (1 ng/dL = 0.03467 nmol/L)
  S  <- shbg_nmol[ok] * 1e-9                # nmol/L -> mol/L
  Ca <- alb_gL[ok] / 66430                  # g/L -> mol/L
  N  <- 1 + Ka * Ca
  a  <- N * Kt
  b  <- N + Kt * S - Tt * Kt
  FT <- (-b + sqrt(b^2 + 4 * a * Tt)) / (2 * a)   # mol/L
  out[ok] <- FT * 1e9                        # nmol/L
  out
}

H <- list()
## 2003-2004 surplus sex steroid (men)
d <- rd("SSCHL_C")
if (!is.null(d)) H[["C"]] <- data.table(cycle = "C", SEQN = as.numeric(d$SEQN),
  tt_ngdl = gv(d, "SSTESTO"), shbg = gv(d, "SSSHBG"), e2_pgml = gv(d, "SSSE2"), wt_h = NA_real_)
## 2011-2012 total T only
d <- rd("TST_G")
if (!is.null(d)) H[["G"]] <- data.table(cycle = "G", SEQN = as.numeric(d$SEQN),
  tt_ngdl = gv(d, "LBXTST"), shbg = NA_real_, e2_pgml = NA_real_, wt_h = NA_real_)
## 2013-2014 / 2015-2016: T + SHBG + E2
for (s in c("H", "I")) {
  d <- rd(paste0("TST_", s))
  if (!is.null(d)) H[[s]] <- data.table(cycle = s, SEQN = as.numeric(d$SEQN),
    tt_ngdl = gv(d, "LBXTST"), shbg = gv(d, "LBXSHBG"), e2_pgml = gv(d, "LBXEST"), wt_h = NA_real_)
}
## 2017-2018: SHBG + E2 (无总睾酮)
d <- rd("SSTST_J")
if (!is.null(d)) H[["J"]] <- data.table(cycle = "J", SEQN = as.numeric(d$SEQN),
  tt_ngdl = NA_real_, shbg = gv(d, "SSTSHBG"), e2_pgml = gv(d, "SSTEST"), wt_h = NA_real_)

HORM <- rbindlist(H, fill = TRUE)
CORE <- readRDS(file.path(BASE, "data/derived/core.rds"))
HORM <- merge(HORM, CORE[, .(SEQN, alb_gL, sex, age)], by = "SEQN", all.x = TRUE)
HORM[, ft_nmol := free_T(tt_ngdl, shbg, alb_gL)]
HORM[, `:=`(ft_pgml = ft_nmol * 288.42, ft_ngdl = ft_nmol * 28.84)]  # 1 nmol/L = 28.84 ng/dL
HORM[, pct_free := ft_ngdl / tt_ngdl * 100]

saveRDS(HORM, file.path(BASE, "data/derived/horm.rds"))

## QC
print(HORM[, .(N = .N, T = sum(!is.na(tt_ngdl)), SHBG = sum(!is.na(shbg)),
               E2 = sum(!is.na(e2_pgml)), FT = sum(!is.na(ft_nmol)),
               FT_pct_median = round(median(pct_free, na.rm = TRUE), 2)), by = cycle][order(cycle)])
cat("\n=== 男性总睾酮 (ng/dL) ===\n")
print(summary(HORM[sex == 1]$tt_ngdl))
cat("\n=== 游离睾酮 (ng/dL) ===\n"); print(summary(HORM$ft_ngdl))
cat("\n=== %游离睾酮（应约 2%）===\n"); print(summary(HORM$pct_free))
