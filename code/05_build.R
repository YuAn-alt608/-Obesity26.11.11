## ============================================================================
## 05_build.R — NHANES 多周期标准化分析数据集（按 SEQN 合并）
## 输出: data/derived/core.rds, results/qc_build_core.csv
## ============================================================================
suppressMessages({library(haven); library(data.table)})
setDTthreads(4)
BASE <- Sys.getenv("AGE22_ROOT", ".")
RAW  <- file.path(BASE, "data/raw")
dir.create(file.path(BASE, "data/derived"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(BASE, "results"), showWarnings = FALSE, recursive = TRUE)

rd <- function(stem) {
  if (length(stem) == 0 || is.na(stem)) return(NULL)
  f <- file.path(RAW, paste0(stem, ".xpt"))
  if (!file.exists(f)) return(NULL)
  as.data.table(suppressWarnings(read_xpt(f)))
}
## 从域文件中按映射抽取变量并统一改名; 映射: 新名 = c(候选原名...)
sel <- function(d, mapping) {
  if (is.null(d) || !"SEQN" %in% names(d)) return(NULL)
  cols <- list(SEQN = as.numeric(d$SEQN))
  for (nm in names(mapping)) {
    v <- intersect(mapping[[nm]], names(d))
    cols[[nm]] <- if (length(v)) as.numeric(d[[v[1]]]) else rep(NA_real_, nrow(d))
  }
  as.data.table(cols)
}

CY <- data.table(suf = c("C","D","F","G","H","I","J","L"),
                 yr  = c(2003,2005,2009,2011,2013,2015,2017,2021))
ST <- list(
  bio = c(C="L40_C",   D="BIOPRO_D", F="BIOPRO_F", G="BIOPRO_G", H="BIOPRO_H", I="BIOPRO_I", J="BIOPRO_J", L="BIOPRO_L"),
  cbc = c(C="L25_C",   D="CBC_D",    F="CBC_F",    G="CBC_G",    H="CBC_H",    I="CBC_I",    J="CBC_J",    L="CBC_L"),
  crp = c(C="L11_C",   D="CRP_D",    F="CRP_F",    I="HSCRP_I",  J="HSCRP_J",  L="HSCRP_L"),
  cot = c(C="L06COT_C",D="COT_D",    F="COTNAL_F", G="COTNAL_G", H="COT_H",    I="COT_I",    J="COT_J",    L="COT_L"),
  tc  = c(C="L13_C",   D="TCHOL_D",  F="TCHOL_F",  G="TCHOL_G",  H="TCHOL_H",  I="TCHOL_I",  J="TCHOL_J",  L="TCHOL_L"),
  hdl = c(C="L13_C",   D="HDL_D",    F="HDL_F",    G="HDL_G",    H="HDL_H",    I="HDL_I",    J="HDL_J",    L="HDL_L"),
  tg  = c(C="L13AM_C", D="TRIGLY_D", F="TRIGLY_F", G="TRIGLY_G", H="TRIGLY_H", I="TRIGLY_I", J="TRIGLY_J", L="TRIGLY_L"),
  glu = c(C="L10AM_C", D="GLU_D",    F="GLU_F",    G="GLU_G",    H="GLU_H",    I="GLU_I",    J="GLU_J",    L="GLU_L"),
  ghb = c(C="L10_C",   D="GHB_D",    F="GHB_F",    G="GHB_G",    H="GHB_H",    I="GHB_I",    J="GHB_J",    L="GHB_L")
)

out <- list()
for (i in seq_len(nrow(CY))) {
  s <- CY$suf[i]
  demo <- rd(paste0("DEMO_", s))
  if (is.null(demo)) { cat("no DEMO for", s, "\n"); next }
  parts <- list(
    sel(demo, list(srvyr="SDDSRVYR", age="RIDAGEYR", sex="RIAGENDR", race="RIDRETH1",
                   edu="DMDEDUC2", pir="INDFMPIR", preg="RIDEXPRG", wt_mec="WTMEC2YR",
                   wt_int="WTINT2YR", psu="SDMVPSU", strata="SDMVSTRA")),
    sel(rd(paste0("BMX_", s)), list(wt_meas="BMXWT", ht_meas="BMXHT", bmi="BMXBMI", waist="BMXWAIST")),
    sel(rd(paste0("WHQ_", s)), list(sr_ht_in="WHD010", sr_wt_now="WHD020", sr_wt_1y="WHD050",
                                    sr_wt_10y="WHD110", sr_wt_max="WHD140", wt_intent="WHQ060",
                                    tried_lose="WHQ070")),
    sel(rd(ST$bio[s]), list(alb="LBXSAL", creat="LBXSCR", glu_bio="LBXSGL", alp="LBXSAPSI")),
    sel(rd(ST$cbc[s]), list(wbc="LBXWBCSI", lym_pct="LBXLYPCT", mcv="LBXMCVSI", rdw="LBXRDW")),
    sel(rd(ST$cot[s]), list(cotinine="LBXCOT")),
    sel(rd(ST$tc[s]),  list(tc="LBXTC")),
    sel(rd(ST$hdl[s]), list(hdl="LBDHDD")),
    sel(rd(ST$tg[s]),  list(tg="LBXTR")),
    sel(rd(ST$glu[s]), list(glu_fast="LBXGLU")),
    sel(rd(ST$ghb[s]), list(hba1c="LBXGH")),
    sel(rd(paste0("DIQ_", s)), list(diab="DIQ010")),
    sel(rd(paste0("BPQ_", s)), list(htn="BPQ020")),
    sel(rd(paste0("SMQ_", s)), list(smoke100="SMQ020", smoke_now="SMQ040")),
    sel(rd(paste0("ALQ_", s)), list(drinks="ALQ130")),
    sel(rd(paste0("PAQ_", s)), list(pa_vig=c("PAQ650","PAQ180"), pa_mod=c("PAQ665","PAQ100")))
  )
  ## CRP: 2015+ 为 hs-CRP (mg/L), 之前为 CRP (mg/dL)
  crp <- rd(ST$crp[s])
  if (!is.null(crp)) {
    if ("LBXHSCRP" %in% names(crp))
      parts <- c(parts, list(sel(crp, list(crp_mgl = "LBXHSCRP"))))
    else
      parts <- c(parts, list({z <- sel(crp, list(v = "LBXCRP")); z[, crp_mgl := v * 10][, v := NULL]; z}))
  }
  parts <- Filter(Negate(is.null), parts)
  D <- Reduce(function(a, b) merge(a, b, by = "SEQN", all.x = TRUE), parts)
  D[, cycle := s]
  out[[s]] <- D
  cat("cycle", s, "N =", nrow(D), "\n")
}
CORE <- rbindlist(out, fill = TRUE)
cat("pooled N =", nrow(CORE), "\n")

## ---------- 单位换算 ----------
CORE[, `:=`(
  alb_gL    = alb * 10,
  creat_um  = creat * 88.4,
  glu_mmol  = glu_bio * 0.0555,
  crp_ok    = ifelse(crp_mgl > 0, crp_mgl, NA_real_),
  sr_kg_now = sr_wt_now * 0.453592,
  sr_kg_1y  = sr_wt_1y  * 0.453592,
  sr_kg_10y = sr_wt_10y * 0.453592,
  ht_sr_m   = sr_ht_in * 0.0254
)]

## ---------- PhenoAge (Levine 2018) ----------
phenoage <- function(alb_gL, creat_um, glu_mmol, crp_mgl, lym_pct, mcv, rdw, alp, wbc, age) {
  xb <- -19.9067 - 0.0336 * alb_gL + 0.0095 * creat_um + 0.1953 * glu_mmol +
         0.0954 * log(pmax(crp_mgl, 0.01)) - 0.0120 * lym_pct + 0.0268 * mcv +
         0.3306 * rdw + 0.00188 * alp + 0.0554 * wbc + 0.0804 * age
  gam <- 0.0076927
  m <- 1 - exp(-exp(xb) * (exp(gam * 120) - 1) / gam)
  m <- pmin(pmax(m, 1e-12), 1 - 1e-12)
  141.50225 + log(-0.00553 * log(1 - m)) / 0.090165
}
CORE[, phenoage := phenoage(alb_gL, creat_um, glu_mmol, crp_ok, lym_pct, mcv, rdw, alp, wbc, age)]
CORE[phenoage < 0 | phenoage > 120, phenoage := NA_real_]

## ---------- 体重变化 ----------
CORE[, `:=`(
  ok_now = sr_kg_now >= 32 & sr_kg_now <= 227,
  ok_1y  = sr_kg_1y  >= 32 & sr_kg_1y  <= 227,
  ok_10y = sr_kg_10y >= 32 & sr_kg_10y <= 227
)]
CORE[, dw1  := fifelse(ok_now & ok_1y,  (sr_kg_1y  - sr_kg_now) / sr_kg_1y  * 100, NA_real_)]
CORE[, dw10 := fifelse(ok_now & ok_10y, (sr_kg_10y - sr_kg_now) / sr_kg_10y * 100, NA_real_)]
CORE[abs(dw1) > 50,  dw1  := NA_real_]
CORE[abs(dw10) > 60, dw10 := NA_real_]
CORE[, bmi_prior := fifelse(ok_1y & ht_sr_m > 1.3 & ht_sr_m < 2.1, sr_kg_1y / ht_sr_m^2, NA_real_)]

## ---------- 协变量 ----------
CORE[, `:=`(
  female = as.integer(sex == 2),
  age2   = age^2,
  race_f = factor(race, levels = 1:5, labels = c("MexAm","OthHisp","NHWhite","NHBlack","Other")),
  edu_f  = factor(edu,   levels = 1:5, labels = c("<9y","9-11y","HS","SomeColl","Coll+")),
  smk    = fifelse(is.na(smoke100) | smoke100 == 2, 0L, fifelse(smoke_now %in% c(1,2), 2L, 1L))
)]
CORE[, egfr := {
  k <- ifelse(female == 1, 0.7, 0.9); a <- ifelse(female == 1, -0.241, -0.302)
  142 * pmin(creat / k, 1)^a * pmax(creat / k, 1)^-1.200 *
    ifelse(female == 1, 1.012, 1) * 0.9938^age
}]
CORE[, pa_any := fifelse(!is.na(pa_vig) | !is.na(pa_mod),
                         as.integer(pa_vig == 1 | pa_mod == 1), NA_integer_)]
CORE[, adult := age >= 20 & (is.na(preg) | preg != 1)]

saveRDS(CORE, file.path(BASE, "data/derived/core.rds"))

## ---------- QC ----------
qc <- CORE[adult == TRUE, .(N = .N,
  wt_history = sum(!is.na(dw1)), loss5 = sum(dw1 >= 5, na.rm = TRUE),
  phenoage = sum(!is.na(phenoage)), alb = sum(!is.na(alb)),
  crp = sum(!is.na(crp_mgl)), wbc = sum(!is.na(wbc)),
  lym = sum(!is.na(lym_pct)), cotinine = sum(!is.na(cotinine)),
  bmi = sum(!is.na(bmi))), by = cycle][order(cycle)]
print(qc)
fwrite(qc, file.path(BASE, "results/qc_build_core.csv"))

cat("\n=== PhenoAge 效度 ===\n")
cat("r(PhenoAge, age) =", round(cor(CORE$phenoage, CORE$age, use = "complete.obs"), 4), "\n")
cat("mean(PhenoAge - age) =", round(mean(CORE$phenoage - CORE$age, na.rm = TRUE), 2), "\n")
print(summary(CORE[adult == TRUE]$phenoage))
cat("\n=== dw1 (%) 分布 ===\n"); print(summary(CORE[adult == TRUE]$dw1))
cat("\n=== 变量缺失率(成人) ===\n")
A <- CORE[adult == TRUE]
vars <- c("age","bmi","bmi_prior","dw1","phenoage","alb","creat","glu_bio","alp","crp_mgl",
          "wbc","lym_pct","mcv","rdw","tc","hdl","cotinine","egfr","pir","edu","pa_any","drinks")
print(round(sapply(A[, ..vars], function(x) mean(is.na(x))) * 100, 1))
cat("done\n")
