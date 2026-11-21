## =====================================================================
## V03  第一层：单位一致性 + 缺失编码（哨兵值）审计
## ---------------------------------------------------------------------
## A. 官方单位：从 CDC 码本页抓取每个变量的 SAS 标签（标签里带单位）
## B. 实测分布：逐周期 P25/中位数/P75/max，自动标出"跨周期 10 倍量级跳跃"
## C. 哨兵值：7/9/77/99/777/999/7777/9999/77777/99999 在衍生数据集中的残留
## D. 人体测量学上界：身高/体重/腰围/1年前体重是否混入 7777/9999
## 产出：verification/official_units.csv、unit_distribution.csv、
##       sentinel_audit.csv、anthropometric_bounds.csv
## =====================================================================
suppressMessages({library(foreign); library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
RAW  <- file.path(BASE, "data/raw")
OUT  <- file.path(BASE, "verification"); dir.create(OUT, showWarnings = FALSE)
CACHE <- file.path(OUT, "_codebook"); dir.create(CACHE, showWarnings = FALSE)

YEAR <- c(C = "2003", D = "2005", F = "2009", G = "2011",
          H = "2013", I = "2015", J = "2017", L = "2021")
url_of <- function(f) {
  sfx <- sub("^.*_([A-Za-z]+)$", "\\1", sub("\\.xpt$", "", f))
  sfx <- toupper(sfx)
  if (!sfx %in% names(YEAR)) return(NA_character_)
  sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%s/DataFiles/%s.htm", YEAR[[sfx]], sub("\\.xpt$", "", f))
}

## ---------- A. 官方单位 ------------------------------------------------------
KEY <- c("PFC_D","PFC_F","PFC_G","PFAS_H","PFAS_I","PFAS_J","PFAS_L","L24PFC_C",
         "TST_G","TST_H","TST_I","SSTST_J",
         "BIOPRO_D","BIOPRO_F","BIOPRO_G","BIOPRO_H","BIOPRO_I","BIOPRO_J","BIOPRO_L",
         "CBC_D","CBC_F","CBC_G","CBC_H","CBC_I","CBC_J","CBC_L","L25_C",
         "CRP_D","CRP_F","HSCRP_I","HSCRP_J","HSCRP_L",
         "BMX_C","BMX_D","BMX_F","BMX_G","BMX_H","BMX_I","BMX_J","BMX_L",
         "DEMO_C","DEMO_D","DEMO_F","DEMO_G","DEMO_H","DEMO_I","DEMO_J","DEMO_L",
         "WHQ_C","WHQ_D","WHQ_F","WHQ_G","WHQ_H","WHQ_I","WHQ_J","WHQ_L")
KEY <- KEY[file.exists(file.path(RAW, paste0(KEY, ".xpt")))]

parse_codebook <- function(f) {
  cf <- file.path(CACHE, paste0(f, ".htm"))
  if (!file.exists(cf)) {
    u <- url_of(f)
    if (is.na(u)) return(NULL)
    ok <- try(download.file(u, cf, quiet = TRUE), silent = TRUE)
    if (inherits(ok, "try-error") || !file.exists(cf) || file.size(cf) < 5000) return(NULL)
  }
  h <- paste(readLines(cf, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  h <- gsub("<script.*?</script>", " ", h, perl = TRUE)
  h <- gsub("<[^>]+>", "|", h)
  h <- gsub("&nbsp;|&#160;", " ", h)
  tok <- unlist(strsplit(h, "\\|"))
  tok <- trimws(tok)
  hit <- grep("^[A-Z][A-Z0-9_]{2,}\\s+-\\s+\\S", tok, value = TRUE)
  if (!length(hit)) return(NULL)
  var <- sub("\\s+-\\s+.*$", "", hit)
  lab <- sub("^[A-Z][A-Z0-9_]{2,}\\s+-\\s+", "", hit)
  unit <- ifelse(grepl("\\(([^)]+)\\)\\s*$", lab),
                 sub("^.*\\(([^)]+)\\)\\s*$", "\\1", lab), "")
  data.table(file = f, var = var, sas_label = lab, unit = unit)
}

cat("== A. 抓取官方码本（", length(KEY), " 个文件）==\n")
UO <- rbindlist(lapply(KEY, function(f) {
  x <- parse_codebook(f)
  if (is.null(x)) { cat("   [无码本]", f, "\n"); return(NULL) }
  x
}), fill = TRUE)
if (nrow(UO)) {
  UO <- unique(UO, by = c("file", "var"))
  ## 只保留本文真正用到的变量
  USED <- c("LBXPFOS","LBXPFOA","LBXPFHS","LBXPFNA","LBXPFDA","LBXPFUA","LBXPFDOA",
            "LBXTST","LBXSHBG","LBXEST","LBXE2","LBXTT","LBXALT","LBXSATSI","LBXSASSI",
            "LBXSGTSI","LBXSAPSI","LBXSTB","LBXSAL","LBXSTR","LBXPLTSI",
            "LBXWBCSI","LBXLYPCT","LBXMOPCT","LBXNEPCT","LBDLYMNO","LBDMONO","LBDNENO",
            "LBXCRP","LBXHSCRP","LBXGH","LBXGLU","LBXTC","LBDHDD","LBXVIDMS","LBXVI025",
            "BMXWT","BMXHT","BMXBMI","BMXWAIST","WHD010","WHD020","WHQ030","WHD050")
  UO <- UO[var %in% USED]
  fwrite(UO, file.path(OUT, "official_units.csv"))
  cat("   提取到", nrow(UO), "条带单位的官方变量定义\n")
}

## ---------- B. 实测分布 ------------------------------------------------------
SPEC <- list(
  list(f = c("L24PFC_C","PFC_D","PFC_F","PFC_G","PFAS_H","PFAS_I","PFAS_J","PFAS_L"),
       v = c("LBXPFOS","LBXPFOA","LBXPFHS","LBXPFNA"), tag = "PFAS"),
  list(f = c("TST_G","TST_H","TST_I"), v = c("LBXTST","LBXSHBG","LBXEST"), tag = "Androgen"),
  list(f = c("BIOPRO_D","BIOPRO_F","BIOPRO_G","BIOPRO_H","BIOPRO_I","BIOPRO_J","BIOPRO_L"),
       v = c("LBXSATSI","LBXSASSI","LBXSGTSI","LBXSAPSI","LBXSTB","LBXSAL","LBXSTR"), tag = "Liver"),
  list(f = c("CBC_D","CBC_F","CBC_G","CBC_H","CBC_I","CBC_J","CBC_L","L25_C"),
       v = c("LBXWBCSI","LBXLYPCT","LBXMOPCT","LBDNENO","LBXPLTSI"), tag = "CBC"),
  list(f = c("CRP_D","CRP_F","HSCRP_I","HSCRP_J","HSCRP_L"),
       v = c("LBXCRP","LBXHSCRP"), tag = "CRP")
)

rows <- list()
for (sp in SPEC) {
  for (f in sp$f) {
    p <- file.path(RAW, paste0(f, ".xpt"))
    if (!file.exists(p)) next
    x <- as.data.table(read.xport(p))
    for (v in sp$v) {
      if (!v %in% names(x)) next
      z <- x[[v]]; z <- z[is.finite(z)]
      if (!length(z)) next
      q <- quantile(z, c(.25, .5, .75), na.rm = TRUE, names = FALSE)
      rows[[length(rows) + 1]] <- data.table(
        tag = sp$tag, file = f, var = v, cycle = sub("^.*_", "", f),
        n = length(z), min = min(z), p25 = q[1], median = q[2], p75 = q[3], max = max(z),
        pct_zero = 100 * mean(z == 0), pct_below_loq_proxy = 100 * mean(z < quantile(z, .01)))
    }
  }
}
UD <- rbindlist(rows, fill = TRUE)
## 跨周期跳跃：同一变量相邻周期中位数之比
UD[, ratio_prev := median / shift(median), by = var]
UD[, flag := fifelse(!is.na(ratio_prev) & (ratio_prev > 5 | ratio_prev < 0.2),
                     sprintf("JUMP x%.1f", ratio_prev), "")]
fwrite(UD, file.path(OUT, "unit_distribution.csv"))

## ---------- C. 哨兵值残留 ----------------------------------------------------
## 判据（避免把真实取值误判为哨兵）：
##   ① 跳过设计与分类变量（srvyr/yr/cycle/strata/psu，以及唯一值<=20 的列）
##   ② 值必须远离主体分布：v > 10 x P99，或 v >= 7777
SENT <- c(7, 9, 77, 99, 777, 999, 7777, 9999, 77777, 99999)
SKIP <- c("SEQN", "srvyr", "yr", "cycle", "strata", "psu", "w5", "w10", "w5s",
          "female", "adult", "comorbid", "diab", "htn", "preg", "smk", "drink2",
          "pa_any", "smoke100", "smoke_now", "race_f", "edu_f", "stable", "crp_ok",
          "ok_now", "ok_1y", "ok_10y", "ok_liver", "ok_infl", "w5s\\.1")
RDS <- c("core", "core_raw", "pop", "horm", "liver", "inflam", "rx_flags")
srows <- list()
for (nm in RDS) {
  p <- file.path(BASE, "data/derived", paste0(nm, ".rds"))
  if (!file.exists(p)) next
  d <- as.data.table(readRDS(p))
  for (v in names(d)) {
    if (v %in% SKIP || !is.numeric(d[[v]])) next
    z <- d[[v]][is.finite(d[[v]])]
    if (length(z) < 100 || uniqueN(z) <= 20) next
    hit <- z[z %in% SENT]
    if (!length(hit)) next
    p99 <- quantile(z[!z %in% SENT], .99, names = FALSE)   # 关键：P99 也要排除哨兵
    ## 只保留"远超主体分布"的取值
    bad <- hit[hit > max(10 * p99, 7777)]
    if (!length(bad)) next
    srows[[length(srows) + 1]] <- data.table(
      dataset = nm, var = v, n = length(z), p99 = signif(p99, 4), max = max(z),
      n_sentinel = length(bad),
      values = paste(sort(unique(bad)), collapse = ";"),
      verdict = "SENTINEL_PRESENT")
  }
}
SA <- if (length(srows)) rbindlist(srows) else
  data.table(dataset = character(), var = character(), n = integer(), p99 = numeric(),
             max = numeric(), n_sentinel = integer(), values = character(), verdict = character())
fwrite(SA, file.path(OUT, "sentinel_audit.csv"))

## ---------- D. 人体测量学上界 ------------------------------------------------
d <- as.data.table(readRDS(file.path(BASE, "data/derived/core.rds")))
chk <- function(v, lo, hi) {
  if (!v %in% names(d)) return(NULL)
  z <- d[[v]][is.finite(d[[v]])]
  data.table(var = v, n = length(z), min = min(z), median = round(median(z), 2), max = max(z),
             n_below = sum(z < lo), n_above = sum(z > hi), bound = sprintf("[%g, %g]", lo, hi))
}
AB <- rbindlist(list(
  chk("ht_meas", 120, 220), chk("bmi", 12, 120), chk("waist", 40, 200),
  chk("sr_kg_now", 20, 400), chk("sr_kg_1y", 20, 400), chk("sr_kg_10y", 20, 400),
  chk("sr_wt_now", 20, 880), chk("sr_wt_1y", 20, 880), chk("sr_wt_10y", 20, 880),
  chk("dw1", -60, 60), chk("dw10", -80, 80)
), fill = TRUE)
fwrite(AB, file.path(OUT, "anthropometric_bounds.csv"))

## ---------- E. 同一指标、不同变量名的量级核对（CRP 陷阱） --------------------
crp_src <- list(
  list(f = "CRP_D",  v = "LBXCRP"),  list(f = "CRP_F",  v = "LBXCRP"),
  list(f = "HSCRP_I", v = "LBXHSCRP"), list(f = "HSCRP_J", v = "LBXHSCRP"),
  list(f = "HSCRP_L", v = "LBXHSCRP"))
ec <- rbindlist(lapply(crp_src, function(z) {
  p <- file.path(RAW, paste0(z$f, ".xpt"))
  if (!file.exists(p)) return(NULL)
  x <- as.data.table(read.xport(p))
  if (!z$v %in% names(x)) return(NULL)
  q <- quantile(x[[z$v]], c(.25, .5, .75), na.rm = TRUE, names = FALSE)
  data.table(file = z$f, var = z$v, n = sum(!is.na(x[[z$v]])),
             p25 = q[1], median = q[2], p75 = q[3])
}))
ec[, implied_unit := fifelse(var == "LBXCRP", "mg/dL（1 mg/dL = 10 mg/L）", "mg/L")]
ec[, median_mgL := fifelse(var == "LBXCRP", median * 10, median)]
fwrite(ec, file.path(OUT, "crp_unit_adjudication.csv"))

## 体重史问卷的哨兵拦截核算
W <- as.data.table(readRDS(file.path(BASE, "data/derived/core.rds")))
wq <- rbindlist(lapply(c("sr_wt_now", "sr_wt_1y", "sr_wt_10y", "sr_wt_max", "sr_ht_in"), function(v) {
  if (!v %in% names(W)) return(NULL)
  z <- W[[v]]
  data.table(var = v, n = sum(is.finite(z)),
             n_9999 = sum(z %in% 9999, na.rm = TRUE),
             n_7777 = sum(z %in% 7777, na.rm = TRUE),
             n_999 = sum(z %in% 999, na.rm = TRUE),
             max = max(z, na.rm = TRUE))
}))
fwrite(wq, file.path(OUT, "whq_sentinel_interception.csv"))

## ---------- 结论 -------------------------------------------------------------
cat("\n========== V03 结论 ==========\n")
jump <- UD[flag != ""]
cat("跨周期量级跳跃：", nrow(jump), "处\n")
if (nrow(jump)) print(jump[, .(tag, var, file, median = round(median, 3), ratio = round(ratio_prev, 2), flag)])
cat("\n哨兵值残留（远超主体分布者）：", nrow(SA), "个变量\n")
if (nrow(SA)) print(SA[order(-n_sentinel)]) else cat("   无\n")
cat("\n人体测量学上界：\n"); print(AB[!is.na(n)], row.names = FALSE)
cat("\nCRP 两代变量量级核对：\n"); print(ec[, .(file, var, median, implied_unit, median_mgL = round(median_mgL, 2))])
cat("\nWHQ 体重史哨兵计数：\n"); print(wq)
