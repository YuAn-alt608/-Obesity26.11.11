## V10 第一步：导出三张核心表的"干净中间表"（脱敏，仅保留分析所需列）
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", "."); OUT <- file.path(BASE, "verification")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds")); CONF <- A$CONF

export <- function(D, y, tag, extra = NULL) {
  d <- D[is.finite(get(y)) & is.finite(dw5) & !is.na(strata) & !is.na(psu)]
  cols <- unique(c(y, "dw5", CONF, "strata", "psu", extra))
  cols <- cols[cols %in% names(d)]
  x <- d[, ..cols]
  x[, id := .I]
  fwrite(x, file.path(OUT, paste0("_py_", tag, ".csv")))
  cat(sprintf("%-12s N=%6d  y=%s\n", tag, nrow(x), y))
}

## (1) 雄激素：清洁样本（排除影响雄激素的药物），总睾酮每 5% 减重
H0 <- as.data.table(readRDS(file.path(BASE, "data/derived/horm.rds")))
RX <- as.data.table(readRDS(file.path(BASE, "data/derived/rx_flags.rds")))
D0 <- A$D
H <- D0[H0, on = "SEQN"]          # 取 dw1 与混杂变量
H <- RX[H, on = "SEQN"]
H[, `:=`(castration = fifelse(is.na(castration), 0L, castration),
         suppress   = fifelse(is.na(suppress), 0L, suppress))]
H <- H[female == 0 & age >= 20 & castration == 0 & suppress == 0]
H[, dw5 := dw1 / 5]
H[, y_tt := log(tt_ngdl)]
export(H, "y_tt", "androgen")

## (2) 肝-脂轴：FLI 每 5% 减重
L <- as.data.table(readRDS(file.path(BASE, "data/derived/liver.rds")))
L <- L[ok_liver == TRUE & !is.na(w5) & is.finite(pfas_z)]
L[, dw5 := dw1 / 5]; L[, y_fli := log(pmax(fli, 0.1))]
export(L, "y_fli", "liver")

## (3) 系统炎症：WBC 每 5% 减重
I <- as.data.table(readRDS(file.path(BASE, "data/derived/inflam.rds")))
I <- I[ok_infl == TRUE & !is.na(w5) & is.finite(pfas_z)]
I[, dw5 := dw1 / 5]
export(I, "y_wbc", "inflam")
