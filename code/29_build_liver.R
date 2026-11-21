## ============================================================================
## 29_build_liver.R — 构建肝-脂轴数据集（ALT/AST/GGT/TB/PLT + FLI + FIB-4 + HSI）
##   8 个周期：C 用 L40_C/L25_C，其余用 BIOPRO_*/CBC_*
##   起点：analysis.rds 的 S0（成人、协变量齐备，N=28,342）
##   输出：data/derived/liver.rds
## ============================================================================
suppressMessages({library(foreign); library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))

## ---------- 周期 → 文件 stem -------------------------------------------------
STEM <- data.table(
  cycle = c("C","D","F","G","H","I","J","L"),
  bio   = c("L40_C","BIOPRO_D","BIOPRO_F","BIOPRO_G","BIOPRO_H","BIOPRO_I","BIOPRO_J","BIOPRO_L"),
  cbc   = c("L25_C","CBC_D","CBC_F","CBC_G","CBC_H","CBC_I","CBC_J","CBC_L")
)
VARS <- c(alt = "LBXSATSI", ast = "LBXSASSI", ggt = "LBXSGTSI",
          alp2 = "LBXSAPSI", tb = "LBXSTB", tp = "LBXSTP",
          alb2 = "LBXSAL", tg2 = "LBXSTR")

pick <- function(f, map) {
  p <- file.path(BASE, "data/raw", paste0(f, ".xpt"))
  if (!file.exists(p)) { cat("missing", f, "\n"); return(NULL) }
  x <- as.data.table(read.xport(p))
  keep <- intersect(unname(map), names(x))
  if (!length(keep)) return(NULL)
  out <- x[, c("SEQN", keep), with = FALSE]
  setnames(out, keep, names(map)[match(keep, map)])
  for (v in names(map)) if (!v %in% names(out)) out[, (v) := NA_real_]
  out[, c("SEQN", names(map)), with = FALSE]
}

LB <- rbindlist(lapply(seq_len(nrow(STEM)), function(i) {
  s <- STEM$cycle[i]
  b <- pick(STEM$bio[i], VARS)
  c0 <- pick(STEM$cbc[i], c(plt = "LBXPLTSI"))
  if (is.null(b)) return(NULL)
  d <- if (is.null(c0)) b else merge(b, c0, by = "SEQN", all.x = TRUE)
  d[, cycle := s]
  nas <- c(sum(!is.na(d$alt)), sum(!is.na(d$ast)), sum(!is.na(d$ggt)), sum(!is.na(d$plt)))
  cat(sprintf("%s: n=%5d  ALT=%4d AST=%4d GGT=%4d PLT=%4d\n", s, nrow(d), nas[1], nas[2], nas[3], nas[4]))
  d
}), use.names = TRUE, fill = TRUE)
cat("\n合并后实验室记录:", nrow(LB), "\n")

## ---------- 并入 S0 ----------------------------------------------------------
d <- merge(copy(A$D), LB, by = c("SEQN","cycle"), all.x = TRUE)
cat("D =", nrow(A$D), " 并入后 =", nrow(d), "\n")

## 单位与合理范围（NHANES：ALT/AST/GGT/ALP U/L；TB mg/dL；TG mg/dL；PLT 1000/uL）
d[alt <= 0 | alt > 5000, alt := NA_real_]
d[ast <= 0 | ast > 5000, ast := NA_real_]
d[ggt <= 0 | ggt > 5000, ggt := NA_real_]
d[plt <= 0 | plt > 1000, plt := NA_real_]
d[tb  <  0 | tb  > 30,   tb  := NA_real_]
d[tg2 <= 0 | tg2 > 2000, tg2 := NA_real_]

## ---------- 指数 -------------------------------------------------------------
## FLI (Bedogni 2006) —— 用全体可得的生化 TG（LBXSTR），非空腹亚样本
d[, fli_lin := 0.953*log(tg2) + 0.139*bmi + 0.718*log(ggt) + 0.053*waist - 15.745]
d[, fli := 100 * exp(fli_lin) / (1 + exp(fli_lin))]
## HSI (Lee 2010) —— 不依赖 TG，作为 FLI 的替代
d[, hsi := 8 * (alt/ast) + bmi + 2*female + 2*(diab == 1)]
## FIB-4 (Sterling 2006)
d[, fib4 := age * ast / (plt * sqrt(alt))]
## De Ritis 比值
d[, ast_alt := ast / alt]

for (v in c("alt","ast","ggt","tb","fib4","ast_alt","fli","hsi","plt")) {
  x <- d[[v]]
  cat(sprintf("%-8s nm=%5d  med=%7.2f  q10=%7.2f  q90=%7.2f\n", v, sum(!is.na(x)),
              median(x, na.rm=TRUE), quantile(x,.1,na.rm=TRUE), quantile(x,.9,na.rm=TRUE)))
}
cat("\nFLI>=60（脂肪肝）比例:", round(100*mean(d$fli >= 60, na.rm=TRUE), 1), "%\n")
cat("FIB-4>2.67（高纤维化风险）比例:", round(100*mean(d$fib4 > 2.67, na.rm=TRUE), 2), "%\n")

## 分析子样本：肝指标 + PFAS + 处理齐备
d[, ok_liver := !is.na(alt) & !is.na(ast) & !is.na(ggt) & !is.na(pfas_z) & !is.na(w5)]
cat("\n肝分析样本 N =", sum(d$ok_liver), " 处理组 =", sum(d$w5[d$ok_liver]), "\n")
print(d[ok_liver == TRUE, .(N = .N, trt = sum(w5), alt = round(median(alt),1),
                            ggt = round(median(ggt),1)), by = cycle][order(cycle)])

saveRDS(d, file.path(BASE, "data/derived/liver.rds"))
cat("\n✅ liver.rds 已保存\n")
