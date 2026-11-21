## =====================================================================
## V07  第二层：版本撕裂审计（Audit Trail）+ 运行环境快照
## ---------------------------------------------------------------------
## 解决的问题："改了上游 05 清洗脚本，却忘记重跑下游 22/34 统计脚本"
##   1) 对每个脚本，抓出它写出的产物文件名（fwrite/saveRDS/ggsave 里的字面量）
##   2) 比对产物 mtime 与脚本 mtime：产物更旧 => 版本撕裂，标 STALE
##   3) 记录运行环境（R 版本、关键包版本、随机种子）供复现声明引用
## 产出：verification/audit_trail.csv、environment_snapshot.txt
## =====================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", ".")
OUT  <- file.path(BASE, "verification")
CODE <- file.path(BASE, "code")
DIRS <- c(results = file.path(BASE, "results"), figures = file.path(BASE, "figures"),
          derived = file.path(BASE, "data/derived"))

scripts <- sort(list.files(CODE, pattern = "\\.R$"))
rows <- list()
for (s in scripts) {
  p <- file.path(CODE, s)
  txt <- paste(readLines(p, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  ## 只认"写"操作（fwrite / saveRDS / ggsave），排除 readRDS / fread 的误抓
  ln <- strsplit(txt, "\n")[[1]]
  ln <- ln[grepl("fwrite\\(|saveRDS\\(|ggsave\\(", ln)]
  lits <- unlist(regmatches(ln, gregexpr('"[^"]+\\.(csv|rds|png|pdf)"', ln)))
  lits <- unique(basename(gsub('"', "", lits)))
  if (!length(lits)) next
  for (f in lits) {
    hit <- sapply(DIRS, function(d) {
      q <- file.path(d, f)
      if (file.exists(q)) as.numeric(file.mtime(q)) else NA_real_
    })
    if (all(is.na(hit))) next
    k <- which.max(replace(hit, is.na(hit), -Inf))
    rows[[length(rows) + 1]] <- data.table(
      script = s, script_mtime = as.numeric(file.mtime(p)),
      output = f, dir = names(DIRS)[k], output_mtime = hit[[k]],
      stale = hit[[k]] < as.numeric(file.mtime(p)))
  }
}
AT <- rbindlist(rows, fill = TRUE)
AT[, script_time := format(as.POSIXct(script_mtime, origin = "1970-01-01"), "%Y-%m-%d %H:%M")]
AT[, output_time := format(as.POSIXct(output_mtime, origin = "1970-01-01"), "%Y-%m-%d %H:%M")]
AT[, lag_min := round((output_mtime - script_mtime) / 60, 1)]
setorder(AT, script, output)
fwrite(AT[, .(script, script_time, output, dir, output_time, lag_min, stale)],
       file.path(OUT, "audit_trail.csv"))

## 环境快照
pk <- c("data.table","ggplot2","patchwork","foreign","MASS","digest","grf","survival",
        "matrixStats","scales","RColorBrewer","readxl","openxlsx")
si <- sessionInfo()
env <- c(
  paste("R version      :", R.version.string),
  paste("Platform       :", si$platform),
  paste("Running under  :", si$running),
  paste("Generated at   :", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "Packages:",
  sprintf("  %-14s %s", pk, sapply(pk, function(x) {
    if (requireNamespace(x, quietly = TRUE)) as.character(packageVersion(x)) else "NOT INSTALLED"})),
  "",
  "Python: 3.x with statsmodels / pandas / numpy / scipy（见 python -c 版本输出）",
  "",
  "Random seeds used in pipeline: set.seed(2026) for grf; 20261117 for verification sampling;",
  "  seed-stability diagnostics use seeds 1..5 in 25_seed_stability.R")
writeLines(env, file.path(OUT, "environment_snapshot.txt"))

cat("\n========== V07 版本撕裂 ==========\n")
cat("脚本数:", length(unique(AT$script)), " 产物记录:", nrow(AT), "\n")
st <- AT[stale == TRUE]
cat("产物早于脚本（撕裂）:", nrow(st), " / ", nrow(AT), "\n")
if (nrow(st)) {
  ## 只列出"落后超过 2 分钟"的，避免"改图例后立即重跑"造成的噪声
  bad <- st[lag_min < -2]
  cat("其中落后 >2 分钟（需重跑确认）:", nrow(bad), "\n")
  if (nrow(bad)) print(bad[order(lag_min)][1:min(12, .N), .(script, output, lag_min)])
}
cat("\n各脚本最旧产物落后（分钟）：\n")
print(AT[, .(min_lag = min(lag_min), n_out = .N), by = script][order(min_lag)][1:min(8, .N)])
cat("\n环境快照 -> verification/environment_snapshot.txt\n")
