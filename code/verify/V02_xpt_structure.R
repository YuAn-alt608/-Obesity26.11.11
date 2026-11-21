## =====================================================================
## V02  第一层：XPT 结构完整性与跨文件主键链接
## ---------------------------------------------------------------------
## 排查三件事：
##   (1) 文件能否解析  -> 截断/半包下载的 XPT 会在 read.xport 处直接失败
##   (2) 文件内 SEQN 是否唯一 -> 决定后续聚类稳健 SE 是否失真
##   (3) 每个文件的人是否都在同周期 DEMO 里 -> 检测"下错周期/串周期"
## 产出：verification/xpt_structure.csv、verification/seqn_linkage.csv
## =====================================================================
suppressMessages({library(foreign); library(data.table); library(digest)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
RAW  <- file.path(BASE, "data/raw")
OUT  <- file.path(BASE, "verification"); dir.create(OUT, showWarnings = FALSE)

fs <- sort(list.files(RAW, pattern = "\\.xpt$", ignore.case = TRUE))
cat("发现 XPT 文件:", length(fs), "个\n")

rd <- function(f) {
  x <- try(read.xport(file.path(RAW, f)), silent = TRUE)
  if (inherits(x, "try-error")) return(NULL)
  as.data.table(x)
}

rows <- vector("list", length(fs)); link <- vector("list", length(fs))
for (i in seq_along(fs)) {
  f  <- fs[i]
  sfx <- sub("^.*_([A-Z])\\.xpt$", "\\1", f)
  x  <- rd(f)
  if (is.null(x)) {
    rows[[i]] <- data.table(file = f, suffix = sfx, ok = FALSE, nrow = NA_integer_, ncol = NA_integer_,
                            n_seqn = NA_integer_, dup_seqn = NA_integer_, seqn_min = NA_integer_,
                            seqn_max = NA_integer_, content_sha256 = NA_character_, note = "PARSE_FAIL")
    next
  }
  has <- "SEQN" %in% names(x)
  sq  <- if (has) x$SEQN else integer(0)
  rows[[i]] <- data.table(
    file = f, suffix = sfx, ok = TRUE, nrow = nrow(x), ncol = ncol(x),
    n_seqn = length(unique(sq)), dup_seqn = sum(duplicated(sq)),
    seqn_min = if (has) min(sq, na.rm = TRUE) else NA_integer_,
    seqn_max = if (has) max(sq, na.rm = TRUE) else NA_integer_,
    content_sha256 = digest(x, algo = "sha256"),
    note = if (!has) "NO_SEQN" else "")
  if (has) link[[i]] <- data.table(file = f, suffix = sfx, seqn = sq)
  if (i %% 40 == 0) cat("  已解析", i, "/", length(fs), "\n")
}
ST <- rbindlist(rows, fill = TRUE)
fwrite(ST, file.path(OUT, "xpt_structure.csv"))

## ---------- 跨文件主键链接：同周期 DEMO 为全集 ----------
LK <- rbindlist(link, fill = TRUE)
demo <- ST[grepl("^DEMO_", file), .(suffix, nrow)]
dm <- LK[grepl("^DEMO_", file), .(suffix, seqn)]
dm[, in_demo := TRUE]
LK <- merge(LK, dm, by = c("suffix", "seqn"), all.x = TRUE)
res <- LK[, .(n = uniqueN(seqn), n_orphan = sum(is.na(in_demo))), by = .(file, suffix)]
res <- merge(res, demo[, .(suffix, demo_n = nrow)], by = "suffix", all.x = TRUE)
res[, pct_orphan := 100 * n_orphan / n]
setorder(res, -pct_orphan)
fwrite(res, file.path(OUT, "seqn_linkage.csv"))

## ---------- 断言 ----------
cat("\n========== V02 结论 ==========\n")
bad_parse <- ST[ok == FALSE]
cat("解析失败文件          :", nrow(bad_parse), if (nrow(bad_parse)) paste(bad_parse$file, collapse = ", ") else "", "\n")
MULTI_OK <- grepl("^RXQ_RX_", ST$file)   # 处方文件按处方逐条记录，人均多行是设计而非错误
dup <- ST[!is.na(dup_seqn) & dup_seqn > 0]
dup_bad <- dup[!MULTI_OK[match(dup$file, ST$file)]]
cat("一对多文件（允许）    :", nrow(dup[dup$file %in% ST$file[MULTI_OK]]), "个 RXQ_RX 文件\n")
if (nrow(dup)) {
  d <- ST[file %in% dup$file, .(file, nrow, n_seqn, dup_seqn)]
  d[, per_person := round(nrow / n_seqn, 2)]
  cat("   处方条数/人：均值", round(mean(d$per_person), 2), "最大", max(d$per_person), "\n")
}
cat("非预期一对多文件      :", nrow(dup_bad), if (nrow(dup_bad)) paste(dup_bad$file, collapse = ", ") else "（无）", "\n")
no_seqn <- ST[note == "NO_SEQN"]
cat("无 SEQN 列的文件      :", nrow(no_seqn), if (nrow(no_seqn)) paste(no_seqn$file, collapse = ", ") else "（无）", "\n")
cat("SEQN 不在同周期 DEMO  :", sum(res$n_orphan > 0), "个文件\n")
if (any(res$n_orphan > 0)) print(res[n_orphan > 0][order(-pct_orphan)][1:min(10, .N)])
cat("\nXPT 合计行数:", format(sum(ST$nrow, na.rm = TRUE), big.mark = ","), "\n")
stopifnot(nrow(bad_parse) == 0, nrow(dup_bad) == 0)
cat("✅ 全部 XPT 解析成功；除处方文件外文件内主键唯一\n")
