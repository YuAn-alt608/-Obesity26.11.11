BASE <- Sys.getenv("AGE22_ROOT", ".")
## Dump variable dictionary for every downloaded XPT
suppressMessages({library(haven)})
raw <- file.path(BASE, "data/raw")
out <- file.path(BASE, "data/meta/vars_all.csv")
fs <- list.files(raw, pattern = "\\.xpt$", full.names = TRUE)
res <- list()
for (f in fs) {
  d <- tryCatch(read_xpt(f), error = function(e) NULL)
  if (is.null(d)) { cat("FAIL", basename(f), "\n"); next }
  for (v in names(d)) {
    lab <- attr(d[[v]], "label"); if (is.null(lab)) lab <- ""
    res[[length(res) + 1]] <- data.frame(
      file = sub("\\.xpt$", "", basename(f)), var = v, label = lab,
      cls = class(d[[v]])[1], n_ok = sum(!is.na(d[[v]])), n = nrow(d),
      stringsAsFactors = FALSE)
  }
  cat("OK", basename(f), nrow(d), "x", ncol(d), "\n")
}
D <- do.call(rbind, res)
write.csv(D, out, row.names = FALSE)
cat("TOTAL", nrow(D), "variables ->", out, "\n")
