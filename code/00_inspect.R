BASE <- Sys.getenv("AGE22_ROOT", ".")
suppressMessages({library(haven)})
lp <- function(f) {
  d <- read_xpt(file.path(BASE, "data/raw", paste0(f, ".xpt")))
  cat("=====", f, " N=", nrow(d), " vars=", ncol(d), "\n", sep = "")
  for (v in names(d)) {
    lab <- attr(d[[v]], "label"); if (is.null(lab)) lab <- ""
    cat(sprintf("  %-10s %-10s %s\n", v, class(d[[v]])[1], lab))
  }
  invisible(d)
}
for (f in c("WHQ_C", "SSCHL_C", "L28NPB_C", "L24PFC_C")) lp(f)
