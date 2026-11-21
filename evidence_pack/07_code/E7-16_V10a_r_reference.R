## V10a  R 端参考实现（与 Python 完全相同的输入 CSV）
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", "."); OUT <- file.path(BASE, "verification")
crsV <- function(X, e, cl) {
  k <- ncol(X); n <- nrow(X); bread <- solve(crossprod(X))
  idx <- split(seq_len(n), cl); meat <- matrix(0, k, k)
  for (g in idx) { s <- crossprod(X[g, , drop = FALSE], e[g]); meat <- meat + s %*% t(s) }
  G <- length(idx); bread %*% meat %*% bread * (G/(G-1)) * ((n-1)/(n-k))
}
fits <- list(
  androgen = list(f = "y_tt ~ dw5", tag = "Total testosterone (log), per 5%"),
  liver    = list(f = "y_fli ~ dw5", tag = "FLI (log), per 5%"),
  inflam   = list(f = "y_wbc ~ dw5", tag = "WBC (log), per 5%"))
res <- rbindlist(lapply(names(fits), function(nm) {
  d <- fread(file.path(OUT, paste0("_py_", nm, ".csv")))
  vars <- all.vars(as.formula(fits[[nm]]$f))
  rhs <- paste(fits[[nm]]$f, "+ age + female + pir + bmi_p + log_cot + drink2 + pa_any + comorbid + tc + hdl + egfr + yr")
  mf <- model.frame(as.formula(rhs), data = d, na.action = na.omit)
  X <- model.matrix(as.formula(rhs), data = mf)
  ## 保留截距（V16 审计：初版 sd(v)>0 被判据误删常数列，R 与 Python 同时犯同一错，
  ## 因此上一轮的"完全一致"只证明了实现一致，未证明设定正确）
  kp <- apply(X, 2, function(v) sd(v) > 0) | colnames(X) == "(Intercept)"
  X <- X[, kp, drop = FALSE]
  y <- model.response(mf)
  d <- d[as.integer(rownames(mf))]
  b <- qr.coef(qr(X), y); V <- crsV(X, as.numeric(y - X %*% b), factor(paste(d$strata, d$psu)))
  j <- which(colnames(X) == "dw5")
  data.table(dataset = nm, term = fits[[nm]]$tag, N = nrow(d), n_cluster = uniqueN(paste(d$strata, d$psu)),
             beta = unname(b[j]), se_cr1 = unname(sqrt(V[j, j])),
             p_cr1 = unname(2 * pnorm(-abs(b[j] / sqrt(V[j, j])))),
             pct_per5 = unname(100 * (exp(b[j]) - 1)))
}))
fwrite(res, file.path(OUT, "python_crosscheck_r.csv"))
print(res[, .(dataset, N, n_cluster, beta = signif(beta, 6), se_cr1 = signif(se_cr1, 6),
              p_cr1 = signif(p_cr1, 4), pct_per5 = round(pct_per5, 3))], row.names = FALSE)
