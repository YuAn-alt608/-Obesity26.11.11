BASE <- Sys.getenv("AGE22_ROOT", ".")
## V06 第一步：导出激素原始三元组与 R 端算出的游离睾酮
suppressMessages(library(data.table))
H <- as.data.table(readRDS(file.path(BASE, "data/derived/horm.rds")))
H <- H[!is.na(tt_ngdl) & !is.na(shbg) & !is.na(alb_gL) & !is.na(ft_nmol)]
set.seed(20261117)
H[, shbg_t := cut(shbg, quantile(shbg, c(0, 1/3, 2/3, 1), na.rm = TRUE), labels = c("low","mid","high"), include.lowest = TRUE)]
pick <- H[, .SD[sample(.N, min(.N, 4))], by = shbg_t]
fwrite(H[, .(SEQN, cycle, age, sex, tt_ngdl, shbg, alb_gL, ft_nmol, ft_ngdl, pct_free)],
       file.path(BASE, "verification/_horm_full.csv"))
fwrite(pick[, .(SEQN, cycle, age, sex, tt_ngdl, shbg, alb_gL, ft_nmol, ft_ngdl, pct_free, shbg_t)],
       file.path(BASE, "verification/_horm_sample.csv"))
cat("导出:", nrow(H), "行；抽样", nrow(pick), "行\n")
cat("R 端 %游离睾酮中位数:", round(median(H$pct_free), 2), "%\n")
