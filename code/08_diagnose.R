suppressMessages({library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
CORE <- readRDS(file.path(BASE, "data/derived/core.rds"))
A <- CORE[adult == TRUE]
cat("=== wt_intent / tried_lose 编码 ===\n")
print(A[, .(n = .N, intent_ok = sum(!is.na(wt_intent)), tried_ok = sum(!is.na(tried_lose)),
            intent_yes = sum(wt_intent == 1, na.rm = TRUE),
            tried_yes = sum(tried_lose == 1, na.rm = TRUE)), by = cycle][order(cycle)])
cat("\n=== dw1>=5 中 有意/尝试 比例 ===\n")
print(A[!is.na(dw1) & dw1 >= 5, .(loss5 = .N,
        tried = sum(tried_lose == 1, na.rm = TRUE),
        intent = sum(wt_intent == 1, na.rm = TRUE)), by = cycle][order(cycle)])
cat("\n=== dw1 分布 (分周期) ===\n")
print(round(A[, .(n = sum(!is.na(dw1)), mean = mean(dw1, na.rm = TRUE),
                  p10 = quantile(dw1, .1, na.rm = TRUE), med = median(dw1, na.rm = TRUE),
                  p90 = quantile(dw1, .9, na.rm = TRUE), pct5 = mean(dw1 >= 5, na.rm = TRUE) * 100),
                by = cycle][order(cycle)], 2))
cat("\n=== dw10 (10年) 可用性 ===\n")
print(A[, .(n10 = sum(!is.na(dw10)), mean10 = round(mean(dw10, na.rm = TRUE), 2)), by = cycle][order(cycle)])
cat("\n=== 男性游离睾酮 (仅 sex==1) ===\n")
H <- readRDS(file.path(BASE, "data/derived/horm.rds"))
print(summary(H[sex == 1]$ft_ngdl)); print(summary(H[sex == 1]$tt_ngdl))
cat("\n=== 男性睾酮与体重史交集 ===\n")
M <- merge(H[sex == 1], A[, .(SEQN, dw1, phenoage, tried_lose)], by = "SEQN")
print(M[, .(N = .N, T = sum(!is.na(tt_ngdl)), FT = sum(!is.na(ft_nmol)),
            dw1 = sum(!is.na(dw1)), both = sum(!is.na(tt_ngdl) & !is.na(dw1))), by = cycle][order(cycle)])
