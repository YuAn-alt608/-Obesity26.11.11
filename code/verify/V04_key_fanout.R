## =====================================================================
## V04  第二层：主键唯一性 + 连接扇出（fan-out）审计
## ---------------------------------------------------------------------
## 审稿人最怕的两类静默错误：
##   (1) 主键重复 -> 聚类稳健 SE 把同一 PSU 算成多个独立簇，P 值虚假变小
##   (2) 多对多 join -> 样本量虚增（例如把处方逐条表直接 join 到人级表）
## 产出：verification/key_integrity.csv、verification/join_fanout.csv
## =====================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", ".")
OUT  <- file.path(BASE, "verification")

LD <- function(nm) as.data.table(readRDS(file.path(BASE, "data/derived", paste0(nm, ".rds"))))
A  <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
DAT <- list(core = LD("core"), core_raw = LD("core_raw"), pop = LD("pop"),
            horm = LD("horm"), liver = LD("liver"), inflam = LD("inflam"),
            rx_flags = LD("rx_flags"))
for (nm in c("D", "S0", "S1", "S2", "S2b", "S3", "S4")) DAT[[paste0("analysis$", nm)]] <- as.data.table(A[[nm]])

## ---------- A. 主键断言 ------------------------------------------------------
ki <- rbindlist(lapply(names(DAT), function(nm) {
  d <- DAT[[nm]]
  data.table(dataset = nm, nrow = nrow(d), n_seqn = uniqueN(d$SEQN),
             dup = anyDuplicated(d$SEQN), n_na = sum(is.na(d$SEQN)))
}))
ki[, unique_ok := dup == 0 & n_na == 0]
fwrite(ki, file.path(OUT, "key_integrity.csv"))
cat("========== A. 主键唯一性 ==========\n"); print(ki)

## 硬断言：任何一张主表出现重复即中止
bad <- ki[unique_ok == FALSE]
if (nrow(bad)) { print(bad); stop("主键重复：聚类标准误将失真") }
cat("✅ 7 个衍生数据集 + 7 个分析子集的主键全部唯一\n\n")

## ---------- B. 连接扇出审计 --------------------------------------------------
## 规则：人级表 join 人级表，行数必须不变；一旦行数增加即为多对多扇出
pairs <- list(
  c("inflam", "analysis$D"), c("liver", "analysis$D"), c("core", "analysis$D"),
  c("pop", "core"), c("horm", "core"), c("rx_flags", "analysis$D"),
  c("inflam", "rx_flags"), c("liver", "rx_flags")
)
fo <- rbindlist(lapply(pairs, function(p) {
  l <- DAT[[p[1]]]; r <- DAT[[p[2]]]
  m <- merge(l[, .(SEQN)], r[, .(SEQN)], by = "SEQN")
  data.table(left = p[1], right = p[2], n_left = nrow(l), n_right = nrow(r),
             n_inner = nrow(m), fanout = nrow(m) > nrow(l),
             dup_inner = anyDuplicated(m$SEQN) > 0)
}))
fwrite(fo, file.path(OUT, "join_fanout.csv"))
cat("========== B. 连接扇出 ==========\n"); print(fo)
stopifnot(!any(fo$fanout), !any(fo$dup_inner))
cat("✅ 全部两两连接均未产生行数膨胀或主键重复\n\n")

## ---------- C. 集合嵌套关系 --------------------------------------------------
rel <- data.table(
  check = c("horm ⊆ pop", "pop ⊆ core", "horm ⊆ core", "liver ⊆ analysis$D", "inflam ⊆ analysis$D"),
  ok = c(all(DAT$horm$SEQN %in% DAT$pop$SEQN),
         all(DAT$pop$SEQN %in% DAT$core$SEQN),
         all(DAT$horm$SEQN %in% DAT$core$SEQN),
         all(DAT$liver$SEQN %in% DAT[["analysis$D"]]$SEQN),
         all(DAT$inflam$SEQN %in% DAT[["analysis$D"]]$SEQN)))
print(rel)
stopifnot(all(rel$ok))
cat("✅ 各层样本嵌套关系正确（不存在串周期的游离记录）\n")
