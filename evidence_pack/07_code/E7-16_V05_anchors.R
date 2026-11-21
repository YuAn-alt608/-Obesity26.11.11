## =====================================================================
## V05  第二层：金标准锚点比对（与 NCHS 公布数字对齐）
## ---------------------------------------------------------------------
## 逻辑：如果清洗/合并/加权任何一步有错，加权患病率不可能同时命中
##       官方公布值。这里比对 4 个可直接查证的公开锚点 + 1 个人口总量锚点。
## 锚点来源
##   肥胖 2015-2016 = 39.8%（粗率，20 岁以上）NCHS Data Brief No.288
##   肥胖 2017-2018 年龄调整 = 42.4%（粗率 42.5%）Data Brief No.360
##   重度肥胖 2015-2016 = 7.7%；2017-2018 = 9.2%
##   2011-2012 男性总睾酮 10-90 百分位 150-698 ng/dL（Vesper, Clin Chem 2015）
## 产出：verification/anchor_benchmarks.csv、anchor_detail.csv
## =====================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", ".")
OUT  <- file.path(BASE, "verification")
A <- readRDS(file.path(BASE, "data/derived/analysis.rds"))
D <- as.data.table(A$D)
H <- as.data.table(readRDS(file.path(BASE, "data/derived/horm.rds")))

## ---------- 设计加权比率 + Taylor 线性化 SE --------------------------------
## !! 变量名陷阱：derived 里的 wt_meas 实为 BMXWT（实测体重 kg），
##    抽样权重是 wt_mec（WTMEC2YR）。这里全部用 wt_mec。
wratio <- function(d, y, w = "wt_mec") {
  d <- d[is.finite(get(y)) & is.finite(get(w)) & get(w) > 0]
  if (!nrow(d)) return(c(est = NA, se = NA, n = 0))
  d[, `:=`(yy = get(y), ww = get(w))]
  est <- d[, sum(ww * yy) / sum(ww)]
  ## PSU 层内合计 -> 以"独立 PSU"近似（NHANES 推荐做法）
  psu_tot <- d[, .(num = sum(ww * yy), den = sum(ww)), by = .(strata, psu)]
  psu_tot[, r := num - est * den]
  G <- nrow(psu_tot); n <- nrow(d)
  v <- (G / (G - 1)) * sum(psu_tot$r^2) / (sum(d$ww)^2)
  c(est = est, se = sqrt(v), n = n)
}

## 年龄直接标准化（2000 年美国标准人口，20-39/40-59/60+）
STD <- c(`20-39` = 0.39639, `40-59` = 0.35577, `60+` = 0.24784)
ageadj <- function(d, y, w = "wt_mec") {
  d <- d[is.finite(get(y)) & get(w) > 0]
  d[, ag := cut(age, c(19, 39, 59, 200), c("20-39", "40-59", "60+"))]
  s <- d[, as.list(wratio(.SD, y, w)), by = ag]
  s <- s[!is.na(est)]
  sum(s$est * STD[s$ag]) / sum(STD[s$ag])
}

res <- list()
add <- function(label, est, se, ref, tol, unit, note = "") {
  res[[length(res) + 1]] <<- data.table(anchor = label, estimate = est, se = se,
                                        published = ref, diff = est - ref,
                                        tol = tol, unit = unit,
                                        pass = abs(est - ref) <= tol, note = note)
}

## ---------- A. 2015-2016（周期 I）------------------------------------------
dI <- D[cycle == "I" & age >= 20 & !is.na(bmi) & is.finite(wt_mec)]
dI[, ob := as.numeric(bmi >= 30)]
dI[, sev := as.numeric(bmi >= 40)]
r1 <- wratio(dI, "ob"); add("肥胖率 2015-2016（粗率）", 100 * r1["est"], 100 * r1["se"], 39.8, 0.5, "%")
r2 <- wratio(dI, "sev"); add("重度肥胖 2015-2016（粗率）", 100 * r2["est"], 100 * r2["se"], 7.7, 0.5, "%")
a1 <- ageadj(dI, "ob");   add("肥胖率 2015-2016（年龄调整）", 100 * a1, NA, 39.6, 0.5, "%", "DB360 趋势表")

## ---------- B. 2017-2018（周期 J）------------------------------------------
dJ <- D[cycle == "J" & age >= 20 & !is.na(bmi) & is.finite(wt_mec)]
dJ[, ob := as.numeric(bmi >= 30)]
dJ[, sev := as.numeric(bmi >= 40)]
r3 <- wratio(dJ, "ob"); add("肥胖率 2017-2018（粗率）", 100 * r3["est"], 100 * r3["se"], 42.5, 0.5, "%")
a2 <- ageadj(dJ, "ob"); add("肥胖率 2017-2018（年龄调整）", 100 * a2, NA, 42.4, 0.5, "%")

## ---------- C. 加权人口总量（权重是否可用） --------------------------------
popI <- dI[is.finite(wt_mec), sum(wt_mec) / 1e6]; popJ <- dJ[is.finite(wt_mec), sum(wt_mec) / 1e6]
## 权重和应落在美国 20+ 非机构化人口的合理区间（约 2.3 亿），仅作合理性区间检验
add("加权成人人口 2015-2016", popI, NA, 232, 6, "百万", "权重和；合理区间 226-238")
add("加权成人人口 2017-2018", popJ, NA, 235, 6, "百万", "权重和；合理区间 229-241")

## ---------- D. 男性总睾酮（周期 H = 2013-2014）-----------------------------
h <- H[cycle == "H" & sex == 1 & age >= 20 & !is.na(tt_ngdl)]
tt_med <- median(h$tt_ngdl)
tt_q <- quantile(h$tt_ngdl, c(.1, .9), names = FALSE)
add("男性总睾酮中位数 2013-2014（20 岁以上）", tt_med, NA, 400, 20, "ng/dL",
    sprintf("10-90 百分位 %.0f-%.0f；文献 NHANES 2011-12 为 150-698", tt_q[1], tt_q[2]))
hI <- H[cycle == "I" & sex == 1 & age >= 20 & !is.na(tt_ngdl)]
add("男性总睾酮中位数 2015-2016（20 岁以上）", median(hI$tt_ngdl), NA, 400, 20, "ng/dL")

AB <- rbindlist(res)
AB[, estimate := round(estimate, 3)]; AB[, diff := round(diff, 3)]
fwrite(AB, file.path(OUT, "anchor_benchmarks.csv"))

## 明细：分年龄组肥胖率、睾酮分布
det <- rbind(
  dI[, .(cycle = "I", group = "total", prev = 100 * sum(wt_mec * ob) / sum(wt_mec))],
  dJ[, .(cycle = "J", group = "total", prev = 100 * sum(wt_mec * ob) / sum(wt_mec))],
  rbind(dI[, .(cycle = "I", group = as.character(cut(age, c(19,39,59,200), c("20-39","40-59","60+"))),
               ob = ob, w = wt_mec)][, .(prev = 100 * sum(w * ob) / sum(w)), by = .(cycle, group)],
        dJ[, .(cycle = "J", group = as.character(cut(age, c(19,39,59,200), c("20-39","40-59","60+"))),
               ob = ob, w = wt_mec)][, .(prev = 100 * sum(w * ob) / sum(w)), by = .(cycle, group)])
)
fwrite(det, file.path(OUT, "anchor_detail.csv"))

cat("\n========== V05 金标准锚点 ==========\n")
print(AB[, .(anchor, estimate, published, diff, tol, pass, note)], row.names = FALSE)
cat("\n命中率：", sum(AB$pass), "/", nrow(AB), "\n")
cat("\n睾酮分布（男性，周期 H）：\n"); print(summary(h$tt_ngdl))
cat("\n%游离睾酮中位数（全部男性，周期 H）：", round(median(H[cycle=="H" & sex==1]$pct_free, na.rm=TRUE), 2), "%\n")
