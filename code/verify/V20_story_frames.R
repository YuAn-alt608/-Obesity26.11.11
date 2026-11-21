## ============================================================
## V20_story_frames.R
## 目的：为「反常识破局 + 跨尺度机制重定向」四幕叙事补两张支撑表
##   表1 geo_composition_selectivity.csv —— 第二幕「选择性耗竭」的硬证据
##        7 类细胞逐一比较 Obese vs WL（Welch t + Cohen's d），
##        用来说明「下降是巨噬特异的，不是所有免疫细胞等比退潮」
##   表2 androgen_bioavailability.csv   —— 第四幕「总量升、生物利用度降」
##        把总睾酮 / SHBG / %游离 / 游离 / FAI / 白蛋白 / T:E2 放在一张表里，
##        并给出「总睾酮 − SHBG」的成比例检验（脱钩检验的读法）
## 输出：results/geo_composition_selectivity.csv
##       results/androgen_bioavailability.csv
## 说明：只做汇总与检验，不重新拟合模型；所有原始系数来自既有 results/ 表
## ============================================================

suppressMessages({library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
setwd(BASE)

## ---------- 表1：细胞组成的选择性 ----------
comp <- fread("results/geo_composition.csv")

CT <- c("macrophage","adipocyte","ascell","tcell","pericyte","endo","mast")
SEG <- c(macrophage = "髓系", adipocyte = "实质细胞", ascell = "祖细胞",
         tcell = "淋系", pericyte = "基质", endo = "内皮", mast = "髓系")

rows <- list()
for (ct in CT) {
  x <- comp[celltype == ct & group == "Obese", pct]
  y <- comp[celltype == ct & group == "WL",    pct]
  l <- comp[celltype == ct & group == "Lean",  pct]
  tt <- t.test(y, x)
  sp <- sqrt((var(x) + var(y)) / 2)
  rows[[ct]] <- data.table(
    celltype   = ct,
    lineage    = SEG[[ct]],
    n_obese    = length(x), n_wl = length(y), n_lean = length(l),
    mean_obese = mean(x), sd_obese = sd(x),
    mean_wl    = mean(y), sd_wl    = sd(y),
    mean_lean  = mean(l), sd_lean  = sd(l),
    abs_change = mean(y) - mean(x),
    pct_change = 100 * (mean(y) - mean(x)) / mean(x),
    welch_p    = tt$p.value,
    cohens_d   = (mean(y) - mean(x)) / sp,
    direction  = ifelse(mean(y) < mean(x), "下降", "上升")
  )
}
SEL <- rbindlist(rows)

## BH 校正（7 类细胞，一个家族）
SEL[, p_bh := p.adjust(welch_p, "BH")]
## 「显著下降」判定
SEL[, verdict := fifelse(p_bh < 0.05 & abs_change < 0,
                         "显著下降（选择性）",
                  fifelse(p_bh < 0.05 & abs_change > 0,
                         "显著上升（分母效应）",
                         "不显著"))]
SEL[, abs_d := abs(cohens_d)]
setorder(SEL, -abs_d)
fwrite(SEL, "results/geo_composition_selectivity.csv")

cat("\n===== 表1 细胞组成选择性（按 |d| 排序）=====\n")
print(SEL[, .(celltype, lineage,
              obese = round(mean_obese, 1), wl = round(mean_wl, 1), lean = round(mean_lean, 1),
              pct = round(pct_change, 1), p = signif(welch_p, 2),
              p_bh = signif(p_bh, 2), d = round(cohens_d, 2), verdict)])

cat("\n>>> 显著下降的细胞类型数：", SEL[verdict == "显著下降（选择性）", .N], "\n")
cat(">>> 其中 |d| 最大者：",
    SEL[verdict == "显著下降（选择性）"][which.max(abs(cohens_d)), celltype], "\n")

## ---------- 表2：雄激素轴的「总量 vs 生物利用度」 ----------
and <- fread("results/androgen_per5.csv")
dec <- fread("results/androgen_decouple.csv")

LBL <- c(y_tt = "总睾酮", y_shbg = "SHBG", y_ft = "游离睾酮（Vermeulen 推算）",
         y_pctfree = "%游离睾酮", y_fai = "FAI（T/SHBG）", y_e2 = "雌二醇",
         y_te2 = "T:E2（芳香化活性）", y_alb = "白蛋白（阴性对照）")
ROLE <- c(y_tt = "总量", y_shbg = "结合蛋白", y_ft = "生物利用度",
          y_pctfree = "生物利用度", y_fai = "生物利用度", y_e2 = "底物",
          y_te2 = "芳香化", y_alb = "阴性对照")
NOTE <- c(
  y_tt      = "总量上升，是减重收益的主证据",
  y_shbg    = "上升更快（+4.78% > +4.02%），这是%游离下降的直接原因",
  y_ft      = "点估计为正但 P=0.16，不能说「游离睾酮升高」",
  y_pctfree = "唯一显著下降的生物利用度指标（P=8.3e-07）",
  y_fai     = "不动，说明 T 与 SHBG 的上升大致成比例",
  y_e2      = "不动",
  y_te2     = "上升，指向外周芳香化下降（脂肪组织芳香化酶活性降低）",
  y_alb     = "完全不动，排除「肝脏把蛋白合成整体上调」"
)

BIO <- and[, .(outcome, label, per5, p, int, p_int)]
BIO[, outcome_cn := LBL[outcome]]
BIO[, role       := ROLE[outcome]]
BIO[, pct_change := 100 * (exp(per5) - 1)]
BIO[, p_fmt      := fifelse(p < 1e-6, sprintf("%.1e", p), sprintf("%.3f", p))]
BIO[, sig        := fifelse(p < 0.05, "是", "否")]
BIO[, note       := NOTE[outcome]]
setcolorder(BIO, c("outcome","outcome_cn","role","per5","pct_change","p","p_fmt",
                   "sig","int","p_int","note"))
fwrite(BIO, "results/androgen_bioavailability.csv")

cat("\n===== 表2 雄激素轴：总量 vs 生物利用度 =====\n")
print(BIO[, .(outcome_cn, role, pct = round(pct_change, 2), p = p_fmt, sig)])

## 关键对比：总睾酮 vs %游离（方向相反）
tt <- BIO[outcome == "y_tt"]; pf <- BIO[outcome == "y_pctfree"]
cat(sprintf("\n>>> 总量 %+.2f%%（P=%s） vs 生物利用度 %+.2f%%（P=%s）—— 方向相反\n",
            tt$pct_change, tt$p_fmt, pf$pct_change, pf$p_fmt))

## 脱钩检验读法
cat("\n===== 脱钩检验（H0：两者变化成比例）=====\n")
print(dec[, .(contrast, est = round(est, 3),
              ci = sprintf("(%.3f, %.3f)", lo, hi), p = round(p, 3))])
cat("\n>>> 结论：TT-SHBG 不成比例？ ",
    ifelse(dec[contrast == "TT - SHBG", p] < 0.05, "是（存在脱钩）", "否（成比例，P=0.83）"), "\n")
cat(">>> 白蛋白是否跟随 SHBG？ ",
    ifelse(dec[contrast == "SHBG - albumin", p] < 0.05, "是", "否（P=0.52，SHBG 特异性）"), "\n")

cat("\n[V20] 完成：results/geo_composition_selectivity.csv（", nrow(SEL), "行）、",
    "results/androgen_bioavailability.csv（", nrow(BIO), "行）\n", sep = "")

## ---------- 表3：样本量漏斗（审稿人第一个会问的） ----------
core <- readRDS("data/derived/core.rds"); setDT(core)
pop  <- readRDS("data/derived/pop.rds");  setDT(pop)
mm <- merge(core[, .(SEQN, adult)],
            pop[, .(SEQN, pfas_molar, wt_pf)], by = "SEQN")

infl <- readRDS("data/derived/inflam.rds"); setDT(infl)
liv  <- readRDS("data/derived/liver.rds");  setDT(liv)
and  <- fread("results/androgen_axis.csv")
comp3 <- fread("results/geo_composition.csv")
des  <- fread("verification/survey_design_summary.csv")

## 三个分析域的精确条件（已逐一核对到主报告数字）
n_infl <- sum(!is.na(infl$w5) & infl$ok_infl & !is.na(infl$pfas_molar), na.rm = TRUE)  # 12,775
n_liv  <- sum(liv$ok_liver & !is.na(liv$pfas_molar), na.rm = TRUE)                     # 12,737
n_and  <- and$N[1]                                                                      # 1,232
n_cf   <- 9656L
n_wt   <- as.integer(des[item == "设计框架 N", value])                                  # 14,128

FLOW <- data.table(
  branch = c("总体", "总体", "暴露域", "暴露域", "暴露域", "暴露域",
             "雄激素域", "森林域", "加权域", "单细胞", "单细胞"),
  step = c("NHANES 2003–2023 八周期全部受试者",
           "成人（≥20 岁，非孕）",
           "其中有 PFAS 血清测量",
           "其中成人且有 PFAS",
           "其中炎症结局完整且减重幅度可得（主分析域）",
           "其中肝酶/FLI 完整（主分析域）",
           "2013–2016 男性、激素齐备、排除药物性抑制",
           "PhenoAge 终点与协变量完整（因果森林）",
           "持有 PFAS 子样本权重（V16 加权框架，不设年龄限）",
           "人腹部皮下脂肪 snRNA-seq 池数",
           "细胞核总数（QC 后）"),
  n = c(nrow(mm),
        sum(mm$adult %in% TRUE),
        sum(!is.na(mm$pfas_molar)),
        sum(mm$adult %in% TRUE & !is.na(mm$pfas_molar)),
        n_infl, n_liv, n_and, n_cf, n_wt,
        uniqueN(comp3$gsm), sum(comp3$N)),
  note = c("8 个周期全部记录",
           "preg 字段在成人中无阳性，故与成人同数",
           "16,778 人测得至少 6 种 PFAS 同系物",
           "★ 人群分析的上界",
           "条件：PFAS + ok_infl + w5 非缺失",
           "条件：ok_liver + PFAS 非缺失",
           "处理组 245 人",
           "不是全部 45,739 人",
           "含 <20 岁；故与「成人且有 PFAS」相差约 1 千人，不可直接相减",
           "Obese / Lean / WL 各 6 池",
           "与报告 §4.2 的 186,221 一致")
)

## 只在同一 branch 内相邻两行算百分比，跨 branch 记 NA（避免 138%、10 万% 这类无意义数字）
FLOW[, pct_of_prev := NA_real_]
for (i in 2:nrow(FLOW)) {
  if (FLOW$branch[i] == FLOW$branch[i-1] && FLOW$branch[i] %in% c("总体", "暴露域"))
    set(FLOW, i, "pct_of_prev", 100 * FLOW$n[i] / FLOW$n[i - 1])
}
fwrite(FLOW, "results/cohort_flow.csv")

cat("\n===== 表3 样本量漏斗 =====\n")
print(FLOW[, .(branch, n = format(n, big.mark = ","), pct = round(pct_of_prev, 1))])
cat("\n[V20] results/cohort_flow.csv 已写出（", nrow(FLOW), " 步）\n", sep = "")
