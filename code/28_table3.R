## ============================================================================
## 28_table3.R — 导出论文用表 3（雄激素轴）与表 4（同系物）
## ============================================================================
suppressMessages(library(data.table))
BASE <- Sys.getenv("AGE22_ROOT", "."); RES <- file.path(BASE, "results")

per5 <- fread(file.path(RES, "androgen_per5.csv"))
ax   <- fread(file.path(RES, "androgen_axis.csv"))
ter  <- fread(file.path(RES, "androgen_tertile.csv"))
cg   <- fread(file.path(RES, "androgen_congener.csv"))
dc   <- fread(file.path(RES, "androgen_decouple.csv"))
ss   <- fread(file.path(RES, "forest_seed_stability.csv"))

ORD <- c("Total testosterone","SHBG","Free testosterone","%Free testosterone",
         "FAI (T/SHBG)","Estradiol","T:E2 ratio","Albumin (control)")

T3 <- per5[match(ORD, label)]
T3[, `:=`(se = abs(per5) / qnorm(1 - p/2))]
T3[, `:=`(pct = 100*(exp(per5)-1),
          lo  = 100*(exp(per5 - 1.96*se)-1),
          hi  = 100*(exp(per5 + 1.96*se)-1))]
T3[, N := ax$N[match(outcome, ax$outcome)]]
T3 <- merge(T3[, .(outcome, label, N, pct, lo, hi, p,
                   interaction_P = pint_joint)],
            ter[, .(outcome, tertile_P = p_tert, cate_T1, cate_T2, cate_T3)], by = "outcome")
T3[, `:=`(cate_T1 = 100*(exp(cate_T1)-1), cate_T2 = 100*(exp(cate_T2)-1),
          cate_T3 = 100*(exp(cate_T3)-1))]
setnames(T3, c("pct","lo","hi","p","interaction_P","tertile_P","cate_T1","cate_T2","cate_T3"),
         c("per5_pct","per5_lo","per5_hi","per5_p","PFAS_interaction_P","PFAS_tertile_P",
           "tertile1_pct","tertile2_pct","tertile3_pct"))
T3[, `:=`(per5_pct = round(per5_pct, 2), per5_lo = round(per5_lo, 2),
          per5_hi = round(per5_hi, 2),
          per5_p = signif(per5_p, 3), PFAS_interaction_P = signif(PFAS_interaction_P, 3),
          PFAS_tertile_P = signif(PFAS_tertile_P, 3),
          tertile1_pct = round(tertile1_pct, 1), tertile2_pct = round(tertile2_pct, 1),
          tertile3_pct = round(tertile3_pct, 1))]
fwrite(T3, file.path(RES, "table3_androgen_axis.csv"))

T4 <- cg[outcome %in% c("y_tt","y_shbg","y_ft","y_fai")]
T4[, `:=`(outcome = c(y_tt="Total testosterone", y_shbg="SHBG", y_ft="Free testosterone",
                      y_fai="FAI")[outcome],
          pct = round(100*(exp(b1_std)-1), 1),
          pct_lo = round(100*(exp(b1_std - 1.96*se_b1_std)-1), 1),
          pct_hi = round(100*(exp(b1_std + 1.96*se_b1_std)-1), 1))]
setnames(T4, "p_joint", "joint_P")
T4 <- T4[order(factor(congener, c("PFOS","PFOA","PFHxS","PFNA")), factor(outcome))]
fwrite(T4[, .(congener, outcome, N, pct, pct_lo, pct_hi, joint_P = signif(joint_P, 3))],
       file.path(RES, "table4_congener.csv"))

cat("== Table 3 ==\n"); print(T3)
cat("\n== Table 4 ==\n"); print(T4[, .(congener, outcome, pct, pct_lo, pct_hi, joint_P)])
cat("\n== 脱钩检验 ==\n"); print(dc)
cat("\n== 种子稳定性 ==\n")
print(ss[, .(sample, seeds = .N, p_min = round(min(p_hetero),3), p_max = round(max(p_hetero),3),
             n_sig = sum(p_hetero < .05)), by = sample])
