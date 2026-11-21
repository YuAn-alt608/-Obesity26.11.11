## ============================================================================
## 06_build_pop.R — 持久性有机污染物暴露负荷
##   · ΣPFAS: 6 种主要同系物摩尔和 (nmol/L), 覆盖 8 个周期
##   · ΣPCB / ΣOCP: 2003-2004 个体级, 脂质校正 (ng/g lipid)
##   输出 data/derived/pop.rds, results/qc_pop.csv
## ============================================================================
suppressMessages({library(haven); library(data.table)})
BASE <- Sys.getenv("AGE22_ROOT", "."); RAW <- file.path(BASE, "data/raw")
CORE <- readRDS(file.path(BASE, "data/derived/core.rds"))
rd <- function(s) as.data.table(suppressWarnings(read_xpt(file.path(RAW, paste0(s, ".xpt")))))
gv <- function(d, cands) { for (v in cands) if (!is.null(d) && v %in% names(d)) return(as.numeric(d[[v]])); rep(NA_real_, nrow(d)) }

## ---------------- PFAS ----------------
MW <- c(PFOA = 414.07, PFOS = 500.13, PFHxS = 400.11, PFNA = 464.08, PFDA = 514.11, PFUnDA = 564.14)
PF <- list()
for (s in c("C","D","F","G","H","I","J","L")) {
  f <- switch(s, C = "L24PFC_C", D = "PFC_D", F = "PFC_F", G = "PFC_G",
              H = "PFAS_H", I = "PFAS_I", J = "PFAS_J", L = "PFAS_L")
  d <- rd(f); if (is.null(d)) next
  n <- nrow(d)
  pfoa <- switch(s,
    C = , D = , F = , G = gv(d, "LBXPFOA"),
    H = gv(d, "SSNPFOA") + gv(d, "SSBPFOA"),
    I = , J = , L = gv(d, "LBXNFOA") + gv(d, "LBXBFOA"))
  pfos <- switch(s,
    C = , D = , F = , G = gv(d, "LBXPFOS"),
    H = gv(d, "SSNPFOS") + gv(d, "SSMPFOS"),
    I = , J = , L = gv(d, "LBXNFOS") + gv(d, "LBXMFOS"))
  ## 2013-2014 的 PFOA/PFOS 在 surplus 文件中按 SEQN 补入
  if (s == "H") {
    ss <- rd("SSPFAS_H")
    tp <- data.table(SEQN = as.numeric(ss$SEQN),
                     pfoa_ss = gv(ss, "SSNPFOA") + gv(ss, "SSBPFOA"),
                     pfos_ss = gv(ss, "SSNPFOS") + gv(ss, "SSMPFOS"))
    tmp <- data.table(SEQN = as.numeric(d$SEQN), pfoa_i = pfoa, pfos_i = pfos)
    tmp <- merge(tmp, tp, by = "SEQN", all.x = TRUE)
    pfoa <- ifelse(is.na(tmp$pfoa_i), tmp$pfoa_ss, tmp$pfoa_i)
    pfos <- ifelse(is.na(tmp$pfos_i), tmp$pfos_ss, tmp$pfos_i)
  }
  X <- data.table(SEQN = as.numeric(d$SEQN), cycle = s,
                  PFOA = pfoa, PFOS = pfos,
                  PFHxS = gv(d, "LBXPFHS"), PFNA = gv(d, "LBXPFNA"),
                  PFDA = gv(d, "LBXPFDE"), PFUnDA = gv(d, "LBXPFUA"))
  X[, wt_pf := gv(d, c("WTSPF2YR", "WTSB2YR", "WTSA2YR", "WTSC2YR"))]
  ## 低于 LOD 的数保留 NHANES 填报值; 仅剔除负值
  for (v in names(MW)) X[get(v) <= 0, (v) := NA_real_]
  PF[[s]] <- X
  cat("PFAS", s, "N =", nrow(X), " 完整6项 =", sum(complete.cases(X[, names(MW), with=FALSE])), "\n")
}
PFA <- rbindlist(PF, fill = TRUE)

## 摩尔和 (nmol/L): ng/mL ÷ (g/mol) × 1000
PFA[, n_ok := rowSums(!is.na(.SD)), .SDcols = names(MW)]
PFA[, paste0("chk_", names(MW)) := NULL]
PFA[, pfas_molar := {
  m <- matrix(NA_real_, .N, length(MW))
  for (j in seq_along(MW)) m[, j] <- get(names(MW)[j]) / MW[j] * 1000
  rowSums(m, na.rm = TRUE)
}]
PFA[n_ok < 6, pfas_molar := NA_real_]
PFA[, `:=`(log_pfas = log(pfas_molar), pfas_wet = rowSums(.SD, na.rm = TRUE)),
    .SDcols = names(MW)]
PFA[n_ok < 6, pfas_wet := NA_real_]
## 周期内 z 标准化（跨周期可比）
PFA[, pfas_z := as.numeric(scale(log_pfas)), by = cycle]
PFA[, `:=`(pfos_z = as.numeric(scale(log(PFOS))), pfoa_z = as.numeric(scale(log(PFOA)))), by = cycle]

## ---------------- PCB / OCP (2003-2004) ----------------
## 注意: PCB(子样本C) 与 OCP(子样本B) 为互斥子样本, 不能合并为单一指标;
##       PCB 与共平面 PCB/二噁英(同一批 2285 人)可合并为总 PCB 负荷。
pcb <- rd("L28NPB_C"); ocp <- rd("L28OCP_C"); dfp <- rd("L28DFP_C")
pcb_lip  <- intersect(grep("^LB[XD][0-9]{3}LA$", names(pcb), value = TRUE), names(pcb))
ocp_lip  <- intersect(paste0(c("LBXHCB","LBXBHC","LBXGHC","LBXPDE","LBXPDT","LBXODT",
                               "LBXOXY","LBXTNA","LBXHPE","LBXMIR","LBXALD"), "LA"), names(ocp))
co_lip   <- grep("^LB[XD][0-9]{3}LA$", names(dfp), value = TRUE)  # 单邻位共平面PCB
cat("\nPCB(NDL) n=", length(pcb_lip), " 共平面PCB n=", length(co_lip),
    " OCP n=", length(ocp_lip), "\n")

POS <- data.table(SEQN = as.numeric(pcb$SEQN), wt_pop = as.numeric(pcb$WTSC2YR),
                  pcb_lipid  = rowSums(as.matrix(pcb[, ..pcb_lip]), na.rm = TRUE),
                  pcb_n      = rowSums(!is.na(as.matrix(pcb[, ..pcb_lip]))))
if (length(co_lip)) {
  DP <- data.table(SEQN = as.numeric(dfp$SEQN),
                   co_lipid = rowSums(as.matrix(dfp[, ..co_lip]), na.rm = TRUE))
  POS <- merge(POS, DP, by = "SEQN", all.x = TRUE)
  POS[, pcb_lipid := pcb_lipid + fifelse(is.na(co_lipid), 0, co_lipid)]
}
OC <- data.table(SEQN = as.numeric(ocp$SEQN),
                 ocp_lipid = rowSums(as.matrix(ocp[, ..ocp_lip]), na.rm = TRUE),
                 ocp_n     = rowSums(!is.na(as.matrix(ocp[, ..ocp_lip]))))
POS <- merge(POS, OC, by = "SEQN", all = TRUE)
POS[pcb_n < 10, pcb_lipid := NA_real_]
POS[ocp_n < 5,  ocp_lipid := NA_real_]
POS[, `:=`(log_pcb = log(pcb_lipid), log_ocp = log(ocp_lipid))]

## ---------------- 合并 ----------------
POP <- merge(CORE[, .(SEQN, cycle)], PFA[, !c("n_ok"), with = FALSE], by = c("SEQN","cycle"), all.x = TRUE)
POP <- merge(POP, POS, by = "SEQN", all.x = TRUE)
POP[, `:=`(pcb_z = as.numeric(scale(log_pcb)), ocp_z = as.numeric(scale(log_ocp)))]

saveRDS(POP, file.path(BASE, "data/derived/pop.rds"))

qc <- POP[, .(N = .N, pfas6 = sum(!is.na(pfas_molar)), pcb = sum(!is.na(pcb_lipid)),
              ocp = sum(!is.na(ocp_lipid)),
              pfas_pfos = sum(!is.na(PFOS))), by = cycle][order(cycle)]
print(qc); fwrite(qc, file.path(BASE, "results/qc_pop.csv"))
cat("\n=== ΣPFAS 分布 (nmol/L) ===\n"); print(summary(POP$pfas_molar))
cat("\n=== 与体重史交集 ===\n")
A <- merge(POP, CORE[, .(SEQN, adult, dw1, phenoage)], by = "SEQN")
print(A[adult == TRUE & !is.na(dw1), .(N = .N, pfas = sum(!is.na(pfas_molar)),
      pcb = sum(!is.na(pcb_lipid)), phenoage = sum(!is.na(phenoage)),
      pfas_pheno = sum(!is.na(pfas_molar) & !is.na(phenoage))), by = cycle][order(cycle)])
cat("done\n")
