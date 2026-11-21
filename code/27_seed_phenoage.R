## 27_seed_phenoage.R — PhenoAge 森林异质性检验的种子稳定性（与 §3.8 同一逻辑）
suppressMessages({library(data.table); library(grf)})
BASE <- Sys.getenv("AGE22_ROOT", ".")
A <- readRDS(file.path(BASE,"data/derived/analysis.rds")); CONF <- A$CONF
mkX <- function(dt){ X <- cbind(as.matrix(dt[,..CONF]), model.matrix(~race_f+edu_f+smk-1,data=dt))
  X <- cbind(X, as.matrix(dt[,.(pfas_z)])); mc <- grep("^miss_",names(dt),value=TRUE)
  if(length(mc)) X <- cbind(X, as.matrix(dt[,..mc])); storage.mode(X)<-"double"; X }
one <- function(seed){
  set.seed(seed); dt <- copy(A$S1)[!is.na(phenoage_accel)&!is.na(pfas_z)&!is.na(w5)]
  dt[, mv := as.numeric(pfas_z)]
  X <- mkX(dt); Y <- dt$phenoage_accel; W <- as.numeric(dt$w5)
  cl <- as.integer(factor(paste(dt$strata,dt$psu)))
  e0 <- as.numeric(regression_forest(X,W,num.trees=1500,tune.parameters="all",seed=seed)$predictions)
  k <- e0>=0.05 & e0<=0.95; X<-X[k,];Y<-Y[k];W<-W[k];e0<-e0[k];cl<-cl[k];dt<-dt[k]
  cf <- causal_forest(X,Y,W,W.hat=e0,num.trees=3000,honesty=TRUE,tune.parameters="all",clusters=cl,seed=seed)
  cal <- test_calibration(cf); blp <- best_linear_projection(cf, A=cbind(mv=dt$mv,mv2=dt$mv^2))
  data.table(seed=seed, N=nrow(dt),
             p_hetero=cal["differential.forest.prediction","Pr(>t)"],
             blp_mv=blp["mv","Estimate"], p_blp_mv=blp["mv","Pr(>|t|)"],
             blp_mv2=blp["mv2","Estimate"], p_blp_mv2=blp["mv2","Pr(>|t|)"])
}
R <- rbindlist(lapply(c(20261117,1,2), one)); print(R); fwrite(R, file.path(BASE,"results/forest_seed_phenoage.csv"))
