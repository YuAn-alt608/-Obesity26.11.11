BASE <- Sys.getenv("AGE22_ROOT", ".")
## 21_geo_burden.R — 组织层面净炎症负荷 (细胞比例 × 单细胞炎症评分)
suppressMessages(library(data.table))
OUT <- file.path(BASE, "results")
comp <- fread(file.path(OUT, "geo_composition.csv"))
sub  <- fread(file.path(OUT, "geo_adipocyte_subscores.csv"))
mac  <- fread(file.path(OUT, "geo_macrophage_modules.csv"))

cat("--- macrophage per-cell inflamma (module z) ---\n")
print(mac[, .(mean_inflam = mean(inflamma), sd = sd(inflamma)), by = group], digits = 3)

f <- comp[celltype == "macrophage", .(frac = mean(pct) / 100), by = .(gsm, group)]
m <- mac[, .(gsm, minflam = inflamma)]
b <- merge(f, m, by = "gsm")
b[, burden := frac * minflam]
cat("\n--- macrophage fraction x per-cell inflammation ---\n")
print(b[, .(burden = mean(burden), se = sd(burden) / sqrt(.N)), by = group], digits = 3)
cat("P(WL vs Obese) =",
    signif(t.test(b$burden[b$group == "WL"], b$burden[b$group == "Obese"])$p.value, 3), "\n")

fa <- comp[celltype == "adipocyte", .(frac = mean(pct) / 100), by = .(gsm, group)]
ad <- sub[, .(gsm, ainflam = inflam)]
ba <- merge(fa, ad, by = "gsm")
ba[, burden := frac * ainflam]
cat("\n--- adipocyte fraction x adipocyte inflammation (log10 CPM) ---\n")
print(ba[, .(burden = mean(burden)), by = group], digits = 3)
cat("P(WL vs Obese) =",
    signif(t.test(ba$burden[ba$group == "WL"], ba$burden[ba$group == "Obese"])$p.value, 3), "\n")
