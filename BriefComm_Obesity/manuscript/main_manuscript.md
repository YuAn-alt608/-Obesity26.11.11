# Two central predictions of the weight-loss toxin-rebound narrative fail in NHANES 2003–2023

**Running title:** Testing the toxin-rebound narrative

**Design:** cross-sectional analysis of eight NHANES cycles with one-year weight change and mortality follow-up.

**Author:** Zizheng Yuan¹

**Affiliation:** ¹Doctoral Programme in Planetary Health Studies, University of Lisbon, Lisbon 1649-013, Portugal

**Corresponding author:** Zizheng Yuan, Doctoral Programme in Planetary Health Studies, University of Lisbon, Lisbon 1649-013, Portugal. E-mail: yuan-zizheng@edu.ulisboa.pt

**Funding:** This research received no specific grant from any funding agency in the public, commercial or not-for-profit sectors.

**Disclosure:** The author declares no conflict of interest.

**Ethics:** Secondary analysis of de-identified, publicly available data; exempt from ethics review.

**Data availability:** NHANES data are publicly available from the US National Center for Health Statistics. All analytic code, the derived tables used here, the verification ledger linking each reported estimate to its source file, and the numbered evidence pack are openly available at https://github.com/YuAn-alt608/-Obesity26.11.11.

**Word count:** abstract 198; main text 1,498; 18 references; 2 figures; 1 table.

---

## Abstract

**Objective:** Weight loss is said to mobilise stored pollutants, raising their serum concentration and accelerating biological ageing. I tested both claims in one national dataset, and asked how easily the same data can be made to support them.

**Methods:** Among 12,806 adults with serum PFAS and one-year weight change, and 9,656 with PhenoAge (NHANES 2003–2023, eight cycles), I compared adults losing ≥5% of body weight with those who did not, using covariate-adjusted models with cluster-robust variance, equivalence testing, and four routes to a spurious interaction.

**Results:** Weight losers appeared 1.88 years biologically older (P=3×10⁻¹⁶) and had 10.9% lower serum PFAS (P=1×10⁻⁷). Adjustment removed both: −0.18 years (95% CI −0.58 to 0.22) and −1.2% (−4.2 to 1.7). Equivalence testing excluded harm beyond one year and any rise in PFAS beyond 5%. Confounding, androgen-suppressing prescriptions, body-mass-bearing formulas and a specification error each generated a significant interaction on its own, and they did not share a direction.

**Conclusions:** Neither core prediction of the toxin-rebound narrative holds in NHANES. Alarms of this kind can be generated from real data by routes that leave no consistent directional signature, so the existence of such an alarm carries little evidential weight.

---

## Study Importance

**What is already known about this subject?** Weight loss raises serum concentrations of organochlorines, and this observation has been extended to per- and polyfluoroalkyl substances (PFAS). NHANES analyses report associations between PFAS and sex hormones,⁶⁻¹⁰ but few account for the drugs that suppress the androgen axis.

**What are the new findings?** In 12,806 adults, weight loss was associated with lower, not higher, serum PFAS, and with no difference in biological age once confounding was addressed; equivalence testing excluded a rise in PFAS beyond 5% and harm beyond one year of biological ageing. Four distinct analytic routes each produced a significant pollutant interaction, in inconsistent directions.

**How might these results change research or practice?** A pollutant–weight interaction should not be interpreted without prescription screening, formula decomposition and a check that the signal survives dropping a single survey cycle. For clinicians, total testosterone overstates the androgen gain from weight loss.

---

## Introduction

Weight loss mobilises lipophilic organochlorines stored in adipose tissue, and their serum concentrations rise as weight falls.¹ That observation has been extended to per- and polyfluoroalkyl substances (PFAS), which are not stored in fat but bind serum proteins, and a clinical narrative has grown around it: weight loss releases a toxic burden that accelerates biological ageing and blunts the hormonal benefit of losing weight.

The narrative makes two predictions that can be tested in the same dataset. If released pollutant drives harm, adults who lose weight should have higher serum pollutant concentrations, and they should appear biologically older. Both are falsifiable with a cross-sectional design, because both compare people who lost weight with those who did not.

I tested both predictions in eight NHANES cycles, then asked a second question: how much of an apparent "toxin rebound" can the data-generating process itself produce? I examined four routes — confounding, unmeasured prescribing, formulas that contain body mass, and a specification error in my own first pass — and checked whether they share a direction.

## Methods

**Data.** Eight NHANES cycles between 2003 and 2023 (2007–2008 was not fielded; 82,096 records) were restricted to non-pregnant adults aged ≥20 years (n=45,739). The PFAS domain comprised 12,806 adults with at least six congeners and computable one-year weight change; the PhenoAge domain, 9,656 with complete biomarkers (Table S1). PFAS were measured in 16,778 participants, 13,162 of them adults.

**Exposure and outcomes.** Weight change was the difference between examined BMI and BMI one year earlier, taken from self-reported prior weight, expressed as the percentage lost; ≥5% defined the exposed group. Outcomes were PhenoAge acceleration²,³ and the molar sum of six PFAS congeners (PFOA, PFOS, PFHxS, PFNA, PFDA, PFUnDA), with individual congeners examined separately; free testosterone was calculated as described by Vermeulen et al.¹³

**Estimation.** Linear models with cluster-robust variance at the survey design level (strata × PSU) were fitted for both outcomes, with covariates for age, sex, family income-to-poverty ratio, BMI one year earlier, cotinine, alcohol, physical activity, comorbidity count, total cholesterol, HDL, eGFR, survey year, race, education and smoking. Because prior BMI lies on the path from weight change to current BMI, both models are reported. Survey-weighted re-estimation, E-values¹⁵ and quantile g-computation¹⁶ served as sensitivity analyses, and causal forests¹⁴ were fitted for the androgen outcome. Equivalence was assessed by two one-sided tests (TOST) at pre-specified margins of 0.5 and 1.0 years for PhenoAge and 5% and 10% for PFAS.

**Four routes to an interaction.** Confounding was quantified as the shift between crude and adjusted estimates. Prescribing was addressed by flagging GnRH agonists and antagonists, antiandrogens, oestrogens, systemic glucocorticoids, opioids and ketoconazole in the prescription files, and by re-estimating after excluding users. Formula dependence was addressed by re-estimating the Hepatic Steatosis Index¹² and Fatty Liver Index¹¹ with their BMI terms removed (sections S2–S5). The fourth route was a specification error in the first pass of these analyses, described below and in S6.

Reporting follows STROBE.⁴

## Results

**Both predictions appear to hold, in opposite directions.** Before adjustment, adults who had lost ≥5% of body weight were 1.88 years biologically older than those who had not (P=3×10⁻¹⁶), and their serum PFAS was 10.9% lower (P=1×10⁻⁷). The first reproduces the published alarm, the second contradicts it, and neither survives adjustment.

**One imbalance explains both.** Weight losers had a mean BMI one year earlier of 32.4 kg/m² against 27.7 in those who did not, but a current BMI of 29.8 against 29.1: they are people who fell from a higher weight. They were also more often diabetic (19.6% vs 10.9%), hypertensive (44.5% vs 34.1%) and smokers (22.0% vs 18.1%), and less often educated beyond high school (19.6% vs 26.3%) (Table 1). Adjusted for this profile, the ageing difference became −0.18 years (95% CI −0.58 to 0.22; P=0.38) and the PFAS difference −1.2% (−4.2 to 1.7; P=0.41). The shifts were 2.06 years and 9.7 percentage points, in opposite directions, from the same covariate set.

**The residue is excluded, not merely undetected.** Two one-sided tests rejected harm larger than one year of PhenoAge acceleration (P=3.5×10⁻⁵) and a rise in serum PFAS larger than 5% (P=0.0060). A margin of 0.5 years was not excluded (P=0.060), so the ageing bound is one year rather than half a year. Individual congeners were unremarkable (all P≥0.099), and the dose gradient was not monotonic: adults losing ≥10% had 5.9% lower PFAS (P=0.005), those losing 5–10% showed no difference (+0.5%, P=0.80), and those gaining ≥5% also had lower PFAS (−3.9%, P=0.006).

**Three further routes each generate an interaction.** Excluding 167 of 2,244 men (7.4%) on androgen-suppressing drugs removed an already marginal forest heterogeneity (one-sided P=0.055; 0.55 after excluding the 131 opioid users alone, 0.67 after all 167). The 167 were mostly opioid users (131), with 54 on 5α-reductase inhibitors, 39 on glucocorticoids, three on ketoconazole, two on antiandrogens, one on a GnRH agonist and none on oestrogen. Removing the BMI terms from the Hepatic Steatosis Index moved its interaction P from 1.2×10⁻⁶ to 0.68, and from the Fatty Liver Index from 0.017 to 0.14. A specification error in my first pass — a zero-inflated collinear term with the intercept dropped — produced five nominal interactions among 30; three survived coefficient-level randomisation, and dropping a single cycle moved forest heterogeneity from 0.10 to 0.70.

**The four routes do not share a direction.** The same specification error that inflated apparent harm for most outcomes attenuated it for the systemic immune-inflammation index, which fell from −8.0% to −1.0%. Confounding pushed the ageing estimate towards harm and the PFAS estimate towards benefit. Direction is therefore not a diagnostic of artefact.

**What the androgen data do and do not show.** Per 5% of weight lost, total testosterone rose 4.0% (95% CI 1.7–6.4) and SHBG rose 4.8% (3.2–6.4), so the free fraction fell 2.3% (−3.2 to −1.4) while calculated free testosterone did not change detectably (+1.6%, −0.7 to 3.9; P=0.16). The free fraction is a deterministic function of SHBG and albumin, so the fall is arithmetic; the free-minus-total difference was −2.3 percentage points (P=1.2×10⁻⁶), identical to the change in the free percentage by construction. In adipose tissue from 18 donors, weight loss reduced macrophage abundance by 46.5% relative to obese controls (BH-adjusted P=0.035) but left a senescence module still above the lean group (P=0.024).

## Discussion

The two predictions failed differently, and that difference is informative. Ageing appeared worse in weight losers; pollutant burden appeared better. Both were products of the same imbalance between people who fell from a higher weight and people who did not. A crude estimate in this design carries the direction of the covariate–outcome relationship, not the direction of weight loss.

The four routes share one feature: none requires fabrication. Confounding is inherent to the question, prescribing is recorded in the same files as the exposure, and a formula containing BMI will transmit any weight-related difference into an interaction. The specification error was mine, and I report it because it also shows that the direction of a spurious signal is set by the error, not by biology: the same mistake inflated harm in most outcomes and benefit in one.

This bears on how the literature should be read. If alarms can be produced in either direction by ordinary analytic choices, then the number of published alarms is not a measure of the evidence behind them, and a significant interaction should be treated as a hypothesis that has not yet been screened. In twelve NHANES studies of persistent pollutants and sex hormones, two screened for any androgen-suppressing medication and none for GnRH agonists or opioid analgesics, the two classes with the largest suppressive effect on testosterone (Table S7).

For clinicians the surviving observation is a reading instruction. Weight loss lowers insulin, and lower insulin releases hepatic repression of SHBG,⁵ partly reversing the adiposity-driven suppression of the axis.¹⁷,¹⁸ Because the free fraction is a deterministic function of SHBG and albumin, a 4.8% rise in SHBG against a 4.0% rise in total testosterone must lower the free percentage by about 2.3% with albumin unchanged. The increment in total testosterone is therefore largely bound hormone: a value rising from 300 to 312 ng/dL while the free fraction falls has not demonstrated recovery of gonadal function. Calculated free testosterone itself did not move, so this is a caution about interpretation rather than a newly identified harm.

Five limitations bound these conclusions. PFAS and weight were measured once, so the design describes differences between people, not change within them; a cohort with repeated measures could still detect a release effect too small or too slow to appear here. Prior weight was self-reported. Free testosterone was calculated, not measured. Prescription records exclude discontinued drugs and lack drug names for 2021–2023, which has no testosterone measurement either. Finally, the forest estimates were not reproducible run to run, so they are reported only as a marginal signal; the parametric estimates carry the argument.

---

## Table 1. The imbalance behind both predictions, and what remains after adjustment

| Characteristic | No loss ≥5% (n=7,779) | Loss ≥5% (n=1,877) |
|---|---|---|
| BMI one year earlier, kg/m² | 27.7 (6.2) | 32.4 (7.8) |
| Current BMI, kg/m² | 29.1 (7.0) | 29.8 (6.9) |
| Diabetes, % | 10.9 | 19.6 |
| Hypertension, % | 34.1 | 44.5 |
| Current smoker, % | 18.1 | 22.0 |
| College or above, % | 26.3 | 19.6 |
| Family income-to-poverty ratio | 2.7 (1.6) | 2.4 (1.5) |
| PhenoAge acceleration, years | −0.5 (7.6) | 1.4 (9.0) |
| Sum PFAS, nmol/L, median (IQR) | 28.7 (14.8–53.3) | 25.0 (12.7–48.1) |

| Estimate | Crude | Adjusted | Shift | Equivalence |
|---|---|---|---|---|
| PhenoAge acceleration, years | +1.88 | −0.18 (−0.58 to 0.22) | 2.06 | harm >1 year excluded (P=3.5×10⁻⁵) |
| Serum sum PFAS, % | −10.9 | −1.2 (−4.2 to 1.7) | 9.7 | rise >5% excluded (P=0.0060) |

Values are mean (SD) unless stated. Adjusted models include age, sex, income-to-poverty ratio, BMI one year earlier, cotinine, alcohol, physical activity, comorbidity count, total cholesterol, HDL, eGFR, survey year, race, education and smoking, with cluster-robust variance at the stratum × PSU level. Equivalence by two one-sided tests at the stated margin.

## Figure legends

**Figure 1. Both predictions invert under the same adjustment.** (a) Crude and adjusted difference in PhenoAge acceleration between adults losing ≥5% of body weight and those who did not; the shaded band is the equivalence margin of one year. (b) The same for serum summed PFAS, with a margin of 5%. Point estimates with 95% confidence intervals.

**Figure 2. Four routes to a pollutant interaction, and where they point.** (a) Forest heterogeneity before and after excluding opioid users and all 167 men on androgen-suppressing drugs. (b) Interaction P values for the Hepatic Steatosis Index and Fatty Liver Index before and after removing their body-mass terms, on a logarithmic scale. (c) The direction of change produced by each route, showing that harm and benefit both occur.

## References

1. Chevrier J, Dewailly É, Ayotte P, Mauriège P, Després J-P, Tremblay A. Body weight loss increases plasma and adipose tissue concentrations of potentially toxic pollutants in obese individuals. Int J Obes Relat Metab Disord. 2000;24(10):1272–1278.
2. Levine ME, Lu AT, Quach A, et al. An epigenetic biomarker of aging for lifespan and healthspan. Aging (Albany NY). 2018;10(4):573–591.
3. Liu Z, Kuo P-L, Horvath S, Crimmins E, Ferrucci L, Levine M. A new aging measure captures morbidity and mortality risk across diverse subpopulations from NHANES IV: a cohort study. PLoS Med. 2018;15(12):e1002718.
4. von Elm E, Altman DG, Egger M, Pocock SJ, Gøtzsche PC, Vandenbroucke JP. The Strengthening the Reporting of Observational Studies in Epidemiology (STROBE) statement: guidelines for reporting observational studies. Lancet. 2007;370(9596):1453–1457.
5. Selva DM, Hogeveen KN, Innis SM, Hammond GL. Monosaccharide-induced lipogenesis regulates the human hepatic sex hormone-binding globulin gene. J Clin Invest. 2007;117(12):3979–3987.
6. Lewis RC, Johns LE, Meeker JD. Serum biomarkers of exposure to perfluoroalkyl substances in relation to serum testosterone and measures of thyroid function among adults and adolescents from NHANES 2011–2012. Int J Environ Res Public Health. 2015;12(6):6098–6114.
7. Xie X, Weng X, Liu S, et al. Perfluoroalkyl and polyfluoroalkyl substance exposure and association with sex hormone concentrations: results from the NHANES 2015–2016. Environ Sci Eur. 2021;33(1):69.
8. Wang Y, Aimuzi R, Nian M, Zhang Y, Luo K, Zhang J. Perfluoroalkyl substances and sex hormones in postmenopausal women: NHANES 2013–2016. Environ Int. 2021;149:106408.
9. Li X, Chen M, Zheng Q, Wang Z, Lin D, Peng M. The association between brominated flame retardants and serum testosterone levels in American adult men: NHANES 2013–2016. Front Public Health. 2025;13:1589047.
10. Madrigal JM, Sargis RM, Persky V, Turyk ME. Multiple organochlorine pesticide exposures and measures of sex steroid hormones in adult males: cross-sectional findings from the 1999–2004 National Health and Nutrition Examination Survey. Int J Hyg Environ Health. 2021;231:113609.
11. Bedogni G, Bellentani S, Miglioli L, et al. The Fatty Liver Index: a simple and accurate predictor of hepatic steatosis in the general population. BMC Gastroenterol. 2006;6:33.
12. Lee JH, Kim D, Kim HJ, et al. Hepatic steatosis index: a simple screening tool reflecting nonalcoholic fatty liver disease. Dig Liver Dis. 2010;42(7):503–508.
13. Vermeulen A, Verdonck L, Kaufman JM. A critical evaluation of simple methods for the estimation of free testosterone in serum. J Clin Endocrinol Metab. 1999;84(10):3666–3672.
14. Wager S, Athey S. Estimation and inference of heterogeneous treatment effects using random forests. J Am Stat Assoc. 2018;113(523):1228–1242.
15. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. Ann Intern Med. 2017;167(4):268–274.
16. Keil AP, Buckley JP, O'Brien KM, Ferguson KK, Zhao S, White AJ. A quantile-based g-computation approach to addressing the effects of exposure mixtures. Environ Health Perspect. 2020;128(4):047004.
17. Grossmann M, Ng Tang Fui M, Dupuis P. Lowered testosterone in male obesity: mechanisms, morbidity and management. Asian J Androl. 2014;16(2):223–231.
18. Corona G, Rastrelli G, Monami M, et al. Body weight loss reverts obesity-associated hypogonadotropic hypogonadism: a systematic review and meta-analysis. Eur J Endocrinol. 2013;168(6):829–843.
