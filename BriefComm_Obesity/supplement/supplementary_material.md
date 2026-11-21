# Supplementary Material

**Two central predictions of the weight-loss toxin-rebound narrative fail in NHANES 2003–2023**

All estimates below were recomputed from the derived datasets in this project. Scripts are in `code/revision/R1`–`R7`; the verification ledger links each number to its source file.

---

## Table S1. Sample flow

| Stage | n | Note |
|---|---|---|
| All records, eight cycles 2003–2023 | 82,096 | 2007–2008 was not fielded |
| Non-pregnant adults aged ≥20 years | 45,739 | sampling frame, not an analytic sample |
| With serum PFAS (≥6 congeners) | 16,778 | of whom 13,162 are adults |
| PFAS domain: PFAS and one-year weight change | 12,806 | 2,528 losing ≥5% |
| PhenoAge domain: complete biomarkers | 9,656 | 1,877 losing ≥5% |
| Androgen domain: men, 2013–2016, complete panel | 1,232 | 245 losing ≥5%; SHBG not measured in 2011–2012 |
| Adipose tissue (GSE199063) | 18 libraries | obese, weight-loss and lean, six each |

The 45,739 figure is the sampling frame and is not used as a denominator for any estimate.

---

## S2. Modelling details and equivalence margins

Both outcomes were fitted with linear models and cluster-robust variance at the stratum × PSU level, using the small-sample correction G/(G−1) × (n−1)/(n−k). Covariates: age, sex, family income-to-poverty ratio, BMI one year earlier, log cotinine, alcohol, physical activity, comorbidity count, total cholesterol, HDL, eGFR, survey year, race, education and smoking.

Prior BMI is on the path from weight change to current BMI: current BMI = prior BMI × (1 − loss). Models are therefore reported with and without it, and a model conditioning on current BMI is shown only for transparency.

Equivalence was assessed by two one-sided tests. For an estimate θ̂ with standard error s and margin m:

p_TOST = max{ Φ((θ̂ − m)/s), 1 − Φ((θ̂ + m)/s) }

| Outcome | Estimate | SE | Margin | 90% CI | p_TOST | Equivalent |
|---|---|---|---|---|---|---|
| PhenoAge acceleration, years | −0.180 | 0.206 | 0.5 | −0.519 to 0.159 | 0.060 | No |
| PhenoAge acceleration, years | −0.180 | 0.206 | 1.0 | −0.519 to 0.159 | 3.5×10⁻⁵ | Yes |
| Serum ΣPFAS, % | −1.24 | 1.50 | 5 | −3.70 to 1.22 | 0.0060 | Yes |
| Serum ΣPFAS, % | −1.24 | 1.50 | 10 | −3.70 to 1.22 | 2.4×10⁻⁹ | Yes |

The 0.5-year margin is not excluded (P=0.060); the paper therefore claims exclusion of harm beyond one year, not beyond half a year.

---

## Table S3. The imbalance, and how much of each estimate it accounts for

| Characteristic | No loss ≥5% (n=7,779) | Loss ≥5% (n=1,877) |
|---|---|---|
| Age, years | 51.2 (17.9) | 52.1 (18.2) |
| Female, % | 50.3 | 56.3 |
| NH Black, % | 18.2 | 21.4 |
| College or above, % | 26.3 | 19.6 |
| Family income-to-poverty ratio | 2.7 (1.6) | 2.4 (1.5) |
| BMI one year earlier, kg/m² | 27.7 (6.2) | 32.4 (7.8) |
| Current BMI, kg/m² | 29.1 (7.0) | 29.8 (6.9) |
| Diabetes, % | 10.9 | 19.6 |
| Hypertension, % | 34.1 | 44.5 |
| Current smoker, % | 18.1 | 22.0 |
| PhenoAge acceleration, years | −0.5 (7.6) | 1.4 (9.0) |
| Sum PFAS, nmol/L, median (IQR) | 28.7 (14.8–53.3) | 25.0 (12.7–48.1) |

Stepwise adjustment, serum ΣPFAS (percentage difference per ≥5% weight loss):

| Model | Estimate, % | SE | P |
|---|---|---|---|
| Unadjusted | −10.90 | 2.17 | 1×10⁻⁷ |
| + race | −11.54 | 2.16 | 1.4×10⁻⁸ |
| + education | −11.35 | 2.14 | 1.9×10⁻⁸ |
| + smoking | −11.33 | 2.16 | 2.7×10⁻⁸ |
| + age | −11.74 | 2.12 | 3.6×10⁻⁹ |
| + sex | −8.76 | 2.12 | 1.6×10⁻⁵ |
| + income-to-poverty ratio | −9.93 | 2.15 | 1.2×10⁻⁶ |
| + BMI one year earlier | −5.83 | 2.07 | 0.0037 |
| Full covariate set | −1.24 | 1.50 | 0.41 |

PhenoAge acceleration, alternative specifications:

| Specification | Estimate, years | SE | P |
|---|---|---|---|
| Crude | +1.877 | 0.230 | 3×10⁻¹⁶ |
| Adjusted, year linear | −0.180 | 0.206 | 0.382 |
| Adjusted, cycle fixed effects | −0.197 | 0.204 | 0.333 |
| AIPW forest, year linear | −0.329 | 0.200 | 0.100 |
| AIPW forest, cycle fixed effects | −0.266 | 0.199 | 0.182 |

Direction reversal is stable across specifications; the adjusted estimate is not distinguishable from zero in any of them.

Dose gradient (reference: weight stable within ±5%):

| Band | n | Change in ΣPFAS, % (95% CI) | P |
|---|---|---|---|
| Stable ±5% | 7,853 | reference | — |
| Loss 5–10% | 1,448 | +0.5 (−3.2 to 4.3) | 0.80 |
| Loss ≥10% | 1,080 | −5.9 (−9.8 to −1.8) | 0.005 |
| Gain ≥5% | 2,425 | −3.9 (−6.6 to −1.2) | 0.006 |

The gradient is not monotonic, so it is not presented as a dose–response.

---

## Table S4. Composition of the 167 men on androgen-suppressing drugs

| Class | All (n=2,244) | Losing ≥5% (n=439) | % of all |
|---|---|---|---|
| Opioid analgesic | 131 | 37 | 5.84 |
| 5α-reductase inhibitor | 54 | 7 | 2.41 |
| Systemic glucocorticoid | 39 | 8 | 1.74 |
| GnRH agonist (ADT) | 1 | 1 | 0.04 |
| Antiandrogen | 2 | 0 | 0.09 |
| Ketoconazole | 3 | 2 | 0.13 |
| Oestrogen | 0 | 0 | 0 |
| Any of the above | 167 | 47 | 7.44 |

Forest heterogeneity, one-sided P, 4,000 trees:

| Exclusion | n | P (one-sided) | P (two-sided) |
|---|---|---|---|
| None | 2,243 | 0.0555 | 0.111 |
| Opioid users | 2,113 | 0.5480 | 1.000 |
| Glucocorticoid users | 2,203 | 0.1676 | 0.335 |
| 5α-reductase inhibitor users | 2,183 | 0.0688 | 0.138 |
| All 167 | 2,077 | 0.6728 | 1.000 |

Parametric main effect per 5% weight loss, and the joint 2-df interaction:

| Exclusion | n | Main effect, % | P | Interaction P |
|---|---|---|---|---|
| None | 2,244 | +4.28 | 1.9×10⁻⁵ | 0.671 |
| Opioid users | 2,113 | +4.19 | 1.4×10⁻⁴ | 0.977 |
| Glucocorticoid users | 2,205 | +4.44 | 8.9×10⁻⁷ | 0.720 |
| 5α-reductase inhibitor users | 2,190 | +4.31 | 2.1×10⁻⁵ | 0.663 |
| All 167 | 2,077 | +4.46 | 3.1×10⁻⁶ | 0.617 |
| All seven classes | 2,029 | +4.56 | 2.4×10⁻⁶ | 0.659 |

The main effect is unchanged; only the interaction moves. The pre-exclusion forest P was already 0.055, so the correct description is that an already marginal signal disappeared, not that the 167 carried it.

---

## S5. Formula decomposition

| Index | Interaction P with BMI terms | Interaction P without BMI terms |
|---|---|---|
| Hepatic Steatosis Index (8 × ALT/AST) | 1.2×10⁻⁶ | 0.68 |
| Fatty Liver Index | 0.017 | 0.14 |
| log ALT (no BMI term to begin with) | 0.61 | 0.61 |

Both indices contain BMI as a term. Removing it leaves the liver-specific component, and the interaction disappears. ALT, which never contained BMI, shows none.

---

## S6. The specification error, and why forest P is not reported as evidence

In the first pass of these analyses, the interaction model was written as `y ~ w5*pfas_t + dw5*pfas_t`, where `dw5` is identically zero in the comparison group. This makes the design matrix collinear and zero-inflated; the intercept was also dropped when constant columns were removed, which rescales the interaction. The consequences:

| Quantity | First pass | After correction |
|---|---|---|
| Nominal interactions with P<0.05, of 30 | 5 | 3 |
| Surviving coefficient-level randomisation (1,000 permutations) | — | 3 |
| Systemic immune-inflammation index, weight-loss effect | −8.0% | −1.0% |

The systemic immune-inflammation index shows that the same error attenuated one association while inflating others: the error multiplies the design, it does not push in one direction.

Forest heterogeneity is not used as evidence for three reasons, all reproducible:

| Specification | n | One-sided P |
|---|---|---|
| Archived run (4,000 trees, seed 20261117) | 2,244 | 0.0048 |
| Identical re-run, same seed, design matrix differing in no cell | 2,244 | 0.1030 |
| Drop cycle G (2011–2012) | 1,480 | 0.7019 |
| Drop cycle H (2013–2014) | 1,515 | 0.4785 |
| Drop cycle I (2015–2016) | 1,493 | 0.1884 |
| Drop two highest-Cook-distance men | 2,242 | 0.0926 |

The archived P (0.0048) could not be reproduced with byte-identical inputs, which indicates run-level non-determinism in the tuning step rather than a data change. Dropping one cycle moves the same P from 0.10 to 0.70. The two Cook-distance points do not matter, so the instability is not driven by individuals.

---

## Table S7. Literature audit

PubMed search, November 2026: NHANES AND (per- or polyfluoroalkyl OR organochlorine OR polychlorinated biphenyl) AND (testosterone OR sex hormone OR androgen). Nineteen records; twelve included after screening.

| Study | Cycle | Exposure | Outcome | Screened for androgen-suppressing drugs |
|---|---|---|---|---|
| Lewis 2015 | 2011–2012 | PFAS | Testosterone, thyroid | Not reported |
| Xie 2021 | 2015–2016 | PFAS | Sex hormones | Not reported |
| Wang 2021 | 2013–2016 | PFAS | Sex hormones (women) | Not reported |
| Li 2025 | 2013–2016 | Brominated flame retardants | Testosterone | Not reported |
| Madrigal 2021 | 1999–2004 | Organochlorine pesticides | Sex steroids | Excluded antiandrogens and 5α-reductase inhibitors |
| Eight further NHANES analyses | various | PFAS or organochlorines | Sex hormones | Not reported |

Two of twelve screened for any androgen-suppressing medication; none screened for GnRH agonists or opioid analgesics. Four studies were read in full; the remaining eight were assessed from abstracts only, which is a limitation of this audit.

---

## S8. Tissue-level observation

In GSE199063 (18 libraries: obese, weight-loss, lean; six each), weight loss reduced macrophage abundance from 25.3% to 13.5% of nuclei (−46.5%, Welch P=0.0099, BH-adjusted P=0.035) with a reciprocal rise in adipocyte fraction (+40.7%), consistent with a denominator effect. A senescence module remained above the lean group (difference +0.023, P=0.024), and the inflammatory module did not differ from obese controls (+0.022, P=0.165). These data are descriptive and come from a different population from the NHANES sample.
