# Supplementary Material

**Three artefacts manufacture a "weight-loss toxin rebound" alarm, and total testosterone masks a fall in the free fraction**

---

## Note on scope

This Brief Communication reports the analyses needed to attribute the "toxin rebound" alarm to its sources, and to describe the androgen read-out that remains. Three bodies of work are reported here in full but barely appear in the main text, because they do not bear on the attribution question:

- the pre-specified design deviations forced by data availability (S1);
- the adipose single-nucleus RNA-sequencing arm (S5);
- the structural comparison of SHBG and albumin (S4, one sentence in the main text).

They are included so that the full evidence base is auditable, not because they support the main claim.

---

## S1. Design deviations from the pre-specified protocol

| Pre-specified assumption | What the data showed | Consequence |
|---|---|---|
| NHANES carries a DNA-methylation clock | No DNAm in any cycle; PhenoAge is a composite of nine clinical biomarkers | Biological age is a constructed score, validated against mortality (S6) |
| Weight history is measured longitudinally | Weight one year earlier is self-reported; no examination-based prior weight | Exposure is a one-year difference, not an incident trajectory |
| Single-cell data carry exposure measurements | GSE295708 has no PFAS measurement and no weight-loss metadata | The single-cell arm cannot be stratified by exposure; see S5 |

These deviations were documented before analysis and are listed here rather than in the main text, which is restricted to 1,500 words.

---

## S2. Statistical audit of the interaction models

The version of this analysis that preceded the present one reported significant pollutant-by-weight-change interactions for five outcomes: white cell count, neutrophils, the systemic immune-inflammation index, monocyte percentage and the platelet-to-lymphocyte ratio. Three of the five did not survive re-analysis.

**What was wrong.** The interaction was specified as `w5*pfas + dw5*pfas`, where `dw5` is the change in treatment status. In the control group `dw5` is identically zero, so the term is zero-inflated and collinear with `w5`; the design matrix additionally dropped the intercept. The combination produces interaction estimates that are not identified by the data but by the parameterisation.

**What was done.** Three non-collinear specifications were fitted — a tertile difference-in-differences, a continuous term per SD of pollutant, and a dose–response term (per 5% weight loss × per SD) — and inference was moved to coefficient-level randomisation: pollutant values were permuted within cycle-by-sex strata, so that the null distribution never passes through the variance estimator. One thousand permutations were run for each of the 30 outcome–specification cells.

**What happened.**

| Stage | Tests reaching P<0.05 |
|---|---|
| Original specification | 5 of 30 |
| Coefficient-level randomisation | 3 of 30, all from two outcomes (white cell count P=0.011; neutrophils 0.031–0.066) |
| Survey-weighted re-estimation | 0 of 30 (white cell count P=0.072) |

The three retracted interactions (systemic immune-inflammation index, monocyte percentage, platelet-to-lymphocyte ratio) are reported in Table S2 with their original and revised estimates. The effect size of the surviving signal is roughly one order of magnitude smaller than originally reported.

---

## S3. Survey weighting

All main effects were re-estimated with the NHANES primary sampling unit, strata and examination weights (`svydesign` + `svyglm`). Design effects were 3.0 (Fatty Liver Index), 3.1 (white cell count), 3.8 (C-reactive protein) and 4.6 (SHBG), so the effective sample was 22–33% of the nominal N (979 to 4,452 participants).

Across 27 main-effect estimates the median absolute difference between weighted and unweighted estimates was 0.12 percentage points and the maximum was 0.75 points. No estimate changed direction. One estimate crossed the significance threshold (bilirubin, P=0.061 to 0.002), a coefficient close to zero on a scale where a percentage change is unstable. Ratios of weighted to unweighted estimates ranged from 0.41 to 2.22, but the extremes occurred only where the estimate was near zero; ratios are therefore not a useful summary for these quantities and absolute differences are reported instead.

For the interaction of interest, weighting attenuated the white cell count interaction from P=0.011 to P=0.072, with an E-value of 1.15 (Table S3). It is therefore described as hypothesis-generating.

---

## S4. Structural comparison of the two carrier proteins

The main text states one sentence. The supporting analysis is summarised here.

**Structures.** Human SHBG with 5α-dihydrotestosterone (PDB 1D2S, 1.95 Å) and human serum albumin with stearate (PDB 1E7I). Ligand geometries for PFOA, PFOS, PFHxS and PFBS were taken from PubChem 3D conformers.

**Steric fit.** The narrowest usable radius of the SHBG steroid pocket, measured over 20 slices along the long axis, is 4.17 Å. The maximum cross-sectional radii of the ligands are 3.90 Å (dihydrotestosterone), 3.64 Å (PFOA) and 3.59 Å (PFOS): all fit. The van der Waals volume of PFOA (221 Å³) is smaller than that of dihydrotestosterone (237 Å³). The steric objection to occupancy therefore fails; the obstacle lies elsewhere.

**Electrostatics.** Within 4.5 Å of the bound steroid, the SHBG pocket contains one charged side chain (Asp65, net −1), which makes a bidentate contact with the 17β-hydroxyl of dihydrotestosterone (2.72 and 2.73 Å) and with Asn82 (2.93 Å). The fatty-acid sites of albumin (FA1–FA7) contain 16 arginine and lysine side chains against 2 acidic residues, a net charge of +14. A perfluoroalkyl acid carries −1 at physiological pH. Occupying the steroid site would therefore place a carboxylate against Asp65, giving a repulsive term of +3.1 to +12.3 kcal/mol across dielectric constants of 40 to 10.

**Occupancy.** With measured concentrations (median PFOA 5.0 nM, PFOS 13.6 nM; SHBG 46 nM monomer, approximately 23 nM dimeric steroid sites; albumin 632 μM), occupancy of the SHBG site is 1.6% for PFOA and 4.3% for PFOS at a dissociation constant of 3 nM, 0.50% and 1.3% at 10 nM, and 0.05% and 0.14% at 100 nM. Reaching 4% for PFOA would require a dissociation constant of about 1 nM, which is at or beyond the affinity of testosterone itself and is not attainable for a same-charge anion.

**What this does and does not show.** It shows that the observed population null is the expected result, not a power failure. It does not measure binding; no PFAS–SHBG complex structure or affinity has been determined here.

---

## S5. Adipose single-nucleus RNA sequencing

Reported in full because it belongs to the same project; not summarised in the main text because it cannot be stratified by exposure.

**Data.** GSE295708, human abdominal subcutaneous adipose tissue, 18 pools (6 with obesity, 6 after weight loss, 6 lean), 186,221 nuclei after quality control.

**Composition.** Across seven cell classes, only the macrophage fraction fell significantly: 25.3% with obesity to 13.5% after weight loss (−46.5%, Welch P=0.010, BH-adjusted P=0.035, Cohen's d=−2.07). The post-weight-loss value is close to the lean value (14.4%). Five of the remaining six classes were not significant; the adipocyte fraction rose (+40.7%, BH-adjusted P=0.019) but reflects a denominator effect and should not be read as adipocyte hyperplasia.

**Function.** Per-cell inflammatory score rose (0.133 to 0.291, P=0.0003), but transcript-level output fell: IL6 −86%, NLRP3 −84%, and the IL6→IL6R ligand–receptor pair −66.8% (q=0.0074). Net inflammatory burden (fraction × per-cell score) was unchanged (0.032 vs 0.039, P=0.19). In an independent paired bulk cohort (GSE199063, Roux-en-Y gastric bypass, 2 and 5 years) the macrophage signature fell by 1.18 and 1.58 SD, while net burden fell by only 0.70 and 0.94 SD — 41.0% and 40.2% less than the change in the macrophage signature.

**Why it is not in the main text.** The dataset carries no exposure measurement, so it cannot show what high-pollutant adipose tissue looks like. It can only exclude adipose inflammation as a necessary mediator, and that exclusion is not needed for the attribution argument in the main text.

---

## S6. Verification and data integrity

- Each reported estimate is linked to a source file in a 256-row ledger (`verification/ledger.csv`), with the generating script in a separate column.
- All 200 analysed NHANES files were compared byte-for-byte with the CDC distribution: 200/200 identical.
- Two genuine errors were found and corrected during the audit: a CRP unit conversion, and the dropped intercept described in S2. Both are documented with before/after estimates.
- PhenoAge was validated against mortality in the full adult sample (RR 1.032 per year; AUC 0.857).
- The imputation of values below the limit of detection is documented congener by congener in S7.

---

## S7. Sensitivity analyses

**Limit of detection.** PFAS values below the limit were replaced with LOD/√2. Re-estimating under LOD (upper bound), LOD/2 (lower bound) and a restricted four-congener panel changed the coefficient by −1.1% to +0.9%; using PFOA alone, the least censored congener (≤1% below LOD), changed it by −5.3%. Percentage below LOD by cycle and congener, and the substitution values used, are listed in Table S7.

**Exposure mixtures.** Quantile g-computation with four to six congeners gave mixture weights dominated by PFNA and PFHxS for white cell count and by PFNA for C-reactive protein; the mixture interaction terms were null for the androgen outcomes (P=0.26 and 0.67).

**Seeds and specification.** The original androgen heterogeneity result (P=0.0049) drifted between P=0.012 and P=0.091 across five random seeds in the all-men sample (2,244), and between 0.383 and 0.670 after drug exclusion (2,077). The same code applied to PhenoAge in 9,656 participants returned P=1.2×10⁻⁷ to 1.0×10⁻⁴ across three seeds, so the instability is specific to the small androgen sample, not to the estimator.

**Power.** In the most sensitive subgroup (n=56) the minimum detectable effect was ±35.7%; point estimates smaller than this are not interpreted. All subgroup estimates are reported in Table S7 regardless of significance.

---

## S8. Directionality audit of the artefacts

Table S8 lists each artefact with its uncontrolled and controlled estimates, the size of the shift, and the direction of the shift. All four — the three exposure-related artefacts and the separate specification error — shifted the estimate towards harm. The arithmetic of the "alarm yield" is: 5 of 30 estimates significant under the original specification, 3 of 30 after coefficient-level randomisation, 0 of 30 after survey weighting.

---

## S9. Literature audit: prescription screening in NHANES pollutant–hormone studies

PubMed was searched on 21 November 2026 (Title/Abstract) for NHANES studies of persistent organic pollutants with a sex-hormone outcome; 19 records were returned and 12 met the inclusion criteria (PFAS 10, polybrominated diphenyl ethers 1, organochlorine pesticides 1). Full texts were checked for the four studies that modelled testosterone in adult men; the remaining eight were assessed from abstracts, and this limitation is stated.

| Screening for androgen-suppressing medication | Studies |
|---|---|
| Excluded named classes (antiandrogens, 5α-reductase inhibitors, adrenal corticosteroids) | 1 (Madrigal 2021) |
| Excluded "hormone therapy" without naming classes | 1 (Li 2025, n=34) |
| Adjusted for contraceptives in women only | 1 (Xie 2021) |
| No medication handling reported | 9 |
| **Screened for GnRH agonists** | **0** |
| **Screened for opioid analgesics** | **0** |

This is a scoping audit, not a systematic review: only PubMed was searched, no PRISMA flow was applied, and eight studies were judged from abstracts.

---

## Supplementary tables (separate files)

| Table | File | Content |
|---|---|---|
| S1 | `TableS1_cohort_flow.csv` | Cohort flow and analytic domains |
| S2 | `TableS2_statistical_audit.csv` | 30 interaction tests, original vs randomised vs weighted |
| S3 | `TableS3_survey_weighted.csv` | Weighted vs unweighted estimates, design effects |
| S7 | `TableS7_sensitivity.csv` | LOD substitution, mixtures, seeds, power |
| S8 | `TableS8_artifact_audit.csv` | Directionality audit |
| S9 | `TableS9_literature_audit.csv` | The 12 studies and their medication handling |
