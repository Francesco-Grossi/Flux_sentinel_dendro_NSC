# Phenology - carbon pipeline: results

Generated on 2026-10-01 by `Script/42_results_report.py` from the outputs of the last pipeline run.

Effects are days of shift in the end of season (EOS) per +1 SD of the predictor unless stated otherwise; negative = earlier senescence. EOS90 / EOS50 / EOS10 = day of year when greenness has fallen to 90 / 50 / 10 % of its seasonal amplitude (onset, middle, end of senescence).

## 1. Data

- Flux sites passing quality control: **129**, with **2177** site-years (1991-2026); 2177 of 2597 site-years passed the growing-season QC.
- Sites by vegetation type (IGBP): ENF 37, GRA 31, DBF 24, OSH 12, MF 7, CSH 5, DNF 4, WSA 3, SAV 3, SNO 1, BSV 1, EBF 1.
- Satellite (HLS): 104,978 image records at 141 sites ({'S30': 66330, 'L30': 38648}); median 51 clear images per site-year.
- Season type: {'summer-green': 101, 'dry-summer': 19} (dry-summer = seasonal GPP peaks before 1 June).

## 2. End-of-season dates from four independent sources

Steps 23-25 fit the same curve to satellite NDVI / NIRv, tower broadband NDVI and PhenoCam greenness. Tables: `data/eos_source_summary.csv`, `data/eos_source_agreement.csv`.

**Coverage, typical dates and year-to-year variation within a site (days)**

| source | n_site_years | n_sites | first_year | last_year | median DOY EOS90 | median DOY EOS50 | median DOY EOS10 | within-site SD EOS90 | within-site SD EOS50 | within-site SD EOS10 |
|---|---|---|---|---|---|---|---|---|---|---|
| PhenoCam (GCC) | 223 | 41 | 2008 | 2025 | 245 | 274 | 300 | 10.7 | 7.3 | 11.4 |
| Tower NDVI | 409 | 58 | 1998 | 2026 | 250 | 286 | 325 | 13.7 | 10.9 | 15.0 |
| Satellite NDVI | 516 | 87 | 2013 | 2026 | 264 | 287 | 311 | 9.3 | 6.5 | 10.0 |
| Satellite NIRv | 770 | 107 | 2013 | 2026 | 226 | 268 | 313 | 10.6 | 7.1 | 10.9 |
| GPP-derived | 1858 | 120 | 1992 | 2025 | - | - | 306 | - | - | 13.6 |

**Do the sources agree on early and late years? Within-site correlation of EOS anomalies**

| pair | EOS90 | EOS50 | EOS10 | n_site_years | n_sites |
|---|---|---|---|---|---|
| PhenoCam (GCC) vs GPP-derived | - | - | 0.25 | 177 | 24 |
| Satellite NDVI vs PhenoCam (GCC) | 0.47 | 0.63 | 0.44 | 92 | 15 |
| Satellite NDVI vs GPP-derived | - | - | 0.27 | 431 | 68 |
| Satellite NDVI vs Tower NDVI | 0.26 | 0.39 | 0.34 | 113 | 21 |
| Satellite NDVI vs Satellite NIRv | 0.56 | 0.66 | 0.43 | 479 | 73 |
| Tower NDVI vs PhenoCam (GCC) | 0.41 | 0.39 | 0.30 | 52 | 9 |
| Tower NDVI vs GPP-derived | - | - | 0.14 | 344 | 40 |
| Satellite NIRv vs PhenoCam (GCC) | 0.23 | 0.34 | 0.19 | 121 | 19 |
| Satellite NIRv vs GPP-derived | - | - | 0.08 | 640 | 87 |
| Satellite NIRv vs Tower NDVI | 0.33 | 0.37 | 0.10 | 148 | 27 |

![Year-to-year agreement between EOS sources, EOS50 (step 31)](../figure/eos_source_agreement/r_within_EOS50.png)

*Year-to-year agreement between EOS sources, EOS50 (step 31)*

![EOS anomalies of each source against PhenoCam (step 31)](../figure/eos_source_agreement/anomaly_scatter_vs_GCC.png)

*EOS anomalies of each source against PhenoCam (step 31)*

![EOS10 of each source against the GPP-derived EOS10 (step 30)](../figure/vi_eos_vs_gpp_eos/EOS10.png)

*EOS10 of each source against the GPP-derived EOS10 (step 30)*

## 3. Are Zani and Lu both right? The split-GPP hypothesis (step 41)

- **H1** cumulative GPP from leaf-out (SOS10) to the summer solstice shifts EOS **earlier** (Zani et al. 2020)
- **H2** cumulative GPP from the solstice to EOS shifts EOS **later**
- **H3** the two cancel, so GPP over the whole season has **no effect** (Lu et al. 2022)

Headline model: strictly within sites (each year minus its site's mean; standard errors clustered by site). 'Controlled' adds the year's leaf-out date and the air temperature before and after the solstice. Table: `data/split_gpp_cancellation_test.csv`; written summary: `Output/split_gpp_cancellation_summary.md`.

**In how many tests does each part hold?** (EOS sources x site groups)

| target | window version | tests | H1 | H2 | H3 | all_three |
|---|---|---|---|---|---|---|
| EOS10 | ends at site mean EOS (unbiased) | 20 | 0 | 3 | 1 | 0 |
| EOS10 | ends at site mean EOS90, before senescence | 20 | 0 | 1 | 8 | 0 |
| EOS10 | ends at same-year EOS (window-length effect built in) | 20 | 0 | 13 | 1 | 0 |
| EOS90 | ends at site mean EOS (unbiased) | 20 | 16 | 7 | 6 | 6 |
| EOS90 | ends at same-year EOS (window-length effect built in) | 20 | 15 | 20 | 0 | 0 |

### 3.1 Onset of senescence (EOS90), windows ending at the site mean EOS90

Days per +1 within-site SD of cumulative GPP. Stars: one-sided test in the direction of the hypothesis (* p<0.05, ** p<0.01, *** p<0.001).

| sites | EOS source | n (sites) | pre-solstice GPP | post-solstice GPP | whole season [95% CI] | all three hold | controlled: pre / post / whole |
|---|---|---|---|---|---|---|---|
| all | PhenoCam (GCC) | 131 (17) | -2.4* | -1.6 | -2.8 [-4.9, -0.8] | no | -2.3 / -2.2 / -2.8 |
| all | Tower NDVI | 303 (35) | -2.4* | +0.5 | -1.6 [-4.0, +0.8] | no | -2.1 / +0.4 / -1.4 |
| all | Satellite NDVI | 362 (56) | -2.5*** | +1.1* | -0.9 [-1.9, +0.1] | yes | -1.8** / +0.8 / -0.5 |
| all | Satellite NIRv | 606 (81) | -2.0*** | +1.1* | -0.8 [-1.7, +0.1] | yes | -2.8*** / +1.2* / -1.1 |
| summer-green | PhenoCam (GCC) | 124 (15) | -2.7* | -2.0 | -3.2 [-5.4, -1.0] | no | -2.5 / -2.5 / -3.3 |
| summer-green | Tower NDVI | 300 (34) | -2.4* | +0.5 | -1.7 [-4.2, +0.8] | no | -2.2 / +0.4 / -1.6 |
| summer-green | Satellite NDVI | 351 (54) | -2.6*** | +1.1* | -1.0 [-2.1, +0.0] | yes | -1.7* / +0.6 / -0.6 |
| summer-green | Satellite NIRv | 564 (76) | -1.9*** | +1.1* | -0.7 [-1.7, +0.3] | yes | -2.7*** / +1.3* / -1.0 |
| dry-summer | Satellite NIRv | 42 (5) | -3.5** | +1.3 | -1.3 [-3.3, +0.7] | no | -4.4*** / +0.9 / -1.7 |
| deciduous | PhenoCam (GCC) | 74 (8) | -1.5 | -3.4 | -3.2 [-6.3, -0.1] | no | +0.1 / -4.2 / -3.3 |
| deciduous | Tower NDVI | 170 (17) | -0.2 | +2.0* | +1.4 [-0.9, +3.7] | no | +0.9 / +1.4 / +1.0 |
| deciduous | Satellite NDVI | 170 (23) | -2.1*** | +1.8* | +0.1 [-1.4, +1.6] | yes | -0.3 / +0.8 / +0.6 |
| deciduous | Satellite NIRv | 196 (25) | -1.5* | +1.6* | +0.1 [-1.8, +2.0] | yes | -1.1 / +1.1 / +0.2 |
| evergreen | Tower NDVI | 79 (12) | -2.0 | -2.0 | -4.1 [-6.0, -2.2] | no | -2.3 / -1.5 / -2.9 |
| evergreen | Satellite NDVI | 76 (14) | -3.3*** | +0.8 | -1.8 [-3.8, +0.2] | no | -2.5* / +0.3 / -1.8 |
| evergreen | Satellite NIRv | 208 (26) | -1.8* | +0.3 | -1.4 [-3.0, +0.2] | no | -2.9*** / +0.6 / -1.9 |
| grass/shrub | PhenoCam (GCC) | 49 (7) | -5.2** | +2.5 | -3.6 [-7.5, +0.4] | no | -6.9*** / +0.8 / -3.8 |
| grass/shrub | Tower NDVI | 54 (6) | -10.1*** | +0.5 | -5.9 [-10.3, -1.5] | no | -9.2** / -2.4 / -8.1 |
| grass/shrub | Satellite NDVI | 116 (19) | -2.6 | -0.3 | -2.0 [-4.3, +0.2] | no | -3.3* / -0.3 / -1.9 |
| grass/shrub | Satellite NIRv | 202 (30) | -2.7** | +1.2 | -1.1 [-2.6, +0.4] | no | -4.4*** / +1.7 / -1.8 |

![EOS90, all sites: pre-solstice, post-solstice and whole-season GPP, for each window version](../figure/split_gpp_cancellation/coefficients_EOS90_all.png)

*EOS90, all sites: pre-solstice, post-solstice and whole-season GPP, for each window version*

![EOS90, deciduous sites](../figure/split_gpp_cancellation/coefficients_EOS90_deciduous.png)

*EOS90, deciduous sites*

![Satellite NDVI: EOS90 anomaly against GPP anomaly, within sites (rows: window versions)](../figure/split_gpp_cancellation/scatter_EOS90_NDVI.png)

*Satellite NDVI: EOS90 anomaly against GPP anomaly, within sites (rows: window versions)*

![Satellite NIRv: EOS90 anomaly against GPP anomaly, within sites (rows: window versions)](../figure/split_gpp_cancellation/scatter_EOS90_NIRv.png)

*Satellite NIRv: EOS90 anomaly against GPP anomaly, within sites (rows: window versions)*

![PhenoCam (GCC): EOS90 anomaly against GPP anomaly, within sites (rows: window versions)](../figure/split_gpp_cancellation/scatter_EOS90_GCC.png)

*PhenoCam (GCC): EOS90 anomaly against GPP anomaly, within sites (rows: window versions)*

![Tower NDVI: EOS90 anomaly against GPP anomaly, within sites (rows: window versions)](../figure/split_gpp_cancellation/scatter_EOS90_NDVI_tower.png)

*Tower NDVI: EOS90 anomaly against GPP anomaly, within sites (rows: window versions)*

### 3.2 End of senescence (EOS10), windows ending at the site mean EOS10

| sites | EOS source | n (sites) | pre-solstice GPP | post-solstice GPP | whole season [95% CI] | all three hold | controlled: pre / post / whole |
|---|---|---|---|---|---|---|---|
| all | PhenoCam (GCC) | 157 (20) | +1.1 | +0.3 | +1.0 [-2.3, +4.2] | no | -0.1 / +1.4 / +1.1 |
| all | Tower NDVI | 313 (37) | -0.8 | +1.1 | +0.8 [-1.0, +2.7] | no | -2.2* / +1.5* / +0.1 |
| all | Satellite NDVI | 380 (59) | +0.7 | +0.5 | +1.4 [-0.1, +2.8] | no | -0.5 / +1.0 / +0.6 |
| all | Satellite NIRv | 624 (85) | +1.6 | +0.6 | +1.7 [+0.3, +3.1] | no | +0.2 / +1.2* / +1.1 |
| summer-green | PhenoCam (GCC) | 131 (16) | +1.4 | +0.1 | +1.0 [-2.7, +4.7] | no | -0.6 / +2.3 / +1.5 |
| summer-green | Tower NDVI | 298 (34) | -0.5 | +0.9 | +0.9 [-1.0, +2.7] | no | -2.1 / +1.4* / +0.0 |
| summer-green | Satellite NDVI | 366 (56) | +0.8 | +0.1 | +1.1 [-0.3, +2.6] | no | -0.3 / +0.7 / +0.4 |
| summer-green | Satellite NIRv | 570 (78) | +1.4 | +0.6 | +1.6 [+0.1, +3.0] | no | -0.3 / +1.4* / +1.0 |
| dry-summer | Satellite NIRv | 54 (7) | +3.7 | +0.4 | +2.9 [-1.3, +7.0] | no | +2.7 / +1.0 / +2.9 |
| deciduous | PhenoCam (GCC) | 81 (9) | +0.3 | +4.1* | +3.8 [-0.4, +8.0] | no | -4.4* / +7.8* / +4.0 |
| deciduous | Tower NDVI | 168 (17) | -1.1 | +1.8* | +1.1 [-1.2, +3.3] | no | -2.9* / +2.7** / +0.7 |
| deciduous | Satellite NDVI | 169 (23) | +0.7 | -0.3 | +0.3 [-1.3, +1.8] | no | -1.0* / +0.9 / +0.1 |
| deciduous | Satellite NIRv | 193 (25) | +0.8 | +0.3 | +0.7 [-1.5, +2.9] | no | -0.9 / +1.5 / +0.6 |
| evergreen | Tower NDVI | 79 (12) | +2.0 | -1.1 | +0.4 [-2.5, +3.3] | no | +1.3 / -0.1 / +0.7 |
| evergreen | Satellite NDVI | 75 (14) | +0.4 | +0.2 | +0.4 [-2.9, +3.8] | no | -1.6 / +1.8 / +0.4 |
| evergreen | Satellite NIRv | 206 (26) | +0.9 | +1.0 | +1.5 [-1.0, +4.0] | no | -0.0 / +1.5* / +1.4 |
| grass/shrub | PhenoCam (GCC) | 68 (9) | +0.5 | -5.7 | -3.1 [-7.1, +0.8] | no | -0.4 / -4.9 / -4.3 |
| grass/shrub | Tower NDVI | 66 (8) | -3.1 | +1.5 | +0.8 [-4.6, +6.2] | no | -4.0 / +0.3 / -3.5 |
| grass/shrub | Satellite NDVI | 136 (22) | +1.0 | +2.7* | +3.7 [+1.4, +6.0] | no | +0.6 / +2.5* / +2.8 |
| grass/shrub | Satellite NIRv | 225 (34) | +3.1 | +0.9 | +3.2 [+1.1, +5.3] | no | +0.8 / +1.1 / +1.7 |

![EOS10, all sites](../figure/split_gpp_cancellation/coefficients_EOS10_all.png)

*EOS10, all sites*

### 3.3 Same-year windows and the window-length effect

When the post-solstice window ends at the same year's EOS, a later EOS makes the window longer and its cumulative GPP larger by construction. The last column is the post-solstice slope produced by window length alone (the site's average GPP curve, no year-specific GPP).

| EOS source | pre | post | whole season | post from window length alone |
|---|---|---|---|---|
| PhenoCam (GCC) | -1.4 | +6.2 | +6.2 | +10.8 |
| Tower NDVI | -2.8 | +10.1 | +6.7 | +13.9 |
| Satellite NDVI | -3.1 | +5.0 | +2.7 | +9.4 |
| Satellite NIRv | -2.7 | +7.9 | +5.0 | +10.4 |

## 4. One-to-one plots: pooled vs within-site (step 40)

The same pairs drawn twice. Pooled plots mix differences between sites with year-to-year changes; the within-site plots keep only the latter and correspond to the tests of section 3. All pairs: `figure/predictor_correlations/` and `figure/predictor_correlations_within_site/`.

**Pooled correlation with EOS90**

| predictor | NDVI | NIRv |
|---|---|---|
| gpp_sos10_to_solstice | 0.22 | -0.30 |
| gpp_solstice_to_eos90 | 0.63 | 0.32 |
| gpp_solstice_to_eos90_fixed | 0.59 | 0.19 |
| total_gpp_growing_season | 0.44 | -0.04 |

![Pooled: EOS90 vs GPP from SOS10 to the solstice](../figure/predictor_correlations/EOS90_vs_gpp_sos10_to_solstice.png)

*Pooled: EOS90 vs GPP from SOS10 to the solstice*

![Within sites: EOS90 vs GPP from SOS10 to the solstice](../figure/predictor_correlations_within_site/EOS90_vs_gpp_sos10_to_solstice.png)

*Within sites: EOS90 vs GPP from SOS10 to the solstice*

![Pooled: EOS90 vs GPP from the solstice to the same-year EOS90](../figure/predictor_correlations/EOS90_vs_gpp_solstice_to_eos90.png)

*Pooled: EOS90 vs GPP from the solstice to the same-year EOS90*

![Within sites: EOS90 vs GPP from the solstice to the same-year EOS90](../figure/predictor_correlations_within_site/EOS90_vs_gpp_solstice_to_eos90.png)

*Within sites: EOS90 vs GPP from the solstice to the same-year EOS90*

![Pooled: EOS90 vs GPP from the solstice to the site mean EOS90](../figure/predictor_correlations/EOS90_vs_gpp_solstice_to_eos90_fixed.png)

*Pooled: EOS90 vs GPP from the solstice to the site mean EOS90*

![Within sites: EOS90 vs GPP from the solstice to the site mean EOS90](../figure/predictor_correlations_within_site/EOS90_vs_gpp_solstice_to_eos90_fixed.png)

*Within sites: EOS90 vs GPP from the solstice to the site mean EOS90*

![Pooled: EOS90 vs GPP over the growing season](../figure/predictor_correlations/EOS90_vs_total_gpp_growing_season.png)

*Pooled: EOS90 vs GPP over the growing season*

![Within sites: EOS90 vs GPP over the growing season](../figure/predictor_correlations_within_site/EOS90_vs_total_gpp_growing_season.png)

*Within sites: EOS90 vs GPP over the growing season*

## 5. Carbon use efficiency and NPP (steps 26, 32)

CUE after Luo et al. (annual), extended to 30-day sliding windows (seasonal / daily). NPP = CUE x GPP. Tables: `data/cue_luo2025_site_year.csv`, `data/cue_seasonal_means_by_site_year.csv`.

- Annual CUE: 2161 site-years at 128 sites, mean 0.52 (SD 0.13).
- Mean of the daily CUE before the solstice 0.62, after 0.59.

**How much do the time scales share? (within-site correlation)**

| a | b | n_site_years | r_within | r_pooled |
|---|---|---|---|---|
| CUE_annual | CUE_pre | 2089 | 0.55 | 0.64 |
| CUE_annual | CUE_post | 2089 | 0.56 | 0.70 |
| CUE_pre | CUE_post | 2092 | 0.43 | 0.64 |
| NPP_pre | NPPd_pre | 2089 | 0.96 | 0.97 |
| NPP_post | NPPd_post | 2089 | 0.92 | 0.96 |

![Seasonal course of CUE by plant type (step 32)](../figure/cue_seasonal/seasonal_curve.png)

*Seasonal course of CUE by plant type (step 32)*

![Annual CUE against pre- and post-solstice CUE, within sites (step 32)](../figure/cue_seasonal/annual_vs_seasonal.png)

*Annual CUE against pre- and post-solstice CUE, within sites (step 32)*

## 6. Which carbon window relates to which EOS? (step 33)

Single-predictor mixed models, windows ending at the site mean EOS (fixed anchors). Table: `data/eos_window_scan.csv`.

**EOS10, window leaf-out to solstice** (days per +1 SD; * p<0.05, ** p<0.01, *** p<0.001, uncorrected)

| predictor | PhenoCam (GCC) | Tower NDVI | Satellite NDVI | Satellite NIRv |
|---|---|---|---|---|
| GPP cum | +1.5 | -0.1 | +2.3 | +4.0** |
| NPP cum | +0.1 | -1.6 | +1.8* | +2.6** |
| NPPd cum | +1.3 | -1.3 | +2.2* | +2.7** |
| GPP mean | -0.0 | -4.1 | -2.0 | +1.3 |
| CUEd mean | -4.1* | -0.0 | -0.8 | +0.4 |

**EOS10, window solstice to mean EOS10** (days per +1 SD; * p<0.05, ** p<0.01, *** p<0.001, uncorrected)

| predictor | PhenoCam (GCC) | Tower NDVI | Satellite NDVI | Satellite NIRv |
|---|---|---|---|---|
| GPP cum | +3.3 | +9.6*** | +6.6*** | +6.2*** |
| NPP cum | -1.7 | +0.2 | -0.4 | +0.5 |
| NPPd cum | +4.6 | +0.8 | +0.6 | +0.8 |
| GPP mean | +1.1 | +7.4** | +5.0** | +4.7*** |
| CUEd mean | -0.1 | -1.0 | -0.6 | -0.8 |

**EOS50, window leaf-out to solstice** (days per +1 SD; * p<0.05, ** p<0.01, *** p<0.001, uncorrected)

| predictor | PhenoCam (GCC) | Tower NDVI | Satellite NDVI | Satellite NIRv |
|---|---|---|---|---|
| GPP cum | -2.5 | -2.2 | -1.3 | -0.9 |
| NPP cum | -0.6 | -0.7 | +0.8 | +2.4*** |
| NPPd cum | +1.0 | -0.1 | +0.8 | +1.7** |
| GPP mean | -0.8 | -3.3 | -0.5 | -0.8 |
| CUEd mean | -1.8 | +0.6 | -0.6 | +0.2 |

**EOS50, window solstice to mean EOS10** (days per +1 SD; * p<0.05, ** p<0.01, *** p<0.001, uncorrected)

| predictor | PhenoCam (GCC) | Tower NDVI | Satellite NDVI | Satellite NIRv |
|---|---|---|---|---|
| GPP cum | +3.3 | +5.1* | +5.1*** | +4.5*** |
| NPP cum | -3.2 | +1.1 | -0.8 | +1.3* |
| NPPd cum | +2.3 | +2.2 | -0.2 | +0.5 |
| GPP mean | +2.3 | +4.3* | +3.7** | +3.7*** |
| CUEd mean | +0.6 | -0.4 | -0.3 | +0.3 |

![PhenoCam (GCC): GPP windows against EOS, fixed anchors (step 33)](../figure/eos_window_scan/fixed_GCC_GPP.png)

*PhenoCam (GCC): GPP windows against EOS, fixed anchors (step 33)*

![PhenoCam (GCC): NPP (daily CUE) windows against EOS, fixed anchors (step 33)](../figure/eos_window_scan/fixed_GCC_NPPd.png)

*PhenoCam (GCC): NPP (daily CUE) windows against EOS, fixed anchors (step 33)*

![Satellite NDVI: GPP windows against EOS, fixed anchors (step 33)](../figure/eos_window_scan/fixed_NDVI_GPP.png)

*Satellite NDVI: GPP windows against EOS, fixed anchors (step 33)*

![Satellite NDVI: NPP (daily CUE) windows against EOS, fixed anchors (step 33)](../figure/eos_window_scan/fixed_NDVI_NPPd.png)

*Satellite NDVI: NPP (daily CUE) windows against EOS, fixed anchors (step 33)*

## 7. When around the solstice is the relation strongest? (step 34)

Mean flux in 15- and 30-day windows starting 75 days before to 75 days after the solstice. Table: `data/eos_solstice_sliding_scan.csv`; the most negative window per case:

| carbon | EOS source | window start (days from solstice) | beta_days_per_sd | p_value | n_obs |
|---|---|---|---|---|---|
| CUEd | PhenoCam (GCC) | -20 | -3.9 | 0.024 | 218 |
| CUEd | Satellite NDVI | -75 | -1.6 | 0.040 | 502 |
| CUEd | Tower NDVI | 5 | -1.7 | 0.168 | 407 |
| CUEd | Satellite NIRv | 45 | -1.0 | 0.128 | 742 |
| GPP | PhenoCam (GCC) | -75 | -0.6 | 0.853 | 204 |
| GPP | Satellite NDVI | -75 | -0.2 | 0.910 | 472 |
| GPP | Tower NDVI | -75 | -4.4 | 0.050 | 378 |
| GPP | Satellite NIRv | 55 | +1.5 | 0.193 | 702 |
| NPP | PhenoCam (GCC) | 55 | -3.1 | 0.382 | 223 |
| NPP | Satellite NDVI | 5 | -0.5 | 0.589 | 515 |
| NPP | Tower NDVI | -5 | -2.6 | 0.228 | 409 |
| NPP | Satellite NIRv | 75 | +0.2 | 0.817 | 758 |
| NPPd | PhenoCam (GCC) | 70 | +0.0 | 0.987 | 218 |
| NPPd | Satellite NDVI | 5 | +0.0 | 0.999 | 502 |
| NPPd | Tower NDVI | 30 | -3.5 | 0.159 | 407 |
| NPPd | Satellite NIRv | 75 | +0.4 | 0.667 | 742 |

![GPP: effect on EOS of 30-day windows by start date relative to the solstice (step 34)](../figure/eos_solstice_sliding_scan/GPP_L30.png)

*GPP: effect on EOS of 30-day windows by start date relative to the solstice (step 34)*

![NPPd: effect on EOS of 30-day windows by start date relative to the solstice (step 34)](../figure/eos_solstice_sliding_scan/NPPd_L30.png)

*NPPd: effect on EOS of 30-day windows by start date relative to the solstice (step 34)*

![CUEd: effect on EOS of 30-day windows by start date relative to the solstice (step 34)](../figure/eos_solstice_sliding_scan/CUEd_L30.png)

*CUEd: effect on EOS of 30-day windows by start date relative to the solstice (step 34)*

## 8. Does carbon explain EOS beyond climate? (step 35)

Nested mixed models on identical rows. M0 climate only; + leaf-out date (SOS); + source (GPP); + sink (NPP). Tables: `data/eos_env_vs_carbon_comparison.csv`, `..._cv.csv`, `..._lrt.csv`.

**EOS10**

| EOS source | model | n_obs | aic | R2 (fixed effects) | gain over climate | R2 leave-one-site-out |
|---|---|---|---|---|---|---|
| PhenoCam (GCC) | M0_env | 156 | 1389 | 0.23 | 0.00 | 0.15 |
| PhenoCam (GCC) | M1_env+SOS | 156 | 1384 | 0.32 | 0.09 | 0.21 |
| PhenoCam (GCC) | M2_env+source | 156 | 1392 | 0.22 | -0.01 | 0.07 |
| PhenoCam (GCC) | M3_env+sink | 156 | 1393 | 0.22 | -0.01 | 0.06 |
| PhenoCam (GCC) | M4_env+source+sink | 156 | 1394 | 0.22 | -0.01 | 0.09 |
| PhenoCam (GCC) | M3d_env+sink_daily | 156 | 1391 | 0.25 | 0.02 | 0.16 |
| PhenoCam (GCC) | M4d_env+source+sink_daily | 156 | 1390 | 0.19 | -0.04 | -0.09 |
| PhenoCam (GCC) | M5_env+SOS+source+sink | 156 | 1389 | 0.35 | 0.12 | 0.25 |
| Tower NDVI | M5_env+SOS+source+sink | 318 | 2874 | 0.19 | 0.12 | -2.74 |
| Tower NDVI | M4d_env+source+sink_daily | 318 | 2877 | 0.16 | 0.09 | -6.82 |
| Tower NDVI | M3d_env+sink_daily | 318 | 2883 | 0.07 | 0.00 | -1.09 |
| Tower NDVI | M4_env+source+sink | 318 | 2877 | 0.16 | 0.09 | -3.38 |
| Tower NDVI | M3_env+sink | 318 | 2883 | 0.07 | 0.00 | -0.42 |
| Tower NDVI | M2_env+source | 318 | 2873 | 0.16 | 0.09 | 0.05 |
| Tower NDVI | M1_env+SOS | 318 | 2881 | 0.08 | 0.01 | -0.13 |
| Tower NDVI | M0_env | 318 | 2880 | 0.07 | 0.00 | -0.11 |
| Satellite NDVI | M5_env+SOS+source+sink | 389 | 3236 | 0.22 | 0.09 | 0.27 |
| Satellite NDVI | M4d_env+source+sink_daily | 389 | 3248 | 0.16 | 0.04 | 0.18 |
| Satellite NDVI | M3d_env+sink_daily | 389 | 3247 | 0.13 | 0.01 | 0.15 |
| Satellite NDVI | M4_env+source+sink | 389 | 3247 | 0.16 | 0.04 | 0.18 |
| Satellite NDVI | M3_env+sink | 389 | 3246 | 0.13 | 0.01 | 0.15 |
| Satellite NDVI | M2_env+source | 389 | 3248 | 0.15 | 0.03 | 0.17 |
| Satellite NDVI | M1_env+SOS | 389 | 3239 | 0.16 | 0.04 | 0.22 |
| Satellite NDVI | M0_env | 389 | 3247 | 0.13 | 0.00 | 0.15 |
| Satellite NIRv | M0_env | 626 | 5272 | 0.07 | 0.00 | -0.07 |
| Satellite NIRv | M1_env+SOS | 626 | 5262 | 0.11 | 0.04 | 0.06 |
| Satellite NIRv | M2_env+source | 626 | 5257 | 0.12 | 0.05 | 0.05 |
| Satellite NIRv | M3_env+sink | 626 | 5275 | 0.08 | 0.00 | -0.07 |
| Satellite NIRv | M4_env+source+sink | 626 | 5260 | 0.12 | 0.05 | 0.05 |
| Satellite NIRv | M3d_env+sink_daily | 626 | 5275 | 0.08 | 0.00 | -0.07 |
| Satellite NIRv | M4d_env+source+sink_daily | 626 | 5260 | 0.12 | 0.05 | 0.05 |
| Satellite NIRv | M5_env+SOS+source+sink | 626 | 5250 | 0.17 | 0.10 | 0.16 |

**EOS50**

| EOS source | model | n_obs | aic | R2 (fixed effects) | gain over climate | R2 leave-one-site-out |
|---|---|---|---|---|---|---|
| PhenoCam (GCC) | M0_env | 156 | 1332 | 0.12 | 0.00 | 0.18 |
| PhenoCam (GCC) | M1_env+SOS | 156 | 1334 | 0.11 | -0.01 | 0.14 |
| PhenoCam (GCC) | M2_env+source | 156 | 1336 | 0.12 | 0.01 | 0.17 |
| PhenoCam (GCC) | M3_env+sink | 156 | 1335 | 0.13 | 0.01 | 0.18 |
| PhenoCam (GCC) | M4_env+source+sink | 156 | 1338 | 0.13 | 0.01 | 0.19 |
| PhenoCam (GCC) | M3d_env+sink_daily | 156 | 1333 | 0.10 | -0.01 | 0.14 |
| PhenoCam (GCC) | M4d_env+source+sink_daily | 156 | 1335 | 0.10 | -0.01 | 0.09 |
| PhenoCam (GCC) | M5_env+SOS+source+sink | 156 | 1338 | 0.12 | 0.00 | 0.18 |
| Tower NDVI | M5_env+SOS+source+sink | 318 | 2692 | 0.09 | 0.06 | -0.26 |
| Tower NDVI | M4d_env+source+sink_daily | 318 | 2690 | 0.09 | 0.05 | -3.78 |
| Tower NDVI | M3d_env+sink_daily | 318 | 2699 | 0.04 | 0.00 | -0.26 |
| Tower NDVI | M4_env+source+sink | 318 | 2690 | 0.09 | 0.05 | -0.28 |
| Tower NDVI | M3_env+sink | 318 | 2699 | 0.04 | 0.00 | -0.57 |
| Tower NDVI | M2_env+source | 318 | 2687 | 0.09 | 0.05 | 0.05 |
| Tower NDVI | M1_env+SOS | 318 | 2697 | 0.03 | -0.00 | -0.17 |
| Tower NDVI | M0_env | 318 | 2695 | 0.04 | 0.00 | -0.15 |
| Satellite NDVI | M5_env+SOS+source+sink | 389 | 2957 | 0.09 | 0.05 | 0.26 |
| Satellite NDVI | M4d_env+source+sink_daily | 389 | 2960 | 0.07 | 0.03 | 0.21 |
| Satellite NDVI | M3d_env+sink_daily | 389 | 2967 | 0.05 | 0.00 | 0.14 |
| Satellite NDVI | M4_env+source+sink | 389 | 2960 | 0.07 | 0.03 | 0.21 |
| Satellite NDVI | M3_env+sink | 389 | 2967 | 0.05 | 0.00 | 0.14 |
| Satellite NDVI | M2_env+source | 389 | 2958 | 0.06 | 0.02 | 0.20 |
| Satellite NDVI | M1_env+SOS | 389 | 2966 | 0.05 | 0.00 | 0.15 |
| Satellite NDVI | M0_env | 389 | 2964 | 0.04 | 0.00 | 0.14 |
| Satellite NIRv | M0_env | 626 | 4931 | 0.02 | 0.00 | -0.07 |
| Satellite NIRv | M1_env+SOS | 626 | 4927 | 0.03 | 0.01 | -0.00 |
| Satellite NIRv | M2_env+source | 626 | 4915 | 0.04 | 0.01 | -0.02 |
| Satellite NIRv | M3_env+sink | 626 | 4934 | 0.02 | -0.00 | -0.07 |
| Satellite NIRv | M4_env+source+sink | 626 | 4918 | 0.04 | 0.01 | -0.02 |
| Satellite NIRv | M3d_env+sink_daily | 626 | 4935 | 0.02 | -0.00 | -0.07 |
| Satellite NIRv | M4d_env+source+sink_daily | 626 | 4918 | 0.04 | 0.01 | -0.02 |
| Satellite NIRv | M5_env+SOS+source+sink | 626 | 4897 | 0.06 | 0.04 | 0.14 |

## 9. Pathways and rate vs cumulative uptake (steps 36, 37)

![PhenoCam (GCC): climate -> GPP -> EOS10 path coefficients (step 36)](../figure/eos_path_analysis/GCC_EOS10_GPP.png)

*PhenoCam (GCC): climate -> GPP -> EOS10 path coefficients (step 36)*

![Satellite NDVI: climate -> GPP -> EOS10 path coefficients (step 36)](../figure/eos_path_analysis/NDVI_EOS10_GPP.png)

*Satellite NDVI: climate -> GPP -> EOS10 path coefficients (step 36)*

**GPP rate (mean) vs cumulative GPP, leaf-out to solstice, all sites** (`data/eos_rate_vs_cumulative.csv`)

| EOS source | target | n_obs | beta_rate | p_rate | beta_cum | p_cum |
|---|---|---|---|---|---|---|
| PhenoCam (GCC) | EOS10 | 177 | -0.0 | 0.999 | +1.5 | 0.639 |
| Tower NDVI | EOS10 | 347 | -4.1 | 0.069 | -0.1 | 0.969 |
| Satellite NDVI | EOS10 | 418 | -2.0 | 0.183 | +2.3 | 0.154 |
| Satellite NIRv | EOS10 | 671 | +1.3 | 0.332 | +4.0 | 0.003 |

## 10. Results by plant type, with multiple-testing correction (step 38)

Tables: `data/eos_results_by_leaf_habit.csv`, `data/eos_results_consistency.csv`.

- 780 tests; 130 (16.7%) with uncorrected p < 0.05 (about 5% expected by chance); 8 significant after FDR correction.

**Effects with the same sign in every EOS source, both ground sources agreeing, at least one significant after FDR**

| leaf_habit | target | model | window | predictor | direction | median_beta | sources |
|---|---|---|---|---|---|---|---|
| all | EOS10 | adj_T | post | GPP cumulative (source) | later EOS | +5.1 | GCC:+1.8, NDVI:+4.4, NDVI_tower:+8.2, NIRv:+5.8 |
| all | EOS50 | adj_T | post | GPP cumulative (source) | later EOS | +3.6 | GCC:+3.1, NDVI:+3.0, NDVI_tower:+4.2, NIRv:+4.1 |
| all | EOS50 | adj_T | post | GPP rate (mean) | later EOS | +3.0 | GCC:+2.2, NDVI:+2.6, NDVI_tower:+3.6, NIRv:+3.5 |

![EOS10: every predictor by plant type and EOS source (step 38)](../figure/eos_results_by_leaf_habit/EOS10_raw.png)

*EOS10: every predictor by plant type and EOS source (step 38)*

![EOS50, air temperature controlled (step 38)](../figure/eos_results_by_leaf_habit/EOS50_adj_T.png)

*EOS50, air temperature controlled (step 38)*

## 11. Timing of anomalies and drought years (step 39)

Tables: `data/anomaly_timing_effects.csv`, `data/drought_years.csv`, `data/drought_eos_contrast.csv`.

- Drought years (May-September water balance at least 1 SD below the site mean): 327 of 2172 site-years (15%).

**Mean EOS shift in drought years against the other years**

| EOS source | target | sites | drought / other site-years | EOS shift in drought years (days) | p_value |
|---|---|---|---|---|---|
| PhenoCam (GCC) | EOS10 | all | 34 / 173 | -1.6 | 0.522 |
| PhenoCam (GCC) | EOS10 | deciduous | 16 / 77 | -3.2 | 0.158 |
| PhenoCam (GCC) | EOS50 | all | 34 / 173 | -2.5 | 0.148 |
| PhenoCam (GCC) | EOS50 | deciduous | 16 / 77 | -0.9 | 0.452 |
| Tower NDVI | EOS50 | deciduous | 25 / 163 | -0.4 | 0.803 |
| Tower NDVI | EOS50 | all | 48 / 335 | -1.0 | 0.510 |
| Tower NDVI | EOS10 | deciduous | 25 / 163 | +0.3 | 0.897 |
| Tower NDVI | EOS10 | all | 48 / 335 | -0.5 | 0.777 |
| Satellite NDVI | EOS50 | deciduous | 33 / 160 | +1.0 | 0.308 |
| Satellite NDVI | EOS50 | all | 78 / 418 | -4.0 | 0.002 |
| Satellite NDVI | EOS10 | deciduous | 33 / 160 | -0.2 | 0.930 |
| Satellite NDVI | EOS10 | all | 78 / 418 | -4.8 | 0.005 |
| Satellite NIRv | EOS10 | all | 123 / 625 | +0.4 | 0.782 |
| Satellite NIRv | EOS10 | deciduous | 40 / 182 | +4.1 | 0.058 |
| Satellite NIRv | EOS50 | all | 123 / 625 | -2.3 | 0.040 |
| Satellite NIRv | EOS50 | deciduous | 40 / 182 | +1.4 | 0.345 |

- Anomaly windows: 3964 tests, 11.4% with uncorrected p < 0.05; after FDR: all years 10, drought years 9, non-drought years 12.

![Satellite NDVI: effect on EOS10 of positive (+) and negative (-) anomalies by time of year (step 39)](../figure/anomaly_timing/NDVI_EOS10.png)

*Satellite NDVI: effect on EOS10 of positive (+) and negative (-) anomalies by time of year (step 39)*

![PhenoCam (GCC): effect on EOS10 of positive (+) and negative (-) anomalies by time of year (step 39)](../figure/anomaly_timing/GCC_EOS10.png)

*PhenoCam (GCC): effect on EOS10 of positive (+) and negative (-) anomalies by time of year (step 39)*

---

Per-step logs are in `logs/`. Method notes are in `README.md` and in the docstring of each script.

