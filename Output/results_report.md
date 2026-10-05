# Phenology - carbon pipeline: results

Generated on 2026-10-05 by `Script/46_results_report.py` from the outputs of the last pipeline run. A written interpretation of these results, with the analyses still to do, is in [findings_and_next_steps.md](findings_and_next_steps.md).

Effects are days of shift in the end of season (EOS) per +1 within-site SD of the predictor (within-site models: each year minus its site's mean, standard errors clustered by site) unless stated otherwise; negative = earlier senescence. EOS90 / EOS50 / EOS10 = day of year when greenness has fallen to 90 / 50 / 10 % of its seasonal amplitude (onset, middle, end of senescence).

## Methods: what was done, why, and what makes it reliable

### M1. The question

Zani et al. (2020) found that more photosynthesis early in the season brings senescence forward. Lu et al. (2022) found that productivity over the growing season does not. The hypothesis tested here is that both are right: GPP before the summer solstice advances senescence (H1), GPP after it delays senescence (H2), and the two cancel so that whole-season GPP shows no effect (H3).

### M2. Flux data (steps 11, 21, 22)

- **What:** daily FLUXNET data (GPP and Reco from night-time partitioning, NEE, air temperature, shortwave radiation, VPD, precipitation) for sites north of 30 N with natural vegetation and a long gap-free record.
- **Quality control:** a site-year is dropped when 50% or more of its growing-season days (1 March - 31 October) are low quality or missing for GPP, NEE or Reco, or when 50 or more bad days follow each other.
- **Why:** GPP is measured at the site, every day, independently of any greenness index. Earlier studies of this question mostly used modelled or satellite-derived productivity, which shares its input with the satellite phenology it is compared to.

### M3. End-of-season dates from four sources (steps 12, 13, 23-25)

- **Satellite NDVI and NIRv:** Harmonized Landsat Sentinel-2 (HLS L30 + S30, 30 m), averaged over a 1 km radius around the tower. Only pixels whose ESA WorldCover class matches the site's vegetation type are used, so roads, fields and water inside the radius do not enter. Cloud, shadow and snow pixels are masked; the snow fraction of each image is kept.
- **Tower NDVI:** broadband NDVI from the tower's own incoming and reflected shortwave and PAR sensors, midday records only, measured (not gap-filled) radiation only.
- **PhenoCam:** canopy greenness (GCC, 90th percentile of the day) from the camera at the tower, for the region of interest that matches the site's vegetation type.
- **One fitting routine for all four** (`pheno_fit.py`): a double-logistic curve per site-year. A year needs at least 12 clear observations. Snow-covered and frozen periods are set to the site's dormant-season background (Beck et al. 2006). Gaps longer than 30 days are filled, at low weight, with the site's multi-year curve shifted to the year's level (as in the MSLSP product), so a gap cannot bend the curve but the year's own observations decide the dates.
- **Dates:** EOS90, EOS50 and EOS10 are the days when the fitted curve has fallen to 90, 50 and 10% of its amplitude after the peak; leaf-out dates (SOS10/50/90) are defined the same way on the rising side.
- **Fit quality control:** R2 >= 0.8 on the real observations; amplitude at least 2 times the fit error; dates in the right order; at least 3 observations after the peak; EOS90 within 60 days of the site's own median. The last rule uses the site as its own reference, so dry-summer sites that really senesce in June are kept.
- **Why four sources:** each has a different weakness (satellite: clouds and mixed pixels; tower NDVI: sensor drift; PhenoCam: few sites, one viewing angle). A result that appears in all four is not an artefact of one instrument. Because the same curve and the same definitions are used, the dates can be compared directly.

### M4. Predictors: windows of GPP around the solstice (steps 27, 28, 41)

- **pre** = GPP from leaf-out (SOS10) to the summer solstice; **post** = GPP from the solstice to EOS; **total** = both.
- **Fixed anchors.** Windows end at the site's mean EOS over all years, not at the same year's EOS. Why: a window that ends at the year's own EOS is longer when EOS is later, so its GPP sum rises with EOS by construction. The size of that artefact is measured with a length-only null (the site's average GPP curve summed over the year's window), and it is as large as the apparent effect.
- **Calendar windows** (60 days before, 45 days from the solstice) and the **GPP rate** (mean per day) are used as well. Why: a window that starts at leaf-out is longer in an early spring, so its GPP sum partly measures leaf-out date.

### M5. The statistical model (eos_common.py; steps 33-44)

- **Within-site model.** Every variable is the year's value minus the site's own mean (sites with at least 3 years). Predictors are divided by their within-site standard deviation. Ordinary least squares on these anomalies, with standard errors clustered by site. An effect is the shift of EOS in days when the predictor is one standard deviation above the site's normal.
- **Why not pooled correlations or mixed models.** The hypothesis is about what happens at a site in a productive year. Sites differ in productivity and in senescence date for many reasons (climate, species), and a pooled plot mostly shows those differences. A mixed model with a random site intercept removes them only partly; here it gave effects about twice as large (section 6). Subtracting the site mean removes everything that is constant at a site.
- **Spring temperature as a standing covariate.** The mean air temperature of the 60 days before the solstice is in every model of the main test. Why: a warm spring raises GPP and advances the onset of senescence by itself (section 14), so without it the GPP effect is overstated.
- **H3 as an equivalence test.** 'Not significant' is not evidence of no effect. H3 is accepted only when the whole-season effect lies significantly inside +/-2 days per SD (two one-sided tests).
- **Many tests.** Where many predictor x source x group cells are tested (step 38), p-values are corrected with Benjamini-Hochberg, and an effect counts as robust only if it has the same sign in every EOS source.

### M6. Sink side: carbon use efficiency (step 26)

- CUE per site-year from the method of Luo et al. (two-round MCMC on day-pair differences of Reco and GPP), ported from the authors' MATLAB code, with two indexing errors of that code corrected. Extended here to 30-day sliding windows for a seasonal CUE. NPP = CUE x GPP.
- Limit: NPP and the respiration terms derived from it are GPP multiplied by a factor, so they are not independent of GPP.

### M7. Checks on the main result (steps 42-44)

- **Same site-years** (step 42): PhenoCam and satellite compared on the site-years both have.
- **Leave one site out, other QC thresholds, statistical power** (step 43).
- **Alternative explanations** (step 44): leaf-out date, water balance, spring temperature, and sink variables, each put in the same model as GPP.

### M8. Strong points

1. **Measured GPP** at the tower, independent of the greenness data that give the senescence dates.
2. **Four independent senescence records** processed with one routine and identical definitions.
3. **Within-site inference**: differences between sites cannot produce the result.
4. **The window-length artefact is removed and quantified**, not just mentioned.
5. **'No effect' is tested, not assumed** (equivalence test).
6. **Robustness is shown**: every single-site removal, twelve QC settings, three versions of the predictor.
7. **Competing explanations are tested in the same model** (temperature, leaf-out, water).
8. **Reproducible**: one command reruns everything from the raw downloads; every number in this report is read from the tables of the last run.

### M9. Limits

- Observational data: the models show association within sites, not causation.
- Senescence dates carry an error of about 7-11 days (year-to-year SD within a site), comparable to the signal; effects of 1-2 days per SD need several hundred site-years.
- PhenoCam and tower NDVI have too few site-years to confirm an effect of this size on their own.
- The satellite record starts in 2013; northern temperate and boreal sites only.

## 1. Data

- Flux sites passing quality control: **129**, with **2177** site-years (1991-2026); 2177 of 2597 site-years passed the growing-season QC.
- Sites by vegetation type (IGBP): ENF 37, GRA 31, DBF 24, OSH 12, MF 7, CSH 5, DNF 4, WSA 3, SAV 3, SNO 1, BSV 1, EBF 1.
- Satellite (HLS): 104,978 image records at 141 sites ({'S30': 66330, 'L30': 38648}); median 51 clear images per site-year.
- Season type: {'summer-green': 101, 'dry-summer': 19} (dry-summer = seasonal GPP peaks before 1 June).

**What this shows**

- The sample is 57% forest; the rest is grassland, shrubland and savanna. Results for 'all sites' are therefore mostly about forests, and the non-forest groups are small.
- A site contributes 17 years on average. That is what makes within-site comparisons (an early year against a late year of the same site) possible.
- With a median of 51 clear satellite images per site-year, about one every 7 days, the autumn decline is sampled densely enough to date its onset. Before Sentinel-2 was added the median was about 13.
- 19 sites have their GPP peak before June. There the canopy dries out in early summer, so 'senescence' is a response to water, and those sites are analysed as a separate group.

## 2. End-of-season dates from four independent sources

Steps 23-25 fit the same curve to satellite NDVI / NIRv, tower broadband NDVI and PhenoCam greenness. Tables: `data/eos_source_summary.csv`, `data/eos_source_agreement.csv`.

**Coverage, typical dates and year-to-year variation within a site (days)**

| source | n_site_years | n_sites | first_year | last_year | median DOY EOS90 | median DOY EOS50 | median DOY EOS10 | within-site SD EOS90 | within-site SD EOS50 | within-site SD EOS10 |
|---|---|---|---|---|---|---|---|---|---|---|
| PhenoCam (GCC) | 266 | 41 | 2006 | 2025 | 244 | 270 | 296 | 10.3 | 7.6 | 11.1 |
| Tower NDVI | 409 | 58 | 1998 | 2026 | 250 | 286 | 325 | 13.7 | 10.9 | 15.0 |
| Satellite NDVI | 516 | 87 | 2013 | 2026 | 264 | 287 | 311 | 9.3 | 6.5 | 10.0 |
| Satellite NIRv | 770 | 107 | 2013 | 2026 | 226 | 268 | 313 | 10.6 | 7.1 | 10.9 |
| GPP-derived | 1858 | 120 | 1992 | 2025 | - | - | 306 | - | - | 13.6 |

**Do the sources agree on early and late years? Within-site correlation of EOS anomalies**

| pair | EOS90 | EOS50 | EOS10 | n_site_years | n_sites |
|---|---|---|---|---|---|
| PhenoCam (GCC) vs GPP-derived | - | - | 0.23 | 219 | 28 |
| Satellite NDVI vs PhenoCam (GCC) | 0.45 | 0.62 | 0.32 | 110 | 18 |
| Satellite NDVI vs GPP-derived | - | - | 0.27 | 431 | 68 |
| Satellite NDVI vs Tower NDVI | 0.26 | 0.39 | 0.34 | 113 | 21 |
| Satellite NDVI vs Satellite NIRv | 0.56 | 0.66 | 0.43 | 479 | 73 |
| Tower NDVI vs PhenoCam (GCC) | 0.38 | 0.36 | 0.27 | 67 | 12 |
| Tower NDVI vs GPP-derived | - | - | 0.14 | 344 | 40 |
| Satellite NIRv vs PhenoCam (GCC) | 0.24 | 0.36 | 0.27 | 147 | 24 |
| Satellite NIRv vs GPP-derived | - | - | 0.08 | 640 | 87 |
| Satellite NIRv vs Tower NDVI | 0.33 | 0.37 | 0.10 | 148 | 27 |

![Year-to-year agreement between EOS sources, EOS50 (step 31)](../figure/eos_source_agreement/r_within_EOS50.png)

*Year-to-year agreement between EOS sources, EOS50 (step 31)*

![EOS anomalies of each source against PhenoCam (step 31)](../figure/eos_source_agreement/anomaly_scatter_vs_GCC.png)

*EOS anomalies of each source against PhenoCam (step 31)*

![EOS10 of each source against the GPP-derived EOS10 (step 30)](../figure/vi_eos_vs_gpp_eos/EOS10.png)

*EOS10 of each source against the GPP-derived EOS10 (step 30)*

**What this shows**

- Site-years with a usable EOS90: PhenoCam (GCC) 266 (41 sites), Tower NDVI 409 (58 sites), Satellite NDVI 516 (87 sites), Satellite NIRv 770 (107 sites). The satellite indices have by far the largest samples; the two ground sources cover fewer sites.
- Within a site, EOS90 varies from year to year by 9.3 to 13.7 days (median SD), EOS10 by 10.0 to 15.0 days. This is the size of the signal every later analysis tries to explain: an effect of 1-2 days per SD of a predictor is a small part of it.
- The sources agree only moderately on which years are early or late. The best pair for EOS50 is Satellite NDVI with Satellite NIRv (within-site r = 0.66); for EOS90 the correlations range from 0.24 to 0.56. A correlation of 0.5 between two measures of the same event means that roughly half of the year-to-year variance of each is measurement error or a real difference in what they see. This limits how strong any relation with GPP can appear.
- Agreement is better for the middle of senescence (EOS50) than for its onset (EOS90) or end (EOS10): the start and the tail of the decline are flat parts of the curve and are harder to date.

## 3. Are Zani and Lu both right? The split-GPP hypothesis (step 41)

- **H1** cumulative GPP from leaf-out (SOS10) to the summer solstice shifts EOS **earlier** (Zani et al. 2020)
- **H2** cumulative GPP from the solstice to EOS shifts EOS **later**
- **H3** the two cancel, so GPP over the whole season has **no effect** (Lu et al. 2022)

Headline model: strictly within sites (each year minus its site's mean; standard errors clustered by site), with the air temperature of the 60 days before the solstice held fixed. 'Spring T' is the effect of that temperature; 'without T' is the same model without it; 'controlled' also adds the year's leaf-out date and the air temperature after the solstice. Table: `data/split_gpp_cancellation_test.csv`; written summary: `Output/split_gpp_cancellation_summary.md`.

**In how many tests does each part hold?** (EOS sources x site groups)

| target | window version | tests | H1 | H2 | H3 | all_three |
|---|---|---|---|---|---|---|
| EOS10 | calendar windows around the solstice (independent of leaf-out) | 21 | 1 | 2 | 8 | 0 |
| EOS10 | ends at site mean EOS (unbiased) | 20 | 1 | 2 | 1 | 0 |
| EOS10 | ends at site mean EOS90, before senescence | 20 | 1 | 1 | 9 | 0 |
| EOS10 | mean daily GPP instead of the sum | 20 | 10 | 11 | 5 | 2 |
| EOS10 | ends at same-year EOS (window-length effect built in) | 20 | 2 | 13 | 0 | 0 |
| EOS90 | calendar windows around the solstice (independent of leaf-out) | 21 | 4 | 1 | 5 | 0 |
| EOS90 | ends at site mean EOS (unbiased) | 20 | 12 | 3 | 6 | 3 |
| EOS90 | mean daily GPP instead of the sum | 20 | 4 | 1 | 4 | 0 |
| EOS90 | ends at same-year EOS (window-length effect built in) | 20 | 15 | 20 | 0 | 0 |

### 3.1 Onset of senescence (EOS90), windows ending at the site mean EOS90

Days per +1 within-site SD of cumulative GPP. Stars on GPP: one-sided test in the direction of the hypothesis; on spring T: two-sided (* p<0.05, ** p<0.01, *** p<0.001).

| sites | EOS source | n (sites) | pre-solstice GPP | post-solstice GPP | whole season [95% CI] | all three hold | spring T | without T: pre / post / whole | controlled: pre / post / whole |
|---|---|---|---|---|---|---|---|---|---|
| all | PhenoCam (GCC) | 174 (22) | -0.7 | -1.8 | -2.0 [-3.6, -0.4] | no | -2.6** | -1.3 / -1.0 / -1.8 | -1.2 / -1.8 / -2.3 |
| all | Tower NDVI | 303 (35) | -1.5 | +0.0 | -1.3 [-3.5, +0.9] | no | -2.2 | -2.4* / +0.5 / -1.6 | -1.7 / +0.1 / -2.0 |
| all | Satellite NDVI | 362 (56) | -2.7*** | +1.3* | -1.0 [-1.9, -0.0] | yes | +0.4 | -2.5*** / +1.1* / -0.9 | -2.1** / +1.1 / -1.0 |
| all | Satellite NIRv | 606 (81) | -1.6** | +0.7 | -0.8 [-1.7, +0.1] | no | -1.7** | -2.0*** / +1.1* / -0.8 | -2.6*** / +1.0* / -1.0 |
| summer-green | PhenoCam (GCC) | 167 (20) | -0.7 | -2.0 | -2.1 [-3.8, -0.4] | no | -2.6** | -1.4 / -1.2 / -1.9 | -1.2 / -2.0 / -2.6 |
| summer-green | Tower NDVI | 300 (34) | -1.5 | -0.1 | -1.4 [-3.7, +0.8] | no | -2.2 | -2.4* / +0.5 / -1.7 | -1.7 / +0.1 / -2.1 |
| summer-green | Satellite NDVI | 351 (54) | -3.0*** | +1.4* | -1.1 [-2.1, -0.1] | yes | +0.8 | -2.6*** / +1.1* / -1.0 | -2.2** / +1.1 / -1.1 |
| summer-green | Satellite NIRv | 564 (76) | -1.4** | +0.7 | -0.7 [-1.6, +0.3] | no | -1.9** | -1.9*** / +1.1* / -0.7 | -2.4*** / +1.0* / -0.9 |
| dry-summer | Satellite NIRv | 42 (5) | -3.5** | +1.4 | -1.4 [-3.6, +0.7] | no | +0.2 | -3.5** / +1.3 / -1.3 | -4.3*** / +1.0 / -1.6 |
| deciduous | PhenoCam (GCC) | 98 (8) | -0.3 | -2.5 | -2.0 [-4.2, +0.1] | no | -0.4 | -0.5 / -2.3 / -2.1 | +2.0 / -4.6 / -2.7 |
| deciduous | Tower NDVI | 170 (17) | +0.5 | +1.3 | +1.4 [-0.8, +3.5] | no | -1.9 | -0.2 / +2.0* / +1.4 | +0.9 / +1.4 / +0.9 |
| deciduous | Satellite NDVI | 170 (23) | -2.7*** | +2.4** | +0.1 [-1.5, +1.6] | yes | +1.3 | -2.1*** / +1.8* / +0.1 | -0.6 / +1.2 / +0.7 |
| deciduous | Satellite NIRv | 196 (25) | -0.6 | +0.7 | +0.2 [-1.9, +2.2] | no | -2.2** | -1.5* / +1.6* / +0.1 | -0.7 / +0.7 / +0.0 |
| evergreen | Tower NDVI | 79 (12) | -0.6 | -2.1 | -3.1 [-5.9, -0.3] | no | -3.2 | -2.0 / -2.0 / -4.1 | -1.9 / -1.1 / -3.9 |
| evergreen | Satellite NDVI | 76 (14) | -3.2*** | +0.8 | -1.5 [-3.4, +0.4] | no | -0.2 | -3.3*** / +0.8 / -1.8 | -3.0** / +0.5 / -1.7 |
| evergreen | Satellite NIRv | 208 (26) | -1.6* | +0.1 | -1.3 [-2.8, +0.1] | no | -1.4 | -1.8* / +0.3 / -1.4 | -2.8*** / +0.4 / -1.8 |
| grass/shrub | PhenoCam (GCC) | 65 (11) | -4.0** | +1.2 | -3.1 [-6.2, +0.0] | no | -3.8*** | -4.4** / +2.5* / -2.8 | -5.6*** / +0.9 / -3.7 |
| grass/shrub | Tower NDVI | 54 (6) | -10.4*** | +0.7 | -5.7 [-9.5, -1.9] | no | +0.6 | -10.1*** / +0.5 / -5.9 | -10.0** / -1.6 / -6.5 |
| grass/shrub | Satellite NDVI | 116 (19) | -2.4* | -0.4 | -2.1 [-4.4, +0.2] | no | -0.3 | -2.6 / -0.3 / -2.0 | -3.1* / -0.4 / -2.6 |
| grass/shrub | Satellite NIRv | 202 (30) | -2.7** | +0.9 | -1.4 [-2.8, +0.1] | no | -2.1 | -2.7** / +1.2 / -1.1 | -4.3*** / +1.5 / -1.7 |

**What this shows**

- **Pre-solstice GPP (H1).** With spring temperature held fixed, a year with one SD more GPP between leaf-out and the solstice has an onset of senescence shifted by PhenoCam (GCC) -0.7, Tower NDVI -1.5, Satellite NDVI -2.7***, Satellite NIRv -1.6** days. The sign is negative in 4 of 4 sources and significant in 2 of 4. This is the Zani-type effect, and it is small: the year-to-year SD of EOS90 is 7-10 days.
- **What spring temperature takes away.** Without the temperature covariate the same slopes are PhenoCam (GCC) -1.3, Tower NDVI -2.4*, Satellite NDVI -2.5***, Satellite NIRv -2.0***. Spring temperature itself shifts EOS90 by PhenoCam (GCC) -2.6**, Tower NDVI -2.2, Satellite NDVI +0.4, Satellite NIRv -1.7** days per SD. A warm spring raises GPP and advances senescence on its own, so part of what looks like a GPP effect is a temperature effect.
- **Post-solstice GPP (H2).** PhenoCam (GCC) -1.8, Tower NDVI +0.0, Satellite NDVI +1.3*, Satellite NIRv +0.7 days per SD; the hypothesis needs a positive effect and finds one in 1 of 4 sources. There is no consistent delaying effect of late-season GPP.
- **Whole-season GPP (H3).** PhenoCam (GCC) -2.0, Tower NDVI -1.3, Satellite NDVI -1.0, Satellite NIRv -0.8 days per SD. The 95% interval lies inside +/-2 days (a formal 'no effect') for 2 of 4 sources. Where it does not, the interval is too wide to tell, or the effect is negative.
- **All three together** hold for 1 of 4 sources over all sites.
- **By plant type.** Deciduous forests: pre-solstice PhenoCam (GCC) -0.3, Tower NDVI +0.5, Satellite NDVI -2.7***, Satellite NIRv -0.6; evergreen forests: Tower NDVI -0.6, Satellite NDVI -3.2***, Satellite NIRv -1.6*. Groups are small (76 to 208 site-years), so single-group results are unstable.
- **Reading.** The data support a weak advancing effect of early-season GPP on the onset of senescence and little effect of whole-season GPP. They do not show a late-season effect of opposite sign. Whole-season GPP has little effect because the early effect is small and is diluted by the larger, neutral late-season GPP - not because two opposite effects cancel.

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

### 3.2 EOS90 with predictors that do not depend on leaf-out

The windows above start at the year's leaf-out, so an early spring lengthens the pre-solstice window. **Calendar windows** (60 days before / 45 days from the solstice) and the **GPP rate** (mean daily GPP over the same windows) remove that.

**Calendar windows**

| sites | EOS source | n (sites) | pre-solstice GPP | post-solstice GPP | whole season [95% CI] | all three hold | spring T | without T: pre / post / whole | controlled: pre / post / whole |
|---|---|---|---|---|---|---|---|---|---|
| all | PhenoCam (GCC) | 233 (30) | -0.2 | -0.5 | -0.6 [-2.4, +1.1] | no | -2.5** | -0.8 / +0.3 / -0.3 | -0.3 / -0.6 / -0.7 |
| all | Tower NDVI | 356 (41) | -1.8 | +0.3 | -1.1 [-3.6, +1.4] | no | -2.4* | -2.8** / +0.8 / -1.6 | -2.0* / +0.0 / -1.5 |
| all | Satellite NDVI | 446 (67) | -0.6 | -0.6 | -1.1 [-2.0, -0.1] | no | -1.6* | -1.2* / -0.0 / -1.0 | -0.5 / -0.8 / -1.1 |
| all | Satellite NIRv | 675 (88) | -1.5** | +0.6 | -0.8 [-1.7, +0.1] | no | -1.8** | -1.9*** / +1.1* / -0.7 | -1.6** / +0.5 / -0.8 |
| summer-green | PhenoCam (GCC) | 202 (25) | -0.5 | -0.9 | -1.1 [-2.9, +0.7] | no | -2.8** | -1.2 / +0.1 / -0.8 | -0.5 / -1.0 / -1.2 |
| summer-green | Tower NDVI | 341 (38) | -2.1* | +0.3 | -1.4 [-3.9, +1.2] | no | -2.3* | -3.1** / +0.8 / -1.9 | -2.4* / +0.0 / -1.8 |
| summer-green | Satellite NDVI | 426 (63) | -0.5 | -0.8 | -1.0 [-2.0, -0.1] | no | -1.6* | -1.1* / -0.1 / -1.0 | -0.3 / -0.9 / -1.0 |
| summer-green | Satellite NIRv | 605 (79) | -1.6** | +0.7 | -0.7 [-1.7, +0.2] | no | -1.9** | -2.1*** / +1.3* / -0.7 | -1.7** / +0.6 / -0.8 |
| dry-summer | PhenoCam (GCC) | 31 (5) | +2.2 | +1.9 | +3.5 [-1.5, +8.6] | no | +0.2 | +2.2 / +1.9 / +3.5 | - |
| dry-summer | Satellite NIRv | 70 (9) | -0.6 | -0.2 | -0.6 [-4.3, +3.2] | no | -0.6 | -0.4 / -0.1 / -0.4 | -1.0 / -0.0 / -0.8 |
| deciduous | PhenoCam (GCC) | 104 (9) | +0.4 | -1.7 | -1.1 [-3.6, +1.3] | no | +0.0 | +0.4 / -1.7 / -1.1 | +2.8 / -4.7 / -1.6 |
| deciduous | Tower NDVI | 174 (17) | -0.4 | +2.5* | +1.8 [+0.1, +3.6] | no | -1.8 | -1.2 / +3.1** / +1.6 | -0.5 / +2.1 / +1.4 |
| deciduous | Satellite NDVI | 170 (23) | -1.5 | +1.4 | -0.1 [-1.5, +1.4] | no | +0.7 | -1.1 / +1.0 / -0.0 | -0.2 / +0.6 / +0.3 |
| deciduous | Satellite NIRv | 195 (25) | -1.1 | +1.3 | +0.2 [-2.0, +2.4] | no | -1.7 | -2.1** / +2.3* / +0.2 | -1.2 / +1.2 / +0.1 |
| evergreen | Tower NDVI | 86 (12) | -3.3 | +0.3 | -2.1 [-6.1, +1.8] | no | -3.1 | -5.1* / +0.7 / -3.7 | -3.4 / +0.2 / -2.5 |
| evergreen | Satellite NDVI | 82 (15) | -0.7 | -1.3 | -1.7 [-3.5, +0.1] | no | -2.0 | -1.3 / -0.9 / -1.9 | -0.2 / -2.0 / -1.9 |
| evergreen | Satellite NIRv | 211 (27) | -1.9** | +0.0 | -1.6 [-2.9, -0.3] | no | -1.5 | -2.2*** / +0.3 / -1.5 | -1.7** / -0.2 / -1.6 |
| grass/shrub | PhenoCam (GCC) | 118 (18) | -1.8 | +0.9 | -0.6 [-3.6, +2.3] | no | -4.1** | -1.8 / +2.1 / +0.3 | -1.9 / +0.9 / -0.7 |
| grass/shrub | Tower NDVI | 96 (12) | -2.0 | -5.4 | -6.4 [-9.7, -3.1] | no | -2.3 | -2.7 / -4.9 / -6.6 | -2.6 / -5.4 / -6.9 |
| grass/shrub | Satellite NDVI | 194 (29) | -1.2 | -2.0 | -2.4 [-4.3, -0.4] | no | -2.8* | -1.6 / -1.2 / -2.1 | -1.4 / -2.2 / -2.7 |
| grass/shrub | Satellite NIRv | 269 (36) | -1.6 | +0.2 | -1.1 [-2.7, +0.5] | no | -2.3 | -1.6 / +0.6 / -0.8 | -1.8* / +0.2 / -1.2 |

**GPP rate**

| sites | EOS source | n (sites) | pre-solstice GPP | post-solstice GPP | whole season [95% CI] | all three hold | spring T | without T: pre / post / whole | controlled: pre / post / whole |
|---|---|---|---|---|---|---|---|---|---|
| all | PhenoCam (GCC) | 174 (22) | -0.5 | -0.8 | -1.5 [-3.4, +0.5] | no | -2.4** | -1.3 / +0.2 / -1.2 | -0.2 / -1.0 / -1.6 |
| all | Tower NDVI | 303 (35) | -1.5 | +0.4 | -1.4 [-3.5, +0.7] | no | -2.4* | -2.2* / +1.0 / -1.5 | -2.3* / +0.6 / -1.7 |
| all | Satellite NDVI | 362 (56) | +0.3 | -0.2 | -0.3 [-1.4, +0.8] | no | -0.8 | -0.0 / +0.2 / -0.2 | -0.7 / +0.4 / -0.5 |
| all | Satellite NIRv | 606 (81) | -1.3** | +0.6 | -0.7 [-1.5, +0.1] | no | -1.7** | -1.7*** / +1.1* / -0.7 | -1.2* / +0.5 / -0.8 |
| summer-green | PhenoCam (GCC) | 167 (20) | -0.4 | -1.1 | -1.7 [-3.8, +0.5] | no | -2.5** | -1.3 / +0.1 / -1.3 | -0.1 / -1.2 / -1.8 |
| summer-green | Tower NDVI | 300 (34) | -1.5 | +0.3 | -1.5 [-3.7, +0.7] | no | -2.4* | -2.2* / +0.9 / -1.6 | -2.3 / +0.5 / -1.9 |
| summer-green | Satellite NDVI | 351 (54) | +0.5 | -0.3 | -0.3 [-1.5, +0.9] | no | -0.7 | +0.2 / +0.0 / -0.2 | -0.5 / +0.2 / -0.6 |
| summer-green | Satellite NIRv | 564 (76) | -1.1* | +0.5 | -0.7 [-1.5, +0.2] | no | -2.0** | -1.7*** / +1.1* / -0.8 | -1.0* / +0.4 / -0.7 |
| dry-summer | Satellite NIRv | 42 (5) | -2.1 | +0.8 | -0.6 [-3.3, +2.0] | no | +0.2 | -2.2 / +0.7 / -0.5 | -1.3 / +0.3 / -0.7 |
| deciduous | PhenoCam (GCC) | 98 (8) | +2.2 | -3.6 | -1.0 [-3.1, +1.2] | no | -1.5 | +1.4 / -2.8 / -1.1 | +1.7 / -4.3 / -2.5 |
| deciduous | Tower NDVI | 170 (17) | +0.7 | +1.7 | +1.8 [+0.1, +3.4] | no | -1.9 | +0.0 / +2.3* / +1.8 | +0.5 / +1.9 / +1.2 |
| deciduous | Satellite NDVI | 170 (23) | +0.9 | +0.5 | +1.3 [-0.6, +3.1] | no | -0.2 | +0.7 / +0.7 / +1.3 | -0.7 / +1.3 / +0.7 |
| deciduous | Satellite NIRv | 196 (25) | -0.1 | +0.3 | +0.3 [-1.8, +2.4] | no | -2.5** | -1.4* / +1.7 / +0.3 | -0.1 / +0.2 / +0.1 |
| evergreen | Tower NDVI | 79 (12) | -4.5* | +0.6 | -3.3 [-7.2, +0.7] | no | -3.0** | -5.0* / +0.9 / -3.5 | -8.2* / +1.4 / -4.5 |
| evergreen | Satellite NDVI | 76 (14) | -1.7 | +0.2 | -1.0 [-3.0, +0.9] | no | -1.2 | -2.0 / +0.5 / -1.1 | -1.9 / -0.2 / -1.5 |
| evergreen | Satellite NIRv | 208 (26) | -2.2** | +0.5 | -1.3 [-2.4, -0.3] | no | -1.3 | -2.4** / +0.8 / -1.3 | -2.0** / +0.3 / -1.4 |
| grass/shrub | PhenoCam (GCC) | 65 (11) | -3.4 | +1.8* | -1.9 [-6.2, +2.4] | no | -3.2** | -3.9 / +3.2** / -1.5 | -3.5 / +1.9 / -1.9 |
| grass/shrub | Tower NDVI | 54 (6) | -1.4 | -4.3 | -5.0 [-8.4, -1.6] | no | -3.4 | -2.3 / -3.4 / -5.1 | -7.2 / -1.5 / -4.7 |
| grass/shrub | Satellite NDVI | 116 (19) | +0.2 | -1.5 | -1.3 [-3.2, +0.5] | no | -1.3 | -0.1 / -0.9 / -1.4 | -0.1 / -1.2 / -1.5 |
| grass/shrub | Satellite NIRv | 202 (30) | -1.3 | +0.2 | -1.1 [-2.3, +0.2] | no | -2.0 | -1.3 / +0.7 / -1.3 | -1.2 / +0.4 / -1.2 |

**What this shows**

- **Calendar windows:** pre-solstice PhenoCam (GCC) -0.2, Tower NDVI -1.8, Satellite NDVI -0.6, Satellite NIRv -1.5**. **Rate:** PhenoCam (GCC) -0.5, Tower NDVI -1.5, Satellite NDVI +0.3, Satellite NIRv -1.3**. Compare the sum from leaf-out: PhenoCam (GCC) -0.7, Tower NDVI -1.5, Satellite NDVI -2.7***, Satellite NIRv -1.6**.
- Over the three versions the pre-solstice slope is negative in 11 of 12 source x version cells and significant in 4 of 12. The effect is smaller when the predictor cannot carry leaf-out date, so part of the leaf-out-to-solstice result is an effect of an early spring as such.
- The source that keeps a significant effect in every version is the one with the largest sample; the others agree in sign but cannot confirm an effect of 1-2 days (see the power calculation, section 13).
- Formal verdict, all site groups: all three parts hold in 3 of 20 tests (sum), 0 of 21 tests (calendar), 0 of 20 tests (rate); H1 alone in 12, 4, 4.

### 3.3 End of senescence (EOS10), windows ending at the site mean EOS10

| sites | EOS source | n (sites) | pre-solstice GPP | post-solstice GPP | whole season [95% CI] | all three hold | spring T | without T: pre / post / whole | controlled: pre / post / whole |
|---|---|---|---|---|---|---|---|---|---|
| all | PhenoCam (GCC) | 201 (25) | +0.0 | +0.6 | +1.0 [-1.5, +3.5] | no | +0.9 | +0.2 / +0.3 / +0.8 | -0.6 / +1.0 / +0.8 |
| all | Tower NDVI | 313 (37) | -0.9 | +1.1 | +0.8 [-1.0, +2.6] | no | +0.3 | -0.8 / +1.1 / +0.8 | -2.4* / +1.7* / +0.4 |
| all | Satellite NDVI | 380 (59) | +0.9 | +0.3 | +1.3 [-0.1, +2.7] | no | -0.4 | +0.7 / +0.5 / +1.4 | -0.5 / +1.0 / +0.9 |
| all | Satellite NIRv | 624 (85) | +2.0 | +0.1 | +1.5 [+0.2, +2.9] | no | -1.9** | +1.6 / +0.6 / +1.7 | +0.3 / +1.0 / +1.3 |
| summer-green | PhenoCam (GCC) | 174 (21) | +0.0 | +0.6 | +1.0 [-1.8, +3.8] | no | +1.2 | +0.4 / +0.1 / +0.8 | -1.0 / +1.4 / +1.0 |
| summer-green | Tower NDVI | 298 (34) | -0.6 | +0.9 | +0.8 [-1.0, +2.7] | no | +0.2 | -0.5 / +0.9 / +0.9 | -2.3* / +1.5* / +0.4 |
| summer-green | Satellite NDVI | 366 (56) | +0.9 | -0.0 | +1.1 [-0.3, +2.6] | no | -0.3 | +0.8 / +0.1 / +1.1 | -0.3 / +0.7 / +0.8 |
| summer-green | Satellite NIRv | 570 (78) | +1.8 | +0.1 | +1.5 [+0.0, +2.9] | no | -1.6* | +1.4 / +0.6 / +1.6 | -0.1 / +1.3* / +1.2 |
| dry-summer | Satellite NIRv | 54 (7) | +3.3 | -0.0 | +2.2 [-1.1, +5.5] | no | -4.1 | +3.7 / +0.4 / +2.9 | +2.7 / +0.2 / +2.2 |
| deciduous | PhenoCam (GCC) | 105 (9) | -1.8* | +3.7** | +2.2 [-1.0, +5.4] | no | +1.7 | -1.0 / +2.9** / +2.2 | -3.7** / +5.9* / +2.8 |
| deciduous | Tower NDVI | 168 (17) | -1.2 | +2.0 | +1.1 [-1.3, +3.4] | no | +0.4 | -1.1 / +1.8* / +1.1 | -2.8* / +2.6* / +0.7 |
| deciduous | Satellite NDVI | 169 (23) | +0.7 | -0.3 | +0.3 [-1.2, +1.9] | no | -0.0 | +0.7 / -0.3 / +0.3 | -1.0* / +0.8 / +0.1 |
| deciduous | Satellite NIRv | 193 (25) | -0.1 | +1.2 | +0.8 [-1.3, +2.9] | no | +2.0* | +0.8 / +0.3 / +0.7 | -1.6* / +2.2* / +0.6 |
| evergreen | Tower NDVI | 79 (12) | +2.2 | -1.1 | -0.0 [-3.2, +3.2] | no | -0.5 | +2.0 / -1.1 / +0.4 | +0.6 / +0.2 / -0.7 |
| evergreen | Satellite NDVI | 75 (14) | -0.7 | +0.9 | +0.3 [-3.0, +3.7] | no | +1.7 | +0.4 / +0.2 / +0.4 | -1.9 / +2.1 / +0.1 |
| evergreen | Satellite NIRv | 206 (26) | +1.2 | +0.7 | +1.5 [-0.9, +3.8] | no | -1.8 | +0.9 / +1.0 / +1.5 | +0.1 / +1.5* / +1.4 |
| grass/shrub | PhenoCam (GCC) | 85 (13) | -0.5 | -4.8 | -2.0 [-6.0, +2.0] | no | +0.1 | -0.4 / -4.8 / -2.4 | -0.5 / -4.8 / -2.5 |
| grass/shrub | Tower NDVI | 66 (8) | -3.5 | +1.8 | +0.7 [-4.8, +6.2] | no | +0.9 | -3.1 / +1.5 / +0.8 | -4.0 / +0.1 / +0.4 |
| grass/shrub | Satellite NDVI | 136 (22) | +1.4 | +2.2* | +3.6 [+1.4, +5.7] | no | -1.3 | +1.0 / +2.7* / +3.7 | +0.6 / +2.2* / +3.1 |
| grass/shrub | Satellite NIRv | 225 (34) | +3.2 | -0.1 | +2.6 [+0.8, +4.4] | no | -4.1*** | +3.1 / +0.9 / +3.2 | +1.2 / +0.7 / +2.3 |

**EOS10, calendar windows**

| sites | EOS source | n (sites) | pre-solstice GPP | post-solstice GPP | whole season [95% CI] | all three hold | spring T | without T: pre / post / whole | controlled: pre / post / whole |
|---|---|---|---|---|---|---|---|---|---|
| all | PhenoCam (GCC) | 233 (30) | +0.4 | +0.3 | +0.5 [-1.4, +2.5] | no | +0.8 | +0.6 / -0.0 / +0.4 | +0.2 / +0.4 / +0.5 |
| all | Tower NDVI | 356 (41) | -0.9 | +1.2 | +0.3 [-1.5, +2.1] | no | +0.5 | -0.7 / +1.1 / +0.3 | -1.6* / +1.1 / -0.3 |
| all | Satellite NDVI | 446 (67) | +0.7 | +0.1 | +0.6 [-0.7, +1.9] | no | -1.0 | +0.3 / +0.4 / +0.6 | +0.3 / -0.1 / +0.1 |
| all | Satellite NIRv | 675 (88) | +1.4 | -0.4 | +0.9 [-0.2, +2.0] | no | -2.4*** | +0.9 / +0.3 / +0.9 | +0.8 / +0.0 / +0.7 |
| summer-green | PhenoCam (GCC) | 202 (25) | -0.0 | +0.4 | +0.3 [-1.9, +2.5] | no | +1.2 | +0.3 / -0.1 / +0.2 | -0.5 / +0.7 / +0.2 |
| summer-green | Tower NDVI | 341 (38) | -1.2 | +1.2 | +0.1 [-1.7, +2.0] | no | +0.6 | -0.9 / +1.1 / +0.2 | -1.9* / +1.1 / -0.6 |
| summer-green | Satellite NDVI | 426 (63) | +1.0 | -0.6 | +0.3 [-0.9, +1.6] | no | -1.1 | +0.6 / -0.2 / +0.3 | +0.6 / -0.7 / -0.1 |
| summer-green | Satellite NIRv | 605 (79) | +1.2 | -0.4 | +0.6 [-0.4, +1.7] | no | -2.2** | +0.6 / +0.2 / +0.7 | +0.4 / +0.1 / +0.4 |
| deciduous | PhenoCam (GCC) | 104 (9) | -2.7* | +3.3* | +0.7 [-1.7, +3.1] | no | +1.8 | -1.8** / +2.4* / +0.7 | -4.6* / +5.7 / +1.0 |
| deciduous | Tower NDVI | 174 (17) | -1.6 | +2.6* | +1.0 [-1.2, +3.2] | no | +0.4 | -1.4 / +2.5** / +0.9 | -2.3** / +2.4** / +0.4 |
| deciduous | Satellite NDVI | 170 (23) | +0.7 | -0.4 | +0.2 [-1.3, +1.8] | no | -0.2 | +0.5 / -0.3 / +0.2 | -0.3 / -0.3 / -0.5 |
| deciduous | Satellite NIRv | 195 (25) | -0.9 | +0.6 | -0.2 [-1.3, +0.9] | no | +1.9* | +0.2 / -0.4 / -0.2 | -1.3* / +0.9 / -0.2 |

**What this shows**

- For the END of senescence the early-season effect is absent: pre-solstice GPP gives PhenoCam (GCC) +0.0, Tower NDVI -0.9, Satellite NDVI +0.9, Satellite NIRv +2.0 days per SD, negative in 1 of 4 sources. H1 holds in 1 of 20 tests.
- Post-solstice GPP is PhenoCam (GCC) +0.6, Tower NDVI +1.1, Satellite NDVI +0.3, Satellite NIRv +0.1. A positive value here is not evidence that GPP delays senescence: the window runs into the weeks in which the canopy is senescing, and a canopy that stays green longer photosynthesises more in those weeks. When the window stops at the mean onset of senescence instead, the post-solstice slope is PhenoCam (GCC) +0.4, Tower NDVI +0.1, Satellite NDVI -0.5, Satellite NIRv -1.3.
- So whatever early GPP does to the onset of senescence, it does not carry through to its end. The end of the season is set by something else (autumn temperature and photoperiod in the literature).

![EOS10, all sites](../figure/split_gpp_cancellation/coefficients_EOS10_all.png)

*EOS10, all sites*

### 3.4 Same-year windows and the window-length effect

When the post-solstice window ends at the same year's EOS, a later EOS makes the window longer and its cumulative GPP larger by construction. The last column is the post-solstice slope produced by window length alone (the site's average GPP curve, no year-specific GPP).

| EOS source | pre | post | whole season | post from window length alone |
|---|---|---|---|---|
| PhenoCam (GCC) | -1.3 | +6.2 | +5.7 | +10.5 |
| Tower NDVI | -3.0 | +10.2 | +6.6 | +13.9 |
| Satellite NDVI | -4.0 | +5.8 | +2.6 | +9.4 |
| Satellite NIRv | -2.7 | +8.0 | +4.8 | +10.4 |

**What this shows**

- With windows ending at the same year's EOS90 the post-solstice slope is +5.8 to +10.2 days per SD, far larger than anything above. Window length alone - no year-specific GPP at all - produces +9.4 to +13.9. The apparent delaying effect of late-season GPP in this version is therefore almost entirely arithmetic: a later EOS makes the window longer, and a longer window holds more GPP.
- This is why every other analysis in this report ends its windows at the site's mean EOS. Published positive relations between growing-season or late-season productivity and senescence date that use same-year windows should be read with this in mind.

## 4. One-to-one plots: pooled vs within-site (step 40)

The same pairs drawn twice. Pooled plots mix differences between sites with year-to-year changes; the within-site plots keep only the latter and correspond to the tests of section 3. Growing-season totals and means run from the year's leaf-out to the site's mean EOS10 (not the same-year EOS10). All pairs: `figure/predictor_correlations/` and `figure/predictor_correlations_within_site/`.

**Pooled correlation with EOS90**

| predictor | GCC | NDVI | NDVI_tower | NIRv |
|---|---|---|---|---|
| gpp_sos10_to_solstice | -0.08 | 0.22 | -0.05 | -0.30 |
| gpp_solstice_to_eos90 | 0.59 | 0.63 | 0.60 | 0.32 |
| gpp_solstice_to_eos90_fixed | 0.54 | 0.59 | 0.43 | 0.19 |
| total_gpp_growing_season | 0.30 | 0.44 | 0.24 | -0.04 |

**Within-site correlation with EOS90** (year minus site mean; sites with at least 3 years)

| predictor | GCC | NDVI | NDVI_tower | NIRv |
|---|---|---|---|---|
| gpp_sos10_to_solstice | -0.10 | -0.18 | -0.11 | -0.10 |
| gpp_solstice_to_eos90 | 0.45 | 0.35 | 0.63 | 0.58 |
| gpp_solstice_to_eos90_fixed | -0.11 | 0.00 | 0.03 | 0.02 |
| total_gpp_growing_season | -0.03 | -0.03 | -0.00 | 0.06 |

**What this shows**

- **GPP from leaf-out to the solstice.** Pooled over all site-years the correlation with EOS90 is -0.30 to +0.22 and changes sign between sources. Within sites it is -0.18 to -0.10, negative in 4 of 4 sources. The pooled plot is dominated by differences between sites (productive sites differ from unproductive ones in many ways); only the within-site value speaks to the hypothesis.
- **GPP from the solstice to the same-year EOS90:** within sites +0.35 to +0.63, a strong positive relation in every source. **To the site's mean EOS90:** -0.11 to +0.03. The relation vanishes when the window length is fixed. This pair of plots is the clearest picture of the window-length artefact.
- **Between sites the post-solstice relation is real:** pooled, with the fixed window, r = +0.19 to +0.59. Sites with a productive late summer senesce later. That is a statement about sites, not about what a productive late summer does at a given site.
- **GPP over the whole growing season:** within sites -0.03 to +0.06, i.e. no relation with the onset of senescence in any source. This is the Lu-type result.
- **Mean growing-season temperature** (within sites -0.17 to -0.09) and **leaf-out date** (-0.04 to +0.07) are also weakly related to EOS90. With the window ending at the same year's EOS the temperature correlation was about -0.5; that was the same artefact.

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
| CUE_annual | CUE_pre | 1968 | 0.54 | 0.65 |
| CUE_annual | CUE_post | 2001 | 0.57 | 0.71 |
| CUE_pre | CUE_post | 1961 | 0.40 | 0.62 |
| NPP_pre | NPPd_pre | 1968 | 0.79 | 0.94 |
| NPP_post | NPPd_post | 2001 | 0.78 | 0.93 |

![Seasonal course of CUE by plant type (step 32)](../figure/cue_seasonal/seasonal_curve.png)

*Seasonal course of CUE by plant type (step 32)*

![Annual CUE against pre- and post-solstice CUE, within sites (step 32)](../figure/cue_seasonal/annual_vs_seasonal.png)

*Annual CUE against pre- and post-solstice CUE, within sites (step 32)*

**What this shows**

- Annual CUE averages 0.52: about 52% of GPP ends up as NPP and 48% is respired by the plants. That is in the range reported for forests (0.4-0.6). The spread between site-years is large (SD 0.13), and part of it is estimation noise: the method infers CUE from how respiration follows GPP from day to day.
- The seasonal estimate gives 0.62 before and 0.59 after the solstice: a slightly larger share of the fixed carbon is retained early in the season.
- How much the time scales share (within-site r): CUE_annual vs CUE_pre +0.54; CUE_annual vs CUE_post +0.57; CUE_pre vs CUE_post +0.40; NPP_pre vs NPPd_pre +0.79; NPP_post vs NPPd_post +0.78. Where the correlation is low, the annual and the seasonal CUE carry different information; where NPP series correlate above 0.9, NPP is GPP in other units and cannot tell a sink effect from a source effect.
- Consequence for everything below: NPP and the respiration terms derived from CUE are GPP multiplied by a factor that varies little, so 'sink' variables built this way mostly repeat the GPP result.

## 6. Which carbon window relates to which EOS? (step 33)

Single-predictor within-site models, windows ending at the site mean EOS (fixed anchors). Table: `data/eos_window_scan.csv`.

**EOS10, window leaf-out to solstice** (days per +1 SD; * p<0.05, ** p<0.01, *** p<0.001, uncorrected)

| predictor | PhenoCam (GCC) | Tower NDVI | Satellite NDVI | Satellite NIRv |
|---|---|---|---|---|
| GPP cum | +0.2 | -0.6 | +1.0 | +1.6* |
| NPP cum | -0.4 | -0.4 | +0.1 | +0.9 |
| NPPd cum | -2.2 | +0.0 | +0.3 | +0.8 |
| GPP mean | -0.2 | -2.9** | -1.1* | +0.2 |
| CUEd mean | -2.3 | +0.2 | -0.6 | +0.1 |

**EOS10, window solstice to mean EOS10** (days per +1 SD; * p<0.05, ** p<0.01, *** p<0.001, uncorrected)

| predictor | PhenoCam (GCC) | Tower NDVI | Satellite NDVI | Satellite NIRv |
|---|---|---|---|---|
| GPP cum | +0.8 | +1.4 | +1.3 | +1.3* |
| NPP cum | +0.7 | +0.5 | +0.1 | +0.5 |
| NPPd cum | +0.2 | -0.3 | +0.2 | +0.2 |
| GPP mean | +0.8 | +1.5 | +1.9* | +1.7** |
| CUEd mean | -0.2 | -1.1 | -0.6 | -0.5 |

**EOS50, window leaf-out to solstice** (days per +1 SD; * p<0.05, ** p<0.01, *** p<0.001, uncorrected)

| predictor | PhenoCam (GCC) | Tower NDVI | Satellite NDVI | Satellite NIRv |
|---|---|---|---|---|
| GPP cum | -0.6 | -1.4 | -0.5 | -0.1 |
| NPP cum | -0.9 | -1.1 | -1.0* | -0.4 |
| NPPd cum | -1.9 | -0.6 | -0.6* | -0.3 |
| GPP mean | -0.5 | -2.2* | -0.4 | -0.3 |
| CUEd mean | -1.0 | +0.4 | -0.6 | -0.0 |

**EOS50, window solstice to mean EOS10** (days per +1 SD; * p<0.05, ** p<0.01, *** p<0.001, uncorrected)

| predictor | PhenoCam (GCC) | Tower NDVI | Satellite NDVI | Satellite NIRv |
|---|---|---|---|---|
| GPP cum | +0.7 | +0.3 | +0.8 | +1.3** |
| NPP cum | +0.7 | +0.1 | +0.1 | +0.8* |
| NPPd cum | +0.5 | -0.2 | +0.1 | +1.0* |
| GPP mean | +0.8 | +0.5 | +1.2* | +1.5** |
| CUEd mean | +0.4 | -0.4 | -0.3 | +0.3 |

**Within-site model vs the random-intercept mixed model used before** (1512 cells): median absolute effect 0.77 vs 1.84 days per SD; p < 0.05 in 418 vs 499 cells; same sign in 87% of cells.

**What this shows**

- **Leaf-out to solstice, cumulative GPP.** Effect on EOS90: PhenoCam (GCC) -1.2, Tower NDVI -1.9, Satellite NDVI -1.9***, Satellite NIRv -1.3**; on EOS10: PhenoCam (GCC) +0.2, Tower NDVI -0.6, Satellite NDVI +1.0, Satellite NIRv +1.6*. The early-season effect is on the onset of senescence, not on its end - the same pattern as in section 3, here without covariates.
- **Solstice to mean EOS10, cumulative GPP on EOS10:** PhenoCam (GCC) +0.8, Tower NDVI +1.4, Satellite NDVI +1.3, Satellite NIRv +1.3*. These positive values come from a window that overlaps senescence (late senescence causes late-season GPP), so they are not evidence of a delaying effect.
- **Whole season (leaf-out to EOS10) on EOS10:** with the same-year window PhenoCam (GCC) +3.6, Tower NDVI +3.0**, Satellite NDVI +2.3*, Satellite NIRv +3.0***; with the fixed window PhenoCam (GCC) +0.8, Tower NDVI +0.8, Satellite NDVI +1.4, Satellite NIRv +1.7*. The difference between the two lines is the window-length effect.
- **NPP in place of GPP** (leaf-out to solstice, EOS90): -2.0 to -0.9 days per SD over sources and CUE versions - the same as GPP, as expected from section 5.
- **Mixed model against within-site model:** see the line above. The mixed model, used in the first versions of this project, gives effects about twice as large because a random site intercept does not fully remove differences between sites. All numbers in this report are within-site.

![PhenoCam (GCC): GPP windows against EOS, fixed anchors (step 33)](../figure/eos_window_scan/fixed_GCC_GPP.png)

*PhenoCam (GCC): GPP windows against EOS, fixed anchors (step 33)*

![PhenoCam (GCC): NPP (daily CUE) windows against EOS, fixed anchors (step 33)](../figure/eos_window_scan/fixed_GCC_NPPd.png)

*PhenoCam (GCC): NPP (daily CUE) windows against EOS, fixed anchors (step 33)*

![Satellite NDVI: GPP windows against EOS, fixed anchors (step 33)](../figure/eos_window_scan/fixed_NDVI_GPP.png)

*Satellite NDVI: GPP windows against EOS, fixed anchors (step 33)*

![Satellite NDVI: NPP (daily CUE) windows against EOS, fixed anchors (step 33)](../figure/eos_window_scan/fixed_NDVI_NPPd.png)

*Satellite NDVI: NPP (daily CUE) windows against EOS, fixed anchors (step 33)*

## 7. When around the solstice is the relation strongest? (step 34)

Mean flux in 15- and 30-day windows starting 120 days before to 75 days after the solstice, within-site models. `at edge` = the minimum is on the first or last window of the scan. Table: `data/eos_solstice_sliding_scan.csv`; the most negative window per case:

**EOS90, 30-day windows**

| carbon | EOS source | window start (days from solstice) | beta_days_per_sd | p_value | n_obs | at edge |
|---|---|---|---|---|---|---|
| CUEd | PhenoCam (GCC) | -30 | -0.2 | 0.877 | 232 | False |
| CUEd | Satellite NDVI | -10 | -0.8 | 0.092 | 446 | False |
| CUEd | Tower NDVI | 55 | -0.7 | 0.579 | 362 | False |
| CUEd | Satellite NIRv | -5 | -0.5 | 0.262 | 678 | False |
| GPP | PhenoCam (GCC) | -30 | -0.9 | 0.397 | 234 | False |
| GPP | Satellite NDVI | -55 | -1.5 | 0.006 | 455 | False |
| GPP | Tower NDVI | -60 | -2.6 | 0.029 | 356 | False |
| GPP | Satellite NIRv | -50 | -1.7 | 0.000 | 682 | False |
| NPP | PhenoCam (GCC) | -105 | -0.9 | 0.452 | 233 | False |
| NPP | Satellite NDVI | -75 | -1.6 | 0.002 | 454 | False |
| NPP | Tower NDVI | -85 | -2.4 | 0.037 | 352 | False |
| NPP | Satellite NIRv | -55 | -1.4 | 0.002 | 682 | False |
| NPPd | PhenoCam (GCC) | -95 | -1.1 | 0.364 | 232 | False |
| NPPd | Satellite NDVI | -55 | -1.4 | 0.015 | 450 | False |
| NPPd | Tower NDVI | -85 | -2.2 | 0.011 | 352 | False |
| NPPd | Satellite NIRv | -55 | -1.2 | 0.012 | 680 | False |

**EOS10, 30-day windows**

| carbon | EOS source | window start (days from solstice) | beta_days_per_sd | p_value | n_obs | at edge |
|---|---|---|---|---|---|---|
| CUEd | PhenoCam (GCC) | -20 | -3.1 | 0.152 | 232 | False |
| CUEd | Satellite NDVI | -80 | -1.3 | 0.036 | 448 | False |
| CUEd | Tower NDVI | 5 | -1.7 | 0.042 | 362 | False |
| CUEd | Satellite NIRv | 45 | -0.7 | 0.128 | 681 | False |
| GPP | PhenoCam (GCC) | 20 | -0.9 | 0.537 | 234 | False |
| GPP | Satellite NDVI | -120 | -0.7 | 0.264 | 447 | True |
| GPP | Tower NDVI | -120 | -2.1 | 0.102 | 337 | True |
| GPP | Satellite NIRv | 55 | +0.1 | 0.902 | 684 | False |
| NPP | PhenoCam (GCC) | -90 | -1.4 | 0.419 | 233 | False |
| NPP | Satellite NDVI | -120 | -1.3 | 0.025 | 446 | True |
| NPP | Tower NDVI | -120 | -2.2 | 0.058 | 337 | True |
| NPP | Satellite NIRv | 20 | -0.4 | 0.444 | 685 | False |
| NPPd | PhenoCam (GCC) | -25 | -2.6 | 0.149 | 232 | False |
| NPPd | Satellite NDVI | -120 | -0.9 | 0.186 | 442 | True |
| NPPd | Tower NDVI | -120 | -1.8 | 0.100 | 337 | True |
| NPPd | Satellite NIRv | 45 | -0.6 | 0.250 | 681 | False |

**What this shows**

- **EOS90.** The 30-day window in which GPP has its most negative relation with the onset of senescence starts 30 days before the solstice for PhenoCam (GCC) (-0.9), 60 days before the solstice for Tower NDVI (-2.6*), 55 days before the solstice for Satellite NDVI (-1.5**), 50 days before the solstice for Satellite NIRv (-1.7***). It is significant in 3 of 4 sources.
- Where it is significant the window starts 50 to 60 days before the solstice, i.e. in late April to early May and running into late May or early June. The sensitive period is spring, the weeks after leaf-out - earlier than the weeks around the solstice, where Zohner et al. (2023) place the switch.
- The scan covers 120 days before to 75 days after the solstice. A minimum on the edge of that range would mean the true optimum lies outside it; this is the case for 0 of 4 sources.
- **EOS10.** The most negative window gives -2.1 to +0.1 days per SD and is significant in 0 of 4 sources: no period of the season has a GPP that predicts the END of senescence.
- A caution on reading the minima: taking the most negative of 40 windows overstates the effect (the p-values are not corrected for the search). The location of the minimum is more informative than its size.

![GPP: effect on EOS of 30-day windows by start date relative to the solstice (step 34)](../figure/eos_solstice_sliding_scan/GPP_L30.png)

*GPP: effect on EOS of 30-day windows by start date relative to the solstice (step 34)*

![NPPd: effect on EOS of 30-day windows by start date relative to the solstice (step 34)](../figure/eos_solstice_sliding_scan/NPPd_L30.png)

*NPPd: effect on EOS of 30-day windows by start date relative to the solstice (step 34)*

![CUEd: effect on EOS of 30-day windows by start date relative to the solstice (step 34)](../figure/eos_solstice_sliding_scan/CUEd_L30.png)

*CUEd: effect on EOS of 30-day windows by start date relative to the solstice (step 34)*

## 8. Does carbon explain EOS beyond climate? (step 35)

Nested within-site models on identical rows. M0 climate only; + leaf-out date (SOS); + source (GPP); + sink (NPP). Tables: `data/eos_env_vs_carbon_comparison.csv`, `..._cv.csv`, `..._lrt.csv`.

**EOS10**

| EOS source | model | n_obs | aic | R2 (within sites) | gain over climate | R2 leave-one-site-out |
|---|---|---|---|---|---|---|
| PhenoCam (GCC) | M0_env | 199 | 1651 | 0.02 | 0.00 | -0.11 |
| PhenoCam (GCC) | M1_env+SOS | 199 | 1652 | 0.03 | 0.01 | -0.15 |
| PhenoCam (GCC) | M2_env+source | 199 | 1655 | 0.02 | 0.00 | -0.18 |
| PhenoCam (GCC) | M3_env+sink | 199 | 1655 | 0.02 | 0.00 | -0.20 |
| PhenoCam (GCC) | M4_env+source+sink | 199 | 1656 | 0.04 | 0.02 | -0.26 |
| PhenoCam (GCC) | M3d_env+sink_daily | 199 | 1652 | 0.04 | 0.02 | -0.14 |
| PhenoCam (GCC) | M4d_env+source+sink_daily | 199 | 1651 | 0.06 | 0.04 | -0.20 |
| PhenoCam (GCC) | M5_env+SOS+source+sink | 199 | 1657 | 0.05 | 0.02 | -0.28 |
| Tower NDVI | M5_env+SOS+source+sink | 313 | 2676 | 0.08 | 0.02 | -0.02 |
| Tower NDVI | M4d_env+source+sink_daily | 313 | 2677 | 0.07 | 0.01 | -0.03 |
| Tower NDVI | M3d_env+sink_daily | 313 | 2678 | 0.05 | 0.00 | -0.03 |
| Tower NDVI | M4_env+source+sink | 313 | 2679 | 0.06 | 0.01 | -0.03 |
| Tower NDVI | M3_env+sink | 313 | 2676 | 0.06 | 0.01 | -0.02 |
| Tower NDVI | M2_env+source | 313 | 2676 | 0.06 | 0.01 | -0.03 |
| Tower NDVI | M1_env+SOS | 313 | 2674 | 0.06 | 0.01 | -0.02 |
| Tower NDVI | M0_env | 313 | 2674 | 0.05 | 0.00 | -0.02 |
| Satellite NDVI | M5_env+SOS+source+sink | 375 | 2839 | 0.08 | 0.02 | -0.01 |
| Satellite NDVI | M4d_env+source+sink_daily | 375 | 2840 | 0.08 | 0.01 | -0.02 |
| Satellite NDVI | M3d_env+sink_daily | 375 | 2839 | 0.07 | 0.00 | 0.00 |
| Satellite NDVI | M4_env+source+sink | 375 | 2841 | 0.07 | 0.01 | -0.03 |
| Satellite NDVI | M3_env+sink | 375 | 2840 | 0.07 | 0.00 | -0.00 |
| Satellite NDVI | M2_env+source | 375 | 2839 | 0.07 | 0.00 | -0.00 |
| Satellite NDVI | M1_env+SOS | 375 | 2834 | 0.08 | 0.01 | 0.02 |
| Satellite NDVI | M0_env | 375 | 2836 | 0.07 | 0.00 | 0.02 |
| Satellite NIRv | M0_env | 621 | 4820 | 0.11 | 0.00 | 0.09 |
| Satellite NIRv | M1_env+SOS | 621 | 4820 | 0.11 | 0.00 | 0.09 |
| Satellite NIRv | M2_env+source | 621 | 4816 | 0.12 | 0.01 | 0.09 |
| Satellite NIRv | M3_env+sink | 621 | 4823 | 0.11 | 0.00 | 0.08 |
| Satellite NIRv | M4_env+source+sink | 621 | 4819 | 0.12 | 0.01 | 0.08 |
| Satellite NIRv | M3d_env+sink_daily | 621 | 4823 | 0.11 | 0.00 | 0.08 |
| Satellite NIRv | M4d_env+source+sink_daily | 621 | 4818 | 0.12 | 0.01 | 0.09 |
| Satellite NIRv | M5_env+SOS+source+sink | 621 | 4819 | 0.13 | 0.02 | 0.08 |

**EOS50**

| EOS source | model | n_obs | aic | R2 (within sites) | gain over climate | R2 leave-one-site-out |
|---|---|---|---|---|---|---|
| PhenoCam (GCC) | M0_env | 199 | 1529 | 0.02 | 0.00 | -0.16 |
| PhenoCam (GCC) | M1_env+SOS | 199 | 1530 | 0.02 | 0.00 | -0.19 |
| PhenoCam (GCC) | M2_env+source | 199 | 1532 | 0.02 | 0.00 | -0.18 |
| PhenoCam (GCC) | M3_env+sink | 199 | 1531 | 0.03 | 0.01 | -0.18 |
| PhenoCam (GCC) | M4_env+source+sink | 199 | 1532 | 0.04 | 0.02 | -0.24 |
| PhenoCam (GCC) | M3d_env+sink_daily | 199 | 1528 | 0.04 | 0.03 | -0.16 |
| PhenoCam (GCC) | M4d_env+source+sink_daily | 199 | 1530 | 0.05 | 0.03 | -0.21 |
| PhenoCam (GCC) | M5_env+SOS+source+sink | 199 | 1531 | 0.05 | 0.04 | -0.25 |
| Tower NDVI | M5_env+SOS+source+sink | 313 | 2467 | 0.08 | 0.02 | -0.03 |
| Tower NDVI | M4d_env+source+sink_daily | 313 | 2465 | 0.08 | 0.02 | -0.03 |
| Tower NDVI | M3d_env+sink_daily | 313 | 2467 | 0.07 | 0.00 | -0.03 |
| Tower NDVI | M4_env+source+sink | 313 | 2466 | 0.08 | 0.02 | -0.02 |
| Tower NDVI | M3_env+sink | 313 | 2464 | 0.07 | 0.01 | -0.02 |
| Tower NDVI | M2_env+source | 313 | 2463 | 0.08 | 0.02 | -0.02 |
| Tower NDVI | M1_env+SOS | 313 | 2466 | 0.06 | 0.00 | -0.02 |
| Tower NDVI | M0_env | 313 | 2464 | 0.06 | 0.00 | -0.01 |
| Satellite NDVI | M5_env+SOS+source+sink | 375 | 2498 | 0.07 | 0.04 | -0.03 |
| Satellite NDVI | M4d_env+source+sink_daily | 375 | 2497 | 0.06 | 0.03 | -0.02 |
| Satellite NDVI | M3d_env+sink_daily | 375 | 2503 | 0.04 | 0.01 | -0.02 |
| Satellite NDVI | M4_env+source+sink | 375 | 2496 | 0.07 | 0.04 | -0.02 |
| Satellite NDVI | M3_env+sink | 375 | 2495 | 0.06 | 0.03 | -0.00 |
| Satellite NDVI | M2_env+source | 375 | 2497 | 0.05 | 0.02 | -0.01 |
| Satellite NDVI | M1_env+SOS | 375 | 2504 | 0.03 | 0.00 | -0.03 |
| Satellite NDVI | M0_env | 375 | 2503 | 0.03 | 0.00 | -0.02 |
| Satellite NIRv | M0_env | 621 | 4461 | 0.08 | 0.00 | 0.05 |
| Satellite NIRv | M1_env+SOS | 621 | 4462 | 0.08 | 0.00 | 0.04 |
| Satellite NIRv | M2_env+source | 621 | 4447 | 0.10 | 0.03 | 0.07 |
| Satellite NIRv | M3_env+sink | 621 | 4449 | 0.10 | 0.02 | 0.07 |
| Satellite NIRv | M4_env+source+sink | 621 | 4449 | 0.11 | 0.03 | 0.06 |
| Satellite NIRv | M3d_env+sink_daily | 621 | 4454 | 0.09 | 0.02 | 0.06 |
| Satellite NIRv | M4d_env+source+sink_daily | 621 | 4449 | 0.11 | 0.03 | 0.06 |
| Satellite NIRv | M5_env+SOS+source+sink | 621 | 4444 | 0.12 | 0.04 | 0.07 |

**What this shows**

- **Climate alone** (temperature, radiation and precipitation before and after the solstice) explains 2 to 11% of the year-to-year variance of EOS10 within sites. Most of the variation in senescence date is left unexplained by seasonal climate means.
- **Adding leaf-out date, GPP and NPP** raises this to 5 to 13%, a gain of 2 to 2 percentage points. Carbon uptake adds very little to climate.
- **Prediction for a site the model has not seen** (leave-one-site-out R2): -0.11 to 0.09 for climate alone, -0.28 to 0.08 for the full model. Values near zero or below mean the relations do not transfer between sites: a model fitted on some sites does not predict the early and late years of another. Adding predictors makes this worse, not better (overfitting).
- In the full model, leaf-out date has the most consistent coefficient: PhenoCam (GCC) -1.9, Tower NDVI -3.2*, Satellite NDVI -1.5, Satellite NIRv -1.2 days per SD (a later leaf-out goes with an earlier end of season). No GPP or NPP term is consistent in sign across the four sources.

## 9. Pathways and rate vs cumulative uptake (steps 36, 37)

![PhenoCam (GCC): climate -> GPP -> EOS10 path coefficients (step 36)](../figure/eos_path_analysis/GCC_EOS10_GPP.png)

*PhenoCam (GCC): climate -> GPP -> EOS10 path coefficients (step 36)*

![Satellite NDVI: climate -> GPP -> EOS10 path coefficients (step 36)](../figure/eos_path_analysis/NDVI_EOS10_GPP.png)

*Satellite NDVI: climate -> GPP -> EOS10 path coefficients (step 36)*

**GPP rate (mean) vs cumulative GPP, leaf-out to solstice, all sites** (`data/eos_rate_vs_cumulative.csv`)

| EOS source | target | n_obs | beta_rate | p_rate | beta_cum | p_cum |
|---|---|---|---|---|---|---|
| PhenoCam (GCC) | EOS10 | 205 | -0.2 | 0.744 | +0.2 | 0.869 |
| Tower NDVI | EOS10 | 318 | -2.9 | 0.004 | -0.6 | 0.592 |
| Satellite NDVI | EOS10 | 395 | -1.1 | 0.024 | +1.0 | 0.163 |
| Satellite NIRv | EOS10 | 650 | +0.2 | 0.614 | +1.6 | 0.025 |

**What this shows**

- A cumulative sum is a rate times a duration. Between leaf-out and the solstice the duration is set by the leaf-out date, so 'more cumulative GPP' can mean 'photosynthesised faster' or 'started earlier'.
- **Rate** (mean daily GPP) on EOS10: PhenoCam (GCC) -0.2, Tower NDVI -2.9**, Satellite NDVI -1.1*, Satellite NIRv +0.2. **Cumulative**: PhenoCam (GCC) +0.2, Tower NDVI -0.6, Satellite NDVI +1.0, Satellite NIRv +1.6*. The rate is negative in 3 of 4 sources, the sum in 1 of 4.
- With rate and window length in the same model, the length term is PhenoCam (GCC) +0.7, Tower NDVI +1.7, Satellite NDVI +1.8**, Satellite NIRv +2.7*** days per SD: a longer pre-solstice window, i.e. an earlier leaf-out, goes with a LATER end of season. That positive duration effect is why the cumulative sum shows no negative effect on EOS10 while the rate does in some sources.
- The path figures decompose the same thing: climate acts on EOS10 mostly directly, and the part that passes through GPP is small.

## 10. Results by plant type, with multiple-testing correction (step 38)

Tables: `data/eos_results_by_leaf_habit.csv`, `data/eos_results_consistency.csv`.

- 780 tests; 90 (11.5%) with uncorrected p < 0.05 (about 5% expected by chance); 1 significant after FDR correction.

**Effects with the same sign in every EOS source, both ground sources agreeing, at least one significant after FDR**

_No data._

![EOS10: every predictor by plant type and EOS source (step 38)](../figure/eos_results_by_leaf_habit/EOS10_raw.png)

*EOS10: every predictor by plant type and EOS source (step 38)*

![EOS50, air temperature controlled (step 38)](../figure/eos_results_by_leaf_habit/EOS50_adj_T.png)

*EOS50, air temperature controlled (step 38)*

**What this shows**

- 780 combinations of plant type, EOS source, target (EOS10, EOS50), window and predictor were tested. 11.5% have p < 0.05 before correction, against 5% expected if nothing were going on - so there is some signal, but most single 'significant' cells are chance.
- After correcting for the number of tests, 1 remain: grass/shrub, Satellite NDVI, EOS10, post GPP rate (mean) (+5.1 days per SD).
- 71 of 208 predictor x group combinations have the same sign in every EOS source, and 0 meet the full robustness rule (same sign everywhere, both ground sources agreeing, one significant after correction).
- For the middle and the end of senescence, then, no carbon variable - GPP, NPP, CUE, as a sum or as a rate, before or after the solstice - has an effect that is consistent across plant types and EOS sources. The one consistent carbon effect in this project is the one on the ONSET of senescence (section 3).

## 11. Timing of anomalies and drought years (step 39)

Tables: `data/anomaly_timing_effects.csv`, `data/drought_years.csv`, `data/drought_eos_contrast.csv`.

- Drought years (May-September water balance at least 1 SD below the site mean): 327 of 2172 site-years (15%).

**Mean EOS shift in drought years against the other years**

| EOS source | target | sites | drought / other site-years | EOS shift in drought years (days) | p_value |
|---|---|---|---|---|---|
| PhenoCam (GCC) | EOS10 | all | 34 / 222 | -1.9 | 0.421 |
| PhenoCam (GCC) | EOS10 | deciduous | 18 / 98 | -3.4 | 0.125 |
| PhenoCam (GCC) | EOS50 | all | 34 / 222 | -1.9 | 0.100 |
| PhenoCam (GCC) | EOS50 | deciduous | 18 / 98 | -1.3 | 0.294 |
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

- Anomaly windows: 3964 tests, 11.7% with uncorrected p < 0.05; after FDR: all years 12, drought years 9, non-drought years 14.

**What this shows**

- A drought year is one whose May-September water balance (precipitation minus potential evaporation) is at least one SD below the site's mean: 15% of site-years.
- In drought years senescence is earlier in 7 of 8 source x target combinations (shift -4.8 to +0.4 days), significantly so in 3: Satellite NDVI EOS10 -4.8 days; Satellite NDVI EOS50 -4.0 days; Satellite NIRv EOS50 -2.3 days.
- The effect is clearest in the satellite indices, which have the most site-years; the ground sources agree in sign. Drought advances senescence by a few days - a larger effect than that of GPP.
- Timing of anomalies: of 3964 window x variable x sign tests, 11.7% have uncorrected p < 0.05 and 35 survive correction. No 30-day window of temperature, radiation, VPD, precipitation, water balance or GPP shifts EOS10 or EOS50 consistently across sources. The heat maps show scattered cells, not a band at a particular time of year.

![Satellite NDVI: effect on EOS10 of positive (+) and negative (-) anomalies by time of year (step 39)](../figure/anomaly_timing/NDVI_EOS10.png)

*Satellite NDVI: effect on EOS10 of positive (+) and negative (-) anomalies by time of year (step 39)*

![PhenoCam (GCC): effect on EOS10 of positive (+) and negative (-) anomalies by time of year (step 39)](../figure/anomaly_timing/GCC_EOS10.png)

*PhenoCam (GCC): effect on EOS10 of positive (+) and negative (-) anomalies by time of year (step 39)*

## 12. PhenoCam vs satellite EOS on equal terms (step 42)

The same within-site model on each source's own site-years, on the PhenoCam sites only, and on the shared site-years; predictors are calendar windows, identical for every source. Tables: `data/phenocam_vs_satellite_models.csv`, `..._difference.csv`, `..._levels.csv`.

**EOS90: days per +1 within-site SD of GPP [95% CI]** (two-sided stars)

| pair | sample | EOS source | n (sites) | pre | post | total |
|---|---|---|---|---|---|---|
| GCC vs NDVI | own | PhenoCam (GCC) | 233 (30) | -0.8 [-3.1, +1.5] | +0.3 [-2.0, +2.6] | -0.3 [-2.1, +1.5] |
| GCC vs NDVI | own | Satellite NDVI | 446 (67) | -1.2 [-2.4, +0.0] | -0.0 [-1.3, +1.2] | -1.0* [-2.0, -0.1] |
| GCC vs NDVI | shared | PhenoCam (GCC) | 99 (16) | -1.4 [-3.9, +1.1] | +0.4 [-2.2, +2.9] | -0.8 [-2.1, +0.6] |
| GCC vs NDVI | shared | Satellite NDVI | 99 (16) | -0.9 [-4.1, +2.4] | -1.6 [-4.0, +0.8] | -2.0 [-4.4, +0.5] |
| GCC vs NDVI | sites | Satellite NDVI | 183 (27) | -0.1 [-2.5, +2.2] | -1.1 [-3.1, +0.9] | -1.0 [-2.8, +0.8] |
| GCC vs NIRv | own | Satellite NIRv | 675 (88) | -1.9*** [-2.9, -0.9] | +1.1 [-0.0, +2.2] | -0.7 [-1.6, +0.2] |
| GCC vs NIRv | shared | PhenoCam (GCC) | 137 (23) | -0.8 [-2.9, +1.2] | +0.7 [-1.3, +2.8] | -0.1 [-1.5, +1.4] |
| GCC vs NIRv | shared | Satellite NIRv | 137 (23) | -0.9 [-4.0, +2.2] | +1.1 [-1.1, +3.3] | +0.2 [-1.5, +2.0] |
| GCC vs NIRv | sites | Satellite NIRv | 266 (34) | -0.9 [-3.2, +1.4] | +0.7 [-0.8, +2.3] | -0.1 [-1.8, +1.6] |

**EOS10: days per +1 within-site SD of GPP [95% CI]** (two-sided stars)

| pair | sample | EOS source | n (sites) | pre | post | total |
|---|---|---|---|---|---|---|
| GCC vs NDVI | own | PhenoCam (GCC) | 233 (30) | +0.6 [-1.5, +2.8] | -0.0 [-2.9, +2.9] | +0.4 [-1.6, +2.5] |
| GCC vs NDVI | own | Satellite NDVI | 446 (67) | +0.3 [-1.1, +1.7] | +0.4 [-1.1, +1.9] | +0.6 [-0.7, +1.9] |
| GCC vs NDVI | shared | PhenoCam (GCC) | 99 (16) | +1.1 [-1.4, +3.6] | -0.3 [-2.2, +1.7] | +0.6 [-1.7, +2.8] |
| GCC vs NDVI | shared | Satellite NDVI | 99 (16) | -0.2 [-1.9, +1.5] | +0.8 [-1.1, +2.8] | +0.5 [-1.4, +2.5] |
| GCC vs NDVI | sites | Satellite NDVI | 183 (27) | -0.6 [-2.7, +1.6] | +1.9* [+0.0, +3.8] | +1.2 [-0.7, +3.1] |
| GCC vs NIRv | own | Satellite NIRv | 675 (88) | +0.9 [-0.3, +2.0] | +0.3 [-0.9, +1.4] | +0.9 [-0.2, +2.1] |
| GCC vs NIRv | shared | PhenoCam (GCC) | 137 (23) | +0.6 [-1.5, +2.6] | +0.1 [-1.7, +2.0] | +0.6 [-1.5, +2.6] |
| GCC vs NIRv | shared | Satellite NIRv | 137 (23) | +0.4 [-1.1, +1.9] | -0.1 [-1.2, +1.1] | +0.3 [-0.9, +1.5] |
| GCC vs NIRv | sites | Satellite NIRv | 266 (34) | +1.5 [-0.6, +3.6] | +1.0 [-0.6, +2.6] | +2.0* [+0.4, +3.7] |

**Does GPP shift the two EOS measures differently?** Satellite EOS minus PhenoCam EOS, regressed on GPP (shared site-years). A slope different from zero would mean the sources respond differently.

| difference | level | n (sites) | mean difference (days) | pre | post |
|---|---|---|---|---|---|
| NDVI - GCC | EOS90 | 99 (16) | +16.7 | +0.5 [-2.8, +3.9] | -1.9 [-5.5, +1.6] |
| NDVI - GCC | EOS50 | 99 (16) | +15.9 | +0.0 [-1.2, +1.2] | -0.6 [-2.0, +0.9] |
| NDVI - GCC | EOS10 | 99 (16) | +14.3 | -1.3 [-4.9, +2.4] | +1.1 [-1.0, +3.1] |
| NIRv - GCC | EOS90 | 137 (23) | -8.1 | -0.1 [-2.9, +2.7] | +0.4 [-2.6, +3.4] |
| NIRv - GCC | EOS50 | 137 (23) | -0.0 | -0.0 [-2.1, +2.1] | -0.0 [-1.8, +1.7] |
| NIRv - GCC | EOS10 | 137 (23) | +8.9 | -0.1 [-2.2, +2.0] | -0.2 [-2.3, +1.9] |

**What this shows**

- PhenoCam and the satellite cover different sites and years, so different results could come from the sample or from what each instrument sees. Fitting both on the SAME site-years separates the two.
- **Shared site-years, PhenoCam vs NDVI** (99 site-years, 16 sites): pre-solstice GPP PhenoCam (GCC) -1.4, Satellite NDVI -0.9; whole season PhenoCam (GCC) -0.8, Satellite NDVI -2.0.
- **Shared site-years, PhenoCam vs NIRv** (137 site-years, 23 sites): pre-solstice GPP PhenoCam (GCC) -0.8, Satellite NIRv -0.9; whole season PhenoCam (GCC) -0.1, Satellite NIRv +0.2.
- On the shared samples nothing is significant (0 of 4 pre-solstice slopes): with about 100 site-years an effect of 1-2 days cannot be detected. The comparison can show that the two sources do not contradict each other; it cannot confirm the effect.
- **Difference model.** If GPP moved the satellite EOS90 and the PhenoCam EOS90 differently, the difference between them would depend on GPP. The slopes are -0.1 to +0.5 (pre) and -1.9 to +0.4 (post) days per SD, significant in 0 of 2 cases: the two sources respond to GPP in the same way.
- **Systematic offsets** (satellite minus PhenoCam, mean days): NDVI - PhenoCam EOS90 +17; NDVI - PhenoCam EOS50 +16; NDVI - PhenoCam EOS10 +14; NIRv - PhenoCam EOS90 -8; NIRv - PhenoCam EOS50 -0; NIRv - PhenoCam EOS10 +9. The sources date the same stage at different calendar days, but that offset is constant and drops out of within-site analyses.
- **Do they measure the same event?** Within-site correlation of the same stage in both sources: Satellite NDVI EOS90 r = 0.45; Satellite NDVI EOS50 r = 0.62; Satellite NDVI EOS10 r = 0.32; Satellite NIRv EOS90 r = 0.24; Satellite NIRv EOS50 r = 0.36; Satellite NIRv EOS10 r = 0.27. NDVI follows the PhenoCam closely in the middle of senescence; NIRv agrees poorly with it and should be read as a different measure of canopy state, not a noisier copy.

![EOS90: the same model on three samples (step 42)](../figure/phenocam_vs_satellite/models_EOS90.png)

*EOS90: the same model on three samples (step 42)*

![Which PhenoCam stage corresponds to which satellite NDVI stage (step 42)](../figure/phenocam_vs_satellite/levels_NDVI.png)

*Which PhenoCam stage corresponds to which satellite NDVI stage (step 42)*

![PhenoCam vs satellite NDVI EOS anomalies on shared site-years (step 42)](../figure/phenocam_vs_satellite/scatter_NDVI.png)

*PhenoCam vs satellite NDVI EOS anomalies on shared site-years (step 42)*

## 13. How solid is the pre-solstice effect on EOS90? (step 43)

One-predictor within-site model, two versions of pre-solstice GPP: `sos` = cumulative from leaf-out to the solstice, `cal` = the 60 days before the solstice. Tables: `data/presolstice_loso*.csv`, `data/presolstice_qc_sensitivity.csv`, `data/presolstice_power.csv`.


### 1. Leave one site out

| EOS source | predictor | all sites (p) | range with one site removed | negative in | p < 0.05 in | most influential site | slope without it (p) |
|---|---|---|---|---|---|---|---|
| PhenoCam | sos | -1.2 (0.306) | -1.8 to -0.6 | 100% | 0% | US-Me6 | -1.8 (0.083) |
| PhenoCam | cal | -0.7 (0.509) | -1.0 to -0.1 | 100% | 0% | US-KFS | -0.1 (0.916) |
| tower NDVI | sos | -1.9 (0.151) | -2.2 to -1.1 | 100% | 0% | JP-Tmd | -1.1 (0.353) |
| tower NDVI | cal | -2.5 (0.039) | -2.8 to -1.7 | 100% | 83% | JP-Tmd | -1.7 (0.126) |
| satellite NDVI | sos | -1.9 (0.000) | -2.1 to -1.7 | 100% | 100% | US-NC2 | -1.7 (0.002) |
| satellite NDVI | cal | -1.1 (0.026) | -1.2 to -0.9 | 100% | 96% | US-NC2 | -0.9 (0.085) |
| satellite NIRv | sos | -1.3 (0.005) | -1.5 to -1.1 | 100% | 100% | US-NC2 | -1.1 (0.020) |
| satellite NIRv | cal | -1.5 (0.001) | -1.6 to -1.3 | 100% | 100% | US-NC2 | -1.3 (0.003) |

### 2. QC thresholds

**GPP from leaf-out to the solstice** - slope (site-years)

| min R2 | max EOS90 distance | PhenoCam | tower NDVI | satellite NDVI | satellite NIRv |
|---|---|---|---|---|---|
| 0.7 | 30 d | -0.9 (213) | -0.7 (325) | -1.4** (403) | -1.0* (658) |
| 0.7 | 60 d | -1.0 (225) | -2.7** (372) | -2.0*** (418) | -1.3** (680) |
| 0.7 | 90 d | -1.2 (232) | -2.6* (377) | -2.0*** (418) | -1.4** (684) |
| 0.7 | no rule | -2.1 (236) | -3.0** (381) | -2.2*** (419) | -0.6 (689) |
| 0.8 | 30 d | -1.2 (192) | -0.6 (290) | -1.3** (382) | -1.0* (629) |
| 0.8 (default) | 60 d | -1.2 (205) | -1.9 (318) | -1.9*** (395) | -1.3** (650) |
| 0.8 | 90 d | -1.6 (211) | -1.3 (323) | -1.9*** (395) | -1.4** (653) |
| 0.8 | no rule | -2.2 (213) | -1.6 (324) | -2.2*** (396) | -0.8 (657) |
| 0.9 | 30 d | -1.0 (148) | -0.5 (218) | -1.2** (305) | -1.5*** (519) |
| 0.9 | 60 d | -2.2 (153) | -0.5 (239) | -1.6** (315) | -1.8*** (536) |
| 0.9 | 90 d | -2.3* (154) | -1.0 (240) | -1.6** (315) | -1.9*** (538) |
| 0.9 | no rule | -2.3* (154) | -1.0 (240) | -1.6** (315) | -1.9*** (538) |

**GPP in the 60 days before the solstice** - slope (site-years)

| min R2 | max EOS90 distance | PhenoCam | tower NDVI | satellite NDVI | satellite NIRv |
|---|---|---|---|---|---|
| 0.7 | 30 d | -0.1 (243) | -0.5 (362) | -0.7* (460) | -1.2** (690) |
| 0.7 | 60 d | -0.3 (261) | -2.4* (410) | -1.0* (482) | -1.5** (713) |
| 0.7 | 90 d | -0.1 (268) | -2.8* (416) | -1.0* (482) | -1.6*** (717) |
| 0.7 | no rule | -0.8 (273) | -2.9* (420) | -1.5* (483) | -1.3* (722) |
| 0.8 | 30 d | -0.5 (216) | -1.1 (324) | -0.7* (437) | -1.1** (659) |
| 0.8 (default) | 60 d | -0.7 (233) | -2.5* (357) | -1.1* (456) | -1.5*** (681) |
| 0.8 | 90 d | -0.8 (239) | -2.0 (362) | -1.1* (456) | -1.7*** (684) |
| 0.8 | no rule | -1.5 (244) | -2.1 (363) | -1.6* (457) | -1.2* (688) |
| 0.9 | 30 d | -1.3 (158) | -0.6 (240) | -0.8* (341) | -1.4*** (545) |
| 0.9 | 60 d | -2.2* (163) | -1.2 (261) | -1.1* (353) | -1.6*** (559) |
| 0.9 | 90 d | -2.3* (164) | -1.5 (262) | -1.1* (353) | -1.7*** (561) |
| 0.9 | no rule | -2.3* (164) | -1.5 (262) | -1.1* (353) | -1.7*** (561) |

### 3. Power to detect 1.5 days per SD (80% power, two-sided alpha 0.05)

| sample | EOS source | predictor | site-years (sites) | slope | smallest detectable effect | power now | site-years needed | sites needed |
|---|---|---|---|---|---|---|---|---|
| own site-years | PhenoCam | sos | 205 (26) | -1.2 | 3.2 | 25% | 954 | 121 |
| own site-years | PhenoCam | cal | 233 (30) | -0.7 | 3.0 | 28% | 950 | 122 |
| own site-years | tower NDVI | sos | 318 (37) | -1.9 | 3.7 | 21% | 1889 | 220 |
| own site-years | tower NDVI | cal | 357 (41) | -2.5 | 3.3 | 24% | 1774 | 204 |
| own site-years | satellite NDVI | sos | 395 (62) | -1.9 | 1.5 | 82% | 376 | 59 |
| own site-years | satellite NDVI | cal | 456 (69) | -1.1 | 1.4 | 85% | 400 | 60 |
| own site-years | satellite NIRv | sos | 650 (87) | -1.3 | 1.3 | 88% | 511 | 68 |
| own site-years | satellite NIRv | cal | 681 (88) | -1.5 | 1.2 | 93% | 444 | 57 |
| shared PhenoCam + NDVI | PhenoCam | sos | 81 (13) | -2.0 | 2.8 | 33% | 278 | 45 |
| shared PhenoCam + NDVI | PhenoCam | cal | 99 (16) | -1.3 | 3.0 | 29% | 390 | 63 |
| shared PhenoCam + NDVI | satellite NDVI | sos | 78 (12) | -2.7 | 4.5 | 15% | 717 | 110 |
| shared PhenoCam + NDVI | satellite NDVI | cal | 99 (16) | -1.3 | 4.0 | 18% | 696 | 113 |
| shared PhenoCam + NIRv | PhenoCam | sos | 117 (19) | -1.1 | 2.2 | 48% | 249 | 40 |
| shared PhenoCam + NIRv | PhenoCam | cal | 137 (23) | -0.6 | 2.6 | 37% | 411 | 69 |
| shared PhenoCam + NIRv | satellite NIRv | sos | 128 (22) | -1.1 | 4.2 | 17% | 1024 | 176 |
| shared PhenoCam + NIRv | satellite NIRv | cal | 137 (23) | -0.6 | 4.0 | 18% | 978 | 164 |

**What this shows**

- These checks use the simplest model (EOS90 on pre-solstice GPP alone, no temperature covariate), so the slopes are larger than the headline ones of section 3. The question here is whether the SIGN and rough size depend on a single site or on a threshold choice.
- **Leave one site out.** Over all sources and both predictors the slope stays negative in 100% of the refits. Where the full-sample slope is significant (Tower NDVI cal, Satellite NDVI sos, Satellite NDVI cal, Satellite NIRv sos, Satellite NIRv cal), it stays significant in 83 to 100% of the refits. The most influential sites are JP-Tmd, US-KFS, US-Me6, US-NC2; removing them changes the slope by at most 0.8 days. No single site produces the result.
- **QC thresholds.** Over 12 combinations of minimum fit R2 and the site-median rule, 96 of 96 slopes are negative. By source the range is PhenoCam (GCC) -2.3 to -0.1; Tower NDVI -3.0 to -0.5; Satellite NDVI -2.2 to -0.7; Satellite NIRv -1.9 to -0.6. The satellite slopes move little; the ground-based ones move more because their samples are small and each threshold changes which years are in.
- **Power.** To detect 1.5 days per SD with 80% probability, 376 to 1889 site-years are needed, depending on how noisy the source is. The satellite indices have 395 to 681 and a power of 82 to 93%. PhenoCam and tower NDVI have 205 to 357 site-years and a power of 21 to 28%.
- So a non-significant result from PhenoCam or tower NDVI is not evidence against the effect: with their sample sizes they would miss it three times out of four. Only the satellite indices can confirm or reject an effect of this size; the ground sources can only agree or disagree in sign.

![Slope with each site left out in turn (step 43)](../figure/presolstice_robustness/loso.png)

*Slope with each site left out in turn (step 43)*

![Slope under other phenology QC thresholds (step 43)](../figure/presolstice_robustness/qc_sensitivity.png)

*Slope under other phenology QC thresholds (step 43)*

## 14. Why? Leaf-out, water and sink variables (step 44)

Tables: `data/mechanism_leafout_vs_gpp.csv`, `data/mechanism_water.csv`, `data/mechanism_sink.csv`.


### A. Leaf-out date against early GPP

| sites | EOS source | early GPP | n (sites) | r(GPP, leaf-out) | GPP alone | leaf-out alone | GPP, leaf-out fixed | leaf-out, GPP fixed | R2 only GPP | R2 only leaf-out |
|---|---|---|---|---|---|---|---|---|---|---|
| deciduous | PhenoCam | sos | 105 (9) | -0.45 | -0.4 | +1.6 | +0.4 | +1.7 | 0.001 | 0.014 |
| deciduous | PhenoCam | cal | 104 (9) | -0.28 | -0.2 | +1.6 | +0.3 | +1.6 | 0.000 | 0.014 |
| deciduous | PhenoCam | rate | 105 (9) | +0.31 | +0.5 | +1.6 | +0.0 | +1.6 | 0.000 | 0.012 |
| deciduous | tower NDVI | sos | 171 (17) | -0.35 | +0.3 | +0.8 | +0.7 | +1.0 | 0.003 | 0.006 |
| deciduous | tower NDVI | cal | 174 (17) | -0.16 | -0.4 | +0.9 | -0.3 | +0.8 | 0.000 | 0.004 |
| deciduous | tower NDVI | rate | 171 (17) | +0.51 | +0.8 | +0.8 | +0.5 | +0.5 | 0.001 | 0.001 |
| deciduous | satellite NDVI | sos | 173 (23) | -0.46 | -1.4* | +3.6*** | +0.4 | +3.8*** | 0.001 | 0.131 |
| deciduous | satellite NDVI | cal | 172 (23) | -0.26 | -0.6 | +3.6*** | +0.4 | +3.7*** | 0.001 | 0.147 |
| deciduous | satellite NDVI | rate | 173 (23) | +0.23 | +1.1 | +3.6*** | +0.3 | +3.5*** | 0.001 | 0.137 |
| deciduous | satellite NIRv | sos | 198 (25) | -0.49 | -0.9 | +0.4 | -0.9 | -0.0 | 0.005 | 0.000 |
| deciduous | satellite NIRv | cal | 197 (25) | -0.24 | -1.0 | +0.4 | -1.0 | +0.2 | 0.008 | 0.000 |
| deciduous | satellite NIRv | rate | 198 (25) | +0.21 | -0.5 | +0.4 | -0.7 | +0.6 | 0.004 | 0.003 |
| deciduous | all sources stacked | sos | 647 (29) | -0.41 | -0.7 | +1.4** | -0.1 | +1.4* | 0.000 | 0.013 |
| deciduous | all sources stacked | cal | 647 (29) | -0.21 | -0.6 | +1.4** | -0.3 | +1.4** | 0.001 | 0.014 |
| deciduous | all sources stacked | rate | 647 (29) | +0.32 | +0.4 | +1.4** | -0.0 | +1.4** | 0.000 | 0.014 |
| evergreen | tower NDVI | sos | 81 (12) | -0.48 | -2.7 | -0.9 | -4.1* | -2.9 | 0.039 | 0.019 |
| evergreen | tower NDVI | cal | 87 (12) | -0.10 | -4.6** | +0.9 | -4.5*** | +0.4 | 0.062 | 0.001 |
| evergreen | tower NDVI | rate | 81 (12) | +0.75 | -4.2* | -0.9 | -8.2*** | +5.3* | 0.088 | 0.037 |
| evergreen | satellite NDVI | sos | 82 (16) | -0.21 | -2.5*** | +0.9 | -2.5*** | +0.4 | 0.031 | 0.001 |
| evergreen | satellite NDVI | cal | 88 (17) | -0.05 | -1.4 | +0.8 | -1.4 | +0.8 | 0.008 | 0.003 |
| evergreen | satellite NDVI | rate | 82 (16) | +0.31 | -1.4 | +0.9 | -1.8 | +1.5 | 0.016 | 0.011 |
| evergreen | satellite NIRv | sos | 210 (26) | -0.47 | -1.7* | -0.8 | -2.6*** | -2.0* | 0.043 | 0.025 |
| evergreen | satellite NIRv | cal | 213 (27) | +0.06 | -2.0*** | -0.6 | -2.0** | -0.5 | 0.031 | 0.002 |
| evergreen | satellite NIRv | rate | 210 (26) | +0.34 | -2.0** | -0.8 | -1.9** | -0.1 | 0.026 | 0.000 |
| evergreen | all sources stacked | sos | 384 (28) | -0.39 | -1.7* | -0.9 | -2.5*** | -1.9 | 0.028 | 0.017 |
| evergreen | all sources stacked | cal | 399 (28) | -0.01 | -2.0*** | -0.2 | -2.0*** | -0.2 | 0.021 | 0.000 |
| evergreen | all sources stacked | rate | 384 (28) | +0.41 | -2.2** | -0.9 | -2.1** | -0.1 | 0.021 | 0.000 |
| grass/shrub | PhenoCam | sos | 89 (14) | -0.45 | -3.6** | +0.7 | -4.1* | -1.1 | 0.098 | 0.008 |
| grass/shrub | PhenoCam | cal | 118 (18) | -0.08 | -1.9 | +1.4 | -1.9 | +1.2 | 0.021 | 0.009 |
| grass/shrub | PhenoCam | rate | 89 (14) | +0.27 | -2.7 | +0.7 | -3.2 | +1.6 | 0.069 | 0.017 |
| grass/shrub | tower NDVI | sos | 66 (8) | -0.25 | -6.5 | +3.5 | -6.0 | +2.0 | 0.096 | 0.011 |
| grass/shrub | tower NDVI | cal | 96 (12) | -0.15 | -5.3* | -2.3 | -5.8** | -3.2* | 0.097 | 0.029 |
| grass/shrub | tower NDVI | rate | 66 (8) | +0.30 | -4.0 | +3.5 | -5.6* | +5.2* | 0.082 | 0.071 |
| grass/shrub | satellite NDVI | sos | 140 (23) | -0.60 | -2.4 | +1.4 | -2.4 | -0.1 | 0.039 | 0.000 |
| grass/shrub | satellite NDVI | cal | 196 (29) | -0.17 | -1.7 | -0.7 | -1.9* | -1.0 | 0.024 | 0.007 |
| grass/shrub | satellite NDVI | rate | 140 (23) | +0.36 | +0.0 | +1.4 | -0.5 | +1.5 | 0.002 | 0.021 |
| grass/shrub | satellite NIRv | sos | 242 (36) | -0.48 | -1.4 | -1.4* | -2.7** | -2.7*** | 0.029 | 0.028 |
| grass/shrub | satellite NIRv | cal | 271 (36) | -0.18 | -1.5 | -0.9 | -1.7 | -1.2 | 0.014 | 0.008 |
| grass/shrub | satellite NIRv | rate | 242 (36) | +0.07 | -1.0 | -1.4* | -0.9 | -1.3 | 0.004 | 0.009 |
| grass/shrub | all sources stacked | sos | 537 (37) | -0.46 | -2.7** | -0.1 | -3.5*** | -1.7** | 0.054 | 0.013 |
| grass/shrub | all sources stacked | cal | 681 (37) | -0.16 | -2.2* | -0.8 | -2.4** | -1.2 | 0.029 | 0.007 |
| grass/shrub | all sources stacked | rate | 537 (37) | +0.17 | -1.4* | -0.1 | -1.5* | +0.2 | 0.012 | 0.000 |
| all | PhenoCam | sos | 205 (26) | -0.43 | -1.2 | -0.1 | -1.5 | -0.7 | 0.011 | 0.003 |
| all | PhenoCam | cal | 233 (30) | -0.14 | -0.7 | +0.4 | -0.7 | +0.3 | 0.003 | 0.001 |
| all | PhenoCam | rate | 205 (26) | +0.26 | -0.7 | -0.1 | -0.8 | +0.1 | 0.003 | 0.000 |
| all | tower NDVI | sos | 318 (37) | -0.35 | -1.9 | +0.6 | -1.9 | -0.1 | 0.014 | 0.000 |
| all | tower NDVI | cal | 357 (41) | -0.13 | -2.5* | -0.3 | -2.5* | -0.6 | 0.026 | 0.002 |
| all | tower NDVI | rate | 318 (37) | +0.51 | -1.5 | +0.6 | -2.4 | +1.8 | 0.019 | 0.011 |
| all | satellite NDVI | sos | 395 (62) | -0.44 | -1.9*** | +2.0* | -1.2* | +1.5 | 0.011 | 0.016 |
| all | satellite NDVI | cal | 456 (69) | -0.14 | -1.1* | +0.4 | -1.1* | +0.3 | 0.008 | 0.001 |
| all | satellite NDVI | rate | 395 (62) | +0.28 | +0.2 | +2.0* | -0.4 | +2.1* | 0.001 | 0.038 |
| all | satellite NIRv | sos | 650 (87) | -0.42 | -1.3** | -0.8 | -2.0*** | -1.7*** | 0.023 | 0.016 |
| all | satellite NIRv | cal | 681 (88) | -0.10 | -1.5*** | -0.6 | -1.6*** | -0.8 | 0.016 | 0.004 |
| all | satellite NIRv | rate | 650 (87) | +0.16 | -1.2** | -0.8 | -1.1* | -0.7 | 0.007 | 0.003 |
| all | all sources stacked | sos | 1568 (94) | -0.40 | -1.6*** | +0.1 | -1.8*** | -0.6 | 0.017 | 0.002 |
| all | all sources stacked | cal | 1727 (94) | -0.12 | -1.5*** | -0.2 | -1.5*** | -0.4 | 0.014 | 0.001 |
| all | all sources stacked | rate | 1568 (94) | +0.27 | -0.8* | +0.1 | -0.9* | +0.4 | 0.005 | 0.001 |

### B. Water (90 days before the solstice) against pre-solstice GPP

Plant-available soil water = measured soil water content relative to the site's own dry and wet ends.

| sites | EOS source | water measure | GPP | n (sites) | r(GPP, water) | GPP alone | water alone | GPP, water fixed | water, GPP fixed | GPP, water + T fixed | R2 only GPP | R2 only water |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| dry-summer | PhenoCam | water balance (P - PET) | sos | 31 (5) | +0.30 | -0.9 | +2.4 | -1.7 | +2.9 | -1.7 | 0.032 | 0.089 |
| dry-summer | PhenoCam | water balance (P - PET) | cal | 31 (5) | +0.23 | +3.0 | +2.4 | +2.6 | +1.8 | +2.9 | 0.076 | 0.034 |
| dry-summer | satellite NIRv | water balance (P - PET) | sos | 69 (9) | +0.04 | -1.2 | +1.6*** | -1.3 | +1.7** | -1.3 | 0.015 | 0.026 |
| dry-summer | satellite NIRv | water balance (P - PET) | cal | 70 (9) | +0.05 | -0.5 | +1.6** | -0.5 | +1.6*** | -0.5 | 0.003 | 0.023 |
| dry-summer | satellite NIRv | plant-available soil water | sos | 37 (6) | +0.25 | +0.5 | -0.1 | +0.5 | -0.2 | +0.4 | 0.005 | 0.001 |
| dry-summer | satellite NIRv | plant-available soil water | cal | 37 (6) | +0.41 | +0.2 | -0.1 | +0.3 | -0.2 | -0.1 | 0.002 | 0.000 |
| dry-summer | satellite NIRv | plant-available soil water, after the solstice | sos | 36 (5) | +0.16 | +1.0 | -0.8 | +1.2 | -1.0 | +1.0 | 0.018 | 0.012 |
| dry-summer | satellite NIRv | plant-available soil water, after the solstice | cal | 36 (5) | -0.06 | +0.4 | -0.8 | +0.3 | -0.7 | -0.0 | 0.002 | 0.008 |
| dry-summer | all sources stacked | water balance (P - PET) | sos | 132 (9) | +0.15 | -0.2 | +2.0** | -0.5 | +2.0* | -0.5 | 0.002 | 0.040 |
| dry-summer | all sources stacked | water balance (P - PET) | cal | 136 (9) | +0.09 | +0.5 | +1.8** | +0.3 | +1.8** | +0.4 | 0.001 | 0.033 |
| dry-summer | all sources stacked | plant-available soil water | sos | 81 (6) | +0.23 | +0.9 | +1.4 | +0.6 | +1.3 | +0.5 | 0.004 | 0.020 |
| dry-summer | all sources stacked | plant-available soil water | cal | 81 (6) | +0.50 | +1.8 | +1.4 | +1.5 | +0.7 | +1.2 | 0.020 | 0.005 |
| dry-summer | all sources stacked | plant-available soil water, after the solstice | sos | 79 (6) | +0.13 | +1.0 | -0.1 | +1.0 | -0.3 | +1.0 | 0.012 | 0.001 |
| dry-summer | all sources stacked | plant-available soil water, after the solstice | cal | 79 (6) | +0.08 | +1.9 | -0.1 | +1.9 | -0.3 | +1.7 | 0.040 | 0.001 |
| summer-green | PhenoCam | water balance (P - PET) | sos | 174 (21) | -0.24 | -1.2 | +1.2 | -1.0 | +0.9 | -0.8 | 0.006 | 0.005 |
| summer-green | PhenoCam | water balance (P - PET) | cal | 202 (25) | -0.24 | -1.2 | +1.2 | -1.0 | +1.0 | -0.7 | 0.005 | 0.005 |
| summer-green | PhenoCam | plant-available soil water | sos | 140 (19) | -0.01 | -0.9 | -1.3 | -0.9 | -1.3 | -0.1 | 0.004 | 0.009 |
| summer-green | PhenoCam | plant-available soil water | cal | 160 (22) | -0.02 | -0.8 | -1.0 | -0.8 | -1.0 | -0.2 | 0.004 | 0.005 |
| summer-green | PhenoCam | plant-available soil water, after the solstice | sos | 141 (19) | -0.24 | -0.9 | +1.0 | -0.6 | +0.8 | -0.3 | 0.002 | 0.004 |
| summer-green | PhenoCam | plant-available soil water, after the solstice | cal | 161 (22) | -0.21 | -0.8 | +0.7 | -0.7 | +0.5 | -0.3 | 0.003 | 0.001 |
| summer-green | tower NDVI | water balance (P - PET) | sos | 303 (34) | -0.17 | -2.2 | +0.9 | -2.1 | +0.6 | -1.3 | 0.019 | 0.001 |
| summer-green | tower NDVI | water balance (P - PET) | cal | 342 (38) | -0.19 | -2.8* | +1.0 | -2.7* | +0.5 | -2.0 | 0.029 | 0.001 |
| summer-green | tower NDVI | plant-available soil water | sos | 250 (31) | -0.03 | -2.3 | +0.4 | -2.3 | +0.3 | -1.1 | 0.022 | 0.000 |
| summer-green | tower NDVI | plant-available soil water | cal | 288 (35) | +0.02 | -3.0* | +0.4 | -3.1* | +0.5 | -2.1 | 0.037 | 0.001 |
| summer-green | tower NDVI | plant-available soil water, after the solstice | sos | 254 (31) | -0.11 | -2.4 | +0.4 | -2.4 | +0.1 | -1.3 | 0.024 | 0.000 |
| summer-green | tower NDVI | plant-available soil water, after the solstice | cal | 293 (35) | -0.11 | -3.1* | +0.4 | -3.1* | +0.0 | -2.2 | 0.038 | 0.000 |
| summer-green | satellite NDVI | water balance (P - PET) | sos | 378 (58) | +0.04 | -2.1*** | -0.4 | -2.1*** | -0.3 | -2.0*** | 0.038 | 0.001 |
| summer-green | satellite NDVI | water balance (P - PET) | cal | 436 (65) | +0.01 | -1.1* | -0.4 | -1.1* | -0.3 | -0.7 | 0.008 | 0.001 |
| summer-green | satellite NDVI | plant-available soil water | sos | 309 (49) | -0.09 | -2.2*** | +0.2 | -2.2*** | +0.0 | -2.2*** | 0.047 | 0.000 |
| summer-green | satellite NDVI | plant-available soil water | cal | 346 (54) | -0.06 | -1.1* | +0.4 | -1.1* | +0.3 | -0.8 | 0.010 | 0.001 |
| summer-green | satellite NDVI | plant-available soil water, after the solstice | sos | 310 (49) | -0.15 | -2.2*** | +0.9 | -2.1*** | +0.6 | -2.1*** | 0.043 | 0.003 |
| summer-green | satellite NDVI | plant-available soil water, after the solstice | cal | 347 (54) | -0.17 | -1.1* | +0.7 | -1.0* | +0.6 | -0.8 | 0.008 | 0.002 |
| summer-green | satellite NIRv | water balance (P - PET) | sos | 581 (78) | +0.07 | -1.3* | -0.1 | -1.3** | +0.0 | -1.0 | 0.012 | 0.000 |
| summer-green | satellite NIRv | water balance (P - PET) | cal | 611 (79) | +0.03 | -1.6*** | -0.1 | -1.6*** | -0.0 | -1.2** | 0.017 | 0.000 |
| summer-green | satellite NIRv | plant-available soil water | sos | 456 (68) | +0.01 | -1.8*** | +0.4 | -1.8*** | +0.4 | -1.2* | 0.026 | 0.001 |
| summer-green | satellite NIRv | plant-available soil water | cal | 480 (69) | +0.03 | -1.9*** | +0.5 | -1.9*** | +0.5 | -1.2* | 0.029 | 0.002 |
| summer-green | satellite NIRv | plant-available soil water, after the solstice | sos | 458 (68) | -0.08 | -1.8** | +0.4 | -1.8** | +0.2 | -1.2* | 0.024 | 0.000 |
| summer-green | satellite NIRv | plant-available soil water, after the solstice | cal | 482 (69) | -0.10 | -1.9*** | +0.3 | -1.8*** | +0.2 | -1.2* | 0.026 | 0.000 |
| summer-green | all sources stacked | water balance (P - PET) | sos | 1436 (85) | -0.03 | -1.7*** | +0.3 | -1.7*** | +0.2 | -1.3** | 0.018 | 0.000 |
| summer-green | all sources stacked | water balance (P - PET) | cal | 1591 (85) | -0.06 | -1.6*** | +0.3 | -1.6*** | +0.2 | -1.2** | 0.015 | 0.000 |
| summer-green | all sources stacked | plant-available soil water | sos | 1155 (74) | -0.03 | -1.9*** | +0.2 | -1.9*** | +0.1 | -1.3* | 0.024 | 0.000 |
| summer-green | all sources stacked | plant-available soil water | cal | 1274 (74) | -0.00 | -1.8*** | +0.3 | -1.8*** | +0.3 | -1.2* | 0.020 | 0.000 |
| summer-green | all sources stacked | plant-available soil water, after the solstice | sos | 1163 (74) | -0.13 | -1.9*** | +0.6 | -1.9*** | +0.4 | -1.3* | 0.023 | 0.001 |
| summer-green | all sources stacked | plant-available soil water, after the solstice | cal | 1283 (74) | -0.14 | -1.8*** | +0.5 | -1.8*** | +0.3 | -1.2* | 0.019 | 0.000 |
| all | PhenoCam | water balance (P - PET) | sos | 205 (26) | -0.18 | -1.2 | +1.3 | -1.0 | +1.1 | -0.8 | 0.006 | 0.008 |
| all | PhenoCam | water balance (P - PET) | cal | 233 (30) | -0.19 | -0.7 | +1.3 | -0.5 | +1.2 | -0.3 | 0.001 | 0.008 |
| all | PhenoCam | plant-available soil water | sos | 164 (22) | +0.02 | -1.0 | -0.7 | -0.9 | -0.7 | -0.4 | 0.005 | 0.003 |
| all | PhenoCam | plant-available soil water | cal | 184 (25) | +0.02 | -0.2 | -0.5 | -0.2 | -0.5 | +0.3 | 0.000 | 0.002 |
| all | PhenoCam | plant-available soil water, after the solstice | sos | 165 (22) | -0.22 | -1.0 | +0.9 | -0.8 | +0.7 | -0.5 | 0.004 | 0.003 |
| all | PhenoCam | plant-available soil water, after the solstice | cal | 185 (25) | -0.20 | -0.2 | +0.6 | -0.1 | +0.6 | +0.2 | 0.000 | 0.002 |
| all | tower NDVI | water balance (P - PET) | sos | 318 (37) | -0.14 | -1.9 | +1.0 | -1.8 | +0.7 | -1.0 | 0.013 | 0.002 |
| all | tower NDVI | water balance (P - PET) | cal | 357 (41) | -0.17 | -2.5* | +1.1 | -2.4* | +0.7 | -1.6 | 0.022 | 0.002 |
| all | tower NDVI | plant-available soil water | sos | 262 (33) | -0.02 | -1.9 | +0.5 | -1.9 | +0.4 | -0.8 | 0.016 | 0.001 |
| all | tower NDVI | plant-available soil water | cal | 300 (37) | +0.03 | -2.7* | +0.5 | -2.8 | +0.6 | -1.7 | 0.031 | 0.001 |
| all | tower NDVI | plant-available soil water, after the solstice | sos | 265 (33) | -0.11 | -2.1 | +0.5 | -2.0 | +0.2 | -1.0 | 0.018 | 0.000 |
| all | tower NDVI | plant-available soil water, after the solstice | cal | 304 (37) | -0.11 | -2.9* | +0.4 | -2.8* | +0.1 | -1.9 | 0.032 | 0.000 |
| all | satellite NDVI | water balance (P - PET) | sos | 395 (62) | +0.04 | -1.9*** | -0.3 | -1.9*** | -0.2 | -1.8*** | 0.032 | 0.000 |
| all | satellite NDVI | water balance (P - PET) | cal | 456 (69) | +0.01 | -1.1* | -0.3 | -1.1* | -0.3 | -0.8 | 0.009 | 0.001 |
| all | satellite NDVI | plant-available soil water | sos | 317 (51) | -0.08 | -2.0*** | +0.2 | -2.0*** | +0.0 | -1.9*** | 0.038 | 0.000 |
| all | satellite NDVI | plant-available soil water | cal | 354 (56) | -0.04 | -1.1* | +0.4 | -1.1* | +0.3 | -0.9 | 0.010 | 0.001 |
| all | satellite NDVI | plant-available soil water, after the solstice | sos | 318 (51) | -0.14 | -2.0*** | +0.9 | -1.9*** | +0.6 | -1.9*** | 0.035 | 0.004 |
| all | satellite NDVI | plant-available soil water, after the solstice | cal | 355 (56) | -0.15 | -1.1* | +0.7 | -1.1* | +0.5 | -0.8 | 0.009 | 0.002 |
| all | satellite NIRv | water balance (P - PET) | sos | 650 (87) | +0.07 | -1.3** | +0.1 | -1.3** | +0.2 | -1.1* | 0.012 | 0.000 |
| all | satellite NIRv | water balance (P - PET) | cal | 681 (88) | +0.03 | -1.5*** | +0.1 | -1.5*** | +0.1 | -1.2** | 0.015 | 0.000 |
| all | satellite NIRv | plant-available soil water | sos | 493 (74) | +0.02 | -1.6** | +0.3 | -1.7** | +0.4 | -1.1* | 0.022 | 0.001 |
| all | satellite NIRv | plant-available soil water | cal | 517 (75) | +0.05 | -1.8*** | +0.4 | -1.8*** | +0.5 | -1.1* | 0.026 | 0.002 |
| all | satellite NIRv | plant-available soil water, after the solstice | sos | 494 (73) | -0.07 | -1.6** | +0.3 | -1.6** | +0.2 | -1.0* | 0.019 | 0.000 |
| all | satellite NIRv | plant-available soil water, after the solstice | cal | 518 (74) | -0.10 | -1.7*** | +0.3 | -1.7*** | +0.1 | -1.1* | 0.022 | 0.000 |
| all | all sources stacked | water balance (P - PET) | sos | 1568 (94) | -0.02 | -1.6*** | +0.4 | -1.6*** | +0.3 | -1.2** | 0.016 | 0.001 |
| all | all sources stacked | water balance (P - PET) | cal | 1727 (94) | -0.05 | -1.5*** | +0.4 | -1.5*** | +0.3 | -1.1** | 0.013 | 0.001 |
| all | all sources stacked | plant-available soil water | sos | 1236 (80) | -0.01 | -1.7** | +0.2 | -1.7** | +0.2 | -1.2* | 0.020 | 0.000 |
| all | all sources stacked | plant-available soil water | cal | 1355 (80) | +0.02 | -1.6** | +0.3 | -1.6** | +0.3 | -1.0* | 0.016 | 0.001 |
| all | all sources stacked | plant-available soil water, after the solstice | sos | 1242 (80) | -0.11 | -1.7** | +0.6 | -1.7** | +0.4 | -1.2* | 0.019 | 0.001 |
| all | all sources stacked | plant-available soil water, after the solstice | cal | 1362 (80) | -0.13 | -1.6** | +0.5 | -1.6** | +0.3 | -1.0* | 0.015 | 0.000 |

### C. Sink variables in place of GPP (60 days before / 45 days from the solstice)

| sites | EOS source | variable | n (sites) | pre | post | R2 | R2 of GPP, same rows | r with GPP (pre) | variable next to GPP | GPP next to it |
|---|---|---|---|---|---|---|---|---|---|---|
| all | PhenoCam | GPP | 233 (30) | -0.8 | +0.3 | 0.003 | - | - | - | - |
| all | PhenoCam | NPPd | 231 (30) | -1.0 | +1.1 | 0.010 | 0.004 | +0.67 | -0.5 | -0.5 |
| all | PhenoCam | Ra | 231 (30) | +0.2 | -0.9 | 0.004 | 0.004 | +0.52 | +0.4 | -1.0 |
| all | PhenoCam | Rg | 231 (30) | -1.9 | +1.1 | 0.026 | 0.004 | +0.66 | -2.3 | +0.7 |
| all | PhenoCam | Rm | 231 (30) | +0.6 | -1.1 | 0.007 | 0.004 | +0.36 | +0.6 | -1.0 |
| all | PhenoCam | Rm_model | 233 (30) | +0.2 | -0.4 | 0.000 | 0.003 | +0.46 | +0.2 | -0.8 |
| all | PhenoCam | SI | 256 (33) | +0.6 | -0.1 | 0.002 | 0.003 | -0.06 | +0.6 | -0.7 |
| all | PhenoCam | fT | 256 (33) | -1.7* | +0.4 | 0.016 | 0.003 | +0.18 | -1.7* | -0.4 |
| all | PhenoCam | fW | 256 (33) | +1.6 | -0.5 | 0.012 | 0.003 | -0.19 | +1.4 | -0.5 |
| all | tower NDVI | GPP | 356 (41) | -2.8* | +0.8 | 0.028 | - | - | - | - |
| all | tower NDVI | NPPd | 356 (41) | -1.7 | +0.8 | 0.011 | 0.028 | +0.68 | +0.6 | -2.9* |
| all | tower NDVI | Ra | 356 (41) | -1.3 | +0.1 | 0.007 | 0.028 | +0.37 | -0.5 | -2.3 |
| all | tower NDVI | Rg | 356 (41) | -2.0 | +0.1 | 0.016 | 0.028 | +0.65 | -0.6 | -2.1 |
| all | tower NDVI | Rm | 356 (41) | -0.8 | +0.0 | 0.003 | 0.028 | +0.21 | -0.3 | -2.4* |
| all | tower NDVI | Rm_model | 356 (41) | -3.1 | +1.8 | 0.018 | 0.028 | +0.50 | -0.7 | -2.1 |
| all | tower NDVI | SI | 385 (43) | -0.8 | +0.5 | 0.003 | 0.028 | +0.31 | -0.1 | -2.5* |
| all | tower NDVI | fT | 385 (43) | -2.9** | +0.1 | 0.036 | 0.028 | +0.41 | -2.8** | -1.3 |
| all | tower NDVI | fW | 385 (43) | +1.6 | -0.7 | 0.009 | 0.028 | -0.07 | +1.4 | -2.4* |
| all | satellite NDVI | GPP | 446 (67) | -1.2 | -0.0 | 0.010 | - | - | - | - |
| all | satellite NDVI | NPPd | 441 (66) | -1.0 | -0.4 | 0.010 | 0.012 | +0.78 | -0.2 | -1.1 |
| all | satellite NDVI | Ra | 441 (66) | -0.5 | +0.5 | 0.002 | 0.012 | +0.39 | +0.2 | -1.4* |
| all | satellite NDVI | Rg | 441 (66) | -0.7 | -1.1 | 0.014 | 0.012 | +0.76 | +0.1 | -1.4 |
| all | satellite NDVI | Rm | 441 (66) | -0.2 | +0.6 | 0.002 | 0.012 | +0.17 | +0.1 | -1.3* |
| all | satellite NDVI | Rm_model | 446 (67) | -1.6 | +1.5 | 0.006 | 0.010 | +0.48 | +0.3 | -1.3 |
| all | satellite NDVI | SI | 497 (73) | -1.1* | -0.1 | 0.010 | 0.010 | +0.18 | -1.0* | -1.0* |
| all | satellite NDVI | fT | 497 (73) | -1.3* | +0.4 | 0.014 | 0.010 | +0.24 | -1.0 | -0.9 |
| all | satellite NDVI | fW | 497 (73) | -0.1 | -0.3 | 0.001 | 0.010 | -0.02 | -0.3 | -1.2* |
| all | satellite NIRv | GPP | 675 (88) | -1.9*** | +1.1 | 0.022 | - | - | - | - |
| all | satellite NIRv | NPPd | 672 (88) | -1.2* | +0.8 | 0.011 | 0.023 | +0.69 | +0.1 | -1.6* |
| all | satellite NIRv | Ra | 672 (88) | -0.7 | +0.2 | 0.003 | 0.023 | +0.41 | -0.1 | -1.5*** |
| all | satellite NIRv | Rg | 672 (88) | -1.1* | +0.5 | 0.008 | 0.023 | +0.68 | +0.2 | -1.6* |
| all | satellite NIRv | Rm | 672 (88) | -0.4 | +0.1 | 0.001 | 0.023 | +0.24 | -0.1 | -1.5*** |
| all | satellite NIRv | Rm_model | 675 (88) | -0.7 | +0.5 | 0.001 | 0.022 | +0.42 | +0.4 | -1.7*** |
| all | satellite NIRv | SI | 751 (94) | -1.5** | +0.4 | 0.015 | 0.022 | +0.23 | -1.3* | -1.2** |
| all | satellite NIRv | fT | 751 (94) | -2.0*** | -0.5 | 0.030 | 0.022 | +0.17 | -1.8** | -1.2** |
| all | satellite NIRv | fW | 751 (94) | -0.1 | +0.3 | 0.000 | 0.022 | +0.08 | -0.1 | -1.5*** |
| deciduous | PhenoCam | GPP | 104 (9) | +0.4 | -1.7 | 0.014 | - | - | - | - |
| deciduous | PhenoCam | NPPd | 104 (9) | -1.0 | +1.2 | 0.009 | 0.014 | +0.65 | -0.9 | +0.4 |
| deciduous | PhenoCam | Ra | 104 (9) | +1.3 | -2.5* | 0.034 | 0.014 | +0.46 | +0.8 | -0.6 |
| deciduous | PhenoCam | Rg | 104 (9) | -2.3 | +1.6 | 0.041 | 0.014 | +0.59 | -3.2 | +1.7 |
| deciduous | PhenoCam | Rm | 104 (9) | +1.6 | -2.6** | 0.041 | 0.014 | +0.31 | +1.2 | -0.6 |
| deciduous | PhenoCam | Rm_model | 104 (9) | +3.6 | -5.2 | 0.034 | 0.014 | +0.54 | -1.5 | +0.6 |
| deciduous | PhenoCam | SI | 116 (10) | -1.8 | +1.3 | 0.020 | 0.014 | -0.00 | -1.5 | -0.2 |
| deciduous | PhenoCam | fT | 116 (10) | -0.3 | +1.3 | 0.008 | 0.014 | +0.36 | +0.1 | -0.2 |
| deciduous | PhenoCam | fW | 116 (10) | -1.5 | +1.0 | 0.011 | 0.014 | -0.30 | -1.4 | -0.6 |
| deciduous | tower NDVI | GPP | 174 (17) | -1.2 | +3.1* | 0.059 | - | - | - | - |
| deciduous | tower NDVI | NPPd | 174 (17) | -0.9 | +2.7** | 0.047 | 0.059 | +0.67 | +0.0 | -0.4 |
| deciduous | tower NDVI | Ra | 174 (17) | -0.1 | -0.0 | 0.000 | 0.059 | +0.33 | -0.0 | -0.4 |
| deciduous | tower NDVI | Rg | 174 (17) | -1.8 | +2.2* | 0.040 | 0.059 | +0.63 | -1.6 | +0.6 |
| deciduous | tower NDVI | Rm | 174 (17) | +0.3 | -0.4 | 0.001 | 0.059 | +0.17 | +0.2 | -0.4 |
| deciduous | tower NDVI | Rm_model | 174 (17) | +0.4 | -0.3 | 0.000 | 0.059 | +0.56 | +0.5 | -0.7 |
| deciduous | tower NDVI | SI | 189 (18) | +0.5 | +2.5 | 0.048 | 0.059 | +0.20 | +1.4 | -0.7 |
| deciduous | tower NDVI | fT | 189 (18) | -1.2 | +1.2 | 0.017 | 0.059 | +0.40 | -1.8 | +0.3 |
| deciduous | tower NDVI | fW | 189 (18) | +1.5 | +1.3 | 0.038 | 0.059 | -0.17 | +2.6** | +0.0 |
| deciduous | satellite NDVI | GPP | 170 (23) | -1.1 | +1.0 | 0.013 | - | - | - | - |
| deciduous | satellite NDVI | NPPd | 169 (23) | -0.6 | +1.2 | 0.016 | 0.013 | +0.76 | +0.4 | -0.9 |
| deciduous | satellite NDVI | Ra | 169 (23) | -0.4 | -0.6 | 0.006 | 0.013 | +0.43 | -0.3 | -0.5 |
| deciduous | satellite NDVI | Rg | 169 (23) | -0.2 | -0.1 | 0.001 | 0.013 | +0.72 | +0.3 | -0.8 |
| deciduous | satellite NDVI | Rm | 169 (23) | -0.3 | -0.5 | 0.005 | 0.013 | +0.23 | -0.3 | -0.5 |
| deciduous | satellite NDVI | Rm_model | 170 (23) | -1.7 | +0.6 | 0.015 | 0.013 | +0.61 | -1.2 | +0.1 |
| deciduous | satellite NDVI | SI | 194 (26) | -2.0*** | +0.9 | 0.050 | 0.013 | +0.16 | -1.6** | -0.3 |
| deciduous | satellite NDVI | fT | 194 (26) | -1.0 | +0.3 | 0.013 | 0.013 | +0.38 | -0.4 | -0.4 |
| deciduous | satellite NDVI | fW | 194 (26) | -1.4 | +1.1* | 0.022 | 0.013 | -0.18 | -1.2 | -0.8 |
| deciduous | satellite NIRv | GPP | 195 (25) | -2.1* | +2.3 | 0.044 | - | - | - | - |
| deciduous | satellite NIRv | NPPd | 194 (25) | -1.0 | +2.0** | 0.036 | 0.044 | +0.73 | +0.4 | -1.4 |
| deciduous | satellite NIRv | Ra | 194 (25) | -0.6 | -0.4 | 0.007 | 0.044 | +0.47 | -0.3 | -0.9 |
| deciduous | satellite NIRv | Rg | 194 (25) | -1.6** | +1.5 | 0.030 | 0.044 | +0.72 | -0.9 | -0.4 |
| deciduous | satellite NIRv | Rm | 194 (25) | -0.3 | -0.7 | 0.006 | 0.044 | +0.30 | -0.2 | -1.0 |
| deciduous | satellite NIRv | Rm_model | 195 (25) | -2.0 | +0.9 | 0.015 | 0.044 | +0.57 | -0.9 | -0.5 |
| deciduous | satellite NIRv | SI | 224 (28) | -0.9 | +1.2 | 0.015 | 0.044 | +0.23 | -0.3 | -1.0 |
| deciduous | satellite NIRv | fT | 224 (28) | -2.7*** | -0.3 | 0.064 | 0.044 | +0.39 | -2.7** | -0.0 |
| deciduous | satellite NIRv | fW | 224 (28) | +1.6 | +0.4 | 0.028 | 0.044 | -0.10 | +1.5 | -0.9 |
| evergreen | tower NDVI | GPP | 86 (12) | -5.1* | +0.7 | 0.070 | - | - | - | - |
| evergreen | tower NDVI | NPPd | 86 (12) | -0.9 | -1.6 | 0.014 | 0.070 | +0.52 | +2.0 | -6.0*** |
| evergreen | tower NDVI | Ra | 86 (12) | -3.5 | +1.5 | 0.035 | 0.070 | +0.44 | -1.7 | -4.2* |
| evergreen | tower NDVI | Rg | 86 (12) | -0.1 | -2.3 | 0.016 | 0.070 | +0.47 | +2.7 | -6.3*** |
| evergreen | tower NDVI | Rm | 86 (12) | -3.2 | +1.9 | 0.031 | 0.070 | +0.32 | -1.8 | -4.4* |
| evergreen | tower NDVI | Rm_model | 86 (12) | -8.5* | +3.6 | 0.140 | 0.070 | +0.44 | -5.2* | -3.0 |
| evergreen | tower NDVI | SI | 89 (12) | -3.9 | +0.8 | 0.043 | 0.070 | +0.47 | -2.3 | -3.8*** |
| evergreen | tower NDVI | fT | 89 (12) | -5.6*** | +1.0 | 0.095 | 0.070 | +0.52 | -4.5 | -2.4 |
| evergreen | tower NDVI | fW | 89 (12) | +2.1 | -1.1 | 0.013 | 0.070 | -0.00 | +1.7 | -4.7** |
| evergreen | satellite NDVI | GPP | 82 (15) | -1.3 | -0.9 | 0.016 | - | - | - | - |
| evergreen | satellite NDVI | NPPd | 82 (15) | +1.0 | -4.7* | 0.085 | 0.016 | +0.79 | +1.0 | -2.5 |
| evergreen | satellite NDVI | Ra | 82 (15) | -2.6 | +3.8* | 0.060 | 0.016 | +0.37 | -0.7 | -1.5 |
| evergreen | satellite NDVI | Rg | 82 (15) | +0.5 | -4.9** | 0.095 | 0.016 | +0.81 | -0.3 | -1.6 |
| evergreen | satellite NDVI | Rm | 82 (15) | -2.1 | +4.2* | 0.071 | 0.016 | +0.14 | -0.5 | -1.7 |
| evergreen | satellite NDVI | Rm_model | 82 (15) | -4.2 | +4.6 | 0.042 | 0.016 | +0.41 | -0.1 | -1.7 |
| evergreen | satellite NDVI | SI | 89 (17) | -1.3 | +1.2 | 0.009 | 0.016 | +0.32 | -0.7 | -1.5 |
| evergreen | satellite NDVI | fT | 89 (17) | -1.0 | +1.1 | 0.011 | 0.016 | +0.16 | -0.9 | -1.6 |
| evergreen | satellite NDVI | fW | 89 (17) | -0.4 | +0.3 | 0.001 | 0.016 | +0.23 | -0.1 | -1.7 |
| evergreen | satellite NIRv | GPP | 211 (27) | -2.2** | +0.3 | 0.034 | - | - | - | - |
| evergreen | satellite NIRv | NPPd | 211 (27) | -0.8 | -0.6 | 0.010 | 0.034 | +0.67 | +0.7 | -2.5* |
| evergreen | satellite NIRv | Ra | 211 (27) | -1.4 | +0.8 | 0.015 | 0.034 | +0.31 | -0.5 | -1.9** |
| evergreen | satellite NIRv | Rg | 211 (27) | -0.5 | -0.3 | 0.004 | 0.034 | +0.64 | +1.2 | -2.8** |
| evergreen | satellite NIRv | Rm | 211 (27) | -1.1 | +0.7 | 0.010 | 0.034 | +0.15 | -0.6 | -1.9*** |
| evergreen | satellite NIRv | Rm_model | 211 (27) | -1.6 | +1.6 | 0.009 | 0.034 | +0.35 | +0.5 | -2.2*** |
| evergreen | satellite NIRv | SI | 217 (27) | -0.9 | +0.6 | 0.007 | 0.034 | +0.28 | -0.2 | -2.0** |
| evergreen | satellite NIRv | fT | 217 (27) | -1.1 | -0.9 | 0.019 | 0.034 | +0.11 | -1.1 | -1.9** |
| evergreen | satellite NIRv | fW | 217 (27) | -0.2 | +0.7 | 0.003 | 0.034 | +0.19 | +0.5 | -2.1** |

### D. Spring temperature against GPP (60 days before the solstice)

'all weather' = temperature, radiation and water balance of the same window.

| sites | EOS source | n (sites) | r(GPP, T) | GPP alone | T alone | GPP, T fixed | T, GPP fixed | GPP, all weather fixed | T, all else fixed | R2 only GPP | R2 only T |
|---|---|---|---|---|---|---|---|---|---|---|---|
| all | PhenoCam | 233 (30) | +0.14 | -0.7 | -2.4** | -0.4 | -2.4** | -0.4 | -2.5* | 0.001 | 0.031 |
| all | tower NDVI | 357 (41) | +0.35 | -2.5* | -3.0** | -1.6 | -2.5** | -1.3 | -2.2** | 0.009 | 0.022 |
| all | satellite NDVI | 456 (69) | +0.20 | -1.1* | -1.5* | -0.8 | -1.4* | -0.8 | -1.6* | 0.005 | 0.013 |
| all | satellite NIRv | 681 (88) | +0.13 | -1.5*** | -2.1*** | -1.2** | -1.9** | -1.2** | -2.2*** | 0.010 | 0.025 |
| all | all sources stacked | 1727 (94) | +0.19 | -1.5*** | -2.2*** | -1.1** | -2.0*** | -1.1* | -2.1*** | 0.007 | 0.022 |
| summer-green | PhenoCam | 202 (25) | +0.18 | -1.2 | -2.6** | -0.7 | -2.5* | -0.8 | -2.7* | 0.003 | 0.032 |
| summer-green | tower NDVI | 342 (38) | +0.36 | -2.8* | -3.0** | -2.0 | -2.3* | -1.6 | -2.1** | 0.014 | 0.019 |
| summer-green | satellite NDVI | 436 (65) | +0.22 | -1.1* | -1.6* | -0.8 | -1.4* | -0.7 | -1.7* | 0.004 | 0.013 |
| summer-green | satellite NIRv | 611 (79) | +0.17 | -1.6*** | -2.3*** | -1.2** | -2.1** | -1.2** | -2.4*** | 0.010 | 0.027 |
| summer-green | all sources stacked | 1591 (85) | +0.22 | -1.6*** | -2.3*** | -1.2** | -2.0*** | -1.2** | -2.1*** | 0.008 | 0.022 |
| deciduous | PhenoCam | 104 (9) | +0.36 | -0.2 | +0.7 | -0.5 | +0.9 | -0.6 | +0.4 | 0.001 | 0.004 |
| deciduous | tower NDVI | 174 (17) | +0.38 | -0.4 | -2.4*** | +0.6 | -2.7* | +1.0 | -2.3* | 0.002 | 0.040 |
| deciduous | satellite NDVI | 172 (23) | +0.28 | -0.6 | -0.2 | -0.6 | -0.0 | -0.4 | -1.0 | 0.004 | 0.000 |
| deciduous | satellite NIRv | 197 (25) | +0.31 | -1.0 | -2.4** | -0.3 | -2.3** | -0.4 | -2.4** | 0.001 | 0.042 |
| deciduous | all sources stacked | 647 (29) | +0.32 | -0.6 | -1.3** | -0.2 | -1.3* | -0.1 | -1.7** | 0.000 | 0.011 |
| evergreen | tower NDVI | 87 (12) | +0.53 | -4.6** | -4.7** | -2.8 | -3.2 | -1.9 | -3.1* | 0.018 | 0.023 |
| evergreen | satellite NDVI | 88 (17) | +0.18 | -1.4 | -2.0 | -1.0 | -1.8 | -1.1 | -1.3 | 0.005 | 0.015 |
| evergreen | satellite NIRv | 213 (27) | +0.11 | -2.0*** | -1.7* | -1.8** | -1.5 | -1.9** | -1.2 | 0.027 | 0.018 |
| evergreen | all sources stacked | 399 (28) | +0.17 | -2.0*** | -2.5** | -1.6** | -2.2** | -1.5** | -1.8* | 0.014 | 0.025 |

**Mean EOS90 anomaly (days) by class of year** - the last two classes are the test

| sites | EOS source | years | reading | n (sites) | share of years | mean EOS90 anomaly [95% CI] |
|---|---|---|---|---|---|---|
| all | PhenoCam | warm + high GPP | both say earlier | 62 (24) | 27% | +0.3 [-4.0, +4.5] |
| all | PhenoCam | cool + low GPP | both say later | 69 (25) | 30% | +0.3 [-3.7, +4.3] |
| all | PhenoCam | warm + low GPP | earlier if temperature drives it | 54 (26) | 23% | -1.1 [-3.1, +0.8] |
| all | PhenoCam | cool + high GPP | earlier if GPP drives it | 48 (23) | 21% | +0.5 [-3.5, +4.5] |
| all | tower NDVI | warm + high GPP | both say earlier | 112 (39) | 31% | -3.9** [-6.4, -1.5] |
| all | tower NDVI | cool + low GPP | both say later | 110 (39) | 31% | +2.3 [-0.2, +4.8] |
| all | tower NDVI | warm + low GPP | earlier if temperature drives it | 71 (36) | 20% | +1.8 [-1.0, +4.5] |
| all | tower NDVI | cool + high GPP | earlier if GPP drives it | 64 (33) | 18% | +0.9 [-2.7, +4.6] |
| all | satellite NDVI | warm + high GPP | both say earlier | 119 (63) | 26% | -1.6 [-3.9, +0.7] |
| all | satellite NDVI | cool + low GPP | both say later | 129 (64) | 28% | +0.9 [-0.9, +2.7] |
| all | satellite NDVI | warm + low GPP | earlier if temperature drives it | 100 (53) | 22% | -0.8 [-2.9, +1.2] |
| all | satellite NDVI | cool + high GPP | earlier if GPP drives it | 108 (60) | 24% | +1.6 [-0.2, +3.3] |
| all | satellite NIRv | warm + high GPP | both say earlier | 183 (82) | 27% | -2.9*** [-4.6, -1.2] |
| all | satellite NIRv | cool + low GPP | both say later | 192 (84) | 28% | +2.2** [+0.9, +3.6] |
| all | satellite NIRv | warm + low GPP | earlier if temperature drives it | 152 (72) | 22% | -0.4 [-2.1, +1.4] |
| all | satellite NIRv | cool + high GPP | earlier if GPP drives it | 154 (81) | 23% | +1.0 [-0.8, +2.8] |
| all | all sources stacked | warm + high GPP | both say earlier | 476 (89) | 28% | -2.4** [-3.9, -1.0] |
| all | all sources stacked | cool + low GPP | both say later | 500 (91) | 29% | +1.6* [+0.4, +2.9] |
| all | all sources stacked | warm + low GPP | earlier if temperature drives it | 377 (85) | 22% | -0.2 [-1.5, +1.1] |
| all | all sources stacked | cool + high GPP | earlier if GPP drives it | 374 (89) | 22% | +1.1 [-0.4, +2.6] |
| summer-green | PhenoCam | warm + high GPP | both say earlier | 54 (21) | 27% | +0.0 [-4.8, +4.8] |
| summer-green | PhenoCam | cool + low GPP | both say later | 63 (22) | 31% | +1.1 [-3.1, +5.2] |
| summer-green | PhenoCam | warm + low GPP | earlier if temperature drives it | 46 (22) | 23% | -1.0 [-3.2, +1.2] |
| summer-green | PhenoCam | cool + high GPP | earlier if GPP drives it | 39 (19) | 19% | -0.6 [-5.3, +4.1] |
| summer-green | tower NDVI | warm + high GPP | both say earlier | 108 (36) | 32% | -4.0** [-6.6, -1.5] |
| summer-green | tower NDVI | cool + low GPP | both say later | 108 (37) | 32% | +2.6* [+0.1, +5.1] |
| summer-green | tower NDVI | warm + low GPP | earlier if temperature drives it | 67 (33) | 20% | +2.1 [-0.7, +5.0] |
| summer-green | tower NDVI | cool + high GPP | earlier if GPP drives it | 59 (31) | 17% | +0.2 [-3.6, +4.1] |
| summer-green | satellite NDVI | warm + high GPP | both say earlier | 116 (60) | 27% | -1.5 [-3.9, +0.8] |
| summer-green | satellite NDVI | cool + low GPP | both say later | 126 (61) | 29% | +1.0 [-0.8, +2.8] |
| summer-green | satellite NDVI | warm + low GPP | earlier if temperature drives it | 95 (50) | 22% | -0.9 [-3.1, +1.2] |
| summer-green | satellite NDVI | cool + high GPP | earlier if GPP drives it | 99 (56) | 23% | +1.4 [-0.4, +3.2] |
| summer-green | satellite NIRv | warm + high GPP | both say earlier | 171 (75) | 28% | -3.0*** [-4.7, -1.3] |
| summer-green | satellite NIRv | cool + low GPP | both say later | 177 (76) | 29% | +2.5*** [+1.1, +3.9] |
| summer-green | satellite NIRv | warm + low GPP | earlier if temperature drives it | 132 (64) | 22% | -0.5 [-2.4, +1.5] |
| summer-green | satellite NIRv | cool + high GPP | earlier if GPP drives it | 131 (72) | 21% | +1.0 [-1.0, +3.0] |
| summer-green | all sources stacked | warm + high GPP | both say earlier | 449 (81) | 28% | -2.5** [-4.0, -1.0] |
| summer-green | all sources stacked | cool + low GPP | both say later | 474 (83) | 30% | +1.9** [+0.6, +3.2] |
| summer-green | all sources stacked | warm + low GPP | earlier if temperature drives it | 340 (77) | 21% | -0.1 [-1.6, +1.3] |
| summer-green | all sources stacked | cool + high GPP | earlier if GPP drives it | 328 (80) | 21% | +0.8 [-0.8, +2.4] |
| deciduous | PhenoCam | warm + high GPP | both say earlier | 30 (9) | 29% | +4.1 [-3.2, +11.4] |
| deciduous | PhenoCam | cool + low GPP | both say later | 37 (9) | 36% | -1.4 [-7.2, +4.4] |
| deciduous | PhenoCam | warm + low GPP | earlier if temperature drives it | 17 (8) | 16% | -1.0 [-5.4, +3.4] |
| deciduous | PhenoCam | cool + high GPP | earlier if GPP drives it | 20 (7) | 19% | -2.7 [-9.8, +4.3] |
| deciduous | tower NDVI | warm + high GPP | both say earlier | 59 (17) | 34% | -1.9* [-3.7, -0.0] |
| deciduous | tower NDVI | cool + low GPP | both say later | 55 (17) | 32% | +1.0 [-1.2, +3.1] |
| deciduous | tower NDVI | warm + low GPP | earlier if temperature drives it | 29 (15) | 17% | +0.9 [-4.0, +5.7] |
| deciduous | tower NDVI | cool + high GPP | earlier if GPP drives it | 31 (15) | 18% | +1.0 [-4.3, +6.3] |
| deciduous | satellite NDVI | warm + high GPP | both say earlier | 46 (21) | 27% | -0.3 [-2.8, +2.2] |
| deciduous | satellite NDVI | cool + low GPP | both say later | 55 (23) | 32% | -0.2 [-2.6, +2.1] |
| deciduous | satellite NDVI | warm + low GPP | earlier if temperature drives it | 33 (18) | 19% | -0.4 [-3.4, +2.7] |
| deciduous | satellite NDVI | cool + high GPP | earlier if GPP drives it | 38 (21) | 22% | +1.0 [-1.3, +3.4] |
| deciduous | satellite NIRv | warm + high GPP | both say earlier | 58 (23) | 29% | -1.7 [-4.1, +0.8] |
| deciduous | satellite NIRv | cool + low GPP | both say later | 63 (25) | 32% | +3.0* [+0.7, +5.3] |
| deciduous | satellite NIRv | warm + low GPP | earlier if temperature drives it | 38 (21) | 19% | -2.3 [-5.1, +0.4] |
| deciduous | satellite NIRv | cool + high GPP | earlier if GPP drives it | 38 (22) | 19% | -0.1 [-3.1, +3.0] |
| deciduous | all sources stacked | warm + high GPP | both say earlier | 193 (27) | 30% | -0.5 [-2.3, +1.3] |
| deciduous | all sources stacked | cool + low GPP | both say later | 210 (29) | 32% | +0.8 [-0.9, +2.5] |
| deciduous | all sources stacked | warm + low GPP | earlier if temperature drives it | 117 (27) | 18% | -0.8 [-3.1, +1.5] |
| deciduous | all sources stacked | cool + high GPP | earlier if GPP drives it | 127 (27) | 20% | +0.1 [-1.9, +2.2] |
| evergreen | tower NDVI | warm + high GPP | both say earlier | 28 (12) | 32% | -9.7*** [-13.3, -6.0] |
| evergreen | tower NDVI | cool + low GPP | both say later | 33 (12) | 38% | +6.5** [+2.6, +10.5] |
| evergreen | tower NDVI | warm + low GPP | earlier if temperature drives it | 14 (9) | 16% | +5.1 [-0.4, +10.6] |
| evergreen | tower NDVI | cool + high GPP | earlier if GPP drives it | 12 (7) | 14% | -1.4 [-10.6, +7.8] |
| evergreen | satellite NDVI | warm + high GPP | both say earlier | 27 (17) | 31% | -0.8 [-7.9, +6.2] |
| evergreen | satellite NDVI | cool + low GPP | both say later | 26 (15) | 30% | +0.1 [-4.8, +4.9] |
| evergreen | satellite NDVI | warm + low GPP | earlier if temperature drives it | 14 (9) | 16% | -1.8 [-12.0, +8.4] |
| evergreen | satellite NDVI | cool + high GPP | earlier if GPP drives it | 21 (14) | 24% | +2.2 [-1.5, +5.8] |
| evergreen | satellite NIRv | warm + high GPP | both say earlier | 61 (26) | 29% | -1.9 [-5.0, +1.1] |
| evergreen | satellite NIRv | cool + low GPP | both say later | 63 (27) | 30% | +2.1 [-0.3, +4.5] |
| evergreen | satellite NIRv | warm + low GPP | earlier if temperature drives it | 40 (20) | 19% | -0.6 [-4.7, +3.4] |
| evergreen | satellite NIRv | cool + high GPP | earlier if GPP drives it | 49 (25) | 23% | +0.3 [-2.5, +3.0] |
| evergreen | all sources stacked | warm + high GPP | both say earlier | 118 (27) | 30% | -3.4 [-6.9, +0.1] |
| evergreen | all sources stacked | cool + low GPP | both say later | 124 (27) | 31% | +2.8* [+0.6, +5.0] |
| evergreen | all sources stacked | warm + low GPP | earlier if temperature drives it | 72 (23) | 18% | +0.2 [-3.8, +4.2] |
| evergreen | all sources stacked | cool + high GPP | earlier if GPP drives it | 85 (27) | 21% | +0.5 [-1.9, +2.9] |

### E. Leaf age: does senescence follow leaf-out day for day? (deciduous forests)

Within-site slope of the EOS date on the leaf-out date (SOS50), in days per day. A fixed leaf life span predicts a slope of 1; no link predicts 0. 'PhenoCam leaf-out' rows take the leaf-out date from the camera and the EOS from another instrument, so errors of one curve fit cannot produce the relation.

| EOS from | leaf-out from | EOS level | model | n (sites) | slope [95% CI] | differs from 0 | differs from 1 | SD leaf-out / EOS / life span (days) |
|---|---|---|---|---|---|---|---|---|
| PhenoCam | PhenoCam | EOS90 | leaf-out alone | 116 (10) | +0.15 [-0.21, +0.51] | p = 0.413 | p = 0.000 | 6.5 / 13.0 / 14.1 |
| PhenoCam | PhenoCam | EOS90 | + spring temperature and GPP | 104 (9) | +0.20 [-0.18, +0.57] | p = 0.304 | p = 0.000 | 6.4 / 13.5 / 14.5 |
| PhenoCam | PhenoCam | EOS50 | leaf-out alone | 116 (10) | +0.05 [-0.43, +0.53] | p = 0.833 | p = 0.000 | 6.5 / 10.6 / 12.3 |
| PhenoCam | PhenoCam | EOS50 | + spring temperature and GPP | 104 (9) | -0.00 [-0.56, +0.56] | p = 0.992 | p = 0.000 | 6.4 / 11.1 / 12.8 |
| PhenoCam | PhenoCam | EOS10 | leaf-out alone | 116 (10) | -0.03 [-0.80, +0.74] | p = 0.941 | p = 0.009 | 6.5 / 13.5 / 15.1 |
| PhenoCam | PhenoCam | EOS10 | + spring temperature and GPP | 104 (9) | -0.20 [-1.26, +0.85] | p = 0.707 | p = 0.026 | 6.4 / 14.2 / 15.8 |
| tower NDVI | tower NDVI | EOS90 | leaf-out alone | 189 (18) | -0.05 [-0.35, +0.25] | p = 0.750 | p = 0.000 | 7.7 / 12.3 / 14.7 |
| tower NDVI | tower NDVI | EOS90 | + spring temperature and GPP | 174 (17) | +0.01 [-0.23, +0.26] | p = 0.906 | p = 0.000 | 7.4 / 12.3 / 14.0 |
| tower NDVI | tower NDVI | EOS50 | leaf-out alone | 189 (18) | -0.15 [-0.38, +0.08] | p = 0.192 | p = 0.000 | 7.7 / 9.6 / 13.0 |
| tower NDVI | tower NDVI | EOS50 | + spring temperature and GPP | 174 (17) | -0.15 [-0.39, +0.09] | p = 0.232 | p = 0.000 | 7.4 / 9.6 / 12.4 |
| tower NDVI | tower NDVI | EOS10 | leaf-out alone | 189 (18) | -0.29 [-0.64, +0.05] | p = 0.093 | p = 0.000 | 7.7 / 13.0 / 16.2 |
| tower NDVI | tower NDVI | EOS10 | + spring temperature and GPP | 174 (17) | -0.35 [-0.68, -0.01] | p = 0.044 | p = 0.000 | 7.4 / 12.8 / 15.6 |
| tower NDVI | PhenoCam | EOS90 | leaf-out alone | 37 (5) | +0.24 [-0.52, +1.01] | p = 0.538 | p = 0.052 | 7.2 / 10.2 / 11.5 |
| tower NDVI | PhenoCam | EOS50 | leaf-out alone | 37 (5) | +0.07 [-0.39, +0.53] | p = 0.774 | p = 0.000 | 7.2 / 7.3 / 9.9 |
| tower NDVI | PhenoCam | EOS10 | leaf-out alone | 37 (5) | -0.14 [-0.39, +0.11] | p = 0.279 | p = 0.000 | 7.2 / 10.5 / 13.3 |
| satellite NDVI | satellite NDVI | EOS90 | leaf-out alone | 194 (26) | +0.26 [-0.04, +0.56] | p = 0.093 | p = 0.000 | 6.2 / 9.1 / 10.0 |
| satellite NDVI | satellite NDVI | EOS90 | + spring temperature and GPP | 172 (23) | +0.26 [-0.14, +0.65] | p = 0.201 | p = 0.000 | 6.2 / 9.3 / 10.3 |
| satellite NDVI | satellite NDVI | EOS50 | leaf-out alone | 194 (26) | +0.05 [-0.10, +0.19] | p = 0.522 | p = 0.000 | 6.2 / 5.3 / 7.9 |
| satellite NDVI | satellite NDVI | EOS50 | + spring temperature and GPP | 172 (23) | +0.03 [-0.14, +0.20] | p = 0.753 | p = 0.000 | 6.2 / 5.5 / 8.1 |
| satellite NDVI | satellite NDVI | EOS10 | leaf-out alone | 194 (26) | -0.14 [-0.39, +0.11] | p = 0.267 | p = 0.000 | 6.2 / 9.7 / 12.0 |
| satellite NDVI | satellite NDVI | EOS10 | + spring temperature and GPP | 172 (23) | -0.17 [-0.50, +0.17] | p = 0.325 | p = 0.000 | 6.2 / 10.1 / 12.3 |
| satellite NDVI | PhenoCam | EOS90 | leaf-out alone | 56 (8) | +0.22 [-0.08, +0.52] | p = 0.159 | p = 0.000 | 5.2 / 6.5 / 7.6 |
| satellite NDVI | PhenoCam | EOS50 | leaf-out alone | 56 (8) | +0.03 [-0.26, +0.32] | p = 0.843 | p = 0.000 | 5.2 / 4.3 / 6.6 |
| satellite NDVI | PhenoCam | EOS10 | leaf-out alone | 56 (8) | -0.17 [-0.65, +0.32] | p = 0.499 | p = 0.000 | 5.2 / 7.8 / 9.8 |
| satellite NIRv | satellite NIRv | EOS90 | leaf-out alone | 224 (28) | -0.12 [-0.47, +0.23] | p = 0.494 | p = 0.000 | 5.8 / 10.7 / 12.5 |
| satellite NIRv | satellite NIRv | EOS90 | + spring temperature and GPP | 197 (25) | -0.29 [-0.66, +0.07] | p = 0.119 | p = 0.000 | 5.7 / 10.7 / 12.3 |
| satellite NIRv | satellite NIRv | EOS50 | leaf-out alone | 224 (28) | -0.18 [-0.42, +0.06] | p = 0.141 | p = 0.000 | 5.8 / 6.9 / 9.6 |
| satellite NIRv | satellite NIRv | EOS50 | + spring temperature and GPP | 197 (25) | -0.26 [-0.48, -0.03] | p = 0.025 | p = 0.000 | 5.7 / 6.8 / 9.5 |
| satellite NIRv | satellite NIRv | EOS10 | leaf-out alone | 224 (28) | -0.18 [-0.48, +0.12] | p = 0.240 | p = 0.000 | 5.8 / 10.9 / 12.8 |
| satellite NIRv | satellite NIRv | EOS10 | + spring temperature and GPP | 197 (25) | -0.19 [-0.50, +0.11] | p = 0.215 | p = 0.000 | 5.7 / 11.0 / 12.9 |
| satellite NIRv | PhenoCam | EOS90 | leaf-out alone | 58 (8) | +0.00 [-0.66, +0.66] | p = 0.992 | p = 0.003 | 5.1 / 12.0 / 13.1 |
| satellite NIRv | PhenoCam | EOS50 | leaf-out alone | 58 (8) | +0.14 [-0.21, +0.50] | p = 0.423 | p = 0.000 | 5.1 / 5.5 / 7.0 |
| satellite NIRv | PhenoCam | EOS10 | leaf-out alone | 58 (8) | +0.34 [+0.05, +0.62] | p = 0.020 | p = 0.000 | 5.1 / 7.6 / 8.2 |

For comparison, EOS90 on SOS50 in the other plant types (same source for both dates): evergreen, tower NDVI +0.01; evergreen, satellite NDVI +0.20; evergreen, satellite NIRv +0.22*; grass/shrub, PhenoCam +0.32; grass/shrub, tower NDVI +0.19; grass/shrub, satellite NDVI +0.23*; grass/shrub, satellite NIRv +0.11.

**What this shows**

- Each explanation is tested by putting the competing variable in the same within-site model as GPP. 'Stacked' rows use all four EOS sources at once (each site-year once per source); they have the most power and are quoted here.
- **A. Leaf-out date.** In deciduous forests (647 site-years) leaf-out date predicts EOS90 (+1.4** days per SD with GPP in the model: earlier leaf-out, earlier senescence) and early GPP does not (-0.3). Over all vegetation types it is the reverse: GPP -1.5***, leaf-out -0.4. In deciduous trees the timing of the whole leaf cycle shifts together (a leaf that emerges early ages early); the GPP effect is not a deciduous-forest effect.
- **B. Water.** At dry-summer sites (136 site-years, 9 sites) the spring water balance predicts EOS90 (+1.8**: a wet spring delays senescence) and GPP does not (+0.3). At summer-green sites (1591 site-years) water has no effect (+0.2) and the GPP effect is unchanged by it (-1.6***). Measured plant-available soil water gives the same answer at summer-green sites (water +0.3, GPP -1.8***; 1274 site-years). Water limitation explains senescence where summers are dry; it does not explain the GPP effect elsewhere.
- **C. Sink variables.** NPP and the respiration terms (correlation with GPP within sites 0.17 to 0.78) give pre-solstice effects of -2.0 to +0.6 days per SD; next to GPP in the same model they are significant in 0 of 16 cases. They add nothing to GPP, of which they are a rescaled copy. Maintenance respiration from the fitted model, which does not contain the day's GPP: PhenoCam (GCC) +0.2, Tower NDVI -3.1, Satellite NDVI -1.6, Satellite NIRv -0.7; next to GPP PhenoCam (GCC) +0.2, Tower NDVI -0.7, Satellite NDVI +0.3, Satellite NIRv +0.4. The temperature part of the sink-limitation index does predict EOS90 (PhenoCam (GCC) -1.7*, Tower NDVI -2.9**, Satellite NDVI -1.3*, Satellite NIRv -2.0***) and keeps an effect next to GPP in 3 of 4 sources.
- **D. Spring temperature.** GPP and temperature of the 60 days before the solstice are only weakly correlated within sites (r = 0.19), so their effects can be separated. Together in one model: temperature -2.0***, GPP -1.1** days per SD. With radiation and water balance also held fixed, GPP keeps -1.1*. Temperature alone explains 4 times more variance than GPP adds to it.
- By plant type, GPP net of all weather: evergreen -1.5**, deciduous -0.1. The GPP effect that the weather does not explain is carried by evergreen forests.
- **Years in which temperature and GPP diverge.** Mean EOS90 anomaly: warm + high GPP -2.4** days; cool + low GPP +1.6* days; warm + low GPP -0.2 days; cool + high GPP +1.1 days. Senescence is early when spring is both warm and productive and late when it is both cool and unproductive. When only one of the two is high the shift is small. Neither variable alone drives the onset of senescence; the two add up.
- **E. Leaf age.** If leaves had a fixed life span, a leaf-out one day earlier would bring senescence one day earlier: a slope of 1 day per day. In deciduous forests the within-site slope of EOS90 on leaf-out (SOS50) is PhenoCam (GCC) +0.15, Tower NDVI -0.05, Satellite NDVI +0.26, Satellite NIRv -0.12 days per day; with the PhenoCam, whose leaf-out dates are the most precise, +0.15 (95% CI -0.21 to +0.51, 116 site-years). Every one of these slopes is significantly below 1 (4 of 4), and 0 of 4 differ from 0.
- Taking leaf-out from the PhenoCam and EOS90 from another instrument (so that errors of one curve fit cannot link the two dates) gives Tower NDVI +0.24, Satellite NDVI +0.22, Satellite NIRv +0.00 days per day.
- The leaf life span itself (EOS90 minus leaf-out) varies by 10.0 to 14.7 days between years at a site, more than EOS90 does (9.1 to 13.0 days). A fixed life span would make it vary less. **The leaf-age explanation is rejected:** senescence follows leaf-out by at most a fifth to a quarter of a day per day, and that weak link is not distinguishable from zero in most sources. Leaf-out date matters for the onset of senescence in deciduous forests a little, but not through a fixed leaf life span.
- **Overall.** About half of the uncorrected 'pre-solstice GPP' effect is a spring-temperature effect. What remains (about one day per SD) is not explained by leaf-out date, water or weather, and is found mainly in evergreen forests. Whether it is a sink effect cannot be told from flux-derived variables, because every one of them is tied to GPP.

![Deciduous forests: leaf-out date against early GPP (step 44)](../figure/mechanism/leafout_vs_gpp.png)

*Deciduous forests: leaf-out date against early GPP (step 44)*

![Pre-solstice GPP against the spring water balance, by season type (step 44)](../figure/mechanism/water.png)

*Pre-solstice GPP against the spring water balance, by season type (step 44)*

![GPP and sink variables in the same two-window model (step 44)](../figure/mechanism/sink.png)

*GPP and sink variables in the same two-window model (step 44)*

![Spring temperature against pre-solstice GPP, and the years in which they diverge (step 44)](../figure/mechanism/temperature.png)

*Spring temperature against pre-solstice GPP, and the years in which they diverge (step 44)*

![Deciduous forests: onset of senescence against leaf-out date, with the 1:1 line of a fixed leaf life span (step 44)](../figure/mechanism/leaf_age.png)

*Deciduous forests: onset of senescence against leaf-out date, with the 1:1 line of a fixed leaf life span (step 44)*

## 15. Measured stem growth: the sink itself (steps 29, 45)

Dendrometer data at the few flux sites that have them (`Output/dendrometer_datasets.md`). Growth is counted by the zero-growth concept (only new maxima of stem size), scaled so that 1 is a normal year's growth of a tree, and averaged over trees. Tables: `data/dendro_growth_by_site_year.csv`, `data/growth_vs_gpp.csv`, `data/growth_vs_senescence.csv`.


### Growth data

| site | dataset | years | units (median) | years with flux data | years with timing | share of growth done by the solstice |
|---|---|---|---|---|---|---|
| AT-Zoe | LTER Zoebelboden | 24 (1996-2019) | 4 | 4 | 22 | 40% |
| CH-Dav | DenDrought2018 | 3 (2016-2018) | 9 | 3 | 3 | 34% |
| CZ-BK1 | DenDrought2018 | 3 (2016-2018) | 8 | 3 | 3 | 44% |
| CZ-RAJ | DenDrought2018 | 2 (2017-2018) | 7 | 2 | 2 | 36% |
| CZ-Stn | DenDrought2018 | 1 (2018-2018) | 7 | 1 | 1 | 78% |
| FR-Fon | DenDrought2018 | 3 (2016-2018) | 12 | 0 | 3 | 69% |
| US-Ha1 | Harvard Forest HF069 | 24 (1999-2024) | 620 | 24 | 2 | 39% |
| US-Ton | Rao et al. (Dryad) | 1 (2023-2023) | 5 | 1 | 1 | 93% |

### A. Is growth coupled to GPP?

| relation | within-site r | p | site-years | sites (years) |
|---|---|---|---|---|
| annual growth ~ annual GPP | +0.42 | 0.003 | 34 | US-Ha1 (24), AT-Zoe (4), CH-Dav (3), CZ-BK1 (3) |
| growth by the solstice ~ pre-solstice GPP | +0.07 | 0.731 | 34 | US-Ha1 (24), AT-Zoe (4), CH-Dav (3), CZ-BK1 (3) |
| growth rate before the solstice ~ pre-solstice GPP | +0.05 | 0.794 | 34 | US-Ha1 (24), AT-Zoe (4), CH-Dav (3), CZ-BK1 (3) |
| growth by the solstice ~ spring temperature | +0.41 | 0.019 | 36 | US-Ha1 (24), AT-Zoe (6), CH-Dav (3), CZ-BK1 (3) |

### B. Does growth predict EOS90? (days per +1 within-site SD)

| EOS source | predictor | slope [95% CI] | p | site-years | sites (years) |
|---|---|---|---|---|---|
| PhenoCam | growth by the solstice | -4.9 [-12.0, +2.1] | 0.171 | 15 | US-Ha1 (15) |
| PhenoCam | growth rate before the solstice | -3.3 [-10.6, +4.0] | 0.373 | 15 | US-Ha1 (15) |
| PhenoCam | annual growth | -6.5 [-12.5, -0.5] | 0.034 | 15 | US-Ha1 (15) |
| PhenoCam | share of growth done by the solstice | +0.7 [-7.2, +8.5] | 0.869 | 15 | US-Ha1 (15) |
| PhenoCam | GPP, 60 days before the solstice | -6.0 [-9.6, -2.4] | 0.001 | 15 | US-Ha1 (15) |
| PhenoCam | air temperature, 60 days before the solstice | +2.2 [-4.5, +8.9] | 0.518 | 15 | US-Ha1 (15) |
| satellite NIRv | growth by the solstice | +4.5 [-2.2, +11.2] | 0.190 | 16 | US-Ha1 (9), AT-Zoe (4), CH-Dav (3) |
| satellite NIRv | growth rate before the solstice | +4.5 [-2.2, +11.2] | 0.185 | 16 | US-Ha1 (9), AT-Zoe (4), CH-Dav (3) |
| satellite NIRv | annual growth | +3.5 [-4.0, +11.1] | 0.356 | 16 | US-Ha1 (9), AT-Zoe (4), CH-Dav (3) |
| satellite NIRv | share of growth done by the solstice | +3.7 [-1.3, +8.7] | 0.150 | 16 | US-Ha1 (9), AT-Zoe (4), CH-Dav (3) |
| satellite NIRv | date growth stops | +3.6 [-3.3, +10.5] | 0.308 | 11 | AT-Zoe (4), CH-Dav (3), CZ-BK1 (2), CZ-RAJ (2) |
| satellite NIRv | date half of the growth is done | +4.3 [-1.9, +10.4] | 0.172 | 11 | AT-Zoe (4), CH-Dav (3), CZ-BK1 (2), CZ-RAJ (2) |
| satellite NIRv | date of fastest growth | +6.6 [+2.3, +10.9] | 0.003 | 11 | AT-Zoe (4), CH-Dav (3), CZ-BK1 (2), CZ-RAJ (2) |
| satellite NIRv | fastest growth rate | +9.2 [+6.7, +11.7] | 0.000 | 11 | AT-Zoe (4), CH-Dav (3), CZ-BK1 (2), CZ-RAJ (2) |
| satellite NIRv | GPP, 60 days before the solstice | +1.7 [-6.5, +9.9] | 0.683 | 12 | US-Ha1 (9), CH-Dav (3) |
| satellite NIRv | air temperature, 60 days before the solstice | -0.4 [-5.2, +4.4] | 0.866 | 16 | US-Ha1 (9), AT-Zoe (4), CH-Dav (3) |
| all sources stacked | growth by the solstice | +0.8 [-3.2, +4.9] | 0.690 | 40 | US-Ha1 (33), AT-Zoe (4), CH-Dav (3) |
| all sources stacked | growth rate before the solstice | +1.4 [-2.6, +5.4] | 0.499 | 40 | US-Ha1 (33), AT-Zoe (4), CH-Dav (3) |
| all sources stacked | annual growth | -2.1 [-6.8, +2.6] | 0.386 | 40 | US-Ha1 (33), AT-Zoe (4), CH-Dav (3) |
| all sources stacked | share of growth done by the solstice | +2.5 [-0.8, +5.9] | 0.141 | 40 | US-Ha1 (33), AT-Zoe (4), CH-Dav (3) |
| all sources stacked | date growth stops | +3.2 [-3.0, +9.4] | 0.310 | 13 | AT-Zoe (6), CH-Dav (3), CZ-BK1 (2), CZ-RAJ (2) |
| all sources stacked | date half of the growth is done | +3.9 [-1.6, +9.4] | 0.167 | 13 | AT-Zoe (6), CH-Dav (3), CZ-BK1 (2), CZ-RAJ (2) |
| all sources stacked | date of fastest growth | +6.0 [+2.1, +9.9] | 0.003 | 13 | AT-Zoe (6), CH-Dav (3), CZ-BK1 (2), CZ-RAJ (2) |
| all sources stacked | fastest growth rate | +8.4 [+6.2, +10.7] | 0.000 | 13 | AT-Zoe (6), CH-Dav (3), CZ-BK1 (2), CZ-RAJ (2) |
| all sources stacked | GPP, 60 days before the solstice | -2.4 [-5.8, +1.0] | 0.166 | 36 | US-Ha1 (33), CH-Dav (3) |
| all sources stacked | air temperature, 60 days before the solstice | +0.7 [-2.3, +3.6] | 0.645 | 40 | US-Ha1 (33), AT-Zoe (4), CH-Dav (3) |

### C. Carbon budget: total sink, NPP and what is left (DOY 90-305, gC m-2)

- respiration = autotrophic respiration = (1 - CUE) x GPP, with the daily CUE of step 26
- growth = measured biomass growth, the dendrometer curve scaled to the stand's mean annual growth: US-Ha1: 170 gC m-2 per year, aboveground woody increment (HF069); AT-Zoe: 255 gC m-2 per year, stem + branch + foliage growth, dry mass x 0.5 (LTER Zoebelboden)
- **total sink = growth + respiration**
- **NPP = GPP - respiration**
- **GPP - total sink = NPP - growth**: carbon fixed but used neither for respiration nor for the measured growth. It goes to what the dendrometers do not see (roots; at US-Ha1 also leaves) and to reserves.

| site | years | window | GPP | autotrophic respiration | growth | total sink (growth + respiration) | NPP (GPP - respiration) | GPP - total sink | growth respiration (model) | maintenance respiration (model) | autotrophic respiration (model) | total sink (growth + model respiration) | GPP - total sink (model respiration) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| AT-Zoe | 5 | whole season | 983 | 315 | 205 | 521 | 667 | 462 | 133 | 192 | 325 | 530 | 452 |
| AT-Zoe | 5 | before the solstice | 347 | 88 | 83 | 181 | 260 | 166 | 47 | 22 | 69 | 163 | 184 |
| AT-Zoe | 5 | after the solstice | 623 | 207 | 122 | 329 | 417 | 295 | 84 | 166 | 250 | 372 | 251 |
| US-Ha1 | 24 | whole season | 1566 | 353 | 165 | 518 | 1213 | 1048 | 226 | 260 | 486 | 651 | 914 |
| US-Ha1 | 24 | before the solstice | 398 | 97 | 64 | 161 | 301 | 237 | 57 | 19 | 75 | 140 | 258 |
| US-Ha1 | 24 | after the solstice | 1168 | 256 | 101 | 357 | 912 | 811 | 169 | 242 | 411 | 512 | 656 |

Model respiration (step 26): growth respiration = gR x CUE x GPP; maintenance respiration = temperature response x biomass built since 1 January. The maintenance term does not contain the day's GPP, so a total sink built with it is not GPP times a factor.

Growth as a share of NPP over the season: AT-Zoe 33%, US-Ha1 14%.

**EOS90 against the budget terms** (days per +1 within-site SD)

| EOS source | window | term | slope [95% CI] | p | site-years | sites (years) |
|---|---|---|---|---|---|---|
| PhenoCam | pre-solstice | GPP | -6.9 [-10.9, -3.0] | 0.001 | 15 | US-Ha1 (15) |
| PhenoCam | pre-solstice | autotrophic respiration | -5.5 [-9.9, -1.1] | 0.014 | 15 | US-Ha1 (15) |
| PhenoCam | pre-solstice | growth | -4.9 [-11.9, +2.2] | 0.180 | 15 | US-Ha1 (15) |
| PhenoCam | pre-solstice | total sink (growth + respiration) | -5.9 [-10.6, -1.3] | 0.012 | 15 | US-Ha1 (15) |
| PhenoCam | pre-solstice | NPP (GPP - respiration) | -3.5 [-10.5, +3.5] | 0.327 | 15 | US-Ha1 (15) |
| PhenoCam | pre-solstice | GPP - total sink | -2.5 [-10.0, +4.9] | 0.503 | 15 | US-Ha1 (15) |
| PhenoCam | pre-solstice | growth respiration (model) | -6.6 [-12.4, -0.9] | 0.024 | 15 | US-Ha1 (15) |
| PhenoCam | pre-solstice | maintenance respiration (model) | -3.6 [-6.8, -0.5] | 0.025 | 15 | US-Ha1 (15) |
| PhenoCam | pre-solstice | autotrophic respiration (model) | -6.5 [-10.8, -2.2] | 0.003 | 15 | US-Ha1 (15) |
| PhenoCam | pre-solstice | total sink (growth + model respiration) | -7.3 [-12.9, -1.8] | 0.010 | 15 | US-Ha1 (15) |
| PhenoCam | pre-solstice | GPP - total sink (model respiration) | -6.3 [-9.6, -3.0] | 0.000 | 15 | US-Ha1 (15) |
| PhenoCam | post-solstice | GPP | -7.7 [-11.2, -4.1] | 0.000 | 15 | US-Ha1 (15) |
| PhenoCam | post-solstice | autotrophic respiration | -2.9 [-8.0, +2.3] | 0.272 | 15 | US-Ha1 (15) |
| PhenoCam | post-solstice | growth | -3.3 [-12.8, +6.1] | 0.489 | 15 | US-Ha1 (15) |
| PhenoCam | post-solstice | total sink (growth + respiration) | -3.2 [-8.1, +1.6] | 0.194 | 15 | US-Ha1 (15) |
| PhenoCam | post-solstice | NPP (GPP - respiration) | -4.2 [-12.9, +4.5] | 0.345 | 15 | US-Ha1 (15) |
| PhenoCam | post-solstice | GPP - total sink | -4.1 [-13.0, +4.9] | 0.373 | 15 | US-Ha1 (15) |
| PhenoCam | post-solstice | growth respiration (model) | -5.6 [-13.9, +2.8] | 0.193 | 15 | US-Ha1 (15) |
| PhenoCam | post-solstice | maintenance respiration (model) | -3.9 [-8.4, +0.5] | 0.083 | 15 | US-Ha1 (15) |
| PhenoCam | post-solstice | autotrophic respiration (model) | -5.4 [-10.9, +0.1] | 0.054 | 15 | US-Ha1 (15) |
| PhenoCam | post-solstice | total sink (growth + model respiration) | -5.6 [-11.3, +0.2] | 0.057 | 15 | US-Ha1 (15) |
| PhenoCam | post-solstice | GPP - total sink (model respiration) | -3.7 [-14.2, +6.7] | 0.485 | 15 | US-Ha1 (15) |
| satellite NIRv | pre-solstice | growth | +10.7 [+5.3, +16.0] | 0.000 | 12 | US-Ha1 (9), AT-Zoe (3) |
| satellite NIRv | post-solstice | GPP | +7.5 [+1.1, +14.0] | 0.022 | 12 | US-Ha1 (9), AT-Zoe (3) |
| satellite NIRv | post-solstice | autotrophic respiration | -1.9 [-9.3, +5.5] | 0.611 | 12 | US-Ha1 (9), AT-Zoe (3) |
| satellite NIRv | post-solstice | growth | -1.5 [-10.9, +7.9] | 0.755 | 12 | US-Ha1 (9), AT-Zoe (3) |
| satellite NIRv | post-solstice | total sink (growth + respiration) | -2.2 [-9.0, +4.7] | 0.540 | 12 | US-Ha1 (9), AT-Zoe (3) |
| satellite NIRv | post-solstice | NPP (GPP - respiration) | +7.3 [+1.6, +13.1] | 0.013 | 12 | US-Ha1 (9), AT-Zoe (3) |
| satellite NIRv | post-solstice | GPP - total sink | +7.5 [+2.0, +13.0] | 0.008 | 12 | US-Ha1 (9), AT-Zoe (3) |
| satellite NIRv | post-solstice | growth respiration (model) | +5.7 [-1.5, +12.8] | 0.119 | 12 | US-Ha1 (9), AT-Zoe (3) |
| satellite NIRv | post-solstice | maintenance respiration (model) | +7.6 [+2.2, +13.1] | 0.006 | 12 | US-Ha1 (9), AT-Zoe (3) |
| satellite NIRv | post-solstice | autotrophic respiration (model) | +8.6 [+3.4, +13.7] | 0.001 | 12 | US-Ha1 (9), AT-Zoe (3) |
| satellite NIRv | post-solstice | total sink (growth + model respiration) | +8.1 [+1.1, +15.0] | 0.023 | 12 | US-Ha1 (9), AT-Zoe (3) |
| satellite NIRv | post-solstice | GPP - total sink (model respiration) | +5.6 [-1.7, +12.9] | 0.132 | 12 | US-Ha1 (9), AT-Zoe (3) |
| all sources stacked | pre-solstice | GPP | -2.8 [-6.1, +0.5] | 0.102 | 33 | US-Ha1 (33) |
| all sources stacked | pre-solstice | autotrophic respiration | -2.6 [-5.4, +0.3] | 0.077 | 33 | US-Ha1 (33) |
| all sources stacked | pre-solstice | growth | +1.8 [-3.1, +6.7] | 0.471 | 36 | US-Ha1 (33), AT-Zoe (3) |
| all sources stacked | pre-solstice | total sink (growth + respiration) | -2.2 [-5.2, +0.8] | 0.157 | 33 | US-Ha1 (33) |
| all sources stacked | pre-solstice | NPP (GPP - respiration) | -1.1 [-5.5, +3.2] | 0.610 | 33 | US-Ha1 (33) |
| all sources stacked | pre-solstice | GPP - total sink | -1.3 [-5.5, +2.8] | 0.532 | 33 | US-Ha1 (33) |
| all sources stacked | pre-solstice | growth respiration (model) | -2.6 [-6.8, +1.5] | 0.218 | 33 | US-Ha1 (33) |
| all sources stacked | pre-solstice | maintenance respiration (model) | -1.4 [-4.3, +1.5] | 0.353 | 33 | US-Ha1 (33) |
| all sources stacked | pre-solstice | autotrophic respiration (model) | -2.6 [-5.9, +0.7] | 0.128 | 33 | US-Ha1 (33) |
| all sources stacked | pre-solstice | total sink (growth + model respiration) | -1.4 [-6.0, +3.3] | 0.565 | 33 | US-Ha1 (33) |
| all sources stacked | pre-solstice | GPP - total sink (model respiration) | -3.0 [-6.0, -0.1] | 0.046 | 33 | US-Ha1 (33) |
| all sources stacked | post-solstice | GPP | -2.1 [-7.1, +2.9] | 0.415 | 36 | US-Ha1 (33), AT-Zoe (3) |
| all sources stacked | post-solstice | autotrophic respiration | -1.8 [-5.0, +1.3] | 0.258 | 36 | US-Ha1 (33), AT-Zoe (3) |
| all sources stacked | post-solstice | growth | -2.6 [-7.6, +2.4] | 0.304 | 36 | US-Ha1 (33), AT-Zoe (3) |
| all sources stacked | post-solstice | total sink (growth + respiration) | -2.1 [-5.1, +0.8] | 0.157 | 36 | US-Ha1 (33), AT-Zoe (3) |
| all sources stacked | post-solstice | NPP (GPP - respiration) | -0.5 [-6.1, +5.1] | 0.856 | 36 | US-Ha1 (33), AT-Zoe (3) |
| all sources stacked | post-solstice | GPP - total sink | -0.3 [-6.0, +5.3] | 0.909 | 36 | US-Ha1 (33), AT-Zoe (3) |
| all sources stacked | post-solstice | growth respiration (model) | -1.7 [-7.0, +3.6] | 0.536 | 36 | US-Ha1 (33), AT-Zoe (3) |
| all sources stacked | post-solstice | maintenance respiration (model) | -1.3 [-4.5, +2.0] | 0.448 | 36 | US-Ha1 (33), AT-Zoe (3) |
| all sources stacked | post-solstice | autotrophic respiration (model) | -1.7 [-5.4, +1.9] | 0.350 | 36 | US-Ha1 (33), AT-Zoe (3) |
| all sources stacked | post-solstice | total sink (growth + model respiration) | -2.0 [-5.7, +1.7] | 0.283 | 36 | US-Ha1 (33), AT-Zoe (3) |
| all sources stacked | post-solstice | GPP - total sink (model respiration) | -0.8 [-6.5, +4.8] | 0.774 | 36 | US-Ha1 (33), AT-Zoe (3) |

**What this shows**

- **Data.** 61 site-years of stem growth at 8 sites; 5 sites have at least three years. Only 37 site-years have readings dense enough to date the start and end of growth.
- **Seasonal course.** On average 43% of a year's stem growth is done by the summer solstice. Growth starts around day 143, is fastest around day 165 and is 90% complete by day 233 - weeks before the canopy starts to senesce (EOS90 around day 260-270). Stem growth and leaf senescence are separated in time.
- **A. Growth and GPP.** Annual growth follows annual GPP (within-site r = +0.42, p = 0.003). Growth by the solstice does not follow GPP before the solstice (r = +0.07) but follows spring temperature (r = +0.41, p = 0.019). Early in the season the sink runs on temperature, not on the current carbon supply - the same conclusion as section 14 D, from an independent measurement.
- **B. Growth and the onset of senescence.** All sources stacked (40 site-years, mostly US-Ha1): growth by the solstice +0.8, annual growth -2.1, pre-solstice GPP on the same years -2.4 days per SD; none can be distinguished from zero. Single sources disagree in sign. With one site supplying most of the years and band readings only 4-5 times a year there, this is not a test of the sink hypothesis yet.
- **Timing of growth and EOS90.** Only 13 site-years at 4 sites have both a dated growth curve and an EOS90 (sites with two years admitted). A later end of growth goes with a later onset of senescence (+3.2 days per SD, p = 0.31); so does a later peak of growth (+6.0). The direction is what a sink mechanism predicts (growth finished early, canopy shed early), but with a dozen site-years, two-year sites and mostly one EOS source, this is an indication to follow up, not a result. The p-values of this row are not reliable.
- **C. Carbon budget.** Over the growing season, AT-Zoe: GPP 983, respiration 315 (32% of GPP), measured growth 205, total sink 521, GPP minus total sink 462 gC m-2; US-Ha1: GPP 1566, respiration 353 (23% of GPP), measured growth 165, total sink 518, GPP minus total sink 1048 gC m-2. Measured growth is 33% of NPP at AT-Zoe, 14% of NPP at US-Ha1.
- The residual (GPP minus total sink) is large because the measured growth is aboveground only and because the respiration estimate is low: forests typically respire 45-55% of GPP. The absolute level of the residual should not be interpreted; its year-to-year variation may be.
- With respiration from the fitted model (growth + maintenance; the maintenance term follows temperature and accumulated biomass, not the day's GPP), respiration before the solstice is 69 gC m-2 at AT-Zoe (against 88), 75 gC m-2 at US-Ha1 (against 97): the year's total is the same by construction, but its seasonal distribution differs.
- **Budget terms before the solstice against EOS90** (stacked, days per SD): GPP -2.8, autotrophic respiration -2.6, growth +1.8, total sink (growth + respiration) -2.2, total sink (growth + model respiration) -1.4, NPP (GPP - respiration) -1.1, GPP - total sink -1.3, GPP - total sink (model respiration) -3.0*. Total sink behaves like GPP because respiration is its largest part and is derived from GPP. The terms that are independent of the day's GPP - measured growth, and the model-based sink - are the ones to watch, and they show no clear effect on this small sample.

![US-Ha1: cumulative stem growth, one line per year (step 29)](../figure/dendro_growth/US-Ha1.png)

*US-Ha1: cumulative stem growth, one line per year (step 29)*

![AT-Zoe: cumulative stem growth, one line per year (step 29)](../figure/dendro_growth/AT-Zoe.png)

*AT-Zoe: cumulative stem growth, one line per year (step 29)*

![CH-Dav: cumulative stem growth, one line per year (step 29)](../figure/dendro_growth/CH-Dav.png)

*CH-Dav: cumulative stem growth, one line per year (step 29)*

![Stem growth against GPP, within sites (step 45)](../figure/dendro_growth/growth_vs_gpp.png)

*Stem growth against GPP, within sites (step 45)*

![Onset of senescence against stem growth, within sites (step 45)](../figure/dendro_growth/growth_vs_eos90.png)

*Onset of senescence against stem growth, within sites (step 45)*

![Carbon budget through the season: GPP, respiration, growth, total sink and the residual (step 45)](../figure/dendro_growth/carbon_budget.png)

*Carbon budget through the season: GPP, respiration, growth, total sink and the residual (step 45)*

## 16. What the results add up to

1. **Early-season GPP and the onset of senescence.** A year with more GPP before the summer solstice has a slightly earlier onset of senescence at the same site. The effect is about 1-2.5 days per SD without covariates and about 1 day once spring temperature is held fixed. It is found in every EOS source in sign, survives the removal of any site and every QC setting, and is statistically secured only in the satellite record (sections 3, 13).
2. **Whole-season GPP** has little or no effect on the onset of senescence (sections 3, 4). In that sense Zani et al. and Lu et al. are both right.
3. **But not because two effects cancel.** Late-season GPP has no consistent delaying effect once the window-length artefact is removed (sections 3.4, 4). The early effect is simply small and is diluted in the seasonal total.
4. **The end of senescence (EOS10)** is not related to GPP in any window (sections 3.3, 6, 7, 10).
5. **Spring temperature** is a stronger predictor of the onset of senescence than GPP and accounts for about half of the uncorrected GPP effect (section 14 D). Stem growth before the solstice also follows temperature, not GPP (section 15).
6. **Deciduous forests:** leaf-out date, not GPP, predicts the onset of senescence - but only weakly, and not through a fixed leaf life span (senescence follows leaf-out by about a fifth of a day per day, section 14 E). **Dry-summer sites:** spring water balance, not GPP. The GPP effect that remains is mainly in evergreen forests (section 14).
7. **Sink variables from the flux data** (NPP, respiration) cannot separate sink from source, because they are GPP times a factor. Measured stem growth can, but exists for too few site-years to decide (sections 5, 14 C, 15).
8. **Methodological result.** Same-year windows and mixed models both inflate the apparent effects - the first by construction, the second by a factor of about two (sections 3.4, 6). Within-site models with fixed windows are needed for this question.

---

Per-step logs are in `logs/`. Method notes are in `README.md` and in the docstring of each script.

