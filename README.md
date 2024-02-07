# [**Common infections and neuroimaging markers of dementia in three UK cohort studies - Green et al., 2024**](https://doi.org/10.1002/alz.13613)

_Rebecca E Green, Carole H Sudre, Charlotte Warren-Gash, Julia Butt, Tim Waterboer, Alun D Hughes, Jonathan M Schott, Marcus Richards, Nish Chaturvedi, Dylan M Williams_

---

## Background

This repo contains analysis code for the manuscript [Common infections and neuroimaging markers of dementia in three UK cohort studies](https://doi.org/10.1002/alz.13613) published in _Alzheimer's & Dementia_.

#### Repo structure:

```
.
├── README.md
├── scripts
|   ├── functions.R # user defined functions
│   ├── 1.analysis_nshd.R # main serostatus/pbi/seroreactivity analyses (NSHD)
│   ├── 1.analysis_sabre.R # main serostatus/pbi/seroreactivity analyses (SABRE)
│   ├── 1.analysis_ukb.R # main serostatus/pbi/seroreactivity analyses (UKB)
│   ├── 2.analysis_main_ma.R # meta-analyses of results from step 1
│   ├── 3.apoe_interactions_nshd.R # testing apoe e2&e4 interactions (NSHD)
│   ├── 3.apoe_interactions_sabre.R # testing apoe e2&e4 interactions (SABRE)
│   ├── 3.apoe_interactions_ukb.R # testing apoe e2&e4 interactions (UKB)
│   ├── 4.apoe_interactions_ma.R # meta-analyses of results from step 3
│   ├── 5.apoe_stratified_nshd.R # stratification by e2 carrier status for findings from step 4 (NSHD)
│   ├── 5.apoe_stratified_sabre.R # stratification by e2 carrier status for findings from step 4 (SABRE)
│   ├── 5.apoe_stratified_ukb.R # stratification by e2 carrier status for findings from step 4 (UKB)
│   ├── 6.apoe_stratified_ma.R # meta-analyses of results from step 5
│   ├── 7.figure2.R
│   └── 8.figure3.R
├── docs
│   └── figure1.pdf - analysis workflow figure
└── LICENSE.md

```

--- 

### Abstract 
**INTRODUCTION**
We aimed to investigate associations between common infections and neuroimaging markers of dementia risk (brain volume, hippocampal volume, white matter lesions) across three population-based studies.

**METHODS**
We tested associations between serology measures (pathogen serostatus, cumulative burden, continuous antibody responses) and outcomes using linear regression, including adjustments for total intracranial volume and scanner/clinic information (basic model), age, sex, ethnicity, education, socioeconomic position, alcohol, body mass index, and smoking (fully adjusted model). Interactions between serology measures and apolipoprotein E (APOE) genotype were tested. Findings were meta-analysed across cohorts (N<sub>main</sub> = 2632; N<sub>APOE-interaction</sub> = 1810).

**RESULTS**
Seropositivity to John Cunningham virus associated with smaller brain volumes in basic models (β = −3.89 mL [−5.81, −1.97], P<sub>adjusted</sub> < 0.05); these were largely attenuated in fully adjusted models (β = −1.59 mL [−3.55, 0.36], P = 0.11). No other relationships were robust to multiple testing corrections and sensitivity analyses, but several suggestive associations were observed.

**DISCUSSION**
We did not find clear evidence for relationships between common infections and markers of dementia risk. Some suggestive findings warrant testing for replication.

---

### Project workflow
<img src="docs/figure1.png" width="50%"/>


---
