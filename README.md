# [**Common infections and neuroimaging markers of dementia in three UK cohort studies**](https://doi.org/10.1101/2023.07.12.23292538)

_Rebecca E Green, Carole H Sudre, Charlotte Warren-Gash, Julia Butt, Tim Waterboer, Alun D Hughes, Jonathan M Schott, Marcus Richards, Nish Chaturvedi, Dylan M Williams_

---

## Background

This repo contains analysis code for the preprint [Common infections and neuroimaging markers of dementia in three UK cohort studies](https://doi.org/10.1101/2023.07.12.23292538) - _Green et al_ 2023. 

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

### Project workflow
<img src="docs/figure1.png" width="50%"/>
  
---
