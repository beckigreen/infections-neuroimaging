# ================= #
# 1. Analysis: UKB  #
# ================= #

## Set up ----
filepath <- ""
setwd(filepath)

## Load libraries ----
library(data.table)
library(dplyr)
library(magrittr)
library(mice)
library(purrr)

## Read in data ----
source("scripts/functions.R")

# Serostatus
dat_serostatus <- readRDS("data_derived/ukb_serostatus.RDS")

# Seroreactivity
# Convert seronegatives to NAs for seroreactivity analyses
dat_seroreactivity <-
  readRDS("data_derived/ukb_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  mutate(across(all_of(contains("_seroreac")),
                ~ na_if(., 0))) %>%
  as.mids(.)

## Select list of predictors for analysis ----
# Serostatus
serostatus <- dat_serostatus$data %>%
  select(contains("_serostatus")) %>%
  select(!(contains("_v2"))) %>%
  colnames()

n_distinct(serostatus) # 17 pathogens

# Seroreactivity
seroreac <- dat_seroreactivity$data %>%
  select(
    HSV1_1gG_seroreac,
    HSV2_2mgGunique_seroreac,
    VZV_gE_gI_seroreac,
    EBV_VCAp18_seroreac,
    CMV_pp150NTerm_seroreac,
    HHV6_1E1A_seroreac,
    HHV6_1E1B_seroreac,
    HHV7_U14_seroreac,
    BK_VP1_seroreac,
    JC_VP1_seroreac,
    MCV_VP1_seroreac,
    Tg_p22_seroreac,
    HP0547_1_seroreac,
    Ct_pGP3_seroreac
  ) %>%
  colnames()

n_distinct(seroreac) # 14 antigens

pbi <- c("PBI", "neuro_PBI")

## Set up models ----
m1_covar <- "+ TIV + scanner_x + scanner_y + scanner_z + centre_serology"
m2_covar <- paste0(m1_covar,
                   " + serology_age + imaging_age + ethnicity + sex")
m3_covar <- paste0(m2_covar,
                   " + education + townsend_cat + smok_stat + bmi + alcohol_cat")

## Serostatus ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(serostatus,
             # regression analyses
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(paste(
                 "brainvol ~", .x, m1_covar
               )))
             )))) %>%
  # bind pathogen results together
  rbindlist() %>%
  # select cols of interest
  select(pathogen = term, estimate, std.error, p.value) %>%
  # extract only serostatus res
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  # indicate which model
  mutate(model = "model 1")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_brainvol_lmchecks_m1.pdf")
map(serostatus, ~ mids_diagnostics(model = with(
  dat_serostatus,
  lm(formula = as.formula(paste(
    "brainvol ~", .x, m1_covar
  ))), name = .x
)))
dev.off()

#### Hippocampal volume ----
hv_m1 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(paste(
                 "hippovol ~", .x, m1_covar
               )))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 1")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_hippovol_lmchecks_m1.pdf")
map(serostatus, ~ mids_diagnostics(model = with(
  dat_serostatus,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  ))), name = .x
)))
dev.off()

#### WMHV ----
wmhv_m1 <- map(serostatus,
               ~ summary(pool(with(
                 dat_serostatus,
                 lm(formula = as.formula(paste(
                   "log_wmhv ~", .x, m1_covar
                 )))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 1")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_wmhv_lmchecks_m1.pdf")
map(serostatus, ~ mids_diagnostics(model = with(
  dat_serostatus,
  lm(formula = as.formula(paste(
    "log_wmhv ~", .x, m1_covar
  ))), name = .x
)))
dev.off()

### Model 2 ----
#### Brain volume ----
bv_m2 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(paste(
                 "brainvol ~", .x, m2_covar
               )))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 2")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_brainvol_lmchecks_m2.pdf")
map(serostatus, ~ mids_diagnostics(model = with(
  dat_serostatus,
  lm(formula = as.formula(paste(
    "brainvol ~", .x, m2_covar
  ))), name = .x
)))
dev.off()

#### Hippocampal volume ----
hv_m2 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(paste(
                 "hippovol ~", .x, m2_covar
               )))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 2")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_hippovol_lmchecks_m2.pdf")
map(serostatus, ~ mids_diagnostics(model = with(
  dat_serostatus,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  ))), name = .x
)))
dev.off()

#### WMHV ----
wmhv_m2 <- map(serostatus,
               ~ summary(pool(with(
                 dat_serostatus,
                 lm(formula = as.formula(paste(
                   "log_wmhv ~", .x, m2_covar
                 )))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 2")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_wmhv_lmchecks_m2.pdf")
map(serostatus, ~ mids_diagnostics(model = with(
  dat_serostatus,
  lm(formula = as.formula(paste(
    "log_wmhv ~", .x, m2_covar
  ))), name = .x
)))
dev.off()

### Model 3 ----
#### Brain volume ----
bv_m3 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(paste(
                 "brainvol ~", .x, m3_covar
               )))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 3")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_brainvol_lmchecks_m3.pdf")
map(serostatus, ~ mids_diagnostics(model = with(
  dat_serostatus,
  lm(formula = as.formula(paste(
    "brainvol ~", .x, m3_covar
  ))), name = .x
)))
dev.off()

#### Hippocampal volume ----
hv_m3 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(paste(
                 "hippovol ~", .x, m3_covar
               )))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 3")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_hippovol_lmchecks_m3.pdf")
map(serostatus, ~ mids_diagnostics(model = with(
  dat_serostatus,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  ))), name = .x
)))
dev.off()

#### WMHV ----
wmhv_m3 <- map(serostatus,
               ~ summary(pool(with(
                 dat_serostatus,
                 lm(formula = as.formula(paste(
                   "log_wmhv ~", .x, m3_covar
                 )))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 3")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_wmhv_lmchecks_m3.pdf")
map(serostatus, ~ mids_diagnostics(model = with(
  dat_serostatus,
  lm(formula = as.formula(paste(
    "log_wmhv ~", .x, m3_covar
  ))), name = .x
)))
dev.off()

### Combine & save ----
bv_res <- rbind(bv_m1,
                bv_m2,
                bv_m3) %>%
  mutate(outcome = "brain volume")

hv_res <- rbind(hv_m1,
                hv_m2,
                hv_m3) %>%
  mutate(outcome = "hippocampal volume")

wmhv_res <- rbind(wmhv_m1,
                  wmhv_m2,
                  wmhv_m3) %>%
  mutate(outcome = "WMHV")

rbind(bv_res, hv_res, wmhv_res) %>%
  mutate(pathogen = gsub(pathogen, pattern = "_serostatus1", replacement = "")) %>%
  write.csv("results/ukb_serostatus.csv", row.names = F)

# Remove items for next analysis
rm(
  "wmhv_res",
  "wmhv_m1",
  "wmhv_m2",
  "wmhv_m3",
  "hv_res",
  "hv_m1",
  "hv_m2",
  "hv_m3",
  "bv_res",
  "bv_m1",
  "bv_m2",
  "bv_m3"
)

## PBI ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(paste(
                 "brainvol ~", .x, m1_covar
               )))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 1")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_brainvol_lmchecks_m1_PBI.pdf")
map(pbi, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(paste(
    "brainvol ~", .x, m1_covar
  ))), name = .x
)))
dev.off()

#### Hippocampal volume ----
hv_m1 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(paste(
                 "hippovol ~", .x, m1_covar
               )))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 1")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_hippovol_lmchecks_m1_PBI.pdf")
map(pbi, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  ))), name = .x
)))
dev.off()

#### WMHV ----
wmhv_m1 <- map(pbi,
               ~ summary(pool(with(
                 dat_seroreactivity,
                 lm(formula = as.formula(paste(
                   "log_wmhv ~", .x, m1_covar
                 )))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 1")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_wmhv_lmchecks_m1_PBI.pdf")
map(pbi, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(paste(
    "log_wmhv ~", .x, m1_covar
  ))), name = .x
)))
dev.off()

### Model 2 ----
#### Brain volume ----
bv_m2 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(paste(
                 "brainvol ~", .x, m2_covar
               )))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 2")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_brainvol_lmchecks_m2_PBI.pdf")
map(pbi, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(paste(
    "brainvol ~", .x, m2_covar
  ))), name = .x
)))
dev.off()

#### Hippocampal volume ----
hv_m2 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(paste(
                 "hippovol ~", .x, m2_covar
               )))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 2")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_hippovol_lmchecks_m2_PBI.pdf")
map(pbi, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  ))), name = .x
)))
dev.off()

#### WMHV ----
wmhv_m2 <- map(pbi,
               ~ summary(pool(with(
                 dat_seroreactivity,
                 lm(formula = as.formula(paste(
                   "log_wmhv ~", .x, m2_covar
                 )))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 2")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_wmhv_lmchecks_m2_PBI.pdf")
map(pbi, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(paste(
    "log_wmhv ~", .x, m2_covar
  ))), name = .x
)))
dev.off()

### Model 3 ----
#### Brain volume ----
bv_m3 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(paste(
                 "brainvol ~", .x, m3_covar
               )))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 3")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_brainvol_lmchecks_m3_PBI.pdf")
map(pbi, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(paste(
    "brainvol ~", .x, m3_covar
  ))), name = .x
)))
dev.off()

#### Hippocampal volume ----
hv_m3 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(paste(
                 "hippovol ~", .x, m3_covar
               )))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 3")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_hippovol_lmchecks_m3_PBI.pdf")
map(pbi, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  ))), name = .x
)))
dev.off()

#### WMHV ----
wmhv_m3 <- map(pbi,
               ~ summary(pool(with(
                 dat_seroreactivity,
                 lm(formula = as.formula(paste(
                   "log_wmhv ~", .x, m3_covar
                 )))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 3")

# Diagnostics
pdf("figures/lmplots/serostatus/ukb_wmhv_lmchecks_m3_PBI.pdf")
map(pbi, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(paste(
    "log_wmhv ~", .x, m3_covar
  ))), name = .x
)))
dev.off()

### Combine & save ----
bv_res <- rbind(bv_m1,
                bv_m2,
                bv_m3) %>%
  mutate(outcome = "brain volume")

hv_res <- rbind(hv_m1,
                hv_m2,
                hv_m3) %>%
  mutate(outcome = "hippocampal volume")

wmhv_res <- rbind(wmhv_m1,
                  wmhv_m2,
                  wmhv_m3) %>%
  mutate(outcome = "WMHV")

rbind(bv_res, hv_res, wmhv_res) %>%
  write.csv("results/ukb_pbi.csv", row.names = F)

# Remove items for next analysis
rm(
  "wmhv_res",
  "wmhv_m1",
  "wmhv_m2",
  "wmhv_m3",
  "hv_res",
  "hv_m1",
  "hv_m2",
  "hv_m3",
  "bv_res",
  "bv_m1",
  "bv_m2",
  "bv_m3"
)

## Seroreactivity ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(seroreac,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("brainvol ~ factor(", .x, ", ordered=T)", m1_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ").L")) %>%
  mutate(
    pathogen = gsub(pathogen,
                    pattern = "factor|ordered = T|[(|.|)|,]", replacement = ""),
    model = "model 1"
  )

# Diagnostics
pdf("figures/lmplots/seroreactivity/ukb_brainvol_lmchecks_m1.pdf")
map(seroreac, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(
    paste("brainvol ~ factor(", .x, ", ordered=T)", m1_covar)
  ))
), name = .x))
dev.off()

#### Hippocampal volume ----
hv_m1 <- map(seroreac,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("hippovol ~ factor(", .x, ", ordered=T)", m1_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ").L")) %>%
  mutate(
    pathogen = gsub(pathogen, pattern = "factor|ordered = T|[(|.|)|,]", replacement = ""),
    model = "model 1"
  )

# Diagnostics
pdf("figures/lmplots/seroreactivity/ukb_hippovol_lmchecks_m1.pdf")
map(seroreac, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(
    paste("hippovol ~ factor(", .x, ", ordered=T)", m1_covar)
  ))
), name = .x))
dev.off()

#### WMHV ----
wmhv_m1 <- map(seroreac,
               ~ summary(pool(with(
                 dat_seroreactivity,
                 lm(formula = as.formula(
                   paste("log_wmhv ~ factor(", .x, ", ordered=T)", m1_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ").L")) %>%
  mutate(
    pathogen = gsub(pathogen, pattern = "factor|ordered = T|[(|.|)|,]", replacement = ""),
    model = "model 1"
  )

# Diagnostics
pdf("figures/lmplots/seroreactivity/ukb_wmhv_lmchecks_m1.pdf")
map(seroreac, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(
    paste("log_wmhv ~ factor(", .x, ", ordered=T)", m1_covar)
  ))
), name = .x))
dev.off()

### Model 2 ----
#### Brain volume ----
bv_m2 <- map(seroreac,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("brainvol ~ factor(", .x, ", ordered=T)", m2_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ").L")) %>%
  mutate(
    pathogen = gsub(pathogen, pattern = "factor|ordered = T|[(|.|)|,]", replacement = ""),
    model = "model 2"
  )

# Diagnostics
pdf("figures/lmplots/seroreactivity/ukb_brainvol_lmchecks_m2.pdf")
map(seroreac, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(
    paste("brainvol ~ factor(", .x, ", ordered=T)", m2_covar)
  ))
), name = .x))
dev.off()

#### Hippocampal volume ----
hv_m2 <- map(seroreac,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("hippovol ~ factor(", .x, ", ordered=T)", m2_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ").L")) %>%
  mutate(
    pathogen = gsub(pathogen, pattern = "factor|ordered = T|[(|.|)|,]", replacement = ""),
    model = "model 2"
  )

# Diagnostics
pdf("figures/lmplots/seroreactivity/ukb_hippovol_lmchecks_m2.pdf")
map(seroreac, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(
    paste("hippovol ~ factor(", .x, ", ordered=T)", m2_covar)
  ))
), name = .x))
dev.off()

#### WMHV ----
wmhv_m2 <- map(seroreac,
               ~ summary(pool(with(
                 dat_seroreactivity,
                 lm(formula = as.formula(
                   paste("log_wmhv ~ factor(", .x, ", ordered=T)", m2_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ").L")) %>%
  mutate(
    pathogen = gsub(pathogen, pattern = "factor|ordered = T|[(|.|)|,]", replacement = ""),
    model = "model 2"
  )

# Diagnostics
pdf("figures/lmplots/seroreactivity/ukb_wmhv_lmchecks_m2.pdf")
map(seroreac, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(
    paste("log_wmhv ~ factor(", .x, ", ordered=T)", m2_covar)
  ))
), name = .x))
dev.off()

### Model 3 ----
#### Brain volume ----
bv_m3 <- map(seroreac,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("brainvol ~ factor(", .x, ", ordered=T)", m3_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ").L")) %>%
  mutate(
    pathogen = gsub(pathogen, pattern = "factor|ordered = T|[(|.|)|,]", replacement = ""),
    model = "model 3"
  )

# Diagnostics
pdf("figures/lmplots/seroreactivity/ukb_brainvol_lmchecks_m3.pdf")
map(seroreac, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(
    paste("brainvol ~ factor(", .x, ", ordered=T)", m3_covar)
  ))
), name = .x))
dev.off()

#### Hippocampal volume ----
hv_m3 <- map(seroreac,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("hippovol ~ factor(", .x, ", ordered=T)", m3_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ").L")) %>%
  mutate(
    pathogen = gsub(pathogen, pattern = "factor|ordered = T|[(|.|)|,]", replacement = ""),
    model = "model 3"
  )

# Diagnostics
pdf("figures/lmplots/seroreactivity/ukb_hippovol_lmchecks_m3.pdf")
map(seroreac, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(
    paste("hippovol ~ factor(", .x, ", ordered=T)", m3_covar)
  ))
), name = .x))
dev.off()

#### WMHV ----
wmhv_m3 <- map(seroreac,
               ~ summary(pool(with(
                 dat_seroreactivity,
                 lm(formula = as.formula(
                   paste("log_wmhv ~ factor(", .x, ", ordered=T)", m3_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ").L")) %>%
  mutate(
    pathogen = gsub(pathogen, pattern = "factor|ordered = T|[(|.|)|,]", replacement = ""),
    model = "model 3"
  )

# Diagnostics
pdf("figures/lmplots/seroreactivity/ukb_wmhv_lmchecks_m3.pdf")
map(seroreac, ~ mids_diagnostics(model = with(
  dat_seroreactivity,
  lm(formula = as.formula(
    paste("log_wmhv ~ factor(", .x, ", ordered=T)", m3_covar)
  ))
), name = .x))
dev.off()

### Combine & save ----
bv_res <- rbind(bv_m1,
                bv_m2,
                bv_m3) %>%
  mutate(outcome = "brain volume")

hv_res <- rbind(hv_m1,
                hv_m2,
                hv_m3) %>%
  mutate(outcome = "hippocampal volume")

wmhv_res <- rbind(wmhv_m1,
                  wmhv_m2,
                  wmhv_m3) %>%
  mutate(outcome = "WMHV")

rbind(bv_res, hv_res, wmhv_res) %>%
  mutate(pathogen = gsub(pathogen, pattern = " L", replacement = "")) %>%
  write.csv("results/ukb_seroreactivity.csv", row.names = F)

## End ----