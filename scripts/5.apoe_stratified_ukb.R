# ========================= #
# 5. APOE2 stratified: UKB  #
# ========================= #

## Set up ----
filepath <- ""
setwd(filepath)

## Load libraries ----
library(data.table)
library(dplyr)
library(magrittr)
library(mice)
library(purrr)
library(tidyr)

## Read in data ----
# Principal components
serology_dat_genqc <-
  readRDS("data_derived/serology_dat_genqc.RDS") %>%
  select(n_eid,
         PC1:PC10)

# Serostatus
dat_serostatus_e2 <- readRDS("data_derived/ukb_serostatus.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(n_eid %in% serology_dat_genqc$n_eid) %>%
  inner_join(serology_dat_genqc, using = "n_eid") %>%
  filter(apoe2_stat == 1) %>%
  as.mids(.)

dim(dat_serostatus_e2$data)

dat_serostatus_noe2 <-
  readRDS("data_derived/ukb_serostatus.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(n_eid %in% serology_dat_genqc$n_eid) %>%
  inner_join(serology_dat_genqc, using = "n_eid") %>%
  filter(apoe2_stat == 0) %>%
  as.mids(.)

dim(dat_serostatus_noe2$data)

# PBI
dat_seroreactivity_e2 <-
  readRDS("data_derived/ukb_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(n_eid %in% serology_dat_genqc$n_eid) %>%
  inner_join(serology_dat_genqc, using = "n_eid") %>%
  filter(apoe2_stat == 1) %>%
  as.mids(.)

dim(dat_seroreactivity_e2$data)

dat_seroreactivity_noe2 <-
  readRDS("data_derived/ukb_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(n_eid %in% serology_dat_genqc$n_eid) %>%
  inner_join(serology_dat_genqc, using = "n_eid") %>%
  filter(apoe2_stat == 0) %>%
  as.mids(.)

dim(dat_seroreactivity_noe2$data)

## Select list of predictors for analysis ----
serostatus <- c("Tg_serostatus", "MCV_serostatus")
pbi <- c("PBI", "neuro_PBI")

## Set up models ----
m1_covar <-
  "+ TIV + scanner_x + scanner_y + scanner_z + centre_serology"
m2_covar <- paste0(m1_covar,
                   " + serology_age + imaging_age + sex") # subset all same self-reported ethnicity
m3_covar <- paste0(m2_covar,
                   " + education + townsend_cat + smok_stat + bmi + alcohol_cat")

## APOE e2 carriers ----
### Model 1 ----
# MCV & T.gondii
hv_m1_sero <- map(serostatus,
                  ~ summary(pool(with(
                    dat_serostatus_e2,
                    lm(formula = as.formula(paste(
                      "hippovol ~", .x, m1_covar
                    )))
                  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 1")

hv_m1_pbi <- map(pbi,
                 ~ summary(pool(with(
                   dat_seroreactivity_e2,
                   lm(formula = as.formula(paste(
                     "hippovol ~", .x, m1_covar
                   )))
                 )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 1")

# Combine
hv_m1 <- rbind(hv_m1_sero, hv_m1_pbi)

### Model 2 ----
# MCV & T.gondii
hv_m2_sero <- map(serostatus,
                  ~ summary(pool(with(
                    dat_serostatus_e2,
                    lm(formula = as.formula(paste(
                      "hippovol ~", .x, m2_covar
                    )))
                  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 2")

hv_m2_pbi <- map(pbi,
                 ~ summary(pool(with(
                   dat_seroreactivity_e2,
                   lm(formula = as.formula(paste(
                     "hippovol ~", .x, m2_covar
                   )))
                 )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 2")

# Combine
hv_m2 <- rbind(hv_m2_sero, hv_m2_pbi)

### Model 3 ----
# MCV & T.gondii
hv_m3_sero <- map(serostatus,
                  ~ summary(pool(with(
                    dat_serostatus_e2,
                    lm(formula = as.formula(paste(
                      "hippovol ~", .x, m3_covar
                    )))
                  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 3")

hv_m3_pbi <- map(pbi,
                 ~ summary(pool(with(
                   dat_seroreactivity_e2,
                   lm(formula = as.formula(paste(
                     "hippovol ~", .x, m3_covar
                   )))
                 )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 3")

# Combine
hv_m3 <- rbind(hv_m3_sero, hv_m3_pbi)

### Combine ----
hv_e2 <- rbind(hv_m1, hv_m2, hv_m3) %>%
  mutate(outcome = "hippocampal volume",
         analysis = "e2 carrier",
         study = "UKB")

# Remove items for next analysis
rm(
  "hv_m1",
  "hv_m1_sero",
  "hv_m1_pbi",
  "hv_m2",
  "hv_m2_sero",
  "hv_m2_pbi",
  "hv_m3",
  "hv_m3_sero",
  "hv_m3_pbi"
)

## APOE e2 non-carriers ----
### Model 1 ----
# MCV & T.gondii
hv_m1_sero <- map(serostatus,
                  ~ summary(pool(with(
                    dat_serostatus_noe2,
                    lm(formula = as.formula(paste(
                      "hippovol ~", .x, m1_covar
                    )))
                  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 1")

hv_m1_pbi <- map(pbi,
                 ~ summary(pool(with(
                   dat_seroreactivity_noe2,
                   lm(formula = as.formula(paste(
                     "hippovol ~", .x, m1_covar
                   )))
                 )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 1")

# Combine
hv_m1 <- rbind(hv_m1_sero, hv_m1_pbi)

### Model 2 ----
# MCV & T.gondii
hv_m2_sero <- map(serostatus,
                  ~ summary(pool(with(
                    dat_serostatus_noe2,
                    lm(formula = as.formula(paste(
                      "hippovol ~", .x, m2_covar
                    )))
                  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 2")

hv_m2_pbi <- map(pbi,
                 ~ summary(pool(with(
                   dat_seroreactivity_noe2,
                   lm(formula = as.formula(paste(
                     "hippovol ~", .x, m2_covar
                   )))
                 )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 2")

# Combine
hv_m2 <- rbind(hv_m2_sero, hv_m2_pbi)

### Model 3 ----
# MCV & T.gondii
hv_m3_sero <- map(serostatus,
                  ~ summary(pool(with(
                    dat_serostatus_noe2,
                    lm(formula = as.formula(paste(
                      "hippovol ~", .x, m3_covar
                    )))
                  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 3")

hv_m3_pbi <- map(pbi,
                 ~ summary(pool(with(
                   dat_seroreactivity_noe2,
                   lm(formula = as.formula(paste(
                     "hippovol ~", .x, m3_covar
                   )))
                 )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 3")

# Combine
hv_m3 <- rbind(hv_m3_sero, hv_m3_pbi)

### Combine ----
hv_noe2 <- rbind(hv_m1, hv_m2, hv_m3) %>%
  mutate(outcome = "hippocampal volume",
         analysis = "non-e2 carrier",
         study = "UKB")

### Combine & save ----
rbind(hv_e2, hv_noe2) %>%
  mutate(pathogen = gsub(pathogen, pattern = "_serostatus1", replacement = "")) %>%
  write.csv("results/ukb_e2_strat.csv", row.names = F)

## End ----