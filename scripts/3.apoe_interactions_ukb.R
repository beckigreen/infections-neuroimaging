# ========================== #
# 2. APOE interactions: UKB  #
# ========================== #

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
# Serostatus
serology_dat_genqc <-
  readRDS("data_derived/serology_dat_genqc.RDS") %>%
  select(n_eid,
         PC1:PC10)

dat_serostatus <- readRDS("data_derived/ukb_serostatus.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(n_eid %in% serology_dat_genqc$n_eid) %>%
  inner_join(serology_dat_genqc, using = "n_eid") %>%
  as.mids(.)

dim(dat_serostatus$data)

# PBI
dat_seroreactivity <-
  readRDS("data_derived/ukb_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(n_eid %in% serology_dat_genqc$n_eid) %>%
  inner_join(serology_dat_genqc, using = "n_eid") %>%
  as.mids(.)

dim(dat_serostatus$data)

## Select list of predictors for analysis ----
# Serostatus
serostatus <- dat_serostatus$data %>%
  select(contains("_serostatus")) %>%
  select(!(contains("_v2"))) %>%
  colnames()

n_distinct(serostatus) # 17 pathogens

pbi <- c("PBI", "neuro_PBI")

## Set up models ----
m1_covar <-
  "+ TIV + scanner_x + scanner_y + scanner_z + centre_serology + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
m2_covar <- paste0(m1_covar,
                   " + serology_age + imaging_age + sex")
m3_covar <- paste0(m2_covar,
                   " + education + townsend_cat + smok_stat + bmi + alcohol_cat")

## Serostatus (APOE e2) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(serostatus,
             # regression analyses
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe2_stat", m1_covar)
               ))
             )))) %>%
  # bind pathogen results together
  rbindlist() %>%
  # select cols of interest
  select(pathogen = term, estimate, std.error, p.value) %>%
  # extract only serostatus res
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  # indicate which model
  mutate(model = "model 1")

#### Hippocampal volume ----
hv_m1 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(
                 paste("hippovol ~", .x, "*apoe2_stat", m1_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <- map(serostatus,
               ~ summary(pool(with(
                 dat_serostatus,
                 lm(formula = as.formula(
                   paste("log_wmhv ~", .x, "*apoe2_stat", m1_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe2_stat", m2_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(
                 paste("hippovol ~", .x, "*apoe2_stat", m2_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <- map(serostatus,
               ~ summary(pool(with(
                 dat_serostatus,
                 lm(formula = as.formula(
                   paste("log_wmhv ~", .x, "*apoe2_stat", m2_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe2_stat", m3_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(
                 paste("hippovol ~", .x, "*apoe2_stat", m3_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <- map(serostatus,
               ~ summary(pool(with(
                 dat_serostatus,
                 lm(formula = as.formula(
                   paste("log_wmhv ~", .x, "*apoe2_stat", m3_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 3")

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
  write.csv("results/ukb_serostatus_apoe2.csv", row.names = F)

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

## Serostatus (APOE e4) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe4_stat", m1_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 1")

#### Hippocampal volume ----
hv_m1 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(
                 paste("hippovol ~", .x, "*apoe4_stat", m1_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <- map(serostatus,
               ~ summary(pool(with(
                 dat_serostatus,
                 lm(formula = as.formula(
                   paste("log_wmhv ~", .x, "*apoe4_stat", m1_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe4_stat", m2_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(
                 paste("hippovol ~", .x, "*apoe4_stat", m2_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <- map(serostatus,
               ~ summary(pool(with(
                 dat_serostatus,
                 lm(formula = as.formula(
                   paste("log_wmhv ~", .x, "*apoe4_stat", m2_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe4_stat", m3_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus,
               lm(formula = as.formula(
                 paste("hippovol ~", .x, "*apoe4_stat", m3_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <- map(serostatus,
               ~ summary(pool(with(
                 dat_serostatus,
                 lm(formula = as.formula(
                   paste("log_wmhv ~", .x, "*apoe4_stat", m3_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 3")

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
  write.csv("results/ukb_serostatus_apoe4.csv", row.names = F)

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

## PBI (APOE e2) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe2_stat", m1_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 1")

#### Hippocampal volume ----
hv_m1 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("hippovol ~", .x, "*apoe2_stat", m1_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <- map(pbi,
               ~ summary(pool(with(
                 dat_seroreactivity,
                 lm(formula = as.formula(
                   paste("log_wmhv ~", .x, "*apoe2_stat", m1_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe2_stat", m2_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("hippovol ~", .x, "*apoe2_stat", m2_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <- map(pbi,
               ~ summary(pool(with(
                 dat_seroreactivity,
                 lm(formula = as.formula(
                   paste("log_wmhv ~", .x, "*apoe2_stat", m2_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe2_stat", m3_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("hippovol ~", .x, "*apoe2_stat", m3_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <- map(pbi,
               ~ summary(pool(with(
                 dat_seroreactivity,
                 lm(formula = as.formula(
                   paste("log_wmhv ~", .x, "*apoe2_stat", m3_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_stat")) %>%
  mutate(model = "model 3")

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
  write.csv("results/ukb_pbi_apoe2.csv", row.names = F)

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

## PBI (APOE e4) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe4_stat", m1_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 1")

#### Hippocampal volume ----
hv_m1 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("hippovol ~", .x, "*apoe4_stat", m1_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <- map(pbi,
               ~ summary(pool(with(
                 dat_seroreactivity,
                 lm(formula = as.formula(
                   paste("log_wmhv ~", .x, "*apoe4_stat", m1_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe4_stat", m2_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("hippovol ~", .x, "*apoe4_stat", m2_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <- map(pbi,
               ~ summary(pool(with(
                 dat_seroreactivity,
                 lm(formula = as.formula(
                   paste("log_wmhv ~", .x, "*apoe4_stat", m2_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe4_stat", m3_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("hippovol ~", .x, "*apoe4_stat", m3_covar)
               ))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <- map(pbi,
               ~ summary(pool(with(
                 dat_seroreactivity,
                 lm(formula = as.formula(
                   paste("log_wmhv ~", .x, "*apoe4_stat", m3_covar)
                 ))
               )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_stat")) %>%
  mutate(model = "model 3")

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
  write.csv("results/ukb_pbi_apoe4.csv", row.names = F)

# End ----
