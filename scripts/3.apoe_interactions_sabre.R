# ============================ #
# 3. APOE interactions: SABRE  #
# ============================ #

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
# Principal components (EUR)
eur_pcs <- readRDS("data_derived/pca/sabre_pcs_eur.RDS") %>%
  mutate(sabre_id = as.integer(sabre_id))

# Principal components (SA)
sa_pcs <- readRDS("data_derived/pca/sabre_pcs_sa.RDS")  %>%
  mutate(sabre_id = as.integer(sabre_id))

# Serostatus (EUR)
dat_serostatus_eur <- readRDS("data_derived/sabre_serostatus.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(sabre_id %in% eur_pcs$sabre_id) %>%
  inner_join(eur_pcs, using="sabre_id") %>%
  as.mids(.)

dim(dat_serostatus_eur$data)

# Serostatus (SA)
dat_serostatus_sa <- readRDS("data_derived/sabre_serostatus.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(sabre_id %in% sa_pcs$sabre_id) %>%
  inner_join(sa_pcs, using="sabre_id") %>%
  as.mids(.)

dim(dat_serostatus_sa$data) 

# PBI
dat_seroreactivity_eur <-
  readRDS("data_derived/sabre_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(sabre_id %in% eur_pcs$sabre_id) %>%
  inner_join(eur_pcs, using="sabre_id") %>%
  as.mids(.)

dim(dat_seroreactivity_eur$data) 

dat_seroreactivity_sa <-
  readRDS("data_derived/sabre_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(sabre_id %in% sa_pcs$sabre_id) %>%
  inner_join(sa_pcs, using="sabre_id") %>%
  as.mids(.)

dim(dat_seroreactivity_sa$data) 

## Select list of predictors for analysis ----
# Serostatus
serostatus <- dat_serostatus_sa$data %>%
  select(contains("_serostatus")) %>%
  select(-HHV6_serostatus) %>% #
  colnames()

n_distinct(serostatus) # 17 pathogens

pbi <- c("PBI", "neuro_PBI")

## Set up models ----
m1_covar <-
  "+ tiv + scanner + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
m2_covar <- paste0(m1_covar,
                   " + age + sex")
m3_covar <- paste0(m2_covar,
                   " + education + sep + smoking + alcohol_cat + bmi")

## Serostatus (Eur, APOE e2) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(serostatus,
             # regression analyses
             ~ summary(pool(with(
               dat_serostatus_eur, lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe2_carrier_imp ", m1_covar)
               )))))) %>%
  # bind pathogen results together
  rbindlist() %>%
  # select cols of interest
  select(pathogen = term, estimate, std.error, p.value) %>%
  # extract only interaction res
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 1") # indicate which model

#### Hippocampal volume ----
hv_m1 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <-
  map(serostatus[-11], ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 1")
# 11 won't run for wmhv models in eur

### Model 2 ----
#### Brain volume ----
bv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe2_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <-
  map(serostatus[-11], ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe2_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <-
  map(serostatus[-11], ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
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
  mutate(pathogen = gsub(pathogen, pattern = "_serostatus1|_imp", replacement = "")) %>%
  write.csv("results/sabre_serostatus_apoe2_eur.csv", row.names = F)

# Remove items to recreate in next analysis
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

## Serostatus (Eur, APOE e4) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus_eur, lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe4_carrier_imp ", m1_covar)
               )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 1") 

#### Hippocampal volume ----
hv_m1 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <-
  map(serostatus[-11], ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe4_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <-
  map(serostatus[-11], ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe4_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <-
  map(serostatus[-11], ~ summary(pool(with(
    dat_serostatus_eur, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
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
  mutate(pathogen = gsub(pathogen, pattern = "_serostatus1|_imp", replacement = "")) %>%
  write.csv("results/sabre_serostatus_apoe4_eur.csv", row.names = F)

# Remove items to recreate in next analysis
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

## Serostatus (SA, APOE e2) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus_sa, lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe2_carrier_imp ", m1_covar)
               )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 1") 

#### Hippocampal volume ----
hv_m1 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe2_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe2_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
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
  mutate(pathogen = gsub(pathogen, pattern = "_serostatus1|_imp", replacement = "")) %>%
  write.csv("results/sabre_serostatus_apoe2_sa.csv", row.names = F)

# Remove items to recreate in next analysis
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

## Serostatus (SA, APOE e4) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus_sa, lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe4_carrier_imp ", m1_covar)
               )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 1") 

#### Hippocampal volume ----
hv_m1 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe4_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe4_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus_sa, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
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
  mutate(pathogen = gsub(pathogen, pattern = "_serostatus1|_imp", replacement = "")) %>%
  write.csv("results/sabre_serostatus_apoe4_sa.csv", row.names = F)

# Remove items to recreate in next analysis
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

## PBI (Eur, APOE e2) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity_eur, lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe2_carrier_imp ", m1_covar)
               )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 1") 

#### Hippocampal volume ----
hv_m1 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe2_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe2_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
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
  mutate(pathogen = gsub(pathogen, pattern = "_imp", replacement = "")) %>%
  write.csv("results/sabre_pbi_apoe2_eur.csv", row.names = F)

# Remove items to recreate in next analysis
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

## PBI (Eur, APOE e4) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity_eur, lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe4_carrier_imp ", m1_covar)
               )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 1") 

#### Hippocampal volume ----
hv_m1 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe4_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe4_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_eur, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
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
  mutate(pathogen = gsub(pathogen, pattern = "_imp", replacement = "")) %>%
  write.csv("results/sabre_pbi_apoe4_eur.csv", row.names = F)

# Remove items to recreate in next analysis
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

## PBI (SA, APOE e2) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity_sa, lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe2_carrier_imp ", m1_covar)
               )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 1") 

#### Hippocampal volume ----
hv_m1 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe2_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe2_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier_imp")) %>%
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
  mutate(pathogen = gsub(pathogen, pattern = "_imp", replacement = "")) %>%
  write.csv("results/sabre_pbi_apoe2_sa.csv", row.names = F)

# Remove items to recreate in next analysis
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

## PBI (SA, APOE e4) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity_sa, lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe4_carrier_imp ", m1_covar)
               )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 1") 

#### Hippocampal volume ----
hv_m1 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier_imp ", m1_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe4_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier_imp ", m2_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe4_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity_sa, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier_imp ", m3_covar)
    )))))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier_imp")) %>%
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
  mutate(pathogen = gsub(pathogen, pattern = "_imp", replacement = "")) %>%
  write.csv("results/sabre_pbi_apoe4_sa.csv", row.names = F)

# End ----
