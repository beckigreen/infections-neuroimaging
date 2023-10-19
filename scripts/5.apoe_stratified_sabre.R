# =========================== #
# 5. APOE2 stratified: SABRE  #
# =========================== #

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
eur_pcs <-
  readRDS("data_derived/pca/sabre_pcs_eur.RDS") %>%
  mutate(sabre_id = as.integer(sabre_id))

# Principal components (SA)
sa_pcs <-
  readRDS("data_derived/pca/sabre_pcs_sa.RDS")  %>%
  mutate(sabre_id = as.integer(sabre_id))

### e2 carrier ----
# Serostatus (EUR, e2 carrier)
dat_serostatus_eur_e2 <-
  readRDS("data_derived/sabre_serostatus.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(sabre_id %in% eur_pcs$sabre_id) %>%
  inner_join(eur_pcs, using = "sabre_id") %>%
  filter(apoe2_carrier_imp == 1) %>%
  as.mids(.)

dim(dat_serostatus_eur_e2$data) 

# Serostatus (SA, e2 carrier)
dat_serostatus_sa_e2 <-
  readRDS("data_derived/sabre_serostatus.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(sabre_id %in% sa_pcs$sabre_id) %>%
  inner_join(sa_pcs, using = "sabre_id") %>%
  filter(apoe2_carrier_imp == 1) %>%
  as.mids(.)

dim(dat_serostatus_sa_e2$data) 

# PBI (EUR, e2 carrier)
dat_seroreactivity_eur_e2 <-
  readRDS("data_derived/sabre_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(sabre_id %in% eur_pcs$sabre_id) %>%
  inner_join(eur_pcs, using = "sabre_id") %>%
  filter(apoe2_carrier_imp == 1) %>%
  as.mids(.)

dim(dat_seroreactivity_eur_e2$data) 

# PBI (SA, e2 carrier)
dat_seroreactivity_sa_e2 <-
  readRDS("data_derived/sabre_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(sabre_id %in% sa_pcs$sabre_id) %>%
  inner_join(sa_pcs, using = "sabre_id") %>%
  filter(apoe2_carrier_imp == 1) %>%
  as.mids(.)

dim(dat_seroreactivity_sa_e2$data)

### Non e2 carrier ----
# Serostatus (EUR, non-e2 carrier)
dat_serostatus_eur_noe2 <-
  readRDS("data_derived/sabre_serostatus.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(sabre_id %in% eur_pcs$sabre_id) %>%
  inner_join(eur_pcs, using = "sabre_id") %>%
  filter(apoe2_carrier_imp == 0) %>%
  as.mids(.)

dim(dat_serostatus_eur_noe2$data)

# Serostatus (SA, non-e2 carrier)
dat_serostatus_sa_noe2 <-
  readRDS("data_derived/sabre_serostatus.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(sabre_id %in% sa_pcs$sabre_id) %>%
  inner_join(sa_pcs, using = "sabre_id") %>%
  filter(apoe2_carrier_imp == 0) %>%
  as.mids(.)

dim(dat_serostatus_sa_noe2$data)

# PBI (EUR, non-e2 carrier)
dat_seroreactivity_eur_noe2 <-
  readRDS("data_derived/sabre_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(sabre_id %in% eur_pcs$sabre_id) %>%
  inner_join(eur_pcs, using = "sabre_id") %>%
  filter(apoe2_carrier_imp == 0) %>%
  as.mids(.)

dim(dat_seroreactivity_eur_noe2$data)

# PBI (SA, non-e2 carrier)
dat_seroreactivity_sa_noe2 <-
  readRDS("data_derived/sabre_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(sabre_id %in% sa_pcs$sabre_id) %>%
  inner_join(sa_pcs, using = "sabre_id") %>%
  filter(apoe2_carrier_imp == 0) %>%
  as.mids(.)

dim(dat_seroreactivity_sa_noe2$data) 

## Set up models ----
m1_covar <- "+ tiv + scanner"
m2_covar <- paste0(m1_covar,
                   " + age + sex") # SA & EUR all same self-reported ethnicity
m3_covar <- paste0(m2_covar,
                   " + education + sep + smoking + alcohol_cat + bmi")

## Select list of predictors for analysis ----
serostatus <- c("Tg_serostatus", "MCV_serostatus")
pbi <- c("PBI", "neuro_PBI")

## EUR ----
### APOE e2 carrier ----
#### Model 1 ----
# MCV & T.gondii
hv_m1_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_eur_e2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 1")

# PBI
hv_m1_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_eur_e2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 1")

# Combine
hv_m1 <- rbind(hv_m1_sero, hv_m1_pbi)

#### Model 2 ----
# MCV & T.gondii
hv_m2_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_eur_e2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 2")

# PBI
hv_m2_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_eur_e2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 2")

# Combine
hv_m2 <- rbind(hv_m2_sero, hv_m2_pbi)

#### Model 3 ----
# MCV & T.gondii
hv_m3_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_eur_e2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 3")

# PBI
hv_m3_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_eur_e2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 3")

# Combine
hv_m3 <- rbind(hv_m3_sero, hv_m3_pbi)

#### Combine ----
hv_e2 <- rbind(hv_m1, hv_m2, hv_m3) %>%
  mutate(outcome = "hippocampal volume",
         analysis = "e2 carrier",
         study = "SABRE (EUR)")

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

### APOE e2 non-carriers ----
#### Model 1 ----
# MCV & T.gondii
hv_m1_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_eur_noe2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 1")

# PBI
hv_m1_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_eur_noe2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 1")

# Combine
hv_m1 <- rbind(hv_m1_sero, hv_m1_pbi)

#### Model 2 ----
# MCV & T.gondii
hv_m2_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_eur_noe2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 2")

# PBI
hv_m2_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_eur_noe2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 2")

# Combine
hv_m2 <- rbind(hv_m2_sero, hv_m2_pbi)

#### Model 3 ----
# MCV & T.gondii
hv_m3_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_eur_noe2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 3")

# PBI
hv_m3_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_eur_noe2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 3")

# Combine
hv_m3 <- rbind(hv_m3_sero, hv_m3_pbi)

#### Combine ----
hv_noe2 <- rbind(hv_m1, hv_m2, hv_m3) %>%
  mutate(outcome = "hippocampal volume",
         analysis = "non-e2 carrier",
         study = "SABRE (EUR)")

### Save ----
rbind(hv_e2, hv_noe2) %>%
  mutate(pathogen = gsub(pathogen, pattern = "_serostatus1", rep = "")) %>%
  write.csv("results/sabre_e2_strat_eur.csv", row.names = F)

rm(
  "hv_m1",
  "hv_m1_sero",
  "hv_m1_pbi",
  "hv_m2",
  "hv_m2_sero",
  "hv_m2_pbi",
  "hv_m3",
  "hv_m3_sero",
  "hv_m3_pbi",
  "hv_e2",
  "hv_noe2"
)

## SA ----
### APOE e2 carrier ----
#### Model 1 ----
# MCV & T.gondii
hv_m1_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_sa_e2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 1")

# PBI
hv_m1_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_sa_e2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 1")

# Combine
hv_m1 <- rbind(hv_m1_sero, hv_m1_pbi)

#### Model 2 ----
# MCV & T.gondii
hv_m2_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_sa_e2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 2")

# PBI
hv_m2_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_sa_e2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 2")

# Combine
hv_m2 <- rbind(hv_m2_sero, hv_m2_pbi)

#### Model 3 ----
# MCV & T.gondii
hv_m3_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_sa_e2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 3")

# PBI
hv_m3_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_sa_e2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 3")

# Combine
hv_m3 <- rbind(hv_m3_sero, hv_m3_pbi)

#### Combine ----
hv_e2 <- rbind(hv_m1, hv_m2, hv_m3) %>%
  mutate(outcome = "hippocampal volume",
         analysis = "e2 carrier",
         study = "SABRE (SA)")

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

### APOE e2 non-carriers ----
#### Model 1 ----
# MCV & T.gondii
hv_m1_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_sa_noe2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 1")

# PBI
hv_m1_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_sa_noe2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 1")

# Combine
hv_m1 <- rbind(hv_m1_sero, hv_m1_pbi)

#### Model 2 ----
# MCV & T.gondii
hv_m2_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_sa_noe2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 2")

# PBI
hv_m2_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_sa_noe2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 2")

# Combine
hv_m2 <- rbind(hv_m2_sero, hv_m2_pbi)

#### Model 3 ----
# MCV & T.gondii
hv_m3_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_sa_noe2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 3")

# PBI
hv_m3_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_sa_noe2, lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  )))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 3")

# Combine
hv_m3 <- rbind(hv_m3_sero, hv_m3_pbi)

#### Combine ----
hv_noe2 <- rbind(hv_m1, hv_m2, hv_m3) %>%
  mutate(outcome = "hippocampal volume",
         analysis = "non-e2 carrier",
         study = "SABRE (SA)")

### Save ----
rbind(hv_e2, hv_noe2) %>%
  mutate(pathogen = gsub(pathogen, pattern = "_serostatus1", rep = "")) %>%
  write.csv("results/sabre_e2_strat_sa.csv", row.names = F)

## End ----
