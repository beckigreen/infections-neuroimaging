# ========================== #
# 5. APOE2 stratified: NSHD  #
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
library(tidyr)

## Read in data ----
# Principal components
pcs <- read.table("data_derived/pca/nshd.pca.evec") %>%
  select(
    id = V1,
    PC1 = V2,
    PC2 = V3,
    PC3 = V4,
    PC4 = V5,
    PC5 = V6,
    PC6 = V7,
    PC7 = V8,
    PC8 = V9,
    PC9 = V10,
    PC10 = V11
  ) %>%
  separate(id, into = c("id0", "nshd_id"), ":NSHD") %>% # format ID for merging
  select(-id0) %>%
  mutate(nshd_id = as.integer(nshd_id))

# Serostatus
dat_serostatus_e2 <- readRDS("data_derived/nshd_serostatus.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(nshd_id %in% pcs$nshd_id) %>%
  inner_join(pcs, using = "nshd_id") %>%
  filter(apoe2_carrier == 1) %>%
  filter(!(is.na(apoe2_carrier))) %>%
  as.mids(.)

dim(dat_serostatus_e2$data) 

dat_serostatus_noe2 <-
  readRDS("data_derived/nshd_serostatus.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(nshd_id %in% pcs$nshd_id) %>%
  inner_join(pcs, using = "nshd_id") %>%
  filter(apoe2_carrier == 0) %>%
  filter(!(is.na(apoe2_carrier))) %>%
  as.mids(.)

dim(dat_serostatus_noe2$data) 

# PBI
dat_seroreactivity_e2 <-
  readRDS("data_derived/nshd_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(nshd_id %in% pcs$nshd_id) %>%
  inner_join(pcs, using = "nshd_id") %>%
  filter(apoe2_carrier == 1) %>%
  as.mids(.)

dim(dat_seroreactivity_e2$data)

dat_seroreactivity_noe2 <-
  readRDS("data_derived/nshd_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(nshd_id %in% pcs$nshd_id) %>%
  inner_join(pcs, using = "nshd_id") %>%
  filter(apoe2_carrier == 0) %>%
  as.mids(.)

dim(dat_seroreactivity_noe2$data) 

## Set up models ----
m1_covar <- "+ tiv + serology_clin"
m2_covar <- paste0(m1_covar,
                   " + serology_age + imaging_age + sex")
m3_covar <- paste0(m2_covar,
                   " + education + sep_adult + eversmoker_comb + alcohol_cat + bmi")

## Select list of predictors for analysis ----
serostatus <- c("Tg_serostatus", "MCV_serostatus")
pbi <- c("PBI", "neuro_PBI")

## APOE e2 carriers ----
### Model 1 ----
# MCV & T.gondii
hv_m1_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_e2,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  )),
  weights =
    as.numeric(inf))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 1")

# PBI
hv_m1_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_e2,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  )),
  weights =
    as.numeric(inf))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 1")

# Combine
hv_m1 <- rbind(hv_m1_sero, hv_m1_pbi)

### Model 2 ----
# MCV & T.gondii
hv_m2_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_e2,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  )),
  weights =
    as.numeric(inf))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 2")

# PBI
hv_m2_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_e2,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  )),
  weights =
    as.numeric(inf))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 2")

# Combine
hv_m2 <- rbind(hv_m2_sero, hv_m2_pbi)

### Model 3 ----
# MCV & T.gondii
hv_m3_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_e2,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  )),
  weights =
    as.numeric(inf))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 3")

# PBI
hv_m3_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_e2,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  )),
  weights =
    as.numeric(inf))
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
         study = "NSHD")

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
hv_m1_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_noe2,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  )),
  weights =
    as.numeric(inf))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 1")

# PBI
hv_m1_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_noe2,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m1_covar
  )),
  weights =
    as.numeric(inf))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 1")

# Combine
hv_m1 <- rbind(hv_m1_sero, hv_m1_pbi)

### Model 2 ----
# MCV & T.gondii
hv_m2_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_noe2,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  )),
  weights =
    as.numeric(inf))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 2")

# PBI
hv_m2_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_noe2,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m2_covar
  )),
  weights =
    as.numeric(inf))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "PBI")) %>%
  mutate(model = "model 2")

# Combine
hv_m2 <- rbind(hv_m2_sero, hv_m2_pbi)

### Model 3 ----
# MCV & T.gondii
hv_m3_sero <- map(serostatus, ~ summary(pool(with(
  dat_serostatus_noe2,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  )),
  weights =
    as.numeric(inf))
)))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = "serostatus")) %>%
  mutate(model = "model 3")

# PBI
hv_m3_pbi <- map(pbi, ~ summary(pool(with(
  dat_seroreactivity_noe2,
  lm(formula = as.formula(paste(
    "hippovol ~", .x, m3_covar
  )),
  weights =
    as.numeric(inf))
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
         study = "NSHD")

## Save ----
rbind(hv_e2, hv_noe2) %>%
  mutate(pathogen = gsub(pathogen, pattern = "_serostatus1", rep = "")) %>%
  write.csv("results/nshd_e2_strat.csv", row.names = F)

## End ----
