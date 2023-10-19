# =========================== #
# 3. APOE interactions: NSHD  #
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
dat_serostatus <- readRDS("data_derived/nshd_serostatus.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(nshd_id %in% pcs$nshd_id) %>%
  inner_join(pcs, using = "nshd_id") %>%
  filter(!(is.na(apoe2_carrier))) %>%
  as.mids(.)

dim(dat_serostatus$data) 

# PBI
dat_seroreactivity <-
  readRDS("data_derived/nshd_seroreactivity.RDS") %>%
  complete(action = "long", include = T) %>%
  filter(nshd_id %in% pcs$nshd_id) %>%
  inner_join(pcs, using = "nshd_id") %>%
  filter(!(is.na(apoe2_carrier))) %>%
  as.mids(.)

dim(dat_seroreactivity$data) 

## Select list of predictors for analysis ----
# Serostatus
serostatus <- dat_serostatus$data %>%
  select(contains("_serostatus")) %>%
  select(-HHV6_serostatus) %>% #
  colnames()

n_distinct(serostatus) # 17 pathogens

pbi <- c("PBI", "neuro_PBI")

## Set up models ----
m1_covar <-
  "+ tiv + serology_clin + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10"
m2_covar <- paste0(m1_covar,
                   " + serology_age + imaging_age + sex")
m3_covar <- paste0(m2_covar,
                   " + education + sep_adult + eversmoker_comb + alcohol_cat + bmi")

## Serostatus (APOE e2) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(serostatus,
             # regression analyses
             ~ summary(pool(with(
               dat_serostatus, lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe2_carrier ", m1_covar)
               ),
               weights =
                 as.numeric(inf))
             )))) %>%
  # bind pathogen results together
  rbindlist() %>%
  # select cols of interest
  select(pathogen = term, estimate, std.error, p.value) %>%
  # extract only interaction res
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 1") # indicate which model

#### Hippocampal volume ----
hv_m1 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier ", m1_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier ", m1_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe2_carrier ", m2_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier ", m2_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier ", m2_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe2_carrier ", m3_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier ", m3_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier ", m3_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
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
  write.csv("results/nshd_serostatus_apoe2.csv", row.names = F)

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

## Serostatus (APOE e4) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(serostatus,
             ~ summary(pool(with(
               dat_serostatus, lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe4_carrier ", m1_covar)
               ),
               weights =
                 as.numeric(inf))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 1")

#### Hippocampal volume ----
hv_m1 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier ", m1_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier ", m1_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe4_carrier ", m2_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier ", m2_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier ", m2_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe4_carrier ", m3_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier ", m3_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <-
  map(serostatus, ~ summary(pool(with(
    dat_serostatus, lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier ", m3_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
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
  write.csv("results/nshd_serostatus_apoe4.csv", row.names = F)

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

## PBI (APOE e2) ----
### Model 1 ----
#### Brain volume ----
bv_m1 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe2_carrier ", m1_covar)
               ),
               weights =
                 as.numeric(inf))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 1")

#### Hippocampal volume ----
hv_m1 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier ", m1_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier ", m1_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe2_carrier ", m2_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier ", m2_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier ", m2_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe2_carrier ", m3_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe2_carrier ", m3_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe2_carrier ", m3_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe2_carrier")) %>%
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
  write.csv("results/nshd_pbi_apoe2.csv", row.names = F)

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

## PBI (APOE e4) ----
### Model 1 ----
#### Brain volume ----

bv_m1 <- map(pbi,
             ~ summary(pool(with(
               dat_seroreactivity,
               lm(formula = as.formula(
                 paste("brainvol ~", .x, "*apoe4_carrier ", m1_covar)
               ),
               weights =
                 as.numeric(inf))
             )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 1")

#### Hippocampal volume ----
hv_m1 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier ", m1_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 1")

#### WMHV ----
wmhv_m1 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier ", m1_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 1")

### Model 2 ----
#### Brain volume ----
bv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe4_carrier ", m2_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 2")

#### Hippocampal volume ----
hv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier ", m2_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 2")

#### WMHV ----
wmhv_m2 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier ", m2_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 2")

### Model 3 ----
#### Brain volume ----
bv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("brainvol ~", .x, "*apoe4_carrier ", m3_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 3")

#### Hippocampal volume ----
hv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("hippovol ~", .x, "*apoe4_carrier ", m3_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
  mutate(model = "model 3")

#### WMHV ----
wmhv_m3 <-
  map(pbi, ~ summary(pool(with(
    dat_seroreactivity,
    lm(formula = as.formula(
      paste("log_wmhv ~", .x, "*apoe4_carrier ", m3_covar)
    ),
    weights =
      as.numeric(inf))
  )))) %>%
  rbindlist() %>%
  select(pathogen = term, estimate, std.error, p.value) %>%
  filter(grepl(pathogen, pattern = ":apoe4_carrier")) %>%
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
  write.csv("results/nshd_pbi_apoe4.csv", row.names = F)

# End ----