# ================== #
#  2. Meta-Analyses  #
# ================== #

## Set up ----
filepath <- ""
setwd(filepath)

## Load libraries ----
library(data.table)
library(dplyr)
library(magrittr)
library(metafor)
library(purrr)

## Read in data (serostatus & PBI) ----
nshd <- read.csv("results/nshd_serostatus.csv")
nshd_pbi <- read.csv("results/nshd_pbi.csv")
nshd %<>%
  rbind(nshd_pbi) %>%
  mutate(study = "NSHD")

sabre <- read.csv("results/sabre_serostatus.csv")
sabre_pbi <- read.csv("results/sabre_pbi.csv")
sabre %<>%
  rbind(sabre_pbi) %>%
  mutate(study = "SABRE")

ukb <- read.csv("results/ukb_serostatus.csv")
ukb_pbi <- read.csv("results/ukb_pbi.csv")
ukb %<>%
  rbind(ukb_pbi) %>%
  mutate(study = "UKB")

## Combine data ----
dat <- rbind(nshd, sabre, ukb)

# Check
dim(dat)
n_distinct(dat$pathogen)
table(dat$pathogen) 
table(dat$study) 
table(dat$model) 
table(dat$outcome) 
table(duplicated(dat$estimate)) 
table(duplicated(dat$std.error))
table(duplicated(dat$p.value)) 

## Meta-analyses ----
# Pathogens (& PBI)
pathogen <- dat %>%
  pull(pathogen) %>%
  unique()

models <- c("model 1", "model 2", "model 3")

### Brain volume ----
meta_df_list <- list()

for (i in pathogen) {
  for (j in models) {
    out <- "brain volume"
    print(paste0(out, ": ", i, ", model: ", j))
    res <-
      dat %>% filter(model == j &
                       outcome == out) %>% filter(pathogen == i)
    print(paste("Number of rows of filtered data:", nrow(res)))
    print(paste("Number of studies:", n_distinct(res$study)))
    print(res)
    rma <- rma(
      yi = estimate,
      vi = std.error ^ 2,
      method = "REML",
      data = res
    )
    meta_df <- data.frame(
      i,
      out,
      j,
      rma$beta,
      rma$zval,
      rma$pval,
      rma$ci.lb,
      rma$ci.ub,
      rma$se,
      rma$I2,
      rma$QE,
      rma$QEp,
      row.names = NULL
    )
    names(meta_df) <- c(
      "pathogen",
      "outcome",
      "model",
      "beta",
      "z",
      "p",
      "ci_lower",
      "ci_upper",
      "se",
      "i2",
      "Q",
      "Q-p"
    )
    meta_df_list <- c(meta_df_list, list(meta_df))
  }
}

bv <- rbindlist(meta_df_list) %>%
  mutate(padj = p.adjust(p, "BH"))

### Hippocampal volume ----
meta_df_list <- list()

for (i in pathogen) {
  for (j in models) {
    out <- "hippocampal volume"
    print(paste0(out, ": ", i, ", model: ", j))
    res <-
      dat %>% filter(model == j &
                       outcome == out) %>% filter(pathogen == i)
    print(paste("Number of rows of filtered data:", nrow(res)))
    print(paste("Number of studies:", n_distinct(res$study)))
    print(res)
    rma <- rma(
      yi = estimate,
      vi = std.error ^ 2,
      method = "REML",
      data = res
    )
    meta_df <- data.frame(
      i,
      out,
      j,
      rma$beta,
      rma$zval,
      rma$pval,
      rma$ci.lb,
      rma$ci.ub,
      rma$se,
      rma$I2,
      rma$QE,
      rma$QEp,
      row.names = NULL
    )
    names(meta_df) <- c(
      "pathogen",
      "outcome",
      "model",
      "beta",
      "z",
      "p",
      "ci_lower",
      "ci_upper",
      "se",
      "i2",
      "Q",
      "Q-p"
    )
    meta_df_list <- c(meta_df_list, list(meta_df))
  }
}

hv <- rbindlist(meta_df_list) %>%
  mutate(padj = p.adjust(p, "BH"))

### WMHV ----
meta_df_list <- list()

for (i in pathogen) {
  for (j in models) {
    out <- "WMHV"
    print(paste0(out, ": ", i, ", model: ", j))
    res <-
      dat %>% filter(model == j &
                       outcome == out) %>% filter(pathogen == i)
    print(paste("Number of rows of filtered data:", nrow(res)))
    print(paste("Number of studies:", n_distinct(res$study)))
    print(res)
    rma <- rma(
      yi = estimate,
      vi = std.error ^ 2,
      method = "REML",
      data = res
    )
    meta_df <- data.frame(
      i,
      out,
      j,
      rma$beta,
      rma$zval,
      rma$pval,
      rma$ci.lb,
      rma$ci.ub,
      rma$se,
      rma$I2,
      rma$QE,
      rma$QEp,
      row.names = NULL
    )
    names(meta_df) <- c(
      "pathogen",
      "outcome",
      "model",
      "beta",
      "z",
      "p",
      "ci_lower",
      "ci_upper",
      "se",
      "i2",
      "Q",
      "Q-p"
    )
    meta_df_list <- c(meta_df_list, list(meta_df))
  }
}

wmhv <- rbindlist(meta_df_list) %>%
  mutate(padj = p.adjust(p, "BH"))

## Save ----
rbind(bv, hv, wmhv) %>%
  write.csv("results/serostatus_main_MA.csv", row.names = F)

# Clear environment
rm(list = ls())

## Read in data (seroreactivity) ----
# UKB antigens are named differently to those in NSHD/SABRE
# So we need to match them so that they can be MA
nshd <- read.csv("results/nshd_seroreactivity.csv") %>%
  mutate(
    study = "NSHD",
    pathogen = case_when(
      pathogen == "HHV6_IE1A_truncated_seroreac" ~ "HHV6_1E1A_seroreac",
      pathogen == "HHV6_IE1B_truncated_seroreac" ~ "HHV6_1E1B_seroreac",
      pathogen == "MCV344_VP1_seroreac" ~ "MCV_VP1_seroreac",
      pathogen == "Ct_pGP3DCT_seroreac" ~ "Ct_pGP3_seroreac",
      TRUE ~ as.character(pathogen)
    )
  )

sabre <- read.csv("results/sabre_seroreactivity.csv") %>%
  mutate(
    study = "SABRE",
    pathogen = case_when(
      pathogen == "HHV6_IE1A_truncated_seroreac" ~ "HHV6_1E1A_seroreac",
      pathogen == "HHV6_IE1B_truncated_seroreac" ~ "HHV6_1E1B_seroreac",
      pathogen == "MCV344_VP1_seroreac" ~ "MCV_VP1_seroreac",
      pathogen == "Ct_pGP3DCT_seroreac" ~ "Ct_pGP3_seroreac",
      TRUE ~ as.character(pathogen)
    )
  )

ukb <- read.csv("results/ukb_seroreactivity.csv") %>%
  mutate(study = "UKB")

## Combine data ----
dat <- rbind(nshd, sabre, ukb)

# Check
dim(dat)
n_distinct(dat$pathogen) 
table(dat$pathogen)
table(dat$study)
table(dat$model) 
table(dat$outcome) 
table(duplicated(dat$estimate)) 
table(duplicated(dat$std.error)) 
table(duplicated(dat$p.value)) 

## Meta-analyses ----
pathogen <- dat %>%
  pull(pathogen) %>%
  unique()

models <- c("model 1", "model 2", "model 3")

### Brain volume ----
meta_df_list <- list()

for (i in pathogen) {
  for (j in models) {
    out <- "brain volume"
    print(paste0(out, ": ", i, ", model: ", j))
    res <-
      dat %>% filter(model == j &
                       outcome == out) %>% filter(pathogen == i)
    print(paste("Number of rows of filtered data:", nrow(res)))
    print(paste("Number of studies:", n_distinct(res$study)))
    print(res)
    rma <-
      rma(
        yi = estimate,
        vi = std.error ^ 2,
        method = "REML",
        data = res,
        control = list(maxiter = 1000)
      )
    meta_df <- data.frame(
      i,
      out,
      j,
      rma$beta,
      rma$zval,
      rma$pval,
      rma$ci.lb,
      rma$ci.ub,
      rma$se,
      rma$I2,
      rma$QE,
      rma$QEp,
      row.names = NULL
    )
    names(meta_df) <- c(
      "pathogen",
      "outcome",
      "model",
      "beta",
      "z",
      "p",
      "ci_lower",
      "ci_upper",
      "se",
      "i2",
      "Q",
      "Q-p"
    )
    meta_df_list <- c(meta_df_list, list(meta_df))
  }
}

bv <- rbindlist(meta_df_list)

bv_padj <- bv %>%
  mutate(padj = p.adjust(p, "BH"))

### Hippocampal volume ----
meta_df_list <- list()

for (i in pathogen) {
  for (j in models) {
    out <- "hippocampal volume"
    print(paste0(out, ": ", i, ", model: ", j))
    res <-
      dat %>% filter(model == j &
                       outcome == out) %>% filter(pathogen == i)
    print(paste("Number of rows of filtered data:", nrow(res)))
    print(paste("Number of studies:", n_distinct(res$study)))
    print(res)
    rma <- rma(
      yi = estimate,
      vi = std.error ^ 2,
      method = "REML",
      data = res
    )
    meta_df <- data.frame(
      i,
      out,
      j,
      rma$beta,
      rma$zval,
      rma$pval,
      rma$ci.lb,
      rma$ci.ub,
      rma$se,
      rma$I2,
      rma$QE,
      rma$QEp,
      row.names = NULL
    )
    names(meta_df) <- c(
      "pathogen",
      "outcome",
      "model",
      "beta",
      "z",
      "p",
      "ci_lower",
      "ci_upper",
      "se",
      "i2",
      "Q",
      "Q-p"
    )
    meta_df_list <- c(meta_df_list, list(meta_df))
  }
}

hv <- rbindlist(meta_df_list)

hv_padj <- hv %>%
  mutate(padj = p.adjust(p, "BH"))

### WMHV ----
meta_df_list <- list()

for (i in pathogen) {
  for (j in models) {
    out <- "WMHV"
    print(paste0(out, ": ", i, ", model: ", j))
    res <-
      dat %>% filter(model == j &
                       outcome == out) %>% filter(pathogen == i)
    print(paste("Number of rows of filtered data:", nrow(res)))
    print(paste("Number of studies:", n_distinct(res$study)))
    print(res)
    rma <- rma(
      yi = estimate,
      vi = std.error ^ 2,
      method = "REML",
      data = res
    )
    meta_df <- data.frame(
      i,
      out,
      j,
      rma$beta,
      rma$zval,
      rma$pval,
      rma$ci.lb,
      rma$ci.ub,
      rma$se,
      rma$I2,
      rma$QE,
      rma$QEp,
      row.names = NULL
    )
    names(meta_df) <- c(
      "pathogen",
      "outcome",
      "model",
      "beta",
      "z",
      "p",
      "ci_lower",
      "ci_upper",
      "se",
      "i2",
      "Q",
      "Q-p"
    )
    meta_df_list <- c(meta_df_list, list(meta_df))
  }
}

wmhv <- rbindlist(meta_df_list)

wmhv_padj <- wmhv %>%
  mutate(padj = p.adjust(p, "BH"))

## Save ----
rbind(bv_padj, hv_padj, wmhv_padj)  %>%
  write.csv("results/seroreac_main_MA.csv", row.names = F)

## End ----
