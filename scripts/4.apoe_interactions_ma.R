# ====================== #
#  4. APOE interactions  #  
#      meta-analyses     #
# ====================== #

## Set up ----
filepath <- ""
setwd(filepath)

## Load libraries ----
library(data.table)
library(dplyr)
library(magrittr)
library(metafor)
library(purrr)

## Read in data (APOE e2) ----
nshd <- read.csv("results/nshd_serostatus_apoe2.csv")
nshd_pbi <- read.csv("results/nshd_pbi_apoe2.csv")
nshd %<>%
  rbind(nshd_pbi) %>%
  mutate(study = "NSHD")

sabre_eur <- read.csv("results/sabre_serostatus_apoe2_eur.csv")
sabre_pbi_eur <- read.csv("results/sabre_pbi_apoe2_eur.csv")
sabre_eur %<>%
  rbind(sabre_pbi_eur) %>%
  mutate(study = "SABRE (EUR)")

sabre_sa <- read.csv("results/sabre_serostatus_apoe2_sa.csv")
sabre_pbi_sa <- read.csv("results/sabre_pbi_apoe2_sa.csv")
sabre_sa %<>%
  rbind(sabre_pbi_sa) %>%
  mutate(study = "SABRE (SA)")

ukb <- read.csv("results/ukb_serostatus_apoe2.csv")
ukb_pbi <- read.csv("results/ukb_pbi_apoe2.csv")
ukb %<>%
  rbind(ukb_pbi) %>%
  mutate(study = "UKB")

## Combine data ----
dat_apoe2 <- rbind(nshd, sabre_eur, sabre_sa, ukb) %>%
  mutate(pathogen = gsub(pathogen, pattern=":apoe2_carrier1|:apoe2_stat1", rep=""))

# Check
dim(dat_apoe2)
n_distinct(dat_apoe2$pathogen) 
table(dat_apoe2$pathogen) 
table(dat_apoe2$study) 
table(dat_apoe2$model) 
table(dat_apoe2$outcome) 
table(duplicated(dat_apoe2$estimate)) 
table(duplicated(dat_apoe2$std.error))
table(duplicated(dat_apoe2$p.value)) 

## Meta-analyses (APOE e2) ----
# Pathogens (& PBI)
pathogen <- dat_apoe2 %>%
  pull(pathogen) %>%
  unique()

length(pathogen)

models <- c("model 1", "model 2", "model 3")

### Brain volume ----
meta_df_list <- list()

for (i in pathogen) {
  for (j in models) {
    out <- "brain volume"
    print(paste0(out, ": ", i, ", model: ", j))
    res <-
      dat_apoe2 %>% filter(model == j &
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
      rma$k,
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
      "Q-p",
      "nstudy"
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
      dat_apoe2 %>% filter(model == j &
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
      rma$k,
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
      "Q-p",
      "nstudy"
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
      dat_apoe2 %>% filter(model == j &
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
      rma$k,
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
      "Q-p",
      "nstudy"
    )
    meta_df_list <- c(meta_df_list, list(meta_df))
  }
}

wmhv <- rbindlist(meta_df_list) %>%
  mutate(padj = p.adjust(p, "BH"))

## Save ----
rbind(bv, hv, wmhv) %>%
  write.csv("results/serostatus_apoe2_MA.csv", row.names = F)

# Clear environment
rm(list = ls())

## Read in data (APOE e4) ----
nshd <- read.csv("results/nshd_serostatus_apoe4.csv")
nshd_pbi <- read.csv("results/nshd_pbi_apoe4.csv")
nshd %<>%
  rbind(nshd_pbi) %>%
  mutate(study = "NSHD")

sabre_eur <- read.csv("results/sabre_serostatus_apoe4_eur.csv")
sabre_pbi_eur <- read.csv("results/sabre_pbi_apoe4_eur.csv")
sabre_eur %<>%
  rbind(sabre_pbi_eur) %>%
  mutate(study = "SABRE (EUR)")

sabre_sa <- read.csv("results/sabre_serostatus_apoe4_sa.csv")
sabre_pbi_sa <- read.csv("results/sabre_pbi_apoe4_sa.csv")
sabre_sa %<>%
  rbind(sabre_pbi_sa) %>%
  mutate(study = "SABRE (SA)")

ukb <- read.csv("results/ukb_serostatus_apoe4.csv")
ukb_pbi <- read.csv("results/ukb_pbi_apoe4.csv")
ukb %<>%
  rbind(ukb_pbi) %>%
  mutate(study = "UKB")

## Combine data ----
dat_apoe4 <- rbind(nshd, sabre_eur, sabre_sa, ukb) %>%
  mutate(pathogen = gsub(pathogen, pattern=":apoe4_carrier1|:apoe4_stat1", rep=""))

# Check
dim(dat_apoe4)
n_distinct(dat_apoe4$pathogen) 
table(dat_apoe4$pathogen) 
table(dat_apoe4$study) 
table(dat_apoe4$model) 
table(dat_apoe4$outcome) 
table(duplicated(dat_apoe4$estimate)) 
table(duplicated(dat_apoe4$std.error))
table(duplicated(dat_apoe4$p.value)) 

## Meta-analyses (APOE e4) ----
# Pathogens (& PBI)
pathogen <- dat_apoe4 %>%
  pull(pathogen) %>%
  unique()

length(pathogen) 

models <- c("model 1", "model 2", "model 3")

### Brain volume ----
meta_df_list <- list()

for (i in pathogen) {
  for (j in models) {
    out <- "brain volume"
    print(paste0(out, ": ", i, ", model: ", j))
    res <-
      dat_apoe4 %>% filter(model == j &
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
      rma$k,
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
      "Q-p",
      "nstudy"
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
      dat_apoe4 %>% filter(model == j &
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
      rma$k,
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
      "Q-p",
      "nstudy"
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
      dat_apoe4 %>% filter(model == j &
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
      rma$k,
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
      "Q-p",
      "nstudy"
    )
    meta_df_list <- c(meta_df_list, list(meta_df))
  }
}

wmhv <- rbindlist(meta_df_list) %>%
  mutate(padj = p.adjust(p, "BH"))

## Save ----
rbind(bv, hv, wmhv) %>%
  write.csv("results/serostatus_apoe4_MA.csv", row.names = F)

## End ----
