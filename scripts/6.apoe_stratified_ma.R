# ===================== #
#  6. APOE2 stratified  #  
#     meta-analyses     #
# ===================== #

## Set up ----
filepath <- ""
setwd(filepath)

## Load libraries ----
library(data.table)
library(dplyr)
library(magrittr)
library(metafor)
library(purrr)

## Read in data ----
nshd <- read.csv("results/nshd_e2_strat.csv")
sabre_eur <- read.csv("results/sabre_e2_strat_eur.csv")
sabre_sa <- read.csv("results/sabre_e2_strat_sa.csv")
ukb <- read.csv("results/ukb_e2_strat.csv")

## Combine data ----
dat <- rbind(nshd, sabre_eur, sabre_sa, ukb)

# Check
dim(dat)
n_distinct(dat$pathogen) 
table(dat$pathogen) 
table(dat$study) 
table(dat$model) 
table(dat$outcome) 
table(dat$analysis) 
table(duplicated(dat$estimate))
table(duplicated(dat$std.error)) 
table(duplicated(dat$p.value)) 

## Meta-analyses ----
analysis <- unique(dat$analysis)
pathogen <- unique(dat$pathogen)
models <- c("model 1", "model 2", "model 3")

### APOE e2 carriers ----
meta_df_list <- list()

for (i in pathogen){
  for (j in models){
    res <- dat %>% filter(pathogen == i) %>% 
      filter(model == j) %>% 
      filter(analysis == "e2 carrier")
    print(paste("Number of rows of filtered data:", nrow(res)))
    print(paste("Number of studies:", n_distinct(res$study)))
    print(res)
    rma <- rma(yi=estimate, vi=std.error^2, method="REML", data=res)
    meta_df <- data.frame(i, j, rma$beta, rma$zval,
                          rma$pval, rma$ci.lb, rma$ci.ub,
                          rma$se, rma$I2, rma$QE, rma$QEp, rma$k,
                          row.names = NULL)
    names(meta_df) <- c("pathogen", "model", "beta",
                        "z", "p", "ci_lower", "ci_upper", "se",
                        "i2", "Q", "Q-p", "n_study")
    meta_df_list <- c(meta_df_list, list(meta_df))
  }
}

e2_carrier <- rbindlist(meta_df_list) %>%
  mutate(analysis = "e2 carrier")

### APOE e2 non-carriers ----
meta_df_list <- list()

for (i in pathogen){
  for (j in models){
    res <- dat %>% filter(pathogen == i) %>% 
      filter(model == j) %>% 
      filter(analysis == "non-e2 carrier")
    print(paste("Number of rows of filtered data:", nrow(res)))
    print(paste("Number of studies:", n_distinct(res$study)))
    print(res)
    rma <- rma(yi=estimate, vi=std.error^2, method="REML", data=res)
    meta_df <- data.frame(i, j, rma$beta, rma$zval,
                          rma$pval, rma$ci.lb, rma$ci.ub,
                          rma$se, rma$I2, rma$QE, rma$QEp, rma$k,
                          row.names = NULL)
    names(meta_df) <- c("pathogen", "model", "beta",
                        "z", "p", "ci_lower", "ci_upper", "se",
                        "i2", "Q", "Q-p", "n_study")
    meta_df_list <- c(meta_df_list, list(meta_df))
  }
}

non_e2_carrier <- rbindlist(meta_df_list) %>%
  mutate(analysis = "non-e2 carrier")

## Save ----
write.csv(rbind(e2_carrier, non_e2_carrier), 
          "results/apoe2_strat_ma.csv", row.names=F)

## End ----
