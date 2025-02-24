## ---------------------------
##
## Script name: 04_fit_full_model_whole_data.R
##
## Author: Victor Velasco-Pardo (vvp1@st-andrews.ac.uk)
##
## This script runs the full model (with both demographic and proteomic
## covariates) the complete dataset  
## Before running the script, change the path to your local directory
##
library(tidyverse)
library(parallel)
library(hsstan)
options(mc.cores = detectCores())
setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")

folder_names <- c("AMYLOID", "PAD_MCI", "PAD_CN")
ANALYSIS_AMYLOID <- 1
ANALYSIS_PAD_MCI <- 2
ANALYSIS_PAD_CN  <- 3

args <- commandArgs(trailingOnly = TRUE)
THIS_ANALYSIS <- as.integer(args[1]) # Fold in cross-validation

set.seed(1)

# Load data
df_all <- read_csv("data/df_all.csv")
df_meta <- read_csv("data/proteomics_meta.csv")

scale_cont <- function(x) (x - mean(x))/(2*sd(x))
scale_bin  <- function(x) x - mean(x)

df_all <- df_all |> 
  mutate(
    age = AGE,
    sexMale = ifelse(SEX == "M", TRUE, FALSE), 
    yearsEducation = `Years-Education`,
    raceBlackOrAfricanAmerican = ifelse(`RACE` == "BLACK OR AFRICAN AMERICAN", TRUE, FALSE),
    raceAsian = ifelse(`RACE` == "ASIAN", TRUE, FALSE),
    raceOtherUnknown = ifelse(`RACE` %in% c("UNKOWN", "AMERICAN INDIAN OR ALASKA NATIVE", "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER"), TRUE, FALSE),
    ethnicHispanicOrLatino = ifelse(`ETHNIC` == "HISPANIC OR LATINO", TRUE, FALSE),
    ethnicNotReported = ifelse(`ETHNIC` == "NOT REPORTED", TRUE, FALSE),
    APOEe4 = apoe4 >= 1, 
    site = substring(SUBJID, 1, 5),
    .keep = "unused"
  ) |>
  select(-QVAL, -`Subject-amyloid-classification`, -`VisQ-Amyloid-Classification`, -`Standard-Uptake-Value-Ratio`, -`Safety-QC-Grade`) |> 
  mutate(across(seq.10000.28:seq.9999.1, log)) |>
  relocate(Amyloid, COHORT, age, sexMale, yearsEducation,
           starts_with("APOEe4"), starts_with("race"), starts_with("ethnic"),
           Vascular, Hypertension, Diabetes) |>
  drop_na()

df_all <- df_all |> 
  mutate_if(is.numeric, scale_cont) |>
  mutate_if(is.logical, scale_bin) |> 
  mutate(Amyloid = as.numeric(as.factor(Amyloid))-1)
for (k in 15:7349) {
  df_all[, k] <- unname(residuals(lm(as.formula(paste(names(df_all)[k], "~ site")), data = df_all)))
}

setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")
setwd(folder_names[THIS_ANALYSIS])

if (THIS_ANALYSIS == ANALYSIS_AMYLOID) {
  fit1 <- hsstan(
    df_all |> select(-site),
    as.formula(paste("Amyloid ~ ", paste(names(df_all)[3:14], collapse = " + "))),
    penalized = names(df_all)[15:(ncol(df_all)-1)],
    family = binomial(),
    iter = 6000,
    warmup = 1000,
    seed = 123,
    chains = 2, 
    cores = 2,
    scale.u = 1,
    qr = FALSE
  )
} else if (THIS_ANALYSIS == ANALYSIS_PAD_MCI) {
  df_all <- df_all |>
    select(-site) |>
    filter(COHORT != "Cohort 1 (Healthy)") |>
    mutate(COHORT = as.integer(as.factor(COHORT)) - 1)
  fit1 <- hsstan(
    df_all,
    as.formula(paste("COHORT ~ ", paste(names(df_all)[3:14], collapse = " + "))),
    penalized = names(df_all)[15:(ncol(df_all))],
    family = binomial(),
    iter = 6000,
    warmup = 1000,
    seed = 123,
    chains = 2, 
    cores = 2,
    scale.u = 1,
    qr = FALSE
  )
  
} else if (THIS_ANALYSIS == ANALYSIS_PAD_CN) {
  df_all <- df_all |>
    select(-site) |>
    filter(COHORT != "Cohort 2 (MCI)") |>
    mutate(COHORT = as.integer(as.factor(COHORT)) - 1)
  fit1 <- hsstan(
    df_all,
    as.formula(paste("COHORT ~ ", paste(names(df_all)[3:14], collapse = " + "))),
    penalized = names(df_all)[15:(ncol(df_all))],
    family = binomial(),
    iter = 6000,
    warmup = 1000,
    seed = 123,
    chains = 2, 
    cores = 2,
    scale.u = 1,
    qr = FALSE
  )
}
saveRDS(fit1, "bin/fit_full.rds")
