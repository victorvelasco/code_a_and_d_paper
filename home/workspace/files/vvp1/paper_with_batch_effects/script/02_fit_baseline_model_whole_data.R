## ---------------------------
##
## Script name: 02_fit_baseline_model_whole_data.R
##
## Author: Victor Velasco-Pardo (vvp1@st-andrews.ac.uk)
##
## This script runs the baseline model (using demographic covariates
## only) on the complete dataset 
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

set.seed(1)

# Load data
df_all <- read_csv("data/df_all.csv")
df_meta <- read_csv("data/proteomics_meta.csv")

scale_cont <- function(x) (x - mean(x))/(2*sd(x))
scale_bin  <- function(x) x - mean(x)

data_train <- df_all |> 
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

data_train <- data_train |> 
  mutate_if(is.numeric, scale_cont) |>
  mutate_if(is.logical, scale_bin) 

setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")
THIS_ANALYSIS <- ANALYSIS_AMYLOID
setwd(folder_names[THIS_ANALYSIS])

fit1 <- hsstan(
  data_train |> select(-site) |> mutate(Amyloid = as.integer(as.factor(Amyloid)) - 1),
  as.formula(paste("Amyloid ~ ", paste(names(data_train)[3:14], collapse = " + "))),
  penalized = NULL,
  family = binomial(),
  iter = 6000,
  warmup = 1000,
  seed = 123,
  chains = 2, 
  cores = 2,
  scale.u = 1,
  qr = FALSE
)

saveRDS(fit1, "bin/fit_baseline.rds")

setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")
THIS_ANALYSIS <- ANALYSIS_PAD_MCI
setwd(folder_names[THIS_ANALYSIS])

data_train_pad_mci <- data_train |>
  select(-site) |>
  filter(COHORT != "Cohort 1 (Healthy)") |>
  mutate(COHORT = as.integer(as.factor(COHORT)) - 1)
fit1 <- hsstan(
  data_train_pad_mci,
  as.formula(paste("COHORT ~ ", paste(names(data_train)[3:14], collapse = " + "))),
  penalized = NULL,
  family = binomial(),
  iter = 6000,
  warmup = 1000,
  seed = 123,
  chains = 2, 
  cores = 2,
  scale.u = 1,
  qr = FALSE
)

saveRDS(fit1, "bin/fit_baseline.rds")


setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")
THIS_ANALYSIS <- ANALYSIS_PAD_CN
setwd(folder_names[THIS_ANALYSIS])
data_train_pad_cn <- data_train |>
  select(-site) |>
  filter(COHORT != "Cohort 2 (MCI)") |>
  mutate(COHORT = as.integer(as.factor(COHORT)) - 1)
fit1 <- hsstan(
  data_train_pad_cn,
  as.formula(paste("COHORT ~ ", paste(names(data_train)[3:14], collapse = " + "))),
  penalized = NULL,
  family = binomial(),
  iter = 6000,
  warmup = 1000,
  seed = 123,
  chains = 2, 
  cores = 2,
  scale.u = 1,
  qr = FALSE
)

saveRDS(fit1, "bin/fit_baseline.rds")
