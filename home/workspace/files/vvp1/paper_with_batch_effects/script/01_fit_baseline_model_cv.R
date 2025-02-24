## ---------------------------
##
## Script name: 01_fit_baseline_model_cv.R
##
## Author: Victor Velasco-Pardo (vvp1@st-andrews.ac.uk)
##
## This script runs the baseline model (using demographic covariates
## only) in each of the 10 cross validation folds
## Before running the script, change the path to your local directory
## Run the script from the command line using 
##      Rscript script/01_fit_baseline_model_cv.R <cv-fold>
## replacing <cv-fold> with the index of the CV fold to be used (0 to 9)
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
v <- as.integer(args[1]) # Fold in cross-validation
set.seed(1)

setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")
THIS_ANALYSIS <- ANALYSIS_AMYLOID
# THIS_ANALYSIS <- ANALYSIS_PAD_MCI
# THIS_ANALYSIS <- ANALYSIS_PAD_CN
setwd(folder_names[THIS_ANALYSIS])

data_train <- readRDS(paste0("bin/", v, "/data_train.rds"))

fit1 <- hsstan(
  data_train |> select(-site) |> mutate(Amyloid = as.integer(as.factor(Amyloid)) - 1),
  as.formula(paste("Amyloid ~ ", paste(names(data_train)[3:14], collapse = " + "))),
  penalized = NULL,
  family = binomial(),
  seed = 123,
  chains = 2, 
  cores = 2,
  scale.u = 1,
  qr = FALSE
)

saveRDS(fit1, paste0("bin/", v, "/fit_baseline.rds"))

setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")
THIS_ANALYSIS <- ANALYSIS_PAD_MCI
setwd(folder_names[THIS_ANALYSIS])

data_train <- readRDS(paste0("bin/", v, "/data_train.rds"))

fit1 <- hsstan(
  data_train |> select(-site),
  as.formula(paste("COHORT ~ ", paste(names(data_train)[3:14], collapse = " + "))),
  penalized = NULL,
  family = binomial(),
  seed = 123,
  chains = 2, 
  cores = 2,
  scale.u = 1,
  qr = FALSE
)

saveRDS(fit1, paste0("bin/", v, "/fit_baseline.rds"))

setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")
THIS_ANALYSIS <- ANALYSIS_PAD_CN
setwd(folder_names[THIS_ANALYSIS])

data_train <- readRDS(paste0("bin/", v, "/data_train.rds"))

fit1 <- hsstan(
  data_train |> select(-site),
  as.formula(paste("COHORT ~ ", paste(names(data_train)[3:14], collapse = " + "))),
  penalized = NULL,
  family = binomial(),
  seed = 123,
  chains = 2, 
  cores = 2,
  scale.u = 1,
  qr = FALSE
)

saveRDS(fit1, paste0("bin/", v, "/fit_baseline.rds"))
