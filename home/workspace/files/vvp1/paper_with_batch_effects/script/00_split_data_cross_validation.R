## ---------------------------
##
## Script name: 00_split_data_cross_validation.R
##
## Author: Victor Velasco-Pardo (vvp1@st-andrews.ac.uk)
##
## This script creates 10 cross-validation folds from the data 
## from the data frame of demographic covariates and proteomics
## stored in df_all.csv
#
library(tidyverse)
library(tidymodels)
setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")

# Cross validation folds for three analysis are created.
# CV folds and results from downstream analyses are stored
# in directories listed in "folder_names"
folder_names <- c("AMYLOID", "PAD_MCI", "PAD_CN")
ANALYSIS_AMYLOID <- 1
ANALYSIS_PAD_MCI <- 2
ANALYSIS_PAD_CN  <- 3

set.seed(1)
V <- 10
for (THIS_ANALYSIS in 1:length(folder_names)) {
  # Load data
  df_all <- read_csv("data/df_all.csv")
  df_meta <- read_csv("data/proteomics_meta.csv")
  
  # Data split
  df_all <- df_all |> 
    #   select(-RACE, -ETHNIC) |>
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
  
  if (THIS_ANALYSIS == ANALYSIS_AMYLOID) {
    data_folds <- vfold_cv(df_all, v = V, strata = Amyloid)
  } else if (THIS_ANALYSIS == ANALYSIS_PAD_MCI) {
    df_all <- df_all |> filter(COHORT != "Cohort 1 (Healthy)")
    data_folds <- vfold_cv(df_all, v = V, strata = COHORT)
  } else {
    df_all <- df_all |> filter(COHORT != "Cohort 2 (MCI)")
    data_folds <- vfold_cv(df_all, v = V, strata = COHORT)
  }
  
  setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")
  setwd(folder_names[THIS_ANALYSIS])

  cv_indices <- list()
  cv_train_df <- list()
  cv_test_df <- list()
  
  for (v in 1:V) {
    cv_indices[[v]] <- data_folds[v,1][[1]][[1]]$in_id
    print(v)
  }
  for (v in 1:V) {
    my_df <- df_all[cv_indices[[v]], ] |> 
      mutate(Amyloid = as.numeric(as.factor(Amyloid))-1) |>
      mutate(COHORT = as.numeric(as.factor(COHORT))-1)
    for (k in 15:7349) {
      my_df[, k] <- unname(residuals(lm(as.formula(paste(names(df_all)[k], "~ site")), data = my_df)))
    }
      
    saveRDS(my_df, paste0("bin/", v-1, "/data_train.rds"))
    print(v)
  }
  for (v in 1:V) {
    my_df <- df_all[-cv_indices[[v]], ] |> 
      mutate(Amyloid = as.numeric(as.factor(Amyloid))-1) |>
      mutate(COHORT = as.numeric(as.factor(COHORT))-1)
    for (k in 15:7349) {
      my_df[, k] <- unname(residuals(lm(as.formula(paste(names(df_all)[k], "~ site")), data = my_df)))
    }
    
    saveRDS(my_df, paste0("bin/", v-1, "/data_test.rds"))
    print(v)
  }
}
