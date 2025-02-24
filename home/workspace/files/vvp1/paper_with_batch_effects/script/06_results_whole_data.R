## ---------------------------
##
## Script name: 06_results_whole_data.R
##
## Author: Victor Velasco-Pardo (vvp1@st-andrews.ac.uk)
##
## This script produces the plots summarising the results of the multivariate analyses
#
library(tidyverse)
library(parallel)
library(hsstan)
library(bayesplot)
library(rstan)
library(gridExtra)
library(pROC)
library(ggrepel)
options(mc.cores = detectCores())
setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")

set.seed(1)

folder_names <- c("AMYLOID", "PAD_MCI", "PAD_CN")
ANALYSIS_AMYLOID <- 1
ANALYSIS_PAD_MCI <- 2
ANALYSIS_PAD_CN  <- 3

for (THIS_ANALYSIS in c(ANALYSIS_AMYLOID, ANALYSIS_PAD_MCI, ANALYSIS_PAD_CN)) {
  setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")
  setwd(folder_names[THIS_ANALYSIS])
  
  fits <- readRDS("bin/fit_baseline.rds")
  fits_full <- readRDS("bin/fit_full.rds")
  
  # Load SomaScan metadata
  proteomics_meta <- read_csv("../data/proteomics_meta.csv")
  
  # Most relevant covariates, and most relevant proteins
  summary(fits_full) -> my_summary
  my_summary |>
    as.data.frame() |> 
    round(2) |> 
    cbind(data.frame(parameter = rownames(my_summary)))  |> 
    head(-2) |> tail(-1) |> 
    as_tibble() |> 
    arrange(desc(abs(mean))) |>
    left_join(proteomics_meta, by = c("parameter" = "SeqId")) |>
    relocate(where(is.character)) -> my_summary
  
  # Write csv file with results
  write_csv(my_summary, paste0("../supplementary/results_multivariate_", folder_names[THIS_ANALYSIS], ".csv"))
  
  my_summary |> 
    filter(str_detect(parameter, "^seq")) |> 
    head(25) |> pull(parameter) -> most_relevant_proteins
  
  my_summary |> 
    filter(str_detect(parameter, "^seq")) |> 
    head(25) |> pull(Target) -> most_relevant_proteins_target
  
  my_summary |> 
    filter(str_detect(parameter, "^seq")) |> 
    head(25) |> pull(TargetFullName) -> most_relevant_proteins_target_fn
  
  fit1.mat <- fits_full$stanfit |> as.matrix()
  fit1.mat <- fit1.mat[, most_relevant_proteins]
  colnames(fit1.mat) <- most_relevant_proteins_target
  
  color_scheme_set(scheme = "red")
  
  mcmc_intervals(fit1.mat, point_est = "mean", prob = 0.90, prob_outer = 0.95) + 
    xaxis_text(size = 16) + yaxis_text(size = 16) + 
    scale_x_continuous(breaks = seq(-3, 4, by = 0.5)) 
    ggsave(
    paste0("../figures/multivariate_results_", folder_names[THIS_ANALYSIS],".jpeg"), 
    dpi = 300, 
    units = "px", 
    width =  2880, #ifelse(THIS_ANALYSIS == ANALYSIS_AMYLOID, 2880, 1440), 
    height = 1614 # ifelse(THIS_ANALYSIS == ANALYSIS_AMYLOID, 1614, 807)
  )
}

my_summary |>
  filter(parameter != "raceOtherUnknown") |>
#  head() |>
  ggplot(aes(x = parameter, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`)) +
  theme(axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank()) + 
  geom_text_repel(aes(label = ifelse(abs(mean)>0.1, ifelse(is.na(Target), parameter, Target), "")), , hjust = -0.1, vjust = 0.5, size = 4, color = "black")



