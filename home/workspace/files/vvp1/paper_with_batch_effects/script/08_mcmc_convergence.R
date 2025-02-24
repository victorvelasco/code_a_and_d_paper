## ---------------------------
##
## Script name: 08_mcmc_convergence.R
##
## Author: Victor Velasco-Pardo (vvp1@st-andrews.ac.uk)
##
## This script produces plots to check convergence of the MCMC algorithm
##
library(tidyverse)
library(parallel)
library(hsstan)
library(bayesplot)
library(rstan)
library(gridExtra)
library(pROC)
options(mc.cores = detectCores())
setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")

set.seed(1)

fits_bs <- list()
fits_full <- list()
analyses <- c("AMYLOID", "PAD_MCI", "PAD_CN")
analyses2 <- c("PET Amyloid", "Probable AD vs Mild Cognitive Impairment", "Probable AD vs Cognitively normal controls")
for (i in 1:3) {
  fits_bs[[i]] <- readRDS(paste0(analyses[i], "/bin/fit_baseline.rds"))
  fits_full[[i]] <- readRDS(paste0(analyses[i], "/bin/fit_full.rds"))
}

color_scheme_set("red")
plots <- list()
for (i in 1:3) {
  plots[[2*i-1]] <- mcmc_trace(fits_bs[[i]]$stanfit, pars = "lp__") + 
    xaxis_text(size = 18) + yaxis_text(size = 18) +
    ggtitle(paste0(analyses2[i], " (Baseline)"))
  plots[[2*i]]   <- mcmc_trace(fits_full[[i]]$stanfit, pars = "lp__") + 
    xaxis_text(size = 18) + yaxis_text(size = 18) +
    ggtitle(paste0(analyses2[i], " (Full model)"))
}
do.call("grid.arrange", c(plots, ncol=2))
