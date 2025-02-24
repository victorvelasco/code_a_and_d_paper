## ---------------------------
##
## Script name: 05_model_performance.R
##
## Author: Victor Velasco-Pardo (vvp1@st-andrews.ac.uk)
##
## This script produces the ROC curves and calculates areas under
## the ROC curves (AUC).
##
library(tidyverse)
library(parallel)
library(hsstan)
library(bayesplot)
library(gridExtra)
library(caret)
library(pROC)
library(rstan)
options(mc.cores = detectCores())
setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")

folder_names <- c("AMYLOID", "PAD_MCI", "PAD_CN")
ANALYSIS_AMYLOID <- 1
ANALYSIS_PAD_MCI <- 2
ANALYSIS_PAD_CN  <- 3

# THIS_ANALYSIS <- ANALYSIS_AMYLOID
# THIS_ANALYSIS <- ANALYSIS_PAD_MCI
# THIS_ANALYSIS <- ANALYSIS_PAD_CN

set.seed(1)

df_meta <- read_csv("data/proteomics_meta.csv")
V <- 10
for (THIS_ANALYSIS in c(ANALYSIS_AMYLOID, ANALYSIS_PAD_MCI, ANALYSIS_PAD_CN)) {
  fits <- list()
  fits_baseline <- list()
  setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")
  setwd(folder_names[THIS_ANALYSIS])
  data_test <- list()
  data_train <- list()
  for (v in 1:V) {
    fits_baseline[[v]] <- readRDS(paste0("bin/", v-1, "/fit_baseline.rds"))
    fits[[v]] <- readRDS(paste0("bin/", v-1, "/fit_full.rds"))
    data_train[[v]] <- readRDS(paste0("bin/", v-1, "/data_train.rds")) 
    data_test[[v]] <- readRDS(paste0("bin/", v-1, "/data_test.rds")) 
  }
  
  
  train.pred.bs<- list(); test.pred.bs <- list()
  train.pred <- list(); test.pred <- list()
  train.y <- list(); test.y <- list()
  for (v in 1:V) {
    fit_bs1 <- fits_baseline[[v]]
    fit1 <- fits[[v]]
    print(v)
    
    # Predictions on training set
    train.preds <- posterior_linpred(fit1)
    train.pred[[v]] <- colMeans(train.preds)
    train.y[[v]] <- data_train[[v]] |> pull(Amyloid)
    train.preds.bs <- posterior_linpred(fit_bs1)
    train.pred.bs[[v]] <- colMeans(train.preds.bs)
    
    # Predictions on testing set
    test.preds <- posterior_linpred(fit1, newdata = data_test[[v]])
    test.pred[[v]] <- colMeans(test.preds)
    test.y[[v]] <- data_test[[v]] |> pull(Amyloid)
    test.preds.bs <- posterior_linpred(fit_bs1, newdata = data_test[[v]])
    test.pred.bs[[v]] <- colMeans(test.preds.bs)
    
  }
  p<- grid.arrange(
    ggroc(
      c(
        lapply(1:8, function(v) roc(train.y[[v]], train.pred[[v]])),
        lapply(1:8, function(v) roc(train.y[[v]], train.pred.bs[[v]])),
        list(roc(unlist(train.y), unlist(train.pred))),
        list(roc(unlist(train.y), unlist(train.pred.bs)))
      ), legacy.axes = TRUE
    ) + scale_colour_manual(values = c(rep("#9dcdef",8), rep("#e7a9b0", 8), "#0000FF", "#FF0000")) +
      theme_bw() + theme(legend.position = "none") + 
      theme(strip.text = element_text(size = 20),
        axis.text=element_text(size=20),
        axis.title=element_text(size=20)) +
      annotate("text", x = 0.75, y = 0.35, 
               label = paste0("AUC = ", round(as.numeric(auc(roc(unlist(train.y), unlist(train.pred))), "auc"), 3)), 
               colour = "#0000FF", size = 8) +
      annotate("text", x = 0.75, y = 0.25, 
               label = paste0("AUC = ", round(as.numeric(auc(roc(unlist(train.y), unlist(train.pred.bs))), "auc"), 3)), 
               colour = "#FF0000", size = 8) +
      ggtitle("Training sets") + theme(plot.title = element_text(size=18)),
    
    ggroc(
      c(
        lapply(1:8, function(v) roc(test.y[[v]], test.pred[[v]])),
        lapply(1:8, function(v) roc(test.y[[v]], test.pred.bs[[v]])),
        list(roc(unlist(test.y), unlist(test.pred))),
        list(roc(unlist(test.y), unlist(test.pred.bs)))
      ), legacy.axes = TRUE
    ) + scale_colour_manual(values = c(rep("#9dcdef",8), rep("#e7a9b0", 8), "#0000FF", "#FF0000")) +
      theme_bw() + theme(legend.position = "none") + 
      theme(strip.text = element_text(size = 20),
            axis.text=element_text(size=20),
            axis.title=element_text(size=20)) +
      annotate("text", x = 0.75, y = 0.35, 
               label = paste0("AUC = ", round(as.numeric(auc(roc(unlist(test.y), unlist(test.pred))), "auc"), 3)), 
               colour = "#0000FF", size = 8) +
      annotate("text", x = 0.75, y = 0.25, 
               label = paste0("AUC = ", round(as.numeric(auc(roc(unlist(test.y), unlist(test.pred.bs))), "auc"), 3)), 
               colour = "#FF0000", size = 8) +
      ggtitle("Test sets") + theme(plot.title = element_text(size = 18)),
    ncol = 2
  )
  ggsave(paste0("../figures/roc_plot_", folder_names[THIS_ANALYSIS], ".jpeg"), 
         p,
         dpi = 300, 
         units = "px", 
         width =  2880, #ifelse(THIS_ANALYSIS == ANALYSIS_AMYLOID, 2880, 1440), 
         height = 1614 #ifelse(THIS_ANALYSIS == ANALYSIS_AMYLOID, 1614, 807)
  )
}

