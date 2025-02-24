## ---------------------------
##
## Script name: 07_univariate_analyses.R
##
## Author: Victor Velasco-Pardo (vvp1@st-andrews.ac.uk)
##
## This script runs univariate analyses as described in the paper
##
library(tidyverse)
library(tidymodels)
library(parallel)
library(ggrepel)
setwd("/home/workspace/files/vvp1/paper_with_batch_effects/")

folder_names <- c("AMYLOID", "PAD_MCI", "PAD_CN")
ANALYSIS_AMYLOID <- 1
ANALYSIS_PAD_MCI <- 2
ANALYSIS_PAD_CN  <- 3

set.seed(1)

scale_cont <- function(x) (x - mean(x))/(2*sd(x))
scale_bin  <- function(x) x - mean(x)

for (THIS_ANALYSIS in c(ANALYSIS_PAD_MCI, ANALYSIS_PAD_CN)) {
  # Load data
  df_all <- read_csv("data/df_all.csv")
  df_meta <- read_csv("data/proteomics_meta.csv")
  
  # Data split
  df_all <- df_all %>% 
    mutate(
      age = AGE,
      sexMale = ifelse(SEX == "M", TRUE, FALSE), 
      yearsEducation = `Years-Education`,
      raceBlackOrAfricanAmerican = ifelse(`RACE` == "BLACK OR AFRICAN AMERICAN", TRUE, FALSE),
      raceAsian = ifelse(`RACE` == "ASIAN", TRUE, FALSE),
      raceUnknown = ifelse(`RACE` == "UNKNOWN", TRUE, FALSE),
      raceAmericanIndianOrAlaskan = ifelse(`RACE` == "AMERICAN INDIAN OR ALASKA NATIVE", TRUE, FALSE),
      raceNativeHawaiianOrIslander = ifelse(`RACE` == "NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER", TRUE, FALSE),
      ethnicHispanicOrLatino = ifelse(`ETHNIC` == "HISPANIC OR LATINO", TRUE, FALSE),
      ethnicNotReported = ifelse(`ETHNIC` == "NOT REPORTED", TRUE, FALSE),
      APOEe4 = apoe4 >= 1, APOEe4Homozygous = apoe4 >= 2,
      .keep = "unused"
    ) %>%
    select(-SUBJID, -QVAL, -`Subject-amyloid-classification`, -`VisQ-Amyloid-Classification`, -`Standard-Uptake-Value-Ratio`, -`Safety-QC-Grade`) %>% 
    mutate(across(seq.10000.28:seq.9999.1, log)) %>%
    drop_na()  %>% 
    relocate(Amyloid, age, sexMale, yearsEducation,
             starts_with("APOEe4"), starts_with("race"), starts_with("ethnic"),
             Vascular, Hypertension, Diabetes) %>%
    mutate_if(is.numeric, scale_cont) %>%
    mutate_if(is.logical, scale_bin) %>% 
    mutate(Amyloid = as.numeric(as.factor(Amyloid))-1)

  if (THIS_ANALYSIS == ANALYSIS_AMYLOID) {
    df_all <- df_all |> 
      select(-COHORT) |>
      mutate(outcome = Amyloid, .keep = "unused") |>
      relocate(outcome)
  } else if (THIS_ANALYSIS == ANALYSIS_PAD_MCI) {
    df_all <- df_all |> 
      select(-Amyloid) |>
      filter(COHORT != "Cohort 1 (Healthy)") |>
      mutate(outcome = ifelse(COHORT == "Cohort 3 (Probable AD)", 1, 0), .keep = "unused") |>
      relocate(outcome)
  } else {
    df_all <- df_all |> 
      select(-Amyloid) |>
      filter(COHORT != "Cohort 2 (MCI)") |>
      mutate(outcome = ifelse(COHORT == "Cohort 3 (Probable AD)", 1, 0), .keep = "unused") |>
      relocate(outcome)
  }
  colnames(df_all) |> head(20)
  str_c(paste0(colnames(df_all)[2:16], sep =  " + "), collapse = "")
  
  D <- sum(startsWith(colnames(df_all), "seq"))
  p_values <- rep(0, D)
  betas <- rep(0, D)
  names(p_values) <- colnames(df_all)[which(startsWith(colnames(df_all), "seq"))]
  names(betas)    <- colnames(df_all)[which(startsWith(colnames(df_all), "seq"))]
  for (d in 1:D) {
    if (d %% 100 == 0) 
      print(d)
    my_formula <- as.formula(
      str_c(
        str_c(
          "outcome ~ ", 
          str_c(colnames(df_all)[2:6], 
          collapse = " + ")
        ), 
        " + ",
        colnames(df_all)[16+d]
      )
    )
    my_fit <- glm(my_formula, family = binomial(), data = df_all)
    betas[d] <- coef(summary(my_fit))[7, 1]
    p_values[d] <- coef(summary(my_fit))[7, 4]
  }
  
  df_all |> 
    select(outcome, starts_with("seq")) |>
    group_by(outcome) |>
    summarise_all(mean) |>
    ungroup() |> 
    select(-outcome) |>
    as.matrix() |> exp() |> t() -> logFC
  logFC <- log2(logFC[,1]) - log2(logFC[,2])
  
  volcano_df <- as.data.frame(cbind(logFC, beta=betas, p.val = p_values))
  volcano_df$SeqId <- rownames(volcano_df)
  volcano_df <- volcano_df |> left_join(df_meta, by = join_by(SeqId)) 
  # volcano_df$adj.P.Val <- p_values*D
  volcano_df$adj.P.Val <- p.adjust(p_values, method = "bonferroni")
  padj <- 0.05
  p <- volcano_df |>
    mutate( effect = ifelse( adj.P.Val < 0.05 & logFC > 0, "POSITIVE", 
                                 ifelse( adj.P.Val < 0.05 & logFC < 0, "NEGATIVE", "NS")) ) |>
    ggplot( aes( x=logFC, y=-log10(p.val), color=effect) ) +
    geom_point() +
    scale_color_manual(breaks = c("POSITIVE", "NEGATIVE", "NS"),
                       values=c("red", "blue", "gray")) +
    xlab( expression(log[2](fold-change)) ) +
    ylab( expression(-log[10](p.val)) ) +
    geom_hline(yintercept = -log10(0.05/D), colour="black", linetype = "dotted") +
    geom_vline(xintercept = 0, colour="black", linetype = "dotted") +
  #  facet_wrap(~group) +
    theme_classic() +
    theme(strip.text = element_text(size = ifelse(THIS_ANALYSIS==ANALYSIS_AMYLOID, 18, 14)),
          axis.text=element_text(size=ifelse(THIS_ANALYSIS==ANALYSIS_AMYLOID, 18, 14)),
          axis.title=element_text(size=ifelse(THIS_ANALYSIS==ANALYSIS_AMYLOID, 18, 14)),face="bold") +
    geom_text_repel(aes(label = ifelse(adj.P.Val < padj, Target, "")), , hjust = -0.1, vjust = 0.5, size = ifelse(THIS_ANALYSIS==ANALYSIS_AMYLOID, 6, 4), color = "black")
  ggsave(paste0("figures/univariate_", folder_names[THIS_ANALYSIS], ".jpeg"), 
         p,
         dpi = 300, 
         units = "px", 
         width = ifelse(THIS_ANALYSIS == ANALYSIS_AMYLOID, 2880, 1440), 
         height = 1614 #ifelse(THIS_ANALYSIS == ANALYSIS_AMYLOID, 1614, 807)
  )
  volcano_df |>
    relocate(SeqId, SomaId, TargetFullName, Target, UniProt, logFC, beta, p.val, adj.P.Val) |>
    arrange(adj.P.Val) |>
    write_csv(paste0("supplementary/results_univariate_", folder_names[THIS_ANALYSIS], ".csv"))
}

proteomics_meta <- read_csv("data/proteomics_meta.csv")
relevant_proteins <- proteomics_meta |> 
  filter(SeqId %in% names(which(p_values < 0.05/D)))

volcano_df |> 
  relocate(SeqId, SomaId, TargetFullName, Target, UniProt, beta, p.val, adj.P.Val, logFC) |>
  arrange(adj.P.Val) |>
  tail()

