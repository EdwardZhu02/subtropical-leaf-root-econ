rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggsignif)  # For adding significance markers
library(broom)     # For tidying the t-test results
library(ggpubr)
library(patchwork) # for merging plots
library(plotly)
library(ggsignif)  # For adding significance markers
library(rstatix)   # For pairwise wilcox test

# Perform across-site comparison using community-mean trait values
load("traitDataFujian-spavg-phylo-step2.RData")

# use data frame 'traitDataIndv', the original, not species averaged or normalized data for analysis
traitDataIndv_touse = traitDataIndv # backup original df, while copy it for this script

plasCalc_spIDCounts = as.data.frame(table(traitDataIndv_touse$SpeciesID)) %>%
  dplyr::rename(SpeciesID = Var1, Count = Freq) %>%
  dplyr::filter(Count > 2)

plasCalc_rawdata = traitDataIndv_touse %>%
  dplyr::filter(SpeciesID %in% plasCalc_spIDCounts$SpeciesID) %>%
  dplyr::select(where(is.numeric), GrowthForm, SampleID) %>%   # Select numeric columns and SiteLoc
  pivot_longer(cols = -c(GrowthForm, SampleID), names_to = "trait", values_to = "value") %>% 
  group_by(trait) %>%
  filter(!is.na(value))