# Prepare data statistics for supplementary table S1 and S2.
# Created 7/3-2025

rm(list=ls())
library(dplyr)
library(tidyr)
load("2502-indv-level-code/traitDataFujian-Ind-step1.RData") # from step 1

# For manual confirmation
# write.csv(traitDataIndv_SelectedTraits %>% dplyr::select(-SpeciesCH), file="2502-indv-level-code/traitDataIndv-statistics250307.csv", row.names=FALSE)

# Sum how many species and individuals are sampled for each trait
traitDataIndv_long = traitDataIndv_SelectedTraits %>% 
  dplyr::select(SampleID, SpeciesFullName, GrowthForm,
                RTD,SRL,RD,RNC,SRA,RPC,RCC,RDMC,SRR25,
                LMA,LNC,LCC,LPC,Asat,Vcmax25,Rdark25P) %>%
  tidyr::gather(key = "trait", value = "value", -SampleID, -SpeciesFullName, -GrowthForm) %>%
  dplyr::filter(!is.na(value))


traitDataIndv_spindvsum = traitDataIndv_long %>% 
  group_by(trait) %>% 
  summarise(nspecies = n_distinct(SpeciesFullName), 
            nindv = n_distinct(SampleID),
            traitvaluemean = mean(value, na.rm=TRUE),
            traitvalueSD = sd(value, na.rm=TRUE))

traitDataIndv_spGFsum = traitDataIndv_long %>% 
  group_by(trait, GrowthForm) %>% 
  summarise(nspecies = n_distinct(SpeciesFullName), 
            nindv = n_distinct(SampleID),
            traitvaluemean = mean(value, na.rm=TRUE),
            traitvalueSD = sd(value, na.rm=TRUE))

# TODO: ADDED 13/3-25: examine LMA and LNC difference between trees and shrubs
traitSubData_LMA = traitDataIndv_long %>% dplyr::filter(trait %in% c("LMA")) %>%
  dplyr::filter(GrowthForm %in% c("tree", "shrub"))

traitSet_LMA_tree = as.numeric(traitSubData_LMA %>% dplyr::filter(GrowthForm == "tree") %>% dplyr::pull(value))
traitSet_LMA_shrub = as.numeric(traitSubData_LMA %>% dplyr::filter(GrowthForm == "shrub") %>% dplyr::pull(value))
wilcox_paired_result_LMA = wilcox.test(traitSet_LMA_tree, traitSet_LMA_shrub, paired = F)


traitSubData_LNC = traitDataIndv_long %>% dplyr::filter(trait %in% c("LNC")) %>%
  dplyr::filter(GrowthForm %in% c("tree", "shrub"))

traitSet_LNC_tree = as.numeric(traitSubData_LNC %>% dplyr::filter(GrowthForm == "tree") %>% dplyr::pull(value))
traitSet_LNC_shrub = as.numeric(traitSubData_LNC %>% dplyr::filter(GrowthForm == "shrub") %>% dplyr::pull(value))
wilcox_paired_result_LNC = wilcox.test(traitSet_LNC_tree, traitSet_LNC_shrub, paired = F)

