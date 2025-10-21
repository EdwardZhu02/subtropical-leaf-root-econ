rm(list=ls())
library(tidyverse)
library(dplyr)
library(tidyr)
# devtools::install_github("jinyizju/V.PhyloMaker2")
library(V.PhyloMaker2) # phylo-tree construction

load(file = "traitDataFujian-normbysp-step1.RData")

# Use plantlist package to normalize species names
# Perform round 1 checking
# ------------------------------------------------------------------------------
# devtools::install_github("helixcn/plantlist", build_vignettes = TRUE)
# library(plantlist) # species check
# spTaxaNames = traitDataIndv_normalized %>% dplyr::select(SpeciesFullName) %>% distinct() %>%
#   dplyr::rename(speciesFullName = SpeciesFullName) %>%
#   dplyr::mutate(
#     speciesCheckResult1 = "",
#     speciesFullNameAutoFix1 = "",  # for 'unresolved' or 'synonym' species, replace to accepted species
#   )
#
# for (rowidx in 1:nrow(spTaxaNames)) {
#   indv_spname = spTaxaNames$speciesFullName[rowidx]
#   spname_check_result = status(indv_spname, exact=T)
# 
#   if ( !is.na(spname_check_result[1,3]) ){
# 
#     if ( spname_check_result[1,3] == "Accepted" ) {
#       spTaxaNames$speciesCheckResult1[rowidx] = "PASS-accepted"
# 
#     } else if ( spname_check_result[1,3] %in% c("Unresolved", "Synonym") ) {
#       spTaxaNames$speciesFullNameAutoFix1[rowidx] = spname_check_result[1,5]
#       spTaxaNames$speciesCheckResult1[rowidx] = "PASS-autofixed"
#       print(paste0("Autofixed: ", indv_spname, " -> ", spname_check_result[1,5]))
# 
#     } else {
#       spTaxaNames$speciesCheckResult1[rowidx] = "FAIL"
#       print(paste0("Failed! original name case not in \'synonym\' or \'unresolved\' ", indv_spname))
#     }
# 
#   } else {
#     spTaxaNames$speciesCheckResult1[rowidx] = "FAIL"
#     print(paste0("Failed! no match: ", indv_spname))
#   }
# 
#   # break
# }
# ------------------------------------------------------------------------------

# Step1: Curate spTaxaNames and seperate into genus and species names ----------
# 1. Camphora parthenoxylon -> Cinnamomum porrectum
# Ref: https://asianplant.net/Lauraceae/Cinnamomum_parthenoxylon.htm
#
# 2. Elaeocarpus sylvestis -> Elaeocarpus sylvestris (already solved in raw data)
spTaxaNames_corrected = traitDataIndv_spavg_log %>% dplyr::select(SpeciesFullName) %>% distinct() %>% 
  mutate(SpeciesFullNameNew = gsub("Camphora parthenoxylon", "Cinnamomum porrectum", traitDataIndv_spavg_log$SpeciesFullName)) %>%
  dplyr::rename(speciesFullName_ori = SpeciesFullName, speciesFullName = SpeciesFullNameNew) %>%
  dplyr::mutate(speciesFullName_tosep = speciesFullName) # for separation
# Add genus and species name separately
spTaxaNames_corrected = spTaxaNames_corrected %>% separate(speciesFullName_tosep, into = c("genus_name", "sp_name"), sep=" ")

# Step2: Add family information ------------------------------------------------
genusFamilyTable = read.table("Phylomaker_genus_family.csv", sep=",", header=T)
spTaxaNames_withfam = spTaxaNames_corrected %>%
  left_join(genusFamilyTable, by = c("genus_name" = "Genus")) %>%
  dplyr::rename(family_name = Family) %>% 
  dplyr::select(-Group, -Source) %>% # remove extra columns
  dplyr::select(speciesFullName_ori, speciesFullName, family_name, genus_name, sp_name) # re-order columns
rm(spTaxaNames_corrected) # remove useless vars


# Step3: Also replace species names in the original dataset --------------------
# 3.1 - dataset normalized by sp only
traitDataIndv_spavg = traitDataIndv_spavg %>% dplyr::rename(speciesFullName_ori = SpeciesFullName) %>%
  left_join(spTaxaNames_withfam, by = "speciesFullName_ori") %>%
  dplyr::select(speciesFullName, family_name, genus_name, sp_name, everything()) # re-order columns

traitDataIndv_spavg_log = traitDataIndv_spavg_log %>% dplyr::rename(speciesFullName_ori = SpeciesFullName) %>%
  left_join(spTaxaNames_withfam, by = "speciesFullName_ori") %>%
  dplyr::select(speciesFullName, family_name, genus_name, sp_name, everything()) # re-order columns

traitDataIndv_spavg_log_ztransform = traitDataIndv_spavg_log_ztransform %>% dplyr::rename(speciesFullName_ori = SpeciesFullName) %>%
  left_join(spTaxaNames_withfam, by = "speciesFullName_ori") %>%
  dplyr::select(speciesFullName, family_name, genus_name, sp_name, everything()) # re-order columns

# 3.2 - dataset normalized by sp and gf
traitDataIndv_spgfavg = traitDataIndv_spgfavg %>% dplyr::rename(speciesFullName_ori = SpeciesFullName) %>%
  left_join(spTaxaNames_withfam, by = "speciesFullName_ori") %>%
  dplyr::select(speciesFullName, family_name, genus_name, sp_name, everything()) # re-order columns

traitDataIndv_spgfavg_log = traitDataIndv_spgfavg_log %>% dplyr::rename(speciesFullName_ori = SpeciesFullName) %>%
  left_join(spTaxaNames_withfam, by = "speciesFullName_ori") %>%
  dplyr::select(speciesFullName, family_name, genus_name, sp_name, everything()) # re-order columns

traitDataIndv_spgfavg_log_ztransform = traitDataIndv_spgfavg_log_ztransform %>% dplyr::rename(speciesFullName_ori = SpeciesFullName) %>%
  left_join(spTaxaNames_withfam, by = "speciesFullName_ori") %>%
  dplyr::select(speciesFullName, family_name, genus_name, sp_name, everything()) # re-order columns

# 3.3 - the original, not species averaged dataset
traitDataIndv = traitDataIndv %>% dplyr::rename(speciesFullName_ori = SpeciesFullName) %>%
  left_join(spTaxaNames_withfam, by = "speciesFullName_ori") %>%
  dplyr::select(speciesFullName, family_name, genus_name, sp_name, everything()) # re-order columns

# Step4: Build phylogenetic tree using PhyloMaker2 -----------------------------
# TODO: Update 20/11-24: all species are now bound to the tree after the addition of
# family names from genus names
spTaxaNames_tobuild = spTaxaNames_withfam %>%
  dplyr::rename(species = speciesFullName, genus = genus_name, family = family_name) %>%
  dplyr::select(-sp_name) %>%
  dplyr::select(species,genus,family) # re-order

phylotree_result = phylo.maker(sp.list = spTaxaNames_tobuild, tree = GBOTB.extended.TPL, nodes = nodes.info.1.TPL)
# phylotree_result[[1]] is the resulting phylogenetic tree.
# summary(phylotree_result[[1]])

# Save data
save(traitDataIndv,
     traitDataIndv_spavg, traitDataIndv_spavg_log, traitDataIndv_spavg_log_ztransform,
     traitDataIndv_spgfavg, traitDataIndv_spgfavg_log, traitDataIndv_spgfavg_log_ztransform,
     spTaxaNames_withfam, phylotree_result, 
     file = "traitDataFujian-spavg-phylo-step2.RData")

