# ****************************************************************
# TODO: Modified for use with individual-level analysis workflow
# in response to HW's advice, 2025.2
#
# Keeping the full 87 coupled aboveground and belowground sampled
# individuals for analysis
# 
# 
# UPDATE 2025.03.09 - combined species averaged data
# ****************************************************************

rm(list=ls())
library(tidyverse)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# Load individual- and species- level traits
# Species-level, code is modified, retaining original suffix '_Indv'
traitDataIndv = readxl::read_xlsx("rawdata/traits-individual.xlsx")
traitDataIndv_colcaption = traitDataIndv[1,]
traitDataIndv = traitDataIndv[-1,] 

# Convert trait measurements to as.numeric format
traitDataIndv_numonly = traitDataIndv[,5:ncol(traitDataIndv)] %>% mutate_all(as.numeric)
# Convert all values in df that are equal to 0 to NA
traitDataIndv_numonly[traitDataIndv_numonly == 0] <- NA

# Extract metadata, and add a 'tree/shrub' indicator for subsequent grouping
traitDataIndv_metadata = traitDataIndv[,1:4]
tree_samplenames = c("SH001", "SH002", "SH003", "SH004", "SH005", "SH006", "SH007",
                     "SH009", "SH010", "SH011", "SH012", "SH013", "SH014",
                     "WH901", "WH902", "WH903", "WH904", "WH905", "WH906", "WH907",
                     "WH908", "WH909", "WH910", "WH911", "WH912", "WH913")
liana_samplenames = c("SH033", "SH034", "WH931", "WH932")

traitDataIndv_metadata$GrowthForm = "shrub" # add new column
traitDataIndv_metadata = traitDataIndv_metadata %>%
  mutate(GrowthForm = ifelse(SpeciesID %in% tree_samplenames, "tree", ifelse(SpeciesID %in% liana_samplenames, "liana", "shrub")))
traitDataIndv_metadata$GrowthForm = as.factor(traitDataIndv_metadata$GrowthForm)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Trait value pre-processing
# Add LMA notation (1/SLA), raw value
traitDataIndv_numonly$LMA = 1/traitDataIndv_numonly$SLA
# Convert Rdark25 and Rlight25 to positive value (for equal meaning with other respiration traits)
traitDataIndv_numonly$Rdark25P = -traitDataIndv_numonly$Rdark25
traitDataIndv_numonly$Rlight25P = -traitDataIndv_numonly$Rlight25

# Convert mass-basis traits to area basis
# SLA: cm2 leaf/g leaf, LNC: % - LNC/100 is g N/g leaf
# SRA: cm2 root/g root, RNC: % - RNC/100 is g N/g root
# 
# LNCarea: g N/cm2 leaf -> = LNC/100 / SLA
# RNCarea: g N/cm2 root -> = RNC/100 / SRA
traitDataIndv_numonly = traitDataIndv_numonly %>% 
  dplyr::mutate(RNCarea = (RNC/100)/SRA, LNCarea = (LNC/100)/SLA)
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Curate the final dataframe (individual level, 87 rows)
traitDataIndv = cbind(traitDataIndv_metadata, traitDataIndv_numonly)
rm(traitDataIndv_metadata, traitDataIndv_numonly, traitDataIndv_colcaption)

# Added 7/3-25: fix the species name
# 1. Camphora parthenoxylon -> Cinnamomum porrectum
# Ref: https://asianplant.net/Lauraceae/Cinnamomum_parthenoxylon.htm
#
# 2. Elaeocarpus sylvestis -> Elaeocarpus sylvestris (already solved in raw data)
traitDataIndv = traitDataIndv %>%
  dplyr::mutate(SpeciesFullName = gsub("Camphora parthenoxylon", "Cinnamomum porrectum", SpeciesFullName)) %>%
  dplyr::mutate(SpeciesFullName_tosep = SpeciesFullName) %>%
  separate(SpeciesFullName_tosep, into = c("genus_name", "sp_name"), sep=" ") # will destroy column: SpeciesFullName_tosep

# Add Family information
genusFamilyTable = read.table("rawdata/PhyloMaker_genus_family.csv", sep=",", header=T)
traitDataIndv = traitDataIndv %>% left_join(genusFamilyTable, by = c("genus_name" = "Genus")) %>%
  rename(family_name = Family) %>% dplyr::select(-Group, -Source) # remove extra columns


traitDataIndv_SelectedTraits = traitDataIndv %>%
  dplyr::mutate(SiteID = ifelse(grepl("^SH", SpeciesID), "hilltop", "valley")) %>%
  dplyr::select(SampleID, SpeciesID, family_name, genus_name, sp_name, SpeciesFullName, SpeciesCH, GrowthForm, SiteID, # Metadata
                RTD,SRL,RD,RNC,SRA,RPC,RCC,RDMC,SRR25, # FR traits
                SLA,LMA,LNC,LCC,LPC,Ld13C,Asat,Rdark25P,Vcmax25, # leaf traits
                H,DBH,LA # Size-related traits
  )

# ------------------------------------------------------------------------------
# calculate the species mean values for each trait, discarding NA values.
# group by species only, give duplicated growth form entries `tree+shrub`
traitDataIndv_SelectedTraits_spavg = traitDataIndv_SelectedTraits %>%
  group_by(`SpeciesFullName`) %>%
  dplyr::summarize(
    family_name = first(family_name), genus_name = first(genus_name), 
    sp_name = first(sp_name), GrowthForm = first(GrowthForm),
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)) # Calculate mean for all numeric columns
  )

# give `tree+shrub` notation for duplicated entries
traitDataIndv_SelectedTraits_spavg$GrowthForm = as.character(traitDataIndv_SelectedTraits_spavg$GrowthForm)
traitDataIndv_SelectedTraits_spavg = traitDataIndv_SelectedTraits_spavg %>% dplyr::mutate(
  GrowthForm = ifelse(SpeciesFullName %in% c("Machilus pauhoi", "Castanopsis fordii"), "tree+shrub", GrowthForm)
)
traitDataIndv_SelectedTraits_spavg$GrowthForm = as.factor(traitDataIndv_SelectedTraits_spavg$GrowthForm)


# ------------------------------------------------------------------------------
# Prepare only log-transformed data
# function for data normalization (log-transform only)
func_log_transform <- function(obj_df) {
  numeric_cols <- sapply(obj_df, is.numeric)
  numeric_data <- obj_df[, numeric_cols]
  # Log-transform each trait (signed log transformation, dealing with negative values)
  log_transformed = sign(numeric_data) * log(abs(numeric_data) + 1)
  obj_df[, numeric_cols] <- log_transformed
  return(obj_df)        
}

# TODO: USE - SMA regression (no scale)
traitDataIndv_SelectedTraits_log = func_log_transform(traitDataIndv_SelectedTraits)
# TODO: USE - PGLS regression (no scale)
traitDataIndv_SelectedTraits_spavg_log = func_log_transform(traitDataIndv_SelectedTraits_spavg)


# Normalize data for bivariate analysis and PCA
# Ref: bergmann et al. Sci Adv, 2020
#
#1. Log-transform each trait value.
#2. scale each log-transformed trait to have a mean of 0 and an SD of 1.
func_log_scale <- function(obj_df) {
  numeric_cols <- sapply(obj_df, is.numeric)
  numeric_data <- obj_df[, numeric_cols]
  # Log-transform each trait (signed log transformation, dealing with negative values)
  log_transformed = sign(numeric_data) * log(abs(numeric_data) + 1)
  # Z-transform the log-transformed data
  scaled_data <- scale(log_transformed, center = TRUE, scale = TRUE)
  obj_df[, numeric_cols] <- scaled_data
  return(obj_df)        
}

# TODO: USE - PCA
traitDataIndv_SelectedTraits_log_scale = func_log_scale(traitDataIndv_SelectedTraits)
# TODO: USE - Phylo-informed PCA
traitDataIndv_SelectedTraits_log_scale = func_log_scale(traitDataIndv_SelectedTraits_spavg)


# # TODO: 6/3/25, Statistics for individual-level data ----
# traitDataIndv_uniqspecies = traitDataIndv %>% dplyr::select(SpeciesFullName) %>% distinct()
# traitDataIndv_site1 = traitDataIndv %>% filter(grepl("^SH", SpeciesID)) # hilltop
# traitDataIndv_site2 = traitDataIndv %>% filter(grepl("^WH", SpeciesID)) # valley

# Save data
save(traitDataIndv, traitDataIndv_SelectedTraits, 
     traitDataIndv_SelectedTraits_log, traitDataIndv_SelectedTraits_log_scale,
     file = "individual-level-code/traitDataFujian-Ind-step1.RData")

save(traitDataIndv, traitDataIndv_SelectedTraits_spavg,
     traitDataIndv_SelectedTraits_spavg_log, traitDataIndv_SelectedTraits_log_scale,
     file = "individual-level-code/traitDataFujian-SpAvg-step1.RData")

