rm(list=ls())
library(tidyverse)
library(dplyr)

# ------------------------------------------------------------------------------
# Load individual- and species- level traits
# ------------------------------------------------------------------------------

# Species-level, code is modified, retaining original suffix '_Indv' ----
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
# Trait value pre-processing
# ------------------------------------------------------------------------------
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
# Curate the final dataframe
traitDataIndv = cbind(traitDataIndv_metadata, traitDataIndv_numonly)
rm(traitDataIndv_metadata, traitDataIndv_numonly, traitDataIndv_colcaption)

# calculate the species mean values for each trait, discarding NA values.
# 1. group by species only, give duplicated growth form entries `tree+shrub`
traitDataIndv_spavg = traitDataIndv %>%
  group_by(`SpeciesFullName`) %>%
  dplyr::summarize(
    GrowthForm = first(GrowthForm),
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)) # Calculate mean for all numeric columns
  )
# 2. group by species and growth form
traitDataIndv_spgfavg = traitDataIndv %>%
  group_by(`SpeciesFullName`, `GrowthForm`) %>%
  dplyr::summarize(
    GrowthForm = first(GrowthForm),
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)) # Calculate mean for all numeric columns
  )

# give `tree+shrub` notation for duplicated entries
traitDataIndv_spavg$GrowthForm = as.character(traitDataIndv_spavg$GrowthForm)
traitDataIndv_spavg = traitDataIndv_spavg %>% dplyr::mutate(
  GrowthForm = ifelse(SpeciesFullName %in% c("Machilus pauhoi", "Castanopsis fordii"), "tree+shrub", GrowthForm)
)
traitDataIndv_spavg$GrowthForm = as.factor(traitDataIndv_spavg$GrowthForm)

# # TODO: [TESTING PURPOSES - sum species, compare the 2 sites
# traitDataIndv_spavg_top = traitDataIndv %>%
#   filter(grepl("^SH", SpeciesID)) %>%
#   group_by(`SpeciesFullName`) %>%
#   dplyr::summarize(
#     GrowthForm = first(GrowthForm),
#     across(where(is.numeric), ~ mean(.x, na.rm = TRUE)) # Calculate mean for all numeric columns
#   )
# 
# traitDataIndv_spavg_bottom = traitDataIndv %>%
#   filter(grepl("^WH", SpeciesID)) %>%
#   group_by(`SpeciesFullName`) %>%
#   dplyr::summarize(
#     GrowthForm = first(GrowthForm),
#     across(where(is.numeric), ~ mean(.x, na.rm = TRUE)) # Calculate mean for all numeric columns
#   )

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
traitDataIndv_spavg_log = func_log_transform(traitDataIndv_spavg)
traitDataIndv_spgfavg_log = func_log_transform(traitDataIndv_spgfavg)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Normalize data for bivariate analysis and PCA
# Ref: bergmann et al. Sci Adv, 2020
#
#1. Log-transform each trait value.
#2. Z-transform each log-transformed trait to have a mean of 0 and an SD of 1.
func_log_z_transform <- function(obj_df) {

  numeric_cols <- sapply(obj_df, is.numeric)
  numeric_data <- obj_df[, numeric_cols]
  
  # Log-transform each trait (signed log transformation, dealing with negative values)
  log_transformed = sign(numeric_data) * log(abs(numeric_data) + 1)
  # Z-transform the log-transformed data
  z_transformed <- scale(log_transformed, center = TRUE, scale = TRUE)

  obj_df[, numeric_cols] <- z_transformed
  return(obj_df)        
}
traitDataIndv_spavg_log_ztransform = func_log_z_transform(traitDataIndv_spavg)
traitDataIndv_spgfavg_log_ztransform = func_log_z_transform(traitDataIndv_spgfavg)
# ------------------------------------------------------------------------------

# testing purposes
# sd(traitDataIndv_normalized$`vH(cm2/cm2)`, na.rm = T)
# mean(traitDataIndv_normalized$`vH(cm2/cm2)`, na.rm = T)

# Save data
save(traitDataIndv, #traitDataIndv_colcaption,
     traitDataIndv_spavg, traitDataIndv_spavg_log, traitDataIndv_spavg_log_ztransform,
     traitDataIndv_spgfavg, traitDataIndv_spgfavg_log, traitDataIndv_spgfavg_log_ztransform,
     file = "species-level-code/traitDataFujian-normbysp-step1.RData")
