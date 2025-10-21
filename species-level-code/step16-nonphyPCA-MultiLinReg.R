rm(list=ls())
library(tidyverse)
library(dplyr)
library(factoextra) # perform PCA analysis
library(FactoMineR) # perform PCA analysis
library(relaimpo) # Variable Importance Plot
library(ggthemes)
library(patchwork)

# ------------------------------------------------------------------------------
# Perform PCA analysis
# modified in response to Sandy Harrison's comment -> Use original values without normalization
# to perform the PCA to show if the conclusion was altered
# ------------------------------------------------------------------------------
load("traitDataFujian-fmtdata-step3.RData")

# TODO: remove LIANA entries in the data due to its minority in PCA analysis
# traitDataIndv_normalized = traitDataIndv_normalized %>% filter(GrowthForm != "liana") # (41->37 species, but only 1 contain FR trait measurements)

# # FILTER OUTLIERS (for tentative analysis, 13/10-24)
# traitDataIndv_normalized = traitDataIndv_normalized %>% 
#   filter(RTD > -2)

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
# ------------------------------------------------------------------------------

### SET DATA TO USE ###
traitDataPCA_touse = traitDataIndv_spavg_log

# Perform PCA
# Curate data for PCA analysis
nonphyloPCAData_numonly = traitDataPCA_touse %>% dplyr::select(RTD,SRL,RD,RNC,LMA,LNC,RDMC,LPC)
rownames(nonphyloPCAData_numonly) = traitDataPCA_touse$speciesFullName
nonphyloPCAData_numonly = nonphyloPCAData_numonly %>% na.omit() # remove lines /w NA values

# used for growth form annotation in PCA biplot
nonphyloPCAData_meta_rmna = traitDataPCA_touse %>% 
  dplyr::select(speciesFullName,GrowthForm,RTD,SRL,RD,RNC,LMA,LNC,RDMC,LPC) %>% na.omit()

nonphyloPCAresult = PCA(nonphyloPCAData_numonly, scale.unit = T, graph = F) 
# ------------------------------------------------------------------------------

# Regression Analysis Using Principal Components
# Extract PCs as predictors
PC_scores = as.data.frame(nonphyloPCAresult$ind$coord)
PC_scores$speciesFullName = nonphyloPCAData_meta_rmna$speciesFullName

# Merge with original data to retain response vars
traitDataPCA_MLR_touse = nonphyloPCAData_meta_rmna %>% left_join(PC_scores, by = "speciesFullName")

# Add response variable, which dont participate PCA
traitDataPCA_MLR_responsevar = traitDataPCA_touse %>% dplyr::select(speciesFullName,SRR25,Rdark25P)
traitDataPCA_MLR_touse = traitDataPCA_MLR_touse %>% left_join(traitDataPCA_MLR_responsevar, by = "speciesFullName")

traitDataPCA_MLR_touse_tree = traitDataPCA_MLR_touse %>% dplyr::filter(GrowthForm == "tree")
traitDataPCA_MLR_touse_shrub = traitDataPCA_MLR_touse %>% dplyr::filter(GrowthForm == "shrub")
traitDataPCA_MLR_touse_liana = traitDataPCA_MLR_touse %>% dplyr::filter(GrowthForm == "liana") # may be useless
# ------------------------------------------------------------------------------

# PART 1: Predicting SRR25 using multiple liear regression (MLR, PC1-3) ----
model.MLR.SRR25_PC1to3_total = lm(SRR25 ~ Dim.1 + Dim.2 + Dim.3, data = traitDataPCA_MLR_touse)
summary(model.MLR.SRR25_PC1to3_total)

imp.MLR.SRR25_PC1to3_total = calc.relimp(model.MLR.SRR25_PC1to3_total, type = "lmg")
imp.MLR.SRR25_PC1to3_total_data = as.data.frame(imp.MLR.SRR25_PC1to3_total$lmg) 
imp.MLR.SRR25_PC1to3_total_data$dims = rownames(imp.MLR.SRR25_PC1to3_total_data)
colnames(imp.MLR.SRR25_PC1to3_total_data) = c("values", "dims")
plt_imp.MLR.SRR25_PC1to3_total = ggplot(imp.MLR.SRR25_PC1to3_total_data, aes(x = dims, y = values)) +
  geom_bar(stat = "identity", fill = "#11999e") +
  theme_minimal() +
  labs(x = "PC", y = "Relative importance")


model.MLR.SRR25_PC1to3_tree = lm(SRR25 ~ Dim.1 + Dim.2 + Dim.3, data = traitDataPCA_MLR_touse_tree)
summary(model.MLR.SRR25_PC1to3_tree)

imp.MLR.SRR25_PC1to3_tree = calc.relimp(model.MLR.SRR25_PC1to3_tree, type = "lmg")
imp.MLR.SRR25_PC1to3_tree_data = as.data.frame(imp.MLR.SRR25_PC1to3_tree$lmg) 
imp.MLR.SRR25_PC1to3_tree_data$dims = rownames(imp.MLR.SRR25_PC1to3_tree_data)
colnames(imp.MLR.SRR25_PC1to3_tree_data) = c("values", "dims")
plt_imp.MLR.SRR25_PC1to3_tree = ggplot(imp.MLR.SRR25_PC1to3_tree_data, aes(x = dims, y = values)) +
  geom_bar(stat = "identity", fill = "#11999e") +
  theme_minimal() +
  labs(x = "PC", y = "Relative importance")


model.MLR.SRR25_PC1to3_shrub = lm(SRR25 ~ Dim.1 + Dim.2 + Dim.3, data = traitDataPCA_MLR_touse_shrub)
summary(model.MLR.SRR25_PC1to3_shrub)

imp.MLR.SRR25_PC1to3_shrub = calc.relimp(model.MLR.SRR25_PC1to3_shrub, type = "lmg")
imp.MLR.SRR25_PC1to3_shrub_data = as.data.frame(imp.MLR.SRR25_PC1to3_shrub$lmg) 
imp.MLR.SRR25_PC1to3_shrub_data$dims = rownames(imp.MLR.SRR25_PC1to3_shrub_data)
colnames(imp.MLR.SRR25_PC1to3_shrub_data) = c("values", "dims")
plt_imp.MLR.SRR25_PC1to3_shrub = ggplot(imp.MLR.SRR25_PC1to3_shrub_data, aes(x = dims, y = values)) +
  geom_bar(stat = "identity", fill = "#11999e") +
  theme_minimal() +
  labs(x = "PC", y = "Relative importance")

plt_imp.MLR.SRR25_PC1to3 = plt_imp.MLR.SRR25_PC1to3_total + plt_imp.MLR.SRR25_PC1to3_tree + plt_imp.MLR.SRR25_PC1to3_shrub
ggsave(filename = "plt_imp.MLR.SRR25_PC1to3.pdf", plot = plt_imp.MLR.SRR25_PC1to3, width=6, height=2)


# PART 2: Predicting Rdark25P using multiple liear regression (MLR, PC1-3) ----
model.MLR.Rdark25P_PC1to3_total = lm(Rdark25P ~ Dim.1 + Dim.2 + Dim.3, data = traitDataPCA_MLR_touse)
summary(model.MLR.Rdark25P_PC1to3_total)

imp.MLR.Rdark25P_PC1to3_total = calc.relimp(model.MLR.Rdark25P_PC1to3_total, type = "lmg")
imp.MLR.Rdark25P_PC1to3_total_data = as.data.frame(imp.MLR.Rdark25P_PC1to3_total$lmg) 
imp.MLR.Rdark25P_PC1to3_total_data$dims = rownames(imp.MLR.Rdark25P_PC1to3_total_data)
colnames(imp.MLR.Rdark25P_PC1to3_total_data) = c("values", "dims")
plt_imp.MLR.Rdark25P_PC1to3_total = ggplot(imp.MLR.Rdark25P_PC1to3_total_data, aes(x = dims, y = values)) +
  geom_bar(stat = "identity", fill = "#f38181") +
  theme_minimal() +
  labs(x = "PC", y = "Relative importance")


model.MLR.Rdark25P_PC1to3_tree = lm(Rdark25P ~ Dim.1 + Dim.2 + Dim.3, data = traitDataPCA_MLR_touse_tree)
summary(model.MLR.Rdark25P_PC1to3_tree)

imp.MLR.Rdark25P_PC1to3_tree = calc.relimp(model.MLR.Rdark25P_PC1to3_tree, type = "lmg")
imp.MLR.Rdark25P_PC1to3_tree_data = as.data.frame(imp.MLR.Rdark25P_PC1to3_tree$lmg) 
imp.MLR.Rdark25P_PC1to3_tree_data$dims = rownames(imp.MLR.Rdark25P_PC1to3_tree_data)
colnames(imp.MLR.Rdark25P_PC1to3_tree_data) = c("values", "dims")
plt_imp.MLR.Rdark25P_PC1to3_tree = ggplot(imp.MLR.Rdark25P_PC1to3_tree_data, aes(x = dims, y = values)) +
  geom_bar(stat = "identity", fill = "#f38181") +
  theme_minimal() +
  labs(x = "PC", y = "Relative importance")


model.MLR.Rdark25P_PC1to3_shrub = lm(Rdark25P ~ Dim.1 + Dim.2 + Dim.3, data = traitDataPCA_MLR_touse_shrub)
summary(model.MLR.Rdark25P_PC1to3_shrub)

imp.MLR.Rdark25P_PC1to3_shrub = calc.relimp(model.MLR.Rdark25P_PC1to3_shrub, type = "lmg")
imp.MLR.Rdark25P_PC1to3_shrub_data = as.data.frame(imp.MLR.Rdark25P_PC1to3_shrub$lmg) 
imp.MLR.Rdark25P_PC1to3_shrub_data$dims = rownames(imp.MLR.Rdark25P_PC1to3_shrub_data)
colnames(imp.MLR.Rdark25P_PC1to3_shrub_data) = c("values", "dims")
plt_imp.MLR.Rdark25P_PC1to3_shrub = ggplot(imp.MLR.Rdark25P_PC1to3_shrub_data, aes(x = dims, y = values)) +
  geom_bar(stat = "identity", fill = "#f38181") +
  theme_minimal() +
  labs(x = "PC", y = "Relative importance")

plt_imp.MLR.Rdark25P_PC1to3 = plt_imp.MLR.Rdark25P_PC1to3_total + plt_imp.MLR.Rdark25P_PC1to3_tree + plt_imp.MLR.Rdark25P_PC1to3_shrub
ggsave(filename = "plt_imp.MLR.Rdark25P_PC1to3.pdf", plot = plt_imp.MLR.Rdark25P_PC1to3, width=6, height=2)
