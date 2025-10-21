rm(list=ls())
library(tidyverse)
library(dplyr)
library(factoextra) # perform PCA analysis
library(FactoMineR) # perform PCA analysis
# library(rgl) # 3-D PCA plot
library(smatr) # SMA regression
library(ggthemes)
library(patchwork) # plot merging
library(cowplot) # plot merging

# ------------------------------------------------------------------------------
# Perform PCA analysis
# modified in response to Sandy Harrison's comment -> Use original values without normalization
# to perform the PCA to show if the conclusion was altered
# ------------------------------------------------------------------------------
load("traitDataFujian-spavg-phylo-step2.RData")

# TODO: remove LIANA entries in the data due to its minority in PCA analysis
# traitDataIndv_normalized = traitDataIndv_normalized %>% filter(GrowthForm != "liana") # (41->37 species, but only 1 contain FR trait measurements)

# # FILTER OUTLIERS (for tentative analysis, 13/10-24)
# traitDataIndv_normalized = traitDataIndv_normalized %>% 
#   filter(RTD > -2)

### SET DATA TO USE ###
traitDataPCA_touse = traitDataIndv_spgfavg_log

# Resolve problem: non-unique values when setting 'row.names': ‘Castanopsis fordii’, ‘Machilus pauhoi’ 
# because of sp-gf averaging
traitDataPCA_touse = traitDataPCA_touse %>% mutate(
  PCAIdentifier_spgf = paste(gsub(" ","_",speciesFullName), GrowthForm, sep="-")
)

# Curate data for PCA analysis
nonphyloPCAData_numonly = traitDataPCA_touse %>% ungroup() %>% 
  dplyr::select(PCAIdentifier_spgf,GrowthForm, RTD,SRL,RD,RNC,LMA,LNC,RDMC,LPC) %>% 
  na.omit()

rownames(nonphyloPCAData_numonly) = nonphyloPCAData_numonly$PCAIdentifier_spgf

# used for growth form annotation in PCA biplot
nonphyloPCAData_meta_rmna = nonphyloPCAData_numonly
nonphyloPCAData_numonly = nonphyloPCAData_numonly %>% dplyr::select(-PCAIdentifier_spgf, -GrowthForm)

# perform PCA
nonphyloPCAresult = PCA(nonphyloPCAData_numonly, scale.unit = T, graph = F) 
# ------------------------------------------------------------------------------

# Regression Analysis Using Principal Components
# Extract PCs as predictors
PC_scores = as.data.frame(nonphyloPCAresult$ind$coord)
PC_scores$PCAIdentifier_spgf = nonphyloPCAData_meta_rmna$PCAIdentifier_spgf

# Merge with original data to retain response vars
traitDataPCA_SMA_touse = nonphyloPCAData_meta_rmna %>% left_join(PC_scores, by = "PCAIdentifier_spgf")

# Add response variable, which dont participate PCA
traitDataPCA_SMA_responsevar = traitDataPCA_touse %>% dplyr::select(PCAIdentifier_spgf,SRR25,Rdark25P)
traitDataPCA_SMA_touse = traitDataPCA_SMA_touse %>% left_join(traitDataPCA_SMA_responsevar, by = "PCAIdentifier_spgf")

traitDataPCA_SMA_touse_tree = traitDataPCA_SMA_touse %>% dplyr::filter(GrowthForm == "tree")
traitDataPCA_SMA_touse_shrub = traitDataPCA_SMA_touse %>% dplyr::filter(GrowthForm == "shrub")
traitDataPCA_SMA_touse_liana = traitDataPCA_SMA_touse %>% dplyr::filter(GrowthForm == "liana") # may be useless

# ------------------------------------------------------------------------------
# PART 1: Predicting SRR25 using PC1 (RES conservation, RES+LES with RDMC and RPC) ----

model.sma.SRR25_Dim.1_total = sma(traitDataPCA_SMA_touse$SRR25 ~ traitDataPCA_SMA_touse$Dim.1)
model.sma.SRR25_Dim.1_tree = sma(traitDataPCA_SMA_touse_tree$SRR25 ~ traitDataPCA_SMA_touse_tree$Dim.1)
model.sma.SRR25_Dim.1_shrub = sma(traitDataPCA_SMA_touse_shrub$SRR25 ~ traitDataPCA_SMA_touse_shrub$Dim.1)

#plot(model.sma.SRR25_Dim.1)
summary(model.sma.SRR25_Dim.1_total)
summary(model.sma.SRR25_Dim.1_tree)
summary(model.sma.SRR25_Dim.1_shrub)

smaplt_total_intercept = model.sma.SRR25_Dim.1_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.SRR25_Dim.1_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.SRR25_Dim.1_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.SRR25_Dim.1_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.SRR25_Dim.1_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.SRR25_Dim.1_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.SRR25_Dim.1_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.SRR25_Dim.1_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.SRR25_Dim.1_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.SRR25_Dim.1_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.SRR25_Dim.1_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.SRR25_Dim.1_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.SRR25_Dim.1_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.SRR25_Dim.1_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.SRR25_Dim.1_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.SRR25_Dim.1_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.SRR25_Dim.1_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.SRR25_Dim.1_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitDataPCA_SMA_touse %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, SRR25, Dim.1) %>%
  dplyr::mutate(Dim.1 = as.numeric(Dim.1), SRR25 = as.numeric(SRR25))
smaplt_total_data$predicted_SRR25 = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$Dim.1
smaplt_total_data$predicted_SRR25_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$Dim.1
smaplt_total_data$predicted_SRR25_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$Dim.1

smaplt_tree_data = traitDataPCA_SMA_touse_tree %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, SRR25, Dim.1) %>%
  dplyr::mutate(Dim.1 = as.numeric(Dim.1), SRR25 = as.numeric(SRR25))
smaplt_tree_data$predicted_SRR25 = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$Dim.1
smaplt_tree_data$predicted_SRR25_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$Dim.1
smaplt_tree_data$predicted_SRR25_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$Dim.1

smaplt_shrub_data = traitDataPCA_SMA_touse_shrub %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, SRR25, Dim.1) %>%
  dplyr::mutate(Dim.1 = as.numeric(Dim.1), SRR25 = as.numeric(SRR25))
smaplt_shrub_data$predicted_SRR25 = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$Dim.1
smaplt_shrub_data$predicted_SRR25_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$Dim.1
smaplt_shrub_data$predicted_SRR25_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$Dim.1

# Annotation, R2 and p-value
model.sma.SRR25_Dim.1_annotlabel = paste0(
  "Merged: R2=", round(as.numeric(model.sma.SRR25_Dim.1_total[["r2"]]), digits = 2), " p",
  ifelse(signif(as.numeric(model.sma.SRR25_Dim.1_total[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.SRR25_Dim.1_total[["pval"]]), digits = 2), "\n")),
  
  "Tree: R2=", round(as.numeric(model.sma.SRR25_Dim.1_tree[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.SRR25_Dim.1_tree[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.SRR25_Dim.1_tree[["pval"]]), digits = 2), "\n")),
  
  "Shrub: R2=", round(as.numeric(model.sma.SRR25_Dim.1_shrub[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.SRR25_Dim.1_shrub[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.SRR25_Dim.1_shrub[["pval"]]), digits = 2), "\n"))
)

smaplt.SRR25_Dim.1 = ggplot(smaplt_total_data, aes(x=Dim.1, y=SRR25)) +
  geom_point(aes(color=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$Dim.1, na.rm=T), max(smaplt_total_data$Dim.1, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$Dim.1, na.rm=T), max(smaplt_tree_data$Dim.1, na.rm=T)),
    color = "#F07673", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$Dim.1, na.rm=T), max(smaplt_shrub_data$Dim.1, na.rm=T)),
    color = "#7998AD", linewidth = 1.0
  ) +
  labs(x = "PC1 (RES+LES+RDMC+LPC)", y = "logSRR25", subtitle = "SRR25-PC1\n(mainly root conservation)") + 
  theme_classic() + theme(legend.position = "none") +
  scale_color_manual(values = c("#EECA40", "#7998AD", "#F07673")) + 
  annotate("text", x = 0.6, y = 4.0, label = model.sma.SRR25_Dim.1_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.SRR25_Dim.1, filename = "smaplt.SRR25_PC1-RES-con.pdf", width = 3, height = 2.9)
ggsave(plot = smaplt.SRR25_Dim.1, filename = "smaplt.SRR25_PC1-RES-con.png", width = 3, height = 2.9, dpi=300)

# ------------------------------------------------------------------------------
# PART 2: Predicting SRR25 using PC2 (LES, RES+LES with RDMC and RPC) ----
 
model.sma.SRR25_Dim.2_total = sma(traitDataPCA_SMA_touse$SRR25 ~ traitDataPCA_SMA_touse$Dim.2)
model.sma.SRR25_Dim.2_tree = sma(traitDataPCA_SMA_touse_tree$SRR25 ~ traitDataPCA_SMA_touse_tree$Dim.2)
model.sma.SRR25_Dim.2_shrub = sma(traitDataPCA_SMA_touse_shrub$SRR25 ~ traitDataPCA_SMA_touse_shrub$Dim.2)

#plot(model.sma.SRR25_Dim.2)
summary(model.sma.SRR25_Dim.2_total)
summary(model.sma.SRR25_Dim.2_tree)
summary(model.sma.SRR25_Dim.2_shrub)

smaplt_total_intercept = model.sma.SRR25_Dim.2_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.SRR25_Dim.2_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.SRR25_Dim.2_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.SRR25_Dim.2_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.SRR25_Dim.2_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.SRR25_Dim.2_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.SRR25_Dim.2_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.SRR25_Dim.2_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.SRR25_Dim.2_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.SRR25_Dim.2_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.SRR25_Dim.2_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.SRR25_Dim.2_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.SRR25_Dim.2_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.SRR25_Dim.2_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.SRR25_Dim.2_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.SRR25_Dim.2_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.SRR25_Dim.2_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.SRR25_Dim.2_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitDataPCA_SMA_touse %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, SRR25, Dim.2) %>%
  dplyr::mutate(Dim.2 = as.numeric(Dim.2), SRR25 = as.numeric(SRR25))
smaplt_total_data$predicted_SRR25 = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$Dim.2
smaplt_total_data$predicted_SRR25_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$Dim.2
smaplt_total_data$predicted_SRR25_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$Dim.2

smaplt_tree_data = traitDataPCA_SMA_touse_tree %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, SRR25, Dim.2) %>%
  dplyr::mutate(Dim.2 = as.numeric(Dim.2), SRR25 = as.numeric(SRR25))
smaplt_tree_data$predicted_SRR25 = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$Dim.2
smaplt_tree_data$predicted_SRR25_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$Dim.2
smaplt_tree_data$predicted_SRR25_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$Dim.2

smaplt_shrub_data = traitDataPCA_SMA_touse_shrub %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, SRR25, Dim.2) %>%
  dplyr::mutate(Dim.2 = as.numeric(Dim.2), SRR25 = as.numeric(SRR25))
smaplt_shrub_data$predicted_SRR25 = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$Dim.2
smaplt_shrub_data$predicted_SRR25_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$Dim.2
smaplt_shrub_data$predicted_SRR25_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$Dim.2

# Annotation, R2 and p-value
model.sma.SRR25_Dim.2_annotlabel = paste0(
  "Merged: R2=", round(as.numeric(model.sma.SRR25_Dim.2_total[["r2"]]), digits = 2), " p",
  ifelse(signif(as.numeric(model.sma.SRR25_Dim.2_total[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.SRR25_Dim.2_total[["pval"]]), digits = 2), "\n")),
  
  "Tree: R2=", round(as.numeric(model.sma.SRR25_Dim.2_tree[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.SRR25_Dim.2_tree[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.SRR25_Dim.2_tree[["pval"]]), digits = 2), "\n")),
  
  "Shrub: R2=", round(as.numeric(model.sma.SRR25_Dim.2_shrub[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.SRR25_Dim.2_shrub[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.SRR25_Dim.2_shrub[["pval"]]), digits = 2), "\n"))
)

smaplt.SRR25_Dim.2 = ggplot(smaplt_total_data, aes(x=Dim.2, y=SRR25)) +
  geom_point(aes(color=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$Dim.2, na.rm=T), max(smaplt_total_data$Dim.2, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$Dim.2, na.rm=T), max(smaplt_tree_data$Dim.2, na.rm=T)),
    color = "#F07673", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$Dim.2, na.rm=T), max(smaplt_shrub_data$Dim.2, na.rm=T)),
    color = "#7998AD", linewidth = 1.0
  ) +
  labs(x = "PC2 (RES+LES+RDMC+LPC)", y = "logSRR25", subtitle = "SRR25-PC2\n(mainly leaf fast-slow)") + 
  theme_classic() + theme(legend.position = "none") +
  scale_color_manual(values = c("#EECA40", "#7998AD", "#F07673")) + 
  annotate("text", x = -0.5, y = 4.5, label = model.sma.SRR25_Dim.2_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.SRR25_Dim.2, filename = "smaplt.SRR25_PC2-LES.pdf", width = 3, height = 2.9)
ggsave(plot = smaplt.SRR25_Dim.2, filename = "smaplt.SRR25_PC2-LES.png", width = 3, height = 2.9, dpi=300)

# ------------------------------------------------------------------------------
# PART 3: Predicting SRR25 using PC3 (RES collaboration, RES+LES with RDMC and RPC) ----

model.sma.SRR25_Dim.3_total = sma(traitDataPCA_SMA_touse$SRR25 ~ traitDataPCA_SMA_touse$Dim.3)
model.sma.SRR25_Dim.3_tree = sma(traitDataPCA_SMA_touse_tree$SRR25 ~ traitDataPCA_SMA_touse_tree$Dim.3)
model.sma.SRR25_Dim.3_shrub = sma(traitDataPCA_SMA_touse_shrub$SRR25 ~ traitDataPCA_SMA_touse_shrub$Dim.3)

#plot(model.sma.SRR25_Dim.3)
summary(model.sma.SRR25_Dim.3_total)
summary(model.sma.SRR25_Dim.3_tree)
summary(model.sma.SRR25_Dim.3_shrub)

smaplt_total_intercept = model.sma.SRR25_Dim.3_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.SRR25_Dim.3_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.SRR25_Dim.3_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.SRR25_Dim.3_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.SRR25_Dim.3_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.SRR25_Dim.3_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.SRR25_Dim.3_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.SRR25_Dim.3_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.SRR25_Dim.3_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.SRR25_Dim.3_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.SRR25_Dim.3_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.SRR25_Dim.3_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.SRR25_Dim.3_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.SRR25_Dim.3_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.SRR25_Dim.3_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.SRR25_Dim.3_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.SRR25_Dim.3_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.SRR25_Dim.3_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitDataPCA_SMA_touse %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, SRR25, Dim.3) %>%
  dplyr::mutate(Dim.3 = as.numeric(Dim.3), SRR25 = as.numeric(SRR25))
smaplt_total_data$predicted_SRR25 = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$Dim.3
smaplt_total_data$predicted_SRR25_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$Dim.3
smaplt_total_data$predicted_SRR25_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$Dim.3

smaplt_tree_data = traitDataPCA_SMA_touse_tree %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, SRR25, Dim.3) %>%
  dplyr::mutate(Dim.3 = as.numeric(Dim.3), SRR25 = as.numeric(SRR25))
smaplt_tree_data$predicted_SRR25 = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$Dim.3
smaplt_tree_data$predicted_SRR25_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$Dim.3
smaplt_tree_data$predicted_SRR25_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$Dim.3

smaplt_shrub_data = traitDataPCA_SMA_touse_shrub %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, SRR25, Dim.3) %>%
  dplyr::mutate(Dim.3 = as.numeric(Dim.3), SRR25 = as.numeric(SRR25))
smaplt_shrub_data$predicted_SRR25 = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$Dim.3
smaplt_shrub_data$predicted_SRR25_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$Dim.3
smaplt_shrub_data$predicted_SRR25_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$Dim.3

# Annotation, R2 and p-value
model.sma.SRR25_Dim.3_annotlabel = paste0(
  "Merged: R2=", round(as.numeric(model.sma.SRR25_Dim.3_total[["r2"]]), digits = 2), " p",
  ifelse(signif(as.numeric(model.sma.SRR25_Dim.3_total[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.SRR25_Dim.3_total[["pval"]]), digits = 2), "\n")),
  
  "Tree: R2=", round(as.numeric(model.sma.SRR25_Dim.3_tree[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.SRR25_Dim.3_tree[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.SRR25_Dim.3_tree[["pval"]]), digits = 2), "\n")),
  
  "Shrub: R2=", round(as.numeric(model.sma.SRR25_Dim.3_shrub[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.SRR25_Dim.3_shrub[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.SRR25_Dim.3_shrub[["pval"]]), digits = 2), "\n"))
)

smaplt.SRR25_Dim.3 = ggplot(smaplt_total_data, aes(x=Dim.3, y=SRR25)) +
  geom_point(aes(color=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$Dim.3, na.rm=T), max(smaplt_total_data$Dim.3, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$Dim.3, na.rm=T), max(smaplt_tree_data$Dim.3, na.rm=T)),
    color = "#F07673", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$Dim.3, na.rm=T), max(smaplt_shrub_data$Dim.3, na.rm=T)),
    color = "#7998AD", linewidth = 1.0
  ) +
  labs(x = "PC3 (RES+LES+RDMC+LPC)", y = "logSRR25", subtitle = "SRR25-PC3\n(mainly root collaboration)") + 
  theme_classic() + theme(legend.position = "none") +
  scale_color_manual(values = c("#EECA40", "#7998AD", "#F07673")) + 
  annotate("text", x = -2, y = 5.5, label = model.sma.SRR25_Dim.3_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.SRR25_Dim.3, filename = "smaplt.SRR25_PC3-RES-collab.pdf", width = 3, height = 2.9)
ggsave(plot = smaplt.SRR25_Dim.3, filename = "smaplt.SRR25_PC3-RES-collab.png", width = 3, height = 2.9, dpi=300)


# ------------------------------------------------------------------------------
# PART 4: Predicting RdarK25P using PC1 (RES conservation, RES+LES with RDMC and RPC) ----
# 
model.sma.Rdark25P_Dim.1_total = sma(traitDataPCA_SMA_touse$Rdark25P ~ traitDataPCA_SMA_touse$Dim.1)
model.sma.Rdark25P_Dim.1_tree = sma(traitDataPCA_SMA_touse_tree$Rdark25P ~ traitDataPCA_SMA_touse_tree$Dim.1)
model.sma.Rdark25P_Dim.1_shrub = sma(traitDataPCA_SMA_touse_shrub$Rdark25P ~ traitDataPCA_SMA_touse_shrub$Dim.1)

#plot(model.sma.Rdark25P_Dim.1)
summary(model.sma.Rdark25P_Dim.1_total)
summary(model.sma.Rdark25P_Dim.1_tree)
summary(model.sma.Rdark25P_Dim.1_shrub)

smaplt_total_intercept = model.sma.Rdark25P_Dim.1_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.Rdark25P_Dim.1_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.Rdark25P_Dim.1_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.Rdark25P_Dim.1_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.Rdark25P_Dim.1_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.Rdark25P_Dim.1_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.Rdark25P_Dim.1_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.Rdark25P_Dim.1_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.Rdark25P_Dim.1_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.Rdark25P_Dim.1_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.Rdark25P_Dim.1_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.Rdark25P_Dim.1_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.Rdark25P_Dim.1_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.Rdark25P_Dim.1_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.Rdark25P_Dim.1_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.Rdark25P_Dim.1_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.Rdark25P_Dim.1_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.Rdark25P_Dim.1_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitDataPCA_SMA_touse %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, Rdark25P, Dim.1) %>%
  dplyr::mutate(Dim.1 = as.numeric(Dim.1), Rdark25P = as.numeric(Rdark25P))
smaplt_total_data$predicted_Rdark25P = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$Dim.1
smaplt_total_data$predicted_Rdark25P_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$Dim.1
smaplt_total_data$predicted_Rdark25P_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$Dim.1

smaplt_tree_data = traitDataPCA_SMA_touse_tree %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, Rdark25P, Dim.1) %>%
  dplyr::mutate(Dim.1 = as.numeric(Dim.1), Rdark25P = as.numeric(Rdark25P))
smaplt_tree_data$predicted_Rdark25P = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$Dim.1
smaplt_tree_data$predicted_Rdark25P_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$Dim.1
smaplt_tree_data$predicted_Rdark25P_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$Dim.1

smaplt_shrub_data = traitDataPCA_SMA_touse_shrub %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, Rdark25P, Dim.1) %>%
  dplyr::mutate(Dim.1 = as.numeric(Dim.1), Rdark25P = as.numeric(Rdark25P))
smaplt_shrub_data$predicted_Rdark25P = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$Dim.1
smaplt_shrub_data$predicted_Rdark25P_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$Dim.1
smaplt_shrub_data$predicted_Rdark25P_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$Dim.1

# Annotation, R2 and p-value
model.sma.Rdark25P_Dim.1_annotlabel = paste0(
  "Merged: R2=", signif(as.numeric(model.sma.Rdark25P_Dim.1_total[["r2"]]), digits = 3), " p=", signif(as.numeric(model.sma.Rdark25P_Dim.1_total[["pval"]]), digits = 3), "\n",
  "Tree: R2=", signif(as.numeric(model.sma.Rdark25P_Dim.1_tree[["r2"]]), digits = 3), " p=", signif(as.numeric(model.sma.Rdark25P_Dim.1_tree[["pval"]]), digits = 3), "\n",
  "Shrub: R2=", signif(as.numeric(model.sma.Rdark25P_Dim.1_shrub[["r2"]]), digits = 3), " p=", signif(as.numeric(model.sma.Rdark25P_Dim.1_shrub[["pval"]]), digits = 3), "\n"
)

smaplt.Rdark25P_Dim.1 = ggplot(smaplt_total_data, aes(x=Dim.1, y=Rdark25P)) +
  geom_point(aes(color=GrowthForm), size=3, alpha=0.7) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$Dim.1, na.rm=T), max(smaplt_total_data$Dim.1, na.rm=T)),
    color = "darkblue", linewidth = 1.2
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$Dim.1, na.rm=T), max(smaplt_tree_data$Dim.1, na.rm=T)),
    color = "lightblue", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$Dim.1, na.rm=T), max(smaplt_shrub_data$Dim.1, na.rm=T)),
    color = "green", linewidth = 1.0
  ) +
  labs(x = "Dim.1", y = "logRdark25P") + 
  theme_clean() + theme(legend.position = "none") +
  annotate("text", x = -1.5, y = 1.15, label = model.sma.Rdark25P_Dim.1_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.Rdark25P_Dim.1, filename = "smaplt.Rdark25P_Dim.1.pdf", width = 4, height = 3.5)


# ------------------------------------------------------------------------------
# PART 5: Predicting Rdark25P using PC2 (LES, RES+LES with RDMC and RPC) ----

model.sma.Rdark25P_Dim.2_total = sma(traitDataPCA_SMA_touse$Rdark25P ~ traitDataPCA_SMA_touse$Dim.2)
model.sma.Rdark25P_Dim.2_tree = sma(traitDataPCA_SMA_touse_tree$Rdark25P ~ traitDataPCA_SMA_touse_tree$Dim.2)
model.sma.Rdark25P_Dim.2_shrub = sma(traitDataPCA_SMA_touse_shrub$Rdark25P ~ traitDataPCA_SMA_touse_shrub$Dim.2)

#plot(model.sma.Rdark25P_Dim.2)
summary(model.sma.Rdark25P_Dim.2_total)
summary(model.sma.Rdark25P_Dim.2_tree)
summary(model.sma.Rdark25P_Dim.2_shrub)

smaplt_total_intercept = model.sma.Rdark25P_Dim.2_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.Rdark25P_Dim.2_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.Rdark25P_Dim.2_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.Rdark25P_Dim.2_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.Rdark25P_Dim.2_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.Rdark25P_Dim.2_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.Rdark25P_Dim.2_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.Rdark25P_Dim.2_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.Rdark25P_Dim.2_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.Rdark25P_Dim.2_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.Rdark25P_Dim.2_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.Rdark25P_Dim.2_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.Rdark25P_Dim.2_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.Rdark25P_Dim.2_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.Rdark25P_Dim.2_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.Rdark25P_Dim.2_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.Rdark25P_Dim.2_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.Rdark25P_Dim.2_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitDataPCA_SMA_touse %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, Rdark25P, Dim.2) %>%
  dplyr::mutate(Dim.2 = as.numeric(Dim.2), Rdark25P = as.numeric(Rdark25P))
smaplt_total_data$predicted_Rdark25P = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$Dim.2
smaplt_total_data$predicted_Rdark25P_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$Dim.2
smaplt_total_data$predicted_Rdark25P_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$Dim.2

smaplt_tree_data = traitDataPCA_SMA_touse_tree %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, Rdark25P, Dim.2) %>%
  dplyr::mutate(Dim.2 = as.numeric(Dim.2), Rdark25P = as.numeric(Rdark25P))
smaplt_tree_data$predicted_Rdark25P = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$Dim.2
smaplt_tree_data$predicted_Rdark25P_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$Dim.2
smaplt_tree_data$predicted_Rdark25P_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$Dim.2

smaplt_shrub_data = traitDataPCA_SMA_touse_shrub %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, Rdark25P, Dim.2) %>%
  dplyr::mutate(Dim.2 = as.numeric(Dim.2), Rdark25P = as.numeric(Rdark25P))
smaplt_shrub_data$predicted_Rdark25P = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$Dim.2
smaplt_shrub_data$predicted_Rdark25P_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$Dim.2
smaplt_shrub_data$predicted_Rdark25P_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$Dim.2

# Annotation, R2 and p-value
model.sma.Rdark25P_Dim.2_annotlabel = paste0(
  "Merged: R2=", signif(as.numeric(model.sma.Rdark25P_Dim.2_total[["r2"]]), digits = 3), " p=", signif(as.numeric(model.sma.Rdark25P_Dim.2_total[["pval"]]), digits = 3), "\n",
  "Tree: R2=", signif(as.numeric(model.sma.Rdark25P_Dim.2_tree[["r2"]]), digits = 3), " p=", signif(as.numeric(model.sma.Rdark25P_Dim.2_tree[["pval"]]), digits = 3), "\n",
  "Shrub: R2=", signif(as.numeric(model.sma.Rdark25P_Dim.2_shrub[["r2"]]), digits = 3), " p=", signif(as.numeric(model.sma.Rdark25P_Dim.2_shrub[["pval"]]), digits = 3), "\n"
)

smaplt.Rdark25P_Dim.2 = ggplot(smaplt_total_data, aes(x=Dim.2, y=Rdark25P)) +
  geom_point(aes(color=GrowthForm), size=3, alpha=0.7) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$Dim.2, na.rm=T), max(smaplt_total_data$Dim.2, na.rm=T)),
    color = "darkblue", linewidth = 1.2
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$Dim.2, na.rm=T), max(smaplt_tree_data$Dim.2, na.rm=T)),
    color = "lightblue", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$Dim.2, na.rm=T), max(smaplt_shrub_data$Dim.2, na.rm=T)),
    color = "green", linewidth = 1.0
  ) +
  labs(x = "Dim.2", y = "logRdark25P") + 
  theme_clean() + theme(legend.position = "none") +
  annotate("text", x = -3, y = 0.15, label = model.sma.Rdark25P_Dim.2_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.Rdark25P_Dim.2, filename = "smaplt.Rdark25P_Dim.2.pdf", width = 4, height = 3.5)


# ------------------------------------------------------------------------------
# PART 6: Predicting Rdark25 using PC3 (RES collaboration, RES+LES with RDMC and RPC) ----

model.sma.Rdark25P_Dim.3_total = sma(traitDataPCA_SMA_touse$Rdark25P ~ traitDataPCA_SMA_touse$Dim.3)
model.sma.Rdark25P_Dim.3_tree = sma(traitDataPCA_SMA_touse_tree$Rdark25P ~ traitDataPCA_SMA_touse_tree$Dim.3)
model.sma.Rdark25P_Dim.3_shrub = sma(traitDataPCA_SMA_touse_shrub$Rdark25P ~ traitDataPCA_SMA_touse_shrub$Dim.3)

#plot(model.sma.Rdark25P_Dim.3)
summary(model.sma.Rdark25P_Dim.3_total)
summary(model.sma.Rdark25P_Dim.3_tree)
summary(model.sma.Rdark25P_Dim.3_shrub)

smaplt_total_intercept = model.sma.Rdark25P_Dim.3_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.Rdark25P_Dim.3_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.Rdark25P_Dim.3_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.Rdark25P_Dim.3_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.Rdark25P_Dim.3_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.Rdark25P_Dim.3_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.Rdark25P_Dim.3_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.Rdark25P_Dim.3_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.Rdark25P_Dim.3_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.Rdark25P_Dim.3_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.Rdark25P_Dim.3_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.Rdark25P_Dim.3_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.Rdark25P_Dim.3_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.Rdark25P_Dim.3_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.Rdark25P_Dim.3_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.Rdark25P_Dim.3_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.Rdark25P_Dim.3_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.Rdark25P_Dim.3_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitDataPCA_SMA_touse %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, Rdark25P, Dim.3) %>%
  dplyr::mutate(Dim.3 = as.numeric(Dim.3), Rdark25P = as.numeric(Rdark25P))
smaplt_total_data$predicted_Rdark25P = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$Dim.3
smaplt_total_data$predicted_Rdark25P_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$Dim.3
smaplt_total_data$predicted_Rdark25P_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$Dim.3

smaplt_tree_data = traitDataPCA_SMA_touse_tree %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, Rdark25P, Dim.3) %>%
  dplyr::mutate(Dim.3 = as.numeric(Dim.3), Rdark25P = as.numeric(Rdark25P))
smaplt_tree_data$predicted_Rdark25P = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$Dim.3
smaplt_tree_data$predicted_Rdark25P_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$Dim.3
smaplt_tree_data$predicted_Rdark25P_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$Dim.3

smaplt_shrub_data = traitDataPCA_SMA_touse_shrub %>% dplyr::select(PCAIdentifier_spgf, GrowthForm, Rdark25P, Dim.3) %>%
  dplyr::mutate(Dim.3 = as.numeric(Dim.3), Rdark25P = as.numeric(Rdark25P))
smaplt_shrub_data$predicted_Rdark25P = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$Dim.3
smaplt_shrub_data$predicted_Rdark25P_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$Dim.3
smaplt_shrub_data$predicted_Rdark25P_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$Dim.3

# Annotation, R2 and p-value
model.sma.Rdark25P_Dim.3_annotlabel = paste0(
  "Merged: R2=", signif(as.numeric(model.sma.Rdark25P_Dim.3_total[["r2"]]), digits = 3), " p=", signif(as.numeric(model.sma.Rdark25P_Dim.3_total[["pval"]]), digits = 3), "\n",
  "Tree: R2=", signif(as.numeric(model.sma.Rdark25P_Dim.3_tree[["r2"]]), digits = 3), " p=", signif(as.numeric(model.sma.Rdark25P_Dim.3_tree[["pval"]]), digits = 3), "\n",
  "Shrub: R2=", signif(as.numeric(model.sma.Rdark25P_Dim.3_shrub[["r2"]]), digits = 3), " p=", signif(as.numeric(model.sma.Rdark25P_Dim.3_shrub[["pval"]]), digits = 3), "\n"
)

smaplt.Rdark25P_Dim.3 = ggplot(smaplt_total_data, aes(x=Dim.3, y=Rdark25P)) +
  geom_point(aes(color=GrowthForm), size=3, alpha=0.7) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$Dim.3, na.rm=T), max(smaplt_total_data$Dim.3, na.rm=T)),
    color = "darkblue", linewidth = 1.2
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$Dim.3, na.rm=T), max(smaplt_tree_data$Dim.3, na.rm=T)),
    color = "lightblue", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$Dim.3, na.rm=T), max(smaplt_shrub_data$Dim.3, na.rm=T)),
    color = "green", linewidth = 1.0
  ) +
  labs(x = "Dim.3", y = "logRdark25P") + 
  theme_clean() + theme(legend.position = "none") +
  annotate("text", x = -1.7, y = 1.3, label = model.sma.Rdark25P_Dim.3_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.Rdark25P_Dim.3, filename = "smaplt.Rdark25P_Dim.3.pdf", width = 4, height = 3.5)

# model <- lm(RDMC ~ Dim.2, data = data_for_regression)
# summary(model)
# 
# # Add predicted values to the dataframe
# data_for_regression$predicted_RTD <- predict(model)
# 
# # Plot actual vs predicted values for RTD
# ggplot(data_for_regression, aes(x = predicted_RTD, y = RTD)) +
#   geom_point() +
#   geom_smooth(method = "lm", se = FALSE, color = "blue") +
#   labs(x = "Predicted RTD", y = "Actual RTD", title = "Predicted vs Actual RTD using PCs") +
#   theme_minimal()
