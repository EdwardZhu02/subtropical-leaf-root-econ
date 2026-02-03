rm(list=ls())
library(tidyverse)
library(dplyr)
library(smatr) # SMA regression
library(ggthemes)
library(patchwork) # plot merging
library(cowplot) # plot merging

load("species-level-code/traitDataFujian-spavg-phylo-step2.RData")

### SET DATA TO USE ###
traitData_touse = traitDataIndv_spgfavg_log

# divide data based on growth forms for separate regression
traitData_touse_tree = traitData_touse %>% dplyr::filter(GrowthForm == "tree")
traitData_touse_shrub = traitData_touse %>% dplyr::filter(GrowthForm == "shrub")
traitData_touse_liana = traitData_touse %>% dplyr::filter(GrowthForm == "liana") # may be useless


# ------------------------------------------------------------------------------
# PART 1: LES (LMA-LNC)
# ------------------------------------------------------------------------------
model.sma.LNC_LMA_total = sma(traitData_touse$LNC ~ traitData_touse$LMA)
model.sma.LNC_LMA_tree = sma(traitData_touse_tree$LNC ~ traitData_touse_tree$LMA)
model.sma.LNC_LMA_shrub = sma(traitData_touse_shrub$LNC ~ traitData_touse_shrub$LMA)

#plot(model.sma.LNC_LMA)
summary(model.sma.LNC_LMA_total)
summary(model.sma.LNC_LMA_tree)
summary(model.sma.LNC_LMA_shrub)

smaplt_total_intercept = model.sma.LNC_LMA_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.LNC_LMA_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.LNC_LMA_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.LNC_LMA_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.LNC_LMA_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.LNC_LMA_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.LNC_LMA_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.LNC_LMA_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.LNC_LMA_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.LNC_LMA_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.LNC_LMA_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.LNC_LMA_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.LNC_LMA_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.LNC_LMA_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.LNC_LMA_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.LNC_LMA_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.LNC_LMA_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.LNC_LMA_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, LNC, LMA) %>%
  dplyr::mutate(LMA = as.numeric(LMA), LNC = as.numeric(LNC))
smaplt_total_data$predicted_LNC = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$LMA
smaplt_total_data$predicted_LNC_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$LMA
smaplt_total_data$predicted_LNC_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$LMA

smaplt_tree_data = traitData_touse_tree %>% dplyr::select(speciesFullName, GrowthForm, LNC, LMA) %>%
  dplyr::mutate(LMA = as.numeric(LMA), LNC = as.numeric(LNC))
smaplt_tree_data$predicted_LNC = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$LMA
smaplt_tree_data$predicted_LNC_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$LMA
smaplt_tree_data$predicted_LNC_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$LMA

smaplt_shrub_data = traitData_touse_shrub %>% dplyr::select(speciesFullName, GrowthForm, LNC, LMA) %>%
  dplyr::mutate(LMA = as.numeric(LMA), LNC = as.numeric(LNC))
smaplt_shrub_data$predicted_LNC = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$LMA
smaplt_shrub_data$predicted_LNC_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$LMA
smaplt_shrub_data$predicted_LNC_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$LMA

# Annotation, R2 and p-value
model.sma.LNC_LMA_annotlabel = paste0(
  "Merged: R2=", round(as.numeric(model.sma.LNC_LMA_total[["r2"]]), digits = 2), " p",
  ifelse(signif(as.numeric(model.sma.LNC_LMA_total[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.LNC_LMA_total[["pval"]]), digits = 2), "\n")),
  
  "Tree: R2=", round(as.numeric(model.sma.LNC_LMA_tree[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.LNC_LMA_tree[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.LNC_LMA_tree[["pval"]]), digits = 2), "\n")),
  
  "Shrub: R2=", round(as.numeric(model.sma.LNC_LMA_shrub[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.LNC_LMA_shrub[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.LNC_LMA_shrub[["pval"]]), digits = 2), "\n"))
  )

smaplt.LNC_LMA = ggplot(smaplt_total_data, aes(x=LMA, y=LNC)) +
  geom_point(aes(color=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$LMA, na.rm=T), max(smaplt_total_data$LMA, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$LMA, na.rm=T), max(smaplt_tree_data$LMA, na.rm=T)),
    color = "#F07673", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$LMA, na.rm=T), max(smaplt_shrub_data$LMA, na.rm=T)),
    color = "#7998AD", linewidth = 1.0
    ) +
  labs(subtitle = "LES fast-slow (LNC-LMA)", x = "lnLMA", y = "lnLNC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#EECA40", "#7998AD", "#F07673")) + 
  annotate("text", x = 0.008, y = 1.18, label = model.sma.LNC_LMA_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.LNC_LMA, filename = "species-level-code/smaplt.LNC_LMA.pdf", width = 3, height = 2.9)

# ------------------------------------------------------------------------------
# PART 2: RES collaboration (RD-SRL)
# ------------------------------------------------------------------------------
model.sma.RD_SRL_total = sma(traitData_touse$RD ~ traitData_touse$SRL)
model.sma.RD_SRL_tree = sma(traitData_touse_tree$RD ~ traitData_touse_tree$SRL)
model.sma.RD_SRL_shrub = sma(traitData_touse_shrub$RD ~ traitData_touse_shrub$SRL)

#plot(model.sma.RD_SRL)
summary(model.sma.RD_SRL_total)
summary(model.sma.RD_SRL_tree)
summary(model.sma.RD_SRL_shrub)

smaplt_total_intercept = model.sma.RD_SRL_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.RD_SRL_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.RD_SRL_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.RD_SRL_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.RD_SRL_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.RD_SRL_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.RD_SRL_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.RD_SRL_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.RD_SRL_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.RD_SRL_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.RD_SRL_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.RD_SRL_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.RD_SRL_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.RD_SRL_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.RD_SRL_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.RD_SRL_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.RD_SRL_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.RD_SRL_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, RD, SRL) %>%
  dplyr::mutate(SRL = as.numeric(SRL), RD = as.numeric(RD))
smaplt_total_data$predicted_RD = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$SRL
smaplt_total_data$predicted_RD_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$SRL
smaplt_total_data$predicted_RD_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$SRL

smaplt_tree_data = traitData_touse_tree %>% dplyr::select(speciesFullName, GrowthForm, RD, SRL) %>%
  dplyr::mutate(SRL = as.numeric(SRL), RD = as.numeric(RD))
smaplt_tree_data$predicted_RD = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$SRL
smaplt_tree_data$predicted_RD_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$SRL
smaplt_tree_data$predicted_RD_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$SRL

smaplt_shrub_data = traitData_touse_shrub %>% dplyr::select(speciesFullName, GrowthForm, RD, SRL) %>%
  dplyr::mutate(SRL = as.numeric(SRL), RD = as.numeric(RD))
smaplt_shrub_data$predicted_RD = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$SRL
smaplt_shrub_data$predicted_RD_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$SRL
smaplt_shrub_data$predicted_RD_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$SRL

# Annotation, R2 and p-value
model.sma.RD_SRL_annotlabel = paste0(
  "Merged: R2=", round(as.numeric(model.sma.RD_SRL_total[["r2"]]), digits = 2), " p",
  ifelse(signif(as.numeric(model.sma.RD_SRL_total[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.RD_SRL_total[["pval"]]), digits = 2), "\n")),
  
  "Tree: R2=", round(as.numeric(model.sma.RD_SRL_tree[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.RD_SRL_tree[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.RD_SRL_tree[["pval"]]), digits = 2), "\n")),
  
  "Shrub: R2=", round(as.numeric(model.sma.RD_SRL_shrub[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.RD_SRL_shrub[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.RD_SRL_shrub[["pval"]]), digits = 2), "\n"))
)

smaplt.RD_SRL = ggplot(smaplt_total_data, aes(x=SRL, y=RD)) +
  geom_point(aes(color=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$SRL, na.rm=T), max(smaplt_total_data$SRL, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$SRL, na.rm=T), max(smaplt_tree_data$SRL, na.rm=T)),
    color = "#F07673", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$SRL, na.rm=T), max(smaplt_shrub_data$SRL, na.rm=T)),
    color = "#7998AD", linewidth = 1.0
  ) +
  labs(subtitle = "RES collaboration (RD-SRL)", x = "lnSRL", y = "lnRD") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#EECA40", "#7998AD", "#F07673")) + 
  annotate("text", x = 3.6, y = 0.53, label = model.sma.RD_SRL_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.RD_SRL, filename = "species-level-code/smaplt.RD_SRL.pdf", width = 3, height = 2.9)

# ------------------------------------------------------------------------------
# PART 3: RES conservation (RNC-RTD)
# ------------------------------------------------------------------------------
model.sma.RNC_RTD_total = sma(traitData_touse$RNC ~ traitData_touse$RTD)
model.sma.RNC_RTD_tree = sma(traitData_touse_tree$RNC ~ traitData_touse_tree$RTD)
model.sma.RNC_RTD_shrub = sma(traitData_touse_shrub$RNC ~ traitData_touse_shrub$RTD)

#plot(model.sma.RNC_RTD)
summary(model.sma.RNC_RTD_total)
summary(model.sma.RNC_RTD_tree)
summary(model.sma.RNC_RTD_shrub)

smaplt_total_intercept = model.sma.RNC_RTD_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.RNC_RTD_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.RNC_RTD_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.RNC_RTD_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.RNC_RTD_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.RNC_RTD_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.RNC_RTD_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.RNC_RTD_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.RNC_RTD_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.RNC_RTD_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.RNC_RTD_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.RNC_RTD_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.RNC_RTD_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.RNC_RTD_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.RNC_RTD_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.RNC_RTD_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.RNC_RTD_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.RNC_RTD_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, RNC, RTD) %>%
  dplyr::mutate(RTD = as.numeric(RTD), RNC = as.numeric(RNC))
smaplt_total_data$predicted_RNC = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$RTD
smaplt_total_data$predicted_RNC_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$RTD
smaplt_total_data$predicted_RNC_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$RTD

smaplt_tree_data = traitData_touse_tree %>% dplyr::select(speciesFullName, GrowthForm, RNC, RTD) %>%
  dplyr::mutate(RTD = as.numeric(RTD), RNC = as.numeric(RNC))
smaplt_tree_data$predicted_RNC = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$RTD
smaplt_tree_data$predicted_RNC_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$RTD
smaplt_tree_data$predicted_RNC_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$RTD

smaplt_shrub_data = traitData_touse_shrub %>% dplyr::select(speciesFullName, GrowthForm, RNC, RTD) %>%
  dplyr::mutate(RTD = as.numeric(RTD), RNC = as.numeric(RNC))
smaplt_shrub_data$predicted_RNC = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$RTD
smaplt_shrub_data$predicted_RNC_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$RTD
smaplt_shrub_data$predicted_RNC_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$RTD

# Annotation, R2 and p-value
model.sma.RNC_RTD_annotlabel = paste0(
  "Merged: R2=", round(as.numeric(model.sma.RNC_RTD_total[["r2"]]), digits = 2), " p",
  ifelse(signif(as.numeric(model.sma.RNC_RTD_total[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.RNC_RTD_total[["pval"]]), digits = 2), "\n")),
  
  "Tree: R2=", round(as.numeric(model.sma.RNC_RTD_tree[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.RNC_RTD_tree[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.RNC_RTD_tree[["pval"]]), digits = 2), "\n")),
  
  "Shrub: R2=", round(as.numeric(model.sma.RNC_RTD_shrub[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.RNC_RTD_shrub[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.RNC_RTD_shrub[["pval"]]), digits = 2), "\n"))
)

smaplt.RNC_RTD = ggplot(smaplt_total_data, aes(x=RTD, y=RNC)) +
  geom_point(aes(color=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$RTD, na.rm=T), max(smaplt_total_data$RTD, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$RTD, na.rm=T), max(smaplt_tree_data$RTD, na.rm=T)),
    color = "#F07673", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$RTD, na.rm=T), max(smaplt_shrub_data$RTD, na.rm=T)),
    color = "#7998AD", linewidth = 1.0
  ) +
  labs(subtitle = "RES conservation (RNC-RTD)", x = "lnRTD", y = "lnRNC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#EECA40", "#7998AD", "#F07673")) + 
  annotate("text", x = 0.165, y = 1.48, label = model.sma.RNC_RTD_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.RNC_RTD, filename = "species-level-code/smaplt.RNC_RTD.pdf", width = 3, height = 2.9)


# ------------------------------------------------------------------------------
# PART 4: RDMC-RTD
# ------------------------------------------------------------------------------
model.sma.RDMC_RTD_total = sma(traitData_touse$RDMC ~ traitData_touse$RTD)
model.sma.RDMC_RTD_tree = sma(traitData_touse_tree$RDMC ~ traitData_touse_tree$RTD)
model.sma.RDMC_RTD_shrub = sma(traitData_touse_shrub$RDMC ~ traitData_touse_shrub$RTD)

#plot(model.sma.RDMC_RTD)
summary(model.sma.RDMC_RTD_total)
summary(model.sma.RDMC_RTD_tree)
summary(model.sma.RDMC_RTD_shrub)

smaplt_total_intercept = model.sma.RDMC_RTD_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.RDMC_RTD_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.RDMC_RTD_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.RDMC_RTD_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.RDMC_RTD_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.RDMC_RTD_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.RDMC_RTD_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.RDMC_RTD_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.RDMC_RTD_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.RDMC_RTD_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.RDMC_RTD_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.RDMC_RTD_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.RDMC_RTD_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.RDMC_RTD_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.RDMC_RTD_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.RDMC_RTD_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.RDMC_RTD_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.RDMC_RTD_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, RDMC, RTD) %>%
  dplyr::mutate(RTD = as.numeric(RTD), RDMC = as.numeric(RDMC))
smaplt_total_data$predicted_RDMC = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$RTD
smaplt_total_data$predicted_RDMC_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$RTD
smaplt_total_data$predicted_RDMC_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$RTD

smaplt_tree_data = traitData_touse_tree %>% dplyr::select(speciesFullName, GrowthForm, RDMC, RTD) %>%
  dplyr::mutate(RTD = as.numeric(RTD), RDMC = as.numeric(RDMC))
smaplt_tree_data$predicted_RDMC = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$RTD
smaplt_tree_data$predicted_RDMC_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$RTD
smaplt_tree_data$predicted_RDMC_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$RTD

smaplt_shrub_data = traitData_touse_shrub %>% dplyr::select(speciesFullName, GrowthForm, RDMC, RTD) %>%
  dplyr::mutate(RTD = as.numeric(RTD), RDMC = as.numeric(RDMC))
smaplt_shrub_data$predicted_RDMC = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$RTD
smaplt_shrub_data$predicted_RDMC_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$RTD
smaplt_shrub_data$predicted_RDMC_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$RTD

# Annotation, R2 and p-value
model.sma.RDMC_RTD_annotlabel = paste0(
  "Merged: R2=", round(as.numeric(model.sma.RDMC_RTD_total[["r2"]]), digits = 2), " p",
  ifelse(signif(as.numeric(model.sma.RDMC_RTD_total[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.RDMC_RTD_total[["pval"]]), digits = 2), "\n")),
  
  "Tree: R2=", round(as.numeric(model.sma.RDMC_RTD_tree[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.RDMC_RTD_tree[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.RDMC_RTD_tree[["pval"]]), digits = 2), "\n")),
  
  "Shrub: R2=", round(as.numeric(model.sma.RDMC_RTD_shrub[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.RDMC_RTD_shrub[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.RDMC_RTD_shrub[["pval"]]), digits = 2), "\n"))
)

smaplt.RDMC_RTD = ggplot(smaplt_total_data, aes(x=RTD, y=RDMC)) +
  geom_point(aes(color=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$RTD, na.rm=T), max(smaplt_total_data$RTD, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$RTD, na.rm=T), max(smaplt_tree_data$RTD, na.rm=T)),
    color = "#F07673", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$RTD, na.rm=T), max(smaplt_shrub_data$RTD, na.rm=T)),
    color = "#7998AD", linewidth = 1.0
  ) +
  labs(subtitle = "RES conservation (RDMC-RTD)", x = "lnRTD", y = "lnRDMC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#EECA40", "#7998AD", "#F07673")) + 
  annotate("text", x = 0.05, y = 0.42, label = model.sma.RDMC_RTD_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.RDMC_RTD, filename = "species-level-code/smaplt.RDMC_RTD.pdf", width = 3, height = 3)


# ------------------------------------------------------------------------------
# PART 5: RDMC-RNC
# ------------------------------------------------------------------------------
model.sma.RDMC_RNC_total = sma(traitData_touse$RDMC ~ traitData_touse$RNC)
model.sma.RDMC_RNC_tree = sma(traitData_touse_tree$RDMC ~ traitData_touse_tree$RNC)
model.sma.RDMC_RNC_shrub = sma(traitData_touse_shrub$RDMC ~ traitData_touse_shrub$RNC)

#plot(model.sma.RDMC_RNC)
summary(model.sma.RDMC_RNC_total)
summary(model.sma.RDMC_RNC_tree)
summary(model.sma.RDMC_RNC_shrub)

smaplt_total_intercept = model.sma.RDMC_RNC_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.RDMC_RNC_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.RDMC_RNC_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.RDMC_RNC_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.RDMC_RNC_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.RDMC_RNC_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.RDMC_RNC_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.RDMC_RNC_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.RDMC_RNC_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.RDMC_RNC_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.RDMC_RNC_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.RDMC_RNC_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.RDMC_RNC_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.RDMC_RNC_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.RDMC_RNC_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.RDMC_RNC_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.RDMC_RNC_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.RDMC_RNC_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, RDMC, RNC) %>%
  dplyr::mutate(RNC = as.numeric(RNC), RDMC = as.numeric(RDMC))
smaplt_total_data$predicted_RDMC = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$RNC
smaplt_total_data$predicted_RDMC_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$RNC
smaplt_total_data$predicted_RDMC_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$RNC

smaplt_tree_data = traitData_touse_tree %>% dplyr::select(speciesFullName, GrowthForm, RDMC, RNC) %>%
  dplyr::mutate(RNC = as.numeric(RNC), RDMC = as.numeric(RDMC))
smaplt_tree_data$predicted_RDMC = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$RNC
smaplt_tree_data$predicted_RDMC_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$RNC
smaplt_tree_data$predicted_RDMC_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$RNC

smaplt_shrub_data = traitData_touse_shrub %>% dplyr::select(speciesFullName, GrowthForm, RDMC, RNC) %>%
  dplyr::mutate(RNC = as.numeric(RNC), RDMC = as.numeric(RDMC))
smaplt_shrub_data$predicted_RDMC = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$RNC
smaplt_shrub_data$predicted_RDMC_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$RNC
smaplt_shrub_data$predicted_RDMC_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$RNC

# Annotation, R2 and p-value
model.sma.RDMC_RNC_annotlabel = paste0(
  "Merged: R2=", round(as.numeric(model.sma.RDMC_RNC_total[["r2"]]), digits = 2), " p",
  ifelse(signif(as.numeric(model.sma.RDMC_RNC_total[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.RDMC_RNC_total[["pval"]]), digits = 2), "\n")),
  
  "Tree: R2=", round(as.numeric(model.sma.RDMC_RNC_tree[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.RDMC_RNC_tree[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.RDMC_RNC_tree[["pval"]]), digits = 2), "\n")),
  
  "Shrub: R2=", round(as.numeric(model.sma.RDMC_RNC_shrub[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.RDMC_RNC_shrub[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.RDMC_RNC_shrub[["pval"]]), digits = 2), "\n"))
)

smaplt.RDMC_RNC = ggplot(smaplt_total_data, aes(x=RNC, y=RDMC)) +
  geom_point(aes(color=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$RNC, na.rm=T), max(smaplt_total_data$RNC, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$RNC, na.rm=T), max(smaplt_tree_data$RNC, na.rm=T)),
    color = "#F07673", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$RNC, na.rm=T), max(smaplt_shrub_data$RNC, na.rm=T)),
    color = "#7998AD", linewidth = 1.0
  ) +
  labs(subtitle = "RES conservation (RDMC-RNC)", x = "lnRNC", y = "lnRDMC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#EECA40", "#7998AD", "#F07673")) + 
  annotate("text", x = 0.87, y = 0.335, label = model.sma.RDMC_RNC_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.RDMC_RNC, filename = "species-level-code/smaplt.RDMC_RNC.pdf", width = 3, height = 3)


# ------------------------------------------------------------------------------
# PART 6: LES (LPC-LNC)
# ------------------------------------------------------------------------------
model.sma.LNC_LPC_total = sma(traitData_touse$LNC ~ traitData_touse$LPC)
model.sma.LNC_LPC_tree = sma(traitData_touse_tree$LNC ~ traitData_touse_tree$LPC)
model.sma.LNC_LPC_shrub = sma(traitData_touse_shrub$LNC ~ traitData_touse_shrub$LPC)

#plot(model.sma.LNC_LPC)
summary(model.sma.LNC_LPC_total)
summary(model.sma.LNC_LPC_tree)
summary(model.sma.LNC_LPC_shrub)

smaplt_total_intercept = model.sma.LNC_LPC_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.LNC_LPC_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.LNC_LPC_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.LNC_LPC_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.LNC_LPC_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.LNC_LPC_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.LNC_LPC_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.LNC_LPC_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.LNC_LPC_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.LNC_LPC_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.LNC_LPC_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.LNC_LPC_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.LNC_LPC_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.LNC_LPC_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.LNC_LPC_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.LNC_LPC_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.LNC_LPC_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.LNC_LPC_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, LNC, LPC) %>%
  dplyr::mutate(LPC = as.numeric(LPC), LNC = as.numeric(LNC))
smaplt_total_data$predicted_LNC = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$LPC
smaplt_total_data$predicted_LNC_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$LPC
smaplt_total_data$predicted_LNC_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$LPC

smaplt_tree_data = traitData_touse_tree %>% dplyr::select(speciesFullName, GrowthForm, LNC, LPC) %>%
  dplyr::mutate(LPC = as.numeric(LPC), LNC = as.numeric(LNC))
smaplt_tree_data$predicted_LNC = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$LPC
smaplt_tree_data$predicted_LNC_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$LPC
smaplt_tree_data$predicted_LNC_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$LPC

smaplt_shrub_data = traitData_touse_shrub %>% dplyr::select(speciesFullName, GrowthForm, LNC, LPC) %>%
  dplyr::mutate(LPC = as.numeric(LPC), LNC = as.numeric(LNC))
smaplt_shrub_data$predicted_LNC = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$LPC
smaplt_shrub_data$predicted_LNC_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$LPC
smaplt_shrub_data$predicted_LNC_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$LPC

# Annotation, R2 and p-value
model.sma.LNC_LPC_annotlabel = paste0(
  "Merged: R2=", round(as.numeric(model.sma.LNC_LPC_total[["r2"]]), digits = 2), " p",
  ifelse(signif(as.numeric(model.sma.LNC_LPC_total[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.LNC_LPC_total[["pval"]]), digits = 2), "\n")),
  
  "Tree: R2=", round(as.numeric(model.sma.LNC_LPC_tree[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.LNC_LPC_tree[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.LNC_LPC_tree[["pval"]]), digits = 2), "\n")),
  
  "Shrub: R2=", round(as.numeric(model.sma.LNC_LPC_shrub[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.LNC_LPC_shrub[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.LNC_LPC_shrub[["pval"]]), digits = 2), "\n"))
)

smaplt.LNC_LPC = ggplot(smaplt_total_data, aes(x=LPC, y=LNC)) +
  geom_point(aes(color=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$LPC, na.rm=T), max(smaplt_total_data$LPC, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$LPC, na.rm=T), max(smaplt_tree_data$LPC, na.rm=T)),
    color = "#F07673", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$LPC, na.rm=T), max(smaplt_shrub_data$LPC, na.rm=T)),
    color = "#7998AD", linewidth = 1.0
  ) +
  labs(subtitle = "LES fast-slow (LNC-LPC)", x = "lnLPC", y = "lnLNC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#EECA40", "#7998AD", "#F07673")) + 
  annotate("text", x = 0.37, y = 1.35, label = model.sma.LNC_LPC_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.LNC_LPC, filename = "species-level-code/smaplt.LNC_LPC.pdf", width = 3, height = 3)


# ------------------------------------------------------------------------------
# PART 7: LES (LPC-LMA)
# ------------------------------------------------------------------------------
model.sma.LMA_LPC_total = sma(traitData_touse$LMA ~ traitData_touse$LPC)
model.sma.LMA_LPC_tree = sma(traitData_touse_tree$LMA ~ traitData_touse_tree$LPC)
model.sma.LMA_LPC_shrub = sma(traitData_touse_shrub$LMA ~ traitData_touse_shrub$LPC)

#plot(model.sma.LMA_LPC)
summary(model.sma.LMA_LPC_total)
summary(model.sma.LMA_LPC_tree)
summary(model.sma.LMA_LPC_shrub)

smaplt_total_intercept = model.sma.LMA_LPC_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.LMA_LPC_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.LMA_LPC_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.LMA_LPC_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.LMA_LPC_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.LMA_LPC_total[["coef"]][[1]][2,3]

smaplt_tree_intercept = model.sma.LMA_LPC_tree[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_tree_intercept_lower = model.sma.LMA_LPC_tree[["coef"]][[1]][1,2]
smaplt_tree_intercept_upper = model.sma.LMA_LPC_tree[["coef"]][[1]][1,3] 
smaplt_tree_slope = model.sma.LMA_LPC_tree[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_tree_slope_lower = model.sma.LMA_LPC_tree[["coef"]][[1]][2,2]
smaplt_tree_slope_upper = model.sma.LMA_LPC_tree[["coef"]][[1]][2,3]

smaplt_shrub_intercept = model.sma.LMA_LPC_shrub[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_shrub_intercept_lower = model.sma.LMA_LPC_shrub[["coef"]][[1]][1,2]
smaplt_shrub_intercept_upper = model.sma.LMA_LPC_shrub[["coef"]][[1]][1,3] 
smaplt_shrub_slope = model.sma.LMA_LPC_shrub[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_shrub_slope_lower = model.sma.LMA_LPC_shrub[["coef"]][[1]][2,2]
smaplt_shrub_slope_upper = model.sma.LMA_LPC_shrub[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, LMA, LPC) %>%
  dplyr::mutate(LPC = as.numeric(LPC), LMA = as.numeric(LMA))
smaplt_total_data$predicted_LMA = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$LPC
smaplt_total_data$predicted_LMA_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$LPC
smaplt_total_data$predicted_LMA_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$LPC

smaplt_tree_data = traitData_touse_tree %>% dplyr::select(speciesFullName, GrowthForm, LMA, LPC) %>%
  dplyr::mutate(LPC = as.numeric(LPC), LMA = as.numeric(LMA))
smaplt_tree_data$predicted_LMA = smaplt_tree_intercept + smaplt_tree_slope * smaplt_tree_data$LPC
smaplt_tree_data$predicted_LMA_upper_tree = smaplt_tree_intercept_upper + smaplt_tree_slope_upper * smaplt_tree_data$LPC
smaplt_tree_data$predicted_LMA_lower_tree = smaplt_tree_intercept_lower + smaplt_tree_slope_lower * smaplt_tree_data$LPC

smaplt_shrub_data = traitData_touse_shrub %>% dplyr::select(speciesFullName, GrowthForm, LMA, LPC) %>%
  dplyr::mutate(LPC = as.numeric(LPC), LMA = as.numeric(LMA))
smaplt_shrub_data$predicted_LMA = smaplt_shrub_intercept + smaplt_shrub_slope * smaplt_shrub_data$LPC
smaplt_shrub_data$predicted_LMA_upper_shrub = smaplt_shrub_intercept_upper + smaplt_shrub_slope_upper * smaplt_shrub_data$LPC
smaplt_shrub_data$predicted_LMA_lower_shrub = smaplt_shrub_intercept_lower + smaplt_shrub_slope_lower * smaplt_shrub_data$LPC

# Annotation, R2 and p-value
model.sma.LMA_LPC_annotlabel = paste0(
  "Merged: R2=", round(as.numeric(model.sma.LMA_LPC_total[["r2"]]), digits = 2), " p",
  ifelse(signif(as.numeric(model.sma.LMA_LPC_total[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.LMA_LPC_total[["pval"]]), digits = 2), "\n")),
  
  "Tree: R2=", round(as.numeric(model.sma.LMA_LPC_tree[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.LMA_LPC_tree[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.LMA_LPC_tree[["pval"]]), digits = 2), "\n")),
  
  "Shrub: R2=", round(as.numeric(model.sma.LMA_LPC_shrub[["r2"]]), digits = 2), " p", 
  ifelse(signif(as.numeric(model.sma.LMA_LPC_shrub[["pval"]]), digits = 2)<0.01, "<0.01\n", paste0("=",signif(as.numeric(model.sma.LMA_LPC_shrub[["pval"]]), digits = 2), "\n"))
)

smaplt.LMA_LPC = ggplot(smaplt_total_data, aes(x=LPC, y=LMA)) +
  geom_point(aes(color=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$LPC, na.rm=T), max(smaplt_total_data$LPC, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # tree
    fun = function(x) smaplt_tree_slope * x + smaplt_tree_intercept, 
    xlim = c(min(smaplt_tree_data$LPC, na.rm=T), max(smaplt_tree_data$LPC, na.rm=T)),
    color = "#F07673", linewidth = 1.0
  ) +
  stat_function( # shrub
    fun = function(x) smaplt_shrub_slope * x + smaplt_shrub_intercept, 
    xlim = c(min(smaplt_shrub_data$LPC, na.rm=T), max(smaplt_shrub_data$LPC, na.rm=T)),
    color = "#7998AD", linewidth = 1.0
  ) +
  labs(subtitle = "LES fast-slow (LMA-LPC)", x = "lnLPC", y = "lnLMA") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#EECA40", "#7998AD", "#F07673")) + 
  annotate("text", x = 0.52, y = 0.016, label = model.sma.LMA_LPC_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.LMA_LPC, filename = "species-level-code/smaplt.LMA_LPC.pdf", width = 3, height = 3)


# PART 8: Try to merge all these plots (failed) --------------------------------
# smaplt_auxpredictor_all = cowplot::plot_grid(
#   smaplt.RDMC_RTD, smaplt.RDMC_RNC, smaplt.LNC_LPC, smaplt.LMA_LPC,
#   nrow = 2, ncol = 2, labels = c("A", "B", "C", "D")
# )
