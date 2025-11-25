# ****************************************************************
# TODO: Modified for use with individual-level analysis workflow
# in response to HW's advice, 2025.2
#
# Keeping the full 87 coupled aboveground and belowground sampled
# individuals for analysis
# ****************************************************************

rm(list=ls())
library(tidyverse)
library(dplyr)
library(smatr) # SMA regression
library(ggthemes)
library(patchwork) # plot merging
library(cowplot) # plot merging

load("2502-indv-level-code/traitDataFujian-Ind-step1.RData")

### SET DATA TO USE ###
traitData_touse = traitDataIndv_SelectedTraits_log

# TODO: (25.2.21) Resolve historical problem: nomination between "speciesFullName" and "SpeciesFullName"
traitData_touse = traitData_touse %>% rename(speciesFullName=SpeciesFullName)

# divide data based on growth forms for separate regression
traitData_touse_hilltop = traitData_touse %>% dplyr::filter(SiteID == "hilltop")
traitData_touse_valley = traitData_touse %>% dplyr::filter(SiteID == "valley")

# ------------------------------------------------------------------------------
# PART 1: LES (LMA-LNC)
# ------------------------------------------------------------------------------
model.sma.LNC_LMA_total = sma(traitData_touse$LNC ~ traitData_touse$LMA)
model.sma.LNC_LMA_hilltop = sma(traitData_touse_hilltop$LNC ~ traitData_touse_hilltop$LMA)
model.sma.LNC_LMA_valley = sma(traitData_touse_valley$LNC ~ traitData_touse_valley$LMA)

#plot(model.sma.LNC_LMA)
summary(model.sma.LNC_LMA_total)
summary(model.sma.LNC_LMA_hilltop)
summary(model.sma.LNC_LMA_valley)

smaplt_total_intercept = model.sma.LNC_LMA_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.LNC_LMA_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.LNC_LMA_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.LNC_LMA_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.LNC_LMA_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.LNC_LMA_total[["coef"]][[1]][2,3]

smaplt_hilltop_intercept = model.sma.LNC_LMA_hilltop[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_hilltop_intercept_lower = model.sma.LNC_LMA_hilltop[["coef"]][[1]][1,2]
smaplt_hilltop_intercept_upper = model.sma.LNC_LMA_hilltop[["coef"]][[1]][1,3] 
smaplt_hilltop_slope = model.sma.LNC_LMA_hilltop[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_hilltop_slope_lower = model.sma.LNC_LMA_hilltop[["coef"]][[1]][2,2]
smaplt_hilltop_slope_upper = model.sma.LNC_LMA_hilltop[["coef"]][[1]][2,3]

smaplt_valley_intercept = model.sma.LNC_LMA_valley[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_valley_intercept_lower = model.sma.LNC_LMA_valley[["coef"]][[1]][1,2]
smaplt_valley_intercept_upper = model.sma.LNC_LMA_valley[["coef"]][[1]][1,3] 
smaplt_valley_slope = model.sma.LNC_LMA_valley[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_valley_slope_lower = model.sma.LNC_LMA_valley[["coef"]][[1]][2,2]
smaplt_valley_slope_upper = model.sma.LNC_LMA_valley[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, SiteID, LNC, LMA) %>%
  dplyr::mutate(LMA = as.numeric(LMA), LNC = as.numeric(LNC))
smaplt_total_data$predicted_LNC = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$LMA
smaplt_total_data$predicted_LNC_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$LMA
smaplt_total_data$predicted_LNC_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$LMA

smaplt_hilltop_data = traitData_touse_hilltop %>% dplyr::select(speciesFullName, GrowthForm, SiteID, LNC, LMA) %>%
  dplyr::mutate(LMA = as.numeric(LMA), LNC = as.numeric(LNC))
smaplt_hilltop_data$predicted_LNC = smaplt_hilltop_intercept + smaplt_hilltop_slope * smaplt_hilltop_data$LMA
smaplt_hilltop_data$predicted_LNC_upper_hilltop = smaplt_hilltop_intercept_upper + smaplt_hilltop_slope_upper * smaplt_hilltop_data$LMA
smaplt_hilltop_data$predicted_LNC_lower_hilltop = smaplt_hilltop_intercept_lower + smaplt_hilltop_slope_lower * smaplt_hilltop_data$LMA

smaplt_valley_data = traitData_touse_valley %>% dplyr::select(speciesFullName, GrowthForm, SiteID,SiteID, LNC, LMA) %>%
  dplyr::mutate(LMA = as.numeric(LMA), LNC = as.numeric(LNC))
smaplt_valley_data$predicted_LNC = smaplt_valley_intercept + smaplt_valley_slope * smaplt_valley_data$LMA
smaplt_valley_data$predicted_LNC_upper_valley = smaplt_valley_intercept_upper + smaplt_valley_slope_upper * smaplt_valley_data$LMA
smaplt_valley_data$predicted_LNC_lower_valley = smaplt_valley_intercept_lower + smaplt_valley_slope_lower * smaplt_valley_data$LMA

# Annotation, R2 and p-value
model.sma.LNC_LMA_annotlabel = paste0(
  "Merged:R=", 
  round(ifelse(
    as.numeric(model.sma.LNC_LMA_total[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.LNC_LMA_total[["r2"]])),
    -sqrt(as.numeric(model.sma.LNC_LMA_total[["r2"]]))
  ), digits = 2), 
  " P",
  ifelse(
    signif(as.numeric(model.sma.LNC_LMA_total[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.LNC_LMA_total[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.LNC_LMA_total[["pval"]]), digits = 2), "\n"))),
  
  "hilltop:R=",
  round(ifelse(
    as.numeric(model.sma.LNC_LMA_hilltop[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.LNC_LMA_hilltop[["r2"]])),
    -sqrt(as.numeric(model.sma.LNC_LMA_hilltop[["r2"]]))
  ), digits = 2), 
  " P", 
  ifelse(
    signif(as.numeric(model.sma.LNC_LMA_hilltop[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.LNC_LMA_hilltop[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.LNC_LMA_hilltop[["pval"]]), digits = 2), "\n"))),
  
  "valley:R=",
  round(ifelse(
    as.numeric(model.sma.LNC_LMA_valley[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.LNC_LMA_valley[["r2"]])),
    -sqrt(as.numeric(model.sma.LNC_LMA_valley[["r2"]]))
  ), digits = 2),
  " P", 
  ifelse(
    signif(as.numeric(model.sma.LNC_LMA_valley[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.LNC_LMA_valley[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.LNC_LMA_valley[["pval"]]), digits = 2), "\n")))
)

smaplt.LNC_LMA = ggplot(smaplt_total_data, aes(x=LMA, y=LNC)) +
  geom_point(aes(color=SiteID, shape=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$LMA, na.rm=T), max(smaplt_total_data$LMA, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # hilltop
    fun = function(x) smaplt_hilltop_slope * x + smaplt_hilltop_intercept, 
    xlim = c(min(smaplt_hilltop_data$LMA, na.rm=T), max(smaplt_hilltop_data$LMA, na.rm=T)),
    color = "#94c6cd", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$LMA, na.rm=T), max(smaplt_valley_data$LMA, na.rm=T)),
    color = "#cb9475", linewidth = 1.0
    ) +
  labs(subtitle = "LES (LNC-LMA)", x = "lnLMA", y = "lnLNC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#94c6cd", "#cb9475")) + 
  annotate("text", x = 0.008, y = 1.18, label = model.sma.LNC_LMA_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.LNC_LMA, filename = "2502-indv-level-code/IndSites_sma.LNC_LMA.pdf", width = 2.8, height = 2.8)

# ------------------------------------------------------------------------------
# PART 2: RES collaboration (RD-SRL)
# ------------------------------------------------------------------------------
model.sma.RD_SRL_total = sma(traitData_touse$RD ~ traitData_touse$SRL)
model.sma.RD_SRL_hilltop = sma(traitData_touse_hilltop$RD ~ traitData_touse_hilltop$SRL)
model.sma.RD_SRL_valley = sma(traitData_touse_valley$RD ~ traitData_touse_valley$SRL)

#plot(model.sma.RD_SRL)
summary(model.sma.RD_SRL_total)
summary(model.sma.RD_SRL_hilltop)
summary(model.sma.RD_SRL_valley)

smaplt_total_intercept = model.sma.RD_SRL_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.RD_SRL_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.RD_SRL_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.RD_SRL_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.RD_SRL_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.RD_SRL_total[["coef"]][[1]][2,3]

smaplt_hilltop_intercept = model.sma.RD_SRL_hilltop[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_hilltop_intercept_lower = model.sma.RD_SRL_hilltop[["coef"]][[1]][1,2]
smaplt_hilltop_intercept_upper = model.sma.RD_SRL_hilltop[["coef"]][[1]][1,3] 
smaplt_hilltop_slope = model.sma.RD_SRL_hilltop[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_hilltop_slope_lower = model.sma.RD_SRL_hilltop[["coef"]][[1]][2,2]
smaplt_hilltop_slope_upper = model.sma.RD_SRL_hilltop[["coef"]][[1]][2,3]

smaplt_valley_intercept = model.sma.RD_SRL_valley[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_valley_intercept_lower = model.sma.RD_SRL_valley[["coef"]][[1]][1,2]
smaplt_valley_intercept_upper = model.sma.RD_SRL_valley[["coef"]][[1]][1,3] 
smaplt_valley_slope = model.sma.RD_SRL_valley[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_valley_slope_lower = model.sma.RD_SRL_valley[["coef"]][[1]][2,2]
smaplt_valley_slope_upper = model.sma.RD_SRL_valley[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RD, SRL) %>%
  dplyr::mutate(SRL = as.numeric(SRL), RD = as.numeric(RD))
smaplt_total_data$predicted_RD = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$SRL
smaplt_total_data$predicted_RD_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$SRL
smaplt_total_data$predicted_RD_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$SRL

smaplt_hilltop_data = traitData_touse_hilltop %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RD, SRL) %>%
  dplyr::mutate(SRL = as.numeric(SRL), RD = as.numeric(RD))
smaplt_hilltop_data$predicted_RD = smaplt_hilltop_intercept + smaplt_hilltop_slope * smaplt_hilltop_data$SRL
smaplt_hilltop_data$predicted_RD_upper_hilltop = smaplt_hilltop_intercept_upper + smaplt_hilltop_slope_upper * smaplt_hilltop_data$SRL
smaplt_hilltop_data$predicted_RD_lower_hilltop = smaplt_hilltop_intercept_lower + smaplt_hilltop_slope_lower * smaplt_hilltop_data$SRL

smaplt_valley_data = traitData_touse_valley %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RD, SRL) %>%
  dplyr::mutate(SRL = as.numeric(SRL), RD = as.numeric(RD))
smaplt_valley_data$predicted_RD = smaplt_valley_intercept + smaplt_valley_slope * smaplt_valley_data$SRL
smaplt_valley_data$predicted_RD_upper_valley = smaplt_valley_intercept_upper + smaplt_valley_slope_upper * smaplt_valley_data$SRL
smaplt_valley_data$predicted_RD_lower_valley = smaplt_valley_intercept_lower + smaplt_valley_slope_lower * smaplt_valley_data$SRL

# Annotation, R2 and p-value
model.sma.RD_SRL_annotlabel = paste0(
  "Merged:R=", 
  round(ifelse(
    as.numeric(model.sma.RD_SRL_total[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RD_SRL_total[["r2"]])),
    -sqrt(as.numeric(model.sma.RD_SRL_total[["r2"]]))
  ), digits = 2), 
  " P",
  ifelse(
    signif(as.numeric(model.sma.RD_SRL_total[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RD_SRL_total[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RD_SRL_total[["pval"]]), digits = 2), "\n"))),
  
  "hilltop:R=",
  round(ifelse(
    as.numeric(model.sma.RD_SRL_hilltop[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RD_SRL_hilltop[["r2"]])),
    -sqrt(as.numeric(model.sma.RD_SRL_hilltop[["r2"]]))
  ), digits = 2), 
  " P", 
  ifelse(
    signif(as.numeric(model.sma.RD_SRL_hilltop[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RD_SRL_hilltop[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RD_SRL_hilltop[["pval"]]), digits = 2), "\n"))),
  
  "valley:R=",
  round(ifelse(
    as.numeric(model.sma.RD_SRL_valley[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RD_SRL_valley[["r2"]])),
    -sqrt(as.numeric(model.sma.RD_SRL_valley[["r2"]]))
  ), digits = 2),
  " P", 
  ifelse(
    signif(as.numeric(model.sma.RD_SRL_valley[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RD_SRL_valley[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RD_SRL_valley[["pval"]]), digits = 2), "\n")))
)

smaplt.RD_SRL = ggplot(smaplt_total_data, aes(x=SRL, y=RD)) +
  geom_point(aes(color=SiteID, shape=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$SRL, na.rm=T), max(smaplt_total_data$SRL, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # hilltop
    fun = function(x) smaplt_hilltop_slope * x + smaplt_hilltop_intercept, 
    xlim = c(min(smaplt_hilltop_data$SRL, na.rm=T), max(smaplt_hilltop_data$SRL, na.rm=T)),
    color = "#94c6cd", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$SRL, na.rm=T), max(smaplt_valley_data$SRL, na.rm=T)),
    color = "#cb9475", linewidth = 1.0
  ) +
  labs(subtitle = "RES collaboration (RD-SRL)", x = "lnSRL", y = "lnRD") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#94c6cd", "#cb9475")) +
  annotate("text", x = 3.6, y = 0.53, label = model.sma.RD_SRL_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.RD_SRL, filename = "2502-indv-level-code/IndSites_sma.RD_SRL.pdf", width = 2.8, height = 2.8)

# ------------------------------------------------------------------------------
# PART 3: RES conservation (RNC-RTD)
# ------------------------------------------------------------------------------
model.sma.RNC_RTD_total = sma(traitData_touse$RNC ~ traitData_touse$RTD)
model.sma.RNC_RTD_hilltop = sma(traitData_touse_hilltop$RNC ~ traitData_touse_hilltop$RTD)
model.sma.RNC_RTD_valley = sma(traitData_touse_valley$RNC ~ traitData_touse_valley$RTD)

#plot(model.sma.RNC_RTD)
summary(model.sma.RNC_RTD_total)
summary(model.sma.RNC_RTD_hilltop)
summary(model.sma.RNC_RTD_valley)

smaplt_total_intercept = model.sma.RNC_RTD_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.RNC_RTD_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.RNC_RTD_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.RNC_RTD_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.RNC_RTD_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.RNC_RTD_total[["coef"]][[1]][2,3]

smaplt_hilltop_intercept = model.sma.RNC_RTD_hilltop[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_hilltop_intercept_lower = model.sma.RNC_RTD_hilltop[["coef"]][[1]][1,2]
smaplt_hilltop_intercept_upper = model.sma.RNC_RTD_hilltop[["coef"]][[1]][1,3] 
smaplt_hilltop_slope = model.sma.RNC_RTD_hilltop[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_hilltop_slope_lower = model.sma.RNC_RTD_hilltop[["coef"]][[1]][2,2]
smaplt_hilltop_slope_upper = model.sma.RNC_RTD_hilltop[["coef"]][[1]][2,3]

smaplt_valley_intercept = model.sma.RNC_RTD_valley[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_valley_intercept_lower = model.sma.RNC_RTD_valley[["coef"]][[1]][1,2]
smaplt_valley_intercept_upper = model.sma.RNC_RTD_valley[["coef"]][[1]][1,3] 
smaplt_valley_slope = model.sma.RNC_RTD_valley[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_valley_slope_lower = model.sma.RNC_RTD_valley[["coef"]][[1]][2,2]
smaplt_valley_slope_upper = model.sma.RNC_RTD_valley[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RNC, RTD) %>%
  dplyr::mutate(RTD = as.numeric(RTD), RNC = as.numeric(RNC))
smaplt_total_data$predicted_RNC = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$RTD
smaplt_total_data$predicted_RNC_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$RTD
smaplt_total_data$predicted_RNC_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$RTD

smaplt_hilltop_data = traitData_touse_hilltop %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RNC, RTD) %>%
  dplyr::mutate(RTD = as.numeric(RTD), RNC = as.numeric(RNC))
smaplt_hilltop_data$predicted_RNC = smaplt_hilltop_intercept + smaplt_hilltop_slope * smaplt_hilltop_data$RTD
smaplt_hilltop_data$predicted_RNC_upper_hilltop = smaplt_hilltop_intercept_upper + smaplt_hilltop_slope_upper * smaplt_hilltop_data$RTD
smaplt_hilltop_data$predicted_RNC_lower_hilltop = smaplt_hilltop_intercept_lower + smaplt_hilltop_slope_lower * smaplt_hilltop_data$RTD

smaplt_valley_data = traitData_touse_valley %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RNC, RTD) %>%
  dplyr::mutate(RTD = as.numeric(RTD), RNC = as.numeric(RNC))
smaplt_valley_data$predicted_RNC = smaplt_valley_intercept + smaplt_valley_slope * smaplt_valley_data$RTD
smaplt_valley_data$predicted_RNC_upper_valley = smaplt_valley_intercept_upper + smaplt_valley_slope_upper * smaplt_valley_data$RTD
smaplt_valley_data$predicted_RNC_lower_valley = smaplt_valley_intercept_lower + smaplt_valley_slope_lower * smaplt_valley_data$RTD

# Annotation, R2 and p-value
model.sma.RNC_RTD_annotlabel = paste0(
  "Merged:R=", 
  round(ifelse(
    as.numeric(model.sma.RNC_RTD_total[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RNC_RTD_total[["r2"]])),
    -sqrt(as.numeric(model.sma.RNC_RTD_total[["r2"]]))
  ), digits = 2), 
  " P",
  ifelse(
    signif(as.numeric(model.sma.RNC_RTD_total[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RNC_RTD_total[["pval"]]), digits = 2)<0.05,
      "<0.05\n", 
      paste0("=",signif(as.numeric(model.sma.RNC_RTD_total[["pval"]]), digits = 2), "\n"))),
  
  "hilltop:R=",
  round(ifelse(
    as.numeric(model.sma.RNC_RTD_hilltop[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RNC_RTD_hilltop[["r2"]])),
    -sqrt(as.numeric(model.sma.RNC_RTD_hilltop[["r2"]]))
  ), digits = 2), 
  " P", 
  ifelse(
    signif(as.numeric(model.sma.RNC_RTD_hilltop[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RNC_RTD_hilltop[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RNC_RTD_hilltop[["pval"]]), digits = 2), "\n"))),
  
  "valley:R=",
  round(ifelse(
    as.numeric(model.sma.RNC_RTD_valley[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RNC_RTD_valley[["r2"]])),
    -sqrt(as.numeric(model.sma.RNC_RTD_valley[["r2"]]))
  ), digits = 2),
  " P", 
  ifelse(
    signif(as.numeric(model.sma.RNC_RTD_valley[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RNC_RTD_valley[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RNC_RTD_valley[["pval"]]), digits = 2), "\n")))
)

smaplt.RNC_RTD = ggplot(smaplt_total_data, aes(x=RTD, y=RNC)) +
  geom_point(aes(color=SiteID, shape=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$RTD, na.rm=T), max(smaplt_total_data$RTD, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # hilltop
    fun = function(x) smaplt_hilltop_slope * x + smaplt_hilltop_intercept, 
    xlim = c(min(smaplt_hilltop_data$RTD, na.rm=T), max(smaplt_hilltop_data$RTD, na.rm=T)),
    color = "#94c6cd", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$RTD, na.rm=T), max(smaplt_valley_data$RTD, na.rm=T)),
    color = "#cb9475", linewidth = 1.0
  ) +
  labs(subtitle = "RES conservation (RNC-RTD)", x = "lnRTD", y = "lnRNC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#94c6cd", "#cb9475")) +
  annotate("text", x = 0.165, y = 1.48, label = model.sma.RNC_RTD_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.RNC_RTD, filename = "2502-indv-level-code/IndSites_sma.RNC_RTD.pdf", width = 2.8, height = 2.8)


# ------------------------------------------------------------------------------
# PART 4: RDMC-RTD
# ------------------------------------------------------------------------------
model.sma.RDMC_RTD_total = sma(traitData_touse$RDMC ~ traitData_touse$RTD)
model.sma.RDMC_RTD_hilltop = sma(traitData_touse_hilltop$RDMC ~ traitData_touse_hilltop$RTD)
model.sma.RDMC_RTD_valley = sma(traitData_touse_valley$RDMC ~ traitData_touse_valley$RTD)

#plot(model.sma.RDMC_RTD)
summary(model.sma.RDMC_RTD_total)
summary(model.sma.RDMC_RTD_hilltop)
summary(model.sma.RDMC_RTD_valley)

smaplt_total_intercept = model.sma.RDMC_RTD_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.RDMC_RTD_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.RDMC_RTD_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.RDMC_RTD_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.RDMC_RTD_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.RDMC_RTD_total[["coef"]][[1]][2,3]

smaplt_hilltop_intercept = model.sma.RDMC_RTD_hilltop[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_hilltop_intercept_lower = model.sma.RDMC_RTD_hilltop[["coef"]][[1]][1,2]
smaplt_hilltop_intercept_upper = model.sma.RDMC_RTD_hilltop[["coef"]][[1]][1,3] 
smaplt_hilltop_slope = model.sma.RDMC_RTD_hilltop[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_hilltop_slope_lower = model.sma.RDMC_RTD_hilltop[["coef"]][[1]][2,2]
smaplt_hilltop_slope_upper = model.sma.RDMC_RTD_hilltop[["coef"]][[1]][2,3]

smaplt_valley_intercept = model.sma.RDMC_RTD_valley[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_valley_intercept_lower = model.sma.RDMC_RTD_valley[["coef"]][[1]][1,2]
smaplt_valley_intercept_upper = model.sma.RDMC_RTD_valley[["coef"]][[1]][1,3] 
smaplt_valley_slope = model.sma.RDMC_RTD_valley[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_valley_slope_lower = model.sma.RDMC_RTD_valley[["coef"]][[1]][2,2]
smaplt_valley_slope_upper = model.sma.RDMC_RTD_valley[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RDMC, RTD) %>%
  dplyr::mutate(RTD = as.numeric(RTD), RDMC = as.numeric(RDMC))
smaplt_total_data$predicted_RDMC = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$RTD
smaplt_total_data$predicted_RDMC_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$RTD
smaplt_total_data$predicted_RDMC_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$RTD

smaplt_hilltop_data = traitData_touse_hilltop %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RDMC, RTD) %>%
  dplyr::mutate(RTD = as.numeric(RTD), RDMC = as.numeric(RDMC))
smaplt_hilltop_data$predicted_RDMC = smaplt_hilltop_intercept + smaplt_hilltop_slope * smaplt_hilltop_data$RTD
smaplt_hilltop_data$predicted_RDMC_upper_hilltop = smaplt_hilltop_intercept_upper + smaplt_hilltop_slope_upper * smaplt_hilltop_data$RTD
smaplt_hilltop_data$predicted_RDMC_lower_hilltop = smaplt_hilltop_intercept_lower + smaplt_hilltop_slope_lower * smaplt_hilltop_data$RTD

smaplt_valley_data = traitData_touse_valley %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RDMC, RTD) %>%
  dplyr::mutate(RTD = as.numeric(RTD), RDMC = as.numeric(RDMC))
smaplt_valley_data$predicted_RDMC = smaplt_valley_intercept + smaplt_valley_slope * smaplt_valley_data$RTD
smaplt_valley_data$predicted_RDMC_upper_valley = smaplt_valley_intercept_upper + smaplt_valley_slope_upper * smaplt_valley_data$RTD
smaplt_valley_data$predicted_RDMC_lower_valley = smaplt_valley_intercept_lower + smaplt_valley_slope_lower * smaplt_valley_data$RTD

# Annotation, R2 and p-value
model.sma.RDMC_RTD_annotlabel = paste0(
  "Merged:R=", 
  round(ifelse(
    as.numeric(model.sma.RDMC_RTD_total[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RDMC_RTD_total[["r2"]])),
    -sqrt(as.numeric(model.sma.RDMC_RTD_total[["r2"]]))
  ), digits = 2), 
  " P",
  ifelse(
    signif(as.numeric(model.sma.RDMC_RTD_total[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RDMC_RTD_total[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RDMC_RTD_total[["pval"]]), digits = 2), "\n"))),
  
  "hilltop:R=",
  round(ifelse(
    as.numeric(model.sma.RDMC_RTD_hilltop[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RDMC_RTD_hilltop[["r2"]])),
    -sqrt(as.numeric(model.sma.RDMC_RTD_hilltop[["r2"]]))
  ), digits = 2), 
  " P", 
  ifelse(
    signif(as.numeric(model.sma.RDMC_RTD_hilltop[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RDMC_RTD_hilltop[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RDMC_RTD_hilltop[["pval"]]), digits = 2), "\n"))),
  
  "valley:R=",
  round(ifelse(
    as.numeric(model.sma.RDMC_RTD_valley[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RDMC_RTD_valley[["r2"]])),
    -sqrt(as.numeric(model.sma.RDMC_RTD_valley[["r2"]]))
  ), digits = 2),
  " P", 
  ifelse(
    signif(as.numeric(model.sma.RDMC_RTD_valley[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RDMC_RTD_valley[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RDMC_RTD_valley[["pval"]]), digits = 2), "\n")))
)

smaplt.RDMC_RTD = ggplot(smaplt_total_data, aes(x=RTD, y=RDMC)) +
  geom_point(aes(color=SiteID, shape=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$RTD, na.rm=T), max(smaplt_total_data$RTD, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # hilltop
    fun = function(x) smaplt_hilltop_slope * x + smaplt_hilltop_intercept, 
    xlim = c(min(smaplt_hilltop_data$RTD, na.rm=T), max(smaplt_hilltop_data$RTD, na.rm=T)),
    color = "#94c6cd", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$RTD, na.rm=T), max(smaplt_valley_data$RTD, na.rm=T)),
    color = "#cb9475", linewidth = 1.0
  ) +
  labs(subtitle = "RES conservation (RDMC-RTD)", x = "lnRTD", y = "lnRDMC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#94c6cd", "#cb9475")) +
  annotate("text", x = 0.05, y = 0.42, label = model.sma.RDMC_RTD_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.RDMC_RTD, filename = "2502-indv-level-code/IndSites_sma.RDMC_RTD.pdf", width = 2.8, height = 2.8)


# ------------------------------------------------------------------------------
# PART 5: RDMC-RNC
# ------------------------------------------------------------------------------
model.sma.RDMC_RNC_total = sma(traitData_touse$RDMC ~ traitData_touse$RNC)
model.sma.RDMC_RNC_hilltop = sma(traitData_touse_hilltop$RDMC ~ traitData_touse_hilltop$RNC)
model.sma.RDMC_RNC_valley = sma(traitData_touse_valley$RDMC ~ traitData_touse_valley$RNC)

#plot(model.sma.RDMC_RNC)
summary(model.sma.RDMC_RNC_total)
summary(model.sma.RDMC_RNC_hilltop)
summary(model.sma.RDMC_RNC_valley)

smaplt_total_intercept = model.sma.RDMC_RNC_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.RDMC_RNC_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.RDMC_RNC_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.RDMC_RNC_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.RDMC_RNC_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.RDMC_RNC_total[["coef"]][[1]][2,3]

smaplt_hilltop_intercept = model.sma.RDMC_RNC_hilltop[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_hilltop_intercept_lower = model.sma.RDMC_RNC_hilltop[["coef"]][[1]][1,2]
smaplt_hilltop_intercept_upper = model.sma.RDMC_RNC_hilltop[["coef"]][[1]][1,3] 
smaplt_hilltop_slope = model.sma.RDMC_RNC_hilltop[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_hilltop_slope_lower = model.sma.RDMC_RNC_hilltop[["coef"]][[1]][2,2]
smaplt_hilltop_slope_upper = model.sma.RDMC_RNC_hilltop[["coef"]][[1]][2,3]

smaplt_valley_intercept = model.sma.RDMC_RNC_valley[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_valley_intercept_lower = model.sma.RDMC_RNC_valley[["coef"]][[1]][1,2]
smaplt_valley_intercept_upper = model.sma.RDMC_RNC_valley[["coef"]][[1]][1,3] 
smaplt_valley_slope = model.sma.RDMC_RNC_valley[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_valley_slope_lower = model.sma.RDMC_RNC_valley[["coef"]][[1]][2,2]
smaplt_valley_slope_upper = model.sma.RDMC_RNC_valley[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RDMC, RNC) %>%
  dplyr::mutate(RNC = as.numeric(RNC), RDMC = as.numeric(RDMC))
smaplt_total_data$predicted_RDMC = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$RNC
smaplt_total_data$predicted_RDMC_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$RNC
smaplt_total_data$predicted_RDMC_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$RNC

smaplt_hilltop_data = traitData_touse_hilltop %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RDMC, RNC) %>%
  dplyr::mutate(RNC = as.numeric(RNC), RDMC = as.numeric(RDMC))
smaplt_hilltop_data$predicted_RDMC = smaplt_hilltop_intercept + smaplt_hilltop_slope * smaplt_hilltop_data$RNC
smaplt_hilltop_data$predicted_RDMC_upper_hilltop = smaplt_hilltop_intercept_upper + smaplt_hilltop_slope_upper * smaplt_hilltop_data$RNC
smaplt_hilltop_data$predicted_RDMC_lower_hilltop = smaplt_hilltop_intercept_lower + smaplt_hilltop_slope_lower * smaplt_hilltop_data$RNC

smaplt_valley_data = traitData_touse_valley %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RDMC, RNC) %>%
  dplyr::mutate(RNC = as.numeric(RNC), RDMC = as.numeric(RDMC))
smaplt_valley_data$predicted_RDMC = smaplt_valley_intercept + smaplt_valley_slope * smaplt_valley_data$RNC
smaplt_valley_data$predicted_RDMC_upper_valley = smaplt_valley_intercept_upper + smaplt_valley_slope_upper * smaplt_valley_data$RNC
smaplt_valley_data$predicted_RDMC_lower_valley = smaplt_valley_intercept_lower + smaplt_valley_slope_lower * smaplt_valley_data$RNC

# Annotation, R2 and p-value
model.sma.RDMC_RNC_annotlabel = paste0(
  "Merged:R=", 
  round(ifelse(
    as.numeric(model.sma.RDMC_RNC_total[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RDMC_RNC_total[["r2"]])),
    -sqrt(as.numeric(model.sma.RDMC_RNC_total[["r2"]]))
  ), digits = 2), 
  " P",
  ifelse(
    signif(as.numeric(model.sma.RDMC_RNC_total[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RDMC_RNC_total[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RDMC_RNC_total[["pval"]]), digits = 2), "\n"))),
  
  "hilltop:R=",
  round(ifelse(
    as.numeric(model.sma.RDMC_RNC_hilltop[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RDMC_RNC_hilltop[["r2"]])),
    -sqrt(as.numeric(model.sma.RDMC_RNC_hilltop[["r2"]]))
  ), digits = 2), 
  " P", 
  ifelse(
    signif(as.numeric(model.sma.RDMC_RNC_hilltop[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RDMC_RNC_hilltop[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RDMC_RNC_hilltop[["pval"]]), digits = 2), "\n"))),
  
  "valley:R=",
  round(ifelse(
    as.numeric(model.sma.RDMC_RNC_valley[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RDMC_RNC_valley[["r2"]])),
    -sqrt(as.numeric(model.sma.RDMC_RNC_valley[["r2"]]))
  ), digits = 2),
  " P", 
  ifelse(
    signif(as.numeric(model.sma.RDMC_RNC_valley[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RDMC_RNC_valley[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RDMC_RNC_valley[["pval"]]), digits = 2), "\n")))
)

smaplt.RDMC_RNC = ggplot(smaplt_total_data, aes(x=RNC, y=RDMC)) +
  geom_point(aes(color=SiteID, shape=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$RNC, na.rm=T), max(smaplt_total_data$RNC, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # hilltop
    fun = function(x) smaplt_hilltop_slope * x + smaplt_hilltop_intercept, 
    xlim = c(min(smaplt_hilltop_data$RNC, na.rm=T), max(smaplt_hilltop_data$RNC, na.rm=T)),
    color = "#94c6cd", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$RNC, na.rm=T), max(smaplt_valley_data$RNC, na.rm=T)),
    color = "#cb9475", linewidth = 1.0
  ) +
  labs(subtitle = "RES conservation (RDMC-RNC)", x = "lnRNC", y = "lnRDMC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#94c6cd", "#cb9475")) +
  annotate("text", x = 0.87, y = 0.335, label = model.sma.RDMC_RNC_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.RDMC_RNC, filename = "2502-indv-level-code/IndSites_sma.RDMC_RNC.pdf", width = 2.8, height = 2.8)


# ------------------------------------------------------------------------------
# PART 6: LES (LPC-LNC)
# ------------------------------------------------------------------------------
model.sma.LNC_LPC_total = sma(traitData_touse$LNC ~ traitData_touse$LPC)
model.sma.LNC_LPC_hilltop = sma(traitData_touse_hilltop$LNC ~ traitData_touse_hilltop$LPC)
model.sma.LNC_LPC_valley = sma(traitData_touse_valley$LNC ~ traitData_touse_valley$LPC)

#plot(model.sma.LNC_LPC)
summary(model.sma.LNC_LPC_total)
summary(model.sma.LNC_LPC_hilltop)
summary(model.sma.LNC_LPC_valley)

smaplt_total_intercept = model.sma.LNC_LPC_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.LNC_LPC_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.LNC_LPC_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.LNC_LPC_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.LNC_LPC_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.LNC_LPC_total[["coef"]][[1]][2,3]

smaplt_hilltop_intercept = model.sma.LNC_LPC_hilltop[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_hilltop_intercept_lower = model.sma.LNC_LPC_hilltop[["coef"]][[1]][1,2]
smaplt_hilltop_intercept_upper = model.sma.LNC_LPC_hilltop[["coef"]][[1]][1,3] 
smaplt_hilltop_slope = model.sma.LNC_LPC_hilltop[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_hilltop_slope_lower = model.sma.LNC_LPC_hilltop[["coef"]][[1]][2,2]
smaplt_hilltop_slope_upper = model.sma.LNC_LPC_hilltop[["coef"]][[1]][2,3]

smaplt_valley_intercept = model.sma.LNC_LPC_valley[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_valley_intercept_lower = model.sma.LNC_LPC_valley[["coef"]][[1]][1,2]
smaplt_valley_intercept_upper = model.sma.LNC_LPC_valley[["coef"]][[1]][1,3] 
smaplt_valley_slope = model.sma.LNC_LPC_valley[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_valley_slope_lower = model.sma.LNC_LPC_valley[["coef"]][[1]][2,2]
smaplt_valley_slope_upper = model.sma.LNC_LPC_valley[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, SiteID, LNC, LPC) %>%
  dplyr::mutate(LPC = as.numeric(LPC), LNC = as.numeric(LNC))
smaplt_total_data$predicted_LNC = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$LPC
smaplt_total_data$predicted_LNC_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$LPC
smaplt_total_data$predicted_LNC_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$LPC

smaplt_hilltop_data = traitData_touse_hilltop %>% dplyr::select(speciesFullName, GrowthForm, SiteID, LNC, LPC) %>%
  dplyr::mutate(LPC = as.numeric(LPC), LNC = as.numeric(LNC))
smaplt_hilltop_data$predicted_LNC = smaplt_hilltop_intercept + smaplt_hilltop_slope * smaplt_hilltop_data$LPC
smaplt_hilltop_data$predicted_LNC_upper_hilltop = smaplt_hilltop_intercept_upper + smaplt_hilltop_slope_upper * smaplt_hilltop_data$LPC
smaplt_hilltop_data$predicted_LNC_lower_hilltop = smaplt_hilltop_intercept_lower + smaplt_hilltop_slope_lower * smaplt_hilltop_data$LPC

smaplt_valley_data = traitData_touse_valley %>% dplyr::select(speciesFullName, GrowthForm, SiteID, LNC, LPC) %>%
  dplyr::mutate(LPC = as.numeric(LPC), LNC = as.numeric(LNC))
smaplt_valley_data$predicted_LNC = smaplt_valley_intercept + smaplt_valley_slope * smaplt_valley_data$LPC
smaplt_valley_data$predicted_LNC_upper_valley = smaplt_valley_intercept_upper + smaplt_valley_slope_upper * smaplt_valley_data$LPC
smaplt_valley_data$predicted_LNC_lower_valley = smaplt_valley_intercept_lower + smaplt_valley_slope_lower * smaplt_valley_data$LPC

# Annotation, R2 and p-value
model.sma.LNC_LPC_annotlabel = paste0(
  "Merged:R=", 
  round(ifelse(
    as.numeric(model.sma.LNC_LPC_total[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.LNC_LPC_total[["r2"]])),
    -sqrt(as.numeric(model.sma.LNC_LPC_total[["r2"]]))
  ), digits = 2), 
  " P",
  ifelse(
    signif(as.numeric(model.sma.LNC_LPC_total[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.LNC_LPC_total[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.LNC_LPC_total[["pval"]]), digits = 2), "\n"))),
  
  "hilltop:R=",
  round(ifelse(
    as.numeric(model.sma.LNC_LPC_hilltop[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.LNC_LPC_hilltop[["r2"]])),
    -sqrt(as.numeric(model.sma.LNC_LPC_hilltop[["r2"]]))
  ), digits = 2), 
  " P", 
  ifelse(
    signif(as.numeric(model.sma.LNC_LPC_hilltop[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.LNC_LPC_hilltop[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.LNC_LPC_hilltop[["pval"]]), digits = 2), "\n"))),
  
  "valley:R=",
  round(ifelse(
    as.numeric(model.sma.LNC_LPC_valley[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.LNC_LPC_valley[["r2"]])),
    -sqrt(as.numeric(model.sma.LNC_LPC_valley[["r2"]]))
  ), digits = 2),
  " P", 
  ifelse(
    signif(as.numeric(model.sma.LNC_LPC_valley[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.LNC_LPC_valley[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.LNC_LPC_valley[["pval"]]), digits = 2), "\n")))
)

smaplt.LNC_LPC = ggplot(smaplt_total_data, aes(x=LPC, y=LNC)) +
  geom_point(aes(color=SiteID, shape=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$LPC, na.rm=T), max(smaplt_total_data$LPC, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # hilltop
    fun = function(x) smaplt_hilltop_slope * x + smaplt_hilltop_intercept, 
    xlim = c(min(smaplt_hilltop_data$LPC, na.rm=T), max(smaplt_hilltop_data$LPC, na.rm=T)),
    color = "#94c6cd", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$LPC, na.rm=T), max(smaplt_valley_data$LPC, na.rm=T)),
    color = "#cb9475", linewidth = 1.0
  ) +
  labs(subtitle = "LES (LNC-LPC)", x = "lnLPC", y = "lnLNC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#94c6cd", "#cb9475")) +
  annotate("text", x = 0.37, y = 1.35, label = model.sma.LNC_LPC_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.LNC_LPC, filename = "2502-indv-level-code/IndSites_sma.LNC_LPC.pdf", width = 2.8, height = 2.8)


# ------------------------------------------------------------------------------
# PART 7: LES (LPC-LMA)
# ------------------------------------------------------------------------------
model.sma.LMA_LPC_total = sma(traitData_touse$LMA ~ traitData_touse$LPC)
model.sma.LMA_LPC_hilltop = sma(traitData_touse_hilltop$LMA ~ traitData_touse_hilltop$LPC)
model.sma.LMA_LPC_valley = sma(traitData_touse_valley$LMA ~ traitData_touse_valley$LPC)

#plot(model.sma.LMA_LPC)
summary(model.sma.LMA_LPC_total)
summary(model.sma.LMA_LPC_hilltop)
summary(model.sma.LMA_LPC_valley)

smaplt_total_intercept = model.sma.LMA_LPC_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.LMA_LPC_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.LMA_LPC_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.LMA_LPC_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.LMA_LPC_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.LMA_LPC_total[["coef"]][[1]][2,3]

smaplt_hilltop_intercept = model.sma.LMA_LPC_hilltop[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_hilltop_intercept_lower = model.sma.LMA_LPC_hilltop[["coef"]][[1]][1,2]
smaplt_hilltop_intercept_upper = model.sma.LMA_LPC_hilltop[["coef"]][[1]][1,3] 
smaplt_hilltop_slope = model.sma.LMA_LPC_hilltop[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_hilltop_slope_lower = model.sma.LMA_LPC_hilltop[["coef"]][[1]][2,2]
smaplt_hilltop_slope_upper = model.sma.LMA_LPC_hilltop[["coef"]][[1]][2,3]

smaplt_valley_intercept = model.sma.LMA_LPC_valley[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_valley_intercept_lower = model.sma.LMA_LPC_valley[["coef"]][[1]][1,2]
smaplt_valley_intercept_upper = model.sma.LMA_LPC_valley[["coef"]][[1]][1,3] 
smaplt_valley_slope = model.sma.LMA_LPC_valley[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_valley_slope_lower = model.sma.LMA_LPC_valley[["coef"]][[1]][2,2]
smaplt_valley_slope_upper = model.sma.LMA_LPC_valley[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, SiteID, LMA, LPC) %>%
  dplyr::mutate(LPC = as.numeric(LPC), LMA = as.numeric(LMA))
smaplt_total_data$predicted_LMA = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$LPC
smaplt_total_data$predicted_LMA_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$LPC
smaplt_total_data$predicted_LMA_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$LPC

smaplt_hilltop_data = traitData_touse_hilltop %>% dplyr::select(speciesFullName, GrowthForm, SiteID, LMA, LPC) %>%
  dplyr::mutate(LPC = as.numeric(LPC), LMA = as.numeric(LMA))
smaplt_hilltop_data$predicted_LMA = smaplt_hilltop_intercept + smaplt_hilltop_slope * smaplt_hilltop_data$LPC
smaplt_hilltop_data$predicted_LMA_upper_hilltop = smaplt_hilltop_intercept_upper + smaplt_hilltop_slope_upper * smaplt_hilltop_data$LPC
smaplt_hilltop_data$predicted_LMA_lower_hilltop = smaplt_hilltop_intercept_lower + smaplt_hilltop_slope_lower * smaplt_hilltop_data$LPC

smaplt_valley_data = traitData_touse_valley %>% dplyr::select(speciesFullName, GrowthForm, SiteID, LMA, LPC) %>%
  dplyr::mutate(LPC = as.numeric(LPC), LMA = as.numeric(LMA))
smaplt_valley_data$predicted_LMA = smaplt_valley_intercept + smaplt_valley_slope * smaplt_valley_data$LPC
smaplt_valley_data$predicted_LMA_upper_valley = smaplt_valley_intercept_upper + smaplt_valley_slope_upper * smaplt_valley_data$LPC
smaplt_valley_data$predicted_LMA_lower_valley = smaplt_valley_intercept_lower + smaplt_valley_slope_lower * smaplt_valley_data$LPC

# Annotation, R2 and p-value
model.sma.LMA_LPC_annotlabel = paste0(
  "Merged:R=", 
  round(ifelse(
    as.numeric(model.sma.LMA_LPC_total[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.LMA_LPC_total[["r2"]])),
    -sqrt(as.numeric(model.sma.LMA_LPC_total[["r2"]]))
  ), digits = 2), 
  " P",
  ifelse(
    signif(as.numeric(model.sma.LMA_LPC_total[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.LMA_LPC_total[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.LMA_LPC_total[["pval"]]), digits = 2), "\n"))),
  
  "hilltop:R=",
  round(ifelse(
    as.numeric(model.sma.LMA_LPC_hilltop[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.LMA_LPC_hilltop[["r2"]])),
    -sqrt(as.numeric(model.sma.LMA_LPC_hilltop[["r2"]]))
  ), digits = 2), 
  " P", 
  ifelse(
    signif(as.numeric(model.sma.LMA_LPC_hilltop[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.LMA_LPC_hilltop[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.LMA_LPC_hilltop[["pval"]]), digits = 2), "\n"))),
  
  "valley:R=",
  round(ifelse(
    as.numeric(model.sma.LMA_LPC_valley[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.LMA_LPC_valley[["r2"]])),
    -sqrt(as.numeric(model.sma.LMA_LPC_valley[["r2"]]))
  ), digits = 2),
  " P", 
  ifelse(
    signif(as.numeric(model.sma.LMA_LPC_valley[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.LMA_LPC_valley[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.LMA_LPC_valley[["pval"]]), digits = 2), "\n")))
)

smaplt.LMA_LPC = ggplot(smaplt_total_data, aes(x=LPC, y=LMA)) +
  geom_point(aes(color=SiteID, shape=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$LPC, na.rm=T), max(smaplt_total_data$LPC, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # hilltop
    fun = function(x) smaplt_hilltop_slope * x + smaplt_hilltop_intercept, 
    xlim = c(min(smaplt_hilltop_data$LPC, na.rm=T), max(smaplt_hilltop_data$LPC, na.rm=T)),
    color = "#94c6cd", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$LPC, na.rm=T), max(smaplt_valley_data$LPC, na.rm=T)),
    color = "#cb9475", linewidth = 1.0
  ) +
  labs(subtitle = "LES (LMA-LPC)", x = "lnLPC", y = "lnLMA") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#94c6cd", "#cb9475")) +
  annotate("text", x = 0.52, y = 0.016, label = model.sma.LMA_LPC_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.LMA_LPC, filename = "2502-indv-level-code/IndSites_sma.LMA_LPC.pdf", width = 2.8, height = 2.8)




# NEWLY ADDED: LES (LA-Ld13C) ----
#
model.sma.Ld13C_LA_total = sma(traitData_touse$Ld13C ~ traitData_touse$LA)
model.sma.Ld13C_LA_hilltop = sma(traitData_touse_hilltop$Ld13C ~ traitData_touse_hilltop$LA)
model.sma.Ld13C_LA_valley = sma(traitData_touse_valley$Ld13C ~ traitData_touse_valley$LA)

#plot(model.sma.Ld13C_LA)
summary(model.sma.Ld13C_LA_total)
summary(model.sma.Ld13C_LA_hilltop)
summary(model.sma.Ld13C_LA_valley)

smaplt_total_intercept = model.sma.Ld13C_LA_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.Ld13C_LA_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.Ld13C_LA_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.Ld13C_LA_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.Ld13C_LA_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.Ld13C_LA_total[["coef"]][[1]][2,3]

smaplt_hilltop_intercept = model.sma.Ld13C_LA_hilltop[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_hilltop_intercept_lower = model.sma.Ld13C_LA_hilltop[["coef"]][[1]][1,2]
smaplt_hilltop_intercept_upper = model.sma.Ld13C_LA_hilltop[["coef"]][[1]][1,3] 
smaplt_hilltop_slope = model.sma.Ld13C_LA_hilltop[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_hilltop_slope_lower = model.sma.Ld13C_LA_hilltop[["coef"]][[1]][2,2]
smaplt_hilltop_slope_upper = model.sma.Ld13C_LA_hilltop[["coef"]][[1]][2,3]

smaplt_valley_intercept = model.sma.Ld13C_LA_valley[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_valley_intercept_lower = model.sma.Ld13C_LA_valley[["coef"]][[1]][1,2]
smaplt_valley_intercept_upper = model.sma.Ld13C_LA_valley[["coef"]][[1]][1,3] 
smaplt_valley_slope = model.sma.Ld13C_LA_valley[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_valley_slope_lower = model.sma.Ld13C_LA_valley[["coef"]][[1]][2,2]
smaplt_valley_slope_upper = model.sma.Ld13C_LA_valley[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, SiteID, Ld13C, LA) %>%
  dplyr::mutate(LA = as.numeric(LA), Ld13C = as.numeric(Ld13C))
smaplt_total_data$predicted_Ld13C = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$LA
smaplt_total_data$predicted_Ld13C_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$LA
smaplt_total_data$predicted_Ld13C_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$LA

smaplt_hilltop_data = traitData_touse_hilltop %>% dplyr::select(speciesFullName, GrowthForm, SiteID, Ld13C, LA) %>%
  dplyr::mutate(LA = as.numeric(LA), Ld13C = as.numeric(Ld13C))
smaplt_hilltop_data$predicted_Ld13C = smaplt_hilltop_intercept + smaplt_hilltop_slope * smaplt_hilltop_data$LA
smaplt_hilltop_data$predicted_Ld13C_upper_hilltop = smaplt_hilltop_intercept_upper + smaplt_hilltop_slope_upper * smaplt_hilltop_data$LA
smaplt_hilltop_data$predicted_Ld13C_lower_hilltop = smaplt_hilltop_intercept_lower + smaplt_hilltop_slope_lower * smaplt_hilltop_data$LA

smaplt_valley_data = traitData_touse_valley %>% dplyr::select(speciesFullName, GrowthForm, SiteID, Ld13C, LA) %>%
  dplyr::mutate(LA = as.numeric(LA), Ld13C = as.numeric(Ld13C))
smaplt_valley_data$predicted_Ld13C = smaplt_valley_intercept + smaplt_valley_slope * smaplt_valley_data$LA
smaplt_valley_data$predicted_Ld13C_upper_valley = smaplt_valley_intercept_upper + smaplt_valley_slope_upper * smaplt_valley_data$LA
smaplt_valley_data$predicted_Ld13C_lower_valley = smaplt_valley_intercept_lower + smaplt_valley_slope_lower * smaplt_valley_data$LA

# Annotation, R2 and p-value
model.sma.Ld13C_LA_annotlabel = paste0(
  "Merged:R=", 
  round(ifelse(
    as.numeric(model.sma.Ld13C_LA_total[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.Ld13C_LA_total[["r2"]])),
    -sqrt(as.numeric(model.sma.Ld13C_LA_total[["r2"]]))
  ), digits = 2), 
  " P",
  ifelse(
    signif(as.numeric(model.sma.Ld13C_LA_total[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.Ld13C_LA_total[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.Ld13C_LA_total[["pval"]]), digits = 2), "\n"))),
  
  "hilltop:R=",
  round(ifelse(
    as.numeric(model.sma.Ld13C_LA_hilltop[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.Ld13C_LA_hilltop[["r2"]])),
    -sqrt(as.numeric(model.sma.Ld13C_LA_hilltop[["r2"]]))
  ), digits = 2), 
  " P", 
  ifelse(
    signif(as.numeric(model.sma.Ld13C_LA_hilltop[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.Ld13C_LA_hilltop[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.Ld13C_LA_hilltop[["pval"]]), digits = 2), "\n"))),
  
  "valley:R=",
  round(ifelse(
    as.numeric(model.sma.Ld13C_LA_valley[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.Ld13C_LA_valley[["r2"]])),
    -sqrt(as.numeric(model.sma.Ld13C_LA_valley[["r2"]]))
  ), digits = 2),
  " P", 
  ifelse(
    signif(as.numeric(model.sma.Ld13C_LA_valley[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.Ld13C_LA_valley[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.Ld13C_LA_valley[["pval"]]), digits = 2), "\n")))
)

smaplt.Ld13C_LA = ggplot(smaplt_total_data, aes(x=LA, y=Ld13C)) +
  geom_point(aes(color=SiteID, shape=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$LA, na.rm=T), max(smaplt_total_data$LA, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # hilltop
    fun = function(x) smaplt_hilltop_slope * x + smaplt_hilltop_intercept, 
    xlim = c(min(smaplt_hilltop_data$LA, na.rm=T), max(smaplt_hilltop_data$LA, na.rm=T)),
    color = "#94c6cd", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$LA, na.rm=T), max(smaplt_valley_data$LA, na.rm=T)),
    color = "#cb9475", linewidth = 1.0
  ) +
  labs(subtitle = "LES (Ld13C-LA)", x = "lnLA", y = "lnLd13C") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#94c6cd", "#cb9475")) +
  annotate("text", x = 0.52, y = 0.016, label = model.sma.Ld13C_LA_annotlabel, hjust = 0, vjust = 1, size = 3.2)

ggsave(plot = smaplt.Ld13C_LA, filename = "2502-indv-level-code/IndSites_sma.Ld13C_LA.pdf", width = 2.8, height = 2.8)
