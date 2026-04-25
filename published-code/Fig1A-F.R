rm(list=ls())
library(tidyverse)
library(dplyr)
library(smatr) # SMA regression
library(ggthemes)
library(patchwork) # plot merging
library(cowplot) # plot merging

load("outdata/traitDataFujian-Ind-step1.RData")
if (!dir.exists("outplts/Fig1")) {
dir.create("outplts/Fig1", showWarnings = FALSE, recursive = TRUE)
}

# (25.2.21) Resolve nomination between "speciesFullName" and "SpeciesFullName"
traitData_touse = traitDataIndv_SelectedTraits_log
traitData_touse = traitData_touse %>% rename(speciesFullName=SpeciesFullName)

# divide data based on growth forms for separate regression
traitData_touse_hilltop = traitData_touse %>% dplyr::filter(SiteID == "hilltop")
traitData_touse_valley = traitData_touse %>% dplyr::filter(SiteID == "valley")


# Fig 1A: LES (LMA-LNC) ----
# 
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
  "M:R=", 
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
  
  "H:R=",
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
  
  "V:R=",
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
    color = "#547bb4", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$LMA, na.rm=T), max(smaplt_valley_data$LMA, na.rm=T)),
    color = "#dd7c4f", linewidth = 1.0
  ) +
  labs(subtitle = "LES (LNC-LMA)", x = "lnLMA", y = "lnLNC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#547bb4", "#dd7c4f")) + 
  annotate("text", x = 0.008, y = 1.18, label = model.sma.LNC_LMA_annotlabel, hjust = 0, vjust = 1, size = 3.2)

suppressWarnings(ggsave(plot = smaplt.LNC_LMA, filename = "outplts/Fig1/Fig1A-IndSites_sma.LNC_LMA.pdf", width = 2.8, height = 2.8))


# Fig 1B: RES collaboration (RD-SRL) ----
#
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
  "M:R=", 
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
  
  "H:R=",
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
  
  "V:R=",
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
    color = "#547bb4", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$SRL, na.rm=T), max(smaplt_valley_data$SRL, na.rm=T)),
    color = "#dd7c4f", linewidth = 1.0
  ) +
  labs(subtitle = "RES collaboration (RD-SRL)", x = "lnSRL", y = "lnRD") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#547bb4", "#dd7c4f")) +
  annotate("text", x = 3.6, y = 0.53, label = model.sma.RD_SRL_annotlabel, hjust = 0, vjust = 1, size = 3.2)

suppressWarnings(ggsave(plot = smaplt.RD_SRL, filename = "outplts/Fig1/Fig1B-IndSites_sma.RD_SRL.pdf", width = 2.8, height = 2.8))


# Fig 1C: RES conservation (RNC-RTD) ----
# 
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
  "M:R=", 
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
  
  "H:R=",
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
  
  "V:R=",
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
    color = "#547bb4", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$RTD, na.rm=T), max(smaplt_valley_data$RTD, na.rm=T)),
    color = "#dd7c4f", linewidth = 1.0
  ) +
  labs(subtitle = "RES conservation (RNC-RTD)", x = "lnRTD", y = "lnRNC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#547bb4", "#dd7c4f")) +
  annotate("text", x = 0.165, y = 1.48, label = model.sma.RNC_RTD_annotlabel, hjust = 0, vjust = 1, size = 3.2)

suppressWarnings(ggsave(plot = smaplt.RNC_RTD, filename = "outplts/Fig1/Fig1C-IndSites_sma.RNC_RTD.pdf", width = 2.8, height = 2.8))


# Fig 1D: LES (SRR25-Rdark25P) ----
#
model.sma.Rdark25P_SRR25_total = sma(traitData_touse$Rdark25P ~ traitData_touse$SRR25)
model.sma.Rdark25P_SRR25_hilltop = sma(traitData_touse_hilltop$Rdark25P ~ traitData_touse_hilltop$SRR25)
model.sma.Rdark25P_SRR25_valley = sma(traitData_touse_valley$Rdark25P ~ traitData_touse_valley$SRR25)

#plot(model.sma.Rdark25P_SRR25)
summary(model.sma.Rdark25P_SRR25_total)
summary(model.sma.Rdark25P_SRR25_hilltop)
summary(model.sma.Rdark25P_SRR25_valley)

smaplt_total_intercept = model.sma.Rdark25P_SRR25_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.Rdark25P_SRR25_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.Rdark25P_SRR25_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.Rdark25P_SRR25_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.Rdark25P_SRR25_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.Rdark25P_SRR25_total[["coef"]][[1]][2,3]

smaplt_hilltop_intercept = model.sma.Rdark25P_SRR25_hilltop[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_hilltop_intercept_lower = model.sma.Rdark25P_SRR25_hilltop[["coef"]][[1]][1,2]
smaplt_hilltop_intercept_upper = model.sma.Rdark25P_SRR25_hilltop[["coef"]][[1]][1,3] 
smaplt_hilltop_slope = model.sma.Rdark25P_SRR25_hilltop[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_hilltop_slope_lower = model.sma.Rdark25P_SRR25_hilltop[["coef"]][[1]][2,2]
smaplt_hilltop_slope_upper = model.sma.Rdark25P_SRR25_hilltop[["coef"]][[1]][2,3]

smaplt_valley_intercept = model.sma.Rdark25P_SRR25_valley[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_valley_intercept_lower = model.sma.Rdark25P_SRR25_valley[["coef"]][[1]][1,2]
smaplt_valley_intercept_upper = model.sma.Rdark25P_SRR25_valley[["coef"]][[1]][1,3] 
smaplt_valley_slope = model.sma.Rdark25P_SRR25_valley[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_valley_slope_lower = model.sma.Rdark25P_SRR25_valley[["coef"]][[1]][2,2]
smaplt_valley_slope_upper = model.sma.Rdark25P_SRR25_valley[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, SiteID, Rdark25P, SRR25) %>%
  dplyr::mutate(SRR25 = as.numeric(SRR25), Rdark25P = as.numeric(Rdark25P))
smaplt_total_data$predicted_Rdark25P = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$SRR25
smaplt_total_data$predicted_Rdark25P_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$SRR25
smaplt_total_data$predicted_Rdark25P_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$SRR25

smaplt_hilltop_data = traitData_touse_hilltop %>% dplyr::select(speciesFullName, GrowthForm, SiteID, Rdark25P, SRR25) %>%
  dplyr::mutate(SRR25 = as.numeric(SRR25), Rdark25P = as.numeric(Rdark25P))
smaplt_hilltop_data$predicted_Rdark25P = smaplt_hilltop_intercept + smaplt_hilltop_slope * smaplt_hilltop_data$SRR25
smaplt_hilltop_data$predicted_Rdark25P_upper_hilltop = smaplt_hilltop_intercept_upper + smaplt_hilltop_slope_upper * smaplt_hilltop_data$SRR25
smaplt_hilltop_data$predicted_Rdark25P_lower_hilltop = smaplt_hilltop_intercept_lower + smaplt_hilltop_slope_lower * smaplt_hilltop_data$SRR25

smaplt_valley_data = traitData_touse_valley %>% dplyr::select(speciesFullName, GrowthForm, SiteID, Rdark25P, SRR25) %>%
  dplyr::mutate(SRR25 = as.numeric(SRR25), Rdark25P = as.numeric(Rdark25P))
smaplt_valley_data$predicted_Rdark25P = smaplt_valley_intercept + smaplt_valley_slope * smaplt_valley_data$SRR25
smaplt_valley_data$predicted_Rdark25P_upper_valley = smaplt_valley_intercept_upper + smaplt_valley_slope_upper * smaplt_valley_data$SRR25
smaplt_valley_data$predicted_Rdark25P_lower_valley = smaplt_valley_intercept_lower + smaplt_valley_slope_lower * smaplt_valley_data$SRR25

# Annotation, R2 and p-value
model.sma.Rdark25P_SRR25_annotlabel = paste0(
  "M:R=", 
  round(ifelse(
    as.numeric(model.sma.Rdark25P_SRR25_total[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.Rdark25P_SRR25_total[["r2"]])),
    -sqrt(as.numeric(model.sma.Rdark25P_SRR25_total[["r2"]]))
  ), digits = 2), 
  " P",
  ifelse(
    signif(as.numeric(model.sma.Rdark25P_SRR25_total[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.Rdark25P_SRR25_total[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.Rdark25P_SRR25_total[["pval"]]), digits = 2), "\n"))),
  
  "H:R=",
  round(ifelse(
    as.numeric(model.sma.Rdark25P_SRR25_hilltop[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.Rdark25P_SRR25_hilltop[["r2"]])),
    -sqrt(as.numeric(model.sma.Rdark25P_SRR25_hilltop[["r2"]]))
  ), digits = 2), 
  " P", 
  ifelse(
    signif(as.numeric(model.sma.Rdark25P_SRR25_hilltop[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.Rdark25P_SRR25_hilltop[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.Rdark25P_SRR25_hilltop[["pval"]]), digits = 2), "\n"))),
  
  "V:R=",
  round(ifelse(
    as.numeric(model.sma.Rdark25P_SRR25_valley[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.Rdark25P_SRR25_valley[["r2"]])),
    -sqrt(as.numeric(model.sma.Rdark25P_SRR25_valley[["r2"]]))
  ), digits = 2),
  " P", 
  ifelse(
    signif(as.numeric(model.sma.Rdark25P_SRR25_valley[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.Rdark25P_SRR25_valley[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.Rdark25P_SRR25_valley[["pval"]]), digits = 2), "\n")))
)

smaplt.Rdark25P_SRR25 = ggplot(smaplt_total_data, aes(x=SRR25, y=Rdark25P)) +
  geom_point(aes(color=SiteID, shape=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$SRR25, na.rm=T), max(smaplt_total_data$SRR25, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # hilltop
    fun = function(x) smaplt_hilltop_slope * x + smaplt_hilltop_intercept, 
    xlim = c(min(smaplt_hilltop_data$SRR25, na.rm=T), max(smaplt_hilltop_data$SRR25, na.rm=T)),
    color = "#547bb4", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$SRR25, na.rm=T), max(smaplt_valley_data$SRR25, na.rm=T)),
    color = "#dd7c4f", linewidth = 1.0
  ) +
  labs(subtitle = "Rd25-Rr25", x = "lnRr25", y = "lnRd25") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#547bb4", "#dd7c4f")) +
  annotate("text", x = 1.2, y = 1.4, label = model.sma.Rdark25P_SRR25_annotlabel, hjust = 0, vjust = 1, size = 3.2)

suppressWarnings(ggsave(plot = smaplt.Rdark25P_SRR25, filename = "outplts/Fig1/Fig1D-IndSites_sma.Rdark25P_SRR25.pdf", width = 2.8, height = 2.8))


# Fig 1E: LES (LNC-RNC) ----
#
model.sma.RNC_LNC_total = sma(traitData_touse$RNC ~ traitData_touse$LNC)
model.sma.RNC_LNC_hilltop = sma(traitData_touse_hilltop$RNC ~ traitData_touse_hilltop$LNC)
model.sma.RNC_LNC_valley = sma(traitData_touse_valley$RNC ~ traitData_touse_valley$LNC)

#plot(model.sma.RNC_LNC)
summary(model.sma.RNC_LNC_total)
summary(model.sma.RNC_LNC_hilltop)
summary(model.sma.RNC_LNC_valley)

smaplt_total_intercept = model.sma.RNC_LNC_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.RNC_LNC_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.RNC_LNC_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.RNC_LNC_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.RNC_LNC_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.RNC_LNC_total[["coef"]][[1]][2,3]

smaplt_hilltop_intercept = model.sma.RNC_LNC_hilltop[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_hilltop_intercept_lower = model.sma.RNC_LNC_hilltop[["coef"]][[1]][1,2]
smaplt_hilltop_intercept_upper = model.sma.RNC_LNC_hilltop[["coef"]][[1]][1,3] 
smaplt_hilltop_slope = model.sma.RNC_LNC_hilltop[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_hilltop_slope_lower = model.sma.RNC_LNC_hilltop[["coef"]][[1]][2,2]
smaplt_hilltop_slope_upper = model.sma.RNC_LNC_hilltop[["coef"]][[1]][2,3]

smaplt_valley_intercept = model.sma.RNC_LNC_valley[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_valley_intercept_lower = model.sma.RNC_LNC_valley[["coef"]][[1]][1,2]
smaplt_valley_intercept_upper = model.sma.RNC_LNC_valley[["coef"]][[1]][1,3] 
smaplt_valley_slope = model.sma.RNC_LNC_valley[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_valley_slope_lower = model.sma.RNC_LNC_valley[["coef"]][[1]][2,2]
smaplt_valley_slope_upper = model.sma.RNC_LNC_valley[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RNC, LNC) %>%
  dplyr::mutate(LNC = as.numeric(LNC), RNC = as.numeric(RNC))
smaplt_total_data$predicted_RNC = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$LNC
smaplt_total_data$predicted_RNC_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$LNC
smaplt_total_data$predicted_RNC_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$LNC

smaplt_hilltop_data = traitData_touse_hilltop %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RNC, LNC) %>%
  dplyr::mutate(LNC = as.numeric(LNC), RNC = as.numeric(RNC))
smaplt_hilltop_data$predicted_RNC = smaplt_hilltop_intercept + smaplt_hilltop_slope * smaplt_hilltop_data$LNC
smaplt_hilltop_data$predicted_RNC_upper_hilltop = smaplt_hilltop_intercept_upper + smaplt_hilltop_slope_upper * smaplt_hilltop_data$LNC
smaplt_hilltop_data$predicted_RNC_lower_hilltop = smaplt_hilltop_intercept_lower + smaplt_hilltop_slope_lower * smaplt_hilltop_data$LNC

smaplt_valley_data = traitData_touse_valley %>% dplyr::select(speciesFullName, GrowthForm, SiteID, RNC, LNC) %>%
  dplyr::mutate(LNC = as.numeric(LNC), RNC = as.numeric(RNC))
smaplt_valley_data$predicted_RNC = smaplt_valley_intercept + smaplt_valley_slope * smaplt_valley_data$LNC
smaplt_valley_data$predicted_RNC_upper_valley = smaplt_valley_intercept_upper + smaplt_valley_slope_upper * smaplt_valley_data$LNC
smaplt_valley_data$predicted_RNC_lower_valley = smaplt_valley_intercept_lower + smaplt_valley_slope_lower * smaplt_valley_data$LNC

# Annotation, R2 and p-value
model.sma.RNC_LNC_annotlabel = paste0(
  "M:R=", 
  round(ifelse(
    as.numeric(model.sma.RNC_LNC_total[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RNC_LNC_total[["r2"]])),
    -sqrt(as.numeric(model.sma.RNC_LNC_total[["r2"]]))
  ), digits = 2), 
  " P",
  ifelse(
    signif(as.numeric(model.sma.RNC_LNC_total[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RNC_LNC_total[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RNC_LNC_total[["pval"]]), digits = 2), "\n"))),
  
  "H:R=",
  round(ifelse(
    as.numeric(model.sma.RNC_LNC_hilltop[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RNC_LNC_hilltop[["r2"]])),
    -sqrt(as.numeric(model.sma.RNC_LNC_hilltop[["r2"]]))
  ), digits = 2), 
  " P", 
  ifelse(
    signif(as.numeric(model.sma.RNC_LNC_hilltop[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RNC_LNC_hilltop[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RNC_LNC_hilltop[["pval"]]), digits = 2), "\n"))),
  
  "V:R=",
  round(ifelse(
    as.numeric(model.sma.RNC_LNC_valley[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.RNC_LNC_valley[["r2"]])),
    -sqrt(as.numeric(model.sma.RNC_LNC_valley[["r2"]]))
  ), digits = 2),
  " P", 
  ifelse(
    signif(as.numeric(model.sma.RNC_LNC_valley[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.RNC_LNC_valley[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.RNC_LNC_valley[["pval"]]), digits = 2), "\n")))
)

smaplt.RNC_LNC = ggplot(smaplt_total_data, aes(x=LNC, y=RNC)) +
  geom_point(aes(color=SiteID, shape=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$LNC, na.rm=T), max(smaplt_total_data$LNC, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # hilltop
    fun = function(x) smaplt_hilltop_slope * x + smaplt_hilltop_intercept, 
    xlim = c(min(smaplt_hilltop_data$LNC, na.rm=T), max(smaplt_hilltop_data$LNC, na.rm=T)),
    color = "#547bb4", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$LNC, na.rm=T), max(smaplt_valley_data$LNC, na.rm=T)),
    color = "#dd7c4f", linewidth = 1.0
  ) +
  labs(subtitle = "RNC-LNC", x = "lnLNC", y = "lnRNC") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#547bb4", "#dd7c4f")) +
  annotate("text", x = 0.8, y = 1.4, label = model.sma.RNC_LNC_annotlabel, hjust = 0, vjust = 1, size = 3.2)

suppressWarnings(ggsave(plot = smaplt.RNC_LNC, filename = "outplts/Fig1/Fig1E-IndSites_sma.RNC_LNC.pdf", width = 2.8, height = 2.8))

# Fig 1F: (LMA-SRL) ----
#
model.sma.SRL_LMA_total = sma(traitData_touse$SRL ~ traitData_touse$LMA)
model.sma.SRL_LMA_hilltop = sma(traitData_touse_hilltop$SRL ~ traitData_touse_hilltop$LMA)
model.sma.SRL_LMA_valley = sma(traitData_touse_valley$SRL ~ traitData_touse_valley$LMA)

#plot(model.sma.SRL_LMA)
summary(model.sma.SRL_LMA_total)
summary(model.sma.SRL_LMA_hilltop)
summary(model.sma.SRL_LMA_valley)

smaplt_total_intercept = model.sma.SRL_LMA_total[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_total_intercept_lower = model.sma.SRL_LMA_total[["coef"]][[1]][1,2]
smaplt_total_intercept_upper = model.sma.SRL_LMA_total[["coef"]][[1]][1,3] 
smaplt_total_slope = model.sma.SRL_LMA_total[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_total_slope_lower = model.sma.SRL_LMA_total[["coef"]][[1]][2,2]
smaplt_total_slope_upper = model.sma.SRL_LMA_total[["coef"]][[1]][2,3]

smaplt_hilltop_intercept = model.sma.SRL_LMA_hilltop[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_hilltop_intercept_lower = model.sma.SRL_LMA_hilltop[["coef"]][[1]][1,2]
smaplt_hilltop_intercept_upper = model.sma.SRL_LMA_hilltop[["coef"]][[1]][1,3] 
smaplt_hilltop_slope = model.sma.SRL_LMA_hilltop[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_hilltop_slope_lower = model.sma.SRL_LMA_hilltop[["coef"]][[1]][2,2]
smaplt_hilltop_slope_upper = model.sma.SRL_LMA_hilltop[["coef"]][[1]][2,3]

smaplt_valley_intercept = model.sma.SRL_LMA_valley[["coef"]][[1]][1,1] # coef(SMA), elevation
smaplt_valley_intercept_lower = model.sma.SRL_LMA_valley[["coef"]][[1]][1,2]
smaplt_valley_intercept_upper = model.sma.SRL_LMA_valley[["coef"]][[1]][1,3] 
smaplt_valley_slope = model.sma.SRL_LMA_valley[["coef"]][[1]][2,1] # coef(SMA), slope
smaplt_valley_slope_lower = model.sma.SRL_LMA_valley[["coef"]][[1]][2,2]
smaplt_valley_slope_upper = model.sma.SRL_LMA_valley[["coef"]][[1]][2,3]


smaplt_total_data = traitData_touse %>% dplyr::select(speciesFullName, GrowthForm, SiteID, SRL, LMA) %>%
  dplyr::mutate(LMA = as.numeric(LMA), SRL = as.numeric(SRL))
smaplt_total_data$predicted_SRL = smaplt_total_intercept + smaplt_total_slope * smaplt_total_data$LMA
smaplt_total_data$predicted_SRL_upper_total = smaplt_total_intercept_upper + smaplt_total_slope_upper * smaplt_total_data$LMA
smaplt_total_data$predicted_SRL_lower_total = smaplt_total_intercept_lower + smaplt_total_slope_lower * smaplt_total_data$LMA

smaplt_hilltop_data = traitData_touse_hilltop %>% dplyr::select(speciesFullName, GrowthForm, SiteID, SRL, LMA) %>%
  dplyr::mutate(LMA = as.numeric(LMA), SRL = as.numeric(SRL))
smaplt_hilltop_data$predicted_SRL = smaplt_hilltop_intercept + smaplt_hilltop_slope * smaplt_hilltop_data$LMA
smaplt_hilltop_data$predicted_SRL_upper_hilltop = smaplt_hilltop_intercept_upper + smaplt_hilltop_slope_upper * smaplt_hilltop_data$LMA
smaplt_hilltop_data$predicted_SRL_lower_hilltop = smaplt_hilltop_intercept_lower + smaplt_hilltop_slope_lower * smaplt_hilltop_data$LMA

smaplt_valley_data = traitData_touse_valley %>% dplyr::select(speciesFullName, GrowthForm, SiteID, SRL, LMA) %>%
  dplyr::mutate(LMA = as.numeric(LMA), SRL = as.numeric(SRL))
smaplt_valley_data$predicted_SRL = smaplt_valley_intercept + smaplt_valley_slope * smaplt_valley_data$LMA
smaplt_valley_data$predicted_SRL_upper_valley = smaplt_valley_intercept_upper + smaplt_valley_slope_upper * smaplt_valley_data$LMA
smaplt_valley_data$predicted_SRL_lower_valley = smaplt_valley_intercept_lower + smaplt_valley_slope_lower * smaplt_valley_data$LMA

# Annotation, R2 and p-value
model.sma.SRL_LMA_annotlabel = paste0(
  "M:R=", 
  round(ifelse(
    as.numeric(model.sma.SRL_LMA_total[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.SRL_LMA_total[["r2"]])),
    -sqrt(as.numeric(model.sma.SRL_LMA_total[["r2"]]))
  ), digits = 2), 
  " P",
  ifelse(
    signif(as.numeric(model.sma.SRL_LMA_total[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.SRL_LMA_total[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.SRL_LMA_total[["pval"]]), digits = 2), "\n"))),
  
  "H:R=",
  round(ifelse(
    as.numeric(model.sma.SRL_LMA_hilltop[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.SRL_LMA_hilltop[["r2"]])),
    -sqrt(as.numeric(model.sma.SRL_LMA_hilltop[["r2"]]))
  ), digits = 2), 
  " P", 
  ifelse(
    signif(as.numeric(model.sma.SRL_LMA_hilltop[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.SRL_LMA_hilltop[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.SRL_LMA_hilltop[["pval"]]), digits = 2), "\n"))),
  
  "V:R=",
  round(ifelse(
    as.numeric(model.sma.SRL_LMA_valley[["groupsummary"]][["Slope"]]) > 0,
    sqrt(as.numeric(model.sma.SRL_LMA_valley[["r2"]])),
    -sqrt(as.numeric(model.sma.SRL_LMA_valley[["r2"]]))
  ), digits = 2),
  " P", 
  ifelse(
    signif(as.numeric(model.sma.SRL_LMA_valley[["pval"]]), digits = 2)<0.01,
    "<0.01\n", 
    ifelse(signif(as.numeric(model.sma.SRL_LMA_valley[["pval"]]), digits = 2)<0.05,
           "<0.05\n", 
           paste0("=",signif(as.numeric(model.sma.SRL_LMA_valley[["pval"]]), digits = 2), "\n")))
)

smaplt.SRL_LMA = ggplot(smaplt_total_data, aes(x=LMA, y=SRL)) +
  geom_point(aes(color=SiteID, shape=GrowthForm), size=2, alpha=1.0) +  # data points
  stat_function( # total
    fun = function(x) smaplt_total_slope * x + smaplt_total_intercept, 
    xlim = c(min(smaplt_total_data$LMA, na.rm=T), max(smaplt_total_data$LMA, na.rm=T)),
    color = "black", linewidth = 1.0
  ) +
  stat_function( # hilltop
    fun = function(x) smaplt_hilltop_slope * x + smaplt_hilltop_intercept, 
    xlim = c(min(smaplt_hilltop_data$LMA, na.rm=T), max(smaplt_hilltop_data$LMA, na.rm=T)),
    color = "#547bb4", linewidth = 1.0
  ) +
  stat_function( # valley
    fun = function(x) smaplt_valley_slope * x + smaplt_valley_intercept, 
    xlim = c(min(smaplt_valley_data$LMA, na.rm=T), max(smaplt_valley_data$LMA, na.rm=T)),
    color = "#dd7c4f", linewidth = 1.0
  ) +
  labs(subtitle = "SRL-LMA", x = "lnLMA", y = "lnSRL") + 
  theme_minimal() + 
  theme(axis.line = element_line(colour = "black")) +
  theme(legend.position = "none") +
  scale_color_manual(values = c("#547bb4", "#dd7c4f")) +
  annotate("text", x = 0.004, y = 5.5, label = model.sma.SRL_LMA_annotlabel, hjust = 0, vjust = 1, size = 3.2)

suppressWarnings(ggsave(plot = smaplt.SRL_LMA, filename = "outplts/Fig1/Fig1F-IndSites_sma.SRL_LMA.pdf", width = 2.8, height = 2.8))
