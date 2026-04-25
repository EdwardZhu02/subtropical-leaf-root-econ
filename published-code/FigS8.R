rm(list=ls())
library(tidyverse)
library(dplyr)

load("outdata/traitDataFujian-Ind-step1.RData")
if (!dir.exists("outplts/SIFigures")) {
  dir.create("outplts/SIFigures", showWarnings = FALSE, recursive = TRUE)
}

intrasp_traitdata = traitDataIndv_SelectedTraits %>%
  dplyr::filter(SpeciesFullName %in% c(
    "Castanopsis faberi", "Castanopsis fargesii", "Elaeocarpus sylvestris")) %>%
  dplyr::mutate(SpeciesFullName = as.factor(SpeciesFullName),
                SiteID = as.factor(SiteID)) %>%
  dplyr::select(SampleID,SpeciesFullName,GrowthForm,SiteID,LMA,LCC,LNC,LPC,Ld13C,Asat,Vcmax25)

summary(intrasp_traitdata)

# extract SD, mean and CV for LMA,LCC,LNC,LPC,Asat,Vcmax25, 2 sites
intrasp_traitdata_summary = intrasp_traitdata %>%
  dplyr::group_by(SpeciesFullName, SiteID) %>%
  dplyr::summarise(
    LMA_mean = mean(LMA, na.rm = TRUE), LMA_sd = sd(LMA, na.rm = TRUE),
    LMA_CV = (LMA_sd / LMA_mean) * 100,
    
    LCC_mean = mean(LCC, na.rm = TRUE), LCC_sd = sd(LCC, na.rm = TRUE),
    LCC_CV = (LCC_sd / LCC_mean) * 100,
    
    LNC_mean = mean(LNC, na.rm = TRUE), LNC_sd = sd(LNC, na.rm = TRUE),
    LNC_CV = (LNC_sd / LNC_mean) * 100,
    
    LPC_mean = mean(LPC, na.rm = TRUE), LPC_sd = sd(LPC, na.rm = TRUE),
    LPC_CV = (LPC_sd / LPC_mean) * 100,
    
    Ld13C_mean = mean(Ld13C, na.rm = TRUE), Ld13C_sd = sd(Ld13C, na.rm = TRUE),
    Ld13C_CV_abs = abs((Ld13C_sd / Ld13C_mean) * 100),
    
    Asat_mean = mean(Asat, na.rm = TRUE), Asat_sd = sd(Asat, na.rm = TRUE),
    Asat_CV = (Asat_sd / Asat_mean) * 100,
    
    Vcmax25_mean = mean(Vcmax25, na.rm = TRUE), Vcmax25_sd = sd(Vcmax25, na.rm = TRUE),
    Vcmax25_CV = (Vcmax25_sd / Vcmax25_mean) * 100
  )

intrasp_traitdata_cv = intrasp_traitdata_summary %>%
  dplyr::select(SpeciesFullName, SiteID, LMA_CV, LCC_CV, LNC_CV, LPC_CV, Ld13C_CV_abs, Asat_CV, Vcmax25_CV) %>%
  pivot_longer(cols = -c(SpeciesFullName, SiteID), names_to = "Trait", values_to = "CV")

# grouped barplot for trait CV between sites
plt_intrasp_cv = intrasp_traitdata_cv %>%
  ggplot(aes(x = SpeciesFullName, y = CV, fill = SiteID)) +
  geom_bar(stat = "identity", position = position_dodge(), color="black") +
  facet_wrap(~ Trait, scales = "free_y") +
  labs(#title = "Intra-specific trait variation (CV) between sites",
    x = "Species", y = "Coefficient of Variation (%)") +
  theme_minimal() +
  scale_fill_manual(values = c("hilltop" = "#547bb4", "valley" = "#dd7c4f")) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 0.95))

ggsave(plt_intrasp_cv, filename = "outplts/SIFigures/FigS8-intrasp_leaftrait_cv.pdf", 
       width = 5.2, height = 5.5)
