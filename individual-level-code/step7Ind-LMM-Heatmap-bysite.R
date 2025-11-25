rm(list=ls())

# Figure S2
# devtools::install_github("ggobi/ggally")
library(tidyverse)
library(dplyr)
library(GGally)
library(rstatix)
library(ggsci)


load("individual-level-code/traitDataFujian-Ind-step1.RData") # TODO: individual-level data

### SET DATA TO USE ###
traitData_touse = traitDataIndv_SelectedTraits_log

# TODO: (25.2.21) Resolve historical problem: nomination between "speciesFullName" and "SpeciesFullName"
traitData_touse = traitData_touse %>% dplyr::rename(speciesFullName=SpeciesFullName) %>%
  dplyr::rename(Rr25 = SRR25, Rd25 = Rdark25P) %>%
  dplyr::mutate(SiteIDAbbr = ifelse(SiteID=="hilltop", "H", "V")) # Abbreviate site names for better visualization

# Define the list of variables
variables = c("RTD","SRL","RD","RNC","RDMC","Rr25","SRA","RPC","RCC", # Fine root traits
              "LMA","LNC","LPC","LCC","Ld13C","Rd25","Vcmax25","Asat","LA" # Leaf traits
)


p3 = ggpairs(traitData_touse, columns = variables,
             mapping = ggplot2::aes(color = SiteIDAbbr, shape = GrowthForm),
             upper = NULL, #list(continuous = wrap("cor",method = "spearman")),
             lower = list(na ="na"),
             diag = NULL)+
  ggplot2::scale_color_manual(values = c("#547bb4", "#dd7c4f"))+
  ggplot2::theme_bw()

pdf("individual-level-code/LMM-all-heatmap.pdf",height = 15,width = 15)
print(p3)
dev.off()
