rm(list=ls())
library(tidyverse)
library(dplyr)
library(ggsignif)  # For adding significance markers
library(broom)     # For tidying the t-test results
library(ggpubr)
library(patchwork) # for merging plots
library(plotly)
library(ggsignif)  # For adding significance markers
library(rstatix)   # For pairwise wilcox test

# Perform across-site comparison using community-mean trait values
load("traitDataFujian-spavg-phylo-step2.RData")

# use data frame 'traitDataIndv', the original, not species averaged or normalized data for analysis
traitDataIndv_touse = traitDataIndv # backup original df, while copy it for this script


# Step 1: add SiteLoc column, and filter observations into different DFs -------
traitDataIndv_touse$SiteLoc = NA
traitDataIndv_touse$SiteLoc = ifelse(startsWith(traitDataIndv_touse$SpeciesID, "SH"), "hilltop", 
                                     ifelse(startsWith(traitDataIndv_touse$SpeciesID, "WH"), "valley",
                                            traitDataIndv_touse$SiteLoc)) # NA

# Step 2: reshape data for plotting (box plot, long data) ----------------------
traitDataIndv_gfloc_mean_box = traitDataIndv_touse %>%
  dplyr::select(where(is.numeric), GrowthForm, SiteLoc, SampleID) %>%   # Select numeric columns and SiteLoc
  pivot_longer(cols = -c(GrowthForm, SiteLoc, SampleID), names_to = "trait", values_to = "value") %>% 
  group_by(trait) %>%
  filter(!is.na(value))

# Exclude liana entries
traitDataIndv_gfloc_mean_box = traitDataIndv_gfloc_mean_box %>% dplyr::filter(GrowthForm != "liana")


# Step 3: generate plots -------------------------------------------------------
#
# Filter only data of interest to generate the plot
# traitlist_touse = c(
#   "RTD","SRL","RD","RNC","SRA","RPC","RCC","RDMC","SRR25",
#   "rs25","LMA","LNC","LPC","Rdark25P","Rlight25P","Vcmax25","Asat",
#   "vH","SWCleaf","SWCbranch","WD","Ks","TLP","H","DBH" # hydraulic traits
# )
# 1. photosynthesis and respiration
# traitlist_touse = c("SRR25","rs25","Rdark25P","Rlight25P","Vcmax25","Asat")
# 2. chemical traits
# traitlist_touse = c("RCC","RNC","RPC","LCC","LNC","LPC")
# 3. morphological, leaf+root
# traitlist_touse = c("RTD","SRL","RD","SRA","RDMC","LMA")
# 4. hydraulic traits
traitlist_touse = c( "vH","SWCleaf","SWCbranch","WD","Ks","TLP","H","DBH")

tomerge_plot_list = list() # empty list for plot merging

for (i in 1:length(traitlist_touse)) {
  
  traitDataIndv_gfloc_touse = traitDataIndv_gfloc_mean_box %>% 
    dplyr::filter(trait == traitlist_touse[i]) %>%
    dplyr::mutate(UniqIdentifier = paste(SiteLoc, GrowthForm, sep="_")) %>%
    dplyr::mutate(
      SiteLoc = as.factor(SiteLoc),
      GrowthForm = as.factor(GrowthForm),
      UniqIdentifier = as.factor(UniqIdentifier),
      trait = as.factor(trait)
      )
  
  # perform wilcox test
  traitData_gfloc.test.SiteGF = traitDataIndv_gfloc_touse %>% group_by(trait) %>%
    pairwise_wilcox_test(value ~ UniqIdentifier, paired = F, p.adjust.method = "holm")
  traitData_gfloc.test.Site = traitDataIndv_gfloc_touse %>% group_by(trait) %>%
    pairwise_wilcox_test(value ~ SiteLoc, paired = F, p.adjust.method = "holm")
  traitData_gfloc.test.GF = traitDataIndv_gfloc_touse %>% group_by(trait) %>%
    pairwise_wilcox_test(value ~ GrowthForm, paired = F, p.adjust.method = "holm")
  
  # plot1: site + growth form, individual significance check boxplot ===========
  plt1 = ggboxplot(traitDataIndv_gfloc_touse, x = "GrowthForm", y = "value", fill = "SiteLoc") +
    scale_fill_manual(values = c("#6888F5", "#D77071")) + labs(subtitle = traitlist_touse[i])
  traitData_gfloc.test.SiteGF = traitData_gfloc.test.SiteGF %>% add_xy_position() %>% 
    mutate(xmin = 0.6 * xmin, xmax = 0.6 * xmax)
  plt1 = plt1 + stat_pvalue_manual(traitData_gfloc.test.SiteGF, label = "p.adj.signif", step.increase = 0.04, hide.ns = T) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none",
          axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
  
  # plot2: site, significance check boxplot ====================================
  plt2 = ggboxplot(traitDataIndv_gfloc_touse, x = "SiteLoc", y = "value", fill = "lightgray")
  traitData_gfloc.test.Site = traitData_gfloc.test.Site %>% add_xy_position(x="SiteLoc") %>%
    mutate(xmin = 1, xmax = 2)
  plt2 = plt2 + stat_pvalue_manual(traitData_gfloc.test.Site, label = "p.adj.signif", step.increase = 0, vjust = 1.4, hide.ns = T) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
  
  # plot3: growth form, significance check boxplot ====================================
  plt3 = ggboxplot(traitDataIndv_gfloc_touse, x = "GrowthForm", y = "value", fill = "lightgray")
  traitData_gfloc.test.GF = traitData_gfloc.test.GF %>% add_xy_position(x="GrowthForm") %>%
    mutate(xmin = 1, xmax = 2)
  plt3 = plt3 + stat_pvalue_manual(traitData_gfloc.test.GF, label = "p.adj.signif", step.increase = 0, vjust=1.4, hide.ns = T) +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
  
  pltall = (plt1 + (plt2/plt3)) + plot_layout(widths = c(3, 2))
  
  # add to the total plot list
  tomerge_plot_list[[i]] = pltall
  # break
}

# Combine plots with patchwork
layout <- "
ABCD
EFGH
"

boxplt_GrowthFormCombined = tomerge_plot_list[[1]]
for (i in 2:length(traitlist_touse)) {
  boxplt_GrowthFormCombined = boxplt_GrowthFormCombined / tomerge_plot_list[[i]]
}

boxplt_GrowthFormCombined = boxplt_GrowthFormCombined + 
  # plot_annotation(title = "Part 1: Photosynthesis and respiration traits") +
  # plot_annotation(title = "Part 2: Chemical traits (leaf + fine root)") +
  plot_annotation(title = "Part 4: Hydraulic traits") +
  plot_layout(guides = "collect", design = layout) & theme(legend.position='bottom')

ggsave(plot = boxplt_GrowthFormCombined, filename = "boxplt_CompNew_Part4.pdf",
       width = 12, height = 5.5)
