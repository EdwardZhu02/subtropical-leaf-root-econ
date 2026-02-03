rm(list=ls())
library(tidyverse)
library(dplyr)
library(vegan)
library(ggpmisc) # stat_poly_eq function
library(reshape2)
library(patchwork)

# phylo-PCA dependencies
library(ape) # general phylogenetic analysis
library(geiger) # compare taxa in data and tree
library(treeplyr) # general phylogenetic analysis
library(phytools) # general phylogenetic analysis, phyl.pca

load("species-level-code/traitDataFujian-spavg-phylo-step2.RData")

# Part 1: data preparation (non-phylo) -----------------------------------------
# 
### SET DATA TO USE ###
traitDataPCA_touse = traitDataIndv_spgfavg_log
### SET TRAITS TO USE###
# All traits, need to be log transformed and scaled 
traitName_touse = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","LPC")

# Resolve problem: non-unique values when setting 'row.names': ‘Castanopsis fordii’, ‘Machilus pauhoi’ because of sp-gf averaging
traitDataPCA_touse = traitDataPCA_touse %>% mutate(
  PCAIdentifier_spgf = paste(gsub(" ","_",speciesFullName), GrowthForm, sep="-")
)

# Curate data for PCA analysis
nonphyloPCAData_numonly = traitDataPCA_touse %>% ungroup() %>% 
  dplyr::select(PCAIdentifier_spgf,speciesFullName,family_name,GrowthForm,all_of(traitName_touse)) %>% 
  na.omit()
rownames(nonphyloPCAData_numonly) = nonphyloPCAData_numonly$PCAIdentifier_spgf

# scale traits based on columns (mean=0, SD=1, same as z-transform) 
nonphyloPCAData_numonly[, traitName_touse] = apply(nonphyloPCAData_numonly[, traitName_touse], 2, scale)


# Part 2: data preparation (phylo) ---------------------------------------------
# 
### SET DATA TO USE ###
traitData_touse = traitDataIndv_spavg_log # perform scaling later (line 71, 7/12-24)
### SET TRAITS TO USE###
# All traits, need to be log transformed and scaled 
traitName_touse = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","LPC")


# Convert species Latin name `A B` to `A_B` to match with tree `tip_label`
traitData_touse = traitData_touse %>% dplyr::mutate(tipLabelMatched=gsub(" ", "_", speciesFullName)) %>%
  dplyr::select(tipLabelMatched, everything())

# Extract phylo-tree constructed in step2
# phylotree_result[[1]] is the resulting phylogenetic tree (type: phylo).
phylotree_touse = phylotree_result[[1]] 
rm(phylotree_result)

# Combine tree and data
traitData_touse_PhyloDataCombined = make.treedata(
  tree = phylotree_touse, data = traitData_touse, name_column = "tipLabelMatched")

# Save data and tree separately, facilitating downstream analyses.
traitData_touse_PhyloDataCombined_data = as.data.frame(traitData_touse_PhyloDataCombined$dat) %>% 
  dplyr::mutate(tipLabelMatched=gsub(" ", "_", speciesFullName)) %>%
  dplyr::select(tipLabelMatched, everything())
traitData_touse_PhyloDataCombined_tree = traitData_touse_PhyloDataCombined$phy

# set tree node label to NULL to avoid Error: Labels duplicated between tips and nodes in phylogeny
# Ref: https://stackoverflow.com/questions/51261388/duplicate-tips-label-removal-importing-phylogenetic-tree-in-r-for-comparison
traitData_touse_PhyloDataCombined_tree$node.label = NULL

phyloPCAData_numonly = traitData_touse %>% ungroup() %>% 
  dplyr::select(tipLabelMatched,GrowthForm,all_of(traitName_touse)) %>% 
  na.omit()
phyloPCAData_numonly_rownames = phyloPCAData_numonly$tipLabelMatched

# used for growth form annotation in PCA biplot
phyloPCAData_meta_rmna = phyloPCAData_numonly
phyloPCAData_numonly = phyloPCAData_numonly %>% dplyr::select(-tipLabelMatched, -GrowthForm)

# scale traits based on columns (mean=0, SD=1, same as z-transform) 
phyloPCAData_numonly[, traitName_touse] = apply(phyloPCAData_numonly[, traitName_touse], 2, scale)
rownames(phyloPCAData_numonly) = phyloPCAData_numonly_rownames


# Part 3: Random trait in non-phylo space --------------------------------------
traitData_all_nonPhylo = nonphyloPCAData_numonly
traitData_all_Phylo = phyloPCAData_numonly

traitName_RES  = c("RTD","SRL","RD","RNC")
traitData_RES = traitData_all_nonPhylo[, traitName_RES] %>% na.omit()
rownames(traitData_RES) = traitData_all_nonPhylo$PCAIdentifier_spgf

traitName_RES_phylo = c("RTD","SRL","RD","RNC")
traitData_RES_phylo = traitData_all_Phylo[, traitName_RES_phylo] %>% na.omit()
rownames(traitData_RES_phylo) = phyloPCAData_numonly_rownames

traitName_RLES  = c("RTD","SRL","RD","RNC","LMA","LNC")
traitData_RLES = traitData_all_nonPhylo[, traitName_RLES] %>% na.omit()
rownames(traitData_RLES) = traitData_all_nonPhylo$PCAIdentifier_spgf

traitName_RLES_phylo = c("RTD","SRL","RD","RNC","LMA","LNC")
traitData_RLES_phylo = traitData_all_Phylo[, traitName_RLES_phylo] %>% na.omit()
rownames(traitData_RLES_phylo) = phyloPCAData_numonly_rownames


rep = 1000 #Number of REPETITIONS
traitsRand =  c("Rand_RLES","Rand_RES")
END_PCA_Rand = matrix(NA, nrow = rep, ncol = length(traitsRand), 
                       dimnames = list(paste0("Rep_", 1:rep), traitsRand)) # To store the results

# RLES, non-phylo (every single run has a different EOD, accounting for ramdom vars)
for(i in 1:rep){
  cat(paste("\r", "Rep: ", i, "out of", rep, "\r"))
  
  traitData_RLES[, "Rand"] = rnorm(nrow(traitData_RLES), 0 , sd = 1)
  PCA_aux = princomp(traitData_RLES) # already log-transformed and scaled
  END_PCA_Rand[i , "Rand_RLES"] = vegan::diversity(x = apply(PCA_aux$scores, 2, var), index = "invsimpson")
  
  traitData_RES[, "Rand"] = rnorm(nrow(traitData_RES), 0 , sd = 1)
  PCA_aux = princomp(traitData_RES) # already log-transformed and scaled
  END_PCA_Rand[i , "Rand_RES"] = vegan::diversity(x = apply(PCA_aux$scores, 2, var), index = "invsimpson")
}

# RLES, phylo (every single run has a same EOD, accounting for ramdom vars)
traitData_RLES_phylo[, "Rand"] = rnorm(nrow(traitData_RLES_phylo), 0 , sd = 1)
rownames(traitData_RLES_phylo) = phyloPCAData_numonly_rownames
Phyl_PCA_aux1 = phyl.pca(
  traitData_touse_PhyloDataCombined_tree, traitData_RLES_phylo,
  mode="corr", method='lambda') # already log-transformed and scaled
END_PCA_Rand_RLES_phylo = vegan::diversity(x = diag(Phyl_PCA_aux1$Eval), index = "invsimpson")

# RES, phylo (every single run has a same EOD, accounting for ramdom vars)
traitData_RES_phylo[, "Rand"] = rnorm(nrow(traitData_RES_phylo), 0 , sd = 1)
rownames(traitData_RES_phylo) = phyloPCAData_numonly_rownames
Phyl_PCA_aux2 = phyl.pca(
  traitData_touse_PhyloDataCombined_tree, traitData_RES_phylo,
  mode="corr", method='lambda') # already log-transformed and scaled
END_PCA_Rand_RES_phylo = vegan::diversity(x = diag(Phyl_PCA_aux2$Eval), index = "invsimpson")


# Part 3: Extracting P-values --------------------------------------------------
# Load PCA results, phyl and non-phyl
load("nonphyloPCAresult_RESLES-step31.RData")
load("phyloPCAresult_RESLES-step32.RData")
nonPhyl_PCA = list(
  PCAresult_RESLES=PCAresult_RESLES,
  PCAresult_RESLES_RDMC=PCAresult_RESLES_RDMC,
  PCAresult_RESLES_LPC=PCAresult_RESLES_LPC,
  PCAresult_RESLES_RDMCLPC=PCAresult_RESLES_RDMCLPC,
  
  PCAresult_RES=PCAresult_RES,
  PCAresult_RES_LMA=PCAresult_RES_LMA,
  PCAresult_RES_LNC=PCAresult_RES_LNC,
  PCAresult_RES_LPC=PCAresult_RES_LPC
  ) # named list
# Phyl_PCA = list(
#   phyloPCAresult_RESLES=phyloPCAresult_RESLES,
#   phyloPCAresult_RESLES_RDMC=phyloPCAresult_RESLES_RDMC,
#   phyloPCAresult_RESLES_LPC=phyloPCAresult_RESLES_LPC,
#   phyloPCAresult_RESLES_RDMCLPC=phyloPCAresult_RESLES_RDMCLPC) # named list

names(nonPhyl_PCA)
Comp = c("Random_Obs", names(nonPhyl_PCA))

END_PCA_Mean_andPval = matrix(NA, ncol = length(Comp), nrow = length(traitsRand),
                               dimnames = list(traitsRand, Comp))
for(i in 1:nrow(END_PCA_Mean_andPval)){
  
  END_PCA_Mean_andPval[i, "Random_Obs"] = mean(END_PCA_Rand[,rownames(END_PCA_Mean_andPval)[i]])
  if(i == 1){
    auxnames <- names(nonPhyl_PCA)[1:4]
  } else if (i == 2){
    auxnames <- names(nonPhyl_PCA)[5:8]
  }
  for(j in 1:length(auxnames)){
    aux = nonPhyl_PCA[[auxnames[j]]]
    END_PCA_Mean_andPval[rownames(END_PCA_Mean_andPval)[i], auxnames[j]] = dnorm(
      aux$EffDimens, mean = END_PCA_Mean_andPval[i, "Random_Obs"], 
      sd = sd(END_PCA_Rand[,rownames(END_PCA_Mean_andPval)[i]]))
  }
}
END_PCA_Mean_andPval
## Random_Obs column contains the END expected under the inclusion of Root size traits

# Part 4: Draw the EOD plot ----------------------------------------------------
#
ENDdata_total_list = list(
  nonphy_RLES_Rand = END_PCA_Mean_andPval[1,1],
  nonphy_RES_Rand = END_PCA_Mean_andPval[2,1],
  phy_RLES_Rand = END_PCA_Rand_RLES_phylo,
  phy_RES_Rand = END_PCA_Rand_RES_phylo,
  
  nonphy_RLES = nonPhyl_PCA$PCAresult_RESLES$EffDimens,
  nonphy_RLES_RDMC = nonPhyl_PCA$PCAresult_RESLES_RDMC$EffDimens,
  nonphy_RLES_LPC = nonPhyl_PCA$PCAresult_RESLES_LPC$EffDimens,
  nonphy_RLES_RDMCLPC = nonPhyl_PCA$PCAresult_RESLES_RDMCLPC$EffDimens,
  
  nonphy_RES = nonPhyl_PCA$PCAresult_RES$EffDimens,
  nonphy_RES_LMA = nonPhyl_PCA$PCAresult_RES_LMA$EffDimens,
  nonphy_RES_LNC = nonPhyl_PCA$PCAresult_RES_LNC$EffDimens,
  nonphy_RES_LPC = nonPhyl_PCA$PCAresult_RES_LPC$EffDimens,
  
  phy_RLES = phyloPCAresult_RESLES$EffDimens,
  phy_RLES_RDMC = phyloPCAresult_RESLES_RDMC$EffDimens,
  phy_RLES_LPC = phyloPCAresult_RESLES_LPC$EffDimens,
  phy_RLES_RDMCLPC = phyloPCAresult_RESLES_RDMCLPC$EffDimens,
  
  phy_RES = phyloPCAresult_RES$EffDimens,
  phy_RES_LMA = phyloPCAresult_RES_LMA$EffDimens,
  phy_RES_LNC = phyloPCAresult_RES_LNC$EffDimens,
  phy_RES_LPC = phyloPCAresult_RES_LPC$EffDimens
)
rm(list=setdiff(ls(), c("ENDdata_total_list"))) # clear all workspace

ENDdata_total_df_long = as.data.frame(ENDdata_total_list) %>%
  pivot_longer(cols = everything(), names_to = "combination", values_to = "END") %>%
  mutate(comb_group = ifelse(grepl("nonphy", combination), "non-phylo", "phylo")) %>%
  mutate(combination = gsub("nonphy_", "", combination)) %>%
  mutate(combination = gsub("phy_", "", combination)) %>%
  filter(!(combination %in% c("RLES_Rand", "RES_Rand")))


ENDdata_total_df_long_plt1 = ENDdata_total_df_long %>%
  # exclude combinations containing "RES"
  filter(!(combination %in% c("RES","RES_Rand","RES_LMA","RES_LNC","RES_LPC")))

ENDdata_total_df_long_plt2 = ENDdata_total_df_long %>%
  # select combinations containing "RES"
  filter(combination %in% c("RES","RES_Rand","RES_LMA","RES_LNC","RLES"))

plt_END_RLESaltpred = ggplot(ENDdata_total_df_long_plt1, aes(x = comb_group,y = END,fill = combination))+
  geom_hline(yintercept = 0, linewidth=0.2, color="gray") +
  geom_hline(yintercept = 1, linewidth=0.2, color="gray") +
  geom_hline(yintercept = 2, linewidth=0.2, color="gray") +
  geom_hline(yintercept = 3, linewidth=0.2, color="gray") +
  geom_hline(yintercept = 4, linewidth=0.2, color="gray") +
  
  geom_bar(stat = "identity",width = 0.6,position = "dodge",color="black",linewidth=0.2) +
  geom_segment(aes(x = 0.5, xend = 1.4, y = ENDdata_total_list$nonphy_RLES_Rand, yend = ENDdata_total_list$nonphy_RLES_Rand), linewidth=0.7, color="darkred", linetype = "dashed") +
  geom_segment(aes(x = 1.6, xend = 2.5, y = ENDdata_total_list$phy_RLES_Rand, yend = ENDdata_total_list$phy_RLES_Rand), linewidth=0.7, color="darkblue", linetype = "dashed") +
  geom_text(aes(x=1, y=5.2, label = paste0(
    "END (Backbone+Random) = ",
    signif(ENDdata_total_list$nonphy_RLES_Rand, digits = 3))), 
    fontface = "bold", size = 3.5, angle=0) +
  geom_text(aes(x=2, y=5.2, label = paste0(
    "END (Backbone+Random) = ",
    signif(ENDdata_total_list$phy_RLES_Rand, digits = 3))), 
    fontface = "bold", size = 3.5, angle=0) +
  geom_text(aes(label = signif(END, digits = 3)), fontface = "bold",
            position=position_dodge(width = 0.65), vjust=-0.4, size = 3.5, angle=0) +
  scale_fill_manual(values = c("#1999B2","#95BCE5","#E84445","#F39DA0")) +
  #scale_fill_manual(values = c("#956642", "#678499", "#74070E", "#efe1a9")) +
  labs(y="Effective Number of Dimensions (END)") +
  theme_classic() +
  theme(axis.title.x = element_blank(), legend.title = element_blank())

ggsave(plot=plt_END_RLESaltpred , filename = "plt_END_RLESaltpred.pdf", width = 7.5, height = 3.5)
#ggsave(plot=plt_END_RLESaltpred , filename = "plt_END_RLESaltpred.png", width = 6.5, height = 3, dpi=300)


plt_END_RESLESorthog = ggplot(ENDdata_total_df_long_plt2, aes(x = comb_group,y = END,fill = combination))+
  geom_hline(yintercept = 0, linewidth=0.2, color="gray") +
  geom_hline(yintercept = 1, linewidth=0.2, color="gray") +
  geom_hline(yintercept = 2, linewidth=0.2, color="gray") +
  geom_hline(yintercept = 3, linewidth=0.2, color="gray") +
  
  geom_bar(stat = "identity",width = 0.6,position = "dodge",color="black",linewidth=0.2) +
  geom_segment(aes(x = 0.5, xend = 1.4, y = ENDdata_total_list$nonphy_RES_Rand, yend = ENDdata_total_list$nonphy_RES_Rand), linewidth=0.7, color="darkred", linetype = "dashed") +
  geom_segment(aes(x = 1.6, xend = 2.5, y = ENDdata_total_list$phy_RES_Rand, yend = ENDdata_total_list$phy_RES_Rand), linewidth=0.7, color="darkblue", linetype = "dashed") +
  geom_text(aes(x=1, y=4.8, label = paste0(
    "END (RES+Random) = ",
    signif(ENDdata_total_list$nonphy_RES_Rand, digits = 3))), 
    fontface = "bold", size = 3.5, angle=0) +
  geom_text(aes(x=2, y=4.8, label = paste0(
    "END (RES+Random) = ",
    signif(ENDdata_total_list$phy_RES_Rand, digits = 3))), 
    fontface = "bold", size = 3.5, angle=0) +
  geom_text(aes(label = signif(END, digits = 3)), fontface = "bold",
            position=position_dodge(width = 0.65), vjust=-0.8, size = 3.5, angle=0) +
  scale_fill_manual(values = c("#1999B2","#95BCE5","#E84445","#F39DA0")) +
  #scale_fill_manual(values = c("#956642", "#678499", "#74070E", "#efe1a9")) +
  labs(y="Effective Number of Dimensions (END)") +
  theme_classic() +
  theme(axis.title.x = element_blank(), legend.title = element_blank())

ggsave(plot=plt_END_RESLESorthog, filename = "species-level-code/plt_END_RESLESorthog.pdf", width = 7, height = 3.5)
#ggsave(plot=plt_END_RESLESorthog, filename = "plt_END_RESLESorthog.png", width = 6.5, height = 3, dpi=300)
