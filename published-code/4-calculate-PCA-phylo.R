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

load("outdata/traitDataFujian-SpAvg-phylo-step2.RData")
if (!dir.exists("outplts/SIFigures")) {
  dir.create("outplts/SIFigures", showWarnings = FALSE, recursive = TRUE)
}


### Function to estimate angles between pairs of vectors (2023, comments)
func_calc_angle_pcavars <- function(x, y){
  norm_x <- sqrt(sum(x^2)) #sd of eigenvector
  norm_y <- sqrt(sum(y^2))
  cosXY <- round(as.numeric((x %*% y) / (norm_x * norm_y)), 8) #dot product/magnitude of vectors
  # round introduced to avoid numerical problems with the acos function
  angle <- acos(cosXY) * 180 / pi
  return(angle)
}

# Part 1: combine raw data and tree --------------------------------------------
# 
### SET DATA TO USE ###
traitData_touse = traitDataIndv_spavg_log # perform scaling later (line 71, 7/12-24)
### SET TRAITS TO USE###
# All traits, need to be log transformed and scaled 
traitName_touse = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA","LPC")


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


# Part 2: Curate data for PCA analysis (normalized data) -----------------------
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


# Part 3: define functions for phylo PCA calculation ----------------------
func_phyloPCA_stepwise <- function(var_traitName_touse_subset, var_traitData_touse_PhyloDataCombined_tree, var_phyloPCAData_numonly_total){
  
  var_phyloPCAData_numonly_rownames = rownames(var_phyloPCAData_numonly_total)
  var_phyloPCAData_numonly = var_phyloPCAData_numonly_total[,var_traitName_touse_subset]
  rownames(var_phyloPCAData_numonly) = var_phyloPCAData_numonly_rownames
  
  # Part 3.1: Perform phylo-informed PCA
  phyloPCA.results = phyl.pca(
    var_traitData_touse_PhyloDataCombined_tree, var_phyloPCAData_numonly,
    mode="corr", method='lambda')
  
  phyloPCA.results$lambda # lambda values
  diag(phyloPCA.results$Eval) # eigenvalues
  
  dimensions = 3
  phyloPCA.results$dimensions = dimensions
  
  # Part 3.2: calculate angles between pairs of vectors ------------------------
  phyloPCA.results$angles <- matrix(NA, nrow = ncol(var_phyloPCAData_numonly), ncol = ncol(var_phyloPCAData_numonly),
                                    dimnames = list(colnames(var_phyloPCAData_numonly), 
                                                    colnames(var_phyloPCAData_numonly)))
  for(i in 1:ncol(phyloPCA.results$angles)){
    for(j in 1:ncol(phyloPCA.results$angles)){
      phyloPCA.results$angles[i, j] <-  func_calc_angle_pcavars(phyloPCA.results$L[i, 1:dimensions], phyloPCA.results$L[j, 1:dimensions])
    }
  }
  colSums(phyloPCA.results$L**2)
  
  # Part 3.3: extract the effective number of dimensions (END) -----------------
  phyloPCA.results$EffDimens <- vegan::diversity(
    x = diag(phyloPCA.results$Eval), index = "invsimpson")
  
  # Return curated PCA result
  return(phyloPCA.results)
}

# Part 4: perform stepwise PCA angle calculation (phylo-informed) --------------
# 
# RES + LES only
phyloPCAresult_RESLES = func_phyloPCA_stepwise(
  c("RTD","SRL","RD","RNC","LMA","LNC"), 
  traitData_touse_PhyloDataCombined_tree, phyloPCAData_numonly)
diag(phyloPCAresult_RESLES$Eval) # eigenvalues, determine how much dimensions (>1)

# RES + LES + RDMC
phyloPCAresult_RESLES_RDMC = func_phyloPCA_stepwise(
  c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA"), 
  traitData_touse_PhyloDataCombined_tree, phyloPCAData_numonly)
diag(phyloPCAresult_RESLES_RDMC$Eval) # eigenvalues, determine how much dimensions (>1)

# RES + LES + LPC
phyloPCAresult_RESLES_LPC = func_phyloPCA_stepwise(
  c("RTD","SRL","RD","RNC","LMA","LNC","LPC"), 
  traitData_touse_PhyloDataCombined_tree, phyloPCAData_numonly)
diag(phyloPCAresult_RESLES_LPC$Eval) # eigenvalues, determine how much dimensions (>1)

# RES + LES + RDMC + LPC
phyloPCAresult_RESLES_RDMCLPC = func_phyloPCA_stepwise(
  c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA","LPC"), 
  traitData_touse_PhyloDataCombined_tree, phyloPCAData_numonly)
diag(phyloPCAresult_RESLES_RDMCLPC$Eval) # eigenvalues, determine how much dimensions (>1)

# Aux combinations, for END comparison, RES & RES+LES
phyloPCAresult_RES = func_phyloPCA_stepwise(
  c("RTD","SRL","RD","RNC"), traitData_touse_PhyloDataCombined_tree, phyloPCAData_numonly)
phyloPCAresult_RES_LMA = func_phyloPCA_stepwise(
  c("RTD","SRL","RD","RNC","LMA"), traitData_touse_PhyloDataCombined_tree, phyloPCAData_numonly)
phyloPCAresult_RES_LNC = func_phyloPCA_stepwise(
  c("RTD","SRL","RD","RNC","LNC"), traitData_touse_PhyloDataCombined_tree, phyloPCAData_numonly)
phyloPCAresult_RES_LPC = func_phyloPCA_stepwise(
  c("RTD","SRL","RD","RNC","LPC"), traitData_touse_PhyloDataCombined_tree, phyloPCAData_numonly)

# Remove irrelevant data
rm(list=setdiff(ls(), c("phyloPCAresult_RESLES", "phyloPCAresult_RESLES_RDMC", "phyloPCAresult_RESLES_LPC", "phyloPCAresult_RESLES_RDMCLPC","phyloPCAresult_RES","phyloPCAresult_RES_LMA","phyloPCAresult_RES_LNC","phyloPCAresult_RES_LPC")))
# Save data for downstream analyses
save(phyloPCAresult_RESLES,phyloPCAresult_RESLES_RDMC,phyloPCAresult_RESLES_LPC,phyloPCAresult_RESLES_RDMCLPC,
     phyloPCAresult_RES,phyloPCAresult_RES_LMA,phyloPCAresult_RES_LNC,phyloPCAresult_RES_LPC,
     file = "outdata/phyloPCAresult_RESLES-step3.RData")
