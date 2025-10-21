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

load("traitDataFujian-spavg-phylo-step2.RData")

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
     file = "phyloPCAresult_RESLES-step32.RData")

# Part 5: angle correlation plot, phylo and non-phylo approach -----------------
# Both data underwent log transform and scaling, so the comparison is valid.
# Load data from non-phylo PCA
load("nonphyloPCAresult_RESLES-step31.RData") # PCAresult_RESLES,PCAresult_RESLES_RDMC,PCAresult_RESLES_LPC,PCAresult_RESLES_RDMCLPC

# process angle matrix, extracting only the upper part
func_get_upper_tri = function(cormat){
  cormat[lower.tri(cormat)] = NA
  return(cormat)
}

# - non-phylo
varAngle_uptri_RESLES = func_get_upper_tri(PCAresult_RESLES[["angles"]])
varAngleMelted_uptri_RESLES = reshape2::melt(varAngle_uptri_RESLES, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_LESRES=value)

varAngle_uptri_RESLES_RDMC = func_get_upper_tri(PCAresult_RESLES_RDMC[["angles"]])
varAngleMelted_uptri_RESLES_RDMC = reshape2::melt(varAngle_uptri_RESLES_RDMC, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_LESRES_RDMC=value)

varAngle_uptri_RESLES_LPC = func_get_upper_tri(PCAresult_RESLES_LPC[["angles"]])
varAngleMelted_uptri_RESLES_LPC = reshape2::melt(varAngle_uptri_RESLES_LPC, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_LESRES_LPC=value)

varAngle_uptri_RESLES_RDMCLPC = func_get_upper_tri(PCAresult_RESLES_RDMCLPC[["angles"]])
varAngleMelted_uptri_RESLES_RDMCLPC = reshape2::melt(varAngle_uptri_RESLES_RDMCLPC, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_LESRES_RDMCLPC=value)

# - phylo
varAngle_uptri_phyloRESLES = func_get_upper_tri(phyloPCAresult_RESLES[["angles"]])
varAngleMelted_uptri_phyloRESLES = reshape2::melt(varAngle_uptri_phyloRESLES, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_phyloLESRES=value)

varAngle_uptri_phyloRESLES_RDMC = func_get_upper_tri(phyloPCAresult_RESLES_RDMC[["angles"]])
varAngleMelted_uptri_phyloRESLES_RDMC = reshape2::melt(varAngle_uptri_phyloRESLES_RDMC, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_phyloLESRES_RDMC=value)

varAngle_uptri_phyloRESLES_LPC = func_get_upper_tri(phyloPCAresult_RESLES_LPC[["angles"]])
varAngleMelted_uptri_phyloRESLES_LPC = reshape2::melt(varAngle_uptri_phyloRESLES_LPC, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_phyloLESRES_LPC=value)

varAngle_uptri_phyloRESLES_RDMCLPC = func_get_upper_tri(phyloPCAresult_RESLES_RDMCLPC[["angles"]])
varAngleMelted_uptri_phyloRESLES_RDMCLPC = reshape2::melt(varAngle_uptri_phyloRESLES_RDMCLPC, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_phyloLESRES_RDMCLPC=value)


# Combine the dataframes into a list
varAngleMelted_dflist = list(
  varAngleMelted_uptri_RESLES, 
  varAngleMelted_uptri_RESLES_RDMC, 
  varAngleMelted_uptri_RESLES_LPC, 
  varAngleMelted_uptri_RESLES_RDMCLPC,
  varAngleMelted_uptri_phyloRESLES,
  varAngleMelted_uptri_phyloRESLES_RDMC,
  varAngleMelted_uptri_phyloRESLES_LPC,
  varAngleMelted_uptri_phyloRESLES_RDMCLPC
  )
# Perform an outer join on all dataframes using reduce
varAngleMelted_combined = reduce(varAngleMelted_dflist, function(x, y) merge(x, y, by = c("trait1", "trait2"), all = TRUE))

varAngleMelted_combined = varAngleMelted_combined %>%
  dplyr::mutate(pair_name = paste0(trait1,"-",trait2)) %>%
  dplyr::mutate(pair_annotation = ifelse(
    pair_name %in% c("SRL-RD","RD-SRR25","SRL-SRR25","SRL-SRA","RD-SRA"), "RES collaboration",
    ifelse(pair_name %in% c("RTD-RNC","RNC-RDMC","RTD-RDMC"), "RES conservation",
           ifelse(pair_name %in% c("LMA-LNC","LMA-LPC"), "LES", "other")))
    )


# plot1: raw, nonphylo-phylo ---------------------------------------------------
plt_anglecorr_phy_RLES = ggplot(data=varAngleMelted_combined,
                               aes(x=angle_LESRES, y=angle_phyloLESRES)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone", y="Backbone (phylo)") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

# plot2: raw+RDMC, nonphylo-phylo ----------------------------------------------
plt_anglecorr_phy_RLESRDMC = ggplot(data=varAngleMelted_combined,
                                   aes(x=angle_LESRES_RDMC, y=angle_phyloLESRES_RDMC)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone + RES co-pred", y="Backbone + RES co-pred (phylo)") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

# plot3: raw+LPC, nonphylo-phylo -----------------------------------------------
plt_anglecorr_phy_RLESLPC = ggplot(data=varAngleMelted_combined,
                              aes(x=angle_LESRES_LPC, y=angle_phyloLESRES_LPC)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone + LES co-pred", y="Backbone + LES co-pred (phylo)") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

# plot3: raw+RDMC+LPC, nonphylo-phylo ------------------------------------------
plt_anglecorr_phy_RLESRDMCLPC = ggplot(data=varAngleMelted_combined,
                                  aes(x=angle_LESRES_RDMCLPC, y=angle_phyloLESRES_RDMCLPC)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone + all co-pred", y="Backbone + all co-pred (phylo)") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

plt_anglecorr_phy_total = plt_anglecorr_phy_RLES + plt_anglecorr_phy_RLESRDMC + plt_anglecorr_phy_RLESLPC + plt_anglecorr_phy_RLESRDMCLPC + plot_layout(guides = "collect") & theme(legend.position='bottom', legend.title = element_blank())

ggsave(plot = plt_anglecorr_phy_total, filename = "plt_anglecorr_phycomp_total.pdf", width = 5.5, height= 6)
# ggsave(plot = plt_anglecorr_phy_total, filename = "plt_anglecorr_phycomp_total.png", width = 5.5, height= 6, dpi=300)

