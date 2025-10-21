rm(list=ls())
library(tidyverse)
library(dplyr)
library(vegan)
library(ggpmisc) # stat_poly_eq function
library(reshape2)
library(patchwork)

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

# Part 1: data preparation -----------------------------------------------------
#
### SET DATA TO USE ###
traitDataPCA_touse = traitDataIndv_spgfavg_log
### SET TRAITS TO USE###
# All traits, need to be log transformed and scaled 
traitName_touse = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA","LPC")

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


# Part 2 : define functions for non-phylo PCA calculation ----------------------
func_nonphyPCA_stepwise <- function(var_traitName_touse_subset, var_nonphyloPCAData_numonly){
  
  # Part 2.1: eigen decomposition ----------------------------------------------
  # 
  # For each eigendecomposition :
  # - Extract the main axis of trait variation defining the trait space 
  # (i.e. Reduced trait space in  the text);
  # - Estimated trait loadings, which reflect the correlation between traits and relevant dimension of trait variation;
  # - We estimated the correlations between trait loadings considering all relevant dimensions in the space (i.e. Reduced space) by extracting the angles that each pair of traits form in the reduced space;
  # - Extract the Effective Number of Dimensions (END) 
  
  # the pairwise correlation between traits
  pwcorr = cor(var_nonphyloPCAData_numonly[, var_traitName_touse_subset], use = "pairwise.complete")
  
  AllTraits = eigen(pwcorr)
  dimnames(AllTraits$vectors) = list(rownames(pwcorr), paste0("PC", 1:ncol(AllTraits$vectors)))
  # Variance vector/total variance = cumulative variance explained by dimensions.
  AllTraits$CumVariance = cumsum(AllTraits$values / sum(AllTraits$values)) 
  # see how many dimensions are needed 
  AllTraits$values # 3 dimensions with eigenval > 1
  dimensions = 3 # TODO: modify this according to threshold
  AllTraits$dimensions = dimensions
  
  # Part 2.2: create matrix of loadings (eigenvectors * sqrt(eigenvalues)) -----
  # 
  # The loadings matrix is the product of the eigenvectors and the square root of the eigenvalues.
  Eival = AllTraits$values
  loadingsEig = AllTraits$vectors
  
  for(i in 1:ncol(loadingsEig)){
    loadingsEig[, i] = loadingsEig[, i]*sqrt(Eival[i])
  }
  AllTraits$loadings = loadingsEig
  
  # Part 2.3: calculate angles between pairs of vectors ------------------------
  # 
  loadEigVar = AllTraits$loadings[, 1:dimensions]
  AllTraits$angles <- matrix(NA, nrow = nrow(loadEigVar), ncol = nrow(loadEigVar), 
                             dimnames = list(rownames(loadEigVar), rownames(loadEigVar)))
  for(i in 1:nrow(loadEigVar)){
    for(j in 1:nrow(loadEigVar)){
      AllTraits$angles[i, j] <- func_calc_angle_pcavars(loadEigVar[i, ], loadEigVar[j, ])
    }
  }
  round(AllTraits$angles, 2)
  
  # Part 2.4: extract the effective number of dimensions (END) -----------------
  AllTraits$EffDimens <- vegan::diversity(x = AllTraits$values, index = "invsimpson")
  
  # Return curated PCA result
  return(AllTraits)
}

# Part 3: perform stepwise PCA angle calculation (non-phylo) -------------------
# 
# RES + LES only
PCAresult_RESLES = func_nonphyPCA_stepwise(
  c("RTD","SRL","RD","RNC","LMA","LNC"), nonphyloPCAData_numonly)
# RES + LES + RDMC
PCAresult_RESLES_RDMC = func_nonphyPCA_stepwise(
  c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA"), nonphyloPCAData_numonly)
# RES + LES + LPC
PCAresult_RESLES_LPC = func_nonphyPCA_stepwise(
  c("RTD","SRL","RD","RNC","LMA","LNC","LPC"), nonphyloPCAData_numonly)
# RES + LES + RDMC + LPC
PCAresult_RESLES_RDMCLPC = func_nonphyPCA_stepwise(
  c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA","LPC"), nonphyloPCAData_numonly)

# Aux combinations, for END comparison, RES & RES+LES
PCAresult_RES = func_nonphyPCA_stepwise(
  c("RTD","SRL","RD","RNC"), nonphyloPCAData_numonly)
PCAresult_RES_LMA = func_nonphyPCA_stepwise(
  c("RTD","SRL","RD","RNC","LMA"), nonphyloPCAData_numonly)
PCAresult_RES_LNC = func_nonphyPCA_stepwise(
  c("RTD","SRL","RD","RNC","LNC"), nonphyloPCAData_numonly)
PCAresult_RES_LPC = func_nonphyPCA_stepwise(
  c("RTD","SRL","RD","RNC","LPC"), nonphyloPCAData_numonly)

# Remove irrelevant data
rm(list=setdiff(ls(), c("PCAresult_RESLES", "PCAresult_RESLES_RDMC", "PCAresult_RESLES_LPC", "PCAresult_RESLES_RDMCLPC", "PCAresult_RES", "PCAresult_RES_LMA", "PCAresult_RES_LNC", "PCAresult_RES_LPC")))
# Save data for phylo/nonphylo comparisons
save(PCAresult_RESLES,PCAresult_RESLES_RDMC,PCAresult_RESLES_LPC,PCAresult_RESLES_RDMCLPC,
     PCAresult_RES, PCAresult_RES_LMA, PCAresult_RES_LNC, PCAresult_RES_LPC,
     file = "nonphyloPCAresult_RESLES-step31.RData")


# Part 4: angle correlation plot, non phylo, alternative predictor vars --------
#
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


# varAngleMelted_combined = varAngleMelted_uptri_RESLES %>% 
#   right_join(varAngleMelted_uptri_RESLES_RDMC, by=join_by(trait1,trait2)) %>%
#   right_join(varAngleMelted_uptri_RESLES_LPC, by=join_by(trait1,trait2)) %>%
#   right_join(varAngleMelted_uptri_RESLES_RDMCLPC, by=join_by(trait1,trait2)) %>%
#   na.omit() %>% # exclude variables that can't form a point in scatter plot.
#   dplyr::mutate(pair_name = paste0(trait1,"-",trait2)) %>%
#   dplyr::mutate(pair_annotation = ifelse(
#     pair_name == "SRL-RD", "RES collaboration",
#     ifelse(pair_name == "RTD-RNC", "RES convervation",
#            ifelse(pair_name == "LMA-LNC", "LES", "other")))
#   )

# Combine the dataframes into a list
varAngleMelted_dflist = list(
  varAngleMelted_uptri_RESLES, 
  varAngleMelted_uptri_RESLES_RDMC, 
  varAngleMelted_uptri_RESLES_LPC, 
  varAngleMelted_uptri_RESLES_RDMCLPC
)
# Perform an outer join on all dataframes using reduce
varAngleMelted_combined = reduce(varAngleMelted_dflist, function(x, y) merge(x, y, by = c("trait1", "trait2"), all = TRUE))

varAngleMelted_combined = varAngleMelted_combined %>%
  dplyr::mutate(pair_name = paste0(trait1,"-",trait2)) %>%
  dplyr::mutate(pair_annotation = ifelse(
    pair_name == "SRL-RD", "RES collaboration",
    ifelse(pair_name %in% c("RTD-RNC"), "RES conservation",
           ifelse(pair_name %in% c("LMA-LNC"), "LES", "other")))
  )


# plot1: raw-raw+RDMC ----------------------------------------------------------
plt_anglecorr_addRDMC = ggplot(data=varAngleMelted_combined,
                               aes(x=angle_LESRES, y=angle_LESRES_RDMC)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone", y="Backbone + RES co-pred") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

# plot2: raw-raw+LPC
plt_anglecorr_addLPC = ggplot(data=varAngleMelted_combined,
                              aes(x=angle_LESRES, y=angle_LESRES_LPC)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone", y="Backbone + LES co-pred") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

# plot3: raw-raw+RDMCLPC
plt_anglecorr_addRDMCLPC = ggplot(data=varAngleMelted_combined,
                                  aes(x=angle_LESRES, y=angle_LESRES_RDMCLPC)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone", y="Backbone + all co-pred") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

plt_anglecorr_nonphy_total = plt_anglecorr_addRDMC + plt_anglecorr_addLPC + plt_anglecorr_addRDMCLPC + plot_layout(guides = "collect") & theme(legend.position='bottom', legend.title = element_blank())

ggsave(plot = plt_anglecorr_nonphy_total, filename = "plt_anglecorr_nonphy_total.pdf", width = 8, height= 3.2)
# ggsave(plot = plt_anglecorr_nonphy_total, filename = "plt_anglecorr_nonphy_total.png", width = 8, height= 3.2, dpi=300)

