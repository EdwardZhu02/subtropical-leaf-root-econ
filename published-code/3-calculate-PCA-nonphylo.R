rm(list=ls())
library(tidyverse)
library(dplyr)
library(vegan)
library(ggpmisc) # stat_poly_eq function
library(reshape2)
library(patchwork)

load("outdata/traitDataFujian-SpAvg-phylo-step2.RData")
if (!dir.exists("outplts/SIFigures")) {
  dir.create("outplts/SIFigures", showWarnings = FALSE, recursive = TRUE)
}

# 1: non-phylo informed PCA results --------------------------------------------

### Function to estimate angles between pairs of vectors (2023, comments)
func_calc_angle_pcavars <- function(x, y){
  norm_x <- sqrt(sum(x^2)) #sd of eigenvector
  norm_y <- sqrt(sum(y^2))
  cosXY <- round(as.numeric((x %*% y) / (norm_x * norm_y)), 8) #dot product/magnitude of vectors
  # round introduced to avoid numerical problems with the acos function
  angle <- acos(cosXY) * 180 / pi
  return(angle)
}

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


# define functions for non-phylo PCA calculation
func_nonphyPCA_stepwise <- function(var_traitName_touse_subset, var_nonphyloPCAData_numonly){
  
  # eigen decomposition
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
  
  # create matrix of loadings (eigenvectors * sqrt(eigenvalues))
  # 
  # The loadings matrix is the product of the eigenvectors and the square root of the eigenvalues.
  Eival = AllTraits$values
  loadingsEig = AllTraits$vectors
  
  for(i in 1:ncol(loadingsEig)){
    loadingsEig[, i] = loadingsEig[, i]*sqrt(Eival[i])
  }
  AllTraits$loadings = loadingsEig
  
  # calculate angles between pairs of vectors
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
  
  # extract the effective number of dimensions (END)
  AllTraits$EffDimens <- vegan::diversity(x = AllTraits$values, index = "invsimpson")
  
  # Return curated PCA result
  return(AllTraits)
}

# perform stepwise PCA angle calculation (non-phylo)
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
     file = "outdata/nonphyloPCAresult_RESLES-step3.RData")
