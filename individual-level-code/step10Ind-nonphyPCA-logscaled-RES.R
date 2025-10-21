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
library(factoextra) # perform PCA analysis
library(FactoMineR) # perform PCA analysis
library(patchwork) # plot merging
library(cowplot) # plot merging

load("2502-indv-level-code/traitDataFujian-Ind-step1.RData")

### SET DATA TO USE ### 
# TODO: perform scaling later (line 36, 27/12-24)
traitDataPCA_touse = traitDataIndv_SelectedTraits_log # perform scaling later (line 36, 27/12-24)

### SET TRAITS TO USE###
# All traits, need to be log transformed and scaled 
traitName_touse = c("RTD","SRL","RD","RNC")

# Resolve problem: non-unique values when setting 'row.names': ‘Castanopsis fordii’, ‘Machilus pauhoi’ 
# because of sp-gf averaging
traitDataPCA_touse = traitDataPCA_touse %>% mutate(
  PCAIdentifier_spgf = paste(gsub(" ","_",SpeciesFullName), GrowthForm, SampleID, sep="-")
)

# Curate data for PCA analysis
nonphyloPCAData_numonly = traitDataPCA_touse %>% ungroup() %>% 
  dplyr::select(PCAIdentifier_spgf,GrowthForm,all_of(traitName_touse)) %>% 
  na.omit()

# TODO: Added 27/12-24, concordant with angle calculation
# scale traits based on columns (mean=0, SD=1, same as z-transform) 
nonphyloPCAData_numonly[, traitName_touse] = apply(nonphyloPCAData_numonly[, traitName_touse], 2, scale)
# fix rownames
rownames(nonphyloPCAData_numonly) = nonphyloPCAData_numonly$PCAIdentifier_spgf

# used for growth form annotation in PCA biplot
nonphyloPCAData_meta_rmna = nonphyloPCAData_numonly
nonphyloPCAData_numonly = nonphyloPCAData_numonly %>% dplyr::select(-PCAIdentifier_spgf, -GrowthForm)


#-----------------------------------------------------------
# Perform PCA using FactoMineR (original result)
nonphyloPCAresult = PCA(nonphyloPCAData_numonly, scale.unit = TRUE, graph = FALSE)

#-----------------------------------------------------------
# Varimax Rotation Procedure
# Apply varimax rotation on the variable loadings from the PCA result
varimax_result <- varimax(nonphyloPCAresult$var$coord)
rotated_loadings <- as.matrix(varimax_result$loadings)

# Compute rotated individual scores by applying the rotation matrix
rotated_scores <- nonphyloPCAresult$ind$coord %*% varimax_result$rotmat

# Create a new PCA object by updating the original PCA result with rotated values.
nonphyloPCAresult_varimax_rotated <- nonphyloPCAresult
nonphyloPCAresult_varimax_rotated$var$coord <- rotated_loadings
nonphyloPCAresult_varimax_rotated$ind$coord <- rotated_scores  


# ------------------------------------------------------------------------------
# Scree plot: contributions of axes
# ------------------------------------------------------------------------------
nonphyloPCA.screeplot = fviz_screeplot(nonphyloPCAresult, addlabels = TRUE, ylim = c(0, 51),
                                       barfill="#A9082C", barcolor="black", linecolor="black") +
  labs(title="R-LES") + theme_classic()
ggsave(filename = "2502-indv-level-code/Ind_nonphyPCA_scree-RES.pdf", plot = nonphyloPCA.screeplot, width = 2.2, height = 2.5)

# ------------------------------------------------------------------------------
# 2-D PCA biplot: GrowthForm as grouping factor
# ------------------------------------------------------------------------------
# 
# TODO: Added 27/12-24
# SET GLOBAL LIMIT FOR contribution levels ('contrib') FOR LEGEND MERGING
contrib_limits = c(5, 30)
#
# 
plt_nonphyloPCA_biplot_ax12 = fviz_pca_biplot(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta_rmna$GrowthForm,
  axes=c(1,2),
  repel=T,
  col.var = "gray20",
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#EECA40", "#7998AD", "#F07673")) + # liana, shrub, tree
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank())

plt_nonphyloPCA_contribplot_ax12 = fviz_pca_var(
  nonphyloPCAresult, col.var = "contrib",
  axes=c(1,2),
  repel=T
  #gradient.cols = c("#7998AD", "#EECA40", "#F07673")
  ) +
  theme_classic() +
  theme(legend.direction = 'vertical', legend.position = 'none') + labs(title = "") +
  scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08"), limits = contrib_limits)


# Save plot after merging ------------------------------------------------------
plt_ax12 = plt_nonphyloPCA_biplot_ax12 + plt_nonphyloPCA_contribplot_ax12 + 
  #plot_annotation(title = "PC1-2, biplot and contribution plot") +
  plot_layout(guides = "collect") & theme(legend.position='right')

# TODO: 25-02-01, extract contrib and cos2 in numbers
#plt_nonphyloPCA_contribplot_ax12[["data"]] %>% view()
#plt_nonphyloPCA_contribplot_ax23[["data"]] %>% view()

# Combine the two plots side by side and save
ggsave(plot = plt_ax12, filename = "2502-indv-level-code/Ind_nonphyPCA_ax12_RES.pdf",
       width = 6.5, height = 3.5)

