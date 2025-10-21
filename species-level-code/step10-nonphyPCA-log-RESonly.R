rm(list=ls())
library(tidyverse)
library(dplyr)
library(factoextra) # perform PCA analysis
library(FactoMineR) # perform PCA analysis
# library(rgl) # 3-D PCA plot
library(patchwork) # plot merging
library(cowplot) # plot merging

# ------------------------------------------------------------------------------
# Perform PCA analysis
# modified in response to Sandy Harrison's comment -> Use original values without normalization
# to perform the PCA to show if the conclusion was altered
# ------------------------------------------------------------------------------
load("traitDataFujian-spavg-phylo-step2.RData")

### SET DATA TO USE ###
traitDataPCA_touse = traitDataIndv_spgfavg_log # perform scaling later (line 36, 27/12-24)
### SET TRAITS TO USE###
# All traits, need to be log transformed and scaled 
traitName_touse = c("RTD","SRL","RD","RNC")

# Resolve problem: non-unique values when setting 'row.names': ‘Castanopsis fordii’, ‘Machilus pauhoi’ 
# because of sp-gf averaging
traitDataPCA_touse = traitDataPCA_touse %>% mutate(
  PCAIdentifier_spgf = paste(gsub(" ","_",speciesFullName), GrowthForm, sep="-")
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

# perform PCA
nonphyloPCAresult = PCA(nonphyloPCAData_numonly, scale.unit = T, graph = F) 

#phyloPCA.results = phyl.pca(traitPhyloDataSH_tree, phyloPCAData_numonly, method = 'lambda')
#phyloPCA.results.prcomp = as.prcomp(phyloPCA.results)
#
# ------------------------------------------------------------------------------
# Scree plot: contributions of axes
# ------------------------------------------------------------------------------
nonphyloPCA.screeplot = fviz_screeplot(nonphyloPCAresult, addlabels = TRUE, ylim = c(0, 60),
                                       barfill="#7998AD", barcolor="#7998AD", linecolor="black") +
  labs(title="RES") + theme_classic()
ggsave(filename = "plt_nonphyPCA_scree-RES.pdf", plot = nonphyloPCA.screeplot, width = 3, height = 3)

# ------------------------------------------------------------------------------
# 2-D PCA biplot: GrowthForm as grouping factor
# ------------------------------------------------------------------------------
# 
# TODO: Added 27/12-24
# SET GLOBAL LIMIT FOR contribution levels ('contrib') FOR LEGEND MERGING
# contrib_limits = c(5, 30)
#
# 
plt_nonphyloPCA_biplot_ax12 = fviz_pca_biplot(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta_rmna$GrowthForm,
  axes=c(1,2),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  # select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
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
  scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08"))#, limits = contrib_limits)


# Save plot after merging ------------------------------------------------------
plt_ax12 = plt_nonphyloPCA_biplot_ax12 + plt_nonphyloPCA_contribplot_ax12 + 
  #plot_annotation(title = "PC1-2, biplot and contribution plot") +
  plot_layout(guides = "collect") & theme(legend.position='right')

ggsave(plot = plt_ax12, filename = "plt_nonphyPCA_ax12Combined_RES.pdf",
       width = 6.5, height = 3.5)


