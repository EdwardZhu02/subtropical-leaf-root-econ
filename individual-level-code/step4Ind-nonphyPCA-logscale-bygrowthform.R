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

load("individual-level-code/traitDataFujian-Ind-step1.RData")

.is_grouping_var <- function(x) {length(x) > 1 & (is.character(x) | is.factor(x))}
.is_continuous_var <- function(x) {x[1] %in% c("cos2", "contrib", "x", "y") | is.numeric(x)}

fviz_pca_biplot_MODIFIED <- function (X, axes = c(1, 2), geom = c("point", "text"), geom.ind = geom, 
          geom.var = c("arrow", "text"), col.ind = "black", fill.ind = "white", 
          col.var = "steelblue", fill.var = "white", gradient.cols = NULL, 
          label = "all", invisible = "none", repel = FALSE, habillage = "none", 
          palette = NULL, addEllipses = FALSE, title = "", 
          ...) 
{
  is.individuals.colored.by.variable <- .is_grouping_var(fill.ind) | 
    .is_grouping_var(col.ind)
  is.variables.colored.by.variable <- .is_continuous_var(col.var) | 
    .is_grouping_var(col.var)
  is.gradient.color <- .is_continuous_var(col.ind) | .is_continuous_var(col.var)
  is.gradient.fill <- .is_continuous_var(fill.ind) | .is_continuous_var(fill.var)
  is.discrete.color <- .is_grouping_var(col.ind) | .is_grouping_var(habillage) | 
    .is_grouping_var(col.var)
  is.discrete.fill <- .is_grouping_var(fill.ind) | .is_grouping_var(fill.var) | 
    .is_grouping_var(habillage) | (.is_grouping_var(col.ind) & 
                                     addEllipses)
  var <- facto_summarize(X, element = "var", result = c("coord", 
                                                        "contrib", "cos2"), axes = axes)
  colnames(var)[2:3] <- c("x", "y")
  pca.ind <- get_pca_ind(X)
  ind <- data.frame(pca.ind$coord[, axes, drop = FALSE], stringsAsFactors = TRUE)
  colnames(ind) <- c("x", "y")
  r <- min((max(ind[, "x"]) - min(ind[, "x"])/(max(var[, "x"]) - 
                                                 min(var[, "x"]))), (max(ind[, "y"]) - min(ind[, "y"])/(max(var[, 
                                                                                                                "y"]) - min(var[, "y"]))))
  ellipse.border.remove <- T
  # if (is.individuals.colored.by.variable & is.variables.colored.by.variable) 
  #   ellipse.border.remove <- TRUE
  p <- fviz_pca_ind(X, axes = axes, geom = geom.ind, repel = repel, 
                    col.ind = col.ind, fill.ind = fill.ind, label = label, 
                    invisible = invisible, habillage = habillage, addEllipses = addEllipses, 
                    ellipse.border.remove = ellipse.border.remove, ...)
  p <- fviz_pca_var(X, axes = axes, geom = geom.var, repel = repel, 
                    col.var = col.var, fill.var = fill.var, label = label, 
                    invisible = invisible, scale. = r * 0.7, ggp = p, ...)
  if (!is.null(gradient.cols)) {
    if (is.gradient.color) 
      p <- p + ggpubr::gradient_color(gradient.cols)
    if (is.gradient.fill) 
      p <- p + ggpubr::gradient_fill(gradient.cols)
  }
  if (!is.null(palette)) {
    if (is.discrete.color) 
      p <- p + ggpubr::color_palette(palette)
    if (is.discrete.fill) 
      p <- p + ggpubr::fill_palette(palette)
  }
  p
}


### SET DATA TO USE ### 
# TODO: perform scaling later (line 36, 27/12-24)
traitDataPCA_touse = traitDataIndv_SelectedTraits_log # perform scaling later (line 36, 27/12-24)

### SET TRAITS TO USE###

## Fine root trait space (Figure 2)
#traitName_touse = c("RTD","SRL","RD","RNC","RDMC","SRR25","SRA","RPC","RCC")

## Leaf trait space (Figure 2)
#traitName_touse = c("LMA","LNC","LPC","LCC","Ld13C","Rdark25P","Vcmax25","Asat","LA")

# Markers and identified co-predictors, leaf+root (Figure 2)
traitName_touse = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA","LPC")


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
nonphyloPCA.screeplot = fviz_screeplot(nonphyloPCAresult, addlabels = TRUE, ylim = c(0, 40),
                                       barfill="#7998AD", barcolor="#7998AD", linecolor="black") +
  labs(title="R-LES") + theme_classic()
ggsave(filename = "2502-indv-level-code/Ind_nonphyPCA_scree-LESall.pdf", plot = nonphyloPCA.screeplot, width = 3, height = 3)

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
  #select.var = list(cos2=0.3),
  select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#EECA40", "#7998AD", "#F07673")) + # liana, shrub, tree
  annotate("text", x = -2.5, y = 3.8, label=paste0("Top ", length(traitName_touse)-2), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank())

# plt_nonphyloPCA_contribplot_ax12 = fviz_pca_var(
#   nonphyloPCAresult, col.var = "contrib",
#   axes=c(1,2),
#   repel=T
#   #gradient.cols = c("#7998AD", "#EECA40", "#F07673")
# ) +
#   theme_classic() +
#   #annotate("text", x = -0.9, y = 1, label="All", color = "black") + 
#   theme(legend.direction = 'vertical', legend.position = 'none') + labs(title = "") +
#   scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08"))#, limits = contrib_limits)


plt_nonphyloPCA_biplot_ax23 = fviz_pca_biplot(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta_rmna$GrowthForm, 
  axes=c(2,3),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  #select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#EECA40", "#7998AD", "#F07673")) + # liana, shrub, tree
  annotate("text", x = -2.5, y = 3.8, label=paste0("Top ", length(traitName_touse)-2), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank())

# plt_nonphyloPCA_contribplot_ax23 = fviz_pca_var(
#   nonphyloPCAresult, col.var = "contrib",
#   axes=c(2,3),
#   repel=T
#   #gradient.cols = c("#7998AD", "#EECA40", "#F07673")
# ) +
#   theme_classic() +
#   #annotate("text", x = -0.9, y = 1, label="All", color = "black") + 
#   theme(legend.direction = 'vertical', legend.position = "none") + labs(title = "") +
#   scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08"))#, limits = contrib_limits)


plt_nonphyloPCA_biplot_ax12.varimax = fviz_pca_biplot_MODIFIED(
  nonphyloPCAresult_varimax_rotated, label = "var", habillage=nonphyloPCAData_meta_rmna$GrowthForm,
  axes=c(1,2),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#EECA40", "#7998AD", "#F07673")) + # liana, shrub, tree
  #annotate("text", x = -4, y = 4.9, label=paste0("PC1-PC2"), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank()) + xlim(-5, 5) + ylim(-5, 5)

# plt_nonphyloPCA_contribplot_ax12.varimax = fviz_pca_var(
#   nonphyloPCAresult_varimax_rotated, col.var = "contrib",
#   axes=c(1,2),
#   repel=T
#   #gradient.cols = c("#7998AD", "#EECA40", "#F07673")
# ) +
#   theme_classic() +
#   #annotate("text", x = -0.9, y = 1, label="All", color = "black") + 
#   theme(legend.direction = 'vertical', legend.position = 'none') + labs(title = "") +
#   scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08"))#, limits = contrib_limits)


plt_nonphyloPCA_biplot_ax23.varimax = fviz_pca_biplot_MODIFIED(
  nonphyloPCAresult_varimax_rotated, label = "var", habillage=nonphyloPCAData_meta_rmna$GrowthForm, 
  axes=c(2,3),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#EECA40", "#7998AD", "#F07673")) + # liana, shrub, tree
  #annotate("text", x = -4, y = 4.9, label=paste0("PC2-PC3"), color = "black") +  
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank()) + xlim(-5, 5) + ylim(-5, 5)

# plt_nonphyloPCA_contribplot_ax23.varimax = fviz_pca_var(
#   nonphyloPCAresult_varimax_rotated, col.var = "contrib",
#   axes=c(2,3),
#   repel=T
#   #gradient.cols = c("#7998AD", "#EECA40", "#F07673")
# ) +
#   theme_classic() +
#   annotate("text", x = -0.9, y = 1, label="All", color = "black") + 
#   theme(legend.direction = 'vertical', legend.position = "none") + labs(title = "") +
#   scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08"))#, limits = contrib_limits)


# TODO: 25-02-01, extract contrib and cos2 in numbers
#plt_nonphyloPCA_contribplot_ax12[["data"]] %>% view()
#plt_nonphyloPCA_contribplot_ax23[["data"]] %>% view()
nonphyloPCAresult_varimax_rotated[["var"]][["coord"]] %>% view()


# Make the output plot
plt_ax123_bi_varimax = plt_nonphyloPCA_biplot_ax12.varimax + plt_nonphyloPCA_biplot_ax23.varimax +
  plot_layout(guides = "collect") & theme(legend.position='bottom')
# ggsave(plot = plt_ax123_bi_varimax, filename = "2502-indv-level-code/Ind_new_BiPlot_RLES_varimax.pdf",
#        width = 5.8, height = 3.6)
ggsave(plot = plt_ax123_bi_varimax, filename = "2502-indv-level-code/Ind_new_BiPlot_RLES_varimax.pdf",
       width = 6.4, height = 4)

# Combine the two plots side by side and save
# plt_ax12_ax23_combined = plot_grid(plt_ax12, plt_ax23, ncol = 1)
# ggsave(plot = plt_ax12_ax23_combined, filename = "2502-indv-level-code/Ind_nonphyPCA_ax123_LESall.pdf",
#        width = 7.5, height = 8.1)
# 
# plt_ax12_ax23_combined.varimax = plot_grid(plt_ax12.varimax, plt_ax23.varimax, ncol = 1)
# ggsave(plot = plt_ax12_ax23_combined.varimax, filename = "2502-indv-level-code/Ind_nonphyPCA_ax123_LESall_varimax.pdf",
#        width = 7.5, height = 8.1)

