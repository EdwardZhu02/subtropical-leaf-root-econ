rm(list=ls())
library(tidyverse)
library(dplyr)
library(factoextra) # perform PCA analysis
library(FactoMineR) # perform PCA analysis
library(patchwork) # plot merging
library(cowplot) # plot merging

load("outdata/traitDataFujian-Ind-step1.RData")
if (!dir.exists("outplts/SIFigures")) {
  dir.create("outplts/SIFigures", showWarnings = FALSE, recursive = TRUE)
}

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

# Fig S4A-B: Site H, contrib plots ----------------------------------------------
traitDataPCA_touse = traitDataIndv_SelectedTraits_log %>% dplyr::filter(SiteID %in% c("hilltop"))
traitName_touse = c("RTD","SRL","RD","RNC","LMA","LNC")

# Resolve problem: non-unique values when setting 'row.names': ‘Castanopsis fordii’, ‘Machilus pauhoi’ 
# because of sp-gf averaging
traitDataPCA_touse = traitDataPCA_touse %>% mutate(
  PCAIdentifier_spgf = paste(gsub(" ","_",SpeciesFullName), SiteID, SampleID, sep="-")
)

# Curate data for PCA analysis
nonphyloPCAData_numonly = traitDataPCA_touse %>% ungroup() %>% 
  dplyr::select(PCAIdentifier_spgf,SiteID,all_of(traitName_touse)) %>% 
  dplyr::mutate(SiteID = as.factor(SiteID)) %>%
  na.omit()

# Added 27/12-24, concordant with angle calculation
# scale traits based on columns (mean=0, SD=1, same as z-transform) 
nonphyloPCAData_numonly[, traitName_touse] = apply(nonphyloPCAData_numonly[, traitName_touse], 2, scale)
# fix rownames
rownames(nonphyloPCAData_numonly) = nonphyloPCAData_numonly$PCAIdentifier_spgf

# used for growth form annotation in PCA biplot
nonphyloPCAData_meta_rmna = nonphyloPCAData_numonly
nonphyloPCAData_numonly = nonphyloPCAData_numonly %>% dplyr::select(-PCAIdentifier_spgf, -SiteID)

nonphyloPCAresult = PCA(nonphyloPCAData_numonly, scale.unit = TRUE, graph = FALSE)

# Apply varimax rotation on the variable loadings from the PCA result
varimax_result <- varimax(nonphyloPCAresult$var$coord)
rotated_loadings <- as.matrix(varimax_result$loadings)

# Compute rotated individual scores by applying the rotation matrix
rotated_scores <- nonphyloPCAresult$ind$coord %*% varimax_result$rotmat

# Create a new PCA object by updating the original PCA result with rotated values.
nonphyloPCAresult_varimax_rotated <- nonphyloPCAresult
nonphyloPCAresult_varimax_rotated$var$coord <- rotated_loadings
nonphyloPCAresult_varimax_rotated$ind$coord <- rotated_scores  

# configurable: use raw or varimax rotated data
# nonphyloPCAresult / nonphyloPCAresult_varimax_rotated
nonphyloPCAresult_touse = nonphyloPCAresult_varimax_rotated

# Loadings and eigenvalues
loadings_mat = as.matrix(nonphyloPCAresult_touse$var$coord)
loadings_full = loadings_mat[, , drop = FALSE]
print(round(loadings_full, 4))

eig = nonphyloPCAresult_touse$eig
print(round(eig, 4))

# 2-D PCA biplot: SiteID as grouping factor
# (24.12.27) SET GLOBAL LIMIT FOR contribution levels ('contrib') FOR LEGEND MERGING
#contrib_limits = c(5, 30)

plt_nonphyloPCA_contrib_ax12 = fviz_pca_var(
  nonphyloPCAresult, col.var = "contrib", axes=c(1,2), repel=T) +
  theme_classic() +
  theme(legend.direction = 'vertical', legend.position = "right") + labs(title = "") +
  scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08"))#, limits = contrib_limits)

plt_nonphyloPCA_biplot_ax12 = fviz_pca_biplot_MODIFIED(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta_rmna$SiteID,
  axes=c(1,2),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  #select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#547bb4", "#dd7c4f")) + # hilltop,valley
  #annotate("text", x = -2.5, y = 3.8, label=paste0("Top ", length(traitName_touse)-2), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank()) + xlim(-5, 5) + ylim(-5, 5)

plt_nonphyloPCA_contrib_ax23 = fviz_pca_var(
  nonphyloPCAresult, col.var = "contrib", axes=c(2,3), repel=T) +
  theme_classic() +
  theme(legend.direction = 'vertical', legend.position = "right") + labs(title = "") +
  scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08"))#, limits = contrib_limits)

plt_nonphyloPCA_biplot_ax23 = fviz_pca_biplot_MODIFIED(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta_rmna$SiteID, 
  axes=c(2,3),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  #select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#547bb4", "#dd7c4f")) + # hilltop, valley
  #annotate("text", x = -2.5, y = 3.8, label=paste0("Top ", length(traitName_touse)-2), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank()) + xlim(-5, 5) + ylim(-5, 5)

# plt_ax123_bi = plt_nonphyloPCA_biplot_ax12 + plt_nonphyloPCA_biplot_ax23 +
#   plot_layout(guides = "collect") & theme(legend.position='bottom')
# ggsave(plot = plt_ax123_bi_varimax, filename = "outplts/SIFigures/Fig4AB-Ind_PCAbi_siteH.pdf",
#        width = 7.5, height = 4.5)

plt_ax123_contrib = plt_nonphyloPCA_contrib_ax12 + plt_nonphyloPCA_contrib_ax23
ggsave(plot = plt_ax123_contrib, filename = "outplts/SIFigures/FigS4AB-Ind_PCAcontrib_siteH.pdf",
       width = 9, height = 5)

# Fig S4C-D: Site V, contrib plots ---------------------------------------------
traitDataPCA_touse = traitDataIndv_SelectedTraits_log %>% dplyr::filter(SiteID %in% c("valley"))
traitName_touse = c("RTD","SRL","RD","RNC","LMA","LNC")

# Resolve problem: non-unique values when setting 'row.names': ‘Castanopsis fordii’, ‘Machilus pauhoi’ 
# because of sp-gf averaging
traitDataPCA_touse = traitDataPCA_touse %>% mutate(
  PCAIdentifier_spgf = paste(gsub(" ","_",SpeciesFullName), SiteID, SampleID, sep="-")
)

# Curate data for PCA analysis
nonphyloPCAData_numonly = traitDataPCA_touse %>% ungroup() %>% 
  dplyr::select(PCAIdentifier_spgf,SiteID,all_of(traitName_touse)) %>% 
  dplyr::mutate(SiteID = as.factor(SiteID)) %>%
  na.omit()

# Added 27/12-24, concordant with angle calculation
# scale traits based on columns (mean=0, SD=1, same as z-transform) 
nonphyloPCAData_numonly[, traitName_touse] = apply(nonphyloPCAData_numonly[, traitName_touse], 2, scale)
# fix rownames
rownames(nonphyloPCAData_numonly) = nonphyloPCAData_numonly$PCAIdentifier_spgf

# used for growth form annotation in PCA biplot
nonphyloPCAData_meta_rmna = nonphyloPCAData_numonly
nonphyloPCAData_numonly = nonphyloPCAData_numonly %>% dplyr::select(-PCAIdentifier_spgf, -SiteID)

nonphyloPCAresult = PCA(nonphyloPCAData_numonly, scale.unit = TRUE, graph = FALSE)

# Apply varimax rotation on the variable loadings from the PCA result
varimax_result <- varimax(nonphyloPCAresult$var$coord)
rotated_loadings <- as.matrix(varimax_result$loadings)

# Compute rotated individual scores by applying the rotation matrix
rotated_scores <- nonphyloPCAresult$ind$coord %*% varimax_result$rotmat

# Create a new PCA object by updating the original PCA result with rotated values.
nonphyloPCAresult_varimax_rotated <- nonphyloPCAresult
nonphyloPCAresult_varimax_rotated$var$coord <- rotated_loadings
nonphyloPCAresult_varimax_rotated$ind$coord <- rotated_scores  

# configurable: use raw or varimax rotated data
# nonphyloPCAresult / nonphyloPCAresult_varimax_rotated
nonphyloPCAresult_touse = nonphyloPCAresult_varimax_rotated

# Loadings and eigenvalues
loadings_mat = as.matrix(nonphyloPCAresult_touse$var$coord)
loadings_full = loadings_mat[, , drop = FALSE]
print(round(loadings_full, 4))

eig = nonphyloPCAresult_touse$eig
print(round(eig, 4))

# 2-D PCA biplot: SiteID as grouping factor
# (24.12.27) SET GLOBAL LIMIT FOR contribution levels ('contrib') FOR LEGEND MERGING
#contrib_limits = c(5, 30)

plt_nonphyloPCA_contrib_ax12 = fviz_pca_var(
  nonphyloPCAresult, col.var = "contrib", axes=c(1,2), repel=T) +
  theme_classic() +
  theme(legend.direction = 'vertical', legend.position = "right") + labs(title = "") +
  scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08"))#, limits = contrib_limits)

plt_nonphyloPCA_biplot_ax12 = fviz_pca_biplot_MODIFIED(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta_rmna$SiteID,
  axes=c(1,2),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  #select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#547bb4", "#dd7c4f")) + # hilltop,valley
  #annotate("text", x = -2.5, y = 3.8, label=paste0("Top ", length(traitName_touse)-2), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank()) + xlim(-5, 5) + ylim(-5, 5)

plt_nonphyloPCA_contrib_ax23 = fviz_pca_var(
  nonphyloPCAresult, col.var = "contrib", axes=c(2,3), repel=T) +
  theme_classic() +
  theme(legend.direction = 'vertical', legend.position = "right") + labs(title = "") +
  scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08"))#, limits = contrib_limits)

plt_nonphyloPCA_biplot_ax23 = fviz_pca_biplot_MODIFIED(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta_rmna$SiteID, 
  axes=c(2,3),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  #select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#547bb4", "#dd7c4f")) + # hilltop, valley
  #annotate("text", x = -2.5, y = 3.8, label=paste0("Top ", length(traitName_touse)-2), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank()) + xlim(-5, 5) + ylim(-5, 5)

# plt_ax123_bi = plt_nonphyloPCA_biplot_ax12 + plt_nonphyloPCA_biplot_ax23 +
#   plot_layout(guides = "collect") & theme(legend.position='bottom')
# ggsave(plot = plt_ax123_bi_varimax, filename = "outplts/SIFigures/Fig4AB-Ind_PCAbi_siteH.pdf",
#        width = 7.5, height = 4.5)

plt_ax123_contrib = plt_nonphyloPCA_contrib_ax12 + plt_nonphyloPCA_contrib_ax23
ggsave(plot = plt_ax123_contrib, filename = "outplts/SIFigures/FigS4CD-Ind_PCAcontrib_siteV.pdf",
       width = 9, height = 5)

# Fig S4E-F: Combined biplot ---------------------------------------------------
traitDataPCA_touse = traitDataIndv_SelectedTraits_log
traitName_touse = c("RTD","SRL","RD","RNC","LMA","LNC")

# Resolve problem: non-unique values when setting 'row.names': ‘Castanopsis fordii’, ‘Machilus pauhoi’ 
# because of sp-gf averaging
traitDataPCA_touse = traitDataPCA_touse %>% mutate(
  PCAIdentifier_spgf = paste(gsub(" ","_",SpeciesFullName), SiteID, SampleID, sep="-")
)

# Curate data for PCA analysis
nonphyloPCAData_numonly = traitDataPCA_touse %>% ungroup() %>% 
  dplyr::select(PCAIdentifier_spgf,SiteID,all_of(traitName_touse)) %>% 
  dplyr::mutate(SiteID = as.factor(SiteID)) %>%
  na.omit()

# Added 27/12-24, concordant with angle calculation
# scale traits based on columns (mean=0, SD=1, same as z-transform) 
nonphyloPCAData_numonly[, traitName_touse] = apply(nonphyloPCAData_numonly[, traitName_touse], 2, scale)
# fix rownames
rownames(nonphyloPCAData_numonly) = nonphyloPCAData_numonly$PCAIdentifier_spgf

# used for growth form annotation in PCA biplot
nonphyloPCAData_meta_rmna = nonphyloPCAData_numonly
nonphyloPCAData_numonly = nonphyloPCAData_numonly %>% dplyr::select(-PCAIdentifier_spgf, -SiteID)

nonphyloPCAresult = PCA(nonphyloPCAData_numonly, scale.unit = TRUE, graph = FALSE)

# Apply varimax rotation on the variable loadings from the PCA result
varimax_result <- varimax(nonphyloPCAresult$var$coord)
rotated_loadings <- as.matrix(varimax_result$loadings)

# Compute rotated individual scores by applying the rotation matrix
rotated_scores <- nonphyloPCAresult$ind$coord %*% varimax_result$rotmat

# Create a new PCA object by updating the original PCA result with rotated values.
nonphyloPCAresult_varimax_rotated <- nonphyloPCAresult
nonphyloPCAresult_varimax_rotated$var$coord <- rotated_loadings
nonphyloPCAresult_varimax_rotated$ind$coord <- rotated_scores  

# configurable: use raw or varimax rotated data
# nonphyloPCAresult / nonphyloPCAresult_varimax_rotated
nonphyloPCAresult_touse = nonphyloPCAresult_varimax_rotated

# Loadings and eigenvalues
loadings_mat = as.matrix(nonphyloPCAresult_touse$var$coord)
loadings_full = loadings_mat[, , drop = FALSE]
print(round(loadings_full, 4))

eig = nonphyloPCAresult_touse$eig
print(round(eig, 4))

# 2-D PCA biplot: SiteID as grouping factor
# (24.12.27) SET GLOBAL LIMIT FOR contribution levels ('contrib') FOR LEGEND MERGING
#contrib_limits = c(5, 30)

plt_nonphyloPCA_contrib_ax12 = fviz_pca_var(
  nonphyloPCAresult, col.var = "contrib", axes=c(1,2), repel=T) +
  theme_classic() +
  theme(legend.direction = 'vertical', legend.position = "right") + labs(title = "") +
  scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08"))#, limits = contrib_limits)

plt_nonphyloPCA_biplot_ax12 = fviz_pca_biplot_MODIFIED(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta_rmna$SiteID,
  axes=c(1,2),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  #select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#547bb4", "#dd7c4f")) + # hilltop,valley
  #annotate("text", x = -2.5, y = 3.8, label=paste0("Top ", length(traitName_touse)-2), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank()) + xlim(-5, 5) + ylim(-5, 5)

plt_nonphyloPCA_contrib_ax23 = fviz_pca_var(
  nonphyloPCAresult, col.var = "contrib", axes=c(2,3), repel=T) +
  theme_classic() +
  theme(legend.direction = 'vertical', legend.position = "right") + labs(title = "") +
  scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08"))#, limits = contrib_limits)

plt_nonphyloPCA_biplot_ax23 = fviz_pca_biplot_MODIFIED(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta_rmna$SiteID, 
  axes=c(2,3),
  repel=T,
  col.var = "gray20",
  #select.var = list(cos2=0.3),
  #select.var = list(contrib=length(traitName_touse)-2),
  # cos2: if cos2 is in [0, 1], ex: 0.6, then individuals/variables with a cos2 > 0.6 are drawn.
  # contrib: if contrib > 1, ex: 5, then the top 5 individuals/variables with the highest contrib are drawn
  # addEllipses=T, palette = c("#8583A9", "#A0BDD5", "#CA8BA8")) + # liana, shrub, tree
  addEllipses=T, palette = c("#547bb4", "#dd7c4f")) + # hilltop, valley
  #annotate("text", x = -2.5, y = 3.8, label=paste0("Top ", length(traitName_touse)-2), color = "black") + 
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'none', legend.title = element_blank()) + xlim(-5, 5) + ylim(-5, 5)

plt_ax123_bi = plt_nonphyloPCA_biplot_ax12 + plt_nonphyloPCA_biplot_ax23 +
  plot_layout(guides = "collect") & theme(legend.position='bottom')
ggsave(plot = plt_ax123_bi, filename = "outplts/SIFigures/FigS4EF-Ind_PCAbi_allsites.pdf",
       width = 7.5, height = 4.5)
