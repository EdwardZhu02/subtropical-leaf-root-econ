#+#+#+#+############################################################
# Species-mean PCA with individual projections
# Addresses nesting/replication bias by fitting PCA on species means
# and projecting individuals into that space.
############################################################

rm(list=ls())
library(tidyverse)
library(dplyr)
library(factoextra)
library(FactoMineR)
library(patchwork)
library(cowplot)

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
traitDataPCA_touse = traitDataIndv_SelectedTraits_log

### SET TRAITS TO USE ###
## Fine root trait space (Figure 2)
#traitName_touse = c("RTD","SRL","RD","RNC","RDMC","SRR25","SRA","RPC","RCC")

## Leaf trait space (Figure 2)
traitName_touse = c("LMA","LNC","LPC","LCC","Ld13C","Rdark25P","Vcmax25","Asat","LA")

# Markers and identified co-predictors, leaf+root (Figure 2)
#traitName_touse = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA","LPC")

# Create unique ID for individuals
traitDataPCA_touse = traitDataPCA_touse %>% mutate(
	PCAIdentifier_spgf = paste(gsub(" ","_",SpeciesFullName), SiteID, SampleID, sep="-")
)

# Individual-level data for projection
indv_data <- traitDataPCA_touse %>% ungroup() %>%
	dplyr::select(PCAIdentifier_spgf, SpeciesFullName, SiteID, all_of(traitName_touse)) %>%
	dplyr::mutate(SiteID = as.factor(SiteID)) %>%
	tidyr::drop_na(all_of(traitName_touse))

# Species means (equal weight per species)
spmean_data <- indv_data %>%
	dplyr::group_by(SpeciesFullName) %>%
	dplyr::summarise(across(all_of(traitName_touse), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
	tidyr::drop_na(all_of(traitName_touse))

spmean_matrix <- spmean_data %>% dplyr::select(all_of(traitName_touse))
indv_matrix <- indv_data %>% dplyr::select(all_of(traitName_touse))

rownames(spmean_matrix) <- spmean_data$SpeciesFullName
rownames(indv_matrix) <- indv_data$PCAIdentifier_spgf

#-----------------------------------------------------------
# PCA on species means
pca_spmean <- PCA(spmean_matrix, scale.unit = TRUE, graph = FALSE)

# Project individuals into species-mean PCA space
pred_indv <- predict.PCA(pca_spmean, newdata = indv_matrix)

# Build a PCA object with projected individuals for plotting
pca_indv_proj <- pca_spmean
pca_indv_proj$ind$coord <- pred_indv$coord
pca_indv_proj$ind$cos2 <- pred_indv$cos2
if (is.null(pred_indv$contrib)) {
  pca_indv_proj$ind$contrib <- matrix(
    0,
    nrow = nrow(pred_indv$coord),
    ncol = ncol(pred_indv$coord),
    dimnames = dimnames(pred_indv$coord)
  )
} else {
  pca_indv_proj$ind$contrib <- pred_indv$contrib
}

#-----------------------------------------------------------
# Varimax rotation (optional, for comparison with original workflow)
varimax_result <- varimax(pca_spmean$var$coord)
rotated_loadings <- as.matrix(varimax_result$loadings)

pca_spmean_varimax <- pca_spmean
pca_spmean_varimax$var$coord <- rotated_loadings
pca_spmean_varimax$ind$coord <- pca_spmean$ind$coord %*% varimax_result$rotmat

pca_indv_proj_varimax <- pca_indv_proj
pca_indv_proj_varimax$var$coord <- rotated_loadings
pca_indv_proj_varimax$ind$coord <- pred_indv$coord %*% varimax_result$rotmat
if (is.null(pca_indv_proj_varimax$ind$contrib)) {
  pca_indv_proj_varimax$ind$contrib <- matrix(
    0,
    nrow = nrow(pca_indv_proj_varimax$ind$coord),
    ncol = ncol(pca_indv_proj_varimax$ind$coord),
    dimnames = dimnames(pca_indv_proj_varimax$ind$coord)
  )
}

#-----------------------------------------------------------
# Output tables
write.csv(pca_spmean$var$coord, "revision-code/step5-spmeanPCA-var-coord.csv", row.names = TRUE)
write.csv(pca_spmean$ind$coord, "revision-code/step5-spmeanPCA-species-scores.csv", row.names = TRUE)
write.csv(pred_indv$coord, "revision-code/step5-spmeanPCA-individual-proj-scores.csv", row.names = TRUE)

write.csv(pca_spmean_varimax$var$coord, "revision-code/step5-spmeanPCA-varimax-var-coord.csv", row.names = TRUE)
write.csv(pca_spmean_varimax$ind$coord, "revision-code/step5-spmeanPCA-varimax-species-scores.csv", row.names = TRUE)
write.csv(pca_indv_proj_varimax$ind$coord, "revision-code/step5-spmeanPCA-varimax-individual-proj-scores.csv", row.names = TRUE)

#-----------------------------------------------------------
# Biplots: projected individuals colored by SiteID
plt_indv_proj_ax12 <- fviz_pca_biplot_MODIFIED(
	pca_indv_proj, label = "var", habillage = indv_data$SiteID,
	axes = c(1,2),
	repel = TRUE,
	col.var = "gray20",
	addEllipses = TRUE, palette = c("#547bb4", "#dd7c4f")
) +
	labs(title = "Species-mean PCA + individual projection") +
	theme_classic() +
	theme(legend.direction = 'horizontal', legend.position = 'bottom', legend.title = element_blank()) +
	xlim(-5, 5) + ylim(-5, 5)

plt_indv_proj_ax23 <- fviz_pca_biplot_MODIFIED(
	pca_indv_proj, label = "var", habillage = indv_data$SiteID,
	axes = c(2,3),
	repel = TRUE,
	col.var = "gray20",
	addEllipses = TRUE, palette = c("#547bb4", "#dd7c4f")
) +
	labs(title = "") +
	theme_classic() +
	theme(legend.direction = 'horizontal', legend.position = 'bottom', legend.title = element_blank()) +
	xlim(-5, 5) + ylim(-5, 5)

plt_indv_proj_ax12_varimax <- fviz_pca_biplot_MODIFIED(
	pca_indv_proj_varimax, label = "var", habillage = indv_data$SiteID,
	axes = c(1,2),
	repel = TRUE,
	col.var = "gray20",
	addEllipses = TRUE, palette = c("#547bb4", "#dd7c4f")
) +
	labs(title = "Species-mean PCA + individual projection (varimax)") +
	theme_classic() +
	theme(legend.direction = 'horizontal', legend.position = 'bottom', legend.title = element_blank()) +
	xlim(-5, 5) + ylim(-5, 5)

plt_indv_proj_ax23_varimax <- fviz_pca_biplot_MODIFIED(
	pca_indv_proj_varimax, label = "var", habillage = indv_data$SiteID,
	axes = c(2,3),
	repel = TRUE,
	col.var = "gray20",
	addEllipses = TRUE, palette = c("#547bb4", "#dd7c4f")
) +
	labs(title = "") +
	theme_classic() +
	theme(legend.direction = 'horizontal', legend.position = 'bottom', legend.title = element_blank()) +
	xlim(-5, 5) + ylim(-5, 5)

plt_spmean_proj_biplots <- plt_indv_proj_ax12 + plt_indv_proj_ax23 +
	plot_layout(guides = "collect") & theme(legend.position='bottom')

plt_spmean_proj_biplots_varimax <- plt_indv_proj_ax12_varimax + plt_indv_proj_ax23_varimax +
	plot_layout(guides = "collect") & theme(legend.position='bottom')

ggsave(plot = plt_spmean_proj_biplots, filename = "revision-code/step5-spmean-indvproj-bi-LES.pdf",
			 width = 5.8, height = 3.6)
ggsave(plot = plt_spmean_proj_biplots_varimax, filename = "revision-code/step5-spmean-indvproj-bi-varimax-LES.pdf",
			 width = 5.8, height = 3.6)
