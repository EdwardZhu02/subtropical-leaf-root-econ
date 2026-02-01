# ****************************************************************
# One-individual-per-species resampling for individual-level PCA
# Repeatedly sample one individual per species within each community
# and recompute PCA (with optional varimax rotation)
# ****************************************************************

rm(list=ls())
library(tidyverse)
library(dplyr)
library(factoextra) # perform PCA analysis
library(FactoMineR) # perform PCA analysis
library(patchwork) # plot merging

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
## Fine root trait space (Figure 2)
#traitName_touse = c("RTD","SRL","RD","RNC","RDMC","SRR25","SRA","RPC","RCC")
## Leaf trait space (Figure 2)
traitName_touse = c("LMA","LNC","LPC","LCC","Ld13C","Rdark25P","Vcmax25","Asat","LA")
# Markers and identified co-predictors, leaf+root (Figure 2)
#traitName_touse = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA","LPC")

### RESAMPLING SETTINGS ###
n_iter = 200
set.seed(1)

sample_one_ind_per_sp_site <- function(data, trait_cols) {
	data %>%
		dplyr::filter(if_all(all_of(trait_cols), ~ !is.na(.))) %>%
		group_by(SiteID, SpeciesFullName) %>%
		slice_sample(n = 1) %>%
		ungroup()
}

resampling_varcoords = vector("list", n_iter)
resampling_varcoords_varimax = vector("list", n_iter)
resampling_eig = vector("list", n_iter)
resampling_meta = vector("list", n_iter)
last_nonphyloPCAresult = NULL
last_nonphyloPCAresult_varimax = NULL
last_nonphyloPCAData_meta_rmna = NULL

for (i in seq_len(n_iter)) {
	message("Resampling iteration: ", i, "/", n_iter)

	# Resample one individual per species within each site
	traitDataPCA_iter = traitDataPCA_touse %>%
		sample_one_ind_per_sp_site(trait_cols = traitName_touse) %>%
		mutate(PCAIdentifier_spgf = paste(gsub(" ","_",SpeciesFullName), SiteID, SampleID, sep="-"))

	# Curate data for PCA analysis
	nonphyloPCAData_numonly = traitDataPCA_iter %>% ungroup() %>%
		dplyr::select(PCAIdentifier_spgf, SiteID, all_of(traitName_touse)) %>%
		dplyr::mutate(SiteID = as.factor(SiteID)) %>%
		na.omit()

	# scale traits based on columns (mean=0, SD=1)
	nonphyloPCAData_numonly[, traitName_touse] = apply(nonphyloPCAData_numonly[, traitName_touse], 2, scale)
	rownames(nonphyloPCAData_numonly) = nonphyloPCAData_numonly$PCAIdentifier_spgf

	# used for site annotation in PCA biplot
	nonphyloPCAData_meta_rmna = nonphyloPCAData_numonly
	nonphyloPCAData_numonly = nonphyloPCAData_numonly %>% dplyr::select(-PCAIdentifier_spgf, -SiteID)

	#-----------------------------------------------------------
	# Perform PCA using FactoMineR
	nonphyloPCAresult = PCA(nonphyloPCAData_numonly, scale.unit = TRUE, graph = FALSE)

	#-----------------------------------------------------------
	# Varimax Rotation 
	varimax_result <- varimax(nonphyloPCAresult$var$coord)
	rotated_loadings <- as.matrix(varimax_result$loadings)
	rotated_scores <- nonphyloPCAresult$ind$coord %*% varimax_result$rotmat
	nonphyloPCAresult_varimax_rotated <- nonphyloPCAresult
	nonphyloPCAresult_varimax_rotated$var$coord <- rotated_loadings
	nonphyloPCAresult_varimax_rotated$ind$coord <- rotated_scores

	#-----------------------------------------------------------
	# Store outputs
	var_coords = as.data.frame(nonphyloPCAresult$var$coord[, 1:5, drop = FALSE])
	var_coords$trait = rownames(var_coords)
	var_coords$iteration = i

	var_coords_varimax = as.data.frame(nonphyloPCAresult_varimax_rotated$var$coord[, 1:5, drop = FALSE])
	var_coords_varimax$trait = rownames(var_coords_varimax)
	var_coords_varimax$iteration = i

	eig_df = as.data.frame(nonphyloPCAresult$eig)
	colnames(eig_df) = c("eigenvalue", "variance_percent", "cumulative_percent")
	eig_df$axis = seq_len(nrow(eig_df))
	eig_df$iteration = i

	resampling_varcoords[[i]] = var_coords
	resampling_varcoords_varimax[[i]] = var_coords_varimax
	resampling_eig[[i]] = eig_df
	resampling_meta[[i]] = nonphyloPCAData_meta_rmna %>%
		dplyr::select(PCAIdentifier_spgf, SiteID) %>%
		dplyr::mutate(iteration = i)
}

### COMBINE RESULTS ###
resampling_varcoords_df = dplyr::bind_rows(resampling_varcoords)
resampling_varcoords_varimax_df = dplyr::bind_rows(resampling_varcoords_varimax)
resampling_eig_df = dplyr::bind_rows(resampling_eig)
resampling_meta_df = dplyr::bind_rows(resampling_meta)

### MEAN LOADINGS (PC1-PC5) ###
mean_loadings_pc1_5 = resampling_varcoords_varimax_df %>%
	dplyr::select(trait, Dim.1, Dim.2, Dim.3, Dim.4, Dim.5) %>%
	dplyr::group_by(trait) %>%
	dplyr::summarise(
		PC1 = mean(Dim.1, na.rm = TRUE),
		PC2 = mean(Dim.2, na.rm = TRUE),
		PC3 = mean(Dim.3, na.rm = TRUE),
		PC4 = mean(Dim.4, na.rm = TRUE),
		PC5 = mean(Dim.5, na.rm = TRUE),
		.n = dplyr::n(),
		.groups = "drop"
	)

# SAVE MEAN LOADINGS AFTER VARIMAX ROTATION, PC1-5
write.csv(
	mean_loadings_pc1_5,
	file = "revision-code/step1-1indv-nonphyPCA-PC1to5-LES.csv",
	row.names = FALSE
)
