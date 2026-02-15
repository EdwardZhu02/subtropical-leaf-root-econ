###########################################################################
# One-individual-per-species resampling for PCA angle calculation
# Repeatedly sample one individual per species within each community,
# recompute nonphylo- and phylo-informed PCA, and summarize trait-pair angles
###########################################################################

rm(list = ls())
library(tidyverse)
library(dplyr)
library(vegan)
library(reshape2)
library(ggpmisc) # stat_poly_eq
library(patchwork)

# phylo-PCA dependencies
library(ape)
library(geiger)
library(treeplyr)
library(phytools)

load("individual-level-code/traitDataFujian-Ind-step1.RData")
load("species-level-code/traitDataFujian-spavg-phylo-step2.RData")

### SET DATA TO USE ###
traitData_touse = traitDataIndv_SelectedTraits_log

### TRAITS TO USE ###
trait_sets = list(
	RESLES = c("RTD","SRL","RD","RNC","LMA","LNC"),
	RESLES_RDMC = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA"),
	RESLES_LPC = c("RTD","SRL","RD","RNC","LMA","LNC","LPC"),
	RESLES_RDMCLPC = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA","LPC")
)

### RESAMPLING SETTINGS ###
n_iter = 200
set.seed(1)

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------
func_calc_angle_pcavars <- function(x, y){
	norm_x <- sqrt(sum(x^2))
	norm_y <- sqrt(sum(y^2))
	cosXY <- round(as.numeric((x %*% y) / (norm_x * norm_y)), 8)
	angle <- acos(cosXY) * 180 / pi
	return(angle)
}

func_get_upper_tri <- function(cormat){
	cormat[lower.tri(cormat)] = NA
	return(cormat)
}

sample_one_ind_per_sp_site <- function(data, trait_cols) {
	data %>%
		dplyr::filter(if_all(all_of(trait_cols), ~ !is.na(.))) %>%
		group_by(SiteID, SpeciesFullName) %>%
		slice_sample(n = 1) %>%
		ungroup()
}

calc_angles_nonphy <- function(df, trait_cols, dimensions = 3) {
	aux = df %>% dplyr::select(all_of(trait_cols)) %>% na.omit()
	if (nrow(aux) < 3) return(NULL)
	aux_scaled = scale(aux)
	PCA_aux = tryCatch(princomp(aux_scaled), error = function(e) NULL)
	if (is.null(PCA_aux)) return(NULL)

	load_mat = as.matrix(PCA_aux$loadings)
	load_mat = load_mat[, seq_len(min(dimensions, ncol(load_mat))), drop = FALSE]
	angles = matrix(NA, nrow = ncol(aux), ncol = ncol(aux),
									dimnames = list(colnames(aux), colnames(aux)))
	for (i in seq_len(ncol(angles))) {
		for (j in seq_len(ncol(angles))) {
			angles[i, j] <- func_calc_angle_pcavars(load_mat[i, ], load_mat[j, ])
		}
	}
	angles
}

calc_angles_phylo <- function(df_species, trait_cols, tree, dimensions = 3) {
	aux = df_species %>%
		dplyr::select(SpeciesFullName, all_of(trait_cols)) %>%
		dplyr::filter(if_all(all_of(trait_cols), ~ !is.na(.))) %>%
		dplyr::mutate(tipLabelMatched = gsub(" ", "_", SpeciesFullName))
	if (nrow(aux) < 3) return(NULL)

	aux_tree = tree
	aux_tree$node.label = NULL
	aux_treedata = tryCatch(
		make.treedata(tree = aux_tree, data = aux, name_column = "tipLabelMatched"),
		error = function(e) NULL
	)
	if (is.null(aux_treedata) || nrow(aux_treedata$dat) < 3) return(NULL)

	trait_mat = aux_treedata$dat %>% dplyr::select(all_of(trait_cols))
	trait_mat = scale(trait_mat)
	rownames(trait_mat) = aux_treedata$dat$tipLabelMatched

	Phyl_PCA_aux = tryCatch(
		phyl.pca(aux_treedata$phy, trait_mat, mode = "corr", method = "lambda"),
		error = function(e) NULL
	)
	if (is.null(Phyl_PCA_aux)) return(NULL)

	load_mat = Phyl_PCA_aux$L
	load_mat = load_mat[, seq_len(min(dimensions, ncol(load_mat))), drop = FALSE]
	angles = matrix(NA, nrow = ncol(trait_mat), ncol = ncol(trait_mat),
									dimnames = list(colnames(trait_mat), colnames(trait_mat)))
	for (i in seq_len(ncol(angles))) {
		for (j in seq_len(ncol(angles))) {
			angles[i, j] <- func_calc_angle_pcavars(load_mat[i, ], load_mat[j, ])
		}
	}
	angles
}

melt_angles <- function(angle_mat, angle_name) {
	if (is.null(angle_mat)) return(NULL)
	angle_uptri = func_get_upper_tri(angle_mat)
	reshape2::melt(angle_uptri, na.rm = TRUE) %>%
		dplyr::filter(value != 0) %>%
		dplyr::rename(trait1 = Var1, trait2 = Var2, !!angle_name := value)
}

annotate_pairs <- function(df) {
	df %>%
		dplyr::mutate(pair_name = paste0(trait1, "-", trait2)) %>%
		dplyr::mutate(pair_annotation = ifelse(
			pair_name %in% c("SRL-RD","RD-SRR25","SRL-SRR25","SRL-SRA","RD-SRA"), "RES collaboration",
			ifelse(pair_name %in% c("RTD-RNC","RNC-RDMC","RTD-RDMC"), "RES conservation",
						 ifelse(pair_name %in% c("LMA-LNC","LMA-LPC"), "LES", "other")))
		)
}

# ---------------------------------------------------------------------------
# Resampling loop
# ---------------------------------------------------------------------------
tree_touse = phylotree_result[[1]]

angle_results = vector("list", n_iter)

for (i in seq_len(n_iter)) {
	message("Resampling iteration: ", i, "/", n_iter)

	iter_list = list()
	idx = 1

	for (trait_name in names(trait_sets)) {
		trait_cols = trait_sets[[trait_name]]

		data_iter = traitData_touse %>%
			sample_one_ind_per_sp_site(trait_cols = trait_cols)

		# non-phylo angles from sampled individuals
		angles_nonphy = calc_angles_nonphy(data_iter, trait_cols)
		df_nonphy = melt_angles(angles_nonphy, "angle")
		if (!is.null(df_nonphy)) {
			iter_list[[idx]] = df_nonphy %>%
				dplyr::mutate(iteration = i, method = "nonphy", trait_set = trait_name)
			idx = idx + 1
		}

		# phylo angles from species means of sampled individuals
		data_iter_sp = data_iter %>%
			group_by(SpeciesFullName) %>%
			summarise(across(all_of(trait_cols), mean, na.rm = TRUE), .groups = "drop")
		angles_phylo = calc_angles_phylo(data_iter_sp, trait_cols, tree_touse)
		df_phylo = melt_angles(angles_phylo, "angle")
		if (!is.null(df_phylo)) {
			iter_list[[idx]] = df_phylo %>%
				dplyr::mutate(iteration = i, method = "phylo", trait_set = trait_name)
			idx = idx + 1
		}
	}

	angle_results[[i]] = dplyr::bind_rows(iter_list)
}

angle_results_raw = dplyr::bind_rows(angle_results) %>%
	annotate_pairs()

# ---------------------------------------------------------------------------
# Summaries (mean and 5/95% CI for each group)
# ---------------------------------------------------------------------------
angle_group_iter = angle_results_raw %>%
	group_by(iteration, method, trait_set, pair_annotation) %>%
	summarise(mean_angle = mean(angle, na.rm = TRUE), .groups = "drop")

angle_group_summary = angle_results_raw %>%
	group_by(method, trait_set, pair_annotation) %>%
	summarise(
		mean_angle = mean(angle, na.rm = TRUE),
		CI_5 = quantile(angle, 0.05, na.rm = TRUE),
		CI_95 = quantile(angle, 0.95, na.rm = TRUE),
		n = n(),
		.groups = "drop"
	)

write.csv(
	angle_results_raw,
	file = "revision-code/step42-angle-raw.csv",
	row.names = FALSE
)

write.csv(
	angle_group_summary,
	file = "revision-code/step42-angle-summary.csv",
	row.names = FALSE
)

# Pair-level summary for plotting nonphy vs phylo
angle_pair_summary = angle_results_raw %>%
	group_by(method, trait_set, trait1, trait2, pair_annotation) %>%
	summarise(
		mean_angle = mean(angle, na.rm = TRUE),
		CI_5 = quantile(angle, 0.05, na.rm = TRUE),
		CI_95 = quantile(angle, 0.95, na.rm = TRUE),
		.groups = "drop"
	)

angle_pair_wide = angle_pair_summary %>%
	tidyr::pivot_wider(
		names_from = method,
		values_from = c(mean_angle, CI_5, CI_95),
		names_glue = "{method}_{.value}"
	)

# ---------------------------------------------------------------------------
# Visualization: nonphy vs phylo (scatter) for each trait set
# ---------------------------------------------------------------------------
plot_scatter_template <- function(df, xlab, ylab) {
	ggplot(data = df, aes(x = nonphy_mean_angle, y = phylo_mean_angle)) +
		geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
		geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
		geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") +
		geom_point(aes(fill = pair_annotation), size = 2.7, alpha = 0.9, shape = 21, color = "gray15", na.rm = TRUE) +
		geom_segment(aes(x = nonphy_CI_5, xend = nonphy_CI_95,
								 y = phylo_mean_angle, yend = phylo_mean_angle),
				 linewidth = 0.4, color = "gray15", alpha = 0.8, na.rm = TRUE) +
		geom_segment(aes(x = nonphy_mean_angle, xend = nonphy_mean_angle,
								 y = phylo_CI_5, yend = phylo_CI_95),
				 linewidth = 0.4, color = "gray15", alpha = 0.8, na.rm = TRUE) +
		stat_poly_eq(formula = y ~ x,
								 aes(label = paste0("P", paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=", ..p.value..)), ..rr.label.., sep = "~~"))),
								 parse = TRUE, rr.digits = 4) +
		scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) +
		theme_classic() + theme(plot.background = element_blank()) +
		labs(x = xlab, y = ylab) +
		coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))
}

plt_anglecorr_RESLES = plot_scatter_template(
	angle_pair_wide %>% filter(trait_set == "RESLES"),
	"Backbone", "Backbone (phylo)"
)

plt_anglecorr_RESLES_RDMC = plot_scatter_template(
	angle_pair_wide %>% filter(trait_set == "RESLES_RDMC"),
	"Backbone + RES co-pred", "Backbone + RES co-pred (phylo)"
)

plt_anglecorr_RESLES_LPC = plot_scatter_template(
	angle_pair_wide %>% filter(trait_set == "RESLES_LPC"),
	"Backbone + LES co-pred", "Backbone + LES co-pred (phylo)"
)

plt_anglecorr_RESLES_RDMCLPC = plot_scatter_template(
	angle_pair_wide %>% filter(trait_set == "RESLES_RDMCLPC"),
	"Backbone + all co-pred", "Backbone + all co-pred (phylo)"
)

plt_anglecorr_total = plt_anglecorr_RESLES + plt_anglecorr_RESLES_RDMC +
	plt_anglecorr_RESLES_LPC + plt_anglecorr_RESLES_RDMCLPC +
	plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.title = element_blank())

ggsave(plot = plt_anglecorr_total, filename = "revision-code/step42plt_anglecorr.pdf", width = 5.5, height = 6)

# ---------------------------------------------------------------------------
# Visualization: group mean angles with 5/95% CI
# ---------------------------------------------------------------------------
angle_group_summary$pair_annotation = factor(
	angle_group_summary$pair_annotation,
	levels = c("RES collaboration", "RES conservation", "LES", "other")
)
