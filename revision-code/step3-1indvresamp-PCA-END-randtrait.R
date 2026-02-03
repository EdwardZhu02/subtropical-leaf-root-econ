#+#+#+#+######################################################################
# One-individual-per-species resampling for PCA END (random trait)
# Repeatedly sample one individual per species within each community
# and recompute END using both phylo and non-phylo methods
#+#+#+#+######################################################################

rm(list=ls())
library(tidyverse)
library(dplyr)
library(vegan)

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
# Backbone (RES + LES) and co-predictor trait sets follow Fig. S6 logic
trait_sets = list(
	RES = c("RTD","SRL","RD","RNC"),
	RES_LMA = c("RTD","SRL","RD","RNC","LMA"),
	RES_LNC = c("RTD","SRL","RD","RNC","LNC"),
	RLES = c("RTD","SRL","RD","RNC","LMA","LNC"),
	RLES_LPC = c("RTD","SRL","RD","RNC","LMA","LNC","LPC"),
	RLES_RESco = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRA","SRR25"),
	RLES_allco = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRA","SRR25","LPC")
)

# Random traits only for RES and RES+LES backbone (validation)
random_trait_sets = list(
	RES_Rand = c("RTD","SRL","RD","RNC"),
	RLES_Rand = c("RTD","SRL","RD","RNC","LMA","LNC")
)

### RESAMPLING SETTINGS ###
n_iter = 200
set.seed(1)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
sample_one_ind_per_sp_site <- function(data, trait_cols) {
	data %>%
		dplyr::filter(if_all(all_of(trait_cols), ~ !is.na(.))) %>%
		group_by(SiteID, SpeciesFullName) %>%
		slice_sample(n = 1) %>%
		ungroup()
}

calc_end_nonphy <- function(df, trait_cols) {
	aux = df %>% dplyr::select(all_of(trait_cols)) %>% na.omit()
	if (nrow(aux) < 3) return(NA_real_)
	aux_scaled = scale(aux)
	PCA_aux = tryCatch(princomp(aux_scaled), error = function(e) NULL)
	if (is.null(PCA_aux)) return(NA_real_)
	vegan::diversity(x = apply(PCA_aux$scores, 2, var), index = "invsimpson")
}

calc_end_phylo <- function(df_species, trait_cols, tree) {
	aux = df_species %>%
		dplyr::select(SpeciesFullName, all_of(trait_cols)) %>%
		dplyr::filter(if_all(all_of(trait_cols), ~ !is.na(.))) %>%
		dplyr::mutate(tipLabelMatched = gsub(" ", "_", SpeciesFullName))
	if (nrow(aux) < 3) return(NA_real_)

	aux_tree = tree
	aux_tree$node.label = NULL
	aux_treedata = tryCatch(
		make.treedata(tree = aux_tree, data = aux, name_column = "tipLabelMatched"),
		error = function(e) NULL
	)
	if (is.null(aux_treedata) || nrow(aux_treedata$dat) < 3) return(NA_real_)

	trait_mat = aux_treedata$dat %>% dplyr::select(all_of(trait_cols))
	trait_mat = scale(trait_mat)
	rownames(trait_mat) = aux_treedata$dat$tipLabelMatched

	Phyl_PCA_aux = tryCatch(
		phyl.pca(aux_treedata$phy, trait_mat, mode = "corr", method = "lambda"),
		error = function(e) NULL
	)
	if (is.null(Phyl_PCA_aux)) return(NA_real_)
	vegan::diversity(x = diag(Phyl_PCA_aux$Eval), index = "invsimpson")
}

# -----------------------------------------------------------------------------
# Resampling loop
# -----------------------------------------------------------------------------
tree_touse = phylotree_result[[1]]

results = vector("list", n_iter)

for (i in seq_len(n_iter)) {
	message("Resampling iteration: ", i, "/", n_iter)

	iter_results = vector("list", (length(trait_sets) + length(random_trait_sets)) * 2)
	idx = 1

	for (trait_name in names(trait_sets)) {
		trait_cols = trait_sets[[trait_name]]
		trait_label = trait_name

		data_iter = traitData_touse %>%
			sample_one_ind_per_sp_site(trait_cols = trait_cols)
		end_nonphy = calc_end_nonphy(data_iter, trait_cols)

		data_iter_sp = data_iter %>%
			group_by(SpeciesFullName) %>%
			summarise(across(all_of(trait_cols), mean, na.rm = TRUE), .groups = "drop")
		end_phy = calc_end_phylo(data_iter_sp, trait_cols, tree_touse)

		iter_results[[idx]] = tibble::tibble(
			iteration = i,
			method = "nonphy",
			trait_set = trait_label,
			END = end_nonphy
		)
		idx = idx + 1
		iter_results[[idx]] = tibble::tibble(
			iteration = i,
			method = "phylo",
			trait_set = trait_label,
			END = end_phy
		)
		idx = idx + 1
	}

	for (trait_name in names(random_trait_sets)) {
		trait_cols = random_trait_sets[[trait_name]]
		trait_label = trait_name

		data_iter = traitData_touse %>%
			sample_one_ind_per_sp_site(trait_cols = trait_cols) %>%
			mutate(Rand = rnorm(nrow(.)))
		end_nonphy = calc_end_nonphy(data_iter, c(trait_cols, "Rand"))

		data_iter_sp = data_iter %>%
			group_by(SpeciesFullName) %>%
			summarise(across(all_of(trait_cols), mean, na.rm = TRUE), .groups = "drop") %>%
			mutate(Rand = rnorm(nrow(.)))
		end_phy = calc_end_phylo(data_iter_sp, c(trait_cols, "Rand"), tree_touse)

		iter_results[[idx]] = tibble::tibble(
			iteration = i,
			method = "nonphy",
			trait_set = trait_label,
			END = end_nonphy
		)
		idx = idx + 1
		iter_results[[idx]] = tibble::tibble(
			iteration = i,
			method = "phylo",
			trait_set = trait_label,
			END = end_phy
		)
		idx = idx + 1
	}

	results[[i]] = dplyr::bind_rows(iter_results)
}

END_results_raw = dplyr::bind_rows(results)

END_results_summary = END_results_raw %>%
	group_by(method, trait_set) %>%
	summarise(
		mean_END = mean(END, na.rm = TRUE),
		CI_5 = quantile(END, 0.05, na.rm = TRUE),
		CI_95 = quantile(END, 0.95, na.rm = TRUE),
		n = sum(!is.na(END)),
		.groups = "drop"
	)

write.csv(
	END_results_raw,
	file = "revision-code/step3-PCA-END-randtrait_raw.csv",
	row.names = FALSE
)

write.csv(
	END_results_summary,
	file = "revision-code/step3-PCA-END-randtrait_summary.csv",
	row.names = FALSE
)

# ----------------------------------------------------------------------------
# Visualization (similar to step33)
# ----------------------------------------------------------------------------
plot_df = END_results_summary %>%
	mutate(
		comb_group = ifelse(method == "nonphy", "non-phylo", "phylo")
	)

# Panel A: RES backbone + LES traits
panelA_sets = c("RES","RES_LMA","RES_LNC","RLES")
panelA_labels = c(
	RES = "RES (4)",
	RES_LMA = "RES + LMA (5)",
	RES_LNC = "RES + LNC (5)",
	RLES = "Core RES + LES backbone (6)"
)

plot_df_A = plot_df %>%
	filter(trait_set %in% panelA_sets) %>%
	mutate(combination = factor(panelA_labels[trait_set], levels = panelA_labels[panelA_sets]))

rand_RES_nonphy = plot_df %>% filter(trait_set == "RES_Rand", method == "nonphy") %>% pull(mean_END)
rand_RES_phy = plot_df %>% filter(trait_set == "RES_Rand", method == "phylo") %>% pull(mean_END)

plt_END_RES_backbone = ggplot(plot_df_A, aes(x = comb_group, y = mean_END, fill = combination)) +
	geom_hline(yintercept = 0:5, linewidth = 0.2, color = "gray") +
	geom_bar(stat = "identity", width = 0.6, position = position_dodge(width = 0.65), color = "black", linewidth = 0.2) +
	geom_errorbar(aes(ymin = CI_5, ymax = CI_95), width = 0.2, position = position_dodge(width = 0.65), linewidth = 0.3) +
	geom_segment(aes(x = 0.5, xend = 1.4, y = rand_RES_nonphy, yend = rand_RES_nonphy),
		linewidth = 0.7, color = "darkred", linetype = "dashed") +
	geom_segment(aes(x = 1.6, xend = 2.5, y = rand_RES_phy, yend = rand_RES_phy),
		linewidth = 0.7, color = "darkblue", linetype = "dashed") +
	geom_text(aes(label = signif(mean_END, digits = 3)), fontface = "bold",
		position = position_dodge(width = 0.65), vjust = -0.4, size = 3.2) +
	scale_fill_manual(values = c("#1999B2", "#95BCE5", "#E84445", "#F39DA0")) +
	labs(y = "Effective Number of Dimensions (END)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.title = element_blank())

ggsave(plot = plt_END_RES_backbone, filename = "revision-code/step3plt_END_RES_backbone_1ind.pdf", width = 7.5, height = 3.5)

# Panel B: Backbone + co-predictors
panelB_sets = c("RLES","RLES_LPC","RLES_RESco","RLES_allco")
panelB_labels = c(
	RLES = "Backbone (6)",
	RLES_LPC = "Backbone + LES co-predictors (7)",
	RLES_RESco = "Backbone + RES co-predictors (9)",
	RLES_allco = "Backbone + all co-predictors (10)"
)

plot_df_B = plot_df %>%
	filter(trait_set %in% panelB_sets) %>%
	mutate(combination = factor(panelB_labels[trait_set], levels = panelB_labels[panelB_sets]))

rand_RLES_nonphy = plot_df %>% filter(trait_set == "RLES_Rand", method == "nonphy") %>% pull(mean_END)
rand_RLES_phy = plot_df %>% filter(trait_set == "RLES_Rand", method == "phylo") %>% pull(mean_END)

plt_END_RLES_copred = ggplot(plot_df_B, aes(x = comb_group, y = mean_END, fill = combination)) +
	geom_hline(yintercept = 0:5, linewidth = 0.2, color = "gray") +
	geom_bar(stat = "identity", width = 0.6, position = position_dodge(width = 0.65), color = "black", linewidth = 0.2) +
	geom_errorbar(aes(ymin = CI_5, ymax = CI_95), width = 0.2, position = position_dodge(width = 0.65), linewidth = 0.3) +
	geom_segment(aes(x = 0.5, xend = 1.4, y = rand_RLES_nonphy, yend = rand_RLES_nonphy),
		linewidth = 0.7, color = "darkred", linetype = "dashed") +
	geom_segment(aes(x = 1.6, xend = 2.5, y = rand_RLES_phy, yend = rand_RLES_phy),
		linewidth = 0.7, color = "darkblue", linetype = "dashed") +
	geom_text(aes(label = signif(mean_END, digits = 3)), fontface = "bold",
		position = position_dodge(width = 0.65), vjust = -0.4, size = 3.2) +
	scale_fill_manual(values = c("#1999B2", "#95BCE5", "#E84445", "#F39DA0")) +
	labs(y = "Effective Number of Dimensions (END)") +
	theme_classic() +
	theme(axis.title.x = element_blank(), legend.title = element_blank())

ggsave(plot = plt_END_RLES_copred, filename = "revision-code/step3plt_END_RLES_copred_1ind.pdf", width = 7, height = 3.5)