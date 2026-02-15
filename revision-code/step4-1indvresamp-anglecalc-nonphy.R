################################################################################
# One-individual-per-species resampling for non-phylo PCA angles
# Repeatedly sample one individual per species (within each community)
# and recompute non-phylo PCA and trait-pair angles
################################################################################

rm(list=ls())
library(tidyverse)
library(dplyr)
library(vegan)
library(ggpmisc) # stat_poly_eq
library(reshape2)
library(patchwork)

load("individual-level-code/traitDataFujian-Ind-step1.RData")

### SET DATA TO USE ###
traitData_touse = traitDataIndv_SelectedTraits_log

# harmonize species name column if needed
if (!"SpeciesFullName" %in% names(traitData_touse) && "speciesFullName" %in% names(traitData_touse)) {
	traitData_touse = traitData_touse %>% dplyr::rename(SpeciesFullName = speciesFullName)
}

### TRAITS TO USE ###
# Backbone (RES + LES) and co-predictor trait sets follow Fig. S6 logic
trait_sets = list(
	Backbone = c("RTD","SRL","RD","RNC","LMA","LNC"),
	Backbone_RESco = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA"),
	Backbone_LESco = c("RTD","SRL","RD","RNC","LMA","LNC","LPC"),
	Backbone_Allco = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA","LPC")
)

### RESAMPLING SETTINGS ###
n_iter = 200
set.seed(1)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
func_calc_angle_pcavars <- function(x, y){
	norm_x <- sqrt(sum(x^2))
	norm_y <- sqrt(sum(y^2))
	cosXY <- round(as.numeric((x %*% y) / (norm_x * norm_y)), 8)
	angle <- acos(cosXY) * 180 / pi
	return(angle)
}

sample_one_ind_per_sp_site <- function(data, trait_cols) {
	data %>%
		dplyr::filter(if_all(all_of(trait_cols), ~ !is.na(.))) %>%
		group_by(SiteID, SpeciesFullName) %>%
		slice_sample(n = 1) %>%
		ungroup()
}

calc_nonphy_angles <- function(df, trait_cols, dimensions = 3) {
	aux = df %>% dplyr::select(all_of(trait_cols)) %>% na.omit()
	if (nrow(aux) < 3) return(NULL)
	aux_scaled = scale(aux)

	pwcorr = cor(aux_scaled, use = "pairwise.complete")
	AllTraits = eigen(pwcorr)
	dimnames(AllTraits$vectors) = list(rownames(pwcorr), paste0("PC", 1:ncol(AllTraits$vectors)))
	AllTraits$CumVariance = cumsum(AllTraits$values / sum(AllTraits$values))
	AllTraits$dimensions = dimensions

	Eival = AllTraits$values
	loadingsEig = AllTraits$vectors
	for(i in 1:ncol(loadingsEig)){
		loadingsEig[, i] = loadingsEig[, i]*sqrt(Eival[i])
	}
	AllTraits$loadings = loadingsEig

	loadEigVar = AllTraits$loadings[, 1:dimensions]
	angles <- matrix(NA, nrow = nrow(loadEigVar), ncol = nrow(loadEigVar),
		dimnames = list(rownames(loadEigVar), rownames(loadEigVar)))
	for(i in 1:nrow(loadEigVar)){
		for(j in 1:nrow(loadEigVar)){
			angles[i, j] <- func_calc_angle_pcavars(loadEigVar[i, ], loadEigVar[j, ])
		}
	}
	return(angles)
}

func_get_upper_tri = function(cormat){
	cormat[lower.tri(cormat)] = NA
	return(cormat)
}

annotate_pair = function(trait1, trait2) {
	key = paste(sort(c(trait1, trait2)), collapse = "-")
	if (key %in% c("RD-SRL", "SRL-RD")) return("RES collaboration")
	if (key %in% c("RNC-RTD", "RTD-RNC")) return("RES conservation")
	if (key %in% c("LMA-LNC", "LNC-LMA")) return("LES")
	return("other")
}

# ------------------------------------------------------------------------------
# Resampling loop
# ------------------------------------------------------------------------------
results = vector("list", n_iter)

for (i in seq_len(n_iter)) {
	message("Resampling iteration: ", i, "/", n_iter)
	iter_results = vector("list", length(trait_sets))
	idx = 1

	for (trait_name in names(trait_sets)) {
		trait_cols = trait_sets[[trait_name]]

		data_iter = traitData_touse %>%
			sample_one_ind_per_sp_site(trait_cols = trait_cols)

		angle_mat = calc_nonphy_angles(data_iter, trait_cols)
		if (is.null(angle_mat)) {
			iter_results[[idx]] = NULL
			idx = idx + 1
			next
		}

		angle_uptri = func_get_upper_tri(angle_mat)
		angle_long = reshape2::melt(angle_uptri, na.rm = TRUE) %>%
			dplyr::filter(value != 0) %>%
			dplyr::rename(trait1 = Var1, trait2 = Var2, angle = value) %>%
			dplyr::mutate(
				iteration = i,
				trait_set = trait_name
			)

		iter_results[[idx]] = angle_long
		idx = idx + 1
	}

	results[[i]] = dplyr::bind_rows(iter_results)
}

pair_palette = c(
	"LES" = "#71BFB2",
	"other" = "gray70",
	"RES collaboration" = "#237B9F",
	"RES conservation" = "#AD0B08"
)

angles_raw = dplyr::bind_rows(results) %>%
	dplyr::mutate(
		pair_name = paste(trait1, trait2, sep = "-"),
		pair_annotation = pmap_chr(list(trait1, trait2), annotate_pair),
		pair_annotation = factor(pair_annotation, levels = names(pair_palette))
	)

# ------------------------------------------------------------------------------
# Summaries: by trait pair and by group
# ------------------------------------------------------------------------------
angles_pair_summary = angles_raw %>%
	group_by(trait_set, trait1, trait2, pair_name, pair_annotation) %>%
	summarise(
		mean_angle = mean(angle, na.rm = TRUE),
		CI_5 = quantile(angle, 0.05, na.rm = TRUE),
		CI_95 = quantile(angle, 0.95, na.rm = TRUE),
		n = sum(!is.na(angle)),
		.groups = "drop"
	)

angles_group_summary = angles_raw %>%
	group_by(trait_set, pair_annotation) %>%
	summarise(
		mean_angle = mean(angle, na.rm = TRUE),
		CI_5 = quantile(angle, 0.05, na.rm = TRUE),
		CI_95 = quantile(angle, 0.95, na.rm = TRUE),
		n = sum(!is.na(angle)),
		.groups = "drop"
	)

write.csv(
	angles_raw,
	file = "revision-code/step41-angle-raw.csv",
	row.names = FALSE
)

write.csv(
	angles_pair_summary,
	file = "revision-code/step41-angle-summary.csv",
	row.names = FALSE
)

write.csv(
	angles_group_summary,
	file = "revision-code/step41-angle-summary_group_1ind.csv",
	row.names = FALSE
)

# ------------------------------------------------------------------------------
# Visualization: angle correlation plots (mean angles)
# ------------------------------------------------------------------------------
angles_wide = angles_pair_summary %>%
	dplyr::select(trait1, trait2, pair_name, pair_annotation, trait_set, mean_angle, CI_5, CI_95) %>%
	tidyr::pivot_wider(
		names_from = trait_set,
		values_from = c(mean_angle, CI_5, CI_95),
		names_glue = "{trait_set}_{.value}"
	)

angles_wide_resco = angles_wide %>%
	dplyr::filter(is.finite(Backbone_mean_angle), is.finite(Backbone_RESco_mean_angle))

angles_wide_lpc = angles_wide %>%
	dplyr::filter(is.finite(Backbone_mean_angle), is.finite(Backbone_LESco_mean_angle))

angles_wide_all = angles_wide %>%
	dplyr::filter(is.finite(Backbone_mean_angle), is.finite(Backbone_Allco_mean_angle))

plt_anglecorr_addRESco = ggplot(data = angles_wide_resco,
		aes(x = Backbone_mean_angle, y = Backbone_RESco_mean_angle)) +
	geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
	geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
	geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") +
  geom_point(aes(fill = pair_annotation), size = 2.7, alpha = 0.9, shape = 21, color = "gray15", na.rm = TRUE) +
	geom_segment(aes(x = Backbone_CI_5, xend = Backbone_CI_95,
		y = Backbone_RESco_mean_angle, yend = Backbone_RESco_mean_angle),
		linewidth = 0.4, color = "gray15", alpha = 0.8, na.rm = TRUE) +
	geom_segment(aes(x = Backbone_mean_angle, xend = Backbone_mean_angle,
		y = Backbone_RESco_CI_5, yend = Backbone_RESco_CI_95),
		linewidth = 0.4, color = "gray15", alpha = 0.8, na.rm = TRUE) +
	stat_poly_eq(formula = y ~ x,
		aes(label = paste0("P", paste(ifelse(after_stat(p.value) < 0.0001, "<0.0001", paste0("=", signif(after_stat(p.value), 4))), after_stat(rr.label), sep = "~~"))),
		parse = TRUE, rr.digits = 4, na.rm = TRUE) +
	scale_fill_manual(values = pair_palette, drop = FALSE) +
	theme_classic() + theme(plot.background = element_blank()) +
	labs(x = "Backbone", y = "Backbone + RES co-pred") +
	coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

plt_anglecorr_addLPC = ggplot(data = angles_wide_lpc,
		aes(x = Backbone_mean_angle, y = Backbone_LESco_mean_angle)) +
	geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
	geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
	geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") +
  geom_point(aes(fill = pair_annotation), size = 2.7, alpha = 0.9, shape = 21, color = "gray15", na.rm = TRUE) +
	geom_segment(aes(x = Backbone_CI_5, xend = Backbone_CI_95,
		y = Backbone_LESco_mean_angle, yend = Backbone_LESco_mean_angle),
		linewidth = 0.4, color = "gray15", alpha = 0.8, na.rm = TRUE) +
	geom_segment(aes(x = Backbone_mean_angle, xend = Backbone_mean_angle,
		y = Backbone_LESco_CI_5, yend = Backbone_LESco_CI_95),
		linewidth = 0.4, color = "gray15", alpha = 0.8, na.rm = TRUE) +
	stat_poly_eq(formula = y ~ x,
		aes(label = paste0("P", paste(ifelse(after_stat(p.value) < 0.0001, "<0.0001", paste0("=", signif(after_stat(p.value), 4))), after_stat(rr.label), sep = "~~"))),
		parse = TRUE, rr.digits = 4, na.rm = TRUE) +
	scale_fill_manual(values = pair_palette, drop = FALSE) +
	theme_classic() + theme(plot.background = element_blank()) +
	labs(x = "Backbone", y = "Backbone + LES co-pred") +
	coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

plt_anglecorr_addAll = ggplot(data = angles_wide_all,
		aes(x = Backbone_mean_angle, y = Backbone_Allco_mean_angle)) +
	geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
	geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
	geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") +
  geom_point(aes(fill = pair_annotation), size = 2.7, alpha = 0.9, shape = 21, color = "gray15", na.rm = TRUE) +
	geom_segment(aes(x = Backbone_CI_5, xend = Backbone_CI_95,
		y = Backbone_Allco_mean_angle, yend = Backbone_Allco_mean_angle),
		linewidth = 0.4, color = "gray15", alpha = 0.8, na.rm = TRUE) +
	geom_segment(aes(x = Backbone_mean_angle, xend = Backbone_mean_angle,
		y = Backbone_Allco_CI_5, yend = Backbone_Allco_CI_95),
		linewidth = 0.4, color = "gray15", alpha = 0.8, na.rm = TRUE) +
	stat_poly_eq(formula = y ~ x,
		aes(label = paste0("P", paste(ifelse(after_stat(p.value) < 0.0001, "<0.0001", paste0("=", signif(after_stat(p.value), 4))), after_stat(rr.label), sep = "~~"))),
		parse = TRUE, rr.digits = 4, na.rm = TRUE) +
	scale_fill_manual(values = pair_palette, drop = FALSE) +
	theme_classic() + theme(plot.background = element_blank()) +
	labs(x = "Backbone", y = "Backbone + all co-pred") +
	coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

plt_anglecorr_nonphy_total = plt_anglecorr_addRESco + plt_anglecorr_addLPC + plt_anglecorr_addAll +
	plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.title = element_blank())

ggsave(plot = plt_anglecorr_nonphy_total,
	filename = "revision-code/step42plt_anglecorr.pdf",
	width = 8, height = 3.2)
