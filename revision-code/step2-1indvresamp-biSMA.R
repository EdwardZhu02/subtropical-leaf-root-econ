# ****************************************************************
# One-individual-per-species resampling for individual-level SMA
# Repeatedly sample one individual per species within each community
# and recompute SMA regression indices (P, adjP, R)
# ****************************************************************

rm(list=ls())
library(tidyverse)
library(dplyr)
library(smatr) # SMA regression

load("individual-level-code/traitDataFujian-Ind-step1.RData")

### SET DATA TO USE ###
traitData_touse = traitDataIndv_SelectedTraits_log

### RESAMPLING SETTINGS ###
n_iter = 200
set.seed(1)

trait_pairs = tibble::tribble(
	~pair_id, ~x, ~y, ~label,
	"LNC_LMA", "LMA", "LNC", "LES (LNC-LMA)",
	"RD_SRL", "SRL", "RD", "RES collaboration (RD-SRL)",
	"RNC_RTD", "RTD", "RNC", "RES conservation (RNC-RTD)",
	"Rdark25P_SRR25", "SRR25", "Rdark25P", "Rd25-Rr25",
	"RNC_LNC", "LNC", "RNC", "RNC-LNC",
	"SRL_SLA", "SLA", "SRL", "SRL-SLA"
)

sample_one_ind_per_sp_site <- function(data, trait_cols) {
	data %>%
		dplyr::filter(if_all(all_of(trait_cols), ~ !is.na(.))) %>%
		group_by(SiteID, SpeciesFullName) %>%
		slice_sample(n = 1) %>%
		ungroup()
}

run_sma_stats <- function(data, x, y) {
	data_xy = data %>% dplyr::select(all_of(c(x, y))) %>%
		dplyr::filter(if_all(everything(), ~ !is.na(.)))

	if (nrow(data_xy) < 3) {
		return(tibble(P = NA_real_, R = NA_real_, n = nrow(data_xy)))
	}

	fit = tryCatch(
		smatr::sma(data_xy[[y]] ~ data_xy[[x]]),
		error = function(e) NULL
	)

	if (is.null(fit)) {
		return(tibble(P = NA_real_, R = NA_real_, n = nrow(data_xy)))
	}

	pval = as.numeric(fit[["pval"]])
	r2 = as.numeric(fit[["r2"]])
	slope = as.numeric(fit[["coef"]][[1]][2, 1])
	R = ifelse(is.na(r2) | is.na(slope), NA_real_, ifelse(slope >= 0, sqrt(r2), -sqrt(r2)))

	tibble(P = pval, R = R, n = nrow(data_xy))
}

resampling_stats = vector("list", n_iter)

for (i in seq_len(n_iter)) {
	message("Resampling iteration: ", i, "/", n_iter)

	iter_results = vector("list", nrow(trait_pairs) * 3)
	idx = 1

	for (j in seq_len(nrow(trait_pairs))) {
		pair = trait_pairs[j, ]

		# Resample one individual per species within each site for this trait pair
		sampled = traitData_touse %>%
			sample_one_ind_per_sp_site(trait_cols = c(pair$x, pair$y))

		# Total (all sites combined)
		stats_total = run_sma_stats(sampled, pair$x, pair$y)
		iter_results[[idx]] = tibble(
			iteration = i,
			pair_id = pair$pair_id,
			x = pair$x,
			y = pair$y,
			label = pair$label,
			site = "total"
		) %>% bind_cols(stats_total)
		idx = idx + 1

		# Hilltop
		stats_hilltop = run_sma_stats(dplyr::filter(sampled, SiteID == "hilltop"), pair$x, pair$y)
		iter_results[[idx]] = tibble(
			iteration = i,
			pair_id = pair$pair_id,
			x = pair$x,
			y = pair$y,
			label = pair$label,
			site = "hilltop"
		) %>% bind_cols(stats_hilltop)
		idx = idx + 1

		# Valley
		stats_valley = run_sma_stats(dplyr::filter(sampled, SiteID == "valley"), pair$x, pair$y)
		iter_results[[idx]] = tibble(
			iteration = i,
			pair_id = pair$pair_id,
			x = pair$x,
			y = pair$y,
			label = pair$label,
			site = "valley"
		) %>% bind_cols(stats_valley)
		idx = idx + 1
	}

	resampling_stats[[i]] = dplyr::bind_rows(iter_results)
}

### COMBINE RESULTS ###
resampling_stats_df = dplyr::bind_rows(resampling_stats)

# adjP per iteration, all trait pairs (not used here)
resampling_stats_df = resampling_stats_df %>%
	group_by(iteration) %>%
	mutate(adjP_pair = p.adjust(P, method = "fdr")) %>%
	ungroup()

# adjP per iteration per trait pair
# resampling_stats_df = resampling_stats_df %>%
#   group_by(iteration, pair_id, site) %>%
#   mutate(adjP_pair = p.adjust(P, method = "fdr")) %>%
#   ungroup()

### SUMMARY: MEAN AND 5%/95% CI ###
summary_stats = resampling_stats_df %>%
	group_by(pair_id, x, y, label, site) %>%
	summarise(
		P_mean = mean(P, na.rm = TRUE),
		P_prop_lt_0.05 = mean(P < 0.05, na.rm = TRUE),
		P_prop_lt_0.01 = mean(P < 0.01, na.rm = TRUE),
		adjP_mean = mean(adjP_pair, na.rm = TRUE),
		adjP_prop_lt_0.05 = mean(adjP_pair < 0.05, na.rm = TRUE),
		adjP_prop_lt_0.01 = mean(adjP_pair < 0.01, na.rm = TRUE),
		R_mean = mean(R, na.rm = TRUE),
		R_q05 = quantile(R, 0.05, na.rm = TRUE),
		R_q95 = quantile(R, 0.95, na.rm = TRUE),
		n_mean = mean(n, na.rm = TRUE),
		.n = dplyr::n(),
		.groups = "drop"
	)

### SAVE RESULTS ###
write.csv(
	resampling_stats_df,
	file = "revision-code/step2-biSMA-iter.csv",
	row.names = FALSE
)

write.csv(
	summary_stats,
	file = "revision-code/step2-biSMA-summary.csv",
	row.names = FALSE
)
