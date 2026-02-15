rm(list=ls())
library(tidyverse)
library(dplyr)
library(smatr) # SMA regression
library(ggthemes)
library(patchwork) # plot merging
library(cowplot) # plot merging

# ------------------------------------------------------------------------------
# Species-level SMA by site
# Using species means within each community
# ------------------------------------------------------------------------------

load("individual-level-code/traitDataFujian-Ind-step1.RData")

# Compute species means within each community (SiteID)
# Include all traits needed
traitName_all = c("LMA", "LNC", "RD", "SRL", "RNC", "RTD", "SRR25", "Rdark25P", "SLA")

traitData_bysite_mean = traitDataIndv_SelectedTraits %>%
  group_by(SpeciesFullName, SiteID) %>%
  summarise(
    GrowthForm = first(GrowthForm),
    across(any_of(traitName_all), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Log transform
func_log_transform <- function(obj_df) {
  numeric_cols <- sapply(obj_df, is.numeric)
  numeric_data <- obj_df[, numeric_cols]
  log_transformed = sign(numeric_data) * log(abs(numeric_data) + 1)
  obj_df[, numeric_cols] <- log_transformed
  return(obj_df)        
}

traitData_touse = func_log_transform(traitData_bysite_mean)

# Define trait pairs (from revision-code/step2-1indvresamp-biSMA.R)
trait_pairs = tibble::tribble(
  ~pair_id, ~x, ~y, ~label,
  "LNC_LMA", "LMA", "LNC", "LES (LNC-LMA)",
  "RD_SRL", "SRL", "RD", "RES collaboration (RD-SRL)",
  "RNC_RTD", "RTD", "RNC", "RES conservation (RNC-RTD)",
  "Rdark25P_SRR25", "SRR25", "Rdark25P", "Rd25-Rr25",
  "RNC_LNC", "LNC", "RNC", "RNC-LNC",
  "SRL_SLA", "SLA", "SRL", "SRL-SLA"
)

# Function to perform SMA and Plot
run_and_plot_sma <- function(data, x_var, y_var, label_text) {
  
  # Ensure columns are numeric
  data[[x_var]] = as.numeric(data[[x_var]])
  data[[y_var]] = as.numeric(data[[y_var]])
  
  # Remove NA
  data_clean = data %>% filter(!is.na(.data[[x_var]]) & !is.na(.data[[y_var]]))
  data_hilltop = data_clean %>% filter(SiteID == "hilltop")
  data_valley = data_clean %>% filter(SiteID == "valley")
  
  # SMA Models
  valid_total = nrow(data_clean) > 2
  valid_hilltop = nrow(data_hilltop) > 2
  valid_valley = nrow(data_valley) > 2
  
  sma_total = if(valid_total) sma(data_clean[[y_var]] ~ data_clean[[x_var]]) else NULL
  sma_hilltop = if(valid_hilltop) sma(data_hilltop[[y_var]] ~ data_hilltop[[x_var]]) else NULL
  sma_valley = if(valid_valley) sma(data_valley[[y_var]] ~ data_valley[[x_var]]) else NULL
  
  # Extract stats for annotation
  get_stats <- function(model, name) {
    if(is.null(model)) return(paste0(name, ": NA"))
    r2 = as.numeric(model$r2)
    p = as.numeric(model$pval)
    slope = as.numeric(model$coef[[1]][2,1])
    R_val = ifelse(slope>0, sqrt(r2), -sqrt(r2))
    
    p_str = ifelse(p < 0.01, "<0.01", ifelse(p < 0.05, "<0.05", paste0("=", signif(p, 2))))
    return(paste0(name, ": R=", round(R_val, 2), " P", p_str))
  }
  
  annot_total = get_stats(sma_total, "M")
  annot_hilltop = get_stats(sma_hilltop, "H")
  annot_valley = get_stats(sma_valley, "V")
  
  full_annot = paste(annot_total, annot_hilltop, annot_valley, sep="\n")
  
  # Prediction for plotting lines
  # Helper to add line
  add_sma_line <- function(plt, model, color, linetype="solid") {
    if(is.null(model)) return(plt)
    
    intercept = model$coef[[1]][1,1]
    slope = model$coef[[1]][2,1]
    
    # We can use geom_abline, but it extends to infinity. 
    # Better to predict on range, but for simplicity abline is often used in SMA plots or segment
    plt + geom_abline(intercept = intercept, slope = slope, color = color, linetype = linetype, linewidth = 1)
  }
  
  p = ggplot(data_clean, aes(x = .data[[x_var]], y = .data[[y_var]])) +
    geom_point(aes(color = SiteID, shape = GrowthForm), size = 2, alpha=0.95) +
    scale_color_manual(values = c("hilltop" = "#547bb4", "valley" = "#dd7c4f")) +
    theme_classic() +
    theme(legend.position = "none") +
    labs(title = label_text, x = x_var, y = y_var) +
    annotate("text", x = -Inf, y = Inf, label = full_annot, hjust = -0.1, vjust = 1.1, size = 3)
  
  # Add lines
  p = add_sma_line(p, sma_total, "black", "dashed")
  p = add_sma_line(p, sma_hilltop, "#547bb4", "solid")
  p = add_sma_line(p, sma_valley, "#dd7c4f", "solid")
  
  return(p)
}

# Pearson correlation stats for each site
get_cor_stats <- function(data, x_var, y_var, site_label, site_filter = NULL) {
  data_sub = data
  if(!is.null(site_filter)) {
    data_sub = data_sub %>% filter(SiteID == site_filter)
  }
  data_sub = data_sub %>% filter(!is.na(.data[[x_var]]) & !is.na(.data[[y_var]]))
  n_obs = nrow(data_sub)
  if(n_obs > 2) {
    ct = cor.test(data_sub[[x_var]], data_sub[[y_var]], method = "pearson")
    r_val = as.numeric(ct$estimate)
    p_val = as.numeric(ct$p.value)
  } else {
    r_val = NA_real_
    p_val = NA_real_
  }
  tibble(
    site = site_label,
    n = n_obs,
    pearson_r = r_val,
    p_value = p_val
  )
}

# Run loop
plot_list = list()
stats_list = list()

for(i in 1:nrow(trait_pairs)) {
  pair = trait_pairs[i,]
  p = run_and_plot_sma(traitData_touse, pair$x, pair$y, pair$label)
  plot_list[[i]] = p

  stats_list[[i]] = bind_rows(
    get_cor_stats(traitData_touse, pair$x, pair$y, "total"),
    get_cor_stats(traitData_touse, pair$x, pair$y, "hilltop", "hilltop"),
    get_cor_stats(traitData_touse, pair$x, pair$y, "valley", "valley")
  ) %>%
    mutate(
      pair_id = pair$pair_id,
      x = pair$x,
      y = pair$y,
      label = pair$label
    )
  
  # Save individual plots
  ggsave(plot = p, filename = paste0("species-level-code/smaplt_bysite_", pair$pair_id, ".pdf"), width = 2.8, height = 2.8)
}

cor_stats = bind_rows(stats_list) %>%
  group_by(site) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

write_csv(cor_stats, "species-level-code/step5-sma_spmean_bysite_stats.csv")

# ------------------------------------------------------------------------------
# Individual-level SMA by site
# ------------------------------------------------------------------------------

traitData_indv_touse = func_log_transform(traitDataIndv_SelectedTraits)

indv_plot_list = list()
indv_stats_list = list()

for(i in 1:nrow(trait_pairs)) {
  pair = trait_pairs[i,]
  p_indv = run_and_plot_sma(traitData_indv_touse, pair$x, pair$y, paste0(pair$label, " (indv)"))
  indv_plot_list[[i]] = p_indv

  indv_stats_list[[i]] = bind_rows(
    get_cor_stats(traitData_indv_touse, pair$x, pair$y, "total"),
    get_cor_stats(traitData_indv_touse, pair$x, pair$y, "hilltop", "hilltop"),
    get_cor_stats(traitData_indv_touse, pair$x, pair$y, "valley", "valley")
  ) %>%
    mutate(
      pair_id = pair$pair_id,
      x = pair$x,
      y = pair$y,
      label = pair$label
    )

  ggsave(plot = p_indv, filename = paste0("individual-level-code/smaplt_indv_bysite_", pair$pair_id, ".pdf"), width = 2.8, height = 2.8)
}

cor_stats_indv = bind_rows(indv_stats_list) %>%
  group_by(site) %>%
  mutate(p_adjusted = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

write_csv(cor_stats_indv, "individual-level-code/step5-sma_indv_bysite_stats.csv")

# Combine all into one
# combined_plot = wrap_plots(plot_list, ncol = 3)
# ggsave(plot = combined_plot, filename = "species-level-code/smaplt_bysite_AllPairs.pdf", width = 10, height = 7)
