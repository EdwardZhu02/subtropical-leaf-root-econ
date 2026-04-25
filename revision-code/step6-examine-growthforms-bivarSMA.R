rm(list=ls())
library(tidyverse)
library(smatr)

# ------------------------------------------------------------------------------
# Examine whether growth forms (tree vs shrub) influence trait coordination
# Uses species means within each site, comparing SMA slopes between groups
# ------------------------------------------------------------------------------

load("individual-level-code/traitDataFujian-Ind-step1.RData")

# Keep only tree and shrub (exclude lianas, n=4 species)
# Collapse the "tree+shrub" species into their primary forms
traitDataGF = traitDataIndv_SelectedTraits %>%
  filter(GrowthForm %in% c("tree", "shrub")) %>%
  mutate(GrowthForm = droplevels(GrowthForm))

# Compute species means within each community
traitName_all = c("LMA", "LNC", "RD", "SRL", "RNC", "RTD", "SRR25", "Rdark25P", "SLA")

# Log transform (signed log, as used throughout the pipeline)
func_log_transform <- function(obj_df) {
  numeric_cols <- sapply(obj_df, is.numeric)
  numeric_data <- obj_df[, numeric_cols]
  log_transformed = sign(numeric_data) * log(abs(numeric_data) + 1)
  obj_df[, numeric_cols] <- log_transformed
  return(obj_df)
}

traitData_touse = func_log_transform(traitDataGF)

# Trait pairs from step5
trait_pairs = tibble::tribble(
  ~pair_id, ~x, ~y, ~label,
  "LNC_LMA", "LMA", "LNC", "LES (LNC-LMA)",
  "RD_SRL", "SRL", "RD", "RES collaboration (RD-SRL)",
  "RNC_RTD", "RTD", "RNC", "RES conservation (RNC-RTD)",
  "Rdark25P_SRR25", "SRR25", "Rdark25P", "Rd25-Rr25",
  "RNC_LNC", "LNC", "RNC", "RNC-LNC",
  "SRL_LMA", "LMA", "SRL", "SRL-LMA"
)

# Helper: fit SMA and return slope, R2, p-value
fit_sma <- function(data, x_var, y_var) {
  data = data %>% filter(!is.na(.data[[x_var]]) & !is.na(.data[[y_var]]))
  if (nrow(data) <= 2) return(NULL)
  tryCatch(
    sma(data[[y_var]] ~ data[[x_var]]),
    error = function(e) NULL
  )
}

# Helper: extract stats from a model
model_stats <- function(model, group_name) {
  if (is.null(model)) {
    return(tibble(
      group = group_name, n = NA_real_, slope = NA_real_,
      intercept = NA_real_, r2 = NA_real_, pval = NA_real_,
      R_signed = NA_real_
    ))
  }
  slope = as.numeric(model$coef[[1]][2, 1])
  intercept = as.numeric(model$coef[[1]][1, 1])
  r2 = as.numeric(model$r2)
  p = as.numeric(model$pval)
  n = as.numeric(model$n)
  R_signed = ifelse(slope > 0, sqrt(r2), -sqrt(r2))
  tibble(
    group = group_name, n = n, slope = slope,
    intercept = intercept, r2 = r2, pval = p, R_signed = R_signed
  )
}

# Format p-value for annotation
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("<0.001")
  if (p < 0.01) return("<0.01")
  if (p < 0.05) return("<0.05")
  return(paste0("=", signif(p, 2)))
}

# Annotation builder
annot_text <- function(model, group_label) {
  if (is.null(model)) return(paste0(group_label, ": NA"))
  slope = as.numeric(model$coef[[1]][2, 1])
  r2 = as.numeric(model$r2)
  p = as.numeric(model$pval)
  R_signed = ifelse(slope > 0, sqrt(r2), -sqrt(r2))
  paste0(group_label, ": slope=", round(slope, 2), " R=", round(R_signed, 2), " P", fmt_p(p))
}

# Pooled test - does growth form influence covariation?
# Uses only tree vs shrub (lianas excluded above)

extract_common_slope_test <- function(model_obj) {
  out = list(chi_sq = NA_real_, df = NA_real_, pval = NA_real_)
  if (is.null(model_obj) || is.null(model_obj$commoncoef)) return(out)

  cc = model_obj$commoncoef
  if (is.list(cc)) {
    if (!is.null(cc$LR)) out$chi_sq = as.numeric(cc$LR[1])
    if (!is.null(cc$df)) out$df = as.numeric(cc$df[1])
    if (!is.null(cc$p)) out$pval = as.numeric(cc$p[1])
  }

  out
}

extract_single_slope <- function(model_obj) {
  if (is.null(model_obj) || is.null(model_obj$coef)) return(NA_real_)
  as.numeric(model_obj$coef[[1]][2, 1])
}

extract_corr_stats <- function(data, x_var, y_var) {
  data = data %>% filter(!is.na(.data[[x_var]]) & !is.na(.data[[y_var]]))
  n = nrow(data)
  if (n <= 3) {
    return(list(n = n, r = NA_real_, p = NA_real_))
  }
  ct = tryCatch(
    cor.test(data[[x_var]], data[[y_var]], method = "pearson"),
    error = function(e) NULL
  )
  if (is.null(ct)) return(list(n = n, r = NA_real_, p = NA_real_))
  list(n = n, r = as.numeric(ct$estimate), p = as.numeric(ct$p.value))
}

compare_corr_strength <- function(r1, n1, r2, n2) {
  out = list(z = NA_real_, p = NA_real_)
  if (is.na(r1) || is.na(r2) || is.na(n1) || is.na(n2) || n1 <= 3 || n2 <= 3) return(out)

  # Clamp away from +/-1 to avoid infinite Fisher transform in edge cases.
  r1 = max(min(r1, 0.999999), -0.999999)
  r2 = max(min(r2, 0.999999), -0.999999)

  z_val = (atanh(r1) - atanh(r2)) / sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
  p_val = 2 * (1 - pnorm(abs(z_val)))
  out$z = z_val
  out$p = p_val
  out
}

# Keep pooled coordination metrics as context, but test slope differences for all pairs
combined_results = list()

for (i in seq_len(nrow(trait_pairs))) {
  pair = trait_pairs[i, ]
  x_var = pair$x
  y_var = pair$y

  data_sub = traitData_touse %>%
    filter(GrowthForm %in% c("tree", "shrub")) %>%
    filter(!is.na(.data[[x_var]]) & !is.na(.data[[y_var]])) %>%
    mutate(GrowthForm = factor(GrowthForm, levels = c("tree", "shrub")))

  n_tree = sum(data_sub$GrowthForm == "tree", na.rm = TRUE)
  n_shrub = sum(data_sub$GrowthForm == "shrub", na.rm = TRUE)

  model_all = fit_sma(data_sub, x_var, y_var)
  model_tree = fit_sma(data_sub %>% filter(GrowthForm == "tree"), x_var, y_var)
  model_shrub = fit_sma(data_sub %>% filter(GrowthForm == "shrub"), x_var, y_var)

  overall_r2 = if (is.null(model_all)) NA_real_ else as.numeric(model_all$r2)
  overall_p = if (is.null(model_all)) NA_real_ else as.numeric(model_all$pval)
  slope_tree = extract_single_slope(model_tree)
  slope_shrub = extract_single_slope(model_shrub)
  slope_diff = slope_tree - slope_shrub
  sign_change = !is.na(slope_tree) && !is.na(slope_shrub) && (slope_tree * slope_shrub < 0)

  corr_tree_stats = extract_corr_stats(data_sub %>% filter(GrowthForm == "tree"), x_var, y_var)
  corr_shrub_stats = extract_corr_stats(data_sub %>% filter(GrowthForm == "shrub"), x_var, y_var)
  corr_comp = compare_corr_strength(
    corr_tree_stats$r, corr_tree_stats$n,
    corr_shrub_stats$r, corr_shrub_stats$n
  )

  p_slope = NA_real_
  chi_sq = NA_real_
  df_test = NA_real_

  if (n_tree > 2 && n_shrub > 2) {
    formula_grp = as.formula(paste0(y_var, " ~ ", x_var, " * GrowthForm"))
    model_multi = tryCatch(
      sma(formula_grp, data = data_sub),
      error = function(e) NULL
    )

    test_vals = extract_common_slope_test(model_multi)
    chi_sq = test_vals$chi_sq
    df_test = test_vals$df
    p_slope = test_vals$pval
  }

  combined_results[[length(combined_results) + 1]] = tibble(
    pair_id = pair$pair_id,
    trait_pair = pair$label,
    n_tree = n_tree,
    n_shrub = n_shrub,
    overall_r2 = overall_r2,
    overall_p = overall_p,
    slope_tree = slope_tree,
    slope_shrub = slope_shrub,
    slope_difference = slope_diff,
    sign_change = sign_change,
    slope_test_chi_sq = chi_sq,
    slope_test_df = df_test,
    slope_test_p = p_slope,
    corr_tree = corr_tree_stats$r,
    corr_shrub = corr_shrub_stats$r,
    corr_tree_p = corr_tree_stats$p,
    corr_shrub_p = corr_shrub_stats$p,
    corr_strength_tree = abs(corr_tree_stats$r),
    corr_strength_shrub = abs(corr_shrub_stats$r),
    corr_strength_difference = abs(corr_tree_stats$r) - abs(corr_shrub_stats$r),
    corr_diff_z = corr_comp$z,
    corr_diff_p = corr_comp$p
  )
}

growthform_effect_table_full = bind_rows(combined_results) %>%
  mutate(
    slope_test_p_fdr = p.adjust(slope_test_p, method = "fdr"),
    corr_diff_p_fdr = p.adjust(corr_diff_p, method = "fdr"),
    influenced_by_growth_form = case_when(
      !is.na(slope_test_p_fdr) & slope_test_p_fdr < 0.05 ~ "Yes",
      !is.na(slope_test_p_fdr) ~ "No",
      TRUE ~ "Not assessed (insufficient sample size)"
    ),
    corr_strength_differs_by_growth_form = case_when(
      !is.na(corr_diff_p_fdr) & corr_diff_p_fdr < 0.05 ~ "Yes",
      !is.na(corr_diff_p_fdr) ~ "No",
      TRUE ~ "Not assessed (insufficient sample size)"
    )
  )

growthform_effect_table = growthform_effect_table_full %>%
  transmute(
    pair_id,
    trait_pair,
    n_tree,
    n_shrub,
    slope_tree,
    slope_shrub,
    slope_test_df,
    slope_test_p,
    slope_test_p_fdr,
    influenced_by_growth_form,
    corr_tree,
    corr_shrub,
    corr_diff_p,
    corr_diff_p_fdr,
    corr_strength_differs_by_growth_form,
    p_value_basis_slope = "SMA LR slope-equality test compares signed slopes (not absolute values)",
    p_value_basis_correlation = "Fisher r-to-z test compares signed Pearson r (absolute R is descriptive only)"
  )

growthform_effect_table
write_csv(growthform_effect_table, "revision-code/step6-slope-comparison-GF-combined.csv")