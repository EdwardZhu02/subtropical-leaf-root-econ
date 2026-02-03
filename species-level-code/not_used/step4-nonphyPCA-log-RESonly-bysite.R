rm(list=ls())
library(tidyverse)
library(dplyr)
library(factoextra) # perform PCA analysis
library(FactoMineR) # perform PCA analysis
library(patchwork) # plot merging
library(cowplot) # plot merging

# ------------------------------------------------------------------------------
# Perform PCA analysis (RES traits only)
# Species means within each community (Hilltop/Valley)
# ------------------------------------------------------------------------------

# Load individual level data to compute site-specific species means
load("individual-level-code/traitDataFujian-Ind-step1.RData")

# traitDataIndv_SelectedTraits contains the raw data (with some 1/SLA conversions etc)
# Select traits
traitName_res = c("RTD","SRL","RD","RNC")
traitName_all = c("RTD","SRL","RD","RNC","LMA","LNC","RDMC","SRR25","SRA","LPC")

# Compute species means within each community (SiteID)
traitData_bysite_mean = traitDataIndv_SelectedTraits %>%
  group_by(SpeciesFullName, SiteID) %>%
  summarise(
    GrowthForm = first(GrowthForm), # Assume constant within species-site
    across(any_of(traitName_all), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )

# Log transform
func_log_transform <- function(obj_df) {
  numeric_cols <- sapply(obj_df, is.numeric)
  numeric_data <- obj_df[, numeric_cols]
  # Log-transform each trait (signed log transformation, dealing with negative values)
  log_transformed = sign(numeric_data) * log(abs(numeric_data) + 1)
  obj_df[, numeric_cols] <- log_transformed
  return(obj_df)        
}

traitData_bysite_log = func_log_transform(traitData_bysite_mean)

### SET DATA TO USE ###
traitDataPCA_touse = traitData_bysite_log
traitName_touse = traitName_res

# Create identifier for PCA
traitDataPCA_touse = traitDataPCA_touse %>% mutate(
  PCAIdentifier = paste(gsub(" ","_",SpeciesFullName), SiteID, sep="-")
)

# Curate data for PCA analysis
nonphyloPCAData_numonly = traitDataPCA_touse %>% ungroup() %>% 
  dplyr::select(PCAIdentifier, SiteID, all_of(traitName_touse)) %>% 
  na.omit() # Remove rows with NA in used traits

# Scale traits (normalization)
nonphyloPCAData_numonly[, traitName_touse] = apply(nonphyloPCAData_numonly[, traitName_touse], 2, scale)

# Set rownames
rownames(nonphyloPCAData_numonly) = nonphyloPCAData_numonly$PCAIdentifier

# used for annotation
nonphyloPCAData_meta = nonphyloPCAData_numonly
nonphyloPCAData_matrix = nonphyloPCAData_numonly %>% dplyr::select(all_of(traitName_touse))

# perform PCA
nonphyloPCAresult = PCA(nonphyloPCAData_matrix, scale.unit = T, graph = F) 

# ------------------------------------------------------------------------------
# Scree plot
# ------------------------------------------------------------------------------
nonphyloPCA.screeplot = fviz_screeplot(nonphyloPCAresult, addlabels = TRUE, ylim = c(0, 100),
                                       barfill="#7998AD", barcolor="#7998AD", linecolor="black") +
  labs(title="PCA, RES only (Site-Species Mean)") + theme_classic()
ggsave(filename = "species-level-code/plt_nonphyPCA_scree-RESonly-bysite.pdf", plot = nonphyloPCA.screeplot, width = 3, height = 3)

# ------------------------------------------------------------------------------
# 2-D PCA biplot: SiteID as grouping factor
# ------------------------------------------------------------------------------

# Define palette for sites (Hilltop, Valley)
# Hilltop="#547bb4", Valley="#dd7c4f" (observed in individual level analysis)
# Levels are "hilltop", "valley" (alphabetical)
site_colors = c("#547bb4", "#dd7c4f")

plt_nonphyloPCA_biplot_ax12 = fviz_pca_biplot(
  nonphyloPCAresult, label = "var", habillage=nonphyloPCAData_meta$SiteID,
  axes=c(1,2),
  repel=T,
  col.var = "gray20",
  palette = site_colors,
  addEllipses=T,ellipse.level = 0.95) +
  labs(title = "") + theme_classic() + theme(legend.direction = 'horizontal', legend.position = 'top', legend.title = element_blank())

plt_nonphyloPCA_contribplot_ax12 = fviz_pca_var(
  nonphyloPCAresult, col.var = "contrib",
  axes=c(1,2),
  repel=T
) +
  theme_classic() +
  scale_color_gradientn(colors = c("#71BFB2", "#F2BA2F", "#AD0B08")) +
  labs(title = "")

# Save plot
plt_ax12 = plt_nonphyloPCA_biplot_ax12 + plt_nonphyloPCA_contribplot_ax12 + 
  plot_layout(guides = "collect") & theme(legend.position='top')

ggsave(plot = plt_ax12, filename = "species-level-code/plt_nonphyPCA_ax12Combined_RESonly_bysite.pdf",
       width = 6.5, height = 4)

