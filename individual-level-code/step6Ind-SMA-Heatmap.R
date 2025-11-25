rm(list=ls())
library(tidyverse)
library(dplyr)
library(smatr) # SMA regression
library(ggthemes)
library(pheatmap)
library(reshape2) # for drawing heatmap using ggplot2

load("individual-level-code/traitDataFujian-Ind-step1.RData") # TODO: individual-level data

# == Already done in preparation of the data, updated 7/11-24
# # ------------------------------------------------------------------------------
# # Prepare only log-transformed data
# # function for data normalization (log-transform only)
# func_log_transform <- function(obj_df) {
#   numeric_cols <- sapply(obj_df, is.numeric)
#   numeric_data <- obj_df[, numeric_cols]
#   # Log-transform each trait (signed log transformation, dealing with negative values)
#   log_transformed = sign(numeric_data) * log(abs(numeric_data) + 1)
#   obj_df[, numeric_cols] <- log_transformed
#   return(obj_df)        
# }
# traitDataIndv_spavg_log = func_log_transform(traitDataIndv_spavg)
# # ------------------------------------------------------------------------------

### SET DATA TO USE ###
traitData_touse = traitDataIndv_SelectedTraits_log

# TODO: (25.2.21) Resolve historical problem: nomination between "speciesFullName" and "SpeciesFullName"
traitData_touse = traitData_touse %>% rename(speciesFullName=SpeciesFullName)

# divide data based on growth forms for separate regression
traitData_touse_tree = traitData_touse %>% dplyr::filter(GrowthForm == "tree")
traitData_touse_shrub = traitData_touse %>% dplyr::filter(GrowthForm == "shrub")
traitData_touse_liana = traitData_touse %>% dplyr::filter(GrowthForm == "liana") # may be useless

# ------------------------------------------------------------------------------
# Batch SMA implementation to test variable correlation
# ------------------------------------------------------------------------------
# Define the list of variables
variables = c("RTD","SRL","RD","RNC","RDMC","SRR25","SRA","RPC","RCC", # Fine root traits
              "LMA","LNC","LPC","LCC","Ld13C","Rdark25P","Vcmax25","Asat","LA" # Leaf traits
)

# Initialize an empty matrix to store the P values and R square values
pvalue_matrix = matrix(NA, nrow = length(variables), ncol = length(variables),
                        dimnames = list(variables, variables))
rvalue_matrix = matrix(NA, nrow = length(variables), ncol = length(variables), 
                        dimnames = list(variables, variables))

# Loop through each pair of variables
for (i in 1:length(variables)) {
  for (j in 1:length(variables)) {
    if (i != j) {
      # Fit the SMA model for the variable pair
      formula <- as.formula(paste0("traitData_touse$", variables[i], " ~ ", "traitData_touse$", variables[j]))
      model.sma <- tryCatch({ sma(formula) }, error = function(e) return(NULL)) # Error handling
      
      if (!is.null(model.sma)) {
        p_value <- as.numeric(model.sma[["pval"]])
        r_value <- ifelse(
          as.numeric(model.sma[["groupsummary"]][["Slope"]]) > 0,
          sqrt(as.numeric(model.sma[["r2"]])),
          -sqrt(as.numeric(model.sma[["r2"]]))
        )
        
        pvalue_matrix[i,j] <- p_value
        rvalue_matrix[i,j] <- r_value
      }
    }
    #break # for testing
  }
}

# Plot the heatmap (method 1: ggplot)
get_upper_tri = function(cormat){ # extracting only the lower part
  cormat[lower.tri(cormat)] = NA
  return(cormat)
}

pvalue_matrix_uptri = get_upper_tri(pvalue_matrix)
pvalue_matrix_uptri = round(pvalue_matrix_uptri,2)
pvalue_matrix_uptri = ifelse(pvalue_matrix_uptri>0.1,0.1,pvalue_matrix_uptri) # filter later, just for coloring
melted_pvalue_matrix = melt(pvalue_matrix_uptri, na.rm = TRUE)

rvalue_matrix_uptri = get_upper_tri(rvalue_matrix)
rvalue_matrix_uptri = round(rvalue_matrix_uptri,2)
melted_rvalue_matrix = melt(rvalue_matrix_uptri, na.rm = TRUE)


plt_SMA_htmap_pvalue = ggplot(data = melted_pvalue_matrix, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#CE403F", high = "white",
                       midpoint = 0.085, limit = c(0,0.1), space = "Lab",
                       name="SMA - P value\nCutoff: P > 0.1") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11), axis.text.y = element_text(size = 10)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = ifelse(value<0.1, value, " ")), color = "black", size = 2.8) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))
ggsave(plot=plt_SMA_htmap_pvalue, filename="individual-level-code/Ind_SMAhtmap_Pvalue.pdf", width=8, height=4.5)


plt_SMA_htmap_rvalue = ggplot(data = melted_rvalue_matrix, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#08519C", mid = "white", high = "#CB181C",
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="SMA - R value\nCutoff: abs(R) < 0.25") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11), axis.text.y = element_text(size = 10)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = ifelse(abs(value)>0.25, value, " ")), color = "black", size = 2.7) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))
ggsave(plot=plt_SMA_htmap_rvalue, filename="individual-level-code/Ind_SMAhtmap_Rvalue.pdf", width=8, height=4.5)
