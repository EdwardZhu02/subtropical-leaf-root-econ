rm(list=ls())
library(tidyverse)
library(dplyr)
library(tidyr)
library(ape) # general phylogenetic analysis
library(caper) # PGLS
library(geiger) # compare taxa in data and tree
# remotes::install_github("uyedaj/treeplyr")
library(treeplyr) # general phylogenetic analysis
library(phytools) # general phylogenetic analysis
library(ggtree) # phylogenetic tree visualization
library(aplot) # combine subplots
library(ggheatmap) #devtools::install_github("XiaoLuo-boy/ggheatmap")
library(pheatmap)
library(reshape2) # for drawing heatmap using ggplot2
# library(rr2) # calculate r2 for PGLS models


# TODO: Updated 21/11-24, use original and perform binding in the latter script
# Load: traitDataIndv, traitDataIndv_spavg, traitDataIndv_spavg_log, traitDataIndv_spavg_log_ztransform, spTaxaNames_withfam, phylotree_result
load("traitDataFujian-spavg-phylo-step2.RData") 

# Refs (different approaches for conducting the lambda-model in phylo regression):
# http://blog.phytools.org/2012/11/fitting-model-in-phylogenetic.html

# Part 1: combine raw data and tree -------------------------------------
# 
### SET DATA TO USE ###
traitData_touse = traitDataIndv_spavg_log

# Convert species Latin name `A B` to `A_B` to match with tree `tip_label`
traitData_touse = traitData_touse %>% dplyr::mutate(tipLabelMatched=gsub(" ", "_", speciesFullName)) %>%
  dplyr::select(tipLabelMatched, everything())

# Extract phylo-tree constructed in step2
# phylotree_result[[1]] is the resulting phylogenetic tree (type: phylo).
phylotree_touse = phylotree_result[[1]] 
rm(phylotree_result)

# Combine tree and data
traitData_touse_PhyloDataCombined = make.treedata(
  tree = phylotree_touse, data = traitData_touse, name_column = "tipLabelMatched")

# Save data and tree separately, facilitating downstream analyses.
traitData_touse_PhyloDataCombined_data = as.data.frame(traitData_touse_PhyloDataCombined$dat) %>% 
  dplyr::mutate(tipLabelMatched=gsub(" ", "_", speciesFullName)) %>%
  dplyr::select(tipLabelMatched, everything())
traitData_touse_PhyloDataCombined_tree = traitData_touse_PhyloDataCombined$phy

# set tree node label to NULL to avoid Error: Labels duplicated between tips and nodes in phylogeny
# Ref: https://stackoverflow.com/questions/51261388/duplicate-tips-label-removal-importing-phylogenetic-tree-in-r-for-comparison
traitData_touse_PhyloDataCombined_tree$node.label = NULL


# divide data based on growth forms for separate regression
traitData_touse_tree = traitData_touse %>% dplyr::filter(GrowthForm == "tree")
traitData_touse_shrub = traitData_touse %>% dplyr::filter(GrowthForm == "shrub")
traitData_touse_liana = traitData_touse %>% dplyr::filter(GrowthForm == "liana") # may be useless


# Step2: Batch PGLS implementation to test variable correlation ----------------
#
# Define the list of variables
variables = c("RTD","SRL","RD","RNC","RDMC","SRR25","SRA","RPC","RCC", # Fine root traits
              "LMA","LNC","LPC","LCC","Ld13C","Rdark25P","Vcmax25","Asat","LA" # Leaf traits
)

# Initialize an empty matrix to store the P values and R square values
pvalue_matrix = matrix(NA, nrow = length(variables), ncol = length(variables),
                       dimnames = list(variables, variables))
rsq_matrix = matrix(NA, nrow = length(variables), ncol = length(variables), 
                    dimnames = list(variables, variables))
lambda_matrix = matrix(NA, nrow = length(variables), ncol = length(variables), 
                    dimnames = list(variables, variables))

PGLSdata = comparative.data(phy = traitData_touse_PhyloDataCombined_tree,
                            data = traitData_touse_PhyloDataCombined_data,
                            names.col = tipLabelMatched, vcv = T, na.omit = F, warn.dropped = T)

# Loop through each pair of variables ------------------------------------------
for (i in 1:length(variables)) { # row
  for (j in 1:length(variables)) { # column
    if (i != j) {
      # Fit the SMA model for the variable pair
      formula <- as.formula(paste(variables[i], "~", variables[j])) # i-Y (dependent), j-X (independent)
      model.pgls <- tryCatch({ pgls(formula, data = PGLSdata, lambda = "ML") },
                             error = function(e) return(NULL)) # Error handling
      
      if (!is.null(model.pgls)) {
        anova_result <- anova.pgls(model.pgls)
        p_value <- anova_result[1, 5] # Extract ANOVA Pr(>F)
        model_summary <- summary(model.pgls)
        r_value <- ifelse(
          as.numeric(model_summary[["adj.r.squared"]]) > 0,
          ifelse(
            as.numeric(model.pgls[["model"]][["coef"]][2]) > 0,
            sqrt(as.numeric(model_summary[["adj.r.squared"]])),
            -sqrt(as.numeric(model_summary[["adj.r.squared"]]))
          ),
          0 # if r2 < 0 in PGLS
        )
        lambda_value <- as.numeric(model_summary[["param"]][["lambda"]]) # Extract lambda
        # lambda estimated via max likelihood method
        
        pvalue_matrix[i, j] <- p_value
        rsq_matrix[i,j] <- r_value
        lambda_matrix[i,j] <- lambda_value
      } else {
        # Error in pgls(as.formula("SRR25 ~ RDMC"), data = PGLSdata, lambda = "ML"): 
        # Problem with optim:52ERROR: ABNORMAL_TERMINATION_IN_LNSRCH
        # 
        # Solution (12/11-24): give all lambda=1 for matches failed in lambda search
        model.pgls <- pgls(formula, data = PGLSdata, lambda = 1)
        anova_result <- anova.pgls(model.pgls)
        p_value <- anova_result[1, 5] # Extract ANOVA Pr(>F)
        model_summary <- summary(model.pgls)
        r_value <- ifelse(
          as.numeric(model_summary[["adj.r.squared"]]) > 0,
          ifelse(
            as.numeric(model.pgls[["model"]][["coef"]][2]) > 0,
            sqrt(as.numeric(model_summary[["adj.r.squared"]])),
            -sqrt(as.numeric(model_summary[["adj.r.squared"]]))
          ),
          0 # if r2 < 0 in PGLS
        )
        lambda_value <- as.numeric(model_summary[["param"]][["lambda"]]) # Extract lambda
        # lambda estimated via max likelihood method
        
        pvalue_matrix[i, j] <- p_value
        rsq_matrix[i,j] <- r_value
        lambda_matrix[i,j] <- lambda_value
      }
    }
  }
}

# # Extra part: examining matches that dont have a fit ---------------------------
# model.testing = pgls(as.formula("SRR25 ~ LNC"), data = PGLSdata, lambda = "ML")
# model_summary = summary(model.testing)
# # plot lambda optimization map for SRR25.RDMC model
# # Create a likelihood profile of the lambda estimate
# model_lambda.profile = pgls.profile(model.testing, "lambda")
# # Plot the likelihood profile
# plot(model_lambda.profile)
# # Conclusion: give all lambda=1 for matches failed in lambda search
# # ---------------------------

# Plot the heatmap (method: ggplot)
# PGLS is about response (Y) - predictor (X) relationship, so both half of the
# heatmap is meaningful. Display the whole plot for reference.

# convert non-significant P to 0.1, making color gradients more readily displayed
pvalue_matrix = round(pvalue_matrix,2)
pvalue_matrix = ifelse(pvalue_matrix>0.1,0.1,pvalue_matrix)
melted_pvalue_matrix = melt(pvalue_matrix, na.rm = F)

rsq_matrix = round(rsq_matrix,2)
# rsq_matrix = ifelse(rsq_matrix<0,0,rsq_matrix)
melted_rsq_matrix = melt(rsq_matrix, na.rm = F)

lambda_matrix = round(lambda_matrix,2)
#lambda_matrix = ifelse(lambda_matrix<0,0,lambda_matrix)
melted_lambda_matrix = melt(lambda_matrix, na.rm = F)

# RLES_vars = c("Asat","LCC","LMA","LNC","LPC","Rdark25P","Vcmax25","RCC","RD","RDMC","RNC","RPC","RTD","SRA","SRL","SRR25")
# melted_lambda_matrix = melted_lambda_matrix %>% dplyr::filter(!(Var1 %in% RLES_vars))

# The lambda plot --------------------------------------------------------------
plt_htmap_lambda = ggplot(data = melted_lambda_matrix, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "white", high = "#719aac",
                       na.value = "grey50",
                       midpoint = 0.001, limit = c(0,1), space = "Lab", 
                       name="lambda") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = ifelse(value>0.1, value, " ")), color = "black", size = 2.6) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    #legend.justification = c(1, 0),
    #legend.position = c(0.6, 0.7),
    legend.position = "right",
    legend.direction = "vertical")
#guides(fill = guide_colorbar(barwidth = 1, barheight = 7, title.position = "top", title.hjust = 0.5))
#ggsave(plot=plt_htmap_lambda, filename="plt_PGLS_htmap_lambda.pdf", width=9.5, height=10)
ggsave(plot=plt_htmap_lambda, filename="plt_PGLS_htmap_lambda_all.pdf", width=5.3, height=4.3)
# ggsave(plot=plt_htmap_lambda, filename="plt_PGLS_htmap_lambda_all.png", width=7, height=5.5, dpi=300)

# Filter only significant dependent variables (Y) that have a lambda > 0 across almost all independent variables (X) ------------------------------------------------------------------
#signif_dependent_vars = c("rs25", "LCC", "LMA", "SRR25", "RPC", "RD", "SRL")
# melted_pvalue_matrix_lambdaSignif = melted_pvalue_matrix %>% 
#   dplyr::filter(Var1 %in% signif_dependent_vars)
# melted_rsq_matrix_lambdaSignif = melted_rsq_matrix %>% 
#   dplyr::filter(Var1 %in% signif_dependent_vars)


signif_dependent_vars = c("SRR25", "SRL", "RPC", "RD", "LMA", "LCC")
melted_pvalue_matrix_lambdaSignif = melted_pvalue_matrix %>%
  dplyr::filter(Var1 %in% signif_dependent_vars)
melted_rsq_matrix_lambdaSignif = melted_rsq_matrix %>%
  dplyr::filter(Var1 %in% signif_dependent_vars)

# The ANOVA p value plot -------------------------------------------------------
plt_htmap_pvalue = ggplot(data = melted_pvalue_matrix_lambdaSignif, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#CE403F", high = "white",
                       na.value = "grey50",
                       midpoint = 0.085, limit = c(0,0.1), space = "Lab",
                       name="PGLS - P value (Pr>F)\nCutoff: P>0.1") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = ifelse(value<0.099, value, " ")), color = "black", size = 2.6) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    # legend.justification = c(1, 0),
    # legend.position = c(0.6, 0.7),
    legend.position = "right",
    legend.direction = "vertical")
#guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))
#ggsave(plot=plt_htmap_pvalue, filename="plt_PGLS_htmap_pranova.pdf", width=9.5, height=10)
ggsave(plot=plt_htmap_pvalue, filename="plt_PGLS_htmap_pranova_all.pdf", width=6.3, height=3)
# ggsave(plot=plt_htmap_pvalue, filename="plt_PGLS_htmap_pranova_all.png", width=7, height=3, dpi = 300)

# The adjusted r value plot ---------------------------------------------------------
plt_htmap_rsq = ggplot(data = melted_rsq_matrix_lambdaSignif, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "#08519C", mid = "white", high = "#CB181C",
                       na.value = "grey50",
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="PGLS - R value\nCutoff: abs(R)<0.25") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  coord_fixed() +
  geom_text(aes(Var2, Var1, label = ifelse(abs(value)>0.25, value, " ")), color = "black", size = 2.6) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    # legend.justification = c(1, 0),
    # legend.position = c(0.6, 0.7),
    legend.position = "right",
    legend.direction = "vertical")
#guides(fill = guide_colorbar(barwidth = 7, barheight = 1, title.position = "top", title.hjust = 0.5))
#ggsave(plot=plt_htmap_rsq, filename="plt_PGLS_htmap_rsq.pdf", width=9.5, height=10)
ggsave(plot=plt_htmap_rsq, filename="plt_PGLS_htmap_rsq_all.pdf", width=6.3, height=3)
# ggsave(plot=plt_htmap_rsq, filename="plt_PGLS_htmap_rsq_all.png", width=7, height=3, dpi = 300)
