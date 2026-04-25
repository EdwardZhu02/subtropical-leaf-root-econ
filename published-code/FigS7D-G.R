rm(list=ls())
library(tidyverse)
library(dplyr)
library(vegan)
library(ggpmisc) # stat_poly_eq function
library(reshape2)
library(patchwork)

# phylo-PCA dependencies
library(ape) # general phylogenetic analysis
library(geiger) # compare taxa in data and tree
library(treeplyr) # general phylogenetic analysis
library(phytools) # general phylogenetic analysis, phyl.pca

load("outdata/traitDataFujian-SpAvg-phylo-step2.RData")
load("outdata/nonphyloPCAresult_RESLES-step3.RData")
load("outdata/phyloPCAresult_RESLES-step3.RData")
if (!dir.exists("outplts/SIFigures")) {
  dir.create("outplts/SIFigures", showWarnings = FALSE, recursive = TRUE)
}

### Function to estimate angles between pairs of vectors (2023, comments)
func_calc_angle_pcavars <- function(x, y){
  norm_x <- sqrt(sum(x^2)) #sd of eigenvector
  norm_y <- sqrt(sum(y^2))
  cosXY <- round(as.numeric((x %*% y) / (norm_x * norm_y)), 8) #dot product/magnitude of vectors
  # round introduced to avoid numerical problems with the acos function
  angle <- acos(cosXY) * 180 / pi
  return(angle)
}

# process angle matrix, extracting only the upper part
func_get_upper_tri = function(cormat){
  cormat[lower.tri(cormat)] = NA
  return(cormat)
}

# - non-phylo
varAngle_uptri_RESLES = func_get_upper_tri(PCAresult_RESLES[["angles"]])
varAngleMelted_uptri_RESLES = reshape2::melt(varAngle_uptri_RESLES, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_LESRES=value)

varAngle_uptri_RESLES_RDMC = func_get_upper_tri(PCAresult_RESLES_RDMC[["angles"]])
varAngleMelted_uptri_RESLES_RDMC = reshape2::melt(varAngle_uptri_RESLES_RDMC, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_LESRES_RDMC=value)

varAngle_uptri_RESLES_LPC = func_get_upper_tri(PCAresult_RESLES_LPC[["angles"]])
varAngleMelted_uptri_RESLES_LPC = reshape2::melt(varAngle_uptri_RESLES_LPC, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_LESRES_LPC=value)

varAngle_uptri_RESLES_RDMCLPC = func_get_upper_tri(PCAresult_RESLES_RDMCLPC[["angles"]])
varAngleMelted_uptri_RESLES_RDMCLPC = reshape2::melt(varAngle_uptri_RESLES_RDMCLPC, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_LESRES_RDMCLPC=value)

# - phylo
varAngle_uptri_phyloRESLES = func_get_upper_tri(phyloPCAresult_RESLES[["angles"]])
varAngleMelted_uptri_phyloRESLES = reshape2::melt(varAngle_uptri_phyloRESLES, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_phyloLESRES=value)

varAngle_uptri_phyloRESLES_RDMC = func_get_upper_tri(phyloPCAresult_RESLES_RDMC[["angles"]])
varAngleMelted_uptri_phyloRESLES_RDMC = reshape2::melt(varAngle_uptri_phyloRESLES_RDMC, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_phyloLESRES_RDMC=value)

varAngle_uptri_phyloRESLES_LPC = func_get_upper_tri(phyloPCAresult_RESLES_LPC[["angles"]])
varAngleMelted_uptri_phyloRESLES_LPC = reshape2::melt(varAngle_uptri_phyloRESLES_LPC, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_phyloLESRES_LPC=value)

varAngle_uptri_phyloRESLES_RDMCLPC = func_get_upper_tri(phyloPCAresult_RESLES_RDMCLPC[["angles"]])
varAngleMelted_uptri_phyloRESLES_RDMCLPC = reshape2::melt(varAngle_uptri_phyloRESLES_RDMCLPC, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_phyloLESRES_RDMCLPC=value)


# Combine the dataframes into a list
varAngleMelted_dflist = list(
  varAngleMelted_uptri_RESLES, 
  varAngleMelted_uptri_RESLES_RDMC, 
  varAngleMelted_uptri_RESLES_LPC, 
  varAngleMelted_uptri_RESLES_RDMCLPC,
  varAngleMelted_uptri_phyloRESLES,
  varAngleMelted_uptri_phyloRESLES_RDMC,
  varAngleMelted_uptri_phyloRESLES_LPC,
  varAngleMelted_uptri_phyloRESLES_RDMCLPC
)
# Perform an outer join on all dataframes using reduce
varAngleMelted_combined = reduce(varAngleMelted_dflist, function(x, y) merge(x, y, by = c("trait1", "trait2"), all = TRUE))

varAngleMelted_combined = varAngleMelted_combined %>%
  dplyr::mutate(pair_name = paste0(trait1,"-",trait2)) %>%
  dplyr::mutate(pair_annotation = ifelse(
    pair_name %in% c("SRL-RD","RD-SRR25","SRL-SRR25","SRL-SRA","RD-SRA"), "RES collaboration",
    ifelse(pair_name %in% c("RTD-RNC","RNC-RDMC","RTD-RDMC"), "RES conservation",
           ifelse(pair_name %in% c("LMA-LNC","LMA-LPC"), "LES", "other")))
  )


# plot1: raw, nonphylo-phylo ---------------------------------------------------
plt_anglecorr_phy_RLES = ggplot(data=varAngleMelted_combined,
                                aes(x=angle_LESRES, y=angle_phyloLESRES)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone", y="Backbone (phylo)") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

# plot2: raw+RDMC, nonphylo-phylo ----------------------------------------------
plt_anglecorr_phy_RLESRDMC = ggplot(data=varAngleMelted_combined,
                                    aes(x=angle_LESRES_RDMC, y=angle_phyloLESRES_RDMC)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone + RES co-pred", y="Backbone + RES co-pred (phylo)") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

# plot3: raw+LPC, nonphylo-phylo -----------------------------------------------
plt_anglecorr_phy_RLESLPC = ggplot(data=varAngleMelted_combined,
                                   aes(x=angle_LESRES_LPC, y=angle_phyloLESRES_LPC)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone + LES co-pred", y="Backbone + LES co-pred (phylo)") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

# plot3: raw+RDMC+LPC, nonphylo-phylo ------------------------------------------
plt_anglecorr_phy_RLESRDMCLPC = ggplot(data=varAngleMelted_combined,
                                       aes(x=angle_LESRES_RDMCLPC, y=angle_phyloLESRES_RDMCLPC)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone + all co-pred", y="Backbone + all co-pred (phylo)") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

plt_anglecorr_phy_total = plt_anglecorr_phy_RLES + plt_anglecorr_phy_RLESRDMC + plt_anglecorr_phy_RLESLPC + plt_anglecorr_phy_RLESRDMCLPC + plot_layout(guides = "collect") & theme(legend.position='bottom', legend.title = element_blank())

suppressWarnings(ggsave(plot = plt_anglecorr_phy_total, filename = "outplts/SIFigures/FigS7DEFG-anglecorr_phycomp.pdf", width = 5.5, height= 6))