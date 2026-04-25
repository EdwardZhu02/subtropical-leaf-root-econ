rm(list=ls())
library(tidyverse)
library(dplyr)
library(vegan)
library(ggpmisc) # stat_poly_eq function
library(reshape2)
library(patchwork)

load("outdata/traitDataFujian-SpAvg-phylo-step2.RData")
load("outdata/nonphyloPCAresult_RESLES-step3.RData")
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

varAngle_uptri_RESLES_R = func_get_upper_tri(PCAresult_RESLES_RDMC[["angles"]])
varAngleMelted_uptri_RESLES_R = reshape2::melt(varAngle_uptri_RESLES_R, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_LESRES_R=value)

varAngle_uptri_RESLES_L = func_get_upper_tri(PCAresult_RESLES_LPC[["angles"]])
varAngleMelted_uptri_RESLES_L = reshape2::melt(varAngle_uptri_RESLES_L, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_LESRES_L=value)

varAngle_uptri_RESLES_RL = func_get_upper_tri(PCAresult_RESLES_RDMCLPC[["angles"]])
varAngleMelted_uptri_RESLES_RL = reshape2::melt(varAngle_uptri_RESLES_RL, na.rm = TRUE) %>%
  dplyr::filter(value != 0) %>% # exclude self correlations
  dplyr::rename(trait1 = Var1, trait2=Var2, angle_LESRES_RL=value)


# varAngleMelted_combined = varAngleMelted_uptri_RESLES %>% 
#   right_join(varAngleMelted_uptri_RESLES_RDMC, by=join_by(trait1,trait2)) %>%
#   right_join(varAngleMelted_uptri_RESLES_LPC, by=join_by(trait1,trait2)) %>%
#   right_join(varAngleMelted_uptri_RESLES_RDMCLPC, by=join_by(trait1,trait2)) %>%
#   na.omit() %>% # exclude variables that can't form a point in scatter plot.
#   dplyr::mutate(pair_name = paste0(trait1,"-",trait2)) %>%
#   dplyr::mutate(pair_annotation = ifelse(
#     pair_name == "SRL-RD", "RES collaboration",
#     ifelse(pair_name == "RTD-RNC", "RES convervation",
#            ifelse(pair_name == "LMA-LNC", "LES", "other")))
#   )

# Combine the dataframes into a list
varAngleMelted_dflist = list(
  varAngleMelted_uptri_RESLES, 
  varAngleMelted_uptri_RESLES_R, 
  varAngleMelted_uptri_RESLES_L, 
  varAngleMelted_uptri_RESLES_RL
)
# Perform an outer join on all dataframes using reduce
varAngleMelted_combined = reduce(varAngleMelted_dflist, function(x, y) merge(x, y, by = c("trait1", "trait2"), all = TRUE))

varAngleMelted_combined = varAngleMelted_combined %>%
  dplyr::mutate(pair_name = paste0(trait1,"-",trait2)) %>%
  dplyr::mutate(pair_annotation = ifelse(
    pair_name == "SRL-RD", "RES collaboration",
    ifelse(pair_name %in% c("RTD-RNC"), "RES conservation",
           ifelse(pair_name %in% c("LMA-LNC"), "LES", "other")))
  )


# plot1: raw-raw+R ----------------------------------------------------------
plt_anglecorr_addR = ggplot(data=varAngleMelted_combined,
                               aes(x=angle_LESRES, y=angle_LESRES_R)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone", y="Backbone + RES co-pred") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

# plot2: raw-raw+L
plt_anglecorr_addL = ggplot(data=varAngleMelted_combined,
                              aes(x=angle_LESRES, y=angle_LESRES_L)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone", y="Backbone + LES co-pred") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

# plot3: raw-raw+RL
plt_anglecorr_addRL = ggplot(data=varAngleMelted_combined,
                                  aes(x=angle_LESRES, y=angle_LESRES_RL)) +
  geom_hline(yintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") + 
  geom_vline(xintercept = 90, color = "gray", linewidth = 0.4, linetype = "dotted") +
  geom_abline(intercept = 0, slope = 1, color = "gray", linewidth = 0.5, linetype = "dashed") + 
  geom_point(aes(fill=pair_annotation), size=2.7, alpha=0.9, shape=21, color="gray15") +
  #geom_smooth(method = lm, formula = y ~ x, se = T) +
  stat_poly_eq(formula = y ~ x, aes(label = paste0("P",paste(ifelse(..p.value..<0.0001, "<0.0001", paste0("=",..p.value..)), ..rr.label.., sep = "~~"))), parse = TRUE, rr.digits = 4) + 
  scale_fill_manual(values = c("#71BFB2", "gray70", "#237B9F", "#AD0B08")) + 
  theme_classic() + theme(plot.background = element_blank()) +
  labs(x="Backbone", y="Backbone + all co-pred") + 
  coord_cartesian(xlim = c(0, 180), ylim = c(0, 180))

plt_anglecorr_nonphy_total = plt_anglecorr_addR + plt_anglecorr_addL + plt_anglecorr_addRL + plot_layout(guides = "collect") & theme(legend.position='bottom', legend.title = element_blank())

suppressWarnings(ggsave(plot = plt_anglecorr_nonphy_total, filename = "outplts/SIFigures/FigS7ABC-anglecorr_nonphy.pdf", width = 8, height= 3.2))
