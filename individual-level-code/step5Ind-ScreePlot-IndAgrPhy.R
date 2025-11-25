rm(list=ls())

library(dplyr)
library(tidyverse)
library(tidyr)
library(reshape2)

# Adopted from step17-phylo-informed-PCA-log.R
# Step 5: Compare with non-phylo PCA, grouped barplot --------------------------
screeData_compare_RLES = data.frame( # backbone
  Dimensions = c("1","2","3","4","5","6"),
  NonPhyloIndv = c(33.3,26.4,24.2,9.3,6.8,0.2),
  NonPhylo = c(34.4,27.6,21.5,10,6.3,0.2),
  Phylo = c(35,27,22.1,9.5,6.3,0.2))

screeData_compare_RLES_root_copred = data.frame( # RDMC,Rr25,SRA
  Dimensions = c("1","2","3","4","5","6","7","8","9"),
  NonPhyloIndv = c(35.3,25.9,17.8,8.8,6.6,3.4,1.9,0.1,0),
  NonPhylo = c(34.6,27.7,17,9.4,6.5,3.3,1.2,0.1,0),
  Phylo = c(35.8,27.2,17.1,9.7,5.5,3.3,1.3,0.1,0))

screeData_compare_RLES_leaf_copred = data.frame( # LPC
  Dimensions = c("1","2","3","4","5","6","7"),
  NonPhyloIndv = c(30,28,21.3,8.6,6.9,5,0.1),
  NonPhylo = c(32,27.6,21.2,9,7.2,2.8,0.1),
  Phylo = c(32.7,27.2,21.1,9.1,6.8,3,0.1))

screeData_compare_RLES_all_copred = data.frame( 
  Dimensions = c("1","2","3","4","5","6","7","8","9","90"),
  NonPhyloIndv = c(31.8,25.5,19,8.5,6.4,4,3.1,1.6,0.1,0),
  NonPhylo = c(31.3,25.8,20.4,9.1,5.9,4.5,1.8,1.1,0.1,0),
  Phylo = c(32.4,25.7,19.9,9.4,5,4.4,1.9,1.2,0.1,0))

#screeData_compare_touse = screeData_compare_RLES
#screeData_compare_touse = screeData_compare_RLESLPC
#screeData_compare_touse = screeData_compare_RLESRDMC
screeData_compare_touse = screeData_compare_RLES_all_copred
screeData_compare_touse = melt(screeData_compare_touse,variable.name="Category",value.name = "ExplainedVariance")

# Grouped bar plot
plt_scree_compare = ggplot(screeData_compare_touse) +
  geom_bar(aes(x=Dimensions, y=ExplainedVariance, fill=Category), color="black", # border
           stat="identity", width=0.8, position="dodge",alpha=0.78) + 
  
  geom_point(aes(x=Dimensions, y=ExplainedVariance, color=Category), size=2, position=position_dodge(width = 0.6)) +
  geom_line(aes(x=Dimensions, y=ExplainedVariance, color=Category, group = Category), linewidth=1, position=position_dodge(width = 0.5)) +
  
  scale_fill_manual(values = c("#A9082C","#FDB76D","#4B61A8")) + # bar coloring
  scale_color_manual(values = c("#632024","#A35D13","#25344F")) + # point and line coloring


  geom_text(aes(x = Dimensions, y = ExplainedVariance, group = Category, label = paste0(ExplainedVariance, "%")),
            position = position_dodge(width = 0.6),size = 3.5,vjust = 0.8,hjust=-0.7,angle=90) + 
  labs(title="R-LES backbone + all co-predictors", y="Variance Explained (%)") + 
  theme_classic() +
  scale_y_continuous(limits = c(0, 42), expand = c(0, 0)) + # Removes the blank space
  theme(legend.title = element_blank(), 
        legend.position = c(1, 1),
        legend.justification = c(1, 1)
        # legend.background = element_rect(fill = "white", color = "black")
  )

plt_scree_compare

ggsave(filename = "individual-level-code/Ind_PCAscree-RLESall.pdf", plot = plt_scree_compare, width = 3.9, height = 3.5)
#ggsave(filename = "plt_PCAcompare_scree-RESLESLPC.pdf", plot = plt_scree_compare, width = 3.9, height = 3.5)
#ggsave(filename = "plt_PCAcompare_scree-RESLESRDMC.pdf", plot = plt_scree_compare, width = 3.9, height = 3.5)
#ggsave(filename = "plt_PCAcompare_scree-RESLESLPCRDMC.pdf", plot = plt_scree_compare, width = 3.9, height = 3.5)

