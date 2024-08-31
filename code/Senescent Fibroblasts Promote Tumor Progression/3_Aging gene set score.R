###Empty the environment####
rm(list=ls())
###Load the package####
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(tidyverse)
library(stringr)
library(harmony)
library(ggplot2)
library(ggpubr)
library(AUCell) 
library(clusterProfiler)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsignif)
###Set the working path to###
setwd("/path")
###Read the data###
load("/path/sencelltype.Rdata")
###Assigns the value###
scRNA <- CAF_AURKB
###Read the selection scoring gene set####
gene_list <- read.table("/path/senescence.txt", header = T)
table(scRNA@meta.data$CellType)
Idents(scRNA) <- "CellType"
###Convert genes to character form####
genes_vector <- as.character(gene_list$gene)
##AddModuleScore score####
##Set default analysis to#####
DefaultAssay(scRNA) <- "RNAnew"
##Add module score####
scRNA <- AddModuleScore(scRNA,
                        features = gene_list,
                        ctrl = 100, 
                        name = "AddModuleScore")
meta=scRNA@meta.data
colnames(scRNA@meta.data)[29] <- 'senescence_Score' 
table(scRNA@meta.data$CellType)
###Extract the required data###
data <- FetchData(scRNA, vars = c("CellType", "senescence_Score"))
###Check data###
head(data)
###Create violin diagram and add box diagram###
p <- ggplot(data, aes(x = CellType, y = senescence_Score, fill = CellType)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA, color = "black") +
  labs(title = "senescence_Score") +
  xlab("CellType") +
  ylab("senescence_Score") +
  scale_fill_manual(values = c("nonsenCAF" = "#E5D2DD", "senCAF" = "#53A85F")) +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(color = "black"),
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(colour = "black"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
avg_scores <- scRNA@meta.data %>%
  group_by(CellType) %>%
  summarize(avg_score = mean(senescence_Score))

p <- p + geom_signif(comparisons = list(c("nonsenCAF", "senCAF")),
                     map_signif_level = TRUE,
                     textsize = 3,vjust = 0.5, y_position = c(0.8, 0.86, 0.75,1.18,0.75,0.62,0.55))

print(p)
ggsave(filename = "senescence_Score.pdf", p, width = 6, height = 4, dpi = 300)
dev.off()
