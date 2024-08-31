
rm(list=ls())
library(limma)
library(Seurat)
library(magrittr)
library(celldex)
library(SingleR)
library(tidyverse)
library(stringr)
library(harmony)
library(ggpubr)
library(AUCell) 
library(dplyr)
library(ggplot2)
library(ggsignif)
library(clusterProfiler)
table(OB_Tumor_filtered$seurat_clusters)
scRNA <- OB_Tumor_filtered
Idents(scRNA) <- "seurat_clusters"
table(scRNA$seurat_clusters)
gene_list <- read.table("NaturalKiller_Cell_Cytotoxicity.txt", header = T)
Idents(scRNA) <- "seurat_clusters"
gene_sets <- list.files(pattern = "\\.txt$")
for (gene_set_file in gene_sets) {
  genes <- read.table(gene_set_file, header = TRUE, stringsAsFactors = FALSE)$gene
  genes <- intersect(genes, rownames(scRNA@assays$RNA@data))
  score_name <- gsub(".txt", "", gene_set_file) 
  scRNA <- AddModuleScore(scRNA, features = list(genes), name = paste0(score_name, "_Score"))
}

data <- FetchData(scRNA, vars = c("seurat_clusters", "NaturalKiller_Cell_Cytotoxicity_Score1"))
head(data)
library(ggplot2)
library(ggsignif)

cluster_colors <- c("#E78376", 
                    "#81CAD9",
                    "#68B9AA",
                    "#7888AA",
                    "#F0B7A6",
                    "#ABB1C9",
                    "#B1DBD1")

p1 <- ggplot(data, aes(x = seurat_clusters, y = NaturalKiller_Cell_Cytotoxicity_Score1, fill = seurat_clusters)) +
  geom_violin(trim = FALSE) +
  xlab("Subcelltype") +
  ylab("NaturalKiller_Cell_Cytotoxicity_Score1") +
  scale_fill_manual(values = cluster_colors) +
  theme_minimal() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_line(color = "black"),  
        axis.line = element_line(color = "black"), 
        strip.background = element_rect(fill = "white", colour = "white"),
        strip.text = element_text(colour = "black"))

p1

avg_scores <- tapply(scRNA@meta.data$NaturalKiller_Cell_Cytotoxicity_Score1, scRNA@meta.data$seurat_clusters, mean)
avg_scores_df <- data.frame(cluster = names(avg_scores), avg_score = avg_scores)

p2 <- p1 + geom_text(data = avg_scores_df, 
                     aes(x = cluster, y = avg_score, label = round(avg_score, 4)),
                     inherit.aes = FALSE, 
                     position = position_nudge(y = 0.0), 
                     size = 4)
print(p2)
ggsave(filename = "NaturalKiller_Cell_Cytotoxicity_Score1.pdf", plot = p2, width = 6, height = 3)




##########################################

scRNA <- OB_Tumor_filtered
Idents(scRNA) <- "seurat_clusters"
table(scRNA$seurat_clusters)
gene_sets <- list.files(pattern = "\\.txt$")
for (gene_set_file in gene_sets) {
  genes <- read.table(gene_set_file, header = TRUE, stringsAsFactors = FALSE)$gene
  genes <- intersect(genes, rownames(scRNA@assays$RNA@data))
  score_name <- gsub(".txt", "", gene_set_file)  
  scRNA <- AddModuleScore(scRNA, features = list(genes), name = paste0(score_name, "_Score"))
}
data <- FetchData(scRNA, vars = c("seurat_clusters", "Lung_Met1_Score1"))
cluster_colors <- c("#E78376", 
                    "#81CAD9",
                    "#68B9AA",
                    "#7888AA",
                    "#F0B7A6",
                    "#ABB1C9",
p<-ggviolin(data, 
            x="seurat_clusters", y="Lung_Met1_Score1", 
            width = 0.8,color = "black",
            fill="seurat_clusters",
            xlab = F, 
            add = 'mean_sd', 
            bxp.errorbar=T,
            bxp.errorbar.width=0.05, 
            size=0.5,
            palette = "npg", 
            legend = "right")+
  coord_flip() 
ggsave(filename = "Lung_Met1_Score1_violin_plot.pdf", plot = p, width = 5, height = 3, dpi = 600)



