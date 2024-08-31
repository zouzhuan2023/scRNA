library(GSVA)
library(Seurat)
library(dplyr)
library(readr)
library(ggplot2)
gsvatestdata <- OB_Tumor_filtered

GeneList_immPort <- read_delim("path/GeneList.txt", 
                               delim = "\t", 
                               escape_double = FALSE, 
                               trim_ws = TRUE)
gene_set_list <- split(GeneList_immPort$Symbol, GeneList_immPort$Category)
Idents(gsvatestdata) <- "seurat_clusters"

Average <- AverageExpression(gsvatestdata, assays = NULL, features = NULL, return.seurat = FALSE,
                             add.ident = NULL, slot = "data", use.scale = FALSE, use.counts = FALSE, verbose = TRUE)
data <- as.matrix(Average$RNA)
GSVAresult <- gsva(data, gene_set_list, min.sz=10,max.sz=Inf, tau=1, method="gsva", kcdf="Poisson",
                   mx.diff=TRUE, abs.ranking=FALSE, verbose=TRUE, parallel.sz=10)
write.csv(GSVAresult, file = "OB_immPort_GSVA.csv")
t <- t(scale(t(GSVAresult)))
write.csv(t, file = "OB_immPort_GSVA_t.csv")

t1 <- read.csv("OB_immPort_GSVA_t.csv", header = T, row.names = 1)
t1 <- as.matrix(t1)
new_colnames <- c("0", "1", "2", "3", "4", "5", "6")
colnames(t1) <- new_colnames
p <- pheatmap::pheatmap(t1, #
                        cluster_rows = T, 
                        cluster_cols = F, 
                        show_colnames = T,
                        fontsize_row = 8,
                        fontsize_col = 8, 
                        border_color = "white",
                        color = colorRampPalette(rev(brewer.pal(10, "RdBu")))(30),
                        angle_col = 315)
ggsave("GSVA.pdf", p, width = 6, height = 3)
