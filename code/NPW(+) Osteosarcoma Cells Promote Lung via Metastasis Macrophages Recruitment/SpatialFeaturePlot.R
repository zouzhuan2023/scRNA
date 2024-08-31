######load packages and setting
library(spacexr)
library(Matrix)
library(doParallel)
library(ggplot2)
library(Seurat)
library(data.table)
library(cowplot)
library(STdeconvolve)
library(ggplot2)
library(ggsci)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  cache = TRUE,
  out.width = "100%"
)

st_data <- Load10X_Spatial(data.dir = "/path",
                           filename = "filtered_feature_bc_matrix.h5",
                           assay = "Spatial")

st_data[["percent.mt"]] <- PercentageFeatureSet(st_data, pattern = "^MT-")
minFeature = 400  
minCount = 200   
maxMT = 20        
minCell = 20 

counts <- GetAssayData(st_data, slot = "counts", assay = "Spatial")
genes_to_keep <- rowSums(counts > 0) >= minCell
st_data <- st_data[genes_to_keep, ]

st_data <- subset(st_data, subset = nFeature_Spatial >= minFeature & percent.mt <= maxMT & nCount_Spatial >= minCount)
pdf(file = "SpatialPlot.pdf", width = 8,height = 8)
SpatialPlot(st_data,pt.size.factor=2)
dev.off()
library(SeuratDisk)
str(st_data@meta.data)
str(st_data@reductions)

st_data <- NormalizeData(st_data, normalization.method = "LogNormalize", scale.factor = 10000,seed.use =2024)
st_data <- FindVariableFeatures(st_data, selection.method = "vst", nfeatures = 2000)

p <- SpatialFeaturePlot(st_data, 
                        features = c("NPW","CD163","C1QA"),
                        pt.size.factor=2.5)
ggsave("NPW_Mar.pdf", plot = p, width = 10, height = 8, dpi = 600)


npw_expr <- GetAssayData(st_data, slot = "counts")["NPW", ]
cd163_expr <- GetAssayData(st_data, slot = "counts")["CD163", ]
c1qa_expr <- GetAssayData(st_data, slot = "counts")["C1QA", ]


cor_npwc1qa <- cor(npw_expr, c1qa_expr)
cor_npwcd163 <- cor(npw_expr, cd163_expr)

library(ggplot2)
##############
df_npwc1qa <- data.frame(NPW = npw_expr, C1QA = c1qa_expr)
cor_test <- cor.test(df_npwc1qa$NPW, df_npwc1qa$C1QA)
cor_npwc1qa <- cor_test$estimate
p_value <- cor_test$p.value
p1 <- ggplot(df_npwc1qa, aes(x = NPW, y = C1QA)) +
  geom_point(color = "#123168", size = 3, alpha = 0.6) + 
  geom_smooth(method = "lm", color = "red", se = FALSE) + 
  annotate("text", x = 0.1, y = max(df_npwc1qa$C1QA) * 0.9, label = paste("R =", round(cor_npwc1qa, 2), "\nP =", format.pval(p_value, digits = 2)), color = "red", size = 5, hjust = 0, vjust = 1) +
  theme_classic(base_size = 15) + 
  ggtitle(paste("Correlation: ", round(cor_npwc1qa, 2))) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        plot.margin = margin(15, 15, 15, 15))  

print(p1)
ggsave(filename = "NPW_C1QA.pdf", plot = p1, width = 5, height = 5, units = "in", dpi = 600)


###################
df_npwcd163 <- data.frame(NPW = npw_expr, CD163 = cd163_expr)
cor_test <- cor.test(df_npwcd163$NPW, df_npwcd163$CD163)
cor_npwcd163 <- cor_test$estimate
p_value <- cor_test$p.value
p2 <- ggplot(df_npwcd163, aes(x = NPW, y = CD163)) +
  geom_point(color = "#123168", size = 3, alpha = 0.6) + 
  geom_smooth(method = "lm", color = "red", se = FALSE) + 
  annotate("text", x = 0.1, y = max(df_npwcd163$CD163) * 0.9, label = paste("R =", round(cor_npwcd163, 2), "\nP =", format.pval(p_value, digits = 2)), color = "red", size = 5, hjust = 0, vjust = 1) +
  theme_classic(base_size = 15) + 
  ggtitle(paste("Correlation: ", round(cor_npwcd163, 2))) + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        plot.margin = margin(15, 15, 15, 15)) 

print(p2)
ggsave(filename = "NPW_CD163.pdf", plot = p2, width = 5, height = 5, units = "in", dpi = 600)

