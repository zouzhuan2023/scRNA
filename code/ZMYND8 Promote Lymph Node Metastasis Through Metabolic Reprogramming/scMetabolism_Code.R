

###############################scMetabolism
library(scMetabolism)
library(tidyverse)
library(rsvd)
library(Seurat)
library(pheatmap)
library(ComplexHeatmap)
library(dplyr)
library(ggsci)


load("OB_Tumor_filtered011.Rdata")
sce2<- OB_Tumor_filtered
rm(OB_Tumor_filtered)
sce2@meta.data$CB <- rownames(sce2@meta.data)
sample_CB <- sce2@meta.data %>% 
  group_by(seurat_clusters) %>% 
  sample_frac(1)
sce3 <- subset(sce2,CB %in% sample_CB$CB) 
head(sce3,2)
Idents(sce3) <- "seurat_clusters"
table(sce3$seurat_clusters)
countexp.Seurat <- sc.metabolism.Seurat(obj = sce3,  
                                        method = "AUCell", 
                                        imputation = F, 
                                        ncores = 5, 
                                        metabolism.type = "KEGG")

save(countexp.Seurat,file="countexp.Seurat.Rdata")
signatures_KEGG_metab <- system.file("data", "KEGG_metabolism_nc.gmt", 
                                     package = "scMetabolism")
signatures_KEGG_metab
score <- countexp.Seurat@assays$METABOLISM$score
score[1:4,1:4]
score_change <- score %>% 
  select_all(~str_replace_all(., "\\.", "-"))  
identical(colnames(score_change) , rownames(countexp.Seurat@meta.data))
countexp.Seurat@meta.data <- cbind(countexp.Seurat@meta.data,t(score_change) )

library(RColorBrewer)
library(patchwork)

input.pathway<-c(
  "Pantothenate and CoA biosynthesis",
  "Oxidative phosphorylation",
  "Glycine, serine and threonine metabolism",
  "Citrate cycle (TCA cycle)",
  "Biosynthesis of unsaturated fatty acids",
  "Arginine and proline metabolism",
  "Arginine biosynthesis",
  "Glycolysis / Gluconeogenesis",
  "Purine metabolism",
  "Pyruvate metabolism"
)

p <- DotPlot.metabolism(obj = countexp.Seurat, 
                        pathway = input.pathway, 
                        phenotype = "seurat_clusters", 
                        #            split.by = "group", 
                        norm = "y") +
  theme(
    axis.text.x = element_text(size = 15, angle = 0, hjust = 1), 
    axis.text.y = element_text(size = 12)                        
  )


ggsave("metabolism.pdf", plot = p, width = 7, height = 5)

