library(DESeq2)
library(limma)
library(pheatmap)
library(magrittr)
library(celldex)
library(tidyverse)
library(stringr)
library(harmony)
library(dplyr)
library(Seurat)
library(patchwork)
library(bigmemory)
library(ggplot2)
library(ggsci)
library(ggsignif)
library(ggpubr)
library(ggalluvial)
library(reshape2)
table(OB_Tumor_filtered$CNVtype)
load("OB_Tumor_filtered011.Rdata")
Idents(OB_Tumor_filtered) <- "seurat_clusters"
table(OB_Tumor_filtered$seurat_clusters)
p <- VlnPlot(OB_Tumor_filtered, features = c("SLPI"), ncol = 1,pt.size = 0,
             cols = c("#E78376", 
                      "#81CAD9",
                      "#68B9AA",
                      "#7888AA",
                      "#F0B7A6",
                      "#ABB1C9",
                      "#B1DBD1"))
p
ggsave("SLPI.pdf", p, width = 5, height = 4)
###########################
Idents(OB_Tumor_filtered) <- "Secgroup"
table(OB_Tumor_filtered$Secgroup)
p <- VlnPlot(OB_Tumor_filtered, features = c("SLPI"), ncol = 1,pt.size = 0,
             cols = c("#1F77B4",
                      "#7ABFE2",
                      "#E377C2",
                      "#fad1e0"))
ggsave("SLPI_group.pdf", p, width = 5, height = 4)


################################ZMYND8
Idents(OB_Tumor_filtered) <- "seurat_clusters"
p <- VlnPlot(OB_Tumor_filtered, features = c("ZMYND8"), ncol = 1,pt.size = 0,
             cols = c("#E78376", 
                      "#81CAD9",
                      "#68B9AA",
                      "#7888AA",
                      "#F0B7A6",
                      "#ABB1C9",
                      "#B1DBD1"))
p
ggsave("ZMYND8.pdf", p, width = 5, height = 4)

Idents(OB_Tumor_filtered) <- "Secgroup"
table(OB_Tumor_filtered$Secgroup)
p <- VlnPlot(OB_Tumor_filtered, features = c("ZMYND8"), ncol = 1,pt.size = 0,
             cols = c("#1F77B4",
                      "#7ABFE2",
                      "#E377C2",
                      "#fad1e0"))
p
ggsave("ZMYND8_group.pdf", p, width = 5, height = 4)



##########################NPW
p <- VlnPlot(OB_Tumor_filtered, features = c("NPW"), ncol = 1,pt.size = 0,
             cols = c("#E78376", 
                      "#81CAD9",
                      "#68B9AA",
                      "#7888AA",
                      "#F0B7A6",
                      "#ABB1C9",
                      "#B1DBD1"))
p
ggsave("NPW.pdf", p, width = 5, height = 4)
####################
Idents(OB_Tumor_filtered) <- "Secgroup"
table(OB_Tumor_filtered$Secgroup)
p <- VlnPlot(OB_Tumor_filtered, features = c("NPW"), ncol = 1,pt.size = 0,
             cols = c("#1F77B4",
                      "#7ABFE2",
                      "#E377C2",
                      "#fad1e0"))
p
ggsave("NPW_group.pdf", p, width = 5, height = 4)
################################################

dir <- dir("path")
dir <- paste0("path/", dir)
samples_name = c("BC2","BC3","BC5","BC6","BC16","BC21","BC22","S1","S2","S3","S4")
scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project = samples_name[i], min.cells = 3, min.features = 200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id = samples_name[i])   
  if(T){scRNAlist[[i]][["percent.mt"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern = "^MT-")}
}

names(scRNAlist) <- samples_name
scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])
scRNA$proj <- rep("10x", ncol(scRNA))
table(scRNA$orig.ident)
head(scRNA@meta.data)
save(scRNA, file = 'scRNA_merge.Rdata')
dir.create("QC")

theme.set1 = theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")
group = "orig.ident"
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i],raster=FALSE) + theme.set2 + NoLegend()}
p1 <- wrap_plots(plots = plots, nrow = 2)    
ggsave("QC/vlnplot_before_qc.pdf", plot = p1, width = 14, height = 7) 
ggsave("QC/vlnplot_before_qc.png", plot = p1, width = 14, height = 7)

p2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt",raster=FALSE)
p3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster=FALSE)
plotc <- p2 + p3 + plot_layout(guides = "collect")
ggsave("QC/FeatureScatter_before_qc.pdf", plotc, width = 15, height = 7) 
ggsave("QC/FeatureScatter_before_qc.png", plotc, width = 15, height = 7)

scRNA<- subset(scRNA, subset = nFeature_RNA > 500 & nFeature_RNA < 4500 & percent.mt < 10)
scRNA<- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA<- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
dim(scRNA)

save(scRNA, file = 'QC/scRNA_QC.Rdata')

scRNA1 <- scRNA
rm(scRNA)
plots = list()
for(i in seq_along(plot.featrures)){plots[[i]] = VlnPlot(scRNA1, group.by=group, pt.size = 0, 
                                                         features = plot.featrures[i],raster=FALSE) + theme.set2 + NoLegend()}
p4 <- wrap_plots(plots = plots, nrow = 2)     
ggsave("QC/vlnplot_after_qc.pdf", plot = p4, width = 14, height = 10) 
ggsave("QC/vlnplot_after_qc.png", plot = p4, width = 14, height = 10)
p5 <- FeatureScatter(scRNA1, feature1 = "nCount_RNA", feature2 = "percent.mt",raster=FALSE)
p6 <- FeatureScatter(scRNA1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",raster=FALSE)
plotc <- p5 + p6 + plot_layout(guides = "collect")
ggsave("QC/FeatureScatter_after_qc.pdf", plotc, width = 15, height = 7) 
ggsave("QC/FeatureScatter_after_qc.png", plotc, width = 15, height = 7)

###########################
all.genes <- rownames(scRNA1)
scRNA1 <- ScaleData(scRNA1, features = all.genes)
scRNA1 <- RunPCA(scRNA1, features = VariableFeatures(object = scRNA1))
VizDimLoadings(object = scRNA1, dims = 1:4, reduction = "pca",nfeatures = 20)
pdf(file = "Harmony/PCA.pdf")
DimPlot(object = scRNA1, reduction = "pca",raster=FALSE)
dev.off()
print(scRNA1[["pca"]], dims = 1:5, nfeatures = 5)

pdf(file = "Harmony/ElbowPlot.pdf")
ElbowPlot(scRNA1)
dev.off()


library(harmony)
dir.create("Harmony")
scRNA <- scRNA1
rm(list=setdiff(ls(), "scRNA"))
gc(verbose = FALSE)


scRNA_harmony <- scRNA %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

scRNA_harmony <- scRNA_harmony %>%
  RunPCA(npcs = 14, verbose = FALSE) %>%
  RunTSNE(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  RunUMAP(reduction="harmony", dims=1:15)
save(scRNA_harmony, file = "Harmony/scRNA_harmony.Rdata")

table(scRNA_harmony$orig.ident)

Idents(scRNA_harmony) <- "orig.ident"
scRNA_harmony@meta.data$group = ifelse(scRNA_harmony@meta.data$orig.ident %in%
                                         c("S1","S2","S3","S4"),"NT","PT")
scRNA_harmony$group <- as.factor(scRNA_harmony$group) 
table(scRNA_harmony$group)


rm(scRNA)
p <- ElbowPlot(scRNA_harmony, ndims = 20)
ggsave("Harmony/ElbowPlot50.pdf", p, width = 8, height = 4)
ggsave("Harmony/ElbowPlot50.png", p, width = 8, height = 4)


scRNA_harmony <-FindClusters(scRNA_harmony,resolution = 0.05)
p<-VlnPlot(scRNA_harmony,features = c("ALPL","RUNX2","IBSP","CLEC11A","CTSK","CDH11"),pt.size = 0, ncol = 3,raster=FALSE)

dir.create("003")
ggsave("OB.pdf",width = 8, height = 6)
save(scRNA_harmony, file = "scRNA_harmony005.Rdata")

Idents(scRNA_harmony) <- "seurat_clusters"
scRNA_harmony$seurat_clusters<-as.factor(scRNA_harmony$seurat_clusters)
table(scRNA_harmony$seurat_clusters)

Idents(scRNA_harmony) <- "seurat_clusters"
OB_3 <- subset(scRNA_harmony,idents=("3"))
table(OB_3$seurat_clusters)
OB_3$seurat_clusters <- OB_3@active.ident
OB_3$seurat_clusters <- as.factor(OB_3$seurat_clusters)


Idents(OB_3) <- "group"
p<-VlnPlot(OB_3,features = c("NPW"),pt.size = 0, ncol = 1,raster=FALSE)
p
ggsave("NPW_validation.pdf",width = 4, height = 4)



############################
###############################
load("NK_harmony.Rdata")
markers <- c("KLRD1", "NKG7", "GZMK", "GZMA", "KLRG1", "KLRF1", "PRF1", "FCGR3A")
Idents(NK)="Secgroup"
p <- VlnPlot(object = NK, features = markers, pt.size = 0, combine = TRUE,ncol = 4,
             cols = c("#1F77B4",
                      "#7ABFE2",
                      "#E377C2",
                      "#fad1e0",
                      "#FB6346",
                      "#71c3aa")) +#group.by = "Secgroup", 
  theme(plot.title = element_text(hjust = 0.5))

ggsave("NKgene.pdf", plot = p, width = 12, height = 6)