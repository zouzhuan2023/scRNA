###Empty the environment####
rm(list=ls())
###Load the package####
library(limma)
library(magrittr)
library(SingleR)
library(monocle)
library(tidyverse)
library(stringr)
library(harmony)
library(Seurat)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
table(scRNA$leiden)
###Set the working path to###
setwd("/path")
###Read the data###
load("/path/st_adata_seurat_qc_gene.Rdata")
table(scRNA@meta.data$celltype)
dir.create("Fibroblast")
table(scRNA$leiden)
Idents(scRNA)="celltype"
Fibroblast <- subset(scRNA,ident = c("Fibroblast"))
save(Fibroblast,file="Fibroblast.Rdata")
load("Fibroblast.Rdata")
table(Fibroblast$leiden)
table(scRNA$leiden)
###Data normalization###
Fibroblast <- NormalizeData(Fibroblast, normalization.method = "LogNormalize", scale.factor = 10000)
Fibroblast <- FindVariableFeatures(Fibroblast, selection.method = "vst", nfeatures = 2000)
#####FibroblastPCA###
all.genes <- rownames(Fibroblast)
Fibroblast <- ScaleData(Fibroblast, features = all.genes)
Fibroblast <- RunPCA(Fibroblast, features = VariableFeatures(object = Fibroblast))
pdf(file = "FibroblastPCA.pdf")
DimPlot(object = Fibroblast, reduction = "pca")
dev.off()
pdf(file = "Fibroblast ElbowPlot.pdf")
ElbowPlot(Fibroblast)
dev.off()
Fibroblast <- FindNeighbors(Fibroblast, dims = 1:5)
Fibroblast <- RunUMAP(Fibroblast, dims = 1:5)
Fibroblast <- RunTSNE(Fibroblast, dims = 1:5)
save(Fibroblast, file = "Fibroblast_harmony.Rdata")
###Check suitable resolution###
resolution <- c(seq(0.01,0.1,0.01))
Fibroblast <- FindClusters(Fibroblast,
                           resolution = resolution,
                           verbose = TRUE)

head(Fibroblast@meta.data)
###Count the number of cells in each cluster within each resolution###
apply(Fibroblast@meta.data[, grep("RNA_snn", colnames(Fibroblast@meta.data))], 2, table)
library(clustree)
pdf(file = "tree.pdf")
clustree(Fibroblast@meta.data, prefix ="RNAnew_snn_res.")
dev.off()
save(Fibroblast, file = "Fibroblast_harmony_tree.Rdata")
###Determine the resolution####
load("Fibroblast_harmony_tree.Rdata")
Fibroblast <-FindClusters(Fibroblast,resolution = 0.05)
table(Fibroblast$seurat_clusters)
p <- DimPlot(Fibroblast, reduction = "tsne",
             group.by ="seurat_clusters"
             ,label = T,raster=FALSE)+
  scale_color_manual(values = c('#00766D',
                                '#FF8066',
                                '#FF6F91',
                                '#EFB04B',
                                "#FBD51A",
                                '#845EC0',
                                '#C68EBB',
                                '#84852B',
                                '#CBCFD3',
                                '#C9B9A0',
                                '#7CD4C7',
                                '#819ABB',
                                '#639494',
                                '#965A4F',
                                "#28507D",
                                '#3C5866',
                                '#4B4453'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  labs(title = "Clusters")

ggsave(filename = "fig1 Cluster_tsne.pdf",p,height = 6,width = 8)
###Visualize cell clustering###
p<-DimPlot(Fibroblast, reduction = "umap",
           group.by ="seurat_clusters"
           ,label = T,raster=FALSE)+
  scale_color_manual(values = c('#00766D',
                                '#FF8066',
                                '#FF6F91',
                                '#EFB04B',
                                "#FBD51A",
                                '#845EC0',
                                '#C68EBB',
                                '#84852B',
                                '#CBCFD3',
                                '#C9B9A0',
                                '#7CD4C7',
                                '#819ABB',
                                '#639494',
                                '#965A4F',
                                "#28507D",
                                '#3C5866',
                                '#4B4453'))+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+
  labs(title = "Clusters")

ggsave(filename = "fig2 Cluster_umap.pdf",height = 6,width = 8)
table(Fibroblast$seurat_clusters)  
##Identify differentially expressed genes for each cluster###
logFCfilter =0.25
adjPvalFilter=0.05
Fibroblast.markers <- FindAllMarkers(object = Fibroblast,
                                     only.pos = FALSE,
                                     min.pct = 0.25,
                                     logfc.threshold = logFCfilter)
sig.markers=Fibroblast.markers[(abs(as.numeric(as.vector(Fibroblast.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(Fibroblast.markers$p_val_adj))<adjPvalFilter),]
write.csv(sig.markers,file="clusterMarkers.csv",sep="\t",row.names=F,quote=F)
save(Fibroblast, file = "Fibroblast_harmony_Undefined.Rdata")
load("/path/Fibroblast_harmony_Undefined.Rdata")
dev.new()
pdf(file = "markers_Vlinplot.pdf")
VlnPlot(Fibroblast, c("DCN","TAGLN","POSTN","LUM","GSN"),group.by = "seurat_clusters")
dev.off()
table(Fibroblast@meta.data$seurat_clusters)
markers <- c("CYP4B1","MYOC","ADH1B","CFD","SCN7A","AKAP6", "CADPS","DLGAP2","FGF14","PDZRN4","SPP1","SCG2","RAMP1","RARRES1","CPXM1",
              "ZNF560","IGF2BP3", "AURKB", "TDO2","WDR62","COL2A1","COL9A2","COL9A3","S100A1","COL11A2","MYF6","MYF5","FITM1","DES","RAPSN"
)
p<- DotPlot(Fibroblast, features = markers) +
  scale_colour_gradient(low = "#71B3C2", high = "red") +
  coord_flip() + # Flip the axes
  theme(
    axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 12),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  ) +
  xlab("markers1") + # Set the x-axis label
  ylab("seurat_clusters") # Set the y-axis label)
print(p)
ggsave("Bubble chart.pdf", plot = p, width = 10, height = 8)
###Define cells###
rm(list=ls())
setwd("/path")
load("/path/Fibroblast_harmony_Undefined.Rdata")
View(Fibroblast@meta.data$seurat_clusters)
Idents(Fibroblast)<-"seurat_clusters"
new.cluster.ids <- c("CAF_CFD","CAF_AKAP6","CAF_SPP1","CAF_AURKB","mCAF","CAF_MYF6")
names(new.cluster.ids) <- levels(Fibroblast)
Fibroblast<- RenameIdents(Fibroblast, new.cluster.ids)
Fibroblast@meta.data$CellType <- Fibroblast@active.ident
table(Fibroblast$CellType)
table(Fibroblast$seurat_clusters)
save(Fibroblast, file = "Fibroblast_harmony_Annotation.Rdata")
Idents(Fibroblast)="CellType"
p=DimPlot(Fibroblast, reduction = "tsne",
          group.by ="CellType"
          ,label = F,raster=FALSE)+
  scale_color_manual(values = c(
    my36colors <-c('#E59CC4', '#AB3282','#E63863',
                   '#58A4C3', '#8C549C','#3A6963',
                   '#9FA3A8', '#E0D4CA', '#5F3D69', '#585658', '#BD956A', '#E4C755', '#F7F398',
                   '#AA9A59', #'#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', #'#B53E2B',
                   '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#C5DEBA',
                   '#968175'
    )))+
  theme(panel.border = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  labs(title = "CellType")
p
ggsave(filename = "Fibroblast_Annotation_tsne.pdf",p,height =6,width =8)
###Add grouping information###
Idents(Fibroblast) <- "batch"
table(Fibroblast$batch)
Fibroblast@meta.data$Secgroup = ifelse(Fibroblast@meta.data$batch %in%
                                         c("MLN07","MLN09"),"MLN",
                                       ifelse(Fibroblast@meta.data$batch %in%
                                                c("NLN08","NLN09","NLN10","NLN11"),"NLN",
                                              ifelse(Fibroblast@meta.data$batch %in%
                                                       c("NT01","NT03","NT07","NT08","NTfan","NTli"),"NT",
                                                     ifelse(Fibroblast@meta.data$batch %in%
                                                              c("PT07","PTyi"),"Lymph-MPT",      
                                                            ifelse(Fibroblast@meta.data$batch %in%
                                                                     c("PTli","PT04","PT03","PT02"),"Lung-MPT", "No-MPT")))))
Fibroblast$Secgroup <- as.factor(Fibroblast$Secgroup)
Idents(Fibroblast) <- "Secgroup"
table(Fibroblast$Secgroup)
save(Fibroblast, file = "Fibroblast_Grouping.Rdata")
###Marker gene bubble chart###
Markers <- c("CYP4B1","MYOC","ADH1B","CFD","SCN7A","AKAP6", "CADPS","DLGAP2","FGF14","PDZRN4","SPP1","SCG2","RAMP1","RARRES1","CPXM1",
              "ZNF560","IGF2BP3", "AURKB", "TDO2","WDR62","COL2A1","COL9A2","COL9A3","S100A1","COL11A2","MYF6","MYF5","FITM1","DES","RAPSN"
)
p<- DotPlot(Fibroblast, features = Markers,dot.scale = 10) +
  scale_colour_gradient(low = "#71B3C2", high = "red") +
  coord_flip() + # Flip the axes
  theme(
    axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9),
    plot.title = element_text(size = 12),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
    panel.spacing = unit(0.1, "lines"),
    strip.spacing = unit(0.05, "lines")
  ) +
  xlab("markers") + # Set the x-axis label
  ylab("CellType") # Set the y-axis label)
print(p)
ggsave("Marker_bubble.pdf", plot = p, width = 8, height = 8)
###Evaluate the expression of CDKN2A in different subgroups###
pdf(file = "Subgroup_CDKN2A.pdf")
VlnPlot(Fibroblast, c("CDKN2A"),pt.size = 0,group.by = "CellType")
dev.off()
###Plot a UMAP showing the expression of a specific gene###
umap_plot <- FeaturePlot(Fibroblast, features = "CDKN2A", reduction = "umap", pt.size = 0.5) + 
  scale_colour_gradient(low = "lightgrey", high = "blue") +
  theme_minimal() +
  labs(x = "UMAP_1", y = "UMAP_2", title = "CDKN2A") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        plot.title = element_text(hjust = 0.5))
print(umap_plot)
ggsave(filename = "UMAP_CDKN2A.pdf", plot = umap_plot, width = 8, height = 6)
