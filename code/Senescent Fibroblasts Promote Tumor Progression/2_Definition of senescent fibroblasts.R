###Empty the environment###
rm(list=ls())
###Load the package###
library(Seurat)
library(RColorBrewer) 
library(viridis)
library(wesanderson)
library(ggplot2)
###Set the working path to###
setwd("/path")
load("/path/Fibroblast_Grouping.Rdata")
table(Fibroblast@meta.data$CellType)
Idents(Fibroblast)<-"CellType"
CAF_AURKB<-subset(Fibroblast,ident="CAF_AURKB")
save(CAF_AURKB,file="CAF_AURKB.Rdata")
load("CAF_AURKB.Rdata")
###Determine whether each cell expresses a specific gene###
express_gene <- FetchData(object = CAF_AURKB, vars = "CDKN2A") > 0
###Convert logical values to the desired categories###
cell_type <- ifelse(express_gene, "senCAF", "nonsenCAF")
###Ensure that the length of the generated vector matches the number of cells in the Seurat object###
if (length(cell_type) == ncol(CAF_AURKB)) {
  CAF_AURKB <- AddMetaData(object = CAF_AURKB, metadata = cell_type, col.name = "CellType")
} else {
  stop("Length of the metadata does not match the number of cells in the Seurat object.")
}
p<-DimPlot(CAF_AURKB, reduction = "umap", group.by = "CellType")
ggsave("senCAF_clustering.pdf",p,width=8,height = 6)
table(CAF_AURKB$CellType)
save(CAF_AURKB,file = "sencelltype.Rdata")
load("/path/sencelltype.Rdata")
scRNAident<-CAF_AURKB
markers<- c("CDKN2A","CDKN2B","GLB1","COL1A1","COL12A1","MMP10",
            "MMP13","CXCL9","TNFSF11","CD276","TGFB1","CDKN1A",
            "CDKN1B","TP53","RB1","FOXO3","HMGB1","TNFRSF10D","LMNB1",
            "SPP1","VEGFA"
)
pal <- wes_palette("Zissou1", 10, type = "continuous")
Idents(scRNAident) <- "CellType"
p <- DotPlot(scRNAident, features = markers, dot.scale = 8) + RotatedAxis() + 
  theme(axis.text.x = element_text(angle = 45, face="italic", hjust=1), axis.text.y = element_text(face="bold")) + 
  scale_colour_gradientn(colours = pal)+ theme(legend.position="top")  + labs(title = "", y = "", x="")
ggsave("Senescence-associated_bubble.pdf", p, width = 10, height = 6)
