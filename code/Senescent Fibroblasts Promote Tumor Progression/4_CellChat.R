###Empty the environment###
rm(list=ls())
###Load the package###
library(data.table)
library(dplyr)
library(tidyverse)
library(patchwork)
library(CellChat)
library(NMF)
library(tidyverse)
library(ggalluvial)
library(Seurat)
library(SeuratObject)
library(data.table)
library(ggsci)
###Set the working path to###
setwd("/path")
load("/path/NK.Rdata")
load("/path/DNT.Rdata")
load("/path/CD8T.Rdata")
load("/path/CD4T.Rdata")
load("/path/sencelltype.Rdata")
###Check the column names of the current meta.data###
colnames(CAF_AURKB@meta.data)
colnames(CD4T@meta.data)
colnames(CD8T@meta.data)
colnames(DNT@meta.data)
colnames(NK@meta.data)
###Delete the specified column###
columns_to_remove <- c("RNAnew_snn_res.0.01", "RNAnew_snn_res.0.02","RNAnew_snn_res.0.03","RNAnew_snn_res.0.04","RNAnew_snn_res.0.05",
                       "RNAnew_snn_res.0.06","RNAnew_snn_res.0.07","RNAnew_snn_res.0.08","RNAnew_snn_res.0.09","RNAnew_snn_res.0.1","seurat_clusters",
                       "Secgroup")
CAF_AURKB@meta.data <- CAF_AURKB@meta.data[, !colnames(CAF_AURKB@meta.data) %in% columns_to_remove]

####Merge###
seurat_obj <- merge(CAF_AURKB, y = list(CD4T,CD8T,DNT,NK), add.cell.ids = c("CAF_AURKB","CD4T","CD8T","DNT","NK"), project = "seurat_obj")
seurat_obj <- RenameAssays(seurat_obj, `RNAnew` = "RNA")
save(seurat_obj, file="Merge_NKT_senCAF.Rdata")
load("/path/Merge_NKT_senCAF.Rdata")
###Create a cellchat object###
scRNA_harmony<-seurat_obj
###Extract the expression matrix and cell classification information###
data.input <- GetAssayData(scRNA_harmony, assay = "RNA", slot = "data")
identity <- subset(scRNA_harmony@meta.data, select = "CellType")
cellchat <- createCellChat(object = data.input, meta = identity,  group.by = "CellType")
####CellChatDB.human###
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
unique(CellChatDB$interaction$annotation)
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor")
cellchat@DB <- CellChatDB.use # set the used database in the object

###Preprocess the expression data###
##Subset the expression data of signaling genes to save computational cost###
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 10)
###Identify overexpressed genes###
cellchat <- identifyOverExpressedGenes(cellchat)
###Identify ligand-receptor pairs###
cellchat <- identifyOverExpressedInteractions(cellchat)
###Map ligands and receptors onto the PPI (protein-protein interaction) network###
cellchat <- projectData(cellchat, PPI.human)
###Calculate communication probabilities to infer the cell interaction communication network###
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
df.net.1 <- subsetCommunication(cellchat,slot.name = "netP")
df.net.2 <- subsetCommunication(cellchat )
write.csv(df.net.2,file="netP.csv")
###Infer intercellular communication at the signaling pathway level###
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, "cellchat.rds")
###Visualize the results###
dev.off()
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
###Examine pathways###
levels(cellchat@idents)
vertex.receiver = c(4, 6)
cellchat@netP$pathways
pathways.show<-df.net.2$pathway_name
###Specify the pathways to be displayed###
pathways.show <- "COLLAGEN"
##Circos plot###
pdf("COLLAGEN_circle.pdf")
netVisual_aggregate(cellchat, signaling = "COLLAGEN", layout = "circle")
dev.off()
##Heatmap###
pdf("COLLAGEN_Heatmap.pdf")
netVisual_heatmap(cellchat, signaling = "COLLAGEN", color.heatmap = "Reds")
dev.off()
###Calculate the contribution of each ligand-receptor pair to the signaling pathway###
pdf("COLLAGEN_contribution .pdf")
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()
###Bubble chart###
levels(cellchat@idents)
pdf("Bubble.pdf")
netVisual_bubble(cellchat, sources.use = 6, targets.use = c(1:4), remove.isolate = FALSE)
dev.off()
###Calculate network centrality scores###
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
pdf("centrality_scores.pdf")
netAnalysis_signalingRole_network(cellchat, signaling = "COLLAGEN", width = 8, height = 4, font.size = 10)
dev.off()