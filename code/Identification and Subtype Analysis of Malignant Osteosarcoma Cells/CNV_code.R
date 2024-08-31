
#################################OB CNV
library(rjags)
library(Seurat)
library(AnnoProbe)
library(infercnv)
library(tidyverse)
library(locfit)
library(edgeR)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(scRNAtoolVis)
library(dplyr)
library(tidydr)
library(stringr)
library(viridis)
library(scCustomize)
library(patchwork)
library(sscVis)
library(RColorBrewer)


load("OS.combinedmerge.Rdata")

####
counts <- GetAssayData(OB.combined, slot = 'counts')

anno <- data.frame(Idents(OB.combined))

geneInfor=annoGene(rownames(counts),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)




infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts,
                                    annotations_file = anno,
                                    delim="\t",
                                    gene_order_file = geneFile,
                                    #min_max_counts_per_cell = c(100, +Inf),
                                    ref_group_names = c("NK/T_cell"))


infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir="path",
                             cluster_by_groups=TRUE,  
                             denoise=TRUE, 
                             HMM=F,
                             num_threads=20)#


save(infercnv_obj, file = "infercnv_obj.Rdata")

###############save PDF
print(infercnv_obj)

infercnv::plot_cnv(infercnv_obj,
                   out_dir="path/pdf",
                   obs_title="Observations (Cells)",
                   ref_title="References (Cells)",
                   cluster_by_groups=TRUE,
                   x.center=1,
                   x.range="auto",
                   hclust_method='ward.D',
                   #custom_color_pal = color.palette(c("#0071B2", "white", "#C3250A"), c(2, 2)),
                   custom_color_pal = color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2)),
                   color_safe_pal=FALSE,
                   output_filename="infercnv_pdf",
                   output_format="pdf",
                   #png_res=300,
                   dynamic_resize=0)

##########################################
#########################################CNV score
CNVthreshold <- function(ref){
  cnv_ref <- colMeans((ref-1)^2)
  score_threshold <- mean(cnv_ref)+2*sd(cnv_ref)
  return(score_threshold)
}
estimateCNV <- function(obs, ref, score_threshold, cor_threshold){
  cnv_obs <- colMeans((obs-1)^2)
  cnv_ref <- colMeans((ref-1)^2)
  
  length <- length(cnv_obs)*0.05
  cell_top <- names(sort(cnv_obs, decreasing=T))[1:length]
  cnv_top <- rowMeans(obs[, cell_top])
  
  cor_obs <- apply(obs, 2, function(x)cor(x, cnv_top))
  cor_ref <- apply(ref, 2, function(x)cor(x, cnv_top))
  
  cnv <- data.frame(score=c(cnv_obs, cnv_ref), cor=c(cor_obs, cor_ref), barcode=c(colnames(obs), colnames(ref)))
  
  cnv$type <- 'Other'
  cnv$type[cnv$score>score_threshold & cnv$cor>cor_threshold] <- 'Malignant'
  cnv$type[cnv$score<score_threshold & cnv$cor<cor_threshold] <- 'Not Malignant'
  return(cnv)
}

obs <- read.table("infercnv.observations.txt", header=T, check.names=F)
ref <- read.table("infercnv.references.txt", header=T, check.names=F)
cnv_ref <- colMeans((ref-1)^2)
score_threshold <- CNVthreshold(ref)
result <- estimateCNV(obs, ref, 0.003382432, 0.4)




cellAnnota <- anno
cellAnnota$barcode <- rownames(cellAnnota)
result <- merge(result, cellAnnota, by = "barcode")
write.table(result,file = "CNV_result.txt", row.names = F, sep = '\t')


OB.combined@meta.data$barcode <- rownames(OB.combined@meta.data)
OB.combined@meta.data <- merge(OB.combined@meta.data, result[c(1,4)], by = "barcode")
rownames(OB.combined@meta.data) <- OB.combined@meta.data$barcode
OB.combined@meta.data <- rename(OB.combined@meta.data, CNVtype = type)
view(OB.combined@meta.data)
write.csv(OB.combined@meta.data, file = "OB_combined_metadata11.csv", row.names = TRUE)

save(OB.combined,file="OB.combinedresult.Rdata")



########Umap
load(file = "OB.combinedresult.Rdata")
all.genes <- rownames(OB.combined)
OB.combined <- ScaleData(OB.combined, features = all.genes)
OB.combined <- FindVariableFeatures(OB.combined, selection.method = "vst", nfeatures = 5000)
OB.combined <- RunPCA(OB.combined, features = VariableFeatures(object = OB.combined))
OB.combined <- FindNeighbors(OB.combined, dims = 1:10)
OB.combined <- RunUMAP(OB.combined, dims = 1:10)
OB.combined <- RunTSNE(OB.combined, dims = 1:10)

##################################3

#####
output <- read.table(file = "CNV_result.txt",header = T)

output <- output[complete.cases(output), ]

if (nrow(output) == 0) {
  stop("Data frame is empty after removing rows with missing values")
}


change <- ifelse((output$score > 0.003382432 & output$cor > 0.4), "Malignant", ifelse((output$score < 0.003382432 & output$cor < 0.4),"Not Malignant","Other"))  

library(dplyr)
library(limma)
library(ggplot2)
p <- ggplot(output, aes(x = cor, y = score)) +
  geom_point(alpha=0.8, size=1.5) +
  geom_point(aes(col=change)) +
  scale_color_manual(values=c("#EBAEA9", "#AECDE1", "grey")) +
  geom_vline(xintercept=0.4, lty=3, col="black", lwd=0.8) +
  geom_hline(yintercept = 0.003382432, lty=3, col="black", lwd=0.8) +
  labs(x="cnv_correlation", y="cnv_score") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(plot.title = element_text(hjust = 0.5, size=24), 
        legend.position="right", 
        legend.title = element_blank(),
        legend.text=element_text(size=18),
        legend.key.size = unit(1, 'cm'),
        axis.title.x = element_text(size=18), 
        axis.title.y = element_text(size=18),
        axis.text = element_text(size=14, face = "bold"))

ggsave(filename = "scatterplot.pdf", plot = p, width = 8, height = 5)



OB <- grep("^Tcell", rownames(OB.combined@meta.data), value = TRUE, invert = TRUE)

View(OB.combined@meta.data)
OB.combined <- subset(OB.combined, cells = OB)
table(OB.combined$CNVtype)
table(OB.combined$seurat_clusters)
table(OB.combined$celltype1)
table(OB.combined$batch)
Idents(OB.combined) <- "CNVtype"
Malig <- subset(OB.combined, subset = CNVtype == "Malignant")
table(Malig$CNVtype)
table(Malig$batch)
table(Malig$celltype1)
View(Malig@meta.data)

########################
OB_Tumor <- Malig
OB_Tumor <- FindVariableFeatures(OB_Tumor, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(OB_Tumor)
OB_Tumor <- ScaleData(OB_Tumor, features = all.genes)
OB_Tumor <- RunPCA(OB_Tumor, features = VariableFeatures(object = OB_Tumor))
OB_Tumor <- FindNeighbors(OB_Tumor, dims = 1:10)
#Myeloid <- FindClusters(Myeloid, resolution = 0.1)
OB_Tumor <- RunUMAP(OB_Tumor, dims = 1:10)
OB_Tumor <- RunTSNE(OB_Tumor, dims = 1:10)
save(OB_Tumor,file="OB_Tumor_PCA.Rdata")
load(OB_Tumor,file="OB_Tumor_PCA.Rdata")

OB_Tumor <-FindClusters(OB_Tumor,resolution = 0.11)
dir.create("011")
table(OB_Tumor$seurat_clusters)
table(OB_Tumor$batch)
table(OB_Tumor$Secgroup)

Idents(OB_Tumor) <- "seurat_clusters"
logFCfilter =0.25
adjPvalFilter=0.05
OB_Tumor.markers <- FindAllMarkers(object = OB_Tumor,
                                   only.pos = FALSE,
                                   min.pct = 0.25,
                                   logfc.threshold = logFCfilter)
sig.markers=OB_Tumor.markers[(abs(as.numeric(as.vector(OB_Tumor.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(OB_Tumor.markers$p_val_adj))<adjPvalFilter),]
write.csv(sig.markers,file="011/clusterMarkers.csv",sep="\t",row.names=F,quote=F)
save(OB_Tumor,file="011/OB_Tumor011.Rdata")
table(OB_Tumor$Secgroup)
table(OB_Tumor$CNVtype)
Idents(OB_Tumor) <- "Secgroup"
OB_Tumor$Secgroup <- OB_Tumor@active.ident
table(OB_Tumor$batch)
table(OB_Tumor$Secgroup)
Idents(OB_Tumor) <- "Secgroup"
OB_Tumor$Secgroup <- OB_Tumor@active.ident
table(OB_Tumor$Secgroup)
OB_Tumor_filtered <- subset(OB_Tumor, subset = Secgroup != "NT")
table(OB_Tumor_filtered$Secgroup)
Idents(OB_Tumor_filtered) <- "Secgroup"
OB_Tumor_filtered$Secgroup <- OB_Tumor_filtered@active.ident
save(OB_Tumor_filtered,file="011/OB_Tumor_filtered011.Rdata")


#############################################
view(OB_Tumor_filtered@meta.data)

table(OB_Tumor_filtered$seurat_clusters)
###########################################################3
Idents(OB_Tumor_filtered) <-"seurat_clusters"
p <- DimPlot(OB_Tumor_filtered, reduction = "tsne",
             group.by ="seurat_clusters"
             ,label = T,raster=FALSE)+
  scale_color_manual(values = c("#E78376", 
                                "#81CAD9",
                                "#68B9AA",
                                "#7888AA",
                                "#F0B7A6",
                                "#ABB1C9",
                                "#B1DBD1"))+
  theme(
    panel.border = element_blank(),  
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    axis.line.x.bottom = element_line(color = "black", size = 1),  
    axis.line.y.left = element_line(color = "black", size = 1)  
  ) +
  labs(title = "Clusters")

ggsave(filename = "fig1 Cluster_tsne.pdf",p,height = 7,width = 8)


#########################################
cell.group <- OB_Tumor_filtered
Idents(cell.group) <-"Secgroup"  
cell.group@meta.data$Secgroup <- cell.group@active.ident
cellnum <- table(cell.group$Secgroup,cell.group$seurat_clusters)
cell.prop<-as.data.frame(prop.table(cellnum,1))  
colnames(cell.prop)<-c("Group","Subcelltype","Proportion")  

colors <- c("#E78376", 
            "#81CAD9",
            "#68B9AA",
            "#7888AA",
            "#F0B7A6",
            "#ABB1C9",
            "#B1DBD1")


p <- ggplot(cell.prop, aes(x = Group, y = Proportion, fill = Subcelltype)) +
  geom_bar(stat = "identity", position = "stack", color = "white", size = 0.11) +
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    panel.border = element_blank(),  
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1),
    axis.text.y = element_text(size = rel(1)),
    axis.line.x.bottom = element_line(color = "black", size = 1), 
    axis.line.y.left = element_line(color = "black", size = 1)     
  )
ggsave("011/scalegraph.pdf",width = 6, height = 5)
dev.off()
