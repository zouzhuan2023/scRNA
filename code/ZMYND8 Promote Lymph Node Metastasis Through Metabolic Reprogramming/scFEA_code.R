
library(magrittr)
library(RColorBrewer)
library(tidyverse)
library(stringr)
library(patchwork)
library(Seurat)
library(dplyr)
library(readxl)
table(OB_Tumor_filtered$seurat_clusters)
Idents(OB_Tumor_filtered) <-"seurat_clusters"


##########R
count_matrix <- OB_Tumor_filtered@assays$RNA@counts
write.csv(count_matrix, file='OB_Tumor_filtered_count_matrix.csv', row.names = T)


######################################python
cd /path/
##############
python src/scFEA.py --data_dir data \
--test_file OB_Tumor_filtered_count_matrix.csv \
--moduleGene_file module_gene_m168.csv \
--cName_file cName_c70_m168.csv \
--stoichiometry_matrix cmMat_c70_m168.csv \
--res_dir out \
--sc_imputation True
##########################################End
adj_flux <- read.csv("OB_Tumor_filtered_count_matrix.csv",row.names = 1)


OB_Tumor_filtered$group_cells <- paste0(OB_Tumor_filtered$Secgroup,"_",OB_Tumor_filtered$seurat_clusters)
table(OB_Tumor_filtered$group_cells)
save(OB_Tumor_filtered,file = "OB_Tumor_filtered_group_cells.Rdata")
cell_anno <- data.frame(cellid=rownames(OB_Tumor_filtered@meta.data),
                        
                        group_cells=OB_Tumor_filtered$seurat_clusters)

cell_anno <- cell_anno[order(cell_anno$group_cells),] 
rownames(adj_flux) <- gsub("\\.", "-", rownames(adj_flux))
adj_flux <- adj_flux[cell_anno$cellid,] 

df_averages <- adj_flux %>%
  group_by(group = cell_anno$group_cells) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  dplyr::select(-group)

rownames(df_averages) <- unique(cell_anno$group_cells)
df_averages <- t(df_averages)%>% as.data.frame()

table(OB_Tumor_filtered$seurat_clusters)
df_flux <- df_averages[,c("0",
                          "1","2","3","4","5","6")]
colnames(df_flux) <- c("0",
                       "1","2","3","4","5","6")

df_flux = df_flux[apply(df_flux, 1, function(x) sd(x)!=0),]

library(ComplexHeatmap)
df_flux[is.na(df_flux)] <- 0


annotation_col = data.frame(
  #Secgroup = c("Lymph-MPT", "MLN"),  
  seurat_clusters = c("0",
                      "1","2","3","4","5","6") 
)
row.names(annotation_col) <- colnames(df_flux)


human_moduleInfo <- read.csv("scFEA.human.moduleinfo.csv", header = T, row.names = 1)

annotation_row = human_moduleInfo[rownames(df_flux),]
annotation_row  = as.data.frame(annotation_row[,c("SM_anno")])
rownames(annotation_row) = rownames(df_flux)
colnames(annotation_row) = c("SM_anno")


cellcolor <- c("#E78376" , 
               "#81CAD9",
               "#68B9AA",
               "#7888AA",
               "#F0B7A6",
               "#ABB1C9",
               "#B1DBD1")
table(OB_Tumor_filtered$group_cells)
names(cellcolor) <- c("0",
                      "1","2","3","4","5","6")

table(OB_Tumor_filtered$Secgroup)
groupcolor <- c("#1F77B4",
                "#7ABFE2",
                "#E377C2",
                "#fad1e0",
                "#FFBB78")

modulecolor <- c("#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#FFFF00",
                 "#808000","#FF00FF","#FA8072","#800080","#87CEEB","#40E0D0","#5F9EA0",
                 "#008B8B","#FFE4B5","#228B22","#4682B4","#32CD32","#F0E68C","#FFFFE0",
                 "#FF6347")

names(modulecolor) <- unique(annotation_row$SM_anno)

ann_colors <- list(seurat_clusters=cellcolor, SM_anno=modulecolor) 

#################################
adj_flux <- read.csv("OB_Tumor_filtered_count_matrix.csv",row.names = 1)

colnames(adj_flux) <- gsub("_", "-", colnames(adj_flux))
rownames(adj_flux) <- gsub("\\.", "-", rownames(adj_flux))
predFlux <- t(data.matrix(adj_flux))
OB_Tumor_filtered[["FLUX"]] <- CreateAssayObject(counts = predFlux)

DefaultAssay(OB_Tumor_filtered) <- 'FLUX'
library(ggplot2)
load("OB_Tumor_filtered_group_cells.Rdata")


##########################################
table(OB_Tumor_filtered$seurat_clusters)
Idents(OB_Tumor_filtered)=OB_Tumor_filtered$Secgroup
table(OB_Tumor_filtered$Secgroup)
Idents(OB_Tumor_filtered) <- "seurat_clusters"
table(OB_Tumor_filtered$Secgroup)
Module <- c("M-1","M-2","M-3","M-4","M-5","M-6","M-7","M-8","M-9","M-10","M-11","M-12","M-13","M-14",
            "M-60","M-53","M-21","M-47","M-148")

human_moduleInfo <- read.csv("scFEA.human.moduleinfo.csv", header = T, row.names = 1)
OB_Tumor_filtered_module <- human_moduleInfo[c("M_1","M_2","M_3","M_4","M_5","M_6","M_7","M_8","M_9","M_10","M_11","M_12","M_13","M_14",
                                               "M_60","M_53","M_21","M_47","M_148"),]
OB_Tumor_filtered_module$M_name <- paste0(OB_Tumor_filtered_module$Module_id, ": ",OB_Tumor_filtered_module$Compound_IN_name,"_",OB_Tumor_filtered_module$Compound_OUT_name)
DotPlot(OB_Tumor_filtered, features = Module,cols = c("Spectral")) + 
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"))) + 
  RotatedAxis() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.title = element_blank())+
  scale_x_discrete(labels=OB_Tumor_filtered_module$M_name)+
  scale_y_discrete(labels=c("0","1","2","3","4","5","6"))+
  coord_flip()
ggsave(filename = "DotPlot.pdf",height = 7,width = 8)
################################################################
Idents(OB_Tumor_filtered)="seurat_clusters"
table(OB_Tumor_filtered$seurat_clusters)


######################(PT与NT）###########################
data_1 <- as.data.frame(OB_Tumor_filtered$RNAnew@data)
data_1=as.data.frame(t(data_1))

genelist=read.table(file="geneset.txt",sep = "",header = F)
genelist=genelist$V1

common_genes <- intersect(genelist, colnames(data_1))

data_2 <- data_1[, common_genes]

META<-OB_Tumor_filtered@meta.data
CELLNAMES<-intersect(row.names(META),rownames(data_2))
META<-META[CELLNAMES,]
data_2<-data_2[CELLNAMES,]
META$ID<-row.names(META)
META1<-META[,c("ID","seurat_clusters")]
data<-cbind(data_2,META1)
data2=data[,-(ncol(data)-1)]

data3=melt(data2,id.vars=c("seurat_clusters"))
colnames(data3)=c("Type","Gene","Expression")
head(data)

library(ggpubr)
p=ggboxplot(data3, x="Gene", y="Expression", color = "Type", 
            ylab="Gene expression",
            xlab="",
            legend.title="Type",
            palette = c(  "#1F77B4",
                          "green1",
                          "#E377C2",
                          '#B53E2B',
                          "darkblue",
                          "red",
                          'gold'),
            width=0.6, add = "none")
p=p+rotate_x_text(60)
p
p1 = p + stat_compare_means(aes(group = Type), method = "kruskal.test", label = "p.signif") +
  stat_compare_means(aes(group = Type), method = "dunn.test", label = "p.signif",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")))
p1
ggsave(filename = "gene_ggboxplot.pdf",p1,height = 5,width = 8)