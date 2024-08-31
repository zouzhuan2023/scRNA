##prepare the environment####
setwd("/your/own/path")
library(Seurat)
library(dplyr)
library(tidyverse)
library(patchwork)
library(magrittr)
library(ggsci)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)

scRNA<- get(load("st_adata.Rdata"))

#Sort the cell counts in descending order
table(scRNA@meta.data$celltype)
cell_counts <- table(scRNA@meta.data$celltype)
sorted_cell_counts <- sort(cell_counts, decreasing = TRUE)
sorted_cell_types <- names(sorted_cell_counts)
scRNA@meta.data$celltype <- factor(scRNA@meta.data$celltype, levels = sorted_cell_types)
##Plot a UMAP graph####
Idents(scRNA)="celltype"
p=DimPlot(scRNA, reduction = "umap",
          group.by ="celltype"
          ,label = F,raster=FALSE)+
  scale_color_manual(values = c(
    my36colors <-c("#E78376", 
                   "#81CAD9",
                   "#68B9AA",
                   "#7888AA",
                   "#F0B7A6",
                   "#ABB1C9",
                   "#B1DBD1",
                   "#fad1e0"
    )))+
  theme(panel.border = element_blank(), 
        
        axis.ticks = element_blank())+ 
  labs(title = "Celltype")

ggsave(filename = "UMAP.pdf",p,height =7,width =10)


##Plot a cell proportion bar chart####
Idents(scRNA) <- "batch"
table(scRNA$batch)
scRNA@meta.data$Secgroup <- ifelse(scRNA@meta.data$batch %in% c("MLN07", "MLN09"), "MLN",#lymphatic node transfer(lymphatic node)

                                   ifelse(scRNA@meta.data$batch %in% c("MLN08", "MLN10", "NLN09", "NLN11"), "NLN", #Lymph nodes without metastases

                                          ifelse(scRNA@meta.data$batch %in% c("NT01", "NT03", "NT07", "NT08", "NTfan", "NTli"), "NT",#paraneoplastic bone tissue

                                                 ifelse(scRNA@meta.data$batch %in% c("PT07", "PTyi"), "Lymph-MPT",  #Osteosarcoma with lymph node metastasis

                                                        ifelse(scRNA@meta.data$batch %in% c("PTli", "PT04", "PT03", "PT02"), "Lung-MPT", "No-MPT"))))) #Osteosarcoma with lung metastases and Osteosarcoma tissue without metastases     

scRNA$Secgroup <- as.factor(scRNA$Secgroup)
Idents(scRNA) <- "Secgroup"
table(scRNA$Secgroup)

cell.group <-scRNA
Idents(cell.group) <-"Secgroup" 
cell.group@meta.data$Secgroup <- cell.group@active.ident
cellnum <- table(cell.group$Secgroup,cell.group$celltype)
cell.prop<-as.data.frame(prop.table(cellnum,1)) ##Proportioning cellnum with the prop.table function
#prop.table(data): convert data to percentage; prop.table(data,1): data by row for percentage; prop.table(data,2): data by column for percentage
View(cell.prop)
colnames(cell.prop)<-c("Group","Celltype","Proportion")  #Rename the three rows of variables in cell.prop

table(cell.group$Secgroup)
new_order <- c("No-MPT","Lung-MPT","MLN","Lymph-MPT","NLN","NT")
cell.prop$Group <- factor(cell.prop$Group,levels = new_order)

colors <- c("#E78376", 
            "#81CAD9",
            "#68B9AA",
            "#7888AA",
            "#F0B7A6",
            "#ABB1C9",
            "#B1DBD1",
            "#fad1e0"
)



p <- ggplot(cell.prop, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "stack", color = "white", size = 0.11) + 
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1),
    axis.text.y = element_text(size = rel(1))
  )
p
ggsave("proportion bar chart.pdf",p,width = 6, height = 5)
dev.off()


##Plot a bubble chart####
final.markers <- c(
  'CD14', 'CD68','C1QA','LYZ',
  'CD3D','CD3G','NKG7','GZMK','CD8A',
  'RUNX2', 'ALPL','PAGE2B','IBSP',
  'CD79A','MS4A1','CD79B','IGHM',
  'COL6A3','FBLN1','DCN','COL5A2', 
  'HSPG2','EGFL7','CLDN5','FLT1',
  'MCAM','RGS5','SOX5', 'PDGFRB',
  'ACAN','SOX9','PTH1R','COL2A1')



# Convert to a factor
scRNA$celltype <- as.factor(scRNA$celltype)

levels(scRNA$celltype)[levels(scRNA$celltype) == "NK/T  cell"] <- "NKT cell"


colors <- c("#E78376", 
            "#81CAD9",
            "#68B9AA",
            "#7888AA",
            "#F0B7A6",
            "#ABB1C9",
            "#B1DBD1",
            "#fad1e0")

Idents(scRNA) <- "celltype"
table(scRNA@meta.data$celltype)

pal <- c("white", "red")

p <- DotPlot(scRNA, features = final.markers, dot.scale = 10) +
  RotatedAxis() +
  theme(
    axis.text.x = element_text(angle = 45, face = "italic", hjust = 1),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right",
    panel.border = element_rect(color = "black", size = 1.5, linetype = "solid"),
    axis.title = element_blank()
  ) +
  scale_colour_gradientn(colours = pal) +
  labs(title = "celltype markers", y = "", x = "") +
  geom_point(shape = 21, stroke = 0.5, color = "black", aes(fill = avg.exp.scaled, size = pct.exp)) + 
  scale_size_continuous(range = c(1, 10)) + 
  scale_fill_gradient(low = "white", high = "red") +
  guides(fill = guide_colourbar(title = "Average Expression"), size = guide_legend(title = "Percent Expressed"))

ggsave(filename = "bubble chart.pdf", plot = p, width = 14, height = 8)
dev.off()


##Plot a pie chart####
install.packages("ggplot2")
library(ggplot2)
table(scRNA@meta.data$celltype)

celltype_table <- table(scRNA@meta.data$celltype)
data<- as.data.frame(celltype_table)

colors <-c("#E78376", 
           "#81CAD9",
           "#68B9AA",
           "#7888AA",
           "#F0B7A6",
           "#ABB1C9",
           "#B1DBD1",
           "#fad1e0")

# Calculate the percentage
data$percentage <- data$count / sum(data$count) * 100

ggplot(data, aes(x = "", y = percentage, fill = celltype)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(aes(label = paste0(round(percentage, 1), "%")), 
            position = position_stack(vjust = 0.5)) +
  labs(fill = "Cell Type", title = "Cell Type Distribution")+
  scale_fill_manual(values = colors)
ggsave("pie chart.pdf")



##Plot the proportion chart of cell distribution among the 28 clusters####
Idents(scRNA) <- "batch"
table(scRNA$batch)

scRNA$batch<- as.factor(scRNA$batch)
Idents(scRNA) <- "batch"
table(scRNA$batch)

cell.group <-scRNA
Idents(cell.group) <-"batch"  
cell.group@meta.data$batch <- cell.group@active.ident
cellnum <- table(cell.group$batch,cell.group$celltype)
cell.prop<-as.data.frame(prop.table(cellnum,1)) #Proportioning cellnum with the prop.table function
#prop.table(data): convert data to percentage; prop.table(data,1): data by row for percentage; prop.table(data,2): data by column for percentage
View(cell.prop)
# Calculate the proportion of each cell type within each group
cell.prop <- cell.group@meta.data %>%
  group_by(batch, celltype) %>%
  summarise(count = n()) %>%
  mutate(Proportion = count / sum(count)) %>%
  ungroup()

colnames(cell.prop) <- c("Group", "Celltype", "Count", "Proportion")


colors <- c("#E78376", 
            "#81CAD9",
            "#68B9AA",
            "#7888AA",
            "#F0B7A6",
            "#ABB1C9",
            "#B1DBD1",
            "#fad1e0")

p <- ggplot(cell.prop, aes(x = Group, y = Proportion, fill = Celltype)) +
  geom_bar(stat = "identity", position = "stack", color = "white", size = 0.11) + 
  scale_fill_manual(values = colors) +
  theme_minimal() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"),
    plot.margin = unit(c(1, 1, 1, 1), "lines"),
    axis.text.x = element_text(size = rel(1), angle = 45, hjust = 1),
    axis.text.y = element_text(size = rel(1))
  )
p
CairoPDF("proportion chart of cell distribution .pdf",p,width = 16, height = 8)
dev.off()
