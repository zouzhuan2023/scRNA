##prepare the environment####
setwd("/your/own/path")
set.seed(2024)
library(dplyr)
library(spacexr)
library(ggplot2)
library(Seurat)
library(harmony)
library(rjags)
library(AnnoProbe)
library(infercnv)
library(tidyverse)

##sc_data QC####
scRNA <- get(load("patient_sc.rdata"))
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA<- subset(scRNA, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 20)
#sc_data processing
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))
scRNA <- scRNA %>%
  FindNeighbors(reduction = "pca", dims = 1:20) %>% 
  FindClusters(resolution = 0.05)

marker_all <-c("ALPL","RUNX2","IBSP","BGLAP",
               "CD14","MPO","C1QA","C1QB",
               "NCAM1","CD3E","CD4","NKG7",
               "MS4A1","CD19","CD79A","IGHM",
               "PECAM1","CDH5","VEGFA","VWF",
               "FAP","PDPN","ACTA2","COL1A1")
VlnPlot(scRNA,features = marker_all,pt.size = 0, ncol = 8)
cell_types <- c("all","cell","type")
scRNA@meta.data$cell_type <- factor(scRNA@meta.data$seurat_clusters, 
                                    levels = 0:c(length(cell_types)-1), 
                                    labels = cell_types)
type <- scRNA@meta.data[,"cell_type", drop=F]
celltype_counts <- table(type[,1])
#RCTD could not allow cell types with fewer than 25 cells
valid_celltypes <- names(celltype_counts[celltype_counts >= 25])
type <- type[type[,1] %in% valid_celltypes, ,drop=FALSE]
cell_list <- rownames(type)
type$cell_type <- droplevels(type$cell_type)
type <- setNames(type[[1]], rownames(type))
type <- as.factor(type) 
sc_nUMI <- scRNA@meta.data[,"nCount_RNA", drop=F]
sc_nUMI_selected <- sc_nUMI[cell_list, , drop=F]
sc_nUMI <- setNames(as.integer(sc_nUMI[[1]]), rownames(sc_nUMI))
sc_count <- GetAssayData(scRNA, slot = "counts") %>% as.data.frame()
sc_count <- sc_count[,cell_list]
reference <- Reference(sc_count, type, sc_nUMI)
#st_data processing
st_data <- Load10X_Spatial(data.dir = "/your/own/path",
                           filename = "filtered_feature_bc_matrix.h5",
                           assay = "Spatial")
st_data <- RenameCells(st_data, add.cell.id = "patient_ID")
st_data$orig.ident <- "patient_ID"

#QC parameter, for example:
minFeature = 400
minCount = 800
maxMT = 20
minCell = 20
st_data[["percent.mt"]] <- PercentageFeatureSet(st_data, pattern = "^MT-")
counts <- GetAssayData(st_data, slot = "counts", assay = "Spatial")
genes_to_keep <- rowSums(counts > 0) >= minCell
st_data <- st_data[genes_to_keep, ]
st_data <- subset(st_data, subset = nFeature_Spatial >= minFeature & percent.mt <= maxMT & nCount_Spatial >= minCount)
img<- GetTissueCoordinates(st_data)
write.csv(img[,1:2], file = paste0("/your/own/path/", patient,"/coords.csv"))
sparse_matrix <- GetAssayData(object = st_data, layer = "counts")
sparse_matrix <- as.data.frame(sparse_matrix)
write.csv(sparse_matrix, file = paste0("/your/own/path/", patient,"/counts.csv"))
counts <- sparse_matrix
coords <- img[,1:2]
nUMI <- colSums(counts)
puck <- SpatialRNA(coords, counts, nUMI)

##Run RCTD####
myRCTD <- create.RCTD(puck, reference, max_cores = 32)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
saveRDS(myRCTD,file=paste0("/your/own/path/", patient,"/RCTD.rds"))
barcodes <- colnames(myRCTD@spatialRNA@counts)
weights <- myRCTD@results$weights
norm_weights <- normalize_weights(weights)
columns <- colnames(norm_weights)
#Pie chart of cell proportion
m <- as.matrix(norm_weights)
p <- coords
colors <- c("#DF727C","#E7D176","#E7B976","#BC6091","#73C162","#685BA0")
desired_order <- c( "Fibroblasts",
                    "OB",
                    "B_cells",
                    "Macrophages",
                    "NKT_cells",
                    "Endothelial_cells")

col_indices <- match(desired_order, colnames(m))
plt <- vizAllTopics(theta = m,
                    pos = p,
                    topicOrder=col_indices,
                    topicCols=colors,
                    groups = NA,
                    group_cols = NA,
                    r = 80, # size of scatterpies; adjust depending on the coordinates of the pixels
                    lwd = 0,
                    showLegend = TRUE,
                    plotTitle = "scatterpies")
ggsave(paste0("/your/own/path/", patient,"/Spaital_scatterpies.pdf"), width=12, height=6, plot=plt)
write.csv(m,file = paste0("/your/own/path/", patient,"/ratio.csv"))

#Overall cell proportion
csv_files <- list.files(path = "/your/own/path/", pattern = "ratio.csv", full.names = TRUE, recursive = TRUE)
get_folder_name <- function(file_path) {
  basename(dirname(file_path))
}
dataframes <- list()
for (file_path in csv_files) {
  folder_name <- get_folder_name(file_path)
  dataframe <- read.csv(file_path)
  
  if ("X" %in% colnames(dataframe)) {
    dataframe$X <- paste0(folder_name, "_", dataframe$X)
  }
  dataframes[[folder_name]] <- dataframe
}
all_columns <- unique(unlist(lapply(dataframes, colnames)))
dataframes <- lapply(dataframes, function(df) {
  missing_columns <- setdiff(all_columns, colnames(df))
  df[missing_columns] <- 0
  return(df)
})
merge_data <- bind_rows(dataframes)
write_csv(merge_data, file = "merge_ratio.csv")
#Plot a cell proportion bar chart
colors <- c("#E7D176","#73C162","#DF727C","#685BA0","#BC6091","#E7B976")
desired_order <- c("OB", "NKT_cells", "Fibroblasts","Macrophages" , "Endothelial_cells","B_cells" )
data <- read.csv("merge_ratio_processed.csv")
data$Total <- rowSums(data[, -1])
data[, -c(1, ncol(data))] <- data[, -c(1, ncol(data))] / data$Total
data_long <- melt(data, id.vars = "X")
data_long <- data_long[data_long$variable != "Total", ]
data_long$variable <- factor(data_long$variable, levels = desired_order)
#Fig 2B
pdf("barplot.pdf",width = 6,height = 10)
ggplot(data_long, aes(x = X, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Samples", y = "Proportion", title = "Proportional Stacked Bar Chart") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black")) +
  scale_fill_manual(name = "Cell Types", values = colors) +
  coord_flip()
dev.off()




##pseudotime analysis####
#prepare the st_data
all_st_data <- merge(st_data,y=c(other_st_data))
all_st_data <- SCTransform(all_st_data, assay = "Spatial")
all.genes <- rownames(all_st_data)
all_st_data <- ScaleData(all_st_data, features = all.genes)
all_st_data <- RunPCA(all_st_data,npcs = 30, verbose = FALSE) %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE) %>%
  RunUMAP(reduction = "harmony",dims = 1:15)
#Using CytoTracy to determine the starting point
exp <- as.matrix(GetAssayData(all_st_data, assay = "SCT", slot = "data"))
result <- CytoTRACE(exp, ncores = 1)
ori <- all_st_data@meta.data$orig.ident %>% as.character()
emb <- all_st_data@reductions[["umap"]]@cell.embeddings
names(ori) <- colnames(exp)
plotCytoTRACE(cyto_re,phenotype = ori, emb = emb, outputDir = getwd())
#monocle 3
data <- GetAssayData(all_st_data, assay = "SCT", slot = "data")
cell_metadata <- all_st_data@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 10)
plot_pc_variance_explained(cds)

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(data.sample,reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
cds <- cluster_cells(cds)
cds <- learn_graph(cds)
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime")
pseudotime_values <- cds@principal_graph_aux$UMAP$pseudotime %>% as.data.frame()
write.csv(pseudotime_values, "pseudotime_values.csv", row.names = T)
#Fig2C
type <- read.csv("merge_ratio.csv")
pseudotime_values <- read.csv("pseudotime_values.csv")
merged_data <- merge(pseudotime_values, type, by = "barcode")
merged_data_sub <- merged_data["Specific cell type"]
merged_data<- filter(merged_data, is.finite(pseudotime))
range_pseudotime <- range(merged_data_sub$pseudotime)
cat("Range of pseudotime:", range_pseudotime, "\n")
merged_data_sub$bin <- cut(merged_data$pseudotime, breaks = 50, labels = FALSE)
binned_data <- merged_data_sub %>%
  group_by(bin) %>%
  summarize(
    median_celltype = median(OB, na.rm = TRUE),
    lower_quartile = quantile(OB, 0.25, na.rm = TRUE),
    upper_quartile = quantile(OB, 0.75, na.rm = TRUE),
    bin_center = median(bin, na.rm = TRUE)
  )
p <- ggplot(binned_data, aes(x = bin_center, y = median_celltype)) +
  geom_ribbon(aes(ymin = lower_quartile, ymax = upper_quartile), fill = "#E78376", alpha = 0.1) +
  geom_line(color = "black", linetype = "dotted", size = 1,alpha = 0.5) + 
  geom_smooth(method = "loess", color ="#E78376", size = 1, se = FALSE) +  
  theme_minimal() +
  labs(
    x = "Pseudotime",
    y = "Celltype Median",
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)
  )
ggsave(filename = "celltype_Pseudotime.pdf",plot = p,width = 6,height = 3)

#####inferCNV#####
exp <- as.matrix(GetAssayData(all_st_data, assay = "SCT", slot = "data"))
annotation <- data.frame(
  Cell = rownames(st_data@meta.data),
  CellType = st_data@meta.data$orig.ident
)
rownames(annotation) <- annotation[,1]
annotation <- annotation[,-1,drop=FALSE]
geneInfor=annoGene(rownames(exp),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
geneFile='geneFile.txt'
write.table(geneInfor,file = geneFile,sep = '\t',quote = F,col.names = F,row.names = F)

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = exp,
  annotations_file = annotation,
  delim = "\t",
  gene_order_file = 'geneFile.txt',
  ref_group_names = "NLN11"
)

infercnv_obj <- infercnv::run(infercnv_obj,
                              cutoff=0.1, 
                              out_dir="/your/own/path",
                              cluster_by_groups=TRUE,
                              denoise=TRUE,
                              HMM=F)
infercnv::plot_cnv(infercnv_obj,
                   out_dir="/your/own/path",
                   obs_title="Observations (Cells)",
                   ref_title="References (Cells)",
                   cluster_by_groups=TRUE,
                   x.center=1,
                   x.range="auto",
                   hclust_method='ward.D',
                   custom_color_pal = color.palette(c("#8DD3C7","white","#BC80BD"), c(2, 2)),
                   color_safe_pal=FALSE,
                   output_filename="infercnv_pdf",
                   output_format="pdf",
                   dynamic_resize=0)






