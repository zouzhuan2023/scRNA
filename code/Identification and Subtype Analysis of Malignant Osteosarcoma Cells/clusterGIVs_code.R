#################################################clusterGIVs
library(ClusterGVis)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(devtools)
library(Seurat)
Idents(OB_Tumor_filtered) <- "seurat_clusters"
pbmc.markers.all <- Seurat::FindAllMarkers(OB_Tumor_filtered,
                                           only.pos = TRUE,
                                           min.pct = 0.25,
                                           logfc.threshold = 0.25)
pbmc.markers <- pbmc.markers.all %>%
  dplyr::group_by(cluster) %>%
  dplyr::top_n(n = 20, wt = avg_log2FC)

# check
head(pbmc.markers)
st.data <- prepareDataFromscRNA(object = OB_Tumor_filtered,
                                diffData = pbmc.markers,
                                showAverage = TRUE,
                                assay="RNAnew",
                                keep.uniqGene = FALSE)
str(st.data)
enrich <- enrichCluster(object = st.data,
                        OrgDb = org.Hs.eg.db,
                        type = "BP",
                        organism = "hsa",
                        pvalueCutoff = 0.5,
                        topn = 5,
                        seed = 5201314)
head(enrich)
set.seed(1000)  #5000  1000  4000
markGenes = unique(pbmc.markers$gene)[sample(1:length(unique(pbmc.markers$gene)),50,

visCluster(object = st.data,
           plot.type = "line")
################################################
#go.col <- rep(jjAnno::useMyCol("stallion", n = 20), each = 5)[1:35]  


my_colors <- c(
  "#FF0000", "darkblue","#008080","#800080", "#FF7F50","#0000FF", "#6495ED","#4682B4", "#FF00FF", "#00FFFF",
  "#800000", "#808000", "#008000",  "#808080", 
  "#C0C0C0", "#FFA500", "#A52A2A", "#DEB887", "#5F9EA0", "#7FFF00",
  "#D2691E",   "#DC143C", "#00FA9A"
)

repeated_colors <- rep(my_colors, each = 5)
go.col <- my_colors
go.col <- repeated_colors[1:35] 
sample_colors <- data.frame(
  sample = as.character(0:6), 
  color = c("#EFE2AA", "#AECDE1", "#EBAEA9", "#D2EBC8", "#FCED82", "#8FA4AE", "#F5D2A8")  
)
sample_color_list <- setNames(sample_colors$color, sample_colors$sample)
pdf('clusterGVIs2.pdf',height = 8,width = 10,onefile = F)
visCluster(object = st.data,
           plot.type = "both",
           column_names_rot = 45,
           show_row_dend = F,
           # genes.gp = c('italic',fontsize = 12,col = "black"),
           markGenes = markGenes,
           markGenes.side = "left",
           annoTerm.data = enrich,
           line.side = "left",
           cluster.order = c(1:7),
           go.col = go.col,#rep(jjAnno::useMyCol("stallion",n = 20),each = 5),
           add.bar = F,
           sample.col = sample_color_list)
dev.off()


