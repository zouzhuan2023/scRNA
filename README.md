# scRNA and ST
The ***scRNA and ST*** used for building **Single-cell and spatial transcriptomic analyses reveals the metastasis mechanisms of osteosarcoma**.

## Table of Contents
1. System requirements
2. Installation guide
3. Codes
4. Dataset
5. Contact

## System requirements
The following are the version numbers of the software or algorithms used in this study.

	AUCell 1.12.0
	GSVA 1.38.2
        SEDR  1.0.0
	spacexr 2.2.1
        CellChat 1.6.1
        TESLA  1.2.4
        stlearn 0.4.12
	Seurat 4.3.0.1
	scFEA 1.1
	scanpy 1.10.1
        stlearn 0.4.12
	monocle3 1.3.1
	Python 3.9
 	Ubuntu 18.04
	R 4.0.5, 4.1.0 
 
## Installation guide

Python libraries can be installed in a shell environment using the "pip install" command. 

	pip install "library_name"

R packages can be installed in the R environment using the "install.packege()" or "BiocManager::install()" commands.

	install.packege("packege_name")

	if(!"BiocManager"%in%installed.packages()){ 
	install.packages("BiocManager")}
 	if(!"packege_name"%in%installed.packages()){ 
	BiocManager::install("packege_name")}

	if (!"devtools" %in% installed.packages()) {
  	install.packages("devtools")}
   	devtools::install_github("packege_name")



## Codes
Specific descriptions of the codes can be found in the corresponding documents.

### 1. Quality control and subgroup
"1_Data_preprocessing_code.pyt" is used to import the processed data into the python environment and perform quality control to eliminate low quality cells and genes and finally perform dimensionality reduction clustering on the data.

"2_Single cell huge group analysis.R" was used to define overall cellular subpopulations and plot subpopulation umap plots, cell scale plots, and other plots.

### 2. Spatial Transcriptional Landscape of Osteosarcoma and Metastatic Lymph Nodes

"Spatial Transcriptional Landscape of Osteosarcoma and Metastatic Lymph Nodes.R" was used to construct Overall landscape of spatial transcriptomics for 8 samples.

### 3.Identification and Subtype Analysis of Malignant Osteosarcoma Cells

"clusterGIVs_code.R" was used to analyze the functions of tumor cell subpopulations.

"CNV_code.R" was used to analyze the copy number variations in osteoblasts.

"GSVA_code.R" was used to analyze the functions of tumor cell subpopulations.

"Score_code.R" was used to evaluate the enrichment scores of tumor cell subpopulations across different gene sets.

### 4. NPW(+) Osteosarcoma Cells Promote Lung via Metastasis Macrophages Recruitment 

"Gene_expression.R" was used to analyze the expression of specific genes in tumor cell subpopulations.

"SpatialFeaturePlot.R" was used to analyze the spatial localization of genes.

"stLearn_code.py" was used to analyze ligand-receptor pairs in spatial contexts.

"TESLA_code.py" was used to analyze the spatial localization of genes.

### 5. ZMYND8 Promote Lymph Node Metastasis Through Metabolic Reprogramming
“scFEA_code.R” was used to analyze metabolic pathway activity in tumor cell subpopulations
“scMetabolism_Code.R” was used to analyze metabolic pathway activity of tumor cell subsets
“SEDR_Code.py” was used to analyze the spatial distribution of expression of specific genes.
### 6. SLPI(+) Tumor Cells Inhibit the Function of NK Cell and Promote Survival In Lymph Nodes.
“CIBERSORT_code.R” was used to estimate the abundance of immune cells in the sample.

### 7  Senescent Fibroblasts Promote Tumor Progression
“1_Analysis of fibroblasts.R” was used for the annotation of fibroblasts and the expression of CDKN2A in each subpopulation.
“2_Definition of senescent fibroblasts.R” was used for defining senescent fibroblasts.
“3_Aging gene set score.R” was used for scoring of 21 Senescence-Associated Genes in senCAF and nonsenCAF.
“4_CellChat.R” was used for the analysis of intercellular communication between senCAF Cells and tumor cell subgroups.
“5_SEDR.py” was used for predicting the spatial localization of CDKN1A, IBSP, and CDKN2A.
“6_Colocalization of the CDKN1A with fibroblasts.R” was used for the analysis of spatial localization of fibroblasts and the CDKN1A.

## Dataset
The dataset for this study is maintained in the [_zenodo_](https://zenodo.org/) database under the registration number: [_11577432_](https://zenodo.org/records/11577432). 
The dataset includes processed single-cell sequencing data from 25 types of solid tumors and 3 normal tissues.

## Contact
If you have any questions or feedback, feel free to contact me at zouzhuan@sr.gxmu.edu.cn.
# scRNA
Single-cell and spatial transcriptomic analyses reveals the metastasis mechanisms of osteosarcoma
