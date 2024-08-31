conda activate SEDR_Env
cd /path
## Load the package ####
python
import scanpy as sc
import pandas as pd
from sklearn import metrics
import torch
import matplotlib.pyplot as plt
import seaborn as sns
import os
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')
import SEDR
import anndata as ad
import numpy as np
from PIL import Image
Image.MAX_IMAGE_PIXELS = None
from pathlib import Path
from matplotlib.colors import LinearSegmentedColormap
###gpu##
device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
###path###
data_root = Path('/path/MLN09')
n_clusters = 10
adata = sc.read_visium(data_root)
adata.var_names_make_unique()

adata.layers['count'] = adata.X.toarray()
### Set quality control parameters###
minFeature = 400 
minCount = 800    
maxMT = 20       
minCell = 20      
### Gene filtering###
sc.pp.filter_genes(adata, min_cells=minCell)
### Calculate the total counts for each cell###
adata.obs['total_counts'] = adata.X.sum(axis=1)
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)

### Calculate the proportion of mitochondrial genes###
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / adata.obs['total_counts'] * 100

### Apply quality control criteria for filtering###
adata = adata[(adata.obs['n_genes'] >= minFeature) & 
              (adata.obs['total_counts'] >= minCount) &
              (adata.obs['percent_mito'] <= maxMT), :]
### Data normalization###
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.scale(adata)
### Construct a spatial plot###
graph_dict = SEDR.graph_construction(adata, 12)
##### Export the spatial matrix###
adj_norm = graph_dict['adj_norm']
assert adj_norm is not None, "adj_norm is not defined"
adj_norm_indices = adj_norm.coalesce().indices().numpy()
adj_norm_values = adj_norm.coalesce().values().numpy()
# ´´½¨ adj_norm µÄ DataFrame
df_adj_norm = pd.DataFrame({
    'row': adj_norm_indices[0],
    'col': adj_norm_indices[1],
    'value': adj_norm_values
})
###Export adj_norm to CSV###
df_adj_norm.to_csv('adj_norm.csv', index=False)
##########
### Set a random seed###
random_seed = 2023
SEDR.fix_seed(random_seed)

### Select the execution device ###
device = 'cpu'

### Initialize the SEDR model using the dataset's expression matrix (adata.X) and the graph dictionary (graph_dict)###
### Set the mode to 'imputation'
sedr_net = SEDR.Sedr(adata.X, graph_dict, mode='imputation', device=device)
###Model training###
using_dec = True
if using_dec:
    sedr_net.train_with_dec()  
else:
    sedr_net.train_without_dec()  

### Process data and extract features###
sedr_feat, _, _, _ = sedr_net.process()
adata.obsm['sedr'] = sedr_feat

### Reconstruct data###
de_feat = sedr_net.recon()
adata.obsm['de_feat'] = de_feat

newcmp = LinearSegmentedColormap.from_list('new', ['#EEEEEE','#E7481B'], N=1000)
matrix_name = 'de_feat'

###Check if the matrix name exists in adata.obsm keys###
if matrix_name in adata.obsm.keys():
    matrix = adata.obsm[matrix_name]
    
    gene_names = adata.var.index.tolist()
    
    df = pd.DataFrame(matrix, columns=gene_names, index=adata.obs.index)
    
    df.to_csv(f'{matrix_name}_export.csv')

list_genes = ['CDKN1A','IBSP','CDKN2A']

###Retrieve the corresponding denoised gene expression data###
for gene in list_genes:
    idx = adata.var.index.tolist().index(gene)
    adata.obs[f'{gene}(denoised)'] = adata.obsm['de_feat'][:, idx]


plt.figure(figsize=(8, 8))
fig, axes = plt.subplots(1,len(list_genes),figsize=(1*(len(list_genes)),1*1))
axes = axes.ravel()

if not np.issubdtype(adata.obsm['spatial'].dtype, np.number):
    adata.obsm['spatial'] = adata.obsm['spatial'].astype(float)

for i, gene in enumerate(list_genes):
    if f'{gene}(denoised)' in adata.obs.columns:
        sc.pl.spatial(
            adata, 
            color=f'{gene}(denoised)', 
            ax=axes[i], 
            vmax='p99', 
            vmin='p1', 
            alpha_img=0, 
            cmap=newcmp, 
            colorbar_loc=None,
            size=1.4,  
            show=False
        )
    else:
        axes[i].set_visible(False)

for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('')
    ax.set_ylabel('')

plt.subplots_adjust(wspace=0, hspace=0)
output_file = "/path/gene.pdf"
plt.savefig(output_file, format='pdf')

## Drawing with a scale£º
fig, axes = plt.subplots(1, len(list_genes), figsize=(4 * len(list_genes), 4 * 1))
axes = axes.ravel()

if not np.issubdtype(adata.obsm['spatial'].dtype, np.number):
    adata.obsm['spatial'] = adata.obsm['spatial'].astype(float)

for i, gene in enumerate(list_genes):
    if f'{gene}(denoised)' in adata.obs.columns:
        sc.pl.spatial(
            adata,
            color=f'{gene}(denoised)',
            ax=axes[i],
            vmax='p99',
            vmin='p1',
            alpha_img=0,
            cmap=newcmp,
            colorbar_loc="right",  
            size=1.4,
            show=False
        )
    else:
        axes[i].set_visible(False)

###Adjust the plot style###
for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('')
    ax.set_ylabel('')

plt.subplots_adjust(wspace=0.1, hspace=0)
output_file = "/path/gene.pdf"
plt.savefig(output_file, format='pdf')
plt.show()
### Save the matrix information of CDKN2A###
import pandas as pd
import scanpy as sc

gene_name = 'CDKN2A'
if gene_name not in adata.var.index:
    raise ValueError(f"Gene {gene_name} not found in adata.var.index")

idx = adata.var.index.tolist().index(gene_name)
CDKN2A_matrix = adata.obsm['de_feat'][:, idx]

CDKN2A_df = pd.DataFrame(CDKN2A_matrix, index=adata.obs.index, columns=[f'{gene_name}(denoised)'])

CDKN2A_df.to_csv(f'{gene_name}_matrix_export.csv')



