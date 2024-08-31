#######
conda activate SEDR_Env
cd /path/
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
# gpu
device = 'cuda:0' if torch.cuda.is_available() else 'cpu'
# path
data_root = Path('/path')
n_clusters = 10
adata = sc.read_visium(data_root)
adata.var_names_make_unique()
adata.layers['count'] = adata.X.toarray()
sc.pp.filter_genes(adata, min_cells=50)
adata.obs['total_counts'] = adata.X.sum(axis=1) 
adata.obs['n_genes'] = (adata.X > 0).sum(axis=1) 
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / adata.obs['total_counts'] * 100 
min_count = 1500
min_feature = 500
max_mt = 25

adata = adata[(adata.obs['total_counts'] >= min_count) & 
              (adata.obs['n_genes'] >= min_feature) & 
              (adata.obs['percent_mito'] <= max_mt), :]
if 'ZMYND8' in adata.var.index:
    print("ZMYND8 is present in the dataset.")
else:
    print("ZMYND8 is not found in the dataset.")
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.scale(adata)
graph_dict = SEDR.graph_construction(adata, 12)
adj_norm = graph_dict['adj_norm']
assert adj_norm is not None, "adj_norm is not defined"
adj_norm_indices = adj_norm.coalesce().indices().numpy()
adj_norm_values = adj_norm.coalesce().values().numpy()
df_adj_norm = pd.DataFrame({
    'row': adj_norm_indices[0],
    'col': adj_norm_indices[1],
    'value': adj_norm_values
})
df_adj_norm.to_csv('adj_norm.csv', index=False)
random_seed = 2023
SEDR.fix_seed(random_seed)
device = 'cpu'
sedr_net = SEDR.Sedr(adata.X, graph_dict, mode='imputation', device=device)

using_dec = True
if using_dec:
    sedr_net.train_with_dec()
else:
    sedr_net.train_without_dec()

sedr_feat, _, _, _ = sedr_net.process()
adata.obsm['sedr'] = sedr_feat

# reconstruction
de_feat = sedr_net.recon()
adata.obsm['de_feat'] = de_feat

from matplotlib.colors import ListedColormap, LinearSegmentedColormap
newcmp = LinearSegmentedColormap.from_list('new', ['#EEEEEE','#E7481B'], N=1000)
####################
matrix_name = 'de_feat'

if matrix_name in adata.obsm.keys():
    matrix = adata.obsm[matrix_name]
    gene_names = adata.var.index.tolist()
    df = pd.DataFrame(matrix, columns=gene_names, index=adata.obs.index)
    df.to_csv(f'{matrix_name}_export.csv')
    print(f"Matrix '{matrix_name}' has been exported to '{matrix_name}_export.csv'.")
else:
    print(f"Matrix '{matrix_name}' not found in adata.obsm.")
####################
list_genes = ['ZMYND8','GPI','PGAM1','SUCLG1','MDH1','MDH2']
for gene in list_genes:
    idx = adata.var.index.tolist().index(gene)
    adata.obs[f'{gene}(denoised)'] = adata.obsm['de_feat'][:, idx]
plt.figure(figsize=(8, 8))
fig, axes = plt.subplots(1,len(list_genes),figsize=(4*(len(list_genes)),4*1))
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

for ax in axes:
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_xlabel('')
    ax.set_ylabel('')

plt.subplots_adjust(wspace=0, hspace=0)
output_file = "out.pdf"
plt.savefig(output_file, format='pdf')



