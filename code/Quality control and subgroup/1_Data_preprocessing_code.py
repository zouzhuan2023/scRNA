
conda activate scanpy_env
vim /path/script.py
import os
import scanpy as sc
from tqdm import tqdm
import numpy as np

def preprocess_data(folder_path):
    # Load the data
    adata = sc.read_10x_mtx(folder_path, var_names='gene_symbols', cache=False)
    
    # Predict doublets using scrublet
    sc.external.pp.scrublet(adata, random_state=112,n_prin_comps=14)
    
    # Filter out predicted doublets
    if 'predicted_doublet' in adata.obs.columns:
        adata = adata[adata.obs['predicted_doublet'] == False].copy()
    else:
        print("Warning: 'predicted_doublet' column not found in adata.obs. Doublet filtering skipped.")


    # Continue with further preprocessing
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=10)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    upper_lim = np.quantile(adata.obs.n_genes_by_counts.values, .98)
    adata = adata[adata.obs.n_genes_by_counts < upper_lim].copy()
    adata = adata[adata.obs.pct_counts_mt < 10].copy()
    adata = adata[adata.obs['total_counts'] <= 20000, :]
    
    return adata

data_path = '/path/'
output_path = '/path'

# Ensure output directory exists
os.makedirs(output_path, exist_ok=True)

cancer_types = os.listdir(data_path)
adata_objects = []

for cancer_type in tqdm(cancer_types, desc="Cancer Types"):
    folder_path = os.path.join(data_path, cancer_type)
    
    # Check if the path is a directory to avoid processing files.
    if not os.path.isdir(folder_path):
        continue
    
    try:
        adata = preprocess_data(folder_path)
        prefix = f"{cancer_type}_"
        adata.obs_names = [prefix + barcode for barcode in adata.obs_names]
        adata_objects.append(adata)
    except Exception as e:
        print(f"Error processing {folder_path}: {e}")

if adata_objects:
    all_data = sc.concat(adata_objects, join='outer')
    all_data.obs['predicted_doublet'] = all_data.obs['predicted_doublet'].astype(str)
    # Save the combined data
    save_file_path = os.path.join(output_path, 'st_qc.h5ad')
    all_data.write_h5ad(save_file_path)
    print(f"All data written to {save_file_path}")
else:
    print("No data processed.")

python /path/script.py
python

import scanpy as sc
import scvi
import os
import math
import itertools
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import bbknn
adata = sc.read_h5ad('/path/st_qc.h5ad')
dir = '/path/'


adata.obs['cell_identifier'] = adata.obs.index
def extract_batch(cell_identifier):
    parts = cell_identifier.split('_')
    return parts[0] if len(parts) > 0 else 'Unknown'
adata.obs['batch'] = adata.obs['cell_identifier'].apply(extract_batch)
batch_counts = adata.obs['batch'].value_counts()

unique_batches = adata.obs['batch'].unique()
number_of_batches = len(unique_batches)


#batches_to_keep = adata.obs['batch'].value_counts()[lambda x: x >= 200].index
#adata = adata[adata.obs['batch'].isin(batches_to_keep)].copy()

adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

import scipy.sparse
df = pd.DataFrame(adata.X[:5, :5].todense() if scipy.sparse.issparse(adata.X) else adata.X[:5, :5], 
                  index=adata.obs_names[:5], 
                  columns=adata.var_names[:5])
print(df)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable]    
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])  
sc.pp.scale(adata, max_value=10)

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
plt.savefig(dir+"pca_variance.pdf")
adata.write(dir + 'st_PCA.h5ad')
adata = sc.read_h5ad('path/st_PCA.h5ad')
dir = 'path/'


sc.external.pp.bbknn(adata, batch_key= "batch")

sc.pp.neighbors(adata, n_neighbors=10, n_pcs=8)
sc.tl.umap(adata)
sc.tl.tsne(adata)
sc.tl.leiden(adata, resolution = 0.5)
leiden_counts = adata.obs['leiden'].value_counts()
print(leiden_counts) 
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
plt.close('all')
plt.rcParams['savefig.dpi'] = 300
plt.figure(figsize=(20, 10))
sc.pl.tsne(adata, color='leiden', legend_loc='right margin',show=False)
plt.tight_layout()
plt.savefig(dir+"tsne_leiden.png",bbox_inches='tight')



plt.close('all')
plt.rcParams['savefig.dpi'] = 300
plt.figure(figsize=(20, 10))
sc.pl.umap(adata, color='leiden', legend_loc='right margin',show=False)
plt.tight_layout()
plt.savefig(dir+"umap_leiden.png",bbox_inches='tight')
adata.write(dir + 'st_adata_leiden.h5ad')
###############################################3

import scanpy as sc
import matplotlib.pyplot as plt

plt.rcParams['savefig.dpi'] = 600
plt.figure(figsize=(30, 10))

sc.pl.umap(adata, color='batch', legend_loc='right margin', show=False)

plt.tight_layout()
output_dir = 'path/'
plt.savefig(output_dir + "umap_Batch.png", bbox_inches='tight')
plt.close()
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
adata.write(dir + 'st_adata_leiden.h5ad')
dir = 'path'
adata = sc.read_h5ad(dir + 'st_adata_leiden.h5ad')
result = adata.uns['rank_genes_groups']

groups = result['names'].dtype.names
all_groups_df = pd.DataFrame()
for group in groups:
    group_data = pd.DataFrame({
        f"{group}_names": result['names'][group],
        f"{group}_pvals": result['pvals'][group],
        f"{group}_logfoldchanges": result['logfoldchanges'][group],
        f"{group}_scores": result['scores'][group]
    })
    all_groups_df = pd.concat([all_groups_df, group_data], axis=1)
all_groups_df.to_csv(dir + 'different_gene.csv', index=False)


