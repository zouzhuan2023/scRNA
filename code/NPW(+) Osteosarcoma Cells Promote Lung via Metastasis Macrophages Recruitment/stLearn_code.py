conda activate stlearn 
python
import stlearn as st
import scanpy as sc
import numpy as np
import os
import random
import pandas as pd
import matplotlib.pyplot as plt
random.seed(1000000)
np.random.seed(1000000)
data_dir = "path/"
data = st.Read10X(data_dir) 
data.var_names_make_unique()  
print(data.uns["spatial"].keys())
st.add.image(adata=data,
             imgpath=data_dir+"spatial/tissue_hires_image.png",
             library_id="PT_ST_Visium", visium=True)
minFeature =400
maxMT = 20
minCount = 200
minCell = 20
st.pp.filter_genes(data, min_cells=minCell ) 
data.obs['nFeature_Spatial'] = np.sum(data.X > 0, axis=1) 
data.obs['nCount_Spatial'] = np.sum(data.X, axis=1)  
data.obs['percent_mt'] = np.sum(data[:, data.var_names.str.startswith('MT-')].X, axis=1).A1 / np.sum(data.X, axis=1).A1 * 100
#data.obs['percent_rb'] = np.sum(data[:, data.var_names.str.startswith('RP[LS]')].X, axis=1).A1 / np.sum(data.X, axis=1).A1 * 100
data = data[(data.obs['nFeature_Spatial'] >= minFeature) & 
                         (data.obs['percent_mt'] <= maxMT) &
                         (data.obs['nCount_Spatial'] >= minCount), :]
st.pp.normalize_total(data) 
lrs = st.tl.cci.load_lrs(['connectomeDB2020_lit'], species='human')
print(len(lrs))
np.random.seed(1000000) 
st.tl.cci.run(data,
                        lrs, min_spots = 20, 
                        distance=None,  
                        n_pairs=200,   
                         n_cpus=80)

lr_info = data.uns['lr_summary'] 
print('\n', lr_info) 
st.tl.cci.adj_pvals(data, correct_axis='spot', pval_adj_cutoff=0.05, adj_method='fdr_bh')
st.pl.lr_summary(data, n_top=30, figsize=(9,6))
plt.savefig(data_dir+'s2.lr_summary30.pdf',dpi=600)
plt.close()
st.pl.lr_n_spots(data, n_top=30, figsize=(20,12),max_text=100)
plt.savefig(data_dir+'s2.lr_n_spots30.pdf',dpi=600)
plt.close()
stats = ['lr_scores', 'p_vals', 'p_adjs', '-log10(p_adjs)']
selected_lrs = ["COL1A1_ITGA2" ]
stats = ['lr_scores', 'p_adjs', '-log10(p_adjs)','lr_sig_scores']          
for lr in selected_lrs:
    fig, axes = plt.subplots(ncols=len(stats), figsize=(4 * len(stats), 6))  
    for j, stat in enumerate(stats):
           st.pl.lr_result_plot(data, use_result=stat, use_lr=lr,size=4.5, show_color_bar=False, ax=axes[j])
           axes[j].set_title(f'{lr} {stat}')
           handles, labels = axes[j].get_legend_handles_labels()        
           if handles:            
                  axes[j].legend(handles=handles, labels=labels, loc='best')
    plt.tight_layout()  
    plt.savefig(f'{data_dir}{lr}_result_plots.pdf', dpi=600)  
    plt.close()  

