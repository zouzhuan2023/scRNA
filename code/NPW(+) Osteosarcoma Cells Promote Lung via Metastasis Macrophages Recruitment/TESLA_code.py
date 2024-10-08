conda activate scanpy_env
python
import os,csv,re, time
import pickle
import random
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
from scipy import stats
from scipy.sparse import issparse
import scanpy as sc
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import cv2
import TESLA as tesla
from IPython.display import Image
from scanpy import read_10x_h5
dir = '/path/'
adata = read_10x_h5(dir+"filtered_feature_bc_matrix.h5")

spatial=pd.read_csv(dir+"tissue_positions_list.csv",sep=",",header=None,na_filter=False,index_col=0) 
adata.obs["x1"] = spatial[1]
adata.obs["x2"] = spatial[2]
adata.obs["x3"] = spatial[3]
adata.obs["x4"] = spatial[4]
adata.obs["x5"] = spatial[5]
adata.obs["array_x"] = adata.obs["x2"]
adata.obs["array_y"] = adata.obs["x3"]
adata.obs["pixel_x"] = adata.obs["x4"]
adata.obs["pixel_y"] = adata.obs["x5"]

scale=float( pd.read_json(dir+'scalefactors_json.json',orient='index').loc['tissue_hires_scalef'] )

adata.obs['pixel_x'] = pd.to_numeric(adata.obs['pixel_x'], errors='coerce')
pixel_x = np.array(adata.obs['pixel_x']).astype(int) 
adata.obs['pixel_y'] = pd.to_numeric(adata.obs['pixel_y'], errors='coerce')

pixel_y=np.array(adata.obs['pixel_y']).astype(int)

adata.obs['pixel_x'] = pixel_x
adata.obs['pixel_y'] = pixel_y

adata.obs['x1'] = pd.to_numeric(adata.obs['x1'], errors='coerce')
adata = adata[adata.obs['x1'] == 1]


adata.var_names=[i.upper() for i in list(adata.var_names)]
adata.var["genename"]=adata.var.index.astype("str")

counts=adata
img=cv2.imread(dir+"cytassist_image.tiff")

resize_factor=1000/np.min(img.shape[0:2])
resize_width=int(img.shape[1]*resize_factor)
resize_height=int(img.shape[0]*resize_factor)
counts.var.index=[i.upper() for i in counts.var.index]
counts.var_names_make_unique()
counts.raw=counts
sc.pp.log1p(counts) 
if issparse(counts.X):counts.X=counts.X.A.copy()


cnt=tesla.cv2_detect_contour(img, apertureSize=5,L2gradient = True)

binary=np.zeros((img.shape[0:2]), dtype=np.uint8)
cv2.drawContours(binary, [cnt], -1, (1), thickness=-1)
#Enlarged filter
cnt_enlarged = tesla.scale_contour(cnt, 1.05)
binary_enlarged = np.zeros(img.shape[0:2])
cv2.drawContours(binary_enlarged, [cnt_enlarged], -1, (1), thickness=-1)
img_new = img.copy()
cv2.drawContours(img_new, [cnt], -1, (255), thickness=50)
img_new=cv2.resize(img_new, ((resize_width, resize_height)))



res=30
enhanced_exp_adata=tesla.imputation(img=img, raw=counts, cnt=cnt, genes=counts.var.index.tolist(), shape="None", res=res, s=1, k=2, num_nbs=10)

enhanced_exp_adata.X

cnt_color = clr.LinearSegmentedColormap.from_list('magma', ["#3b0f6f","#8c2980","#f66e5b", "#fd9f6c", "#fbfcbf"], N=100)
g="NPW"
enhanced_exp_adata.obs[g]=enhanced_exp_adata.X[:,enhanced_exp_adata.var.index==g]
fig=sc.pl.scatter(enhanced_exp_adata,alpha=1,x="y",y="x",color=g,color_map=cnt_color,show=False,size=22) 
fig.set_aspect('equal', 'box')
fig.invert_yaxis()
plt.gcf().set_dpi(600)
fig.figure.show()
plt.savefig(dir+"TESELA_NPW_fig.pdf")


genes=['CD163','C1QA']
genes=list(set([i for i in genes if i in enhanced_exp_adata.var.index ]))
#target_size can be set to "small" or "large".
pred_refined, target_clusters, c_m=tesla.annotation(img=img, 
                                                    binary=binary,
                                                    sudo_adata=enhanced_exp_adata, 
                                                    genes=genes, 
                                                    resize_factor=resize_factor,
                                                    num_required=2, 
                                                    target_size="small")

ret_img=tesla.visualize_annotation(img=img, 
                                   binary=binary, 
                                   resize_factor=resize_factor,
                                   pred_refined=pred_refined, 
                                   target_clusters=target_clusters, 
                                   c_m=c_m)
cv2.imwrite(dir+'TESELA_CD163_C1QA.jpg', ret_img)

genes=['CD163','C1QA','NPW']
genes=list(set([i for i in genes if i in enhanced_exp_adata.var.index ]))
#target_size can be set to "small" or "large".
pred_refined, target_clusters, c_m=tesla.annotation(img=img, 
                                                    binary=binary,
                                                    sudo_adata=enhanced_exp_adata, 
                                                    genes=genes,
                                                    resize_factor=resize_factor,
                                                    num_required=3,
                                                    target_size="small")
#Plot
ret_img=tesla.visualize_annotation(img=img, 
                                   binary=binary, 
                                   resize_factor=resize_factor,
                                   pred_refined=pred_refined, 
                                   target_clusters=target_clusters, 
                                   c_m=c_m)
cv2.imwrite(dir+'TESELA_NPW_Mar_tumor.jpg', ret_img)


