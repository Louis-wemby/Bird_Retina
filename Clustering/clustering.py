import os
import anndata
import warnings
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import scrublet as scr
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.sparse import issparse
from scipy.stats import median_abs_deviation

h5ad_in = "4_AfterNL_h5ad/bs_md_.AfterNL.h5ad"  # Select h5ad file with parameter bs & md after preprocessing
h5ad_out = "5_Cluster_h5ad/bs_md_.Seurat-3000_npcs30.h5ad"  
prefix = os.path.basename(h5ad_out).split(".h5ad")[0]

adata = sc.read_h5ad(h5ad_in)
## Highly variable genes
batch_key = "batch"
adata.X = adata.layers["log1p"].copy()
sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key= batch_key)  # Parameters n_top_genes and the following n_comps/n_pcs usually need tuning for better clustering results.

## Integrate maps
adata_hvg = adata[:, adata.var.highly_variable].copy()
## Dimensionality reduction
sc.tl.pca(adata_hvg, n_comps=50)
sc.pl.pca_variance_ratio(adata_hvg, n_pcs=30)
variance_ratio = adata_hvg.uns['pca']['variance_ratio']
cum_var = np.cumsum(variance_ratio)
print(cum_var)

# Visualization
sc.pp.neighbors(adata_hvg, n_pcs=30)
sc.tl.umap(adata_hvg)

for resolution in [0.3,0.5,0.6,0.7,0.8,1.0]:
    key = f'leiden_{resolution}'
    sc.tl.leiden(adata_hvg, resolution=resolution, key_added=key)
cell_clusters = [f"leiden_{res}" for res in [0.3,0.5,0.6,0.7,0.8,1.0]]
sc.pl.embedding(adata_hvg, basis=f"X_umap", color=cell_clusters, ncols=3, legend_loc="on data")

for annot in ['obs', 'uns', 'obsm', 'obsp']:
    for key in eval(f"adata_hvg.{annot}"):
        if key not in eval(f"adata.{annot}"):
            try:
                exec(f"adata.{annot}['{key}'] = adata_hvg.{annot}['{key}'].copy()")
            except:
                exec(f"adata.{annot}['{key}'] = adata_hvg.{annot}['{key}']")
adata.X = adata.layers['counts'].copy()
adata.write_h5ad(h5ad_out, compression='gzip')

from matplotlib import rc_context
with rc_context({'figure.figsize': (4, 6)}):
    sc.pl.embedding(adata, basis=f"spatial", color='leiden_0.8', save=('.' + prefix + '.res0.8.pdf'), size=0.6)
