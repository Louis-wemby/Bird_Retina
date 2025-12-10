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
# from harmony import harmonize

def umapRun(adata, batch, prefix):
    batch_key = batch
    # Optional -- whether to find new hvgs in subclass or use hvgs in all cells
    # sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key=batch_key)
    adata_hvg = adata[:, adata.var.highly_variable].copy()

    # Dimensionality reduction
    sc.tl.pca(adata_hvg, n_comps=50)
    sc.pl.pca_variance_ratio(adata_hvg, n_pcs=50)
    variance_ratio = adata_hvg.uns['pca']['variance_ratio']
    cum_var = np.cumsum(variance_ratio)
    print(cum_var)
    
    # Integration with harmony
    # adata_hvg.obsm["Harmony"] = harmonize(adata_hvg.obsm["X_pca"], adata_hvg.obs, batch_key= batch_key)

    # Visualization
    sc.pp.neighbors(adata_hvg, use_rep='X_pca')
    sc.tl.umap(adata_hvg)

    for resolution in [0.3,0.5,0.6,0.7,0.8,1.0]:
        key = f'leiden_subclass_{resolution}'
        sc.tl.leiden(adata_hvg, resolution=resolution, key_added=key)
    cell_clusters = [f"leiden_subclass_{res}" for res in [0.3,0.5,0.6,0.7,0.8,1.0]]
    sc.pl.embedding(adata_hvg, basis=f"X_umap", color=cell_clusters, ncols=3, legend_loc="on data", save="." + prefix + ".clusters.pdf")
    
    for annot in ['obs', 'uns', 'obsm', 'obsp']:
        for key in eval(f"adata_hvg.{annot}"):
            if key not in eval(f"adata.{annot}"):
                try:
                    exec(f"adata.{annot}['{key}'] = adata_hvg.{annot}['{key}'].copy()")
                except:
                    exec(f"adata.{annot}['{key}'] = adata_hvg.{annot}['{key}']")
    adata.obsm['X_subclass_umap'] = adata_hvg.obsm['X_umap'].copy()
    return adata

## An example
adata = sc.read_h5ad("5_Cluster_h5ad/bs_md_.Seurat-3000_npcs30.h5ad")
adata.X = adata.layers["log1p"].copy()
bdata = adata[adata.obs['celltype'] == "Photoreceptors"].copy()
batch_key = "batch"
prefix = "PR"
umapRun(bdata, batch_key, prefix)
