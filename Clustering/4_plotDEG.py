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

def degRun(adata, prefix, r):
    ## DEGs
    os.makedirs("DEGs", exist_ok=True)
    sc.tl.rank_genes_groups(adata,
                            groupby=f"leiden_{r}",
                            method='wilcoxon', key_added=f"DEGs_res{r}")
    sc.pl.rank_genes_groups_dotplot(adata,
                                        groupby=f"leiden_{r}", dendrogram=False,
                                        standard_scale="var", n_genes=5, key=f"DEGs_res{r}",
                                        save='DEGs_'+prefix+'_res'+str(r)+'.png')
    degr = f"DEGs_res{r}"
    result = adata.uns[degr]
    groups = result['names'].dtype.names
    degs_df = pd.DataFrame({group + '_' + key[:1]: result[key][group] 
                                for group in groups for key in ['names', 'pvals']}).head(100)
    degs_df.to_csv(f"DEGs/{prefix}_res{str(r)}_DEG-top100.csv", index=False, sep="\t")


    sc.tl.filter_rank_genes_groups(adata,
                               min_in_group_fraction=0.2, max_out_group_fraction=0.3,
                               key=f"DEGs_res{r}", key_added=f"DEGsFilt_res{r}",
                               use_raw=False)
    sc.pl.rank_genes_groups_dotplot(adata,
                                        groupby=f"leiden_{r}", dendrogram=False,
                                        standard_scale="var", n_genes=5, key=f"DEGsFilt_res{r}",
                                        save='DEGsFilt_'+prefix+'_res'+str(r)+'.png')
    degr = f"DEGsFilt_res{r}"
    result = adata.uns[degr]
    groups = result['names'].dtype.names
    degs_df = pd.DataFrame({group + '_' + key[:1]: result[key][group] 
                                for group in groups for key in ['names', 'pvals']}).head(100)
    degs_df.to_csv(f"DEGs/{prefix}_res{str(r)}_DEGfilt-top100.csv", index=False, sep="\t")

def QC_violin(adata, r, stage, output_dir):
    qc_metrics = ['n_genes', 'n_counts', 'percent_mito']
    num_metrics = len(qc_metrics)
    num_cols = 3
    num_rows = (num_metrics + num_cols - 1) // num_cols
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows))
    axes = axes.flatten()
    for i, metric in enumerate(qc_metrics):
        sc.pl.violin(adata, keys=metric, groupby=r, jitter=0, rotation=90, stripplot=False, palette='Set2', ax=axes[i], show=False)
        axes[i].set_title(f"{stage}: {metric}")
        group_vals = adata.obs.groupby(r)[metric]
        meds = group_vals.median()
        ax = axes[i]
        for j, (g, v) in enumerate(meds.items()):
            x_pos = j
            ax.text(x_pos, v, f'{v:.1f}', ha='center', va='bottom', fontsize=10)
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{stage}_violinPlots.pdf", dpi=50)
    plt.show()

h5ad_in = "5_Cluster_h5ad/bs_md_.Seurat-3000_npcs30.h5ad"
prefix = os.path.basename(h5ad_in).split(".h5ad")[0]

adata = sc.read_h5ad(h5ad_in)
adata.X = adata.layers['log1p'].copy()
degRun(adata, prefix, 0.8)
QC_violin(adata, "leiden_0.8", f"{prefix}_leiden_0.8", "figures/")
