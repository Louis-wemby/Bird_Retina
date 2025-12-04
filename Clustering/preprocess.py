### Preprocess data & Visualization

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

def calculate_qc_metrics(adata):
    adata.obs['n_counts'] = adata.X.sum(axis=1)
    adata.obs['n_genes'] = (adata.X > 0).sum(axis=1)
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.obs['percent_mito'] = adata[:, adata.var['mt']].X.sum(axis=1) / adata.X.sum(axis=1) * 100
    
def QC_violin(adata, stage, output_dir):
    qc_metrics = ['n_genes', 'n_counts', 'percent_mito']
    num_metrics = len(qc_metrics)
    num_cols = 3
    num_rows = (num_metrics + num_cols - 1) // num_cols
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(15, 5 * num_rows))
    axes = axes.flatten()
    for i, metric in enumerate(qc_metrics):
        sc.pl.violin(adata, keys=metric, groupby='batch', jitter=0, rotation=90, stripplot=False, palette='Set2', ax=axes[i], show=False)
        axes[i].set_title(f"{stage}: {metric}")
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{stage}_violinPlots.pdf", dpi=50)
    plt.show()
    
def QC_dot(adata, stage, output_dir):
    nCount_RNA = adata.obs['n_counts']
    nFeature_RNA = adata.obs['n_genes']
    percent_mt = adata.obs['percent_mito']
    mt_genes = adata.var_names[adata.var_names.str.startswith('mt-')]
    if len(mt_genes) > 0:
        correlations = [
            (nCount_RNA, percent_mt, 'percent.mt'),
            (nCount_RNA, nFeature_RNA, 'nFeature_RNA')
        ]
    else:
        correlations = [
            (nCount_RNA, percent_mt, 'percent.mt'),
            (nCount_RNA, nFeature_RNA, 'nFeature_RNA')
        ]
    plt.figure(figsize=(10, 4)) 
    for i, (x, y, ylabel) in enumerate(correlations, 1):
        corr, _ = pearsonr(x, y)
        plt.subplot(1, len(correlations), i)
        sns.scatterplot(x=x, y=y)
        plt.title(f'nCount_RNA vs {ylabel}\nPearson correlation: {corr:.2f}')
        plt.xlabel('nCount_RNA')
        plt.ylabel(ylabel)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{stage}_dotPlots.pdf", dpi=50)
    plt.show()
    
def cell_number(input_path, output_file):
    results = []
    for file_name in sorted(os.listdir(input_path)):
        if file_name.endswith('.h5ad'):
            file_path = os.path.join(input_path, file_name)
            adata = sc.read_h5ad(file_path)
            cell_count = adata.n_obs
            total_gene = adata.n_vars
            median_counts = round(adata.obs['n_counts'].median(), 2)
            median_genes = round(adata.obs['n_genes'].median(), 2)
            results.append({'file_name': file_name, 'cell_count': cell_count, 
                            'total_gene': total_gene, 'median_counts': median_counts, 'median_genes':median_genes})
    df = pd.DataFrame(results)
    df = df.sort_values(by='file_name')
    df.to_csv(output_file, index=False)
    print(f"Cell counts saved to {output_file}")
    
# Define the fucntion to filter "DOUBLETS_PREDICTION"
def doublets_filter(input_h5ad, output_h5ad, output_txt):
    adata = sc.read_h5ad(input_h5ad)
    if not issparse(adata.X):
        adata.X = scipy.sparse.csr_matrix(adata.X)
    if np.isnan(adata.X.data).any() or np.isinf(adata.X.data).any():
        print(f"Warning: adata.X contains NaN or Inf values in {input_h5ad}. Skipping this file.")
    # Doublets Prediction
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.05)
    doublet_scores, predicted_doublets = scrub.scrub_doublets()
    # Checking Prediction Results
    if predicted_doublets is None:
        print(f"Warning: Scrublet failed to predict doublets for {input_h5ad}. Skipping this file.")
    adata.obs['doublet_scores'] = doublet_scores
    adata.obs['predicted_doublets'] = predicted_doublets
    # Remove doublets prediction
    adata_filtered = adata[~predicted_doublets, :]
    # Save as new data after doublets prediction
    adata_filtered.write(output_h5ad)
    print(f"Processed {file_name} and saved as {output_h5ad}")
            
    # Save doublets_prediction information into ".txt" file
    doublets_info = adata.obs[predicted_doublets]
    doublets_info.to_csv(output_txt, sep='\t', index=True)
    print(f"Doublets information saved to {output_txt}")

# h5ad.txt contains h5ad file path, sample name and species name for processing data
h5ad_list = "h5ad.txt"
h5ad_outdir = "./1_BeforeQC_h5ad"
os.makedirs(h5ad_outdir, exist_ok=True)
h5ad_df = pd.read_csv(h5ad_list, sep="\t", names=["h5adpath", "sample", "species"])

adatas = {}
for h5ad, sam, spe in zip(h5ad_df.iloc[:, 0], h5ad_df.iloc[:, 1], h5ad_df.iloc[:, 2]):
    hread = sc.read_h5ad(h5ad)
    hread.obs_names_make_unique()
    hread.var_names_make_unique()
    hread.obs['batch'] = sam
    hread.obs['species'] = spe
    calculate_qc_metrics(hread)
    output_path = os.path.join(h5ad_outdir, f'{sam}.BeforeQC.h5ad')
    hread.write(output_path)
    adatas[sam] = hread
    
# Visualization and calculate cell number before qc
os.makedirs("figures", exist_ok=True)
adata = sc.concat(adatas, join='outer', index_unique='-')
QC_violin(adata, "Retina_BeforeQC", "figures/")
QC_dot(adata, "Retina_BeforeQC", "figures/")
cell_number(input_path=h5ad_outdir, output_file='Retina_CellNumber_BeforeQC.csv')

##### Step2: quality control for each h5ad, and get basic statistic after qc
h5ad_indir = "./1_BeforeQC_h5ad/"
h5ad_outdir = "./2_AfterQC_h5ad/"
os.makedirs(h5ad_outdir, exist_ok=True)

import logging
logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s', level=logging.INFO, datefmt='%Y-%m-%d %H:%M:%S')
def is_outlier(adata, metric: str, nmads: int, max_cutoff=None):
    M = adata.obs[metric]
    if max_cutoff != None:
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (max_cutoff < M)
    else:
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (np.median(M) + nmads * median_abs_deviation(M) < M)
    return outlier


adatas = {}
for file_name in os.listdir(h5ad_indir):
    if file_name.endswith(".h5ad"):
        sample_name = file_name.replace(".BeforeQC.h5ad", "")
        hdata = sc.read_h5ad(os.path.join(h5ad_indir, file_name))
        # QC-1: filter outliters
        sc.pp.calculate_qc_metrics(hdata, inplace=True, percent_top=[20], log1p=True)
        hdata.obs["outlier"] = (
        is_outlier(hdata, "log1p_total_counts", 5) 
                | is_outlier(hdata, "log1p_n_genes_by_counts", 5)
                | is_outlier(hdata, "pct_counts_in_top_20_genes", 5))
        # Logging.info(hdata.obs.outlier.value_counts())
        n0 = hdata.n_obs
        hdata.obs["mt_outlier"] = is_outlier(hdata, 'percent_mito', 3, 5)
        # Logging.info(f"    Total number of cells: {hdata.n_obs}")
        hdata = hdata[(~hdata.obs.outlier) & (~hdata.obs.mt_outlier)].copy()
        # QC-2: mt (,5%), min_cells [3,). min_genes  [200,), n_genes_by_counts (,5000]
        sc.pp.filter_cells(hdata, min_genes=150)
        n1 = hdata.n_obs
        print(f"min_genes=150 {n0 - n1} cells filtered, {n1} cells remaining")
        sc.pp.filter_genes(hdata, min_cells=3)
        hdata = hdata[hdata.obs['percent_mito'] < 5, :]
        n2 = hdata.n_obs
        # Print(f"percent_mito < 5 {n2 - n1} cells filtered, {n2} cells remaining")
        # Save as new h5ad (AfterQC.h5ad)
        output_file_name = f"{sample_name}.AfterQC.h5ad"
        hdata.write(os.path.join(h5ad_outdir, output_file_name))
        adatas[sample_name] = hdata

# Visualization and calculate cell number after qc
os.makedirs("figures", exist_ok=True)
adata = sc.concat(adatas, join='outer', index_unique='-')
QC_violin(adata, "Retina_AfterQC", "figures/")
QC_dot(adata, "Retina_AfterQC", "figures/")
cell_number(input_path=h5ad_outdir, output_file='Retina_CellNumber_AfterQC.csv')

##### Step3: process doublets for each h5ad, and get basic statistic after doublets filtering
h5ad_indir = "./2_AfterQC_h5ad/"
h5ad_outdir = "./3_AfterDP_h5ad/"
os.makedirs(h5ad_outdir, exist_ok=True)   

# Doublets_filtering
for file_name in os.listdir(h5ad_indir):
    if file_name.endswith(".h5ad"):
        sample_name = file_name.replace(".AfterQC.h5ad", "")
        h5ad_infile = os.path.join(h5ad_indir, file_name)
        h5ad_outfile = os.path.join(h5ad_outdir, f"{sample_name}.AfterDP.h5ad")
        doublets_txt = os.path.join(h5ad_outdir, f"{sample_name}.Doublets.txt")
        doublets_filter(h5ad_infile, h5ad_outfile, doublets_txt)
        

# Visualization and calculate cell number after qc
os.makedirs("figures", exist_ok=True)
adatas = {}
for file_name in os.listdir(h5ad_outdir):
    if file_name.endswith(".h5ad"):
        sample_name = file_name.replace(".AfterDP.h5ad", "")
        hdata = sc.read_h5ad(os.path.join(h5ad_outdir, file_name))
        adatas[sample_name] = hdata
adata = sc.concat(adatas, join='outer', index_unique='-')

QC_violin(adata, "Retina_AfterDP", "figures/")
QC_dot(adata, "Retina_AfterDP", "figures/")
cell_number(input_path=h5ad_outdir, output_file='Retina_CellNumber_AfterDP.csv')

##### Step4: process doublets for each h5ad, and get basic statistic after doublets filtering
h5ad_indir = "./3_AfterDP_h5ad/"
h5ad_outdir = "./4_AfterNL_h5ad/"
os.makedirs(h5ad_outdir, exist_ok=True)   

# Doublets_filtering
for file_name in os.listdir(h5ad_indir):
    if file_name.endswith(".h5ad"):
        sample_name = file_name.replace(".AfterDP.h5ad", "")
        h5ad_infile = os.path.join(h5ad_indir, file_name)
        h5ad_outfile = os.path.join(h5ad_outdir, f"{sample_name}.AfterNL.h5ad")
        adata = sc.read_h5ad(h5ad_infile)
        adata.layers["counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.layers["log1p"] = adata.X.copy()
        adata.X = adata.layers["counts"].copy()
        adata.write_h5ad(h5ad_outfile, compression='gzip')
  
