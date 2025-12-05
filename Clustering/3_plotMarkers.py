## Draw dotplots for marker genes

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

# Dictionary of marker genes in retinal cells
retina_marker_genes = {
    "Housekeeping_Gene": ["GAPDH"],
    "Photoreceptor_Rods": ["RHO", "PDE6B", "PDE6G", "MAFA"],
    "Photoreceptor_Cones": ["ARR3", "RS1", "OPN1SW", "OPN2SW", "OPN1LW", "LOC135286153"],
    "Horizontal_Cells": ["ONECUT1", "ONECUT2", "ONECUT3", "LHX1", "ISL1", "IPCEF1", "OXT", "CALB1", "SLC4A3", "WDR72"], 
    "Bipolar_Cells": ["VSX1", "VSX2", "ISL1", "TRPM1", "GRIK1", "FEZF2", "PRKCA", "OTX1", "SOX5", "CABP5"], 
    "Amacrine_Cells": ["PAX6", "PROX11", "SLC32A1", "SLC6A1", "GAD1", "GAD2", "SLC6A9", "TFAP2A", "TFAP2B", "TFAP2C"],
    "Retinal_Ganglion_Cells": ["RBPMS", "RBPMS2", "SLC17A6", "THY1", "POU4F1", "POU4F2", "POU4F3", "SATB1", "RELN", "CHRNB2", "SS2"],
    "Muller_Glia": ["SLC1A3", "RLBP1", "CHRDL1", "WIF1", "PSCA", "TMEM123", "FOXD1", "GLUL"],
    "Oligodendrocytes": ["OLIG2", "PTPRZ1", "PDGFRA", "SOX2", "HES1", "MYT1", "BMI1", "SOX10", "CNP", "PLP1", "MBP", "BACS1"],
    "RetinalPigmentEpithelium": ["RPE65", "PMEL", "BEST1", "LRAT", "TYR"]
}

h5ad = "5_Cluster_h5ad/bs_md_.Seurat-3000_npcs30.h5ad"  # Select h5ad file with parameters bs & md after clustering
prefix = os.path.basename(h5ad).split(".h5ad")[0]
adata = sc.read_h5ad(h5ad)

marker_in_adata = {cell_type: [gene for gene in genes if gene in adata.var_names]
                       for cell_type, genes in retina_marker_genes.items()}
sc.pl.dotplot(adata,
              marker_in_adata,
              groupby='leiden_0.8',
              standard_scale="var",
              swap_axes=False,
              save='markers_'+prefix+'_res0.8.png')
