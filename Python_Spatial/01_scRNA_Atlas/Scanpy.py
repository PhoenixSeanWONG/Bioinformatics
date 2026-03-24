import scanpy as sc
import pandas as pd
import numpy as np
from typing import Optional, Union
from pathlib import Path
from anndata import AnnData

# --- 1. Global Settings ---
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=100, facecolor='white', frameon=False)

def run_scrna_pipeline(
    data_input: Optional[Union[str, Path, AnnData]] = None,
    min_genes: int = 220,
    max_mt_percent: float = 20.0,
    resolution: float = 0.5
) -> AnnData:
    """
    Advanced Pipeline for Melanoma scRNA-seq Atlas Construction.
    
    Args:
        data_input: Path to .h5ad file or an existing AnnData object.
        min_genes: Minimum number of genes expressed required for a cell to pass filtering.
        max_mt_percent: Maximum mitochondrial gene percentage allowed.
        resolution: Clustering resolution for the Leiden algorithm.
    """

    # --- Step 1: Smart Data Loading ---
    if data_input is None:
        print("No input provided. Loading demonstration dataset (PBMC 3k)...")
        adata = sc.datasets.pbmc3k()
    elif isinstance(data_input, AnnData):
        adata = data_input.copy()  # Avoid modifying the original object in memory
    else:
        print(f"Loading dataset from: {data_input}")
        adata = sc.read_h5ad(data_input)

    # --- Step 2: Quality Control & Filtering ---
    # In melanoma, preserving low-quality samples vs. removing dead cells is a delicate balance.
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True) # sc.pp: Preprocessing
    
    # Logging the filtering process (Very important for reproducibility)
    n_cells_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=3)
    adata = adata[adata.obs.pct_counts_mt < max_mt_percent, :].copy()
    
    print(f"QC Filter: {n_cells_before - adata.n_obs} cells removed ({adata.n_obs} cells remaining).")

    # --- Step 3: Normalization & HVG Selection ---
    # We use raw data storage to keep 'Phase I' denoising insights accessible.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    # Identify HVGs - Crucial for capturing tumor heterogeneity (e.g., AXL/MITF programs)
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    adata.raw = adata  # Keeping full transcriptomic profile in .raw
    adata = adata[:, adata.var.highly_variable].copy()

    # --- Step 4: Dimension Reduction & Integration ---
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    
    # Note: If integrating multiple melanoma patients, Harmony or Scanorama should be used here.
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=40)
    sc.tl.umap(adata) # sc.tl: Tools

    # --- Step 5: Clustering & Markers ---
    # Using Leiden (standard) for robust community detection.
    sc.tl.leiden(adata, resolution=resolution)
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

    print("Pipeline completed successfully.")
    return adata

if __name__ == "__main__":
    # Execute the pipeline
    final_adata = run_scrna_pipeline(resolution=0.6)
    
    # Visualization: Comparing Leiden clusters with key Melanoma Markers
    # AXL: Resistance marker; MITF: Differentiation marker
    sc.pl.umap(
        final_adata, 
        color=['leiden'], 
        title='Leiden Clustering of Melanoma Atlas',
        save='_clusters.png'
    ) # sc.pl: Plotting
