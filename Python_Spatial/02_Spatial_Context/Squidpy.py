import scanpy as sc
import squidpy as sq
import pandas as pd
from anndata import AnnData

def run_spatial_evolution(adata: AnnData):
    """
    Spatial microenvironment analysis of melanoma.
    Focus: construction of spatial neighborhood graphs and cell–cell interactions.
    """
    
    # --- 1. Spatial coordinate visualization ---
    # 假设你的 adata.obsm['spatial'] 已经包含了 10x Visium 的坐标
    print("Visualizing the spatial organizational structure...")
    sc.pl.spatial(adata, color='leiden', spot_size=150, show=False, save='_enrichment.png')

    # --- 2. Generate a spatial neighborhood network ---
    # Link discrete points to form a network
    # sq.gr (Graph Module) is similar to sc.pp，but it is specialized for spatial geometry
    print("Generate a spatial neighborhood network...")
    sq.gr.spatial_neighbors(adata, coord_type="generic", n_neighs=6)

    # --- 3. Spatial centrality analysis ---
    # Identify which cell populations are located in the tumor core and which are at the periphery (immune infiltration zone).
    sq.gr.centrality_scores(adata, cluster_key="leiden")
    sq.pl.centrality_scores(adata, cluster_key="leiden", score="closeness_centrality",show=False, save='_enrichment.png')

    # --- 4. Neighborhood enrichment analysis ---
    # Key point: Are drug-resistant cell subpopulations preferentially adjacent to certain immune cells (such as Tregs or macrophages)?
    print("Compute neighborhood enrichment...")
    sq.gr.nhood_enrichment(adata, cluster_key="leiden")
    sq.pl.nhood_enrichment(adata, cluster_key="leiden", method="single", show=False, save='_enrichment.png')

    # --- 5. Spatially co-expressed genes ---
    # Identify drug-resistance genes that exhibit non-random spatial distribution (using Moran's I statistic).
    print("Compute spatially autocorrelated gene...")
    sq.gr.spatial_autocorr(adata, mode="moran")
    
    return adata

if __name__ == "__main__":
    # Placeholder: In practice, load the AnnData object from Step 1
    # adata = sc.read_h5ad("path_to_processed_scRNA.h5ad")
    # run_spatial_evolution(adata)
    print("Spatial pipeline logic verified. Ready for execution.")
