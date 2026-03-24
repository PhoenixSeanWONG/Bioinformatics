import scanpy as sc
import gseapy as gp
import scvelo as scv
import matplotlib.pyplot as plt
from anndata import AnnData

def run_functional_dynamics(adata: AnnData) -> AnnData:
    """
    Melanoma Functional Enrichment & Fate Mapping.
    Focus: Pathway identification and cellular state transitions.
    """

    # --- 1. Functional Enrichment (GSEApy) ---
    print("Running Pathway Enrichment for Spatially Variable Genes...")
    
    # Ensure sq.gr.spatial_autocorr was run in Step 2
    if 'moranI' in adata.uns:
        svg_genes = adata.uns['moranI'].index[:100].tolist() 
        
        enr = gp.enrichr(gene_list=svg_genes,
                         gene_sets=['MSigDB_Hallmark_2020', 'KEGG_2021_Human'],
                         organism='Human',
                         outdir=None)
        
        print(f"Top Pathway identified: {enr.results.iloc[0]['Term']}")
        # Visualize and save
        enrichment_dotplot = gp.plot.dotplot(enr.results, column="Adjusted P-value", x='Gene_set', size=10, top_term=10)
        plt.savefig('enrichment_dotplot.png', bbox_inches='tight')
    else:
        print("Warning: Moran's I results not found in adata.uns. Skipping enrichment.")

    # --- 2. Trajectory Inference (PAGA) ---
    print("Computing PAGA Trajectory...")
    
    sc.tl.paga(adata, groups='leiden')
    sc.pl.paga(adata, plot_rng_seed=42, show=False, save='_paga_connectivity.png')
    
    # Map back to UMAP for a clear 'Fate Map'
    sc.tl.draw_graph(adata, init_pos='paga')
    sc.pl.draw_graph(adata, color=['leiden', 'MITF', 'AXL'], legend_loc='on data', save='_trajectory.png')

    return adata

if __name__ == "__main__":
    print("Dynamics module (Enrichment + PAGA) ready for melanoma fate mapping.")
