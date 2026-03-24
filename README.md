# Bioinformatics
Immuno-Resistance Discovery in Advanced Melanoma: A Multi-modal Integrative Framework

📝 Project Abstract
This repository implements a robust bioinformatics pipeline to identify potential immunotherapy resistance genes in advanced melanoma. The study transitions from macro-scale clinical cohort denoising to micro-scale spatial architecture resolving, bridging two distinct data modalities:

Phase I (R-based): Large-scale bulk RNA-seq analysis of 142 clinical samples (via TIDE cohorts).

Phase II (Python-based): High-resolution Spatial Transcriptomics (10x Visium) and scRNA-seq trajectory inference.

🛠️ Key Methodologies
1. Robust Denoising with Modified RUV-III PRPS(R)
Standard normalization often "squashes" critical biological signals in clinical cohorts.

Innovation: I optimized the RUV-III PRPS (Pseudo-replicate Proprietary Samples) logic to specifically isolate technical batch effects while preserving treatment-response signatures.

Tools: edgeR, RUVseq, limma.

2. Spatial Architecture & Fate Mapping (Python)
Leveraging graph-based computing to resolve the tumor immune microenvironment (TIME).

Spatial Analysis: Neighborhood graph construction and cell-cell interaction modeling using Squidpy.

Dynamics: Probabilistic trajectory inference and RNA velocity to map resistance evolution using Scanpy and scVelo.

📂 Repository Structure
Plaintext
.
├── R_Denoising/                # Bulk RNA-seq & RUV-III Optimization
│   └── ruv3_prps_logic.R       # Custom PRPS construction scripts
├── Python_Spatial/             # Spatial & Single-cell Workflows
│   ├── spatial_neighborhood.py  # Squidpy-based graph analysis
│   └── trajectory_mapping.ipynb # Fate inference notebooks
├── data/                       # Local data storage (Ignored by Git)
├── nextflow.config             # Pipeline orchestration
└── README.md

🚀 Quick Start
Bash
# Clone the repository
git clone https://github.com/YourUsername/Melanoma-Spatial-Integrated.git

# Phase 1: Run RUV-III Denoising (R)
Rscript R_Denoising/run_ruv3.R --input data/counts.csv

# Phase 2: Spatial Analysis (Python)
python Python_Spatial/spatial_neighborhood.py --adata data/spatial_data.h5ad

🔬 Biological Impact
By integrating White-box statistical modeling (RUV-III) with Manifold learning (Trajectory), this project identifies non-obvious resistance biomarkers that are often masked in standard pipelines.

👤 Contact
Phoenix Wong Master of Science (Bioinformatics) | University of Melbourne

Interests: Biological Information Flow, Category Theory in Biology, Multi-omics Integration.
