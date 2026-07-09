# spatial_transcriptomics_analysis_using_Seurat

A introductory bioinformatics project implementing the Seurat spatial transcriptomics workflow, adapted from the Harvard Chan Bioinformatics Core (HBC) "[Introduction to Spatial Transcriptomics](https://hbctraining.github.io/Intro-to-spatial-transcriptomics/schedule/self-learning.html)" self-learning workshop materials.

## Overview

Visium HD generates gene expression data at very high spatial resolution (as
fine as 2um, binned here into 8um and 16um bins), preserving each
transcript's spatial location within the tissue section. This adds a
spatial dimension on top of standard scRNA-seq-style analysis: in addition
to identifying cell types and differentially expressed genes, this workflow
also detects spatially coherent tissue domains, deconvolves cell-type
mixtures within each spatial bin, identifies genes with spatial expression
patterns, and infers cell-cell signaling across the tissue architecture.
 
**Dataset used (10x Genomics Visium HD CRC dataset):**
1. `P5CRC` (Colorectal tumor tissue)
2. `P5NAT` (Normal adjacent tissue (matched)
 
Both samples are analyzed jointly for QC and clustering, then compared
directly in the differential expression, deconvolution, and spatial
domain sections.

## Repository files and folders
1. README.md (this file)
2. spatial_analysis_script.R (fully-annotated DESeq2 analysis script)
4. output (figures)
5. results (results tables)

This repository contains only the annotated analysis code, showcased as a
worked-through example of the Visium HD spatial transcriptomics workflow. It
does **not** include the raw Space Ranger output, intermediate Seurat/RDS
objects, or reference single-cell dataset, since these are large,
project-specific files not suited to version control.
 
The dataset is publicly provided by HBC as part of their training materials;
see [HBC's self-learning workshop page](https://hbctraining.github.io/Intro-to-spatial-transcriptomics/schedule/self-learning.html)
for data download instructions and the associated lesson pages. To reproduce
the analysis locally, download the Space Ranger output folders
(`P5CRC_cropped/`, `P5NAT_cropped/`), the single-cell reference dataset, and
the workshop's pre-computed intermediate `.rds` checkpoints (used at several
compute-intensive steps) into a local `data/` folder, along with an `outputs/` folder for exported figures
a `results/` folder for exported output tables.

## Workflow summary:

The `spatial_analysis_script.R` script is organized into the following sections:

1. **Loading Visium HD data and QC prep** — import Space Ranger output with
   `Load10X_Spatial()` at multiple bin resolutions (8um/16um), load and merge
   multiple samples into one Seurat object.
   
3. **Quality control** — compute and visualize per-bin QC metrics (UMI
   counts, genes detected, complexity/novelty score, mitochondrial ratio)
   both spatially and as density plots, then apply sample-specific and
   global filtering thresholds.
   
5. **Normalization and sketch downsampling** — log-normalize counts,
   identify highly variable genes, and downsample to a leverage-score-based
   representative "sketch" subset for computationally expensive steps.
   
7. **Dimensionality reduction** — standalone PCA/UMAP walkthrough on the
   sketch assay (variable feature selection, scaling, PCA, UMAP).
   
9. **Clustering** — k-nearest-neighbor graph construction and Leiden
   clustering on a downsampled "sketch" assay, then projection of cluster
   labels back onto the full dataset via `ProjectData()`.
   
11. **Integration/batch assessment** — check whether clusters separate by
   sample (batch effect) or reflect genuine tumor-vs-normal biology (e.g.
   CEACAM5/6 tumor marker expression).

12. **BANKSY spatial clustering** — augment each bin's expression profile
   with its spatial neighborhood's average expression to detect spatially
   coherent tissue domains; integrate across samples with Harmony; compare
   BANKSY-derived clusters against purely expression-derived clusters.

13. **Deconvolution (RCTD)** — estimate the cell-type composition of each
   spatial bin using a single-cell reference dataset via `spacexr::RCTD`,
   including the sketch/project-back workflow needed for large datasets.

15. **Differential expression and pathway analysis** — test for DE genes
    within a specific cell type (Myeloid) between tumor and normal tissue
    using `FindMarkers()`; visualize with volcano plots, violin plots, and
    spatial/UMAP feature plots; run GO over-representation analysis (ORA)
    and MSigDB gene set enrichment analysis (GSEA).
    
17. **Spatially variable genes** — identify genes with significant spatial
    autocorrelation (Moran's I) within the tumor region, and relate them to
    BANKSY-defined tumor sub-clusters.
    
19. **CellChat** — infer ligand-receptor cell-cell communication across
    RCTD-annotated cell types, incorporating spatial proximity constraints;
    visualize with circle plots, heatmaps, bubble plots, and chord diagrams.

## How to Run
1. Install the following packages in R or RStudio
   
       ## ---- Install package managers ----
       install.packages("remotes")
       install.packages("devtools")
       install.packages("BiocManager")
       
       ## ---- Install CRAN packages ----
       install.packages("tidyverse")   # data wrangling + ggplot2
       install.packages("scales")      # axis/label formatting for plots
       install.packages("future")      # parallelization backend
       install.packages("R.utils")     # utility functions (file handling, etc.)
       install.packages("arrow")       # fast columnar data I/O (used internally by Seurat/BPCells)
       install.packages("hdf5r")       # read/write HDF5 files (10x .h5 matrices)
       install.packages("qs2")         # fast serialization for large R objects
       install.packages("Rfast2")      # fast statistical routines (Seurat dependency)
       install.packages("leidenbase")  # Leiden clustering algorithm backend
       install.packages("harmony")     # batch correction/integration
       install.packages("Seurat")      # core single-cell/spatial analysis framework
       
       ## ---- Install Bioconductor packages ----
       BiocManager::install("EnhancedVolcano")  # publication-style volcano plots
       BiocManager::install("clusterProfiler")  # GO/GSEA enrichment testing
       BiocManager::install("org.Hs.eg.db")     # human genome-wide annotation database
       BiocManager::install("msigdbr")          # MSigDB gene sets for GSEA
       BiocManager::install("ComplexHeatmap")   # advanced heatmap visualization
       BiocManager::install("BiocNeighbors")    # nearest-neighbor search backend
       
       ## ---- Install GitHub packages ----
       remotes::install_github("prabhakarlab/Banksy")             # BANKSY spatial clustering
       remotes::install_github("satijalab/seurat-wrappers")       # Seurat wrapper functions (incl. RunBanksy)
       devtools::install_github("dmcable/spacexr",
                                build_vignettes = FALSE)           # RCTD cell-type deconvolution
       devtools::install_github("jinworks/CellChat")               # cell-cell communication inference
       devtools::install_github("immunogenomics/presto")           # fast Wilcoxon test (speeds up FindMarkers/FindAllMarkers)
   
2. Download the project directory from [HBC training](https://hbctraining.github.io/Intro-to-spatial-transcriptomics/lessons/03_loading_spatial_data.html#set-up).
3. Source `spatial_analysis_script.R` and run it.
4. Export generated plots to `outputs/` and results tables to `results/`.
