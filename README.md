# spatial_transcriptomics_analysis_using_Seurat

A introductory bioinformatics project implementing the Seurat spatial transcriptomics workflow, adapted from the Harvard Chan Bioinformatics Core (HBC) "[Introduction to Spatial Transcriptomics](https://hbctraining.github.io/Intro-to-spatial-transcriptomics/schedule/self-learning.html)" self-learning workshop materials.



## Experimental design (Mov10 dataset)



Two pairwise contrast are tested against the control: overexpression vs. control and knockdown vs. control.

## Repository files and folders
1. README.md (this file)
2. spatial_analysis_script.R (fully-annotated DESeq2 analysis script)
4. output (figures)
5. results (results tables)

## Workflow summary:

The `spatial_analysis_script.R` script is organized into the following sections:
 


## How to Run
 
1. Create a new project.
2. Download the spatial data into the project directory from [HBC training](https://hbctraining.github.io/Intro-to-spatial-transcriptomics/lessons/03_loading_spatial_data.html#set-up).
3. Source `spatial_analysis_script.R` and run it.
4. Export generated plots to `outputs/` and results tables to `results/`.
