#######################################################################
# Visium HD Spatial Transcriptomics Analysis

# Adapted from the Harvard Chan Bioinformatics Core (HBC) self-learning
# workshop: "Introduction to Spatial Transcriptomics"
# https://hbctraining.github.io/Intro-to-spatial-transcriptomics/schedule/self-learning.html

# Dataset: 10x Genomics Visium HD colorectal cancer (CRC) dataset -
# a tumor sample (P5CRC) and matched normal adjacent tissue (P5NAT),
# analyzed at 8um and 16um bin resolution.

#######################################################################

# ===================== 1. LOADING VISIUM HD DATA AND QC =================

## ---- Load libraries ----
library(tidyverse)
library(Seurat)

## ---- Load a single sample into Seurat ----

## Create a Seurat object directly from 10x Space Ranger output.
## bin.size = c(8, 16) loads both the 8um and 16um binned resolutions
## as separate assays within the same object.
crc <- Load10X_Spatial(data.dir = "data/P5CRC_cropped/",
                       bin.size = c(8, 16),
                       slice = "P5CRC")

## Assess no. of bins according to bin size
## (finer 8um bins => many more "spots" than coarser 16um bins)
ncol(crc[["Spatial.008um"]])
ncol(crc[["Spatial.016um"]])

## Check metadata (per-bin QC metrics, sample ID, etc.)
crc@meta.data %>% View()

## Delete crc variable to save RAM before reloading via the multi-sample loop below
rm(crc)

## ---- Loading multiple samples into Seurat ----

## List of samples (associated with data.dir); loop below loads each one
samples <- c("P5CRC", "P5NAT")

## Empty list to fill with Seurat objects, one per sample
list_seurat <- list()

for (sample in samples) {
  
  # Path to data directory
  data_dir <- paste0("data/", sample, "_cropped")
  print(data_dir)
  
  # Create seurat object (8um bins only here) and set orig.ident to be sample
  seurat <- Load10X_Spatial(data.dir = data_dir,
                            bin.size = 8,
                            slice = sample)
  seurat$orig.ident <- sample
  
  # Store seurat object in our list, keyed by sample name
  list_seurat[[sample]] <- seurat
}

## Confirm whether both samples were successfully loaded
list_seurat

## ---- Merge both samples for joint QC ----
# Create a single Seurat object out of the two per-sample objects.
# add.cell.id prefixes barcodes with the sample name to keep them unique.
seurat_merged <- merge(x = list_seurat[["P5CRC"]],
                       y = list_seurat[["P5NAT"]],
                       add.cell.id = c("P5CRC", "P5NAT"))
seurat_merged

## Join layers to collapse the per-sample count layers into a single
## unified counts matrix (needed before most downstream Seurat functions)
seurat_merged <- JoinLayers(seurat_merged)

## ---- Evaluate merged_seurat ----

## Sum of cells/bins in P5CRC and P5NAT in the pre-merge list (sanity check)
ncol(list_seurat[["P5CRC"]]) + ncol(list_seurat[["P5NAT"]])

## Check the total number of bins in the merged seurat object matches
seurat_merged

## Check that the merged object has the appropriate sample-specific barcode prefixes
seurat_merged@meta.data %>% head()

## Set Idents to sample IDs so downstream comparisons default to grouping by sample
Idents(seurat_merged) <- "orig.ident"

## Save merged Seurat object so QC/downstream steps can restart from here
saveRDS(seurat_merged, "data/01_seurat_merged.RDS")

## Clear everything except the RDS object to save RAM
rm(list = setdiff(ls(), "seurat_merged"))

# ============================= 2. QUALITY CONTROL =========================

## ---- Number of bins before filtration ----

## Barplot: number of bins captured per sample, before any QC filtering
ggplot(seurat_merged@meta.data) +
  geom_bar(aes(x = orig.ident, fill = orig.ident),
           color = "black") +
  geom_text(aes(x = orig.ident, label=after_stat(count)), 
            stat='count', vjust=-1) +
  theme_classic()
## Save to "figures/"

## ---- Quality metrics ----

## Store metadata as a variable for easier repeated access/plotting
meta <- seurat_merged@meta.data
View(meta)

## ---- UMI counts (transcripts) per bin ----

## Spatial overlay: visualize the spatial distribution of total UMIs
## across the tissue section (low-UMI bins often correspond to background/
## non-tissue areas or poor-quality regions)
SpatialFeaturePlot(seurat_merged, 
                   "nCount_Spatial.008um",
                   pt.size.factor = 15,
                   image.alpha = 0,
                   max.cutoff = "q90")
## Save to "figures/"

## Before filtration density: log10-transformed density of UMIs per sample.
## Vertical lines mark the sample-specific filtering thresholds chosen below.
ggplot(meta) +
  geom_density(aes(x = nCount_Spatial.008um, fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 30, color = "pink") +
  geom_vline(xintercept = 10, color = "lightblue") +
  scale_x_log10() +
  theme_classic()
## Save to "figures/"

## After filtration density: apply the chosen thresholds and re-plot to
## confirm the low-count tail (likely empty/background bins) is removed
meta_filt <- subset(meta,
                    ((orig.ident == "P5CRC") & (nCount_Spatial.008um > 30)) |
                      ((orig.ident == "P5NAT") & (nCount_Spatial.008um > 10)))

ggplot(meta_filt) +
  geom_density(aes(x = nCount_Spatial.008um,
                   fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 30, color = "pink") +
  geom_vline(xintercept = 10, color = "lightblue") +
  scale_x_log10() +
  theme_classic()
## Save to "figures/"

## ---- Genes detected per bin ----

## Spatial overlay: distribution of the number of unique genes detected per bin
SpatialFeaturePlot(seurat_merged, 
                   "nFeature_Spatial.008um",
                   pt.size.factor = 15,
                   image.alpha = 0,
                   max.cutoff = "q90")
## Save to "figures/"

## Before filtration density (same filtering thresholds shown for reference)
ggplot(meta) +
  geom_density(aes(x = nFeature_Spatial.008um,
                   fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 30, color = "pink") +
  geom_vline(xintercept = 10, color = "lightblue") +
  scale_x_log10() +
  theme_classic()
## Save to "figures/"

## After filtration density (repeated view of the same metric post-threshold)
ggplot(meta) +
  geom_density(aes(x = nFeature_Spatial.008um,
                   fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 30, color = "pink") +
  geom_vline(xintercept = 10, color = "lightblue") +
  scale_x_log10() +
  theme_classic()
## Save to "figures/"

## ---- Complexity (novelty) score ----
## log10(genes)/log10(UMIs) - flags bins with unusually low gene diversity
## relative to their sequencing depth (a sign of low-complexity/ambient RNA)

## Add number of genes per UMI for each bin to metadata
seurat_merged$log10GenesPerUMI <- log10(seurat_merged$nFeature_Spatial.008um) / 
  log10(seurat_merged$nCount_Spatial.008um)

## Turn NA values (from bins with 0 counts) into 0 for now
seurat_merged$log10GenesPerUMI[is.na(seurat_merged$log10GenesPerUMI)] <- 0

## Spatial overlay: visualize the spatial distribution of the complexity score
SpatialFeaturePlot(seurat_merged, 
                   "log10GenesPerUMI",
                   pt.size.factor = 17,
                   image.alpha = 0,
                   min.cutoff = "q10")
## Save to "figures/"

## Before filtration density: distribution of complexity score, with the
## chosen filtering threshold (0.80) marked
meta <- seurat_merged@meta.data
ggplot(meta) +
  geom_density(aes(x = log10GenesPerUMI,
                   fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 0.80) +
  theme_classic()
## Save to "figures/"

## After filtration density (repeated view of the same metric/threshold)
meta <- seurat_merged@meta.data
ggplot(meta) +
  geom_density(aes(x = log10GenesPerUMI,
                   fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 0.80) +
  theme_classic()
## Save to "figures/"

## ---- Mitochondrial counts ratio ----
## High mitochondrial content is a classic marker of stressed/dying cells
## or ambient RNA contamination.

## Compute percent mito ratio by finding genes that start with "MT-"
seurat_merged$mitoRatio <- PercentageFeatureSet(object = seurat_merged, 
                                                pattern = "^MT-")
seurat_merged$mitoRatio <- seurat_merged@meta.data$mitoRatio / 100

## Turn NA values into 1.00 (treat unmeasurable bins as maximally "bad") for now
seurat_merged$mitoRatio[is.na(seurat_merged$mitoRatio)] <- 1.00

## Spatial overlay: visualize the spatial distribution of mitochondrial ratio
SpatialFeaturePlot(seurat_merged, 
                   "mitoRatio",
                   pt.size.factor = 15,
                   image.alpha = 0,
                   max.cutoff = "q90")
## Save to "figures/"

## Before filtration density: distribution of mito ratio, with the chosen
## filtering threshold (0.25) marked
meta <- seurat_merged@meta.data

ggplot(meta) +
  geom_density(aes(x = mitoRatio,
                   fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 0.25) +
  theme_classic()
## Save to "figures/"

## After filtration density: apply the mito ratio threshold and re-plot
meta_filt <- subset(meta, mitoRatio < 0.25)

ggplot(meta_filt) +
  geom_density(aes(x = mitoRatio,
                   fill = orig.ident),
               alpha = 0.4,
               color = "black") +
  geom_vline(xintercept = 0.25) +
  theme_classic()
## Save to "figures/"

## ---- Filtration ----
## Apply all thresholds explored above sequentially to produce the final,
## analysis-ready filtered object.

## Per-sample nCount thresholds (different depth cutoffs per sample)
seurat_filtered <- subset(seurat_merged,
                          ((orig.ident == "P5CRC") & (nCount_Spatial.008um > 30)) |
                            ((orig.ident == "P5NAT") & (nCount_Spatial.008um > 10)))

## Per-sample nFeature thresholds
seurat_filtered <- subset(seurat_filtered,
                          ((orig.ident == "P5CRC") & (nFeature_Spatial.008um > 30)) |
                            ((orig.ident == "P5NAT") & (nFeature_Spatial.008um > 10)))

## Global thresholds for mitochondrial ratio and complexity (applied to all samples)
seurat_filtered <- subset(seurat_filtered, mitoRatio < 0.25)
seurat_filtered <- subset(seurat_filtered, log10GenesPerUMI > 0.80)

## Print seurat object after filtration (confirm reduced bin count)
seurat_filtered

## ---- Visualising counts data (post-filtering sanity check) ----

## Violin plot of UMIs per sample after filtering
p_ncount <- VlnPlot(seurat_filtered, 
                    features = "nCount_Spatial.008um", 
                    pt.size = 0, group.by = 'orig.ident') +
  NoLegend()

## Violin plot of number of genes per sample after filtering
p_nfeats <- VlnPlot(seurat_filtered, 
                    features = "nFeature_Spatial.008um", 
                    pt.size = 0, group.by = 'orig.ident') + 
  NoLegend()

## Plot UMIs and gene count violin plots side-by-side
p_ncount | p_nfeats
## Save to "figures/"

## ---- Spatial overlay (post-filtering) ----

## Visualize the spatial distribution of total UMIs and number of genes
## after filtration - confirm removed bins were background, not tissue
SpatialFeaturePlot(seurat_filtered, 
                   c("nFeature_Spatial.008um", 
                     "nCount_Spatial.008um"),
                   pt.size.factor = 16,
                   image.alpha = 0)
## Save to "figures/"

## Save filtered Seurat object as the checkpoint for downstream analysis
saveRDS(seurat_filtered, "data/seurat_filtered.RDS")

## Clear everything except the RDS object to save RAM
rm(list = setdiff(ls(), "seurat_filtered"))

# ============= 3. NORMALIZING AND SKETCH DOWNSAMPLING ======================
# Visium HD 8um data can contain hundreds of thousands of bins - "sketching"
# selects a representative, information-dense subset for computationally
# expensive steps (PCA/clustering), whose results are later projected back.

library(Seurat)
library(tidyverse)
library(scales)
set.seed(12345)

## Load filtered dataset (checkpoint from Section 2)
seurat_filtered <- readRDS("data/seurat_filtered.RDS")

## Log-normalize the raw counts (standard library-size + log1p normalization)
seurat_sketch <- NormalizeData(seurat_filtered)
seurat_sketch

## Identify the most variable genes (used for PCA feature selection)
seurat_sketch <- FindVariableFeatures(seurat_sketch, 
                                      selection.method = "vst",
                                      nfeatures = 2000)
seurat_sketch

## Identify the 15 most highly variable genes, just to inspect them
ranked_variable_genes <- VariableFeatures(seurat_sketch)
top_genes <- ranked_variable_genes[1:15]
top_genes

## ---- Downsampling ----

## Select 10,000 bins based on leverage scores (an importance-sampling
## metric that favors information-rich/distinct bins over redundant ones),
## storing them in a new "sketch" assay
seurat_sketch <- SketchData(object = seurat_sketch,
                            assay = 'Spatial.008um',
                            ncells = 10000,
                            method = "LeverageScore",
                            sketched.assay = "sketch")

## Seurat object summary after sketch downsampling
seurat_sketch
seurat_sketch@meta.data %>% View()

## ---- Visualising leverage score ----

## Histogram of leverage scores, split by sample
ggplot(seurat_sketch@meta.data) +
  geom_histogram(aes(x=leverage.score, 
                     fill=orig.ident), 
                 alpha=0.5, bins=100) +
  theme_classic() +
  scale_x_log10(labels = scales::label_number())
## Save to "figures/"

## Spatial feature plot: where on the tissue were high-leverage bins selected?
SpatialFeaturePlot(seurat_sketch,
                   "leverage.score",
                   pt.size.factor = 15,
                   image.alpha = 0,
                   max.cutoff = 2)

## Save Seurat object with sketch assay for use in clustering and BANKSY
saveRDS(seurat_sketch, "data/seurat_sketch.RDS")

## Clear everything except the RDS object to save RAM
rm(list = setdiff(ls(), "seurat_sketch"))

# ==================== 4. DIMENSIONALITY REDUCTION ==========================

library(Seurat)
library(tidyverse)
set.seed(12345)

seurat_sketch <- readRDS("data/seurat_sketch.RDS")

## Identify the 15 most highly variable genes (from the full assay)
ranked_variable_genes <- VariableFeatures(seurat_sketch, 
                                          assay = "Spatial.008um")
top_genes <- ranked_variable_genes[1:15]
top_genes

## Identify the most variable genes within the sketch assay specifically
seurat_processed <- FindVariableFeatures(seurat_sketch, 
                                         selection.method = "vst",
                                         nfeatures = 2000,
                                         assay = "sketch")
seurat_processed

## Identify the 15 most highly variable genes (sketch assay)
ranked_variable_genes <- VariableFeatures(seurat_processed,
                                          assay = "sketch")
top_genes <- ranked_variable_genes[1:15]

## Plot the average expression and variance of variable genes, labeling
## which genes fall in the top 15 (mean-variance/vst diagnostic plot)
p <- VariableFeaturePlot(seurat_processed,
                         assay = "sketch")

## Check for zero-mean genes
hvf_info <- HVFInfo(seurat_processed, assay = "sketch")
## If >0, drop repel

LabelPoints(plot = p, points = top_genes, repel = FALSE)
## Save to "figures/"

## ---- PCA ----

## Scale the log-normalized expression (zero-center/unit-variance per gene,
## required before PCA so highly-expressed genes don't dominate the components)
seurat_processed <- ScaleData(seurat_processed)
seurat_processed

## Run PCA on the sketch assay
seurat_processed <- RunPCA(seurat_processed,
                           assay = "sketch",
                           reduction.name = "pca.sketch")
seurat_processed

## Plot PCA colored by sample
DimPlot(seurat_processed,
        group.by = "orig.ident",
        reduction = "pca.sketch")

## ---- UMAP ----

## Run UMAP on the PCA embedding; return.model = TRUE stores the UMAP model
## itself, which is required later to project full-dataset bins into the
## same UMAP space
seurat_processed <- RunUMAP(seurat_processed, 
                            reduction = "pca.sketch", 
                            reduction.name = "umap.sketch", 
                            return.model = TRUE, 
                            dims = 1:50,
                            seed.use = 12345)
seurat_processed

## Plot UMAP colored by sample
DimPlot(seurat_processed, 
        group.by = "orig.ident",
        reduction = "umap.sketch")
## UMAP has no fixed coordinate system. The algorithm optimizes
## relative distances between points, not absolute x/y position, so the output
## can come out rotated, flipped, or mirrored on any run
## Save to "figures/"

## View all dimensionality reduction results stored in the object
seurat_processed@reductions

# ============================= 5. CLUSTERING =============================
# This section runs FindNeighbors/FindClusters on the "sketch" assay

## ---- k-nearest neighbours (kNN) ----

## Determine the K-nearest neighbor graph in PCA space, using the sketch
## (downsampled) assay/reduction for computational efficiency
seurat_processed <- FindNeighbors(seurat_processed, 
                                  assay = "sketch", 
                                  reduction = "pca.sketch",
                                  dims = 1:50)

seurat_processed@graphs

## ---- Clustering ----
## algorithm = 4 specifies the Leiden algorithm (higher quality, more
## computationally intensive than default Louvain)
seurat_processed <- FindClusters(seurat_processed, 
                                 cluster.name = "seurat_cluster.sketched", 
                                 algorithm = 4,
                                 resolution = 0.65,
                                 random.seed = 12345)

seurat_processed@meta.data %>% 
  head() %>% 
  relocate("seurat_cluster.sketched")
  
## Check updated idents after running FindClusters()
  Idents(seurat_processed) %>% head()
## Clustering can be variable as sketching is randomized
  
## Visualize UMAP with cluster labels overlaid
p <- DimPlot(seurat_processed, 
             reduction = "umap.sketch", label = T) + 
  ggtitle("Sketched clustering")
LabelClusters(p, id = "ident",  fontface = "bold", size = 5, 
              bg.colour = "white", bg.r = .2, force = 0)

## ---- Project back to the entire dataset ----
## Clustering was done on the downsampled "sketch" bins for speed;
## ProjectData() propagates those cluster labels back onto all bins
## in the full dataset via the shared PCA space.

## Increases the size of the default vector (8GB) - large full-dataset
## projections can otherwise exceed R's default in-memory object size limit
options(future.globals.maxSize = 8000 * 1024^2)

## Project from sketch assay onto entire dataset
seurat_processed <- ProjectData(
  object = seurat_processed,
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched"))
relocate("seurat_cluster.sketched")

seurat_processed@reductions

seurat_processed@meta.data %>% head()

## ---- Visualize projected clusters ----

## Switch to full (non-sketch) assay for downstream full-dataset analysis
DefaultAssay(seurat_processed) <- "Spatial.008um"

## Arrange seurat_cluster.projected so that values are in numeric order
cluster_order <- seurat_processed$seurat_cluster.projected %>%
  unique() %>% as.numeric() %>%
  sort() %>% as.character()
cluster_order

## Factorize seurat_cluster.projected with levels in correct numeric order
seurat_processed$seurat_cluster.projected <- seurat_processed$seurat_cluster.projected %>% 
  factor(levels = cluster_order)

## Print head of seurat_cluster.projected to confirm levels are in correct order
seurat_processed$seurat_cluster.projected %>% head()

## Change the idents to the projected cluster assignments
Idents(seurat_processed) <- "seurat_cluster.projected"

## Plot UMAP now projected across ALL bins
p <- DimPlot(seurat_processed, 
             reduction = "full.umap.sketch", label = T) + 
  ggtitle("Projected clustering")
LabelClusters(p, id = "ident",  fontface = "bold", size = 5, 
              bg.colour = "white", bg.r = .2, force = 0)
## Save to "figures/"

## Overlay projected clusters directly onto the tissue image
SpatialDimPlot(seurat_processed,
               pt.size.factor = 16,
               image.alpha = 0)

## ---- Assessing cluster quality ----
## Check known endothelial marker genes (PECAM1/VWF) to confirm at least
## one cluster corresponds to a biologically sensible cell population
VlnPlot(seurat_processed,
        features = c("PECAM1", "VWF"),
        pt.size = 0.1,
        ncol = 1) +
  NoLegend()

## Save Seurat object with clustering results
saveRDS(seurat_processed, "data/seurat_processed.RDS")

## Clear everything except the RDS object to save RAM
rm(list = setdiff(ls(), "seurat_processed"))

# ============================= 6. INTEGRATION ============================
# (Batch assessment step - checking whether clusters are driven by real
# biology or by sample-of-origin/technical batch effects.)

library(Seurat)
library(tidyverse)

seurat_processed <- readRDS("data/seurat_processed.RDS")

## Visualize PCA and UMAP colored by sample, to check for batch separation
p1 <- DimPlot(seurat_processed,
              group.by="orig.ident",
              reduction = "full.pca.sketch") +
  NoLegend() + ggtitle("PCA (Projected)")

p2 <- DimPlot(seurat_processed,
              group.by="orig.ident",
              reduction = "full.umap.sketch") +
  ggtitle("UMAP (Projected)")

p1 | p2
## Save this figure to "outputs/"

## ---- Distribution of sample per cluster ----

## Barplot of the proportion of each cluster's bins coming from each sample -
## clusters composed almost entirely of one sample suggest a batch effect
## (or a genuinely sample-specific biological population, e.g. tumor cells)
ggplot(seurat_processed@meta.data) +
  geom_bar(aes(x = seurat_cluster.projected, 
               fill = orig.ident), 
           position = position_fill())  +
  theme_classic()
## Save this figure to "outputs/"

## ---- Possible explanations for batch (checking tumor marker genes) ----

## CEACAM5: check if the sample-dominated clusters express a known CRC tumor marker
p_feature <- FeaturePlot(seurat_processed, 
                         features = "CEACAM5",
                         label = TRUE)

p_vln <- VlnPlot(seurat_processed, 
                 features = "CEACAM5", 
                 pt.size = 0) +
  NoLegend()

p_feature | p_vln
## Save this figure to "outputs/"

## CEACAM6: a second CRC tumor marker, for corroboration
p_feature <- FeaturePlot(seurat_processed, 
                         features = "CEACAM6",
                         label = TRUE)

p_vln <- VlnPlot(seurat_processed, 
                 features = "CEACAM6", 
                 pt.size = 0) +
  NoLegend()

p_feature | p_vln
## Save this figure to "outputs/"

## Visualize expression of both tumor markers together across clusters
DotPlot(seurat_processed, 
        features = c("CEACAM5", "CEACAM6"))

# ============================= 7. BANKSY ===================================
# BANKSY augments each bin's own expression profile with its neighborhood's
# average expression, enabling detection of spatially coherent tissue
# domains rather than purely transcriptionally-defined cell clusters.

library(Seurat)
library(Banksy)
library(SeuratWrappers)
library(tidyverse)
library(future)
library(harmony)
library(patchwork)
set.seed(12345)

## Increases the size of the default vector (8GB)
options(future.globals.maxSize = 8000 * 1024^2)
## Set parallelization (multithreading) to speed up calculations
library(future)
plan(multisession, workers = parallel::detectCores() - 1)

## Load dataset (clustered/projected checkpoint)
seurat_processed <- readRDS("data/seurat_processed.RDS")

## Store the cell barcode as a column (rownames can get dropped during joins)
seurat_processed$cell <- rownames(seurat_processed@meta.data)

## Extract spatial (x,y) coordinates per sample image
coords_crc <- GetTissueCoordinates(seurat_processed, 
                                   image = "P5CRC.008um")
coords_nat <- GetTissueCoordinates(seurat_processed, 
                                   image = "P5NAT.008um")

## Merge P5NAT and P5CRC coordinates into one dataframe
coords <- rbind(coords_crc, coords_nat)

## Join spatial coordinates into the metadata (required by RunBanksy)
seurat_processed@meta.data <- left_join(x = seurat_processed@meta.data, 
                                        y = coords, 
                                        by = "cell")
seurat_processed@meta.data %>% head()

## (!) laptop may have limited RAM for this step
## load in processed data at the end of this checkpoint

## ---- BANKSY Checkpoint (! time-consuming) (Start) ----

## Run BANKSY: lambda controls the neighborhood-vs-own-expression weighting
## (0.8 favors spatial domain detection over pure cell typing);
## k_geom sets the neighborhood size; group ensures neighbors are only
## computed within the same sample (avoids cross-sample "neighbors")
seurat_banksy <- RunBanksy(seurat_processed, 
                           lambda = 0.8, 
                           k_geom = 15, 
                           use_agf = TRUE,
                           dimx = 'x',
                           dimy = 'y', 
                           group = "orig.ident", 
                           split.scale = TRUE,
                           assay = 'Spatial.008um', 
                           slot = 'data',
                           verbose = TRUE)

seurat_banksy

## ---- BANKSY matrix ----
## The BANKSY assay concatenates own-expression + neighborhood-averaged
## expression features, roughly doubling (or more) the feature count

## First features in BANKSY matrix
Features(seurat_banksy)[1:6]

## Middle features in BANKSY matrix
Features(seurat_banksy)[2000:2006]

## Last features in BANKSY matrix
Features(seurat_banksy)[4000:4006]

## ---- Calculations on BANKSY matrix ----

## ---- Dimensionality reduction ----

## Run PCA directly on the BANKSY feature matrix
seurat_banksy <- RunPCA(seurat_banksy, 
                        assay = "BANKSY", 
                        reduction.name = "pca.banksy", 
                        features = Features(seurat_banksy), 
                        npcs = 30)
Features(seurat_banksy)[2000:2006]

## BANKSY PCA colored by sample (check for batch separation before integration)
p_orig <- DimPlot(seurat_banksy,
                  group.by = "orig.ident",
                  reduction = "pca.banksy") 

## BANKSY PCA colored by IGKC expression (a plasma cell marker, used as a
## biological reference point in PCA space)
p_igkc <- FeaturePlot(seurat_banksy,
                      feature = "IGKC",
                      reduction = "pca.banksy")

## Plot PCA space side-by-side
p_orig + p_igkc

## ---- Integration (Using Harmony) ----

## Harmony batch correction: removes sample-of-origin effects from the
## BANKSY PCA embedding while preserving genuine biological variation
seurat_banksy <- RunHarmony(seurat_banksy,
                            group.by.vars = "orig.ident",
                            reduction = "pca.banksy",
                            reduction.save = "harmony.banksy",
                            assay.use = "BANKSY")

## Harmony BANKSY PCA colored by sample (should show better mixing than before)
p_orig <- DimPlot(seurat_banksy,
                  group.by = "orig.ident",
                  reduction = "harmony.banksy") 

## Harmony BANKSY PCA colored by IGKC expression
p_igkc <- FeaturePlot(seurat_banksy,
                      feature = "IGKC",
                      reduction = "harmony.banksy")

## Plot latent space side-by-side
p_orig + p_igkc

## ---- Clustering ----

## Find k-Nearest Neighbors on the batch-corrected Harmony/BANKSY embedding
seurat_banksy <- FindNeighbors(seurat_banksy, 
                               reduction = "harmony.banksy", 
                               dims = 1:30)

## Leiden clustering on the spatial-domain-aware BANKSY embedding
seurat_banksy <- FindClusters(seurat_banksy, 
                              algorithm = 4,
                              cluster.name = "banksy_cluster",
                              resolution = 0.65,
                              random.seed = 12345)

## ---- Delete BANKSY assay to conserve memory ----
## (the augmented BANKSY feature matrix is large and no longer needed
## once the embedding/clusters have been computed)

## Change default assay away from BANKSY before removing it
DefaultAssay(seurat_banksy) <- "Spatial.008um"

## Delete BANKSY assay by setting it to NULL
seurat_banksy[["BANKSY"]] <- NULL

## ---- BANKSY Checkpoint (End) ----

## Load pre-generated seurat object with relevant BANKSY results
## (checkpoint provided for workshop reproducibility/speed)
seurat_banksy <- readRDS("data/intermediate_seurat/seurat_banksy.rds")

## Clear everything except the RDS object to save RAM
rm(list = setdiff(ls(), "seurat_banksy"))

## Barplot of proportion of bins in each BANKSY cluster by sample
ggplot(seurat_banksy@meta.data) +
  geom_bar(aes(x = banksy_cluster, 
               fill = orig.ident), 
           position = position_fill())  +
  theme_classic()
## Save to "figures/"

## ---- Comparing clustering results ----

## Spatial overlay: visualize BANKSY-derived clusters alongside the
## purely expression-derived clusters, side-by-side
SpatialDimPlot(seurat_banksy, 
               group.by = c("banksy_cluster",
                            "seurat_cluster.projected"), 
               pt.size.factor = 15, 
               image.alpha = 0)
## Save to "figures/"

## Highlight expression-based cluster 9 on the tissue image
Idents(seurat_banksy) <- "seurat_cluster.projected"
SpatialDimPlot(seurat_banksy,
               pt.size.factor = 15,
               cells.highlight = WhichCells(seurat_banksy,
                                            idents = 9)) +
  plot_annotation(title = "Seurat cluster 9 (projected)")
## Save to "figures/"

## Contrast with BANKSY cluster 12 highlighted on the same tissue image -
## illustrates how spatial-domain clusters can differ from expression-only clusters
Idents(seurat_banksy) <- "banksy_cluster"
SpatialDimPlot(seurat_banksy,
               pt.size.factor = 15,
               cells.highlight = WhichCells(seurat_banksy,
                                            idents = 12)) +
  plot_annotation(title = "BANKSY cluster 12")
## Save to "figures/"

## ---- UMAP comparison ----

## Compare the UMAPs of the BANKSY-derived clusters and expression-derived clusters
DimPlot(seurat_banksy, 
        group.by = c("banksy_cluster",
                     "seurat_cluster.projected"),
        reduction = "full.umap.sketch",
        label = TRUE)

## ---- Assess clustering ----

## Violin plots for two Smooth Muscle marker genes, to confirm a BANKSY
## cluster corresponds to a biologically sensible tissue domain
VlnPlot(seurat_banksy,
        features = c("ACTA2", "MYH11"),
        pt.size = 0,
        ncol = 1) +
  NoLegend()

## Save seurat_banksy checkpoint
saveRDS(seurat_banksy, "data/seurat_banksy.RDS")

## Clear everything except the RDS object to save RAM
rm(list = setdiff(ls(), "seurat_banksy"))

# ============================= 7. DECONVOLUTION ============================
# RCTD (Robust Cell Type Decomposition) estimates the cell-type composition
# of each spatial bin using a reference single-cell/single-nucleus dataset
# with known cell-type labels.

library(Seurat)
library(tidyverse)
library(spacexr) # Contains RCTD function

## Set parallelization (multithreading) to speed up calculations
plan("multicore", workers = parallel::detectCores() - 1)

## Increases the size of the default vector (8GB)
options(future.globals.maxSize = 8000 * 1024^2)

## Load dataset (BANKSY-clustered checkpoint)
seurat_banksy <- readRDS("data/intermediate_seurat/seurat_banksy.rds")

## ---- Reference dataset ----

## Load Seurat reference dataset with known cell-type annotations
seurat_ref <- readRDS("data/crc_flex_ref_downsample.RDS")
seurat_ref

## Set idents to the reference cell-type labels
Idents(seurat_ref) <- "Level1"

## UMAP of reference dataset, labeled by cell type
p <- DimPlot(seurat_ref) +
  ggtitle("RCTD Reference Dataset")
### Add cluster labels
LabelClusters(p,
              id = "ident",
              fontface = "bold",
              size = 3, 
              bg.colour = "white",
              bg.r = .2,
              force = 0)
## Save to "figures/"

## ---- Create reference ----

## Raw counts of reference dataset (RCTD requires raw, not normalized, counts)
counts <- seurat_ref[["RNA"]]$counts

## Cell type annotation of reference dataset (must be a factor for RCTD)
cluster <- seurat_ref$Level1
cluster <- as.factor(cluster)

## Total counts (library size) of each reference cell
nUMI <- seurat_ref$nCount_RNA

## Create the RCTD reference object bundling counts, labels, and depth
reference <- Reference(counts = counts, 
                       cell_types = cluster, 
                       nUMI = nUMI)

## ---- Create query (the spatial data to be deconvolved) ----

## Subset to P5CRC sample only (RCTD run per-sample)
crc <- subset(seurat_banksy, 
              subset = (orig.ident == "P5CRC"))

crc

## Create a sketch of the CRC subset (RCTD is run on a downsampled query
## for speed, then results are projected back to all bins)
DefaultAssay(crc) <- 'Spatial.008um'

## HVGs are required input for SketchData
crc <- FindVariableFeatures(object = crc,
                            assay = "Spatial.008um")

## Calculate leverage score and downsample to 5,000 representative bins
crc <- SketchData(
  object = crc,
  ncells = 5000,
  method = "LeverageScore",
  sketched.assay = "sketch_crc",
  var.name = "leverage.score_crc")

## Run through the standard processing workflow, stopping at PCA
## (only need a PCA embedding for the later ProjectData step)
crc <- FindVariableFeatures(crc)
crc <- ScaleData(crc)
crc <- RunPCA(crc,
              reduction.name = "pca.sketch_crc")

## Get raw counts of the downsampled sketch assay (RCTD query input)
DefaultAssay(crc) <- "sketch_crc"
counts_sketch <- crc[["sketch_crc"]]$counts

## Grab barcodes for the sketched bins
cells_sketch <- colnames(crc[["sketch_crc"]])

## Spatial coordinates for each bin in the sketch assay
coords <- GetTissueCoordinates(crc)[cells_sketch, c("x", "y")]

## Grab total UMI counts per sketch bin
nUMI_sketch <- crc@meta.data[cells_sketch, "nCounts_Spatial.008um"]

## Create the RCTD query object bundling coords, counts, and depth
query <- SpatialRNA(coords = coords, 
                    counts = counts_sketch, 
                    nUMI = nUMI_sketch)

## ---- Run RCTD ----

## Set up RCTD by supplying both the query (spatial) and reference (single-cell) data
RCTD <- create.RCTD(spatialRNA = query, 
                    reference = reference, 
                    max_cores = parallel::detectCores() - 1)

## Run deconvolution. doublet_mode = "doublet" allows each bin to be
## assigned up to 2 cell types (appropriate given bins can capture multiple cells)
## This step may take 5-10 minutes to run
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")

## ---- RCTD results ----

## Grab results dataframe and view first few rows
df_rctd_results <- RCTD@results$results_df
df_rctd_results %>%
  head()

## Barplot: number of bins annotated to each cell type in the sketch assay
ggplot(df_rctd_results,
       aes(x = first_type, 
           fill = first_type)) +
  geom_bar() +
  theme_bw() + NoLegend() +
  ggtitle("P5CRC: RCTD First Type Results") +
  theme(axis.text.x = element_text(angle = 90,
                                   vjust = 0.5,
                                   hjust=1))

## ---- Project RCTD results onto all bins ----
## (RCTD was only run on the downsampled sketch bins; ProjectData()
## propagates the cell-type calls to the full dataset via shared PCA space)
nrow(df_rctd_results)

## Rename columns with a "_sketch" suffix before merging, to distinguish
## sketch-only results from the eventual full-dataset projected columns
colnames(df_rctd_results) <- paste0(colnames(df_rctd_results), 
                                    "_sketch")

## Create a "cell" column to merge on, converting factor columns to character
## for easier merging/joining downstream
df_rctd_results <- df_rctd_results %>%
  mutate(cell = rownames(df_rctd_results)) %>%
  mutate(first_type_sketch = as.character(first_type_sketch),
         second_type_sketch = as.character(second_type_sketch),
         spot_class_sketch = as.character(spot_class_sketch))

## Make sure columns were generated as expected
colnames(df_rctd_results)

## Create cell barcode column in metadata to join on
crc$cell <- rownames(crc@meta.data)
meta <- crc@meta.data

## Merge together RCTD results and metadata
meta <- left_join(x = meta, 
                  y = df_rctd_results,
                  by = "cell")
rownames(meta) <- meta$cell

## Put updated metadata back into the crc Seurat object
crc@meta.data <- meta

## Set value of "unassigned" for sketch bins that lack RCTD results
## (this is required because ProjectData will otherwise error on NAs)
cols <- c("spot_class_sketch", 
          "first_type_sketch", 
          "second_type_sketch")

## Iterate through each column, replacing NA sketch-bin values with "unassigned"
for (col in cols) {
  # Find rows that are in cells_sketch and are NA
  na_rows <- cells_sketch[is.na(crc@meta.data[cells_sketch, col])]
  # Assign the value "unassigned"
  crc@meta.data[na_rows, col] <- "unassigned"
}

## Re-run ProjectData now that all sketch bins have non-NA values -
## this call succeeds and propagates cell-type calls to the full dataset
crc <- ProjectData(
  object = crc,
  assay = "Spatial.008um",
  sketched.assay = "sketch_crc",
  full.reduction = "pca_crc",
  sketched.reduction = "pca.sketch_crc",
  dims = 1:30,
  refdata = list(spot_class = "spot_class_sketch",
                 first_type = "first_type_sketch",
                 second_type = "second_type_sketch"))

## ---- Visualize RCTD results ----

## Barplot of first_type labels after projecting to the full dataset
ggplot(crc@meta.data,
       aes(x = first_type, fill = first_type)) +
  geom_bar() +
  theme_bw() + NoLegend() +
  ggtitle("P5CRC: Projected First Type RCTD") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
## Save to "figures/"

## Spatial plot of first_type labels across the full tissue image
SpatialDimPlot(crc,
               group.by = "first_type",
               pt.size.factor = 15,
               image.alpha = 0.5)
## Save to "figures/"

## UMAP colored by projected RCTD cell type
Idents(crc) <- "first_type"
p <- DimPlot(crc,
             reduction = "full.umap.sketch") +
  ggtitle("P5CRC RCTD First Type")
LabelClusters(p, id = "ident",  fontface = "bold", size = 3, 
              bg.colour = "white", bg.r = .2, force = 0)
## Save to "figures/"

## ---- P5CRC and P5NAT combined results ----

## Load pre-computed Seurat object with cell-type annotations for BOTH samples
## (checkpoint provided since running RCTD twice, once per sample, is time-consuming)
seurat_rctd <- readRDS("data/intermediate_seurat/seurat_RCTD.rds")

## Spatial plot of first_type labels for both samples
SpatialDimPlot(seurat_rctd,
               group.by = "first_type",
               pt.size.factor = 15,
               image.alpha = 0)
## Save to "figures/"

# ==================== 8. DGE AND PATHWAY ANALYSIS ===========================
# Now that bins have cell-type labels (from RCTD), test for differential
# gene expression between conditions WITHIN a specific cell type
# (here: Myeloid cells, comparing tumor vs. normal adjacent tissue).

library(Seurat)
library(tidyverse)
library(EnhancedVolcano)

## ORA/GSEA libraries
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

set.seed(12345)
## Load dataset (cell-type-annotated checkpoint)
seurat_rctd <- readRDS("data/intermediate_seurat/seurat_RCTD.rds")

## Clear everything except the RDS object to save RAM
rm(list = setdiff(ls(), "seurat_rctd"))

## Subset to Myeloid cells only, across both samples
seurat_myeloid <- subset(seurat_rctd,
                         first_type == "Myeloid")
Idents(seurat_myeloid) <- "orig.ident"

## ---- FindMarkers ----

## Calculate differentially expressed genes with a Wilcoxon rank-sum test,
## comparing Myeloid cells in the tumor (P5CRC) vs. normal tissue (P5NAT)
dge <- FindMarkers(seurat_myeloid,
                   ident.1 = "P5CRC",
                   ident.2 = "P5NAT")
dge$gene <- rownames(dge)
dge %>% head()

## ---- Visualizing results ----

## Volcano plot: log2FC vs. adjusted p-value for all tested genes
EnhancedVolcano(dge,
                lab = dge$gene,
                x = "avg_log2FC",
                y = "p_val_adj",
                title = "FindMarkers Myeloid Results",
                subtitle = "CRC vs NAT")
## Save to "figures/"

## Filter to only the significant genes (padj < 0.05)
dge_sig <- dge %>% subset(p_val_adj < 0.05)

## Get the gene names for the top 6 significant genes, for spot-checking
genes <- dge_sig %>%
  pull(gene) %>%
  head(6)
genes

## Violin plots: expression of top significant genes, split by sample/condition
VlnPlot(seurat_myeloid, 
        genes,
        group.by = "orig.ident")
## Save to "figures/"

## UMAP of gene expression for the same top genes
FeaturePlot(seurat_myeloid, 
            genes, 
            reduction = "full.umap.sketch",
            ncol = 3)
## Save to "figures/"

## Spatial plot of expression for the top gene, on the tissue image
SpatialFeaturePlot(seurat_myeloid, 
                   genes[1], 
                   image.alpha = 0.1, 
                   pt.size.factor = 20)
## Save to "figures/"

## ---- Functional Analysis ----

## ---- Over-representation analysis (ORA) ----

## ---- Running ORA with clusterProfiler ----

## Background gene universe for hypergeometric testing = all genes tested for
## significance in the DGE results (not just the significant ones)
all_genes <- as.character(dge$gene)

## Extract significant, up-regulated genes only (tumor > normal)
sigUp <- dplyr::filter(dge_sig, 
                       p_val_adj < 0.05, 
                       avg_log2FC > 0)

## Convert the genes to a character vector
sigUp_genes <- as.character(sigUp$gene)

## Inspect the significantly up-regulated genes
sigUp_genes %>% head()

## Run GO Biological Process enrichment on the up-regulated gene set
egoUp <- enrichGO(gene = sigUp_genes, 
                  universe = all_genes,
                  keyType = "SYMBOL",
                  OrgDb = org.Hs.eg.db, 
                  ont = "BP", 
                  pAdjustMethod = "BH", 
                  qvalueCutoff = 0.05, 
                  readable = TRUE)

## ---- Exploring over-representation analysis results ----

## Output results from GO analysis to a data frame
cluster_summaryUp <- data.frame(egoUp)

## View GO up-regulated output
View(cluster_summaryUp)

## ---- Visualizing over-representation analysis results ----

## Dotplot of the top 20 enriched GO terms (size = gene count, color = padj)
dotplot(egoUp,
        showCategory = 20)
## Save to "figures/"

## ---- Gene Set Enrichment Analysis (GSEA) ----
## Unlike ORA above, GSEA uses the FULL ranked gene list (no hard
## significance cutoff), so it can detect coordinated pathway-level shifts.

## ---- Running GSEA with MSigDB gene sets ----

## Use the C5 (GO signatures) collection from MSigDB as the reference gene sets
m_t2g <- msigdbr(species = "Homo sapiens",
                 collection = "C5") %>%
  dplyr::select(gs_name, gene_symbol)

## Extract fold changes for all significant genes
foldchanges <- dge_sig$avg_log2FC

## Name each fold change with its corresponding gene symbol
names(foldchanges) <- dge_sig$gene

## GSEA requires the gene list sorted in decreasing order of fold change
foldchanges <- sort(foldchanges, decreasing = TRUE)

## Inspect the sorted, named fold-change vector
foldchanges %>% head()

## Run GSEA against the MSigDB C5 gene sets
msig_GSEA <- GSEA(foldchanges,
                  TERM2GENE = m_t2g,
                  verbose = FALSE)

## Extract the GSEA results table
msigGSEA_results <- msig_GSEA@result

## Write GSEA results to file
write.csv(msigGSEA_results,
          file = "results/gsea_msigdb_GO_genesets.csv",
          quote = FALSE)

## Look at results ordered by absolute Normalized Enrichment Score (NES)
msigGSEA_results %>%
  arrange(desc(abs(NES))) %>%
  View()

## ---- GSEA visualization ----

## Plot the GSEA running-enrichment-score plot for a single enriched GO term
gseaplot(msig_GSEA,
         geneSetID = 'GOBP_INFLAMMATORY_RESPONSE')
## Save to "figures/"

# ==================== 9. SPATIALLY VARIABLE GENES ===========================
# Identify genes whose expression shows significant spatial autocorrelation
# (i.e. neighboring bins have similar expression) within a tissue region,
# using Moran's I statistic.

library(tidyverse)
library(Seurat)
set.seed(12345)

## Increases the size of the default vector (8GB)
options(future.globals.maxSize = 8000 * 1024^2)
## Set parallelization (multithreading) to speed up calculations
library(future)
plan(multisession, workers = parallel::detectCores() - 1)

## Load dataset (cell-type-annotated checkpoint)
seurat_rctd <- readRDS("data/intermediate_seurat/seurat_RCTD.rds")

## Clear everything except the RDS object to save RAM
rm(list = setdiff(ls(), "seurat_rctd"))

## ---- Moran's I ----

## Subset to tumor bins in the P5CRC sample only
crc_sub <- subset(seurat_rctd,
                  subset = (orig.ident == "P5CRC") &
                    (first_type == "Tumor"))
crc_sub

## Re-calculate highly variable genes on this tumor-only subset
crc_sub <- FindVariableFeatures(crc_sub)

## Scale log-normalized counts
crc_sub <- ScaleData(crc_sub)
crc_sub

## Run Moran's I spatial autocorrelation test on the top 50 highly variable genes
crc_sub <- FindSpatiallyVariableFeatures(crc_sub,
                                         selection.method = "moransi",
                                         verbose = TRUE,
                                         features = VariableFeatures(crc_sub)[1:50],
                                         image = "P5CRC.008um")

## Inspect Moran's I scores for each tested gene, sorted by strongest
## spatial autocorrelation
SVFInfo(crc_sub, method = "moransi") %>%
  subset(!is.na(MoransI_observed)) %>%
  arrange(desc(MoransI_observed)) %>%
  head()

## Get the top 10 spatially variable genes (SVGs) by Moran's I
top10_svg <- SpatiallyVariableFeatures(crc_sub, 
                                       method = "moransi")
top10_svg <- top10_svg[1:10]
top10_svg

## Spatial plot of the top 6 SVGs, to visually confirm spatial patterning
SpatialFeaturePlot(crc_sub,
                   features = top10_svg[1:6],
                   pt.size.factor = 15,
                   image.alpha = 0.1, 
                   ncol = 3)
## Save to "figures/"

## Grab x,y coordinates of spatially constrained subtypes identified by BANKSY,
## along with sample ID and BANKSY cluster assignment for each bin
df <- FetchData(crc_sub,
                vars = c("x", "y", 
                         "banksy_cluster"))

## Subset to specific tumor-associated BANKSY clusters (13, 14, 15) and plot
## each cluster's spatial footprint in its own facet panel
ggplot(df %>% subset(banksy_cluster %in% c(13, 14, 15)),
       aes(x = x, y = -y,
           color = banksy_cluster)) +
  geom_point(size = 0.05) +
  theme_bw() +
  facet_wrap(~banksy_cluster) +
  NoLegend()
## Save to "figures/"

## Dotplot of tumor BANKSY clusters for the top SVGs - shows which cluster(s)
## most strongly express each spatially variable gene
DotPlot(crc_sub,
        top10_svg,
        group.by = "banksy_cluster",
        cluster.idents = TRUE,
        idents = c(13, 14, 15)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
## Save to "figures/"

# ============================= 10. CELLCHAT =================================
# CellChat infers cell-cell communication (ligand-receptor signaling)
# between annotated cell types, incorporating spatial proximity constraints.

library(tidyverse)
library(Seurat)
library(CellChat)
set.seed(12345)

## Increases the size of the default vector (8GB)
options(future.globals.maxSize = 8000 * 1024^2)
## Set parallelization (multithreading) to speed up calculations
library(future)
plan(multisession, workers = parallel::detectCores() - 1)

## Load dataset (cell-type-annotated checkpoint)
seurat_rctd <- readRDS("data/intermediate_seurat/seurat_RCTD.rds")

## Clear everything except the RDS object to save RAM
rm(list = setdiff(ls(), "seurat_rctd"))

## Subset to the CRC (tumor) sample and remove bins with no confident
## cell-type assignment
crc <- subset(seurat_rctd,
              subset = (orig.ident == "P5CRC") &
                (first_type != "unassigned"))

## CellChat expects a "samples" column that is a factor
crc$samples <- factor(crc$orig.ident)
crc

## ---- Create CellChat object ----

## Get normalized (log-transformed) counts matrix, required input for CellChat
counts_norm_crc <- GetAssayData(crc,
                                assay = "Spatial.008um",
                                layer = "data")
counts_norm_crc[1:5, 1:5]

## Get spatial (x,y) coordinates for each bin
coords_crc <- FetchData(crc, vars = c("x", "y"))
coords_crc %>% head()

## Scale factors stored in the Space Ranger JSON output - needed to convert
## between pixel and micron distance units for spatial signaling range constraints
sf_json <- jsonlite::fromJSON("data/P5CRC_cropped/binned_outputs/square_008um/spatial/scalefactors_json.json")

## Create spatial.factors dataframe: pixel-to-micron ratio and neighbor
## tolerance radius, both derived from the Space Ranger scale factors
spatial_factors_crc <- data.frame(
  ratio = sf_json$bin_size_um / sf_json$spot_diameter_fullres,
  tol = sf_json$bin_size_um / 2)
spatial_factors_crc

## Create the CellChat object, grouping bins by their RCTD-assigned cell type
cellchat_crc <- createCellChat(
  object = counts_norm_crc,
  coordinates = coords_crc,
  spatial.factors = spatial_factors_crc,
  meta = crc@meta.data,
  group.by = "first_type",
  datatype = "spatial")

cellchat_crc

## Load the human CellChat ligand-receptor interaction database
cellchat_db <- CellChatDB.human

## Assign the database to the CellChat object's DB slot
cellchat_crc@DB <- cellchat_db

## Subset step is required even if the full database is used (initializes
## internal data structures for downstream steps)
cellchat_crc <- subsetData(cellchat_crc)

## Identify genes that are over-expressed in specific cell types relative
## to all others (used to prioritize likely signaling genes)
cellchat_crc <- identifyOverExpressedGenes(cellchat_crc)

## Print the top 2 over-expressed genes per cell type
cellchat_crc@var.features$features.info %>%
  group_by(clusters) %>%
  arrange(pvalues.adj, desc(logFC)) %>%
  slice_head(n = 2) %>% 
  View()

## ---- CellChat Checkpoint (! time-consuming) (Start) ----

## ---- Run CellChat and filter results ----
## DO NOT RUN (unless you have time/compute available)

## Identify ligand-receptor interactions that are over-expressed relative
## to other cell type pairs
cellchat_crc <- identifyOverExpressedInteractions(cellchat_crc)

## DO NOT RUN
## (Optional) Smooth expression over the protein-protein interaction network
## to reduce dropout noise. If run, set raw.use = FALSE in computeCommunProb()
cellchat_crc <- smoothData(cellchat_crc,
                           adj = PPI.human)
# For older versions of CellChat, the function is:
# cellchat <- projectData(cellchat, PPI.human)

## DO NOT RUN
## Calculate communication probability for every ligand-receptor pair,
## incorporating spatial distance constraints (interaction.range/contact.range)
cellchat_crc <- computeCommunProb(
  cellchat_crc,
  type = "trimean",
  distance.use = TRUE,
  interaction.range = 200,
  scale.distance = 0.2,
  contact.dependent = TRUE,
  contact.range = 10,
  raw.use = FALSE,
  seed.use = 12345)

## ---- CellChat Checkpoint (End) ----

## ---- Load pre-computed RDS object ----
cellchat_crc <- readRDS("data/intermediate_seurat/cellchat_crc_computeCommunProb.rds")

## Remove cell-type pairs supported by too few bins (unreliable estimates)
cellchat_crc <- filterCommunication(cellchat_crc,
                                    min.cells = 5)

## Summarize ligand-receptor pair probabilities into signaling pathway-level
## probabilities (groups related ligand-receptor pairs into named pathways)
cellchat_crc <- computeCommunProbPathway(cellchat_crc)

## Count total interactions and sum interaction weights into a flat network matrix
cellchat_crc <- aggregateNet(cellchat_crc)

## Compute how much each cell type sends vs. receives signal, across all pathways
cellchat_crc <- netAnalysis_computeCentrality(cellchat_crc, slot.name = "netP")

## ---- Visualize CellChat results ----

## Arrange 2 plots side-by-side
par(mfrow = c(1, 2), xpd = TRUE) 

## Circle plot: number of ligand-receptor interactions between cell types
netVisual_circle(
  cellchat_crc@net$count,
  weight.scale = TRUE,
  label.edge = FALSE,
  edge.width.max = 12,
  title.name = "Number of interactions")

## Circle plot: total strength (probability-weighted) of interactions
netVisual_circle(
  cellchat_crc@net$weight,
  weight.scale = TRUE,
  label.edge = FALSE,
  edge.width.max = 12,
  title.name = "Interaction strength")
## Save to "figures/"

## Reset par to avoid affecting later single-panel plots
par(mfrow = c(1, 1))

## Heatmap: number of interactions between cell type pairs
gg_h_count  <- netVisual_heatmap(
  cellchat_crc, 
  measure = "count",
  color.heatmap = "Blues",
  title.name = "Number of interactions")

## Heatmap: strength of interactions between cell type pairs
gg_h_weight <- netVisual_heatmap(
  cellchat_crc, 
  measure = "weight",
  color.heatmap = "Reds",
  title.name = "Interaction strength")

## Plot heatmaps side by side
gg_h_count + gg_h_weight
## Save to "figures/"

## Scatterplot: incoming vs. outgoing signaling strength for each cell type
## (identifies "sender"/"receiver"/"hub" cell types)
netAnalysis_signalingRole_scatter(cellchat_crc)
## Save to "figures/"

## Ranked barplot: which pathways contribute most to overall signaling strength
rankNet(cellchat_crc, 
        measure = "weight",
        mode = "single", 
        stacked = TRUE, 
        do.stat = TRUE) +
  xlab("Pathways")
## Save to "figures/"

## Heatmap of outgoing pathway signaling strength, by cell type
ht_out <- netAnalysis_signalingRole_heatmap(
  cellchat_crc, 
  pattern = "outgoing",
  width = 10,
  height = 20,
  color.heatmap = "YlOrRd",
  title = "Outgoing signaling strength")

## Heatmap of incoming pathway signaling strength, by cell type
ht_in  <- netAnalysis_signalingRole_heatmap(
  cellchat_crc, 
  pattern = "incoming",
  width = 10,
  height = 15,
  color.heatmap = "GnBu",
  title = "Incoming signaling strength")

## Plot heatmaps side by side
ht_out + ht_in
## Save to "figures/"

## Chord diagram: visualize the full CEACAM signaling pathway network
netVisual_aggregate(cellchat_crc,
                    signaling  = "CEACAM",
                    layout = "chord",
                    title.name = "CEACAM signaling")
## Save to "figures/"

## ---- Bubble plot ----

## Restrict to the top 5 most significant signaling pathways
top_pathways <- cellchat_crc@netP$pathways[1:5]

## Count number of unique cell types present, for indexing sources
n_ct <- seurat_rctd@meta.data %>%
  pull(first_type) %>%
  unique() %>% length()

## Bubble plot: ligand-receptor pair probabilities from all cell types
## toward the Tumor cell type specifically
netVisual_bubble(cellchat_crc,
                 signaling = top_pathways,
                 sort.by.source.priority = TRUE,
                 sources.use = 1:n_ct,
                 targets.use = "Tumor",
                 remove.isolate = TRUE,
                 color.heatmap = "viridis",
                 angle.x = 90)

## ---- Zooming in on specific ligand-receptor complexes ----
ligand <- "COL3A1"
receptor <- "SDC4"

## Spatial plot of ligand and receptor expression, on the tissue image
SpatialFeaturePlot(crc, 
                   features = c(ligand, receptor),
                   pt.size.factor = 15,
                   image.alpha = 0.2)
## Save to "figures/"

## Violin plot: ligand/receptor expression across cell types
VlnPlot(crc, 
        features = c(ligand, receptor), 
        group.by = "first_type",
        pt.size = 0,
        ncol = 1)
top_pathways
## Save to "figures/"

# ============================= END OF SCRIPT ================================
# Record session info and outputs for reproducibility
sink("sessionInfo_spatial_transcriptomics.txt")
sessionInfo()
sink()
