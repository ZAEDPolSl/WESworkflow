#!/usr/bin/env Rscript

# This script performs clustering of WES-derived gene-level features using the PARC algorithm:
#
# (1) Loads previously computed UMAP embeddings and gene-level feature matrix
#     generated in earlier pipeline steps
#
# (2) Applies PARC clustering (via Python implementation using reticulate)
#     on the feature space using user-defined parameters (k-nearest neighbors, resolution)
#
# (3) Maps cluster assignments back to UMAP coordinates for visualization
#
# (4) Generates and exports UMAP plots colored by cluster labels
#
# (5) Saves clustering results and sample-to-cluster mappings for downstream analyses
#
# Input:
# - UMAP coordinates file (raw_umap_result.tsv)
# - Long-format feature table (raw_ft_long.tsv)
#
# Output:
# - Clustered UMAP coordinates (umap_parc_k{K}r{R}.tsv)
# - Cluster visualization (PARC_clustering.pdf)
# - Sample-to-cluster mapping (sample_kit_cluster_map.tsv)
#
# Note:
# - PARC clustering is executed via a Python backend (run_parc.py)
# - Input data (features and metadata) are not included and must be provided by the user
# - Clustering parameters (k, resolution) should be adjusted depending on dataset size and structure

library(data.table)
library(ggplot2)
library(httpgd)
library(tidyr)
library(dplyr)
library(stringr)
library(reshape2)
library(ggrepel)
library(RColorBrewer)
library(patchwork)
library(colorspace)


dir <- "/path/to/results/dir" # directory containing sudirectories from Variant_post_processing and Variant_to_gene steps
dir.create(file.path(dir, "Results"), showWarnings = FALSE)

# Load UMAP coordinates and long-format feature table from previous step
umap_result <- fread(file.path(dir, "Results", "raw_umap_result.tsv"))
long_tsv <- file.path(dir, "Results", "raw_ft_long.tsv")

# Load PARC clustering implementation (Python) via reticulate
library(reticulate)
source_python("functions/run_parc.py")

# PARC clustering parameters
k = 30    # number of nearest neighbors
res = 1   # clustering resolution



run_parc <- function(long_tsv, knn, res) {
  k <- knn
  res <- res
  samples_labels <- parc_labels_from_long_tsv(long_tsv, knn = as.integer(k), res = res)
  samples <- unlist(samples_labels[[1]])
  parc   <- as.integer(unlist(samples_labels[[2]]))
  lab <- data.frame(Sample = samples, PARC = parc)

  umap_result_parc <- umap_result
  umap_result_parc$Cluster <- as.factor(lab$PARC[match(umap_result_parc$Sample, lab$Sample)])
  umap_params <- list(
    n_neighbors = 15,
    min_dist = 0.5,
    metric = "cosine",
    n_epochs = 1500, 
    nn_method = "nndescent"
    )


  plot_title = paste0("k=", k, " resolution=", res)


  p <- ggplot(umap_result_parc, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
        geom_point(size = 1)+
        theme_test() +
        labs(
          title = plot_title,
          x = "UMAP 1",
          y = "UMAP 2"
        ) + scale_color_discrete_qualitative(palette = "Dark3")

  print(p)
  return(list("plot" = p, "data" = umap_result_parc))
}


parc <- run_parc(long_tsv, knn = k, res = res)
print(parc$plot)

# Save PARC UMAP
fwrite(parc$data, file.path(dir, "Results", paste0("umap_parc_k", ks, "r", reso,".tsv")), sep = "\t")

# Save sample-to-cluster mapping
map <- parc$data %>% select(-UMAP1, -UMAP2)
fwrite(map, file.path(dir, "Results", "sample_kit_cluster_map.tsv"), sep = "\t")

# Save UMAP colored by clusters
pdf(file.path(dir, "Results", "Figures", "PARC_clustering.pdf"), width = 6, height = 5.5)
print(parc$plot)
dev.off()





