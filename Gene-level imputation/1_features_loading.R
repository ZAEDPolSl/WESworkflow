#!/usr/bin/env Rscript

# This script performs downstream analysis of gene-level features derived from WES data:
#
# (1) Loads per-sample, per-gene features computed in previous pipeline steps
#     (Variant_post_processing and Variant_to_gene), and aggregates them into a single dataset
#
# (2) Merges feature data with sample metadata (Dataset, capture kit information)
#     and saves a long-format feature table
#
# (3) Transforms the data into a sample × gene matrix and computes a UMAP embedding
#     to visualize sample relationships based on CADD-weighted allele frequency features
#
# (4) Computes batch effect metrics (LISI, kBET) using external helper functions
#
# (5) Generates UMAP visualizations colored by dataset and capture kit,
#     and exports both tabular results and figures
#
# Input:
# - Directory with per-chromosome feature files (*.feature.txt)
# - Metadata file with Sample, Dataset, and Kit columns
#
# Output:
# - Long-format feature table (raw_ft_long.tsv)
# - UMAP coordinates (raw_umap_result.tsv) with computed batch metrics
# - UMAP plots (PDF)

library(data.table)
library(ggplot2)
library(httpgd)
library(rnndescent)
library(uwot)
library(tidyr)
library(dplyr)

dir <- "/path/to/results/dir" # directory containing sudirectories from Variant_post_processing and Variant_to_gene steps
dir.create(file.path(dir, "Results"), showWarnings = FALSE)


meta <- fread("Data/example/sample_path_map_example.tsv") %>%
	select(Sample, Dataset, Kit)

# Step 1: Load the data ==========================================================

load_chr_features <- function(feature_files_dir){
  files <- list.files(
    feature_files_dir,
    pattern = "*.feature.txt", recursive = TRUE, full.names = TRUE
  )

  all_features <- rbindlist(lapply(files, function(f) {
    dt <- fread(f, select = c("Gene", "CADD_weighted_avg_AF"))
    sample_name <- sub(".feature.txt$", "", basename(f))

    dt[, Sample := sample_name]
    return(dt)
  }), use.names = TRUE, fill = TRUE)


  return(all_features)
}


features <- load_chr_features(file.path(dir, "Features"))


features <- left_join(features, meta, by = "Sample")

fwrite(features,
       file = file.path(dir, "Results", "raw_ft_long.tsv"),
       sep = "\t")

# Step 2: Run UMAP ====================================================
umap_input <- features[, .(Sample, Gene, CADD_weighted_avg_AF)]

cat("Reshaping data to wide format...\n")
umap_input <- dcast(
  umap_input,
  Sample ~ Gene,
  value.var = "CADD_weighted_avg_AF",
  fill = 0
)
sample <- umap_input$Sample
umap_input <- as.matrix(umap_input[, -1])

rownames(umap_input) <- sample

umap_params <- list(
  n_neighbors = 15,
  min_dist = 0.5,
  metric = "cosine",
  n_epochs = 1500, 
  nn_method = "nndescent"
  )

set.seed(123)
cat("Running UMAP...\n")
umap_result <- uwot::umap(
  umap_input,
  n_neighbors = umap_params$n_neighbors,
  min_dist = umap_params$min_dist,
  metric = umap_params$metric,
  n_epochs = umap_params$n_epochs,
  n_threads = parallel::detectCores(),
  ret_model = FALSE,  # do not return the model,
  nn_method = umap_params$nn_method 
)


umap_result <- as.data.frame(umap_result)
colnames(umap_result) <- c("UMAP1", "UMAP2")

umap_result$Sample <- sample
umap_result <- left_join(umap_result, meta, by = "Sample")
umap_result$Kit_short <- kit_map[umap_result$Kit]



fwrite(umap_result, 
      file.path(dir, "Results", "raw_umap_result.tsv"),
      sep = "\t")




# Step 3: Plot the results ================================================
source("functions/batch_metrics.R")
batch_stats <- compute_batch_metrics_df(features)
batch_stats$lisi_stats
batch_stats$kbet_stats


umap_caption <- paste0("n_neighbors=", umap_params$n_neighbors, 
  ", min_dist=", umap_params$min_dist, 
  ", metric=", umap_params$metric, 
  ", n_epochs=", umap_params$n_epochs
)

plot_title = "Genotype imputation only\n1kG EUR panel"

umap_plot <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = Dataset)) +
                    geom_point(size = 1)+
                    theme_test() +
                    labs(
                      title = plot_title,
                      x = "UMAP 1",
                      y = "UMAP 2",
                      subtitle = paste0("LISI = ", round(batch_stats$lisi_stats$mean,2), 
                                        " | kBET = ", round(batch_stats$kbet_stats$mean,3)),
                      caption = umap_caption) 
 

print(umap_plot)


umap_plot_kit <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = Kit)) +
                    geom_point(size = 1)+
                    theme_test() +
                    labs(
                      title = plot_title,
                      x = "UMAP 1",
                      y = "UMAP 2",
                      color = "Capture kit",
                      subtitle = paste0("LISI = ", round(batch_stats$lisi_stats$mean,2), 
                                        " | kBET = ", round(batch_stats$kbet_stats$mean,3)),
                      caption = umap_caption) 
 

print(umap_plot_kit)

pdf(file.path(dir, "Results", "Figures", "raw_UMAP.pdf"), width = 7, height = 6)
print(umap_plot)
print(umap_plot_kit)
dev.off()

