#!/usr/bin/env Rscript

# This script performs MNAR-aware feature imputation and downstream analysis
# on gene-level WES-derived features:
#
# (1) Loads gene-level feature data and cluster assignments generated in previous steps
#
# (2) Identifies Missing Not At Random (MNAR) patterns at the gene × cluster level
#     based on detection rate thresholds derived from a Gaussian Mixture Model (GMM)
#
# (3) Constructs a sample × gene feature matrix and applies masked KNN imputation
#     (cosine similarity, cross-group donor selection) to impute only MNAR values
#
# (4) Saves imputed feature matrices in both wide and long formats
#
# (5) Evaluates batch effects using LISI and kBET metrics
#
# (6) Recomputes UMAP embeddings on imputed features and generates visualizations
#     to assess the impact of imputation on data integration
#
# Input:
# - Long-format feature table (raw_ft_long.tsv)
# - Sample-to-cluster mapping (sample_kit_cluster_map.tsv)
# - GMM model with detection rate thresholds (GMM_det_rate.RDS)
# - Metadata file (sample_path_map_example.tsv)
#
# Output:
# - Imputed feature matrices (wide and long format)
# - UMAP coordinates after imputation
# - UMAP plots (dataset- and cluster-colored)
# - Batch effect metrics (printed to console)
#
# Notes:
# - Imputation is restricted to MNAR entries only (masked imputation)
# - External resources (metadata, feature files) are not included in the repository
# - Parameters (a, b, k) control MNAR detection and KNN imputation behavior
# - Computationally intensive (parallel KNN, large matrices)

library(data.table)
library(ggplot2)
library(httpgd)
library(tidyr)
library(dplyr)
library(reshape2)
library(patchwork)
library(rnndescent)
library(uwot)
library(ComplexHeatmap)
library(circlize)

dir <- "/path/to/results/dir" # directory containing sudirectories from Variant_post_processing and Variant_to_gene steps

# Load cluster assignments and long-format feature table from previous step
ft <- fread( file.path(dir, "Results", "raw_ft_long.tsv") )
clusters <- fread( file.path(dir, "Results", "sample_kit_cluster_map.tsv") )

ft$Group <- as.factor(clusters$Cluster[match(ft$Sample, clusters$Sample)])
setDT(ft)

gmm <- readRDS( file.path(dir, "Results", "GMM_det_rate.RDS"))
thresholds <- gmm$threshold

a_params <- thresholds[3] # low detection in group
b_params <- thresholds[6] # high detection in some other group
k = 10
min_n <- 15


fname_base <- function(a, b, k, prefix = "ft_imp", digits = 2) {
      fmt <- function(x) sprintf(paste0("%.", digits, "f"), x)   # keeps "0.10"
      paste0(prefix, "_a", fmt(a), "_b", fmt(b), "_k", k)
    }

options(future.globals.maxSize = 4 * 1024^3)  # 4 GB

for (a in a_params){
  for (b in b_params){
    # ft: data.table with columns: Sample, Group, Gene, CADD_weighted_avg_AF
    base <- fname_base(a, b, k)
    # examples
    #saveRDS(X_imp,   paste0(dir, base, "_raw.rds"))
    #saveRDS(df_long, paste0(dir, base, "_long.rds"))

  # STEP 1: flag MNAR ==================================================
    cat("Flagging MNAR...\n")
    # 1) Total samples per group (denominator for detection rate)
    n_total <- unique(ft[, .(Sample, Group)])[, .(n_total = .N), by = Group]
    # 2) Detections per Gene x group (presence = any value > 0 in sample)
    det_tbl <- ft[
      CADD_weighted_avg_AF > 0, 
      .(det_n = uniqueN(Sample)),
                  by = .(Gene, Group)]
    # 3) Complete all Gene x group pairs and fill missing with zeros
    all_pairs <- CJ(Gene = unique(ft$Gene), Group = unique(ft$Group))
    det_tbl <- det_tbl[all_pairs, on = .(Gene, Group)]
    det_tbl[is.na(det_n), det_n := 0]
    # 4) Attach n_total and compute detection rate
    det_tbl <- det_tbl[
      n_total, on = "Group"][, det_rate := det_n / n_total][]

    # 5) MNAR flagging with min_n applied ONLY to OTHER groups
        # For each (Gene, Group): compute max detection rate among OTHER groups with n_total >= min_n
    others <- det_tbl[
      n_total >= min_n, .(Gene, Group_other = Group, det_rate_other = det_rate)]
    max_other_tbl <- others[
      det_tbl[, .(Gene, Group)], on = .(Gene), allow.cartesian = TRUE][Group_other != Group][, .(max_other = max(det_rate_other)), by = .(Gene, Group)]

    # 6) Flags: current group's det_rate < a AND any other group's max_other > b
    flags <- det_tbl[
      max_other_tbl, on = .(Gene, Group)][
      det_rate < a & max_other > b,
      .(Gene, Group, det_rate, n_total, max_other)]

    # flags: data.table of MNAR candidates (Gene×group)

  # STEP 2: Prepare data for imputation ================================
    cat("Preparing the data for feature imputation...\n")
    # 1) Wide with NA (no zeros)
    mat <- dcast(
      ft, Sample + Group ~ Gene,
                value.var = "CADD_weighted_avg_AF",
                fun.aggregate = function(x) if (length(x)) mean(x) else NA_real_,
                fill = NA_real_)
    X <- as.matrix(mat[, -(1:2)])
    rownames(X) <- mat$Sample
    clustering <- mat$Group
    Genes <- colnames(X)

    # 2) MNAR mask from flags: impute only NA where (Gene, group) is flagged, rest is true 0
    flag_list <- split(flags$Group, flags$Gene)
    MNAR <- matrix(FALSE, nrow(X), ncol(X), dimnames = list(rownames(X), Genes))
    for (g in names(flag_list)) if (g %in% Genes) {
      MNAR[clustering %in% flag_list[[g]], g] <- is.na(X[clustering %in% flag_list[[g]], g])
    }

  # STEP 3: Feature imputation =========================================
   
    cat(paste0("Performing feature imputation with parameters:\ta=", a, "\tb=", b, "...\n"))
    source("functions/knn_impute_cosine_parallel.R")
    set.seed(123)
    X_imp <- knn_impute_mnar_masked_parallel(X_raw = X, MNAR = MNAR, 
                                              grp = clustering, 
                                              workers = 48,
                                              k = k, 
                                              min_k = 5)

    X_imp[X_imp < 0] <- 0
    X_imp[X_imp > 1] <- 1



  # STEP 4: Save the results ===========================================
    meta <- fread("Data/example/sample_path_map_example.tsv") %>%
      select(Sample, Dataset, Kit)

    X_imp <- X_imp %>%
          as.data.frame() %>%
          tibble::rownames_to_column("Sample") %>%
          left_join(meta, by = "Sample") %>%
          mutate(Group = clusters$Cluster[match(Sample, clusters$Sample)]) %>%
          select(Sample, Dataset, Group, Kit, dplyr::everything()) %>% setDT()


    df_long <- data.table::melt(X_imp, id.vars = c("Sample", "Dataset", "Group", "Kit"),
                    variable.name = "Gene", value.name = "CADD_weighted_avg_AF")

    data.table::fwrite(X_imp,  file.path(dir, paste0(base, "_wide.tsv")), sep = "\t")
    data.table::fwrite(df_long, file.path(dir, paste0(base, "_long.tsv")), sep = "\t")

  # STEP 5: Postprocessing =============================================
    # batch metrics 
    cat("Computing batch metrics...\n")
    source("functions/batch_metrics.R")
    batch_stats <- compute_batch_metrics_df(df_long %>% select(Sample, Dataset, Gene, CADD_weighted_avg_AF))
    cat(paste("LISI:", batch_stats$lisi_stats$mean,"\n"))
    cat(paste("kBET:", batch_stats$kbet_stats$mean,"\n"))

    # UMAP 
    cat("Running UMAP...\n")

    umap_input <- X_imp %>% select(-Sample, -Dataset, -Group, -Kit) %>% as.matrix()
    metadata <- X_imp %>% select(Sample, Dataset, Group, Kit)
    umap_params <- list(
      n_neighbors = 15,
      min_dist = 0.5,
      metric = "cosine",
      n_epochs = 1500, 
      nn_method = "nndescent"
      )

    set.seed(123)  
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
    umap_result <- cbind(umap_result, metadata)
    umap_result$Group <- factor(umap_result$Group, levels = sort(unique(umap_result$Group)))
    data.table::fwrite(umap_result, file.path(dir, paste0(base, "_umap_result.tsv")), sep = "\t")


  # Plot the UMAP ------------------------------------------------------
    cat("Making the figures...\n")
    umap_caption <- paste0("UMAP params: n_neighbors=", umap_params$n_neighbors, 
      ", min_dist=", umap_params$min_dist, 
      ", metric=", umap_params$metric, 
      ", n_epochs=", umap_params$n_epochs
    )

    plot_title = paste0("Feature Imputation\na=", round(a, 3), " b=", round(b, 3), " kn=", k,
                        "\nLISI=", round(batch_stats$lisi_stats$mean,2), 
                        " | kBET=", round(batch_stats$kbet_stats$mean,3))

    umap_plot <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = Dataset)) +
                        geom_point()+
                        theme_test() +
                        labs(
                          title = plot_title,
                          x = "UMAP 1",
                          y = "UMAP 2",
                          caption = umap_caption
                        ) + scale_color_hue()
    print(umap_plot)
    pdf(file.path(dir, "Results","Figures", paste0(gsub("ft_imp_", "", base), "_UMAP.pdf")), width = 7, height = 6)
    print(umap_plot)
    dev.off()


    umap_plot2 <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = as.factor(Group))) +
                        geom_point()+
                        theme_test() +
                        labs(
                          title = plot_title,
                          x = "UMAP 1",
                          y = "UMAP 2",
                          caption = umap_caption,
                          color = "Cluster"
                        ) + scale_color_hue() 
    print(umap_plot2)
    pdf(file.path(dir, "Results","Figures", paste0(gsub("ft_imp_", "", base), "_UMAP_cluster_colored.pdf")), width = 6, height = 6)
    print(umap_plot2)
    dev.off()
    cat("UMAP: DONE\n")
  }
}