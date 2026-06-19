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
#
# Usage:
# R_LIBS= R_LIBS_USER= R_LIBS_SITE= R_PROFILE_USER=/dev/null \
# Rscript --vanilla "Gene-level imputation/4_feature_imputation.R" config/local_config.yaml

conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library")
if (nzchar(Sys.getenv("CONDA_PREFIX")) && dir.exists(conda_lib)) {
	.libPaths(conda_lib)
}

# ======== ! SELECT IMPUTATION THERSHOLDS BASED ON GMM RESULTS ! ===========
# Detection-rate thresholds used for MNAR flagging.
# A gene is considered under-detected in a group when its detection rate
# in that group is <= threshold_low_value, while its detection rate in
# another sufficiently large group is >= threshold_high_value.
# currently set thresholds are adjusted to example run
threshold_low_value <- 0.44
threshold_high_value <- 0.85

# ================== these parameters can be adjusted ================

knn_k <- 10 # objective number of nearest neighbors
knn_min_k <- 5 # minimal number of nearest neighbors to impute missing value
min_other_group_size <- 15 # minimal size of group with high detection rate to claim MNAR

# ======== these should be consistent across previous steps ==========
# UMAP parameters 
plot_title <- "Data after genotype imputation"
umap_params <- list(
	n_neighbors = 15,
	min_dist = 0.5,
	metric = "cosine",
	n_epochs = 1500,
	nn_method = "nndescent",
	seed = 123
)
# Batch metric parameters
batch_metric_params <- list(
	lisi_perplexity = 30,
	k_kBET = 15,
	test_size_fraction_kBET = 0.1,
	heuristic_kBET = FALSE,
	adapt_kBET = FALSE,
	PCA_kBET = TRUE
)
# ===========================================================


suppressPackageStartupMessages({
	library(data.table)
	library(ggplot2)
	library(tidyr)
	library(dplyr)
	library(reshape2)
	library(patchwork)
	library(rnndescent)
	library(uwot)
})


get_script_path <- function() {
	args <- commandArgs(trailingOnly = FALSE)
	file_arg <- grep("^--file=", args, value = TRUE)

	if (length(file_arg) == 0) {
		stop("Cannot determine script path.")
	}

	script_path <- sub("^--file=", "", file_arg[1])
	script_path <- gsub("~\\+~", " ", script_path)

	normalizePath(script_path, mustWork = TRUE)
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
repo_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
config <- if (length(args) >= 1) args[1] else "config/local_config.yaml"

if (!grepl("^/", config)) {
	config <- file.path(repo_dir, config)
}

config <- normalizePath(config, mustWork = TRUE)

read_config <- function(key) {
	value <- system2(
		"python",
		args = c(file.path(repo_dir, "scripts/read_config.py"), config, key),
		stdout = TRUE
	)

	if (!is.null(attr(value, "status"))) {
		stop(paste("Failed to read config key:", key))
	}

	if (length(value) == 0 || is.na(value[1]) || value[1] == "") {
		stop(paste("Empty config value for key:", key))
	}

	value[1]
}

resolve_path <- function(path) {
	if (length(path) == 0 || is.na(path) || path == "") {
		stop("Cannot resolve an empty path.")
	}

	if (grepl("^/", path)) {
		return(path)
	}

	file.path(repo_dir, path)
}

dir <- resolve_path(read_config("directories.results_dir"))
results_dir <- file.path(dir, "Gene_level_imputation")
figures_dir <- file.path(results_dir, "Figures")
metadata_file <- resolve_path(read_config("directories.sample_metadata"))
future_workers <- as.integer(read_config("parameters.feature_imputation_workers"))
future_max_size_gb <- as.numeric(read_config("parameters.feature_imputation_future_max_size_gb"))

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

ft_file <- file.path(results_dir, "raw_ft_long.tsv")
clusters_file <- file.path(results_dir, "sample_kit_cluster_map.tsv")


knn_function_file <- file.path(script_dir, "functions", "knn_impute_cosine_parallel.R")
batch_metrics_file <- file.path(script_dir, "functions", "batch_metrics.R")

if (is.na(future_workers) || future_workers < 1) {
	stop("parameters.feature_imputation_workers must be a positive integer.")
}

if (is.na(future_max_size_gb) || future_max_size_gb <= 0) {
	stop("parameters.feature_imputation_future_max_size_gb must be a positive number.")
}

if (!file.exists(ft_file)) {
	stop("Long-format feature table not found: ", ft_file)
}

if (!file.exists(clusters_file)) {
	stop("Sample-to-cluster mapping not found: ", clusters_file)
}


if (!file.exists(metadata_file)) {
	stop("Metadata file not found: ", metadata_file)
}

if (!file.exists(knn_function_file)) {
	stop("KNN imputation function file not found: ", knn_function_file)
}

if (!file.exists(batch_metrics_file)) {
	stop("Batch metrics function file not found: ", batch_metrics_file)
}

# Load cluster assignments and long-format feature table from previous step
ft <- fread(ft_file)
clusters <- fread(clusters_file)

required_ft_cols <- c("Sample", "Gene", "CADD_weighted_avg_AF")
missing_ft_cols <- setdiff(required_ft_cols, colnames(ft))

if (length(missing_ft_cols) > 0) {
	stop("Missing feature columns: ", paste(missing_ft_cols, collapse = ", "))
}

required_cluster_cols <- c("Sample", "Cluster")
missing_cluster_cols <- setdiff(required_cluster_cols, colnames(clusters))

if (length(missing_cluster_cols) > 0) {
	stop("Missing cluster columns: ", paste(missing_cluster_cols, collapse = ", "))
}

ft$Group <- as.factor(clusters$Cluster[match(ft$Sample, clusters$Sample)])
setDT(ft)

ft <- ft[!is.na(Group)]

if (nrow(ft) == 0) {
	stop("No feature records matched cluster assignments.")
}

# Use manually defined detection-rate thresholds.
a_params <- threshold_low_value # low detection in group
b_params <- threshold_high_value # high detection in some other sufficiently large group

k <- knn_k
min_n <- min_other_group_size

cat("Using detection-rate thresholds:\n")
cat("  low detection threshold a =", a_params, "\n")
cat("  high detection threshold b =", b_params, "\n")

fname_base <- function(a, b, k, prefix = "ft_imp", digits = 2) {
	fmt <- function(x) sprintf(paste0("%.", digits, "f"), x) # keeps "0.10"
	paste0(prefix, "_a", fmt(a), "_b", fmt(b), "_k", k)
}

options(future.globals.maxSize = future_max_size_gb * 1024^3)

source(knn_function_file)
source(batch_metrics_file)

for (a in a_params) {
	for (b in b_params) {
		# ft: data.table with columns: Sample, Group, Gene, CADD_weighted_avg_AF
		base <- fname_base(a, b, k)
		# examples
		# saveRDS(X_imp,   paste0(dir, base, "_raw.rds"))
		# saveRDS(df_long, paste0(dir, base, "_long.rds"))

		# STEP 1: flag MNAR ==================================================
		cat("Flagging MNAR...\n")

		# 1) Total samples per group (denominator for detection rate)
		n_total <- unique(ft[, .(Sample, Group)])[, .(n_total = .N), by = Group]

		# 2) Detections per Gene x group (presence = any value > 0 in sample)
		det_tbl <- ft[
			CADD_weighted_avg_AF > 0,
			.(det_n = uniqueN(Sample)),
			by = .(Gene, Group)
		]

		# 3) Complete all Gene x group pairs and fill missing with zeros
		all_pairs <- CJ(Gene = unique(ft$Gene), Group = unique(ft$Group))
		det_tbl <- det_tbl[all_pairs, on = .(Gene, Group)]
		det_tbl[is.na(det_n), det_n := 0]

		# 4) Attach n_total and compute detection rate
		det_tbl <- det_tbl[
			n_total,
			on = "Group"
		][, det_rate := det_n / n_total][]

		# 5) MNAR flagging with min_n applied ONLY to OTHER groups
		# For each (Gene, Group): compute max detection rate among OTHER groups with n_total >= min_n
		others <- det_tbl[
			n_total >= min_n,
			.(Gene, Group_other = Group, det_rate_other = det_rate)
		]

		max_other_tbl <- others[
			det_tbl[, .(Gene, Group)],
			on = .(Gene),
			allow.cartesian = TRUE
		][
			Group_other != Group
		][
			,
			.(max_other = max(det_rate_other)),
			by = .(Gene, Group)
		]

		# 6) Flags: current group's det_rate < a AND any other group's max_other > b
		flags <- det_tbl[
			max_other_tbl,
			on = .(Gene, Group)
		][
			det_rate < a & max_other > b,
			.(Gene, Group, det_rate, n_total, max_other)
		]

		# flags: data.table of MNAR candidates (Gene×group)

		# STEP 2: Prepare data for imputation ================================
		cat("Preparing the data for feature imputation...\n")

		# 1) Wide with NA (no zeros)
		mat <- dcast(
			ft,
			Sample + Group ~ Gene,
			value.var = "CADD_weighted_avg_AF",
			fun.aggregate = function(x) if (length(x)) mean(x) else NA_real_,
			fill = NA_real_
		)

		X <- as.matrix(mat[, -(1:2)])
		rownames(X) <- mat$Sample
		clustering <- mat$Group
		Genes <- colnames(X)

		# 2) MNAR mask from flags: impute only NA where (Gene, group) is flagged, rest is true 0
		flag_list <- split(flags$Group, flags$Gene)
		MNAR <- matrix(FALSE, nrow(X), ncol(X), dimnames = list(rownames(X), Genes))

		for (g in names(flag_list)) {
			if (g %in% Genes) {
				rows <- clustering %in% flag_list[[g]]
				MNAR[rows, g] <- is.na(X[rows, g])
			}
		}

		# STEP 3: Feature imputation =========================================
		cat(
			paste0(
				"Performing feature imputation with parameters:\ta=",
				a,
				"\tb=",
				b,
				"\tk=",
				k,
				"...\n"
			)
		)

		set.seed(umap_params$seed)

		X_imp <- knn_impute_mnar_masked_parallel(
			X_raw = X,
			MNAR = MNAR,
			grp = clustering,
			workers = future_workers,
			k = k,
			min_k = knn_min_k
		)

		X_imp[X_imp < 0] <- 0
		X_imp[X_imp > 1] <- 1

		# STEP 4: Save the results ===========================================
		meta <- fread(metadata_file) %>%
			select(Sample, Dataset, Kit)

		X_imp <- X_imp %>%
			as.data.frame() %>%
			tibble::rownames_to_column("Sample") %>%
			left_join(meta, by = "Sample") %>%
			mutate(Group = clusters$Cluster[match(Sample, clusters$Sample)]) %>%
			select(Sample, Dataset, Group, Kit, dplyr::everything()) %>%
			setDT()

		df_long <- data.table::melt(
			X_imp,
			id.vars = c("Sample", "Dataset", "Group", "Kit"),
			variable.name = "Gene",
			value.name = "CADD_weighted_avg_AF"
		)

		data.table::fwrite(
			X_imp,
			file.path(results_dir, paste0(base, "_wide.tsv")),
			sep = "\t"
		)

		data.table::fwrite(
			df_long,
			file.path(results_dir, paste0(base, "_long.tsv")),
			sep = "\t"
		)

		# STEP 5: Postprocessing =============================================
		# batch metrics
		cat("Computing batch metrics...\n")

		batch_metric_test_size <- batch_metric_params$test_size_fraction *
			length(unique(df_long$Sample))

		batch_stats <- compute_batch_metrics_df(
			df_long %>%
				select(Sample, Dataset, Gene, CADD_weighted_avg_AF),
			lisi_perplexity = batch_metric_params$lisi_perplexity,
			k_kBET = batch_metric_params$k_kBET,
			test_size = batch_metric_test_size,
			heuristic_kBET = batch_metric_params$heuristic_kBET,
			adapt_kBET = batch_metric_params$adapt_kBET,
			PCA_kBET = batch_metric_params$PCA_kBET
		)

		cat(paste("LISI:", batch_stats$lisi_stats$mean, "\n"))
		cat(paste("kBET:", batch_stats$kbet_stats$mean, "\n"))

		# UMAP
		cat("Running UMAP...\n")

		umap_input <- X_imp %>%
			select(-Sample, -Dataset, -Group, -Kit) %>%
			as.matrix()

		metadata <- X_imp %>%
			select(Sample, Dataset, Group, Kit)

		set.seed(umap_params$seed)

		umap_result <- uwot::umap(
			umap_input,
			n_neighbors = umap_params$n_neighbors,
			min_dist = umap_params$min_dist,
			metric = umap_params$metric,
			n_epochs = umap_params$n_epochs,
			n_threads = parallel::detectCores(),
			ret_model = FALSE, # do not return the model
			nn_method = umap_params$nn_method
		)

		umap_result <- as.data.frame(umap_result)
		colnames(umap_result) <- c("UMAP1", "UMAP2")
		umap_result <- cbind(umap_result, metadata)
		umap_result$Group <- factor(umap_result$Group, levels = sort(unique(umap_result$Group)))

		data.table::fwrite(
			umap_result,
			file.path(results_dir, paste0(base, "_umap_result.tsv")),
			sep = "\t"
		)

		# Plot the UMAP ------------------------------------------------------
		cat("Making the figures...\n")

		umap_caption <- paste0(
			"UMAP params: n_neighbors=", umap_params$n_neighbors,
			", min_dist=", umap_params$min_dist,
			", metric=", umap_params$metric,
			", n_epochs=", umap_params$n_epochs
		)

		plot_title <- paste0(
			"Feature Imputation\na=", round(a, 3),
			" b=", round(b, 3),
			" kn=", k,
			"\nLISI=", round(batch_stats$lisi_stats$mean, 2),
			" | kBET=", round(batch_stats$kbet_stats$mean, 3)
		)

		umap_plot <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = Dataset)) +
			geom_point() +
			theme_test() +
			labs(
				title = plot_title,
				x = "UMAP 1",
				y = "UMAP 2",
				caption = umap_caption
			) +
			scale_color_hue() + coord_equal()

		print(umap_plot)

		pdf(
			file.path(figures_dir, paste0(gsub("ft_imp_", "", base), "_UMAP.pdf")),
			width = 7,
			height = 6
		)
		print(umap_plot)
		dev.off()

		umap_plot2 <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = as.factor(Group))) +
			geom_point() +
			theme_test() +
			labs(
				title = plot_title,
				x = "UMAP 1",
				y = "UMAP 2",
				caption = umap_caption,
				color = "Cluster"
			) +
			scale_color_hue() + coord_equal()

		print(umap_plot2)

		pdf(
			file.path(figures_dir, paste0(gsub("ft_imp_", "", base), "_UMAP_cluster_colored.pdf")),
			width = 6,
			height = 6
		)
		print(umap_plot2)
		dev.off()

		cat("UMAP: DONE\n")
	}
}
