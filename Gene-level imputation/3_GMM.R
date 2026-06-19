#!/usr/bin/env Rscript

# This script models gene detection rates across sample clusters.
# It loads the long-format feature table and sample-to-cluster assignments,
# computes detection rates per gene and cluster, fits a 1D Gaussian mixture model,
# and saves the fitted model and diagnostic plot.
# Usage:
# R_LIBS= R_LIBS_USER= R_LIBS_SITE= R_PROFILE_USER=/dev/null \
# Rscript --vanilla "Gene-level imputation/3_GMM.R" config/local_config.yaml

conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library")
if (nzchar(Sys.getenv("CONDA_PREFIX")) && dir.exists(conda_lib)) {
	.libPaths(conda_lib)
}

suppressPackageStartupMessages({
	library(data.table)
	library(ggplot2)
	library(dpGMM)
})

threshold_label_y <- 15
threshold_digits <- 3

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

results_dir <- resolve_path(read_config("directories.results_dir"))

output_dir <- file.path(results_dir, "Gene_level_imputation")
figures_dir <- file.path(output_dir, "Figures")

features_file <- file.path(output_dir, "raw_ft_long.tsv")
clusters_file <- file.path(output_dir, "sample_kit_cluster_map.tsv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(features_file)) {
	stop("Long-format feature table not found: ", features_file)
}

if (!file.exists(clusters_file)) {
	stop("Sample-to-cluster mapping not found: ", clusters_file)
}

ft_long <- fread(features_file)
clusters <- fread(clusters_file)

required_feature_cols <- c("Sample", "Gene", "CADD_weighted_avg_AF")
missing_feature_cols <- setdiff(required_feature_cols, colnames(ft_long))

if (length(missing_feature_cols) > 0) {
	stop("Missing feature columns: ", paste(missing_feature_cols, collapse = ", "))
}

required_cluster_cols <- c("Sample", "Cluster")
missing_cluster_cols <- setdiff(required_cluster_cols, colnames(clusters))

if (length(missing_cluster_cols) > 0) {
	stop("Missing cluster columns: ", paste(missing_cluster_cols, collapse = ", "))
}

ft_long[, Group := clusters$Cluster[match(Sample, clusters$Sample)]]
ft_long <- ft_long[!is.na(Group)]
ft_long[, Group := as.factor(Group)]

if (nrow(ft_long) == 0) {
	stop("No feature records matched cluster assignments.")
}

sample_group <- unique(ft_long[, .(Sample, Group)])
n_total <- sample_group[, .(Total = .N), by = Group]

detected <- unique(
	ft_long[CADD_weighted_avg_AF > 0, .(Sample, Group, Gene)]
)

det_rate <- detected[, .(Detected = .N), by = .(Group, Gene)]
det_rate <- merge(det_rate, n_total, by = "Group")
det_rate[, Detection_rate := Detected / Total]

if (nrow(det_rate) == 0) {
	stop("No detected gene entries found for GMM fitting.")
}

custom_settings <- GMM_1D_opts

mix_test <- runGMM(det_rate$Detection_rate, opts = custom_settings)

saveRDS(
	mix_test,
	file.path(output_dir, "GMM_det_rate.RDS")
)

thresholds <- mix_test[["threshold"]]
base_plot <- mix_test[["fig"]]

plot_data <- ggplot_build(base_plot)$data
y_values <- unlist(lapply(plot_data, function(layer) layer$y))
y_max <- max(y_values, na.rm = TRUE)

threshold_df <- data.frame(
	x = thresholds,
	y = 0.8 * y_max
)

p <- base_plot +
	geom_text(
		data = threshold_df,
		aes(x = x, y = y, label = round(x, 3)),
		angle = 90,
		vjust = -0.5,
		hjust = 0,
		color = "red",
		size = 3
	) +
	labs(x = "Detection rate")

pdf(file.path(figures_dir, "dpGMM_detection_rate_plot.pdf"), width = 6, height = 4)
print(p)
dev.off()

cat("Detection-rate GMM completed.\n")

