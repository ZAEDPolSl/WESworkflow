#!/usr/bin/env Rscript

# This script loads gene-level features, merges them with sample metadata,
# computes UMAP embeddings, estimates batch effect metrics, and saves plots.
# Usage:
# R_LIBS= R_LIBS_USER= R_LIBS_SITE= R_PROFILE_USER=/dev/null \
# Rscript --vanilla "Gene-level imputation/1_features_loading.R" config/local_config.yaml

conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library")
if (nzchar(Sys.getenv("CONDA_PREFIX")) && dir.exists(conda_lib)) {
	.libPaths(conda_lib)
}

suppressPackageStartupMessages({
	library(data.table)
	library(ggplot2)
	library(rnndescent)
	library(uwot)
	library(tidyr)
	library(dplyr)
})

# UMAP parameters - these values can be adjusted 
metadata_file_rel <- "repo_path/Data/example/sample_path_map_example.tsv" # example file of metadata

umap_params <- list(
	n_neighbors = 15,
	min_dist = 0.5,
	metric = "cosine",
	n_epochs = 1500,
	nn_method = "nndescent",
	seed = 123
)

plot_title <- "Data after genotype imputation"

get_script_path <- function() {
	args <- commandArgs(trailingOnly = FALSE)
	file_arg <- grep("^--file=", args, value = TRUE)

	if (length(file_arg) == 0) {
		stop("Cannot determine script path.")
	}

	normalizePath(sub("^--file=", "", file_arg[1]))
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
repo_dir <- normalizePath(file.path(script_dir, ".."))

args <- commandArgs(trailingOnly = TRUE)
config <- if (length(args) >= 1) args[1] else "config/local_config.yaml"

if (!grepl("^/", config)) {
	config <- file.path(repo_dir, config)
}

read_config <- function(key) {
	system2(
		"python",
		args = c(file.path(repo_dir, "scripts/read_config.py"), config, key),
		stdout = TRUE
	)
}

resolve_path <- function(path) {
	if (grepl("^/", path)) {
		return(path)
	}

	file.path(repo_dir, path)
}

results_dir <- resolve_path(read_config("directories.results_dir"))
metadata_file <- resolve_path(metadata_file_rel)

features_dir <- file.path(results_dir, "Features")
output_dir <- file.path(results_dir, "Results")
figures_dir <- file.path(output_dir, "Figures")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(features_dir)) {
	stop("Feature directory not found: ", features_dir)
}

if (!file.exists(metadata_file)) {
	stop("Metadata file not found: ", metadata_file)
}

load_chr_features <- function(feature_files_dir) {
	files <- list.files(
		feature_files_dir,
		pattern = "\\.feature\\.txt$",
		recursive = TRUE,
		full.names = TRUE
	)

	if (length(files) == 0) {
		stop("No .feature.txt files found in: ", feature_files_dir)
	}

	all_features <- rbindlist(lapply(files, function(f) {
		dt <- fread(f, select = c("Gene", "CADD_weighted_avg_AF"))
		sample_name <- sub("\\.feature\\.txt$", "", basename(f))

		dt[, Sample := sample_name]
		dt
	}), use.names = TRUE, fill = TRUE)

	all_features
}

meta <- fread(metadata_file)

required_cols <- c("Sample", "Dataset", "Kit")
missing_cols <- setdiff(required_cols, colnames(meta))

if (length(missing_cols) > 0) {
	stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
}

meta <- meta[, ..required_cols]

features <- load_chr_features(features_dir)
features <- merge(features, meta, by = "Sample", all.x = TRUE)

fwrite(
	features,
	file = file.path(output_dir, "raw_ft_long.tsv"),
	sep = "\t"
)

cat("Reshaping data to wide format...\n")

umap_input <- features[, .(Sample, Gene, CADD_weighted_avg_AF)]

umap_input <- dcast(
	umap_input,
	Sample ~ Gene,
	value.var = "CADD_weighted_avg_AF",
	fill = 0,
	fun.aggregate = sum
)

sample_ids <- umap_input$Sample
gene_cols <- setdiff(colnames(umap_input), "Sample")
umap_matrix <- as.matrix(umap_input[, ..gene_cols])
rownames(umap_matrix) <- sample_ids

set.seed(umap_params$seed)

cat("Running UMAP...\n")

umap_result <- uwot::umap(
	umap_matrix,
	n_neighbors = umap_params$n_neighbors,
	min_dist = umap_params$min_dist,
	metric = umap_params$metric,
	n_epochs = umap_params$n_epochs,
	n_threads = parallel::detectCores(),
	ret_model = FALSE,
	nn_method = umap_params$nn_method
)

umap_result <- data.table(
	Sample = sample_ids,
	UMAP1 = umap_result[, 1],
	UMAP2 = umap_result[, 2]
)

umap_result <- merge(umap_result, meta, by = "Sample", all.x = TRUE)

fwrite(
	umap_result,
	file.path(output_dir, "raw_umap_result.tsv"),
	sep = "\t"
)

source(file.path(script_dir, "functions", "batch_metrics.R"))

batch_stats <- compute_batch_metrics_df(features)

umap_caption <- paste0(
	"n_neighbors=", umap_params$n_neighbors,
	", min_dist=", umap_params$min_dist,
	", metric=", umap_params$metric,
	", n_epochs=", umap_params$n_epochs
)

umap_plot <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = Dataset)) +
	geom_point(size = 1) +
	theme_test() +
	labs(
		title = plot_title,
		x = "UMAP 1",
		y = "UMAP 2",
		subtitle = paste0(
			"LISI = ", round(batch_stats$lisi_stats$mean, 2),
			" | kBET = ", round(batch_stats$kbet_stats$mean, 3)
		),
		caption = umap_caption
	)

umap_plot_kit <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = Kit)) +
	geom_point(size = 1) +
	theme_test() +
	labs(
		title = plot_title,
		x = "UMAP 1",
		y = "UMAP 2",
		color = "Capture kit",
		subtitle = paste0(
			"LISI = ", round(batch_stats$lisi_stats$mean, 2),
			" | kBET = ", round(batch_stats$kbet_stats$mean, 3)
		),
		caption = umap_caption
	)

pdf(file.path(figures_dir, "raw_UMAP.pdf"), width = 7, height = 6)
print(umap_plot)
print(umap_plot_kit)
dev.off()

cat("Feature loading and UMAP analysis completed.\n")
