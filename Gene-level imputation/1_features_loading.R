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

# ============ these parameters can be adjusted ============
# ==== but should be consistent across downstream steps ===+
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
	library(rnndescent)
	library(uwot)
	library(tidyr)
	library(dplyr)
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

results_dir <- resolve_path(read_config("directories.results_dir"))
metadata_file <- resolve_path(read_config("directories.sample_metadata"))

features_dir <- file.path(results_dir, "Features")
output_dir <- file.path(results_dir, "Gene_level_imputation")
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

n_samples <- nrow(umap_matrix)
umap_n_neighbors <- min(umap_params$n_neighbors, n_samples - 1)

if (umap_n_neighbors < 2) {
	stop("At least 3 samples are required to run UMAP.")
}

cat(
	"Running UMAP with n_neighbors=",
	umap_n_neighbors,
	" for ",
	n_samples,
	" samples...\n",
	sep = ""
)

umap_result <- uwot::umap(
	umap_matrix,
	n_neighbors = umap_n_neighbors,
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

batch_metric_test_size <- batch_metric_params$test_size_fraction_kBET *
	length(unique(features$Sample))

batch_stats <- tryCatch(
	compute_batch_metrics_df(
		features,
		lisi_perplexity = batch_metric_params$lisi_perplexity,
		k_kBET = batch_metric_params$k_kBET,
		test_size = batch_metric_test_size,
		heuristic_kBET = batch_metric_params$heuristic_kBET,
		adapt_kBET = batch_metric_params$adapt_kBET,
		PCA_kBET = batch_metric_params$PCA_kBET
	),
	error = function(e) {
		warning(
			"Batch metrics could not be computed: ",
			conditionMessage(e)
		)

		list(
			lisi_stats = data.frame(
				mean = NA_real_,
				median = NA_real_,
				sd = NA_real_
			),
			kbet_stats = data.frame(
				mean = NA_real_,
				median = NA_real_,
				sd = NA_real_
			),
			lisi_raw = NULL,
			kbet_raw = NULL
		)
	}
)

format_metric <- function(x, digits) {
	if (is.null(x) || length(x) == 0 || is.na(x)) {
		return("NA")
	}

	round(x, digits)
}

batch_subtitle <- paste0(
	"LISI = ", format_metric(batch_stats$lisi_stats$mean, 2),
	" | kBET = ", format_metric(batch_stats$kbet_stats$mean, 3)
)

umap_caption <- paste0(
	"n_neighbors=", umap_n_neighbors,
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
		subtitle = batch_subtitle,
		caption = umap_caption
	) + coord_equal()

umap_plot_kit <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = Kit)) +
	geom_point(size = 1) +
	theme_test() +
	labs(
		title = plot_title,
		x = "UMAP 1",
		y = "UMAP 2",
		color = "Capture kit",
		subtitle = batch_subtitle,
		caption = umap_caption
	) + coord_equal()

pdf(file.path(figures_dir, "raw_UMAP.pdf"), width = 7, height = 6)
print(umap_plot)
print(umap_plot_kit)
dev.off()



cat("Feature loading and UMAP analysis completed.\n")
