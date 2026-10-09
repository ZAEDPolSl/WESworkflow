
#!/usr/bin/env Rscript

# Generate UMAP and batch metrics from an existing gene-level feature table.
# Intended for the artificial downstream example dataset.
# Usage:
# Rscript "Data/example_downstream/umap_from_raw_ft_long.R" config/local_config.yaml

conda_lib <- file.path(Sys.getenv("CONDA_PREFIX"), "lib/R/library")
if (nzchar(Sys.getenv("CONDA_PREFIX")) && dir.exists(conda_lib)) {
	.libPaths(conda_lib)
}

# ============ Adjustable parameters ============
plot_title <- "Artificial downstream example dataset"

umap_params <- list(
	n_neighbors = 15,
	min_dist = 0.5,
	metric = "cosine",
	n_epochs = 1500,
	nn_method = "nndescent",
	seed = 123
)

batch_metric_params <- list(
	lisi_perplexity = 30,
	k_kBET = 15,
	test_size_fraction_kBET = 0.1,
	heuristic_kBET = FALSE,
	adapt_kBET = FALSE,
	PCA_kBET = TRUE
)
# ===============================================

suppressPackageStartupMessages({
	library(data.table)
	library(ggplot2)
	library(rnndescent)
	library(uwot)
	library(dplyr)
})

get_script_path <- function() {
	args <- commandArgs(trailingOnly = FALSE)
	file_arg <- grep("^--file=", args, value = TRUE)
	if (!length(file_arg)) stop("Cannot determine script path.")
	normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg[1])), mustWork = TRUE)
}

script_dir <- dirname(get_script_path())
repo_dir <- normalizePath(file.path(script_dir, "../.."), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
config <- if (length(args)) args[1] else "config/local_config.yaml"
if (!grepl("^/", config)) config <- file.path(repo_dir, config)
config <- normalizePath(config, mustWork = TRUE)

read_config <- function(key) {
	value <- system2("python", args = c(file.path(repo_dir, "scripts/read_config.py"), config, key), stdout = TRUE)
	if (!is.null(attr(value, "status")) || !length(value) || is.na(value[1]) || !nzchar(value[1])) {
		stop("Failed to read config key: ", key)
	}
	value[1]
}

resolve_path <- function(path) {
	if (!length(path) || is.na(path) || !nzchar(path)) stop("Cannot resolve an empty path.")
	if (grepl("^/", path)) return(path)
	file.path(repo_dir, path)
}

results_dir <- resolve_path(read_config("directories.results_dir"))
metadata_file <- resolve_path(read_config("directories.sample_metadata"))
feature_column <- read_config("parameters.feature_column")

output_dir <- file.path(results_dir, "Gene_level_imputation", feature_column)
results_dir <- file.path(output_dir, "Results")
figures_dir <- file.path(output_dir, "Figures")
raw_ft_long_file <- file.path(results_dir, "raw_ft_long.tsv")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(raw_ft_long_file)) {
	stop("raw_ft_long.tsv not found: ", raw_ft_long_file)
}

features <- fread(raw_ft_long_file)

required_feature_cols <- c("Sample", "Gene", feature_column)
missing_feature_cols <- setdiff(required_feature_cols, names(features))
if (length(missing_feature_cols)) {
	stop("Missing columns in raw_ft_long.tsv: ", paste(missing_feature_cols, collapse = ", "))
}

if (!all(c("Dataset", "Kit") %in% names(features))) {
	if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)

	meta <- fread(metadata_file)
	required_meta_cols <- c("Sample", "Dataset", "Kit")
	missing_meta_cols <- setdiff(required_meta_cols, names(meta))

	if (length(missing_meta_cols)) {
		stop("Missing metadata columns: ", paste(missing_meta_cols, collapse = ", "))
	}

	meta <- unique(meta[, ..required_meta_cols])
	if (anyDuplicated(meta$Sample)) stop("Duplicated sample IDs in metadata.")

	features <- merge(features, meta, by = "Sample", all.x = TRUE)
} else {
	meta <- unique(features[, .(Sample, Dataset, Kit)])
	if (anyDuplicated(meta$Sample)) stop("Conflicting sample metadata in raw_ft_long.tsv.")
}

if (anyNA(features[, .(Sample, Gene, Dataset, Kit)]) ||
	any(features$Dataset == "" | features$Kit == "")) {
	stop("Missing sample, gene, dataset or capture kit information.")
}

duplicates <- features[, .N, by = .(Sample, Gene)][N > 1]
if (nrow(duplicates)) {
	print(head(duplicates, 20))
	stop("Duplicated Sample-Gene combinations found.")
}

if (anyNA(features[[feature_column]])) {
	stop("Missing feature values in column: ", feature_column)
}

cat("Selected feature: ", feature_column, "\n", sep = "")
cat("Loaded raw_ft_long.tsv with:\n")
cat("  samples: ", uniqueN(features$Sample), "\n", sep = "")
cat("  genes: ", uniqueN(features$Gene), "\n", sep = "")
cat("  rows: ", nrow(features), "\n", sep = "")

cat("Reshaping data to wide format...\n")

umap_input <- features[, c("Sample", "Gene", feature_column), with = FALSE]
umap_input <- dcast(umap_input, Sample ~ Gene, value.var = feature_column, fill = 0)

sample_ids <- umap_input$Sample
gene_cols <- setdiff(names(umap_input), "Sample")
umap_matrix <- as.matrix(umap_input[, ..gene_cols])
rownames(umap_matrix) <- sample_ids

n_samples <- nrow(umap_matrix)
umap_n_neighbors <- min(umap_params$n_neighbors, n_samples - 1)
if (umap_n_neighbors < 2) stop("At least 3 samples are required to run UMAP.")

set.seed(umap_params$seed)

cat("Running UMAP with n_neighbors=", umap_n_neighbors,
	" for ", n_samples, " samples...\n", sep = "")

umap_coords <- uwot::umap(
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
	UMAP1 = umap_coords[, 1],
	UMAP2 = umap_coords[, 2]
)

umap_result <- merge(umap_result, meta, by = "Sample", all.x = TRUE)
fwrite(umap_result, file.path(results_dir, "raw_umap_result.tsv"), sep = "\t")

source(file.path(repo_dir, "Gene-level imputation", "functions", "batch_metrics.R"))

batch_metric_test_size <- batch_metric_params$test_size_fraction_kBET * uniqueN(features$Sample)

batch_stats <- tryCatch(
	compute_batch_metrics_df(
		features,
		feature_column = feature_column,
		lisi_perplexity = batch_metric_params$lisi_perplexity,
		k_kBET = batch_metric_params$k_kBET,
		test_size = batch_metric_test_size,
		heuristic_kBET = batch_metric_params$heuristic_kBET,
		adapt_kBET = batch_metric_params$adapt_kBET,
		PCA_kBET = batch_metric_params$PCA_kBET
	),
	error = function(e) {
		warning("Batch metrics could not be computed: ", conditionMessage(e))
		list(
			lisi_stats = data.frame(mean = NA_real_, median = NA_real_, sd = NA_real_),
			kbet_stats = data.frame(mean = NA_real_, median = NA_real_, sd = NA_real_),
			lisi_raw = NULL,
			kbet_raw = NULL
		)
	}
)

format_metric <- function(x, digits) {
	if (is.null(x) || !length(x) || is.na(x)) return("NA")
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

umap_plot_dataset <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = Dataset)) +
	geom_point(size = 1, alpha = 0.5) +
	theme_test() +
	labs(
		title = plot_title,
		x = "UMAP 1",
		y = "UMAP 2",
		color = "Dataset",
		subtitle = batch_subtitle,
		caption = umap_caption
	) +
	colorspace::scale_color_discrete_qualitative(palette = "Dark3") +
	theme(aspect.ratio = 1) +
	guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.6)))

umap_plot_kit <- ggplot(umap_result, aes(x = UMAP1, y = UMAP2, color = Kit)) +
	geom_point(size = 1, alpha = 0.5) +
	theme_test() +
	labs(
		title = plot_title,
		x = "UMAP 1",
		y = "UMAP 2",
		color = "Capture kit",
		subtitle = batch_subtitle,
		caption = umap_caption
	) +
	colorspace::scale_color_discrete_qualitative(palette = "Dark3") +
	theme(aspect.ratio = 1) +
	guides(color = guide_legend(override.aes = list(size = 3, alpha = 0.6)))

pdf(file.path(figures_dir, "raw_UMAP.pdf"), width = 7, height = 6)
print(umap_plot_dataset)
print(umap_plot_kit)
dev.off()

cat("UMAP generation from existing raw_ft_long.tsv completed.\n")
