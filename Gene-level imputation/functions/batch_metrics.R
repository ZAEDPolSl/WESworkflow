
compute_batch_metrics_df <- function(df, feature_column = "CADD_weighted_avg_AF", lisi_perplexity = 30, k_kBET = 15, test_size = 0.1*length(unique(df$Sample)), heuristic_kBET = FALSE, adapt_kBET = FALSE, PCA_kBET = TRUE) {
	stopifnot(all(c(feature_column, "Sample", "Dataset") %in% names(df)))

  if ("Gene" %in% names(df)) {df <- df %>% 
    rename(gene = Gene)}

  wide <- df |>
    dplyr::select(gene, dplyr::all_of(feature_column), Sample, Dataset) |>
    dplyr::group_by(Sample, Dataset, gene) |>
    dplyr::summarise(value = mean(.data[[feature_column]]), .groups = "drop") |>
    tidyr::pivot_wider(names_from = gene, values_from = value, values_fill = 0)

  meta <- wide |> dplyr::select(Sample, Dataset)
  X <- as.matrix(wide |> dplyr::select(-Sample, -Dataset))

  if (any(!is.finite(X))) stop("NA/NaN/Inf in matrix X.")
  n <- nrow(X)
  if (n <= 3 * lisi_perplexity) stop("n too small for given perplexity.")
  if (n <= (k_kBET + 1))           stop("n too small for given k_kBET.")

  # LISI
  lisi_res <- lisi::compute_lisi(X, meta, "Dataset", perplexity = lisi_perplexity)
  lisi_stats <- tibble::as_tibble(lisi_res) |>
    dplyr::summarise(mean = mean(Dataset), median = median(Dataset), sd = sd(Dataset))

  # kBET
  knn <- RANN::nn2(X, k = k_kBET + 1)$nn.idx[, -1, drop = FALSE]
  kb  <- kBET::kBET(X, meta$Dataset, k0 = k_kBET, knn = knn,
                    heuristic = heuristic_kBET, testSize = min(test_size, n), adapt = adapt_kBET, verbose = FALSE, do.pca = PCA_kBET)

 
  obs <- tryCatch(kb$stats$kBET.observed, error = function(e) NULL)

  kbet_stats <- tibble::tibble(
    mean   = mean(obs, na.rm = TRUE),
    median = median(obs, na.rm = TRUE),
    sd     = sd(obs, na.rm = TRUE)
  )
  
  list(lisi_stats = lisi_stats, kbet_stats = kbet_stats, lisi_raw = lisi_res, kbet_raw = kb)
}

