#' Impute MNAR values using masked cosine kNN across groups
#'
#' For each sample and each MNAR entry, find similar donor samples
#' from other groups using masked cosine similarity on all other genes,
#' then impute a weighted average from their observed values.
#'
#' @param X_raw Numeric matrix [samples x genes]. May contain NA.
#' @param MNAR Logical matrix same size as X_raw. TRUE = MNAR flag -> impute this entry.
#' @param grp Factor/character, length = nrow(X_raw). Group or batch label.
#' @param k Integer. Max number of neighbors to use.
#' @param min_k Integer. Min effective neighbors needed for weighted mean.
#' @param min_overlap Integer. Min shared non-missing features required for cosine similarity.
#' @param workers Integer. Number of parallel workers.
#'
#' @return Numeric matrix same dim as X_raw with imputed values where MNAR==TRUE.
#' @export
#'
knn_impute_mnar_masked_parallel <- function(X_raw, MNAR, grp, k = 10, min_k = 5, min_overlap = 50, workers = 8) {
  stopifnot(all(dim(X_raw) == dim(MNAR)), length(grp) == nrow(X_raw))
  Xr <- as.matrix(X_raw)
  n_samples <- nrow(Xr)

  cos_masked <- function(x, Y, min_overlap) {
    sapply(seq_len(nrow(Y)), function(i) {
      ok <- is.finite(x) & is.finite(Y[i, ])
      if (sum(ok) < min_overlap) return(NA_real_)
      xx <- x[ok]; yy <- Y[i, ok]
      den <- sqrt(sum(xx * xx)) * sqrt(sum(yy * yy))
      if (den == 0) return(NA_real_)
      sum(xx * yy) / den
    })
  }

  future::plan(future::multicore, workers = workers)

  rows_out <- future.apply::future_lapply(seq_len(n_samples), function(i) {
    x <- Xr[i, ]
    mnar_mask <- MNAR[i, ]
    if (!any(mnar_mask)) return(list(vec = x, zeros = 0L, no_donors = 0L, too_few_sim = 0L, too_few_vals = 0L))

    donors <- which(grp != grp[i])
    zeros_local <- 0L; no_donors <- 0L; too_few_sim <- 0L; too_few_vals <- 0L

    if (length(donors) < 2) {
      x[mnar_mask] <- 0
      no_donors <- sum(mnar_mask)
      return(list(vec = x, zeros = no_donors, no_donors = no_donors, too_few_sim = 0L, too_few_vals = 0L))
    }

    similarities <- cos_masked(x, Xr[donors, ], min_overlap)
    for (j in which(mnar_mask)) {
      keep <- which(is.finite(similarities) & is.finite(Xr[donors, j]))
      if (length(keep) < min_k) {
        x[j] <- 0
        zeros_local <- zeros_local + 1L
        next
      }
      ord <- order(similarities[keep], decreasing = TRUE)
      top_idx <- donors[keep][ord][seq_len(min(k, length(ord)))]
      w <- similarities[keep][ord][seq_len(min(k, length(ord)))] / 
          sum(similarities[keep][ord][seq_len(min(k, length(ord)))])
      w_norm <- w / sum(w)
      x[j] <- sum(w_norm * Xr[top_idx, j])
      if (!is.finite(x[j])) {
        too_few_vals <- too_few_vals + 1L
        x[j] <- 0
}


    }

    list(vec = x, zeros = zeros_local, no_donors = no_donors,
         too_few_sim = too_few_sim, too_few_vals = too_few_vals)
  })

  future::plan(future::sequential)

  X_out <- do.call(rbind, lapply(rows_out, `[[`, "vec"))
  zero_total <- sum(vapply(rows_out, `[[`, integer(1), "zeros"))
  no_donors_total <- sum(vapply(rows_out, `[[`, integer(1), "no_donors"))
  too_few_sim_total <- sum(vapply(rows_out, `[[`, integer(1), "too_few_sim"))
  too_few_vals_total <- sum(vapply(rows_out, `[[`, integer(1), "too_few_vals"))

  rownames(X_out) <- rownames(X_raw)
  colnames(X_out) <- colnames(X_raw)

  message(sprintf(
    "Total zeros imputed: %d (no donors: %d | too few similar: %d | too few vals: %d)",
    zero_total, no_donors_total, too_few_sim_total, too_few_vals_total
  ))
  X_out[is.na(X_out) & !MNAR] <- 0

  X_out
}

