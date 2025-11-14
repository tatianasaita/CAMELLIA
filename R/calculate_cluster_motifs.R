#' Calculate Cluster Motifs from K-mer Counts
#'
#' Calculates motif profiles for each cluster by summing k-mer counts of all
#' sequences in the cluster. The resulting motif matrix is normalized using
#' min-max normalization applied independently to each cluster column.
#'
#' @param x Object of class 'cluster_dendrogram_result' from cluster_dendrogram().
#' @param data_result Optional. List containing 'kmers' and 'metadata' data.frames.
#'   If NULL, extracted from x$data_result. Output from create_data() with
#'   cluster assignments added by cluster_dendrogram().
#' @param ... Additional arguments passed to methods.
#'
#' @return A data.frame with k-mer motifs by cluster (normalized).
#'   Rows are k-mers, columns are clusters. Values are normalized between 0 and 1
#'   independently for each cluster column.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item \strong{Validation:} Checks that input contains cluster assignments
#'     in metadata
#'   \item \strong{Extraction:} Gets k-mer counts and cluster assignments
#'   \item \strong{Aggregation:} For each cluster, sums k-mer values across
#'     all member sequences
#'   \item \strong{Matrix construction:} Creates a matrix with k-mers as rows
#'     and clusters as columns
#'   \item \strong{Normalization:} Applies min-max normalization (0-1 range)
#'     independently to each cluster column
#' }
#'
#' Min-max normalization formula (applied per column):
#' \code{x_normalized = (x - min(column)) / (max(column) - min(column))}
#'
#' If all values in a column are identical (min == max), they are set to 0.5
#' to avoid division by zero.
#'
#' @examples
#' \dontrun{
#' # After running create_data, create_dendrogram, and cluster_dendrogram:
#' motif_cluster <- calculate_cluster_motifs(cluster_result)
#'
#' # View normalized motif matrix
#' head(motif_cluster)
#'
#' # Visualize motifs
#' heatmap(as.matrix(motif_cluster))
#' }
#'
#' @export
calculate_cluster_motifs <- function(x, ...) {
  UseMethod("calculate_cluster_motifs")
}

#' @export
calculate_cluster_motifs.default <- function(x, ...) {
  stop("'x' must be an object of class 'cluster_dendrogram_result'")
}
#' @export
calculate_cluster_motifs.cluster_dendrogram_result <- function(
    x,
    data_result = NULL,
    ...
) {

  # ===== INPUT VALIDATION =====

  if (!inherits(x, "cluster_dendrogram_result")) {
    stop("'x' must be an object of class 'cluster_dendrogram_result'")
  }

  if (is.null(data_result)) {
    if (!("data_result" %in% names(x))) {
      stop("'x' must contain 'data_result' element from cluster_dendrogram()")
    }
    data_result <- x$data_result
  }

  if (!is.list(data_result)) {
    stop("'data_result' must be a list")
  }

  if (!all(c("kmers", "metadata") %in% names(data_result))) {
    stop("'data_result' must contain 'kmers' and 'metadata' elements")
  }

  if (!is.data.frame(data_result$kmers)) {
    stop("'data_result$kmers' must be a data.frame")
  }

  if (!is.data.frame(data_result$metadata)) {
    stop("'data_result$metadata' must be a data.frame")
  }

  if (!"cluster" %in% colnames(data_result$metadata)) {
    stop("'data_result$metadata' must contain a 'cluster' column")
  }

  if (nrow(data_result$metadata) != nrow(data_result$kmers)) {
    stop("'kmers' and 'metadata' must have the same number of rows")
  }

  # ===== EXTRACT DATA =====

  kmers_df <- data_result$kmers
  metadata <- data_result$metadata
  cluster_assignments <- metadata$cluster

  # ===== PREPARE K-MER DATA =====

  # Find and exclude CLASS column
  class_col_idx <- which(tolower(colnames(kmers_df)) == "class")

  if (length(class_col_idx) > 0) {
    kmers_numeric <- kmers_df[, -class_col_idx, drop = FALSE]
  } else {
    kmers_numeric <- kmers_df
  }

  # Ensure all columns are numeric
  for (i in seq_len(ncol(kmers_numeric))) {
    if (!is.numeric(kmers_numeric[[i]])) {
      stop("All k-mer columns must be numeric")
    }
  }

  # Get numeric column names (k-mers)
  numeric_cols <- colnames(kmers_numeric)

  # Get unique clusters and cluster count
  cluster_ids <- unique(cluster_assignments)
  n_clusters <- length(cluster_ids)

  # ===== AGGREGATE K-MERS BY CLUSTER (BASE R) =====

  motif_raw <- matrix(
    nrow = length(numeric_cols),
    ncol = n_clusters,
    dimnames = list(numeric_cols, paste0("Cluster_", cluster_ids))
  )

  for (i in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[i]
    cluster_mask <- cluster_assignments == cluster_id

    if (sum(cluster_mask) == 0) {
      stop(sprintf("Cluster %d has no sequences", cluster_id))
    }

    cluster_kmers <- kmers_numeric[cluster_mask, , drop = FALSE]
    motif_raw[, i] <- colSums(cluster_kmers)
  }

  # ===== MIN-MAX NORMALIZATION (per column) =====

  col_mins <- apply(motif_raw, 2, min)
  col_maxs <- apply(motif_raw, 2, max)
  col_ranges <- col_maxs - col_mins

  # Identify columns with zero range (all values identical)
  zero_range <- col_ranges == 0
  col_ranges[zero_range] <- 1  # Avoid division by zero

  # Apply min-max normalization using sweep (vectorized)
  motif_normalized <- sweep(
    sweep(motif_raw, 2, col_mins, "-"),
    2,
    col_ranges,
    "/"
  )

  # Set zero-range columns to 0.5
  motif_normalized[, zero_range] <- 0.5

  # ===== CREATE OUTPUT DATA.FRAME =====

  motif_cluster <- as.data.frame(motif_normalized)

  # ===== PRINT SUMMARY =====

  cat("\n=== Cluster Motif Calculation Summary ===\n")
  cat("Total clusters:       ", n_clusters, "\n")
  cat("Total k-mers:         ", length(numeric_cols), "\n")
  cat("Motif matrix dims:    ", nrow(motif_cluster), "x", ncol(motif_cluster), "\n")
  cat("Normalization:        Min-Max (0-1) per cluster\n")
  cat("Zero-range columns:   ", sum(zero_range), "\n")
  cat("Aggregation method:   base R\n")
  cat("\n")

  # ===== RETURN RESULTS =====

  class(motif_cluster) <- c("cluster_motifs", "data.frame")
  return(invisible(motif_cluster))
}


#' Print Method for Cluster Motifs
#'
#' @param x Object of class 'cluster_motifs'
#' @param ... Additional arguments passed to print()
#'
#' @export
print.cluster_motifs <- function(x, ...) {
  cat("Cluster Motif Matrix\n")
  cat("====================\n")
  cat("Dimensions:", nrow(x), "k-mers x", ncol(x), "clusters\n")
  cat("Normalization: Min-Max (0-1 range)\n")
  cat("Value range: [0, 1]\n\n")

  print.data.frame(x, ...)
}


#' Summary Method for Cluster Motifs
#'
#' @param object Object of class 'cluster_motifs'
#' @param ... Additional arguments (unused)
#'
#' @export
summary.cluster_motifs <- function(object, ...) {
  cat("Cluster Motif Summary\n")
  cat("=====================\n\n")

  cat("Matrix dimensions:", nrow(object), "k-mers x", ncol(object), "clusters\n\n")

  cat("Per-cluster statistics:\n")
  cat("-----------------------\n")

  for (col_idx in seq_len(ncol(object))) {
    col_name <- colnames(object)[col_idx]
    col_vals <- object[, col_idx]

    cat(sprintf(
      "%s: min=%.3f, median=%.3f, mean=%.3f, max=%.3f\n",
      col_name,
      min(col_vals),
      median(col_vals),
      mean(col_vals),
      max(col_vals)
    ))
  }

  cat("\n")

  invisible(object)
}
