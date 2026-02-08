#' Calculate Cluster Motifs from K-mer Counts
#'
#' Calculates motif profiles for each cluster by summing k-mer counts of all
#' sequences in the cluster. The resulting motif matrix is normalized using
#' min-max normalization applied independently to each cluster column.
#'
#' @param result_cluster_dendrogram,
#'   returned by cluster_dendrogram().
#' @param ... Additional arguments passed to methods.
#'
#' @return A data.frame of class 'cluster_motifs' with k-mer motifs by cluster (normalized).
#'   Rows are k-mers, columns are clusters (named 'Cluster_1', 'Cluster_2', etc.).
#'   Values are normalized between 0 and 1 independently for each cluster column.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Gets k-mer counts and cluster assignments
#'   \item For each cluster, sums k-mer values across
#'     all member sequences
#'   \item Creates a matrix with k-mers as rows
#'     and clusters as columns
#'   \item Applies min-max normalization
#' }
#'
#' @note
#' \itemize{
#'   \item Requires .normalize_motif_matrix. See \code{internal-functions.R} for details.
#'   \item S3 methods available. See \code{methods.R} for details.
#' }
#'
#' @examples
#' \dontrun{
#' result_cluster_motifs <- calculate_cluster_motifs(data_result)
#' }
#'
#' @export
calculate_cluster_motifs <- function(result_cluster_dendrogram, ...) {
  cluster_assignments <- result_cluster_dendrogram$cluster_assignment_dendro_order
  kmers_df <- result_cluster_dendrogram$data_result$kmers

  # Find and exclude CLASS column
  class_col_idx <- which(tolower(colnames(kmers_df)) == "class")

  if (length(class_col_idx) > 0) {
    kmers_numeric <- kmers_df[, -class_col_idx, drop = FALSE]
  } else {
    kmers_numeric <- kmers_df
  }

  # Get numeric column names (k-mers)
  numeric_cols <- colnames(kmers_numeric)

  # Get unique clusters and cluster count
  cluster_ids <- sort(unique(cluster_assignments))
  n_clusters <- length(cluster_ids)

  motif_raw <- matrix(
    nrow = length(numeric_cols),
    ncol = n_clusters,
    dimnames = list(numeric_cols, paste0("Cluster_", cluster_ids))
  )

  for (i in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[i]
    cluster_mask <- cluster_assignments == cluster_id

    cluster_kmers <- kmers_numeric[cluster_mask, , drop = FALSE]
    motif_raw[, i] <- colSums(cluster_kmers)
  }

  # Min-max normalization per column
  motif_normalized <- .normalize_motif_matrix(motif_raw)

  # Output data.frame
  motif_cluster <- as.data.frame(motif_normalized)

  # Return result
  class(motif_cluster) <- c("cluster_motifs", "data.frame")
  return(invisible(motif_cluster))
}

