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

  method <- result_cluster_dendrogram$method

  kmers_df <- result_cluster_dendrogram$data_result$kmers

  if (is.null(kmers_df)) {
    stop("data_result$kmers is NULL: verifique o objeto retornado por cluster_dendrogram()")
  }

  # Exclude CLASS column if present
  class_col_idx <- which(tolower(colnames(kmers_df)) == "class")
  kmers_numeric <- if (length(class_col_idx) > 0) {
    kmers_df[, -class_col_idx, drop = FALSE]
  } else {
    kmers_df
  }

  if (method == "apcluster") {

    # ── AJUSTE: subsetar apenas as sequências que participaram da clusterização
    selected_indices <- result_cluster_dendrogram$selected_indices

    if (is.null(selected_indices)) {
      stop("selected_indices is NULL: verifique o objeto retornado por cluster_dendrogram()")
    }

    kmers_numeric      <- kmers_numeric[selected_indices, , drop = FALSE]
    cluster_assignments <- result_cluster_dendrogram$cluster_assignment[selected_indices]

    if (nrow(kmers_numeric) != length(cluster_assignments)) {
      stop(sprintf(
        "Mismatch: kmers subset has %d rows but cluster_assignment subset has length %d.",
        nrow(kmers_numeric), length(cluster_assignments)
      ))
    }

  } else {

    cluster_assignments <- if (!is.null(result_cluster_dendrogram$cluster_assignment_dendro_order)) {
      result_cluster_dendrogram$cluster_assignment_dendro_order
    } else {
      result_cluster_dendrogram$cluster_assignment
    }

    if (nrow(kmers_numeric) != length(cluster_assignments)) {
      stop(sprintf(
        "Mismatch: kmers has %d rows but cluster_assignment has length %d.",
        nrow(kmers_numeric), length(cluster_assignments)
      ))
    }
  }

  if (is.null(cluster_assignments)) {
    stop("cluster_assignments is NULL: verifique o objeto retornado por cluster_dendrogram()")
  }

  # Get unique clusters and cluster count
  cluster_ids <- sort(unique(cluster_assignments[!is.na(cluster_assignments)]))
  n_clusters  <- length(cluster_ids)

  if (n_clusters == 0L) {
    stop("No valid clusters found (all assignments are NA).")
  }

  numeric_cols <- colnames(kmers_numeric)

  motif_raw <- matrix(
    nrow     = length(numeric_cols),
    ncol     = n_clusters,
    dimnames = list(numeric_cols, paste0("Cluster_", cluster_ids))
  )

  for (i in seq_along(cluster_ids)) {
    cluster_id <- cluster_ids[i]

    # NA-safe mask: NA == cluster_id retorna NA, tratado como FALSE
    cluster_mask  <- !is.na(cluster_assignments) & cluster_assignments == cluster_id
    cluster_kmers <- kmers_numeric[cluster_mask, , drop = FALSE]

    if (nrow(cluster_kmers) == 0L) {
      warning(sprintf(
        "Cluster %d has no sequences in kmers_numeric — column will be NA.",
        cluster_id
      ))
      motif_raw[, i] <- NA_real_
      next
    }

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

