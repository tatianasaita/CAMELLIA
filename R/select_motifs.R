#' Select Motifs from Cluster Analysis
#'
#' Selects n motifs based on cluster information, distributing them intelligently
#' across classes according to their distribution across all clusters.
#'
#' @param motif_cluster Data.frame with k-mer motifs by cluster (normalized).
#' @param cluster_result List containing cluster_summary from cluster_dendrogram().
#' @param n Integer. Number of motifs to select.
#' @param verbose Logical. If TRUE, prints detailed information. Default is TRUE.
#'
#' @return An object of class "select_motifs" with selected motifs per class.
#'
#' @details
#' The function applies one of three selection strategies:
#' \itemize{
#'   \item Case 1 (n < k): Returns empty result - insufficient motifs
#'   \item Case 2 (k ≤ n < m): Distributes evenly by class
#'   \item Case 3 (n ≥ m): Distributes evenly by cluster
#'}
#' @note
#' \itemize{
#'   \item Requires .select_by_cluster_fast, .select_by_class_fast. See \code{internal-functions.R} for details.
#'   \item S3 methods available. See \code{methods.R} for details.
#' }
#'
#' @examples
#' \dontrun{
#' result_selected_motifs <- select_motifs(
#'   motif_cluster = result_cluster_motifs,
#'   cluster_result = result_cluster_dendrogram,
#'   n = 20
#' )
#' }
#'
#' @export
select_motifs <- function(motif_cluster, cluster_result, n, verbose = TRUE) {

  n <- as.integer(n)
  cluster_summary <- cluster_result$cluster_summary
  dominant_classes <- unique(cluster_result$cluster_summary$dominant_class)
  k <- length(dominant_classes)
  m <- nrow(cluster_summary)
  classe_order <- names(sort(table(dominant_classes), decreasing = TRUE))

  if (verbose) {
    message(sprintf("Classes (k): %d | Clusters (m): %d | Motifs (n): %d", k, m, n))
  }

# Case 1: n < k
  if (n < k) {
    if (verbose) message(sprintf("CASE 1: n < k (%d < %d) - Insufficient motifs", n, k))
    return(structure(
      list(),
      class = c("select_motifs", "empty_result"),
      n = n, k = k, m = m,
      case = "CASE_1",
      reason = "n < k: insufficient motifs for all classes"
    ))
  }

# Select motifs: Case 2 and 3
  cluster_to_class <- setNames(
    cluster_summary$dominant_class,
    paste0("Cluster_", cluster_summary$cluster_id)
  )

  if (n < m) {
    if (verbose) message(sprintf("CASE 2: Distributing by class (%d <= %d < %d)", k, n, m))
    selected_motifs <- .select_by_class_fast(motif_cluster, cluster_to_class, classe_order, n)
    case_used <- "CASE_2"
  } else {
    if (verbose) message(sprintf("CASE 3: Distributing by cluster (%d >= %d)", n, m))
    selected_motifs <- .select_by_cluster_fast(motif_cluster, cluster_to_class, classe_order, n, m)
    case_used <- "CASE_3"
  }

# Output
  if (verbose) {
    message("\nSelected motifs by class:")
    for (classe in classe_order) {
      if (classe %in% names(selected_motifs)) {
        message(sprintf("  %s: %d motifs", classe, length(selected_motifs[[classe]])))
      }
    }
  }

  structure(
    selected_motifs,
    class = c("select_motifs", "list"),
    n = n, k = k, m = m,
    case = case_used,
    classes_order = classe_order,
    motifs_per_class = lengths(selected_motifs)
  )
}

