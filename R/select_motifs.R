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
#' @export
select_motifs <- function(motif_cluster, cluster_result, n, verbose = TRUE) {

  # ===== VALIDATION =====
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)) {
    stop("'n' must be a positive integer", call. = FALSE)
  }

  if (!is.data.frame(motif_cluster) || nrow(motif_cluster) == 0) {
    stop("'motif_cluster' must be a non-empty data.frame", call. = FALSE)
  }

  if (!is.list(cluster_result) || !"cluster_summary" %in% names(cluster_result)) {
    stop("'cluster_result' must contain 'cluster_summary' element", call. = FALSE)
  }

  cluster_summary <- cluster_result$cluster_summary
  required_cols <- c("cluster_id", "dominant_class")

  if (!all(required_cols %in% colnames(cluster_summary))) {
    stop("'cluster_summary' missing required columns", call. = FALSE)
  }

  # ===== SETUP =====
  n <- as.integer(n)
  dominant_classes <- unique(cluster_summary$dominant_class)
  k <- length(dominant_classes)
  m <- nrow(cluster_summary)
  classe_order <- names(sort(table(dominant_classes), decreasing = TRUE))

  if (verbose) {
    message(sprintf("Classes (k): %d | Clusters (m): %d | Motifs (n): %d", k, m, n))
  }

  # ===== CASE 1: n < k =====
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

  # ===== PREPARE CLUSTER MAPPING =====
  cluster_to_class <- setNames(
    cluster_summary$dominant_class,
    paste0("Cluster_", cluster_summary$cluster_id)
  )

  available_cols <- intersect(names(cluster_to_class), colnames(motif_cluster))
  if (length(available_cols) == 0) {
    stop("No matching cluster columns found in motif_cluster", call. = FALSE)
  }

  cluster_to_class <- cluster_to_class[available_cols]

  # ===== SELECT MOTIFS =====
  if (n < m) {
    if (verbose) message(sprintf("CASE 2: Distributing by class (%d <= %d < %d)", k, n, m))
    selected_motifs <- .select_by_class_fast(motif_cluster, cluster_to_class, classe_order, n, verbose)
    case_used <- "CASE_2"
  } else {
    if (verbose) message(sprintf("CASE 3: Distributing by cluster (%d >= %d)", n, m))
    selected_motifs <- .select_by_cluster_fast(motif_cluster, cluster_to_class, classe_order, n, m, verbose)
    case_used <- "CASE_3"
  }

  # ===== OUTPUT =====
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


#' @keywords internal
.select_by_class_fast <- function(motif_cluster, cluster_to_class, classe_order, n, verbose) {

  # Calculate distribution
  n_classes <- length(classe_order)
  motifs_per_classe <- rep(n %/% n_classes, n_classes)
  remainder <- n %% n_classes
  if (remainder > 0) motifs_per_classe[seq_len(remainder)] <- motifs_per_classe[seq_len(remainder)] + 1L
  names(motifs_per_classe) <- classe_order

  if (verbose) {
    message("Distribution: ", paste(sprintf("%s=%d", classe_order, motifs_per_classe), collapse=", "))
  }

  # Group columns by class
  class_cols <- split(names(cluster_to_class), cluster_to_class)
  selected_motifs <- setNames(vector("list", n_classes), classe_order)
  all_motifs <- rownames(motif_cluster)
  used_motifs <- setNames(rep(FALSE, length(all_motifs)), all_motifs)

  for (classe in classe_order) {
    cols <- class_cols[[classe]]
    if (length(cols) == 0) {
      selected_motifs[[classe]] <- character(0)
      next
    }

    scores <- rowSums(motif_cluster[, cols, drop = FALSE])
    top_idx <- order(scores, decreasing = TRUE)
    selected <- character(0)

    for (idx in top_idx) {
      if (length(selected) >= motifs_per_classe[classe]) break
      motif <- all_motifs[idx]
      if (!used_motifs[motif]) {
        selected <- c(selected, motif)
        used_motifs[motif] <- TRUE
      }
    }
    selected_motifs[[classe]] <- selected
  }

  selected_motifs
}


#' @keywords internal
.select_by_cluster_fast <- function(motif_cluster, cluster_to_class, classe_order, n, m, verbose) {

  motifs_per_cluster <- n %/% m
  remainder <- n %% m

  if (verbose) {
    message(sprintf("Base: %d per cluster | Remainder: %d", motifs_per_cluster, remainder))
  }

  classe_motifs <- setNames(vector("list", length(classe_order)), classe_order)
  all_motifs <- rownames(motif_cluster)
  used_motifs <- setNames(rep(FALSE, length(all_motifs)), all_motifs)

  # Process each cluster
  for (col in names(cluster_to_class)) {
    classe <- cluster_to_class[col]
    values <- motif_cluster[[col]]
    top_idx <- order(values, decreasing = TRUE)
    selected <- character(0)

    for (idx in top_idx) {
      if (length(selected) >= motifs_per_cluster) break
      motif <- all_motifs[idx]
      if (!used_motifs[motif]) {
        selected <- c(selected, motif)
        used_motifs[motif] <- TRUE
      }
    }
    classe_motifs[[classe]] <- c(classe_motifs[[classe]], selected)
  }

  # Distribute remainder
  if (remainder > 0) {
    for (i in seq_len(remainder)) {
      classe_to_add <- classe_order[((i - 1L) %% length(classe_order)) + 1L]
      cols <- names(cluster_to_class)[cluster_to_class == classe_to_add]
      if (length(cols) == 0) next

      scores <- rowSums(motif_cluster[, cols, drop = FALSE])
      top_idx <- order(scores, decreasing = TRUE)

      for (idx in top_idx) {
        motif <- all_motifs[idx]
        if (!used_motifs[motif]) {
          classe_motifs[[classe_to_add]] <- c(classe_motifs[[classe_to_add]], motif)
          used_motifs[motif] <- TRUE
          break
        }
      }
    }
  }

  classe_motifs
}


#' @export
print.select_motifs <- function(x, ...) {
  cat("\n=== Selected Motifs Summary ===\n\n")

  if (inherits(x, "empty_result")) {
    cat("Result: Empty (", attr(x, "reason"), ")\n", sep = "")
    return(invisible(x))
  }

  cat(sprintf("Total motifs (n):    %d\n", attr(x, "n")))
  cat(sprintf("Classes (k):         %d\n", attr(x, "k")))
  cat(sprintf("Clusters (m):        %d\n", attr(x, "m")))
  cat(sprintf("Case:                %s\n\n", attr(x, "case")))

  cat("Motifs per class:\n")
  for (classe in names(x)) {
    cat(sprintf("  %s: %d\n", classe, length(x[[classe]])))
  }
  cat("\n")

  invisible(x)
}


#' @export
summary.select_motifs <- function(object, ...) {
  cat("\n=== Selected Motifs Detailed Summary ===\n\n")

  if (inherits(object, "empty_result")) {
    cat("Result: Empty (", attr(object, "reason"), ")\n", sep = "")
    return(invisible(object))
  }

  cat(sprintf("Total motifs: %d | Classes: %d | Clusters: %d | Case: %s\n\n",
              attr(object, "n"), attr(object, "k"), attr(object, "m"), attr(object, "case")))

  cat("Selected motifs by class:\n")
  for (classe in names(object)) {
    motifs <- object[[classe]]
    cat(sprintf("\n%s (%d motifs):\n", classe, length(motifs)))
    if (length(motifs) > 0) {
      for (j in seq_along(motifs)) cat(sprintf("  [%d] %s\n", j, motifs[j]))
    } else {
      cat("  (none)\n")
    }
  }
  cat("\n")

  invisible(object)
}
