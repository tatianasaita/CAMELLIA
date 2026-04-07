#' K-mer Analysis: Cluster and Class-Specific Motifs (Training + Test)
#'
#' @param cluster_result Object from cluster_dendrogram()
#' @param motif_matrix Data frame from calculate_cluster_motifs()
#'
#' @return List with unique motifs, rankings, and frequency matrices for
#'   training and test sets.
#'
#' @details
#' Performs comprehensive k-mer analysis across two levels:
#' \itemize{
#'   \item Cluster-level: Identifies motifs unique to specific clusters and
#'     ranks all motifs by cluster frequency
#'   \item Class-level: Aggregates cluster motifs to class level, identifying
#'     class-specific biomarkers and computing frequency matrices
#' }
#'
#' If external validation data is available in \code{cluster_result}, the same
#' analysis is performed independently on the validation set, enabling comparison
#' of motif patterns between training and validation datasets. Motifs with values
#' above \code{threshold} are considered present.
#'
#' @note
#' \itemize{
#'   \item Requires .create_motifs_rank. See \code{internal-functions.R} for details.
#'   \item S3 methods available. See \code{methods.R} for details.
#' }
#'
#' @examples
#' \dontrun{
#' result_kmer_analysis <- kmer_analysis(
#'   cluster_result = result_dendrogram,
#'   motif_matrix   = result_cluster_motifs
#' )
#' }
#'
#' @export
kmer_analysis <- function(cluster_result, motif_matrix) {

# Setup
  motif_mat <- as.matrix(motif_matrix)

  if (is.null(rownames(motif_mat)) || is.null(colnames(motif_mat))) {
    stop("'motif_matrix' must have both rownames (motifs) and colnames (clusters).")
  }

  motif_names   <- rownames(motif_mat)
  cluster_names <- colnames(motif_mat)
  cluster_summary <- cluster_result$cluster_summary
  if (is.null(cluster_summary) || nrow(cluster_summary) == 0L) {
    stop("'cluster_result$cluster_summary' is NULL or empty.")
  }

  cluster_to_class_full <- setNames(
    cluster_summary$dominant_class,
    paste0("Cluster_", cluster_summary$cluster_id)
  )

  valid_clusters   <- intersect(names(cluster_to_class_full), cluster_names)
  cluster_to_class <- cluster_to_class_full[valid_clusters]

  motif_presence     <- matrix(motif_mat > 0,
                               nrow = nrow(motif_mat),
                               ncol = ncol(motif_mat),
                               dimnames = dimnames(motif_mat))
  clusters_per_motif <- rowSums(motif_presence)

# Unique clusters motifs
  unique_idx <- which(clusters_per_motif == 1L)
  if (length(unique_idx) > 0L) {
    cluster_idx <- unname(sapply(
      seq_along(unique_idx),
      function(i) which(motif_presence[unique_idx[i], , drop = FALSE])[1L]
    ))
    unique_cluster_motifs <- data.frame(
      motif   = motif_names[unique_idx],
      cluster = cluster_names[cluster_idx],
      class   = cluster_to_class[cluster_names[cluster_idx]],
      value   = motif_mat[cbind(unique_idx, cluster_idx)],
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  } else {
    unique_cluster_motifs <- data.frame(
      motif   = character(0),
      cluster = character(0),
      class   = character(0),
      value   = numeric(0)
    )
  }

  # Cluster frequency ranking
  cluster_frequency_ranking <- data.frame(
    motif      = motif_names,
    motif_mat,
    n_clusters = clusters_per_motif,
    stringsAsFactors = FALSE,
    row.names  = NULL
  )
  cluster_frequency_ranking <- cluster_frequency_ranking[
    order(-cluster_frequency_ranking$n_clusters), ]
  rownames(cluster_frequency_ranking) <- NULL

  # Class frequency matrix
  all_classes    <- unique(cluster_to_class)
  class_freq_mat <- matrix(
    0L,
    nrow     = nrow(motif_mat),
    ncol     = length(all_classes),
    dimnames = list(motif_names, all_classes)
  )

  for (class_name in all_classes) {
    class_cluster_names <- names(cluster_to_class)[cluster_to_class == class_name]
    class_cluster_names <- intersect(class_cluster_names, colnames(motif_presence))

    class_freq_mat[, class_name] <- if (length(class_cluster_names) == 0L) {
      0L
    } else {
      rowSums(motif_presence[, class_cluster_names, drop = FALSE])
    }
  }

  class_frequency_matrix       <- as.data.frame(class_freq_mat)
  class_frequency_matrix$Total <- rowSums(class_freq_mat)
  class_frequency_matrix       <- class_frequency_matrix[
    order(-class_frequency_matrix$Total), ]

  # Unique class motifs
  classes_per_motif <- rowSums(class_freq_mat > 0L)
  unique_class_idx  <- which(classes_per_motif == 1L)

  if (length(unique_class_idx) > 0L) {
    class_idx <- unname(sapply(
      seq_along(unique_class_idx),
      function(i) which(class_freq_mat[unique_class_idx[i], ] > 0L)[1L]
    ))
    class_names_vec <- all_classes[class_idx]
    unique_class_motifs <- data.frame(
      motif                 = motif_names[unique_class_idx],
      class                 = class_names_vec,
      n_clusters_with_motif = class_freq_mat[cbind(unique_class_idx, class_idx)],
      total_class_clusters  = as.integer(table(cluster_to_class)[class_names_vec]),
      stringsAsFactors      = FALSE,
      row.names             = NULL
    )
    unique_class_motifs <- unique_class_motifs[
      order(unique_class_motifs$class, -unique_class_motifs$n_clusters_with_motif), ]
    rownames(unique_class_motifs) <- NULL
  } else {
    unique_class_motifs <- data.frame(
      motif                 = character(0),
      class                 = character(0),
      n_clusters_with_motif = integer(0),
      total_class_clusters  = integer(0)
    )
  }

  # Training motifs by class rank
  motifs_by_class_rank <- data.frame(motif = character(0))
  if (!is.null(cluster_result$data_result)       &&
      !is.null(cluster_result$data_result$kmers) &&
      "CLASS" %in% colnames(cluster_result$data_result$kmers)) {
    motifs_by_class_rank <- .create_motifs_rank(cluster_result$data_result$kmers)
  }

  # Test set analysis (external/held-out)
  has_external_test <- !is.null(cluster_result$test_data)  &&
    is.data.frame(cluster_result$test_data)                 &&
    nrow(cluster_result$test_data) > 0L                     &&
    "CLASS" %in% colnames(cluster_result$test_data)

  test_results <- NULL

  if (has_external_test) {
    test_data     <- cluster_result$test_data
    class_col_idx <- which(colnames(test_data) == "CLASS")
    motif_cols    <- setdiff(seq_len(ncol(test_data)), class_col_idx)

    tst_motifs         <- test_data[, motif_cols, drop = FALSE]
    tst_classes        <- test_data$CLASS
    tst_unique_classes <- unique(tst_classes)
    tst_motif_names    <- colnames(tst_motifs)

    tst_class_mat <- matrix(
      0,
      nrow     = length(tst_motif_names),
      ncol     = length(tst_unique_classes),
      dimnames = list(tst_motif_names, tst_unique_classes)
    )

    for (class_name in tst_unique_classes) {
      class_rows <- tst_classes == class_name
      tst_class_mat[, class_name] <- colSums(
        tst_motifs[class_rows, , drop = FALSE], na.rm = TRUE
      )
    }

    tst_presence        <- tst_class_mat > 0
    tst_classes_per_mot <- rowSums(tst_presence)
    tst_unique_idx      <- which(tst_classes_per_mot == 1L)

    if (length(tst_unique_idx) > 0L) {
      tst_class_idx <- unname(sapply(
        seq_along(tst_unique_idx),
        function(i) which(tst_presence[tst_unique_idx[i], ])[1L]
      ))
      tst_unique_class_motifs <- data.frame(
        motif = tst_motif_names[tst_unique_idx],
        class = tst_unique_classes[tst_class_idx],
        value = tst_class_mat[cbind(tst_unique_idx, tst_class_idx)],
        stringsAsFactors = FALSE,
        row.names = NULL
      )
    } else {
      tst_unique_class_motifs <- data.frame(
        motif = character(0),
        class = character(0),
        value = numeric(0)
      )
    }

    # Frequency ranking — safe assembly to avoid column name collisions
    tst_cfr           <- as.data.frame(tst_class_mat, stringsAsFactors = FALSE)
    tst_cfr$motif     <- tst_motif_names
    tst_cfr$n_classes <- tst_classes_per_mot
    tst_cfr           <- tst_cfr[, c("motif", tst_unique_classes, "n_classes")]
    tst_class_frequency_ranking           <- tst_cfr[order(-tst_cfr$n_classes), ]
    rownames(tst_class_frequency_ranking) <- NULL

    # Frequency matrix
    tst_class_frequency_matrix        <- as.data.frame(tst_class_mat)
    tst_class_frequency_matrix$Total  <- rowSums(tst_class_mat)
    tst_class_frequency_matrix        <- tst_class_frequency_matrix[
      order(-tst_class_frequency_matrix$Total), ]

    test_results <- list(
      unique_class_motifs     = tst_unique_class_motifs,
      class_frequency_ranking = tst_class_frequency_ranking,
      class_frequency_matrix  = tst_class_frequency_matrix,
      motifs_by_class_rank    = .create_motifs_rank(test_data),
      n_sequences             = nrow(test_data),
      n_classes               = length(tst_unique_classes),
      class_names             = tst_unique_classes
    )
  }

  # Return results
  structure(
    list(
      unique_cluster_motifs     = unique_cluster_motifs,
      cluster_frequency_ranking = cluster_frequency_ranking,
      unique_class_motifs       = unique_class_motifs,
      class_frequency_matrix    = class_frequency_matrix,
      motifs_by_class_rank      = motifs_by_class_rank,
      cluster_to_class          = cluster_to_class,
      has_external_test         = has_external_test,
      test                      = test_results
    ),
    class = "kmer_analysis_result"
  )
}

