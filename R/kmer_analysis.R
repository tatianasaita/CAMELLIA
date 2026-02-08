#' K-mer Analysis: Cluster and Class-Specific Motifs (Training + Validation)
#'
#' @param cluster_result Object from cluster_dendrogram()
#' @param motif_matrix Data frame from calculate_cluster_motifs()
#' @param threshold Numeric. Threshold for motif presence (default: 0)
#'
#' @return List with unique motifs, rankings, and frequency matrices for training and validation
#'
#' @deatils
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
#'   cluster_result = dendrogram_result,
#'   motif_matrix = motif_freq_matrix,
#'   threshold = 5
#' )
#'}
#'
#' @export
kmer_analysis <- function(cluster_result, motif_matrix, threshold = 0) {

# Setup
  motif_mat <- as.matrix(motif_matrix)
  motif_names <- rownames(motif_mat)
  cluster_names <- colnames(motif_mat)
  cluster_summary <- cluster_result$cluster_summary

  cluster_to_class <- setNames(
    cluster_summary$dominant_class,
    paste0("Cluster_", cluster_summary$cluster_id)
  )

  motif_presence <- motif_mat > threshold
  clusters_per_motif <- rowSums(motif_presence)

# Unique clusters motifs
  unique_idx <- which(clusters_per_motif == 1)
  if (length(unique_idx) > 0) {
    cluster_idx <- apply(motif_presence[unique_idx, , drop = FALSE], 1, which)
    unique_cluster_motifs <- data.frame(
      motif = motif_names[unique_idx],
      cluster = cluster_names[cluster_idx],
      class = cluster_to_class[cluster_names[cluster_idx]],
      value = motif_mat[cbind(unique_idx, cluster_idx)],
      stringsAsFactors = FALSE
    )
  } else {
    unique_cluster_motifs <- data.frame(
      motif = character(0), cluster = character(0),
      class = character(0), value = numeric(0)
    )
  }

  # Cluster frequency ranking
  cluster_frequency_ranking <- data.frame(
    motif = motif_names,
    motif_mat,
    n_clusters = clusters_per_motif,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  cluster_frequency_ranking <- cluster_frequency_ranking[order(-cluster_frequency_ranking$n_clusters), ]
  rownames(cluster_frequency_ranking) <- NULL

  # Class frequency matrix
  all_classes <- unique(cluster_to_class)
  class_freq_mat <- matrix(
    0L,
    nrow = nrow(motif_mat),
    ncol = length(all_classes),
    dimnames = list(motif_names, all_classes)
  )

  for (class_name in all_classes) {
    class_clusters <- cluster_to_class == class_name
    class_freq_mat[, class_name] <- rowSums(motif_presence[, class_clusters, drop = FALSE])
  }

  class_frequency_matrix <- as.data.frame(class_freq_mat)
  class_frequency_matrix$Total <- rowSums(class_freq_mat)
  class_frequency_matrix <- class_frequency_matrix[order(-class_frequency_matrix$Total), ]

  # Unique class motifs
  classes_per_motif <- rowSums(class_freq_mat > 0)
  unique_class_idx <- which(classes_per_motif == 1)

  if (length(unique_class_idx) > 0) {
    class_idx <- apply(class_freq_mat[unique_class_idx, , drop = FALSE] > 0, 1, which)
    class_names <- all_classes[class_idx]
    unique_class_motifs <- data.frame(
      motif = motif_names[unique_class_idx],
      class = class_names,
      n_clusters_with_motif = class_freq_mat[cbind(unique_class_idx, class_idx)],
      total_class_clusters = as.integer(table(cluster_to_class)[class_names]),
      stringsAsFactors = FALSE
    )
    unique_class_motifs <- unique_class_motifs[order(unique_class_motifs$class,
                                                     -unique_class_motifs$n_clusters_with_motif), ]
    rownames(unique_class_motifs) <- NULL
  } else {
    unique_class_motifs <- data.frame(
      motif = character(0), class = character(0),
      n_clusters_with_motif = integer(0), total_class_clusters = integer(0)
    )
  }

  # Training motifs by class rank
  motifs_by_class_rank <- data.frame(motif = character(0))
  if (!is.null(cluster_result$data_result) &&
      !is.null(cluster_result$data_result$kmers) &&
      "CLASS" %in% colnames(cluster_result$data_result$kmers)) {
    motifs_by_class_rank <- .create_motifs_rank(cluster_result$data_result$kmers)
  }

  # Validation analysis
  has_external_validation <- !is.null(cluster_result$classification_result) &&
    !is.null(cluster_result$classification_result$validation_data) &&
    "CLASS" %in% colnames(cluster_result$classification_result$validation_data)

  validation_results <- NULL
  if (has_external_validation) {
    validation_data <- cluster_result$classification_result$validation_data
    class_col_idx <- which(colnames(validation_data) == "CLASS")
    motif_cols <- setdiff(seq_len(ncol(validation_data)), class_col_idx)

    val_motifs <- validation_data[, motif_cols, drop = FALSE]
    val_classes <- validation_data$CLASS
    val_unique_classes <- unique(val_classes)
    val_motif_names <- colnames(val_motifs)

    val_class_mat <- matrix(0, nrow = length(val_motif_names), ncol = length(val_unique_classes),
                            dimnames = list(val_motif_names, val_unique_classes))

    for (class_name in val_unique_classes) {
      class_rows <- val_classes == class_name
      val_class_mat[, class_name] <- colSums(val_motifs[class_rows, , drop = FALSE], na.rm = TRUE)
    }

    val_presence <- val_class_mat > threshold
    val_classes_per_motif <- rowSums(val_presence)
    val_unique_idx <- which(val_classes_per_motif == 1)

    if (length(val_unique_idx) > 0) {
      val_class_idx <- apply(val_presence[val_unique_idx, , drop = FALSE], 1, which)
      val_unique_class_motifs <- data.frame(
        motif = val_motif_names[val_unique_idx],
        class = val_unique_classes[val_class_idx],
        value = val_class_mat[cbind(val_unique_idx, val_class_idx)],
        stringsAsFactors = FALSE
      )
    } else {
      val_unique_class_motifs <- data.frame(
        motif = character(0), class = character(0), value = numeric(0)
      )
    }

    val_class_frequency_ranking <- data.frame(
      motif = val_motif_names,
      val_class_mat,
      n_classes = val_classes_per_motif,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    val_class_frequency_ranking <- val_class_frequency_ranking[order(-val_class_frequency_ranking$n_classes), ]
    rownames(val_class_frequency_ranking) <- NULL

    val_class_frequency_matrix <- as.data.frame(val_class_mat)
    val_class_frequency_matrix$Total <- rowSums(val_class_mat)
    val_class_frequency_matrix <- val_class_frequency_matrix[order(-val_class_frequency_matrix$Total), ]

    validation_results <- list(
      unique_class_motifs = val_unique_class_motifs,
      class_frequency_ranking = val_class_frequency_ranking,
      class_frequency_matrix = val_class_frequency_matrix,
      motifs_by_class_rank = .create_motifs_rank(validation_data),
      n_sequences = nrow(validation_data),
      n_classes = length(val_unique_classes),
      class_names = val_unique_classes
    )
  }

  # Return results
  structure(
    list(
      unique_cluster_motifs = unique_cluster_motifs,
      cluster_frequency_ranking = cluster_frequency_ranking,
      unique_class_motifs = unique_class_motifs,
      class_frequency_matrix = class_frequency_matrix,
      motifs_by_class_rank = motifs_by_class_rank,
      cluster_to_class = cluster_to_class,
      has_external_validation = has_external_validation,
      validation = validation_results
    ),
    class = "kmer_analysis_result"
  )
}

