#' K-mer Analysis: Cluster and Class-Specific Motifs (Training + Validation)
#'
#' @param cluster_result Object from cluster_dendrogram()
#' @param motif_matrix Optional. Data frame from calculate_cluster_motifs()
#' @param threshold Numeric. Threshold for motif presence (default: 0)
#'
#' @return List with unique motifs, rankings, and frequency matrices for training and validation
#' @export
kmer_analysis <- function(cluster_result, motif_matrix = NULL, threshold = 0) {

  # ===== VALIDATION =====
  if (!is.list(cluster_result) || !"cluster_summary" %in% names(cluster_result)) {
    stop("'cluster_result' must contain 'cluster_summary' element", call. = FALSE)
  }

  cluster_summary <- cluster_result$cluster_summary
  if (!is.data.frame(cluster_summary) ||
      !all(c("cluster_id", "dominant_class") %in% colnames(cluster_summary))) {
    stop("'cluster_summary' must be a data.frame with 'cluster_id' and 'dominant_class'",
         call. = FALSE)
  }

  if (!is.numeric(threshold) || length(threshold) != 1) {
    stop("'threshold' must be a single numeric value", call. = FALSE)
  }

  # ===== SETUP =====
  if (is.null(motif_matrix)) motif_matrix <- calculate_cluster_motifs(cluster_result)

  motif_mat <- as.matrix(motif_matrix)
  motif_names <- rownames(motif_mat)
  cluster_names <- colnames(motif_mat)

  cluster_to_class <- setNames(
    cluster_summary$dominant_class,
    paste0("Cluster_", cluster_summary$cluster_id)
  )

  motif_presence <- motif_mat > threshold
  clusters_per_motif <- rowSums(motif_presence)

  # ===== UNIQUE CLUSTER MOTIFS =====
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

  # ===== CLUSTER FREQUENCY RANKING =====
  cluster_frequency_ranking <- data.frame(
    motif = motif_names,
    motif_mat,
    n_clusters = clusters_per_motif,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  cluster_frequency_ranking <- cluster_frequency_ranking[order(-cluster_frequency_ranking$n_clusters), ]
  rownames(cluster_frequency_ranking) <- NULL

  # ===== CLASS FREQUENCY MATRIX =====
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

  # ===== UNIQUE CLASS MOTIFS =====
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

  # ===== HELPER FUNCTION =====
  create_motifs_rank <- function(data, class_col = "CLASS") {
    class_col_idx <- which(colnames(data) == class_col)
    motif_cols <- setdiff(seq_len(ncol(data)), class_col_idx)

    classes <- data[[class_col]]
    motifs <- data[, motif_cols, drop = FALSE]
    motif_names <- colnames(motifs)
    unique_classes <- unique(classes)

    sums_mat <- matrix(0, nrow = length(motif_names), ncol = length(unique_classes),
                       dimnames = list(motif_names, unique_classes))
    counts_mat <- matrix(0L, nrow = length(motif_names), ncol = length(unique_classes),
                         dimnames = list(motif_names, unique_classes))

    for (class_name in unique_classes) {
      class_rows <- classes == class_name
      class_data <- motifs[class_rows, , drop = FALSE]
      sums_mat[, class_name] <- colSums(class_data, na.rm = TRUE)
      counts_mat[, class_name] <- colSums(class_data > 0, na.rm = TRUE)
    }

    value_width <- nchar(as.character(ceiling(max(sums_mat))))
    formatted_mat <- matrix("", nrow = nrow(sums_mat), ncol = ncol(sums_mat),
                            dimnames = dimnames(sums_mat))

    for (i in seq_len(nrow(sums_mat))) {
      for (j in seq_len(ncol(sums_mat))) {
        formatted_mat[i, j] <- sprintf("%0*d|%d", value_width,
                                       round(sums_mat[i, j]), counts_mat[i, j])
      }
    }

    result <- as.data.frame(formatted_mat, stringsAsFactors = FALSE)
    result$motif <- motif_names
    result$mean <- round(rowMeans(sums_mat, na.rm = TRUE), 2)
    result <- result[, c("motif", unique_classes, "mean")]
    rownames(result) <- NULL
    return(result)
  }

  # ===== TRAINING MOTIFS BY CLASS RANK =====
  motifs_by_class_rank <- data.frame(motif = character(0))
  if (!is.null(cluster_result$data_result) &&
      !is.null(cluster_result$data_result$kmers) &&
      "CLASS" %in% colnames(cluster_result$data_result$kmers)) {
    motifs_by_class_rank <- create_motifs_rank(cluster_result$data_result$kmers)
  }

  # ===== VALIDATION ANALYSIS =====
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
      motifs_by_class_rank = create_motifs_rank(validation_data),
      n_sequences = nrow(validation_data),
      n_classes = length(val_unique_classes),
      class_names = val_unique_classes
    )
  }

  # ===== RETURN =====
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


#' @export
print.kmer_analysis_result <- function(x, ...) {
  cat("\nK-mer Analysis Results\n")
  cat("======================\n\n")
  cat("TRAINING DATA:\n")
  cat(sprintf("  Cluster-specific motifs:  %d\n", nrow(x$unique_cluster_motifs)))
  cat(sprintf("  Class-specific motifs:    %d\n", nrow(x$unique_class_motifs)))
  cat(sprintf("  Total motifs analyzed:    %d\n", nrow(x$cluster_frequency_ranking)))
  cat(sprintf("  Total clusters:           %d\n", length(x$cluster_to_class)))
  cat(sprintf("  Total classes:            %d\n", ncol(x$class_frequency_matrix) - 1))

  if (x$has_external_validation) {
    cat("\nVALIDATION DATA:\n")
    cat(sprintf("  Class-specific motifs:    %d\n", nrow(x$validation$unique_class_motifs)))
    cat(sprintf("  Total sequences:          %d\n", x$validation$n_sequences))
    cat(sprintf("  Total classes:            %d\n", x$validation$n_classes))
  }
  cat("\n")
  invisible(x)
}


#' @export
summary.kmer_analysis_result <- function(object, ...) {
  print(object)
  cat("\n=== TRAINING SUMMARY ===\n")

  if (nrow(object$unique_cluster_motifs) > 0) {
    cat("\nTop 10 cluster-specific motifs:\n")
    print(head(object$unique_cluster_motifs[order(-object$unique_cluster_motifs$value), ], 10))
  }

  if (nrow(object$unique_class_motifs) > 0) {
    cat("\nTop 10 class-specific motifs:\n")
    print(head(object$unique_class_motifs, 10))
  }

  cat("\nTop 10 most frequent motifs:\n")
  print(head(object$cluster_frequency_ranking, 10))

  if (!is.null(object$validation)) {
    cat("\n=== VALIDATION SUMMARY ===\n")
    if (nrow(object$validation$unique_class_motifs) > 0) {
      cat("\nTop 10 class-specific motifs:\n")
      print(head(object$validation$unique_class_motifs[order(-object$validation$unique_class_motifs$value), ], 10))
    }
    cat("\nTop 10 most frequent motifs:\n")
    print(head(object$validation$class_frequency_ranking, 10))
  }

  invisible(object)
}
