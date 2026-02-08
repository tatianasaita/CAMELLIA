#' S3 Methods of create_data.R
#'
#' Print Method for kmer_data Objects
#'
#' @param x A kmer_data object created by \code{\link{create_data}}.
#' @param ... Additional arguments (not currently used).
#'
#' @return The input object \code{x}, invisibly.
#'
#' @export
#' @method print kmer_data
print.kmer_data <- function(x, ...) {

  cat("K-mer Data Object\n")
  cat("=================\n\n")

  # K-mer matrix summary
  cat("K-mer Matrix:\n")
  cat("  Dimensions:", nrow(x$kmers), "sequences x",
      ncol(x$kmers) - 1L, "k-mers + 1 CLASS column\n")

  kmer_cols <- which(colnames(x$kmers) != "CLASS")
  min_count <- min(x$kmers[, kmer_cols])
  max_count <- max(x$kmers[, kmer_cols])

  cat("  Range of counts: [", min_count, ", ", max_count, "]\n", sep = "")

  # Metadata summary
  cat("\nMetadata:\n")
  cat("  Sequences:", nrow(x$metadata), "\n")
  cat("  Mean sequence length:", round(mean(x$metadata$length), 2), "bp\n")
  cat("  Min sequence length:", min(x$metadata$length), "bp\n")
  cat("  Max sequence length:", max(x$metadata$length), "bp\n")

  # Class distribution
  cat("\nClass Distribution:\n")
  class_dist <- table(x$metadata$class)

  for (i in seq_along(class_dist)) {
    cat("  ", names(class_dist)[i], ": ", class_dist[i], " sequences\n", sep = "")
  }

  invisible(x)
}

################################################################################
#' S3 Method of create_dendrogram.R
#'
#' Displays a summary of the dendrogram analysis result.
#'
#' @param x A dendrogram_result object created by \code{\link{create_dendrogram}}.
#' @param ... Additional arguments (not currently used).
#'
#' @return The input object \code{x}, invisibly.
#'
#' @export
#' @method print dendrogram_result
print.dendrogram_result <- function(x, ...) {
  cat("Dendrogram Analysis Result\n")
  cat("==========================\n\n")
  cat("Components:\n")
  cat("  - dendrogram: hierarchical clustering dendrogram (with colored leaves)\n")
  cat("  - hclust: hclust object\n")
  cat("  - order: dendrogram leaf order\n")
  cat("  - labels: ordered class labels\n")
  if (!is.null(x$sequence_names)) {
    cat("  - sequence_names: ordered sequence names\n")
  }
  cat("  - colors: color mapping for each sample\n")
  cat("  - base_colors: color mapping for each class\n")
  cat("\n")

  # Calculate derived information from existing components
  n_samples <- length(x$order)
  n_classes <- length(x$base_colors)
  dend_height <- attr(x$dendrogram, "height")

  cat("Number of samples:", n_samples, "\n")
  cat("Number of classes:", n_classes, "\n")
  cat("Tree height:", round(dend_height, 4), "\n\n")

  cat("Color mapping by class:\n")
  for (i in seq_along(x$base_colors)) {
    class_name <- names(x$base_colors)[i]
    class_color <- x$base_colors[i]
    cat(sprintf("  %s: %s\n", class_name, class_color))
  }

  invisible(x)
}

################################################################################
#' S3 Method of cluster_dendrogram.R
#'
#'  Print Method for Dendrogram Clustering Result
#'
#' @param x An object of class \code{cluster_dendrogram_result} returned by the clustering function.
#' @param ... Additional arguments passed to print methods.
#'
#' @return Invisibly returns \code{x}. The function is called
#'   primarily for its side effect of printing a formatted summary to the console.
#'
#' @export
print.cluster_dendrogram_result <- function(x, ...) {
  cat("\nDendrogram Clustering Result\n")
  cat("=============================\n\n")
  cat(sprintf("Total elements:     %d\n", x$n_elements))
  cat(sprintf("Total clusters:     %d\n", nrow(x$cluster_summary)))
  cat(sprintf("Unassigned:         %d\n", x$n_unassigned))
  cat(sprintf("Min size:           %d\n", x$min_size))
  cat(sprintf("Hom threshold:      %.3f\n", x$hom_thresh))

  if (nrow(x$cluster_summary) > 0) {
    cat("\nCluster Summary:\n")
    print(x$cluster_summary, row.names = FALSE)
  }

  if (x$n_unassigned > 0) {
    cat("\n[WARNING] Some elements could not be assigned to clusters\n")
  }

  invisible(x)
}

################################################################################
#' S3 Method of calculate_cluster_motifs.R
#'
#' Print Method for Cluster Motifs
#'
#' @param x Object of class 'cluster_motifs' returned by calculate_cluster_motifs().
#' @param ... Additional arguments passed to print.data.frame().
#'
#' @return Invisibly returns the input object x.
#'
#' @export
print.cluster_motifs_result <- function(x, ...) {
  # Extract dimensions
  n_kmers <- nrow(x)
  n_clusters <- ncol(x)

  cat("Cluster Motif Matrix\n")
  cat("====================\n")
  cat("Dimensions:", nrow(x), "k-mers x", ncol(x), "clusters\n")
  cat("Normalization: Min-Max (0-1 range)\n")
  cat("Value range: [0, 1]\n\n")

  print.data.frame(x, ...)

  invisible(x)
}

################################################################################
#' S3 Method of select_motifs.R
#'
#' Print Method for Selected Motifs
#'
#' @param x An object of class 'select_motifs'
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the input object
#'
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

################################################################################
#' S3 Method of select_sample_train_validation.R

#' @export
print.sample_selection <- function(x, ...) {
  cat("\n=== Sample Selection Summary ===\n\n")
  cat("Classification sequences:", attr(x, "n_classification"), "\n")
  cat("Validation sequences:    ", attr(x, "n_validation"), "\n")
  cat("Validation type:         ", attr(x, "validation_type"), "\n")
  cat("Min sequence size:       ", attr(x, "min_size"), "\n")
  cat("Sequences per class:     ", attr(x, "seq_per_class"), "\n\n")

  cat("Class distribution (classification):\n")
  print(table(x$classification_metadata$class))

  if (attr(x, "n_validation") > 0) {
    cat("\nClass distribution (validation):\n")
    print(table(x$validation_metadata$class))
  }

  cat("\n")
  invisible(x)
}

################################################################################
#' S3 Method of train_models_rf_xgboost.R
#' @export
print.train_models_rf_xgboost <- function(x, ...) {
  cat("\n=== Train Models RF/XGBoost Summary ===\n\n")
  cat("Training configuration:\n")
  cat("  CV folds:", x$cv_folds, "\n")
  cat("  Train proportion:", x$prop_train, "\n")
  cat("  Motifs used:", length(x$motifs_used), "\n")
  cat("  Classes:", nlevels(x$actuals_test), "\n\n")
  cat("Performance comparison:\n")
  print(x$model_comparison)
  cat("\nBest model (test):", x$best_model_test, "\n")
  cat("Best model (validation):", x$best_model_validation, "\n")
  invisible(x)
}

################################################################################
#' S3 Method of kmer_analysis.R
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

################################################################################
#' S3 Method of kmer_in_seq.R
#' @export
print.kmers_in_seq_result <- function(x, ...) {
  cat("\n=== K-mer Search Results ===\n\n")
  cat(sprintf("Training sequences: %d\n", attr(x, "n_train_sequences")))
  if (attr(x, "has_validation")) {
    cat(sprintf("Validation sequences: %d\n", attr(x, "n_validation_sequences")))
  }
  cat(sprintf("Motifs searched: %d\n", attr(x, "n_motifs")))
  cat(sprintf("Total occurrences: %d\n", attr(x, "n_occurrences")))
  cat(sprintf("Elapsed time: %.2f s\n\n", attr(x, "elapsed_time")))

  print(head(as.data.frame(x), 10))
  if (nrow(x) > 10) cat(sprintf("\n... and %d more rows\n", nrow(x) - 10))
  invisible(x)
}

################################################################################
#' S3 Method of seq_classification.R
#' @export
print.seq_classification <- function(x, ...) {
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")
  cat("SEQUENCE CLASSIFICATION RESULTS\n")
  cat(rep("=", 80), "\n\n", sep = "")

  # Pipeline parameters
  cat("PIPELINE PARAMETERS\n")
  cat(rep("-", 80), "\n", sep = "")
  cat("  K-mer size:", x$parameters$k, "\n")
  cat("  Distance method:", x$parameters$dist_method, "\n")
  cat("  Homogeneity threshold:", x$parameters$hom_thresh, "\n")
  cat("  Sequences per class:", x$parameters$seq_per_class, "\n")
  cat("  Minimum sequence length:", x$parameters$min_size, "\n")
  cat("  Top motifs selected:", x$parameters$n_motifs, "\n")
  cat("  Training proportion:", x$parameters$prop_train, "\n")
  cat("  CV folds:", x$parameters$cv_folds, "\n")
  cat("  Validation type:",
      if (is.null(x$parameters$external_validation_fasta_dir)) "Internal" else "External",
      "\n\n")

  # Best model results
  cat("BEST MODEL PERFORMANCE\n")
  cat(rep("-", 80), "\n", sep = "")
  cat("  Model:", x$best_model, "\n")
  cat("  Validation Accuracy:", round(x$validation_accuracy, 4), "\n")
  cat("  Kappa:", round(x$confusion_matrix$overall["Kappa"], 4), "\n")
  cat("  95% CI: (",
      round(x$confusion_matrix$overall["AccuracyLower"], 4), ", ",
      round(x$confusion_matrix$overall["AccuracyUpper"], 4), ")\n\n")

  # Confusion matrix
  cat("CONFUSION MATRIX\n")
  cat(rep("-", 80), "\n", sep = "")
  print(x$confusion_matrix$table)
  cat("\n")

  # Model comparison
  cat("MODEL COMPARISON\n")
  cat(rep("-", 80), "\n", sep = "")
  print(x$model_comparison)
  cat("\n")

  # K-mer analysis summary
  cat("K-MER ANALYSIS SUMMARY\n")
  cat(rep("-", 80), "\n", sep = "")
  cat("  Cluster-specific motifs:", nrow(x$kmer_analysis$unique_cluster_motifs), "\n")
  cat("  Class-specific motifs:", nrow(x$kmer_analysis$unique_class_motifs), "\n")
  cat("  Total motifs analyzed:", nrow(x$kmer_analysis$cluster_frequency_ranking), "\n")

  if (x$kmer_analysis$has_external_validation) {
    cat("  Validation class-specific motifs:",
        nrow(x$kmer_analysis$validation$unique_class_motifs), "\n")
    cat("  Validation motifs analyzed:",
        nrow(x$kmer_analysis$validation$class_frequency_ranking), "\n")
  }
  cat("\n")

  # Processing info
  cat("PROCESSING INFORMATION\n")
  cat(rep("-", 80), "\n", sep = "")
  cat("  Total time:", format(x$processing_time, digits = 2), "\n")
  cat("  Timestamp:", format(x$timestamp, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("\n")
  cat(rep("=", 80), "\n", sep = "")

  invisible(x)
}
