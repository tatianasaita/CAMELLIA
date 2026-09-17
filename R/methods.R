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
#' Print Method for Dendrogram Clustering Result
#'
#' @param x An object of class \code{cluster_dendrogram_result}
#' @param ... Additional arguments passed to print methods
#'
#' @return Invisibly returns \code{x}
#'
#' @export
print.cluster_dendrogram_result <- function(x, ...) {
  cat("\nHierarchical Clustering Result\n")
  cat(paste(rep("=", 30), collapse = ""), "\n\n")
  
  cat(sprintf("Total elements:     %d\n", x$n_elements))
  cat(sprintf("Total clusters:     %d\n", nrow(x$cluster_summary)))
  cat(sprintf("Unassigned:         %d\n", x$n_unassigned))
  cat(sprintf("Min size:           %d\n", x$min_size_input))
  cat(sprintf("Hom threshold:      %.3f\n", x$hom_thresh_input))
  
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
#' @method print cluster_motifs
print.cluster_motifs <- function(x, ...) {
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
#' S3 Method of select_train_test.R
#' @export
#' @method print train_test_selection
print.train_test_selection <- function(x, ...) {
  cat("\n=== Train/Test Selection Summary ===\n\n")
  cat("Mode                :", attr(x, "mode"), "\n")
  cat("Train sequences     :", attr(x, "n_train"), "\n")
  cat("Test sequences      :", attr(x, "n_test"), "\n")
  cat("Min sequence size   :", attr(x, "min_size"), "\n")
  cat("Sequences per class :", attr(x, "seq_per_class"), "\n\n")
  
  cat("Class distribution (train):\n")
  print(table(x$train_metadata$class))
  
  if (attr(x, "n_test") > 0) {
    cat("\nClass distribution (test):\n")
    print(table(x$test_metadata$class))
  }
  
  cat("\n")
  invisible(x)
}

################################################################################
#' S3 Method of train_model_xgboost_rf.R
#' @export
#' @method print train_model_xgboost
print.train_model_xgboost <- function(x, ...) {
  cat("\n=== Model Training Summary ===\n\n")
  cat("Method            :", toupper(x$method), "\n")
  cat("CV folds          :", x$cv_folds, "\n")
  cat("Motifs used       :", length(x$motifs_used), "\n")
  cat("Classes           :", nlevels(x$actuals_train), "\n\n")
  
  cat("Model metrics:\n")
  print(x$model_metrics)
  
  cat("\nTime elapsed (s)  :", x$time_seconds, "\n")
  invisible(x)
}

################################################################################
#' S3 Method of kmer_analysis.R
#' @export
print.kmer_analysis_result <- function(x, ...) {
  cat("\nK-mer Analysis Results\n")
  cat("======================\n\n")
  cat(sprintf("Cluster-specific motifs:  %d\n", nrow(x$unique_cluster_motifs)))
  cat(sprintf("Class-specific motifs:    %d\n", nrow(x$unique_class_motifs)))
  cat(sprintf("Total motifs analyzed:    %d\n", nrow(x$cluster_frequency_ranking)))
  cat(sprintf("Total clusters:           %d\n", length(x$cluster_to_class)))
  cat(sprintf("Total classes:            %d\n", ncol(x$class_frequency_matrix) - 1))
  cat("\n")
  invisible(x)
}

################################################################################
#' S3 Method of kmer_in_seq.R
#' @export
print.kmers_in_seq_result <- function(x, ...) {
  cat("\n=== K-mer Search Results ===\n\n")
  cat(sprintf("Training sequences: %d\n", attr(x, "n_train_sequences")))
  if (attr(x, "has_test")) {
    cat(sprintf("Test sequences: %d\n", attr(x, "n_test_sequences")))
  }
  cat(sprintf("Motifs searched: %d\n",    attr(x, "n_motifs")))
  cat(sprintf("Total occurrences: %d\n",  attr(x, "n_occurrences")))
  cat(sprintf("Elapsed time: %.2f s\n\n", attr(x, "elapsed_time")))

  print(head(as.data.frame(x), 10))
  if (nrow(x) > 10) cat(sprintf("\n... and %d more rows\n", nrow(x) - 10))
  invisible(x)
}

################################################################################
#' S3 Method of seq_classification_cent.R
#' @title Print method for seq_classification
#' @description Displays a concise summary of the sequence classification pipeline results.
#' @param x An object of class \code{seq_classification}.
#' @param ... Additional arguments (ignored).
#' @export
print.seq_classification <- function(x, ...) {
  cat("\n")
  cat(rep("=", 60), "\n", sep = "")
  cat("  Sequence Classification Pipeline\n")
  cat(rep("=", 60), "\n", sep = "")
  
  cat("\nMethod              :", toupper(x$classification_results$method), "\n")
  cat("Model metrics:\n")
  print(x$model_metrics)
  
  cat("\nParameters:\n")
  cat("  K-mer size        :", x$parameters$k, "\n")
  cat("  Homogeneity thresh:", x$parameters$hom_thresh, "\n")
  cat("  Sequences/class   :", x$parameters$seq_per_class, "\n")
  cat("  Min. seq. length  :", x$parameters$min_size, "\n")
  cat("  Top motifs        :", x$parameters$n_motifs, "\n")
  cat("  Training prop.    :", x$parameters$prop_train, "\n")
  cat("  CV folds          :", x$parameters$cv_folds, "\n")
  cat("  Test type         :",
      ifelse(is.null(x$parameters$external_test_fasta_dir),
             "internal", "external"), "\n")
  
  cat("\nK-mer Analysis:\n")
  cat("  Cluster-specific motifs :", nrow(x$kmer_analysis$unique_cluster_motifs), "\n")
  cat("  Class-specific motifs   :", nrow(x$kmer_analysis$unique_class_motifs), "\n")
  cat("  Total motifs analyzed   :", nrow(x$kmer_analysis$cluster_frequency_ranking), "\n")
  
  cat("\nProcessing Time    :", format(x$processing_time, digits = 2), "\n")
  cat("Timestamp          :", format(x$timestamp, "%Y-%m-%d %H:%M:%S"), "\n")
  cat(rep("=", 60), "\n", sep = "")
  
  invisible(x)
} 
