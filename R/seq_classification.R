#' Sequence Classification Pipeline
#'
#' Complete pipeline for sequence classification using k-mer analysis, clustering,
#' motif selection, and machine learning models (Random Forest and XGBoost).
#'
#' @param input_dir Character. Path to directory containing FASTA files for classification.
#' @param k Integer. K-mer size for sequence analysis. Default: 6.
#' @param dist_method A character string specifying the distance method (default: "euclidean"). See method options in \code{\link[parallelDist]{parDist}}.
#' @param hom_thresh Numeric. Homogeneity threshold for cluster validation (0-1). Default: 0.8.
#' @param seq_per_class Integer. Number of sequences to select per class for classification. Default: 3.
#' @param min_size Integer. Minimum sequence length. Default: 800.
#' @param n_motifs Integer. Number of top motifs to select per class. Default: 256.
#' @param prop_train Numeric. Proportion of data for training (0-1). Default: 0.7.
#' @param cv_folds Integer. Number of cross-validation folds. Default: 10.
#' @param kmer_analysis_threshold Numeric. Threshold for k-mer analysis. Default: 0.8.
#' @param external_validation_fasta_dir Character. Optional path to directory with external
#'   validation FASTA files. If NULL, uses internal validation. Default: NULL.
#' @param k_for_external_validation Integer. K-mer size for external validation. Default: 6L.
#' @param verbose Logical. If TRUE, shows intermediate results and messages. Default: TRUE.
#'
#' @return A list of class "seq_classification" containing:
#' \itemize{
#'   \item \code{classification_results}: Model training and validation results
#'   \item \code{best_model}: Name of best performing model
#'   \item \code{validation_predictions}: Predictions on validation set
#'   \item \code{validation_accuracy}: Accuracy on validation set
#'   \item \code{confusion_matrix}: Confusion matrix for best model
#'   \item \code{kmer_analysis}: K-mer analysis results
#'   \item \code{processing_time}: Total pipeline execution time
#'   \item \code{parameters}: List of all input parameters
#' }
#'
#' #' @details
#' This function executes a complete 8-step pipeline:
#' \enumerate{
#'   \item K-mer counting and data creation
#'   \item Hierarchical clustering dendrogram
#'   \item Cluster validation and homogeneity assessment
#'   \item Cluster motif calculation
#'   \item Top discriminative motif selection
#'   \item Train/validation dataset preparation
#'   \item Model training (Random Forest and XGBoost)
#'   \item K-mer analysis for cluster and class-specific motifs
#' }
#'
#'@note
#'
#' @examples
#' \dontrun{
#' # Basic usage with internal validation
#' result_seq_classification <- seq_classification(
#'   input_dir = "path/to/fasta_files",
#'   k = 6,
#'   seq_per_class = 200,
#'   min_size = 800
#'  )
#' }
#'
setwd("E:/TATIANA/CAMELLIA-main/R")
source("create_data.R")
source("create_dendrogram.R")
source("cluster_dendrogram.R")
source("calculate_cluster_motifs.R")
source("select_motifs.R")
source("select_sample_train_validation.R")
source("train_models_rf_xgboost.R")
source("kmer_analysis.R")
source("internal-functions.R")
source("methods.R")

#' @export
seq_classification <- function(input_dir,
                               k = 6,
                               dist_method = "euclidean",
                               hom_thresh = 0.8,
                               n_motifs = 256,
                               seq_per_class = 200,
                               min_size = 800,
                               prop_train = 0.7,
                               cv_folds = 10,
                               external_validation_fasta_dir = NULL,
                               k_for_external_validation = 6L,
                               kmer_analysis_threshold = 0.8,
                               verbose = TRUE) {

# Start timer
  start_time <- Sys.time()

# Start pipeline
  if (verbose) {
    cat(" \n SEQUENCE CLASSIFICATION PIPELINE\n")
  }

  set.seed(123)
  intermediate <- list()

  # Helper for verbose output
  print_step <- function(step, total, title, details = NULL) {
    if (!verbose) return()
    cat("STEP ", step, "/", total, ": ", title, "\n", sep = "")
    cat(rep("-", 80), "\n", sep = "")
    if (!is.null(details)) cat(details)
  }

# Step 1: Create data
  print_step(1, 8, "K-mer counting and data creation")

  data_result <- create_data(input = input_dir, k = k)

  if (verbose) {
    cat(sprintf(
      "  Sequences: %d | K-mers: %d | Classes: %d | K-mer size: %d\n\n",
      nrow(data_result$metadata), ncol(data_result$kmers),
      length(unique(data_result$metadata$class)), k
    ))
  }

  intermediate$data_result <- data_result

# Step 2: Create dendrogram
  print_step(2, 8, "Creating hierarchical dendrogram")

  dend_result <- create_dendrogram(
    data = data_result$kmers,
    sequence_names = data_result$metadata$sequence_name
  )

  if (verbose) {
    cat(sprintf("  Distance: %s | Dendrogram created\n\n", dist_method))}

  intermediate$dend_result <- dend_result

# Step 3: Cluster validation
  print_step(3, 8, "Cluster validation and homogeneity assessment")

  cluster_result <- cluster_dendrogram(
    dendrogram = dend_result$dendrogram,
    class_labels = data_result$metadata$class[dend_result$order],
    hom_thresh = hom_thresh,
    sequence_names = data_result$metadata$sequence_name,
    data_result = data_result
  )

  if (verbose) {
    cat(sprintf(
    "  Threshold: %.2f | Valid clusters: %d\n\n",
    hom_thresh, length(cluster_result$valid_clusters)
  ))
  }

  intermediate$cluster_result <- cluster_result

  # Step 4: Calculate motifs
  print_step(4, 8, "Calculating cluster motifs")

  motif_result <- calculate_cluster_motifs(cluster_result)

  if (verbose) {
    cat(sprintf("  Motifs calculated\n\n"))
  }

  intermediate$motif_result <- motif_result

  # Step 5: Select motifs
  print_step(5, 8, "Selecting top discriminative motifs")

  select_motifs_result <- select_motifs(motif_result, cluster_result, n = n_motifs, verbose = verbose)

  if (verbose) {
    cat(sprintf("  Top motifs per class: %d\n\n", n_motifs))
  }

  intermediate$select_motifs_result <- select_motifs_result

  # Step 6: Prepare datasets
  print_step(6, 8, "Preparing train/validation datasets")

  datasets_traintest <- select_sample_train_validation(
    cluster_result = cluster_result,
    seq_per_class = 200,
    min_size = min_size,
    external_validation_fasta_dir = external_validation_fasta_dir,
    k_for_external_validation = k_for_external_validation
  )

  if (verbose) {
    cat(sprintf(
    "  Classification: %d | Validation: %d | Type: %s\n\n",
    nrow(datasets_traintest$classification_dataset),
    nrow(datasets_traintest$validation_dataset),
    if (is.null(external_validation_fasta_dir)) "internal" else "external"
  ))
  }

  intermediate$datasets_traintest <- datasets_traintest

  # Step 7: Train models
  print_step(7, 8, "Training classification models")

  classification_result <- train_models_rf_xgboost(
    dataset_traintest = datasets_traintest,
    selected_motifs = select_motifs_result,
    prop_train = prop_train,
    cv_folds = cv_folds
  )

  # Step 8: K-mer analysis
  print_step(8, 8, "K-mer analysis (cluster and class-specific motifs)")

  if (!is.null(external_validation_fasta_dir)) {
    cluster_result$classification_result <- list(
      validation_data = datasets_traintest$validation_dataset
    )
  }

  kmer_analysis_result <- kmer_analysis(
    cluster_result = cluster_result,
    motif_matrix = motif_result,
    threshold = kmer_analysis_threshold
  )

  if (verbose) {
    cat(sprintf(
    "  Cluster-specific: %d | Class-specific: %d | Total: %d%s\n\n",
    nrow(kmer_analysis_result$unique_cluster_motifs),
    nrow(kmer_analysis_result$unique_class_motifs),
    nrow(kmer_analysis_result$cluster_frequency_ranking),
    if (kmer_analysis_result$has_external_validation)
      sprintf(" | Val: %d", nrow(kmer_analysis_result$validation$class_frequency_ranking))
    else ""
  ))
}
  intermediate$kmer_analysis <- kmer_analysis_result

  # Extract best model
  best_model <- classification_result$best_model_validation
  best_cm <- if (best_model == "Random Forest") {
    classification_result$confusion_matrix_validation_rf
  } else {
    classification_result$confusion_matrix_validation_xgb
  }
  best_predictions <- if (best_model == "Random Forest") {
    classification_result$predictions_validation_rf
  } else {
    classification_result$predictions_validation_xgb
  }
  best_accuracy <- best_cm$overall["Accuracy"]

  processing_time <- Sys.time() - start_time

  # Final summary
  if (verbose) {
    cat("PIPELINE COMPLETED SUCCESSFULLY\n")

  cat("  FINAL CLASSIFICATION RESULTS\n")
  cat(rep("=", 80), "\n\n", sep = "")
  cat("Best Model:", best_model, "\n")
  cat("Validation Accuracy:", round(best_accuracy, 4), "\n")
  cat("Validation Kappa:", round(best_cm$overall["Kappa"], 4), "\n")
  cat("95% CI: (", round(best_cm$overall["AccuracyLower"], 4), ", ",
      round(best_cm$overall["AccuracyUpper"], 4), ")\n\n")

  cat("Confusion Matrix (Validation Set):\n")
  print(best_cm$table)
  cat("\n")

  cat("Model Performance Comparison:\n")
  print(classification_result$model_comparison)
  cat("\n")

  cat("K-mer Analysis Summary:\n")
  cat("  Cluster-specific motifs:", nrow(kmer_analysis_result$unique_cluster_motifs), "\n")
  cat("  Class-specific motifs:", nrow(kmer_analysis_result$unique_class_motifs), "\n")
  cat("  Total motifs analyzed:", nrow(kmer_analysis_result$cluster_frequency_ranking), "\n")
  if (kmer_analysis_result$has_external_validation) {
    cat("  Validation class-specific motifs:", nrow(kmer_analysis_result$validation$unique_class_motifs), "\n")
    cat("  Validation motifs analyzed:", nrow(kmer_analysis_result$validation$class_frequency_ranking), "\n")
  }
  cat("\n")

  cat("Total Processing Time:", format(processing_time, digits = 2), "\n")
}

  # Return object
  list(
    classification_results = classification_result,
    best_model = best_model,
    validation_predictions = best_predictions,
    validation_actuals = classification_result$actuals_validation,
    validation_accuracy = best_accuracy,
    confusion_matrix = best_cm,
    model_comparison = classification_result$model_comparison,
    kmer_analysis = kmer_analysis_result,
    processing_time = processing_time,
    parameters = list(
      input_dir = input_dir,
      k = k,
      hom_thresh = hom_thresh,
      seq_per_class = seq_per_class,
      min_size = min_size,
      n_motifs = n_motifs,
      prop_train = prop_train,
      cv_folds = cv_folds,
      external_validation_fasta_dir = external_validation_fasta_dir,
      k_for_external_validation = k_for_external_validation,
      dist_method = dist_method,
      kmer_analysis_threshold = kmer_analysis_threshold
    ),
    intermediate_results = if (verbose) intermediate else NULL,
    timestamp = Sys.time()
  )
}
