#' Sequence Classification Pipeline
#'
#' Complete pipeline for sequence classification using k-mer analysis, clustering,
#' motif selection, and machine learning models (Random Forest and XGBoost).
#'
#' @param input_dir Character. Path to directory containing FASTA files for classification.
#' @param k Integer. K-mer size for sequence analysis. Default: 6.
#' @param dist_method Character. Distance method for clustering. Default: "euclidean".
#' @param hom_thresh Numeric. Homogeneity threshold for cluster validation (0-1). Default: 0.8.
#' @param k_per_class Integer. Number of sequences to select per class for classification. Default: 3.
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
#' @return A list containing classification results, performance metrics, and intermediate outputs
#' @export
seq_classification <- function(input_dir,
                               k = 6,
                               dist_method = "euclidean",
                               hom_thresh = 0.8,
                               k_per_class = 3,
                               min_size = 800,
                               n_motifs = 256,
                               prop_train = 0.7,
                               cv_folds = 10,
                               kmer_analysis_threshold = 0.8,
                               external_validation_fasta_dir = NULL,
                               k_for_external_validation = 6L,
                               verbose = TRUE) {

  # Start timer
  start_time <- Sys.time()

  # ===== PARAMETER VALIDATION =====

  # Validate hom_thresh
  if (!is.numeric(hom_thresh) || length(hom_thresh) != 1 ||
      is.na(hom_thresh) || hom_thresh < 0 || hom_thresh > 1) {
    stop("'hom_thresh' must be a numeric value between 0 and 1", call. = FALSE)
  }

  # Validate k
  if (!is.numeric(k) || length(k) != 1 || is.na(k) ||
      k <= 0 || k != as.integer(k)) {
    stop("'k' must be a single positive integer", call. = FALSE)
  }

  # Validate k_per_class
  if (!is.numeric(k_per_class) || length(k_per_class) != 1 ||
      is.na(k_per_class) || k_per_class <= 0 ||
      k_per_class != as.integer(k_per_class)) {
    stop("'k_per_class' must be a single positive integer", call. = FALSE)
  }

  # Validate min_size
  if (!is.numeric(min_size) || length(min_size) != 1 ||
      is.na(min_size) || min_size <= 0 ||
      min_size != as.integer(min_size)) {
    stop("'min_size' must be a single positive integer", call. = FALSE)
  }

  # Validate n_motifs
  if (!is.numeric(n_motifs) || length(n_motifs) != 1 ||
      is.na(n_motifs) || n_motifs <= 0 ||
      n_motifs != as.integer(n_motifs)) {
    stop("'n_motifs' must be a single positive integer", call. = FALSE)
  }

  # Validate prop_train
  if (!is.numeric(prop_train) || length(prop_train) != 1 ||
      is.na(prop_train) || prop_train <= 0 || prop_train >= 1) {
    stop("'prop_train' must be a numeric value between 0 and 1", call. = FALSE)
  }

  # Validate cv_folds
  if (!is.numeric(cv_folds) || length(cv_folds) != 1 ||
      is.na(cv_folds) || cv_folds <= 0 ||
      cv_folds != as.integer(cv_folds)) {
    stop("'cv_folds' must be a single positive integer", call. = FALSE)
  }

  # Validate kmer_analysis_threshold
  if (!is.numeric(kmer_analysis_threshold) || length(kmer_analysis_threshold) != 1 ||
      is.na(kmer_analysis_threshold) || kmer_analysis_threshold < 0 ||
      kmer_analysis_threshold > 1) {
    stop("'kmer_analysis_threshold' must be a numeric value between 0 and 1", call. = FALSE)
  }

  # Validate dist_method
  if (!is.character(dist_method) || length(dist_method) != 1) {
    stop("'dist_method' must be a single character string", call. = FALSE)
  }

  valid_methods <- c("euclidean", "maximum", "manhattan", "canberra",
                     "binary", "minkowski")
  if (!dist_method %in% valid_methods) {
    stop("'dist_method' must be one of: ",
         paste(valid_methods, collapse = ", "), call. = FALSE)
  }

  # Validate input_dir
  if (!is.character(input_dir) || length(input_dir) != 1) {
    stop("'input_dir' must be a single character string", call. = FALSE)
  }

  if (!dir.exists(input_dir)) {
    stop("Directory does not exist: ", input_dir, call. = FALSE)
  }

  # Validate external_validation_fasta_dir
  if (!is.null(external_validation_fasta_dir)) {
    if (!is.character(external_validation_fasta_dir) ||
        length(external_validation_fasta_dir) != 1) {
      stop("'external_validation_fasta_dir' must be a single character string", call. = FALSE)
    }
    if (!dir.exists(external_validation_fasta_dir)) {
      stop("'external_validation_fasta_dir' does not exist: ",
           external_validation_fasta_dir, call. = FALSE)
    }
  }

  # ===== START PIPELINE =====

  if (verbose) {
    cat("\n")
    cat(strrep("=", 80), "\n")
    cat("  SEQUENCE CLASSIFICATION PIPELINE\n")
    cat(strrep("=", 80), "\n\n")
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

  # ===== STEP 1: CREATE DATA =====
  print_step(1, 8, "K-mer counting and data creation")

  data_result <- create_data(input = input_dir, k = k)

  print_step(1, 8, "", sprintf(
    "  Sequences: %d | K-mers: %d | Classes: %d | K-mer size: %d\n\n",
    nrow(data_result$metadata), ncol(data_result$kmers),
    length(unique(data_result$metadata$class)), k
  ))

  intermediate$data_result <- data_result

  # ===== STEP 2: CREATE DENDROGRAM =====
  print_step(2, 8, "Creating hierarchical dendrogram")

  dend_result <- create_dendrogram(
    data = data_result$kmers,
    sequence_names = data_result$metadata$sequence_name,
    dist_method = dist_method,
    plot_title = "Hierarchical Clustering Dendrogram"
  )

  print_step(2, 8, "", sprintf("  Distance: %s | Dendrogram created\n\n", dist_method))

  intermediate$dend_result <- dend_result

  # ===== STEP 3: CLUSTER VALIDATION =====
  print_step(3, 8, "Cluster validation and homogeneity assessment")

  cluster_result <- cluster_dendrogram(
    dendrogram = dend_result$dendrogram,
    class_labels = data_result$metadata$class[dend_result$order],
    hom_thresh = hom_thresh,
    min_size = 1,
    dendro_order = dend_result$order,
    sequence_names = data_result$metadata$sequence_name,
    data_result = data_result
  )

  print_step(3, 8, "", sprintf(
    "  Threshold: %.2f | Valid clusters: %d\n\n",
    hom_thresh, length(cluster_result$valid_clusters)
  ))

  intermediate$cluster_result <- cluster_result

  # ===== STEP 4: CALCULATE MOTIFS =====
  print_step(4, 8, "Calculating cluster motifs")

  motif_result <- calculate_cluster_motifs(cluster_result)

  print_step(4, 8, "", "  Motifs calculated\n\n")

  intermediate$motif_result <- motif_result

  # ===== STEP 5: SELECT MOTIFS =====
  print_step(5, 8, "Selecting top discriminative motifs")

  select_motifs_result <- select_motifs(motif_result, cluster_result, n = n_motifs, verbose = verbose)

  print_step(5, 8, "", sprintf("  Top motifs per class: %d\n\n", n_motifs))

  intermediate$select_motifs_result <- select_motifs_result

  # ===== STEP 6: PREPARE DATASETS =====
  print_step(6, 8, "Preparing train/validation datasets")

  datasets_traintest <- select_sample_train_validation(
    cluster_result = cluster_result,
    k_per_class = k_per_class,
    min_size = min_size,
    external_validation_fasta_dir = external_validation_fasta_dir,
    k_for_external_validation = k_for_external_validation
  )

  print_step(6, 8, "", sprintf(
    "  Classification: %d | Validation: %d | Type: %s\n\n",
    nrow(datasets_traintest$classification_dataset),
    nrow(datasets_traintest$validation_dataset),
    if (is.null(external_validation_fasta_dir)) "internal" else "external"
  ))

  intermediate$datasets_traintest <- datasets_traintest

  # ===== STEP 7: TRAIN MODELS =====
  print_step(7, 8, "Training classification models")
  if (verbose) cat("\n")

  classification_result <- train_models_rf_xgboost(
    dataset_traintest = datasets_traintest,
    selected_motifs = select_motifs_result,
    prop_train = prop_train,
    cv_folds = cv_folds,
    seed = 123
  )

  # ===== STEP 8: K-MER ANALYSIS =====
  if (verbose) cat("\n")
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

  print_step(8, 8, "", sprintf(
    "  Cluster-specific: %d | Class-specific: %d | Total: %d%s\n\n",
    nrow(kmer_analysis_result$unique_cluster_motifs),
    nrow(kmer_analysis_result$unique_class_motifs),
    nrow(kmer_analysis_result$cluster_frequency_ranking),
    if (kmer_analysis_result$has_external_validation)
      sprintf(" | Val: %d", nrow(kmer_analysis_result$validation$class_frequency_ranking))
    else ""
  ))

  intermediate$kmer_analysis <- kmer_analysis_result

  # ===== EXTRACT BEST MODEL =====
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

  # ===== FINAL SUMMARY =====
  if (verbose) {
    cat("\n", rep("=", 80), "\n", sep = "")
    cat("  PIPELINE COMPLETED SUCCESSFULLY\n")
    cat(rep("=", 80), "\n\n", sep = "")
  }

  cat("\n", rep("=", 80), "\n", sep = "")
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

  cat(rep("=", 80), "\n", sep = "")
  cat("Total Processing Time:", format(processing_time, digits = 2), "\n")
  cat(rep("=", 80), "\n\n", sep = "")

  # ===== RETURN OBJECT =====
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
      k_per_class = k_per_class,
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
