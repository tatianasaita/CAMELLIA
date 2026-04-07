#' Sequence Classification Pipeline
#'
#' Complete pipeline for sequence classification using k-mer analysis, clustering,
#' motif selection, and machine learning models (Random Forest and XGBoost).
#'
#' @param input_dir Character. Path to directory containing FASTA files for classification.
#' @param k Integer. K-mer size for sequence analysis. Default: 6.
#' @param dist_method A character string specifying the distance method (default: "euclidean").
#'   See method options in \code{\link[parallelDist]{parDist}}.
#' @param hom_thresh Numeric. Homogeneity threshold for cluster validation (0-1). Default: 0.8.
#' @param seq_per_class Integer. Number of sequences to select per class for classification. Default: 3.
#' @param min_size Integer. Minimum sequence length. Default: 800.
#' @param n_motifs Integer. Number of top motifs to select per class. Default: 256.
#' @param prop_train Numeric. Proportion of data for training (0-1). Default: 0.7.
#' @param cv_folds Integer. Number of cross-validation folds. Default: 10.
#' @param kmer_analysis_threshold Numeric. Threshold for k-mer analysis. Default: 0.8.
#' @param external_test_fasta_dir Character. Optional path to directory with external
#'   test FASTA files. If NULL, uses internal test split. Default: NULL.
#' @param k_for_external_test Integer. K-mer size for external test set. Default: 6L.
#' @param verbose Logical. If TRUE, shows intermediate results and messages. Default: TRUE.
#' @param cluster_method Character. Clustering method: "dendrogram" or "apcluster".
#'   Default: "dendrogram".
#' @param ap_r Numeric. Distance metric parameter for apcluster. Default: 2.
#' @param ml_method Character. ML method: "rf" or "xgb". Default: "xgb"
#'
#' @return A list of class "seq_classification" containing:
#' \itemize{
#'   \item \code{classification_results}: Model training results
#'   \item \code{best_model}: Name of best performing model
#'   \item \code{test_predictions_rf}: Random Forest predictions on test set
#'   \item \code{test_predictions_xgb}: XGBoost predictions on test set
#'   \item \code{model_comparison}: Accuracy comparison between models
#'   \item \code{confusion_matrix}: Confusion matrix for best model
#'   \item \code{kmer_analysis}: K-mer analysis results (cluster and class-specific motifs)
#'   \item \code{processing_time}: Total pipeline execution time
#'   \item \code{parameters}: List of all input parameters
#'   \item \code{intermediate_results}: Intermediate results from all 8 steps (if verbose=TRUE)
#'   \item \code{timestamp}: Execution timestamp
#' }
#'
#' @details
#' This function executes a complete 8-step pipeline:
#' \enumerate{
#'   \item K-mer counting and data creation
#'   \item Hierarchical clustering dendrogram
#'   \item Cluster validation and homogeneity assessment
#'   \item Cluster motif calculation
#'   \item Top discriminative motif selection
#'   \item Train/validation/test dataset preparation
#'   \item Model training (Random Forest and XGBoost)
#'   \item K-mer analysis for cluster and class-specific motifs
#' }
#'
#'@note
#'
#' @examples
#' \dontrun{
#' # Basic usage with internal test split — dendrogram method
#' result_seq_classification <- seq_classification_cent(
#'   input_dir     = "path/to/fasta_files",
#'   k             = 6,
#'   seq_per_class = 200,
#'   min_size      = 800
#' )
#'
#' # Basic usage with internal test split — apcluster method
#' result_seq_classification <- seq_classification_cent(
#'   input_dir     = "path/to/fasta_files",
#'   k             = 6,
#'   seq_per_class = 200,
#'   min_size      = 800,
#'   cluster_method = "apcluster",
#'   ap_r          = 2
#' )
#'
#' # With external test dataset
#' result_seq_classification <- seq_classification_cent(
#'   input_dir               = "path/to/fasta_files",
#'   k                       = 6,
#'   seq_per_class           = 200,
#'   min_size                = 800,
#'   external_test_fasta_dir = "path/to/test/fasta"
#' )
#' }
#'
setwd("E:/TATIANA/CAMELLIA_V1") ##################################Changing directories
source("create_data.R")
source("create_dendrogram_cent.R")
source("cluster_dendrogram_ap_total.R") 
source("calculate_cluster_motifs.R")
source("select_motifs.R")
source("select_train_test1.R")
source("train_model_xgboost_rf.R")
source("kmer_analysis.R")
source("internal-functions.R")

#' @export
seq_classification_cent <- function(input_dir,
                                    k                       = 6,
                                    dist_method             = "euclidean",
                                    hom_thresh              = 0.8,
                                    n_motifs                = 256,
                                    seq_per_class           = 200,
                                    min_size                = 800,
                                    prop_train              = 0.7,
                                    cv_folds                = 10,
                                    external_test_fasta_dir = NULL,
                                    k_for_external_test     = 6L,
                                    kmer_analysis_threshold = 0.8,
                                    cluster_method          = c("dendrogram", "apcluster"),
                                    ap_r                    = 2,
                                    verbose                 = TRUE) {

  cluster_method <- match.arg(cluster_method)

  # Start timer
  start_time <- Sys.time()
  is_external <- !is.null(external_test_fasta_dir)
  set.seed(123)
  intermediate <- list()

  total_steps <- if (cluster_method == "dendrogram") 8L else 7L

  print_step <- function(step, title, details = NULL) {
    if (!verbose) return()
    cat("STEP ", step, "/", total_steps, ": ", title, "\n", sep = "")
    cat(rep("-", 80), "\n", sep = "")
    if (!is.null(details)) cat(details)
  }

  if (verbose) {
    cat("\n SEQUENCE CLASSIFICATION PIPELINE\n")
    cat(sprintf("Clustering method: %s\n", cluster_method))
    cat(sprintf("External test dataset: %s\n\n", if (is_external) "YES" else "NO"))
  }

  # Step 1: Create data
  print_step(1, "K-mer counting and data creation")

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
  if (cluster_method == "dendrogram") {

    print_step(2, "Creating hierarchical dendrogram")

    dend_result <- create_dendrogram_cent(
      data           = data_result$kmers,
      sequence_names = data_result$metadata$sequence_name,
      dist_method    = dist_method,
      hclust_method  = "ward.D2",
      seq_lengths    = data_result$metadata$length,
      min_seq_length = min_size
    )

    if (verbose) {
      n_used <- if (!is.null(dend_result$selected_indices))
        length(dend_result$selected_indices)
      else
        nrow(data_result$metadata)

      cat(sprintf("  Method: Dendrogram | Sequences used: %d\n", n_used))
      cat(sprintf("  Distance: %s\n\n", dist_method))
    }

  } else {
    # Affinity Propagation method
    dend_result <- NULL
    if (verbose) {
      cat(sprintf("  Method: Affinity Propagation (AP) | r parameter: %d\n\n", ap_r))
    }
  }

  intermediate$dend_result <- dend_result

  # Dendrogram Plot and Legend
  # plot(result_seq_classification$intermediate_results$dend_result$dendrogram)
  # legend("topright",
  #       legend  = names(result_seq_classification$intermediate_results$dend_result$base_colors),
  #        col     = result_seq_classification$intermediate_results$dend_result$base_colors,
  #        pch     = 15,
  #        pt.cex  = 2,
  #        cex     = 0.8,
  #        title   = "Classes",
  #        bg      = "white",
  #        box.lty = 1)

  # Step 3: Cluster validation
  print_step(
    if (cluster_method == "dendrogram") 3L else 2L,
    "Cluster validation and homogeneity assessment"
  )

  if (cluster_method == "dendrogram") {
    if (!is.null(dend_result$selected_indices)) {
      data_result_dend          <- data_result
      data_result_dend$metadata <- data_result$metadata[dend_result$selected_indices, ]
      data_result_dend$kmers    <- data_result$kmers[dend_result$selected_indices, ]
    } else {
      data_result_dend <- data_result
    }

    cluster_result <- cluster_dendrogram(
      dendrogram     = dend_result$dendrogram,
      class_labels   = dend_result$labels,
      hom_thresh     = hom_thresh,
      sequence_names = data_result_dend$metadata$sequence_name,
      data_result    = data_result_dend,
      method         = "dendrogram"
    )
  } else {
    # Affinity Propagation clustering
    feature_mat <- as.matrix(
      data_result$kmers[, !colnames(data_result$kmers) %in% "CLASS"]
    )
    rownames(feature_mat) <- data_result$metadata$sequence_name

    cluster_result <- cluster_dendrogram(
      class_labels   = data_result$metadata$class,
      method         = "apcluster",
      feature_matrix = feature_mat,
      ap_r           = ap_r,
      sequence_names = data_result$metadata$sequence_name,
      data_result    = data_result,
#      seq_lengths    = data_result$metadata$length,
#      min_seq_length = min_size,
      verbose        = verbose
    )
  }

  if (verbose) {
    cat(sprintf(
      "  Method: %s | Threshold: %s | Valid clusters: %d\n\n",
      cluster_method,
      ifelse(cluster_method == "dendrogram",
             sprintf("%.2f", hom_thresh),
             "N/A"),
      length(cluster_result$valid_clusters)
    ))
  }

  intermediate$cluster_result <- cluster_result

  # Step 4: Calculate motifs
  print_step(
    if (cluster_method == "dendrogram") 4L else 3L,
    "Calculating cluster motifs"
  )

  motif_result <- calculate_cluster_motifs(cluster_result)

  if (verbose) cat(sprintf("  Motifs calculated for %d clusters\n\n", length(motif_result)))

  intermediate$motif_result <- motif_result

  # Step 5: Select motifs
  print_step(
    if (cluster_method == "dendrogram") 5L else 4L,
    "Selecting top discriminative motifs"
  )

  select_motifs_result <- select_motifs(
    motif_result, cluster_result, n = n_motifs, verbose = verbose
  )

  if (verbose) cat(sprintf("  Top motifs per class: %d\n\n", n_motifs))

  intermediate$select_motifs_result <- select_motifs_result

  # Step 6: Prepare datasets
  print_step(
    if (cluster_method == "dendrogram") 6L else 5L,
    "Preparing train/validation/test datasets"
  )

  datasets_traintest <- select_train_test1(
    data_result             = data_result,
    seq_per_class           = seq_per_class,
    min_size                = min_size,
    external_test_fasta_dir = external_test_fasta_dir,
    verbose                 = verbose
  )

  if (verbose) {
    cat(sprintf(
      "  Training: %d | Test: %d | Type: %s\n\n",
      nrow(datasets_traintest$train_dataset),
      nrow(datasets_traintest$test_dataset),
      if (is_external) "external" else "internal"
    ))
  }

  intermediate$datasets_traintest <- datasets_traintest

  # Step 7: Train models
  print_step(
    if (cluster_method == "dendrogram") 7L else 6L,
    "Training classification models"
  )

  classification_result <- train_model_xgboost_rf(
    result_train_test      = datasets_traintest,
    result_selected_motifs = select_motifs_result,
    verbose                = verbose
  )

  if (verbose) {
    cat(sprintf("  Models trained and evaluated\n\n"))
  }

  intermediate$classification_result <- classification_result

  # Step 8: K-mer analysis
  print_step(
    if (cluster_method == "dendrogram") 8L else 7L,
    "K-mer analysis (cluster and class-specific motifs)"
  )

  kmer_analysis_result <- kmer_analysis(
    cluster_result = cluster_result,
    motif_matrix   = motif_result
  )

  if (verbose) {
    cat(sprintf("  Motifs analyzed in sequences\n"))
    if (!is.null(kmer_analysis_result)) {
      cat(sprintf("  Total occurrences found: %d\n\n", nrow(kmer_analysis_result)))
    }
  }

  intermediate$kmer_analysis <- kmer_analysis_result

  processing_time <- Sys.time() - start_time

  # Final summary
  if (verbose) {
    cat(rep("=", 80), "\n", sep = "")
    cat("PIPELINE COMPLETED SUCCESSFULLY\n")
    cat(rep("=", 80), "\n\n", sep = "")

#    cat("Model Performance Summary:\n")
#    print(classification_result$model_comparison)
#    cat("\n")

    if (is_external) {
      cat("External Test Set Results:\n")
      cat("  Status: UNLABELLED sequences classified by both models\n")
      cat("  Random Forest predictions: ", nrow(classification_result$predictions_test_ext_rf), " sequences\n")
      cat("  XGBoost predictions      : ", nrow(classification_result$predictions_test_ext_xgb), " sequences\n\n")
    }

#    cat("K-mer Analysis Summary:\n")
#    if (!is.null(kmer_analysis_result)) {
#      cat("  Total motif occurrences: ", nrow(kmer_analysis_result), "\n")
#    }

    cat("\nTotal Processing Time: ", format(processing_time, digits = 2), "\n\n", sep = "")
  }

  # Return object
  structure(
    list(
      classification_results = classification_result,
      test_predictions_rf    = classification_result$predictions_test_ext_rf,
      test_predictions_xgb   = classification_result$predictions_test_ext_xgb,
      model_comparison       = classification_result$model_comparison,
      kmer_analysis          = kmer_analysis_result,
      processing_time        = processing_time,
      parameters = list(
        input_dir               = input_dir,
        k                       = k,
        hom_thresh              = hom_thresh,
        seq_per_class           = seq_per_class,
        min_size                = min_size,
        n_motifs                = n_motifs,
        prop_train              = prop_train,
        cv_folds                = cv_folds,
        external_test_fasta_dir = external_test_fasta_dir,
        k_for_external_test     = k_for_external_test,
        dist_method             = dist_method,
        kmer_analysis_threshold = kmer_analysis_threshold,
        cluster_method          = cluster_method,
        ap_r                    = ap_r
      ),
      intermediate_results = if (verbose) intermediate else NULL,
      timestamp            = Sys.time()
    ),
    class = "seq_classification"
  )
}

