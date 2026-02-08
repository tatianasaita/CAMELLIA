#' Select Sample for Training and Validation Datasets
#'
#' Splits data into training and validation sets.
#'
#' @param cluster_result List from cluster_dendrogram() containing data_result.
#' @param seq_per_class Integer. Number of sequences per class for training.
#' @param min_size Integer. Minimum sequence length for training (not applied to validation).
#' @param external_validation_fasta_dir Character or NULL. Path to external FASTA files.
#' @param k_for_external_validation Integer. K-mer size for external validation. Default 6.
#' @param verbose Logical. If TRUE, prints detailed information. Default is TRUE.
#'
#' @return A list of class "sample_selection" with training and validation datasets.
#'
#' @details
#' Splits data into training (classification) and validation sets:
#' \itemize{
#'   \item Training: Selects k sequences per class with length >= min_size
#'   \item Validation (internal): Uses remaining sequences (no length filter)
#'   \item Validation (external): Loads from FASTA directory if provided
#' }
#'
#' @note
#' \itemize{
#'   \item Requires .select_by_class. See \code{internal-functions.R} for details.
#'   \item S3 methods available. See \code{methods.R} for details.
#' }
#'
#' @examples
#' \dontrun{
#' # Internal validation (default)
#' result_sample_internal <- select_sample_train_validation(
#'   cluster_result = result_cluster_dendrogram,
#'   seq_per_class = 200,
#'   min_size = 800
#' )
#'
#' # External validation
#' sample_external_result <- select_sample_train_validation(
#'   cluster_result = result_cluster_dendrogram,
#'   seq_per_class = 200,
#'   min_size = 800,
#'   external_validation_fasta_dir = "path/to/validation/fastas",
#'   k_for_external_validation = 6
#' )
#'}
#'
#' @export
select_sample_train_validation <- function(cluster_result,
                                           seq_per_class,
                                           min_size = 800,
                                           external_validation_fasta_dir = NULL,
                                           k_for_external_validation = 6L,
                                           verbose = TRUE) {

  # Extract data from cluster_result
  metadata <- cluster_result$data_result$metadata
  kmers <- cluster_result$data_result$kmers

  # Classification dataset
  metadata_filtered <- metadata[metadata$length >= min_size, ]

  if (nrow(metadata_filtered) == 0) {
    stop("No sequences available with length >= ", min_size, call. = FALSE)
  }

  if (verbose) {
    cat("Sequences with length >=", min_size, ":", nrow(metadata_filtered), "\n")
    cat("Available classes:", paste(unique(metadata_filtered$class), collapse = ", "), "\n")
  }

  # Select sequences per class
  classification_metadata <- .select_by_class(metadata_filtered, seq_per_class, "classification")

  if (nrow(classification_metadata) == 0) {
    stop("Could not select any sequences for classification", call. = FALSE)
  }

  # Extract k-mers using sequence names
  selected_indices <- match(classification_metadata$sequence_name, metadata$sequence_name)
  classification_dataset <- kmers[selected_indices, , drop = FALSE]
  rownames(classification_dataset) <- classification_metadata$sequence_name

  if (verbose) {
    cat("Total classification sequences:", nrow(classification_metadata), "\n\n")
  }

  # Validation dataset
  if (is.null(external_validation_fasta_dir)) {
    # Internal validation
    if (verbose) cat("Internal Validation Dataset \n")

    # Get sequences NOT selected for classification
    remaining_metadata <- metadata[!(metadata$sequence_name %in% classification_metadata$sequence_name), ]

    if (nrow(remaining_metadata) == 0) {
      warning("No sequences available for internal validation", call. = FALSE)
      validation_metadata <- data.frame(
        sequence_name = character(0),
        length = integer(0),
        class = character(0),
        stringsAsFactors = FALSE
      )
      validation_dataset <- kmers[0, , drop = FALSE]
    } else {
      validation_metadata <- remaining_metadata

      if (verbose) {
        cat("Validation sequences:", nrow(validation_metadata), "\n")
        cat("Validation classes:", paste(unique(validation_metadata$class), collapse = ", "), "\n\n")
      }

      # Extract k-mers
      validation_indices <- match(validation_metadata$sequence_name, metadata$sequence_name)
      validation_dataset <- kmers[validation_indices, , drop = FALSE]
      rownames(validation_dataset) <- validation_metadata$sequence_name

      # Verify no overlap
      overlap <- intersect(classification_metadata$sequence_name, validation_metadata$sequence_name)
      if (length(overlap) > 0) {
        warning("Using ", length(overlap), " sequences from classification set (WARNING: Data leakage!)",
                call. = FALSE)
      }
    }

  } else {
    # External validation
    if (verbose) cat("External Validation Dataset\n")

    # Use first directory if multiple provided
    validation_dir <- external_validation_fasta_dir[1]

    external_data_result <- tryCatch({
      create_data(input = validation_dir, k = k_for_external_validation)
    }, error = function(e) {
      stop("Error creating external validation dataset: ", e$message, call. = FALSE)
    })

    validation_metadata <- external_data_result$metadata
    validation_dataset <- external_data_result$kmers

    if (nrow(validation_metadata) == 0) {
      stop("No validation sequences with length >= ", min_size, call. = FALSE)
    }

    if (verbose) {
      cat("Available sequences by class:\n")
      print(table(validation_metadata$class))
      cat("Selected:", nrow(validation_metadata), "sequences\n\n")
    }

    # Check k-mer compatibility
    n_kmers_class <- ncol(classification_dataset)
    n_kmers_val <- ncol(validation_dataset)

    if (n_kmers_class != n_kmers_val) {
      warning("K-mer count mismatch: classification=", n_kmers_class,
              ", validation=", n_kmers_val, call. = FALSE)

      # Align columns
      common_kmers <- intersect(colnames(classification_dataset), colnames(validation_dataset))

      if (length(common_kmers) < n_kmers_class * 0.9) {
        warning("Only ", length(common_kmers), " common k-mers found (",
                round(length(common_kmers) / n_kmers_class * 100, 1), "%)",
                call. = FALSE)
      }

      validation_dataset <- validation_dataset[, common_kmers, drop = FALSE]
      classification_dataset <- classification_dataset[, common_kmers, drop = FALSE]
    }
  }

  # Save results
  true_labels_classification <- classification_metadata
  save(true_labels_classification, file = "true_labels_classification.RData")

  save(classification_dataset, file = "classification_dataset.RData")

  if (nrow(validation_metadata) > 0) {
    true_labels_validation <- validation_metadata
    save(true_labels_validation, file = "true_labels_validation.RData")
    save(validation_dataset, file = "validation_dataset.RData")
  }

  # Create backup with all data
  backup_data <- list(
    true_labels_classification = classification_metadata,
    classification_dataset = classification_dataset,
    true_labels_validation = validation_metadata,
    validation_dataset = validation_dataset
  )
  save(backup_data, file = "backup_checkpoint1.RData")

  # Return result
  result <- list(
    classification_dataset = classification_dataset,
    validation_dataset = validation_dataset,
    classification_metadata = classification_metadata,
    validation_metadata = validation_metadata,
    # Aliases for backward compatibility
    true_labels_classification = classification_metadata,
    true_labels_validation = validation_metadata
  )

  # Set attributes
  attr(result, "n_classification") <- nrow(classification_metadata)
  attr(result, "n_validation") <- nrow(validation_metadata)
  attr(result, "validation_type") <- if (is.null(external_validation_fasta_dir)) "internal" else "external"
  attr(result, "min_size") <- min_size
  attr(result, "seq_per_class") <- seq_per_class

  # Set class
  class(result) <- c("sample_selection", "list")

  return(result)
}
