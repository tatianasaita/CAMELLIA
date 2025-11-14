#' Select Training and Validation Samples from Cluster Analysis
#'
#' Selects sequences for classification and validation datasets,
#' ensuring stratification by class.
#'
#' @param cluster_result List containing data_result from cluster_dendrogram().
#'   Must include metadata (with sequence information) and kmers (k-mer data).
#' @param k_per_class Integer. Number of sequences to select per class for classification dataset.
#' @param min_size Integer. Minimum sequence length to consider.
#'
#' @return A list with four elements:
#'   \describe{
#'     \item{true_labels_classification}{Metadata for classification sequences}
#'     \item{true_labels_validation}{Metadata for validation sequences}
#'     \item{classification_dataset}{K-mer values for classification sequences}
#'     \item{validation_dataset}{K-mer values for validation sequences}
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Filters sequences by minimum length
#'   \item Selects k_per_class sequences per class for classification dataset
#'   \item Remaining sequences go to validation dataset
#'   \item Ensures no overlap between datasets
#'   \item Matches sequences with k-mer data using rownames
#'   \item Saves all results to RData files and creates backup
#' }
#'
#' @examples
#' \dontrun{
#' datasets <- select_sample_train_validation(
#'   cluster_result = cluster_result,
#'   k_per_class = 50,
#'   min_size = 300
#' )
#'
#' classification_data <- datasets$classification_dataset
#' validation_data <- datasets$validation_dataset
#' classification_labels <- datasets$true_labels_classification
#' validation_labels <- datasets$true_labels_validation
#' }
#'
#' @export
select_sample_train_validation <- function(cluster_result,
                                           k_per_class,
                                           min_size) {

  # === INPUT VALIDATION ===

  if (!is.list(cluster_result)) {
    stop("'cluster_result' must be a list from cluster_dendrogram()")
  }

  if (!("data_result" %in% names(cluster_result))) {
    stop("'cluster_result' must contain 'data_result'")
  }

  data_result <- cluster_result$data_result

  if (!("metadata" %in% names(data_result))) {
    stop("'data_result' must contain 'metadata'")
  }

  if (!("kmers" %in% names(data_result))) {
    stop("'data_result' must contain 'kmers'")
  }

  metadata <- data_result$metadata
  kmers <- data_result$kmers

  if (!is.data.frame(metadata)) {
    stop("'metadata' must be a data.frame")
  }

  if (!is.data.frame(kmers)) {
    stop("'kmers' must be a data.frame")
  }

  if (!is.numeric(k_per_class) || k_per_class <= 0 || k_per_class != as.integer(k_per_class)) {
    stop("'k_per_class' must be a positive integer")
  }

  if (!is.numeric(min_size) || min_size <= 0 || min_size != as.integer(min_size)) {
    stop("'min_size' must be a positive integer")
  }

  # Check if required columns exist
  required_cols_metadata <- c("sequence_name", "length", "class")
  missing_cols <- setdiff(required_cols_metadata, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop(paste("Missing columns in metadata:", paste(missing_cols, collapse = ", ")))
  }

  # === AUXILIARY FUNCTIONS ===

  # Function to verify sequence availability
  verify_availability <- function(data, size = NULL, title = "") {
    if (!is.null(size)) {
      data_filtered <- data[data$length >= size, ]
      cat(title, "Sequence availability by class (length >=", size, "):\n")
    } else {
      data_filtered <- data
      cat(title, "Sequence availability by class (any length):\n")
    }

    # Count by class
    availability <- table(data_filtered$class)
    print(availability)
    cat("\n")

    return(list(
      data_filtered = data_filtered,
      availability = availability
    ))
  }

  # Function to select sequences by class
  select_by_class <- function(data, k, dataset_name = "", sequences_exclude = NULL) {
    # Exclude already selected sequences, if any
    if (!is.null(sequences_exclude)) {
      data <- data[!(data$sequence_name %in% sequences_exclude), ]
    }

    # Group by class
    data_by_class <- split(data, data$class)

    # Select k sequences from each class
    sequences_selected <- lapply(data_by_class, function(class_data) {
      n_available <- nrow(class_data)
      k_adjusted <- min(k, n_available)

      if (k_adjusted > 0) {
        if (k_adjusted < k) {
          warning(paste(dataset_name, "- Class", unique(class_data$class),
                        "has only", n_available, "sequences available.",
                        "Selecting", k_adjusted, "sequences."))
        }

        # Select random sample
        selected_indices <- sample(nrow(class_data), k_adjusted)
        return(class_data[selected_indices, ])
      } else {
        warning(paste(dataset_name, "- Class", unique(class_data$class),
                      "has no sequences available."))
        return(NULL)
      }
    })

    # Combine results, removing NULLs
    sequences_selected <- do.call(rbind, sequences_selected[!sapply(sequences_selected, is.null)])
    rownames(sequences_selected) <- NULL

    return(sequences_selected)
  }

  # === SELECTION OF CLASSIFICATION DATASET ===

  cat("=== CLASSIFICATION DATASET SELECTION ===\n\n")
  verification_classification <- verify_availability(metadata, min_size)
  data_filtered_classification <- verification_classification$data_filtered

  # Select sequences for classification_dataset
  true_labels_classification <- select_by_class(
    data_filtered_classification,
    k_per_class,
    "Classification"
  )

  # Check if sequences were selected
  if (is.null(true_labels_classification) || nrow(true_labels_classification) == 0) {
    stop("Could not select any sequences for classification_dataset with the specified criteria.")
  }

  # Extract names of selected sequences for classification
  names_classification <- true_labels_classification$sequence_name

  # Select in kmers using rownames
  classification_dataset <- kmers[rownames(kmers) %in% names_classification, ]

  cat("CLASSIFICATION DATASET selected:\n")
  cat("Total sequences selected:", nrow(true_labels_classification), "\n")
  cat("Total sequences found in kmers:", nrow(classification_dataset), "\n\n")

  # === SELECTION OF VALIDATION DATASET ===

  cat("=== VALIDATION DATASET SELECTION ===\n\n")

  # Select all sequences not in classification set
  true_labels_validation <- metadata[!(metadata$sequence_name %in% names_classification), ]

  # Check if sequences are available for validation
  if (nrow(true_labels_validation) == 0) {
    warning("No sequences available for validation set after removing classification sequences.")
    true_labels_validation <- NULL
  } else {
    # Print information about validation set
    cat("Validation - Sequence availability by class:\n")
    print(table(true_labels_validation$class))
    cat("\n")
  }

  # Extract names of selected sequences for validation
  names_validation <- if (!is.null(true_labels_validation)) {
    true_labels_validation$sequence_name
  } else {
    character(0)
  }

  # Select in kmers using rownames
  validation_dataset <- kmers[rownames(kmers) %in% names_validation, ]

  # Check how many sequences were found in the feature dataset
  if (length(names_validation) > 0) {
    sequences_found <- sum(names_validation %in% rownames(kmers))
    cat("Validation - Total sequences found in kmers:", sequences_found, "\n\n")

    if (sequences_found < length(names_validation)) {
      warning(paste("Only", sequences_found, "of", length(names_validation),
                    "validation sequences were found in the feature dataset."))
    }
  } else {
    cat("Validation - No sequences available\n\n")
  }

  # === FINAL SUMMARY ===

  cat("=== SELECTION SUMMARY ===\n\n")

  cat("CLASSIFICATION DATASET:\n")
  cat("Total sequences selected:", nrow(true_labels_classification), "\n")
  cat("Distribution by class:\n")
  print(table(true_labels_classification$class))
  cat("Sequences found in kmers:", nrow(classification_dataset), "\n")
  cat("Percentage found:", round(100 * nrow(classification_dataset) / nrow(true_labels_classification), 2), "%\n\n")

  cat("VALIDATION DATASET:\n")
  if (!is.null(true_labels_validation) && nrow(true_labels_validation) > 0) {
    cat("Total sequences selected:", nrow(true_labels_validation), "\n")
    cat("Distribution by class:\n")
    print(table(true_labels_validation$class))
    cat("Sequences found in kmers:", nrow(validation_dataset), "\n")
    cat("Percentage found:", round(100 * nrow(validation_dataset) / nrow(true_labels_validation), 2), "%\n\n")
  } else {
    cat("No sequences available\n\n")
  }

  # === OVERLAP CHECK ===

  overlap <- intersect(names_classification, names_validation)
  if (length(overlap) > 0) {
    warning(paste("There are", length(overlap), "overlapping sequences between datasets!"))
  } else {
    cat("No overlap between datasets.\n\n")
  }

  # === MISSING SEQUENCES CHECK ===

  names_classification_not_found <- setdiff(names_classification, rownames(kmers))
  names_validation_not_found <- setdiff(names_validation, rownames(kmers))

  if (length(names_classification_not_found) > 0) {
    warning(paste("Classification sequences not found in kmers:", length(names_classification_not_found)))
  } else {
    cat("All classification sequences were found in kmers.\n")
  }

  if (length(names_validation_not_found) > 0) {
    warning(paste("Validation sequences not found in kmers:", length(names_validation_not_found)))
  } else if (length(names_validation) > 0) {
    cat("All validation sequences were found in kmers.\n")
  }

  cat("\n")

  # === SAVE RESULTS ===

  cat("=== SAVING RESULTS ===\n\n")

  save(true_labels_classification, file = "true_labels_classification.RData")
  save(true_labels_validation, file = "true_labels_validation.RData")
  save(classification_dataset, file = "classification_dataset.RData")
  save(validation_dataset, file = "validation_dataset.RData")

  # Save checkpoint
  save.image("backup_checkpoint1.RData")

  cat("Files saved successfully:\n")
  cat("  * true_labels_classification.RData\n")
  cat("  * true_labels_validation.RData\n")
  cat("  * classification_dataset.RData\n")
  cat("  * validation_dataset.RData\n")
  cat("  * backup_checkpoint1.RData\n\n")

  # === FINAL STATISTICS ===

  cat("=== FINAL STATISTICS ===\n\n")
  cat("Total sequences in metadata:", nrow(metadata), "\n")
  cat("Total sequences in kmers:", nrow(kmers), "\n")
  cat("Sequences selected for classification:", nrow(true_labels_classification), "\n")
  cat("Sequences selected for validation:", if (!is.null(true_labels_validation)) nrow(true_labels_validation) else 0, "\n")
  cat("Total sequences used:", nrow(true_labels_classification) + if (!is.null(true_labels_validation)) nrow(true_labels_validation) else 0, "\n")
  cat("Sequences not used:", nrow(metadata) - (nrow(true_labels_classification) + if (!is.null(true_labels_validation)) nrow(true_labels_validation) else 0), "\n\n")

  # Return objects for immediate use
  return(list(
    true_labels_classification = true_labels_classification,
    true_labels_validation = true_labels_validation,
    classification_dataset = classification_dataset,
    validation_dataset = validation_dataset
  ))
}
