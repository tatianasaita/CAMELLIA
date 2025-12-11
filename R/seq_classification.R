#' Select Sample for Training and Validation Datasets
#'
#' @param cluster_result List from cluster_dendrogram() containing data_result.
#' @param k_per_class Integer. Number of sequences per class for training.
#' @param min_size Integer. Minimum sequence length for training (not applied to validation).
#' @param external_validation_fasta_dir Character or NULL. Path to external FASTA files.
#' @param k_for_external_validation Integer. K-mer size for external validation. Default 6.
#'
#' @return A list of class "sample_selection" with training and validation datasets.
#'
#' @export
select_sample_train_validation <- function(cluster_result,
                                           k_per_class,
                                           min_size = 800,
                                           external_validation_fasta_dir = NULL,
                                           k_for_external_validation = 6L) {

  # ===== VALIDATION =====

  # Validate cluster_result structure
  if (!is.list(cluster_result) || !("data_result" %in% names(cluster_result))) {
    stop("'cluster_result' must be a list from cluster_dendrogram() with 'data_result'", call. = FALSE)
  }

  data_result <- cluster_result$data_result

  if (!is.list(data_result)) {
    stop("'data_result' must be a list", call. = FALSE)
  }

  if (!all(c("metadata", "kmers") %in% names(data_result))) {
    stop("'data_result' must contain 'metadata' and 'kmers'", call. = FALSE)
  }

  metadata <- data_result$metadata
  kmers <- data_result$kmers

  if (is.null(metadata) || !is.data.frame(metadata)) {
    stop("'metadata' must be a data.frame", call. = FALSE)
  }

  if (is.null(kmers) || !is.data.frame(kmers)) {
    stop("'kmers' must be a data.frame", call. = FALSE)
  }

  if (nrow(metadata) == 0) {
    stop("'metadata' is empty", call. = FALSE)
  }

  if (nrow(kmers) == 0) {
    stop("'kmers' is empty", call. = FALSE)
  }

  if (nrow(metadata) != nrow(kmers)) {
    stop("'metadata' and 'kmers' must have the same number of rows", call. = FALSE)
  }

  required_cols <- c("sequence_name", "length", "class")
  if (!all(required_cols %in% colnames(metadata))) {
    stop("'metadata' missing required columns: ",
         paste(setdiff(required_cols, colnames(metadata)), collapse = ", "), call. = FALSE)
  }

  # Validate k_per_class
  if (!is.numeric(k_per_class)) {
    stop("'k_per_class' must be numeric", call. = FALSE)
  }
  if (length(k_per_class) != 1L) {
    stop("'k_per_class' must be a single value", call. = FALSE)
  }
  if (is.na(k_per_class)) {
    stop("'k_per_class' cannot be NA", call. = FALSE)
  }
  if (k_per_class <= 0) {
    stop("'k_per_class' must be positive", call. = FALSE)
  }
  if ((k_per_class %% 1) != 0) {
    stop("'k_per_class' must be an integer", call. = FALSE)
  }

  # Validate min_size
  if (!is.numeric(min_size)) {
    stop("'min_size' must be numeric", call. = FALSE)
  }
  if (length(min_size) != 1L) {
    stop("'min_size' must be a single value", call. = FALSE)
  }
  if (is.na(min_size)) {
    stop("'min_size' cannot be NA", call. = FALSE)
  }
  if (min_size <= 0) {
    stop("'min_size' must be positive", call. = FALSE)
  }
  if ((min_size %% 1) != 0) {
    stop("'min_size' must be an integer", call. = FALSE)
  }

  # Validate external validation parameters
  if (!is.null(external_validation_fasta_dir)) {
    if (!is.character(external_validation_fasta_dir)) {
      stop("'external_validation_fasta_dir' must be a character string or NULL", call. = FALSE)
    }

    invalid_dirs <- character(0)
    for (dir_path in external_validation_fasta_dir) {
      if (!dir.exists(dir_path)) {
        invalid_dirs <- c(invalid_dirs, dir_path)
      }
    }

    if (length(invalid_dirs) > 0) {
      stop("External validation directory does not exist: ",
           paste(invalid_dirs, collapse = ", "), call. = FALSE)
    }

    if (!is.numeric(k_for_external_validation)) {
      stop("'k_for_external_validation' must be numeric", call. = FALSE)
    }
    if (length(k_for_external_validation) != 1L) {
      stop("'k_for_external_validation' must be a single value", call. = FALSE)
    }
    if (is.na(k_for_external_validation)) {
      stop("'k_for_external_validation' cannot be NA", call. = FALSE)
    }
    if (k_for_external_validation <= 0) {
      stop("'k_for_external_validation' must be positive", call. = FALSE)
    }
    if ((k_for_external_validation %% 1) != 0) {
      stop("'k_for_external_validation' must be an integer", call. = FALSE)
    }
  }

  # ===== CLASSIFICATION DATASET =====

  cat("\n=== CLASSIFICATION DATASET ===\n")

  # Filter by min_size
  metadata_filtered <- metadata[metadata$length >= min_size, ]

  if (nrow(metadata_filtered) == 0) {
    stop("No sequences available with length >= ", min_size, call. = FALSE)
  }

  cat("Sequences with length >=", min_size, ":", nrow(metadata_filtered), "\n")
  cat("Available classes:", paste(unique(metadata_filtered$class), collapse = ", "), "\n")

  # Select sequences per class
  classification_metadata <- .select_by_class(metadata_filtered, k_per_class, "classification")

  if (nrow(classification_metadata) == 0) {
    stop("Could not select any sequences for classification", call. = FALSE)
  }

  # Extract k-mers using sequence names
  selected_indices <- match(classification_metadata$sequence_name, metadata$sequence_name)
  classification_dataset <- kmers[selected_indices, , drop = FALSE]
  rownames(classification_dataset) <- classification_metadata$sequence_name

  cat("Total classification sequences:", nrow(classification_metadata), "\n\n")

  # ===== VALIDATION DATASET =====

  if (is.null(external_validation_fasta_dir)) {
    # Internal validation
    cat("=== INTERNAL VALIDATION DATASET ===\n")

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

      cat("Validation sequences:", nrow(validation_metadata), "\n")
      cat("Validation classes:", paste(unique(validation_metadata$class), collapse = ", "), "\n\n")

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
    cat("=== EXTERNAL VALIDATION DATASET ===\n")

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

    cat("Available sequences by class:\n")
    print(table(validation_metadata$class))
    cat("Selected:", nrow(validation_metadata), "sequences\n\n")

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

  # ===== SAVE RESULTS =====

  cat("=== SAVING RESULTS ===\n")

  # Save with the expected variable names
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

  cat("Files saved successfully\n\n")

  # ===== RETURN =====

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
  attr(result, "k_per_class") <- k_per_class

  # Set class
  class(result) <- c("sample_selection", "list")

  return(result)
}


# ===== INTERNAL FUNCTIONS =====

#' Select k sequences per class
#' @keywords internal
.select_by_class <- function(data, k, dataset_name = "") {
  classes <- unique(data$class)

  sequences_selected <- lapply(classes, function(cls) {
    class_data <- data[data$class == cls, ]
    n_available <- nrow(class_data)
    k_adjusted <- min(k, n_available)

    if (k_adjusted == 0) {
      return(NULL)
    }

    cat("   ", cls, ":", k_adjusted, "selected for", dataset_name, ",",
        n_available - k_adjusted, "remaining\n")

    set.seed(123)  # For reproducibility in tests
    selected_indices <- sample(nrow(class_data), k_adjusted)
    class_data[selected_indices, ]
  })

  sequences_selected <- do.call(rbind, sequences_selected[!sapply(sequences_selected, is.null)])
  rownames(sequences_selected) <- NULL

  sequences_selected
}


# ===== S3 METHODS =====

#' @export
print.sample_selection <- function(x, ...) {
  cat("\n=== Sample Selection Summary ===\n\n")
  cat("Classification sequences:", attr(x, "n_classification"), "\n")
  cat("Validation sequences:    ", attr(x, "n_validation"), "\n")
  cat("Validation type:         ", attr(x, "validation_type"), "\n")
  cat("Min sequence size:       ", attr(x, "min_size"), "\n")
  cat("Sequences per class:     ", attr(x, "k_per_class"), "\n\n")

  cat("Class distribution (classification):\n")
  print(table(x$classification_metadata$class))

  if (attr(x, "n_validation") > 0) {
    cat("\nClass distribution (validation):\n")
    print(table(x$validation_metadata$class))
  }

  cat("\n")
  invisible(x)
}


#' @export
summary.sample_selection <- function(object, ...) {
  cat("\n=== Detailed Sample Selection Summary ===\n\n")

  cat("Classification Dataset:\n")
  cat("  Sequences:", attr(object, "n_classification"), "\n")
  cat("  Features:", ncol(object$classification_dataset), "\n")
  cat("  Classes:\n")
  print(table(object$classification_metadata$class))
  cat("\n")

  if (attr(object, "n_validation") > 0) {
    cat("Validation Dataset (", attr(object, "validation_type"), "):\n", sep = "")
    cat("  Sequences:", attr(object, "n_validation"), "\n")
    cat("  Features:", ncol(object$validation_dataset), "\n")
    cat("  Classes:\n")
    print(table(object$validation_metadata$class))
    cat("\n")
  } else {
    cat("Validation Dataset: None\n\n")
  }

  cat("Sequence Length Statistics:\n")
  cat("  Classification:\n")
  print(summary(object$classification_metadata$length))

  if (attr(object, "n_validation") > 0) {
    cat("  Validation:\n")
    print(summary(object$validation_metadata$length))
  }

  cat("\n")
  invisible(object)
}
