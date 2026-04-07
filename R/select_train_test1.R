#' Select Sample for Training and Validation Datasets
#'
#' Splits data into training and validation sets.
#'
#' @param data_result List from create_data() containing metadata and kmers.
#' @param seq_per_class Integer. Desired number of sequences per class for training.
#'   If the smallest class (after length filtering) has fewer sequences, a
#'   \code{train_prop}/\code{(1-train_prop)} split is applied to that class and
#'   the resulting n_train is used for all classes.
#' @param min_size Integer. Minimum sequence length (bp) to be eligible for
#'   training. Sequences below this threshold go directly to the test set.
#' @param train_prop Numeric in (0, 1). Proportion of sequences to use for
#'   training when the smallest class has fewer sequences than
#'   \code{seq_per_class} (default: 0.7, i.e. 70\% train / 30\% test).
#'   Ignored when all classes have at least \code{seq_per_class} sequences.
#' @param external_test_fasta_dir Character or NULL. Path to external FASTA
#'   files.
#' @param verbose Logical. If TRUE, prints detailed information. Default is
#'   TRUE.
#'
#' @return A list of class \code{"train_test_selection"} with:
#' \itemize{
#'   \item \code{train_dataset}: Labelled dataset used for training.
#'   \item \code{test_dataset}: Held-out (internal) or external sequences for
#'     final evaluation.
#'   \item \code{train_metadata}: Metadata for training sequences.
#'   \item \code{test_metadata}: Metadata for test sequences.
#' }
#'
#' @details
#' Splits data according to ML convention — train / validation / test:
#' \itemize{
#'   \item \strong{Training set} (\code{train_prop}, default 70\%): Used to
#'     fit the model.
#'   \item \strong{Validation set} (remainder of training split): Used to
#'     evaluate the model during development; created inside
#'     \code{train_models_rf_xgboost} via \code{createDataPartition}.
#'   \item \strong{Test set}: Completely independent data, never seen during
#'     training; either the remaining labelled sequences (internal) or an
#'     external unlabelled FASTA directory.
#' }
#'
#' @note
#' \itemize{
#'   \item Requires \code{.select_by_class}. See \code{internal-functions.R}.
#'   \item S3 methods available. See \code{methods.R}.
#' }
#'
#' @examples
#' \dontrun{
#' # Internal test set — default 70/30 split
#' result_train_test <- select_train_test(
#'   data_result   = result,
#'   seq_per_class = 200,
#'   min_size      = 800
#' )
#'
#' # Increase training proportion to 80 %
#' result_train_test <- select_train_test(
#'   data_result   = result,
#'   seq_per_class = 200,
#'   min_size      = 800,
#'   train_prop    = 0.8
#' )
#'
#' # External test set (unlabelled)
#' result_train_test <- select_train_test(
#'   data_result             = result,
#'   seq_per_class           = 200,
#'   min_size                = 800,
#'   external_test_fasta_dir = "path/to/test/fastas"
#' )
#' }
#'
#' Select Sequences for Training (Internal)
#'
#' Randomly samples \code{n_train} rows from \code{data}, which must already
#' be filtered to a single class. The seed is intentionally not set here —
#' reproducibility is the caller's responsibility (set seed before calling
#' \code{select_train_test}).
#'
#' @param data Data frame of metadata for a single class.
#' @param n_train Integer. Number of sequences to select.
#' @param dataset_name Character label written to the \code{dataset} column
#'   (default: \code{"train"}).
#'
#' @return A data frame with \code{n_train} rows (or fewer if \code{data} has
#'   fewer rows) and a \code{dataset} column set to \code{dataset_name}.
#'
#' @keywords internal
.select_by_class <- function(data, n_train, dataset_name = "train") {
  n_available        <- nrow(data)
  n_select           <- min(as.integer(n_train), n_available)
  selected           <- data[sample(n_available, n_select), ]
  selected$dataset   <- dataset_name
  rownames(selected) <- NULL
  selected
}

#' @export
select_train_test1 <- function(data_result,
                               seq_per_class,
                               min_size                = 800,
                               train_prop              = 0.7,
                               external_test_fasta_dir = NULL,
                               verbose                 = TRUE) {

  # Validate train_prop
  if (!is.numeric(train_prop) || length(train_prop) != 1L ||
      train_prop <= 0 || train_prop >= 1) {
    stop("'train_prop' must be a single numeric value in (0, 1).", call. = FALSE)
  }

  # Extract data from data_result
  metadata <- data_result$metadata
  kmers    <- data_result$kmers
  k        <- data_result$k

  # Filter by minimum sequence length
  metadata_filtered <- metadata[metadata$length >= min_size, ]

  if (nrow(metadata_filtered) == 0) {
    stop("No sequences available with length >= ", min_size, call. = FALSE)
  }

  # Determine n_train per class:
  #   - classes with class_size >= seq_per_class -> use seq_per_class
  #   - classes with class_size <  seq_per_class -> use floor(class_size * train_prop)
  classes     <- unique(metadata_filtered$class)
  class_sizes <- sapply(classes, function(cls) {
    nrow(metadata_filtered[metadata_filtered$class == cls, ])
  })
  names(class_sizes) <- classes
  min_class_size     <- min(class_sizes)

  adjusted         <- any(class_sizes < seq_per_class)
  classes_adjusted <- names(class_sizes[class_sizes < seq_per_class])

  # Per-class n_train vector
  n_train_per_class <- ifelse(
    class_sizes >= seq_per_class,
    seq_per_class,
    pmax(1L, floor(class_sizes * train_prop))
  )
  names(n_train_per_class) <- classes

  if (verbose) {
    cat("Sequences with length >=", min_size, ":\n")
    print(class_sizes)
    cat("Smallest class size :", min_class_size, "\n")

    if (adjusted) {
      cat(sprintf(
        "\nWARNING: %d class(es) have fewer sequences than seq_per_class (%d):\n",
        length(classes_adjusted), seq_per_class
      ))
      for (cls in classes_adjusted) {
        n_tr <- n_train_per_class[cls]
        n_te <- class_sizes[cls] - n_tr
        cat(sprintf(
          "  '%s': %d available -> n_train = %d, n_test = %d (%.0f/%.0f split)\n",
          cls, class_sizes[cls], n_tr, n_te,
          train_prop * 100, (1 - train_prop) * 100
        ))
      }
    }

    cat("\nn_train per class:\n")
    print(n_train_per_class)
    cat("\n")
  }

  # Select sequences per class using individual n_train values
  train_metadata <- do.call(rbind, lapply(classes, function(cls) {
    cls_meta <- metadata_filtered[metadata_filtered$class == cls, ]
    n_tr     <- n_train_per_class[cls]
    .select_by_class(cls_meta, n_tr, "train")
  }))

  if (nrow(train_metadata) == 0) {
    stop("Could not select any sequences for training.", call. = FALSE)
  }

  # Extract k-mers for training set
  selected_indices        <- match(train_metadata$sequence_name, metadata$sequence_name)
  train_dataset           <- kmers[selected_indices, , drop = FALSE]
  rownames(train_dataset) <- train_metadata$sequence_name

  if (verbose) {
    cat("Total training sequences selected:", nrow(train_metadata), "\n")
    cat("Sequences per class (train):\n")
    print(table(train_metadata$class))
    cat("\n")
  }

  # ---------------------------------------------------------------------------
  # Test dataset
  # ---------------------------------------------------------------------------
  if (is.null(external_test_fasta_dir)) {

    if (verbose) cat("Internal Test Dataset\n")

    train_sequences <- train_metadata$sequence_name

    # All sequences not used in training go to test (including those below min_size)
    test_metadata         <- metadata[!(metadata$sequence_name %in% train_sequences), ]
    test_metadata$dataset <- "test"

    if (nrow(test_metadata) == 0) {
      warning("No sequences available for internal test set.", call. = FALSE)
      test_metadata <- data.frame(
        sequence_name    = character(0),
        length           = integer(0),
        class            = character(0),
        dataset          = character(0),
        stringsAsFactors = FALSE
      )
      test_dataset <- kmers[0, , drop = FALSE]
    } else {
      if (verbose) {
        cat("Test sequences :", nrow(test_metadata), "\n")
        cat("Test classes   :", paste(unique(test_metadata$class), collapse = ", "), "\n\n")
      }

      test_indices           <- match(test_metadata$sequence_name, metadata$sequence_name)
      test_dataset           <- kmers[test_indices, , drop = FALSE]
      rownames(test_dataset) <- test_metadata$sequence_name

      overlap <- intersect(train_sequences, test_metadata$sequence_name)
      if (length(overlap) > 0) {
        warning("Data leakage detected: ", length(overlap),
                " sequences appear in both train and test sets.",
                call. = FALSE)
      }
    }

  } else {

    # External test set: reuse same k from data_result
    if (verbose) cat("External Test Dataset (k =", k, ")\n")

    test_dir <- external_test_fasta_dir[1]

    external_data_result <- tryCatch({
      create_data(input = test_dir, k = k)
    }, error = function(e) {
      stop("Error creating external test dataset: ", e$message, call. = FALSE)
    })

    test_metadata <- external_data_result$metadata
    test_dataset  <- external_data_result$kmers

    if (nrow(test_metadata) == 0) {
      stop("No test sequences found in: ", test_dir, call. = FALSE)
    }

    if (verbose) {
      cat("Available sequences:\n")
      print(table(test_metadata$class))
      cat("Total:", nrow(test_metadata), "sequences\n\n")
    }

    # Check and align k-mer columns
    n_kmers_train <- ncol(train_dataset)
    n_kmers_test  <- ncol(test_dataset)

    if (n_kmers_train != n_kmers_test) {
      warning("K-mer count mismatch: train=", n_kmers_train,
              ", test=", n_kmers_test, call. = FALSE)

      common_kmers <- intersect(colnames(train_dataset), colnames(test_dataset))

      if (length(common_kmers) < n_kmers_train * 0.7) {
        warning("Only ", length(common_kmers), " common k-mers found (",
                round(length(common_kmers) / n_kmers_train * 100, 1), "%)",
                call. = FALSE)
      }

      test_dataset  <- test_dataset[,  common_kmers, drop = FALSE]
      train_dataset <- train_dataset[, common_kmers, drop = FALSE]
    }
  }

  # ---------------------------------------------------------------------------
  # Save checkpoints
  # ---------------------------------------------------------------------------
  save(train_metadata, file = "train_metadata.RData")
  save(train_dataset,  file = "train_dataset.RData")

  if (nrow(test_metadata) > 0) {
    save(test_metadata, file = "test_metadata.RData")
    save(test_dataset,  file = "test_dataset.RData")
  }

  backup_data <- list(
    train_metadata = train_metadata,
    train_dataset  = train_dataset,
    test_metadata  = test_metadata,
    test_dataset   = test_dataset
  )
  save(backup_data, file = "backup_checkpoint1.RData")

  # ---------------------------------------------------------------------------
  # Return
  # ---------------------------------------------------------------------------
  result <- list(
    train_dataset  = train_dataset,
    test_dataset   = test_dataset,
    train_metadata = train_metadata,
    test_metadata  = test_metadata
  )

  attr(result, "n_train")       <- nrow(train_metadata)
  attr(result, "n_test")        <- nrow(test_metadata)
  attr(result, "test_type")     <- if (is.null(external_test_fasta_dir)) "internal" else "external"
  attr(result, "min_size")      <- min_size
  attr(result, "seq_per_class") <- seq_per_class
  attr(result, "train_prop")    <- train_prop
  attr(result, "n_train_class") <- n_train_per_class
  attr(result, "adjusted")      <- adjusted
  attr(result, "k")             <- k

  class(result) <- "train_test_selection"
  return(result)
}
