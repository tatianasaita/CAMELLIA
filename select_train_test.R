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
#'   Also ignored when \code{external_train_fasta} is supplied.
#' @param k Integer. K-mer size used by \code{create_data()} to build the
#'   frequency matrix (e.g. \code{k = 4} for tetranucleotides). Required in
#'   Mode 3 (fully external). In Modes 1 and 2 it is extracted automatically
#'   from \code{data_result} if not supplied.
#' @param external_train_fasta Character or NULL. Path to a directory of
#'   FASTA files to use as the training set. When supplied together with
#'   \code{external_test_fasta}, both train and test datasets are built
#'   entirely from external files and \code{seq_per_class}, \code{min_size},
#'   and \code{train_prop} are ignored. When supplied alone (without
#'   \code{external_test_fasta}), an error is raised — you must provide
#'   both or neither external directory.
#' @param external_test_fasta Character or NULL. Path to external FASTA
#'   files used as the test set. When supplied alone (without
#'   \code{external_train_fasta}), \code{data_result} is used to build the
#'   training set as usual and only the test set comes from external files
#'   (original behaviour). When supplied together with
#'   \code{external_train_fasta}, both datasets come from external files.
#' @param verbose Logical. If TRUE, prints detailed information. Default is
#'   TRUE.
#'
#' @return A list of class \code{"train_test_selection"} with:
#' \itemize{
#'   \item \code{train_dataset}: Labelled dataset used for training.
#'   \item \code{test_dataset}: Held-out sequences for final evaluation.
#'   \item \code{train_metadata}: Metadata for training sequences.
#'   \item \code{test_metadata}: Metadata for test sequences.
#' }
#'
#' @details
#' Three operating modes depending on the external directory arguments:
#' \enumerate{
#'   \item \strong{Internal split} (default — both external dirs are NULL):
#'     \code{data_result} is split into train / test according to
#'     \code{seq_per_class}, \code{min_size}, and \code{train_prop}.
#'   \item \strong{External test only} (\code{external_test_fasta} supplied,
#'     \code{external_train_fasta} = NULL): training sequences come from
#'     \code{data_result}; the test set is built from the external FASTA
#'     directory via \code{create_data}.
#'   \item \strong{Fully external} (both dirs supplied): both train and test
#'     datasets are built from their respective FASTA directories via
#'     \code{create_data}, using the same \code{k} stored in
#'     \code{data_result}. \code{seq_per_class}, \code{min_size}, and
#'     \code{train_prop} are ignored in this mode.
#' }
#'
#' @note
#' \itemize{
#'   \item Requires \code{.select_by_class}. See \code{internal-functions.R}.
#'   \item S3 methods available. See \code{methods.R}.
#'   \item Supplying \code{external_train_fasta} without
#'     \code{external_test_fasta} raises an error.
#'   \item In Mode 3, \code{k} must be supplied explicitly.
#' }
#'
#' @examples
#' \dontrun{
#' # Mode 1 — Internal split, default 70/30
#' result_train_test <- select_train_test(
#'   data_result   = result,
#'   seq_per_class = 200,
#'   min_size      = 800
#' )
#'
#'
#' # Mode 2 — Internal train + external test
#' result_train_test <- select_train_test(
#'   data_result             = result,
#'   seq_per_class           = 200,
#'   min_size                = 800,
#'   external_test_fasta = "path/to/test/fastas"
#' )
#'
#' # Mode 3 — Fully external (train + test from FASTA dirs)
#' result_train_test <- select_train_test(
#'   k                        = 4,
#'   external_train_fasta = "path/to/train/fastas",
#'   external_test_fasta  = "path/to/test/fastas"
#' )
#' }

# ---------------------------------------------------------------------------
# Internal helper
# ---------------------------------------------------------------------------

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

# ---------------------------------------------------------------------------
# Internal helper — align k-mer columns between two matrices
# ---------------------------------------------------------------------------

#' Align K-mer Columns Between Train and Test Matrices
#'
#' Checks whether the two matrices share the same k-mer columns and, if not,
#' subsets both to their intersection, emitting warnings as appropriate.
#'
#' @param train_dataset Matrix. K-mer counts for the training set.
#' @param test_dataset  Matrix. K-mer counts for the test set.
#'
#' @return A named list with elements \code{train_dataset} and
#'   \code{test_dataset}, both restricted to common columns.
#'
#' @keywords internal
.align_kmer_columns <- function(train_dataset, test_dataset) {
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
    
    train_dataset <- train_dataset[, common_kmers, drop = FALSE]
    test_dataset  <- test_dataset[,  common_kmers, drop = FALSE]
  }
  
  list(train_dataset = train_dataset, test_dataset = test_dataset)
}

#------------------------------------------------------------------------
.class_from_sequence_name <- function(sequence_names) {
  sub("\\..*", "", sequence_names)
}



# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

#' @export
select_train_test <- function(data_result          = NULL,
                               seq_per_class       = NULL,
                               min_size             = 800,
                               train_prop           = 0.7,
                               k                    = NULL,
                               external_train_fasta = NULL,
                               external_test_fasta  = NULL,
                               verbose              = TRUE) {
  
  
  # ---------------------------------------------------------------------------
  # Argument validation
  # ---------------------------------------------------------------------------
  if (!is.numeric(train_prop) || length(train_prop) != 1L ||
      train_prop <= 0 || train_prop >= 1) {
    stop("'train_prop' must be a single numeric value in (0, 1).", call. = FALSE)
  }
  
  if (!is.null(external_train_fasta) && is.null(external_test_fasta)) {
    stop(
      "'external_train_fasta' was supplied but 'external_test_fasta' is NULL.\n",
      "  Please supply both files together (Mode 3), or supply only\n",
      "  'external_test_fasta' to use an internal training set (Mode 2).",
      call. = FALSE
    )
  }
  
  # Determine operating mode
  mode <- if (!is.null(external_train_fasta) && !is.null(external_test_fasta)) {
    "fully_external"       # Mode 3
  } else if (is.null(external_train_fasta) && !is.null(external_test_fasta)) {
    "external_test_only"   # Mode 2
  } else {
    "internal_split"       # Mode 1 (default)
  }
  
  # Resolve k
  if (mode == "fully_external") {
    if (is.null(k)) {
      stop(
        "In fully external mode (Mode 3), 'k' must be supplied explicitly.\n",
        "  Example: k = 4",
        call. = FALSE
      )
    }
    k <- as.integer(k)
  } else {
    if (is.null(data_result)) {
      stop("'data_result' must be supplied in Modes 1 and 2.", call. = FALSE)
    }
    if (is.null(k)) {
      k <- unlist(data_result$k)[[1]]
    }
    k <- as.integer(k)
  }
  
  # ===========================================================================
  # MODE 3 — Fully external: both train and test come from FASTA directories
  # ===========================================================================
  if (mode == "fully_external") {
    
    if (verbose) {
      cat("Mode: Fully External (train + test from FASTA files)\n")
      cat("K =", k, "\n")
      cat("NOTE: 'seq_per_class', 'min_size', and 'train_prop' are ignored",
          "in this mode.\n\n")
    }
    
    # --- Train ---------------------------------------------------------------
    if (verbose) cat("Building TRAIN dataset from:", external_train_fasta, "\n")
    
    train_tmp_dir <- tempfile(pattern = "camellia_train_")
    dir.create(train_tmp_dir)
    on.exit(unlink(train_tmp_dir, recursive = TRUE), add = TRUE)
    
    if (!file.copy(external_train_fasta, train_tmp_dir)) {
      stop("Failed to copy train FASTA file to temporary directory: ",
           external_train_fasta, call. = FALSE)
    }
    
    train_data_result <- tryCatch({
      create_data(input = train_tmp_dir, k = k)
    }, error = function(e) {
      stop("Error creating external train dataset: ", e$message, call. = FALSE)
    })
    
    train_metadata         <- train_data_result$metadata
    train_dataset          <- train_data_result$kmers
    train_metadata$dataset <- "train"
    train_metadata$class   <- .class_from_sequence_name(train_metadata$sequence_name)
    train_dataset$CLASS    <- train_metadata$class                                     
    rownames(train_dataset) <- train_metadata$sequence_name
    
    if (nrow(train_metadata) == 0) {
      stop("No training sequences found in: ", external_train_fasta, call. = FALSE)
    }
    
    if (verbose) {
      cat("Train sequences per class:\n")
      print(table(train_metadata$class))
      cat("Total:", nrow(train_metadata), "sequences\n\n")
    }
    
    # --- Test ----------------------------------------------------------------
    if (verbose) cat("Building TEST dataset from:", external_test_fasta, "\n")
    
    test_tmp_dir <- tempfile(pattern = "camellia_test_")
    dir.create(test_tmp_dir)
    on.exit(unlink(test_tmp_dir, recursive = TRUE), add = TRUE)
    
    if (!file.copy(external_test_fasta, test_tmp_dir)) {
      stop("Failed to copy test FASTA file to temporary directory: ",
           external_test_fasta, call. = FALSE)
    }
    
    test_data_result <- tryCatch({
      create_data(input = test_tmp_dir, k = k)
    }, error = function(e) {
      stop("Error creating external test dataset: ", e$message, call. = FALSE)
    })
    
    test_metadata         <- test_data_result$metadata
    test_dataset          <- test_data_result$kmers
    test_metadata$dataset <- "test"
    test_metadata$class   <- .class_from_sequence_name(test_metadata$sequence_name)
    test_dataset$CLASS    <- test_metadata$class                                      
    rownames(test_dataset) <- test_metadata$sequence_name
    
    if (nrow(test_metadata) == 0) {
      stop("No test sequences found in: ", external_test_fasta, call. = FALSE)
    }
    
    if (verbose) {
      cat("Test sequences per class:\n")
      print(table(test_metadata$class))
      cat("Total:", nrow(test_metadata), "sequences\n\n")
    }
    
    
    
    # --- Align k-mer columns -------------------------------------------------
    aligned       <- .align_kmer_columns(train_dataset, test_dataset)
    train_dataset <- aligned$train_dataset
    test_dataset  <- aligned$test_dataset
    
#    n_train_per_class <- table(train_metadata$class)
#    adjusted          <- FALSE
    
  } else {
   
    # MODES 1 & 2 — training set always built from data_result
    metadata <- data_result$metadata
    kmers    <- data_result$kmers
    
    # --- Filter by minimum sequence length -----------------------------------
    metadata_filtered <- metadata[metadata$length >= min_size, ]
    
    if (nrow(metadata_filtered) == 0) {
      stop("No sequences available with length >= ", min_size, call. = FALSE)
    }
    
    # --- Determine n_train per class -----------------------------------------
    classes     <- unique(metadata_filtered$class)
    class_sizes <- sapply(classes, function(cls) {
      nrow(metadata_filtered[metadata_filtered$class == cls, ])
    })
    names(class_sizes) <- classes
    min_class_size     <- min(class_sizes)
    
    if (is.null(seq_per_class)) {
      # seq_per_class not supplied: apply train_prop split to all classes
      adjusted          <- FALSE
      classes_adjusted  <- character(0)
      n_train_per_class <- pmax(1L, floor(class_sizes * train_prop))
      names(n_train_per_class) <- classes
    } else {
      adjusted         <- any(class_sizes < seq_per_class)
      classes_adjusted <- names(class_sizes[class_sizes < seq_per_class])
      
      n_train_per_class <- ifelse(
        class_sizes >= seq_per_class,
        seq_per_class,
        pmax(1L, floor(class_sizes * train_prop))
      )
      names(n_train_per_class) <- classes
    }
    
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
    
    # --- Select training sequences -------------------------------------------
    train_metadata <- do.call(rbind, lapply(classes, function(cls) {
      cls_meta <- metadata_filtered[metadata_filtered$class == cls, ]
      n_tr <- if (!is.null(seq_per_class)) {
        min(nrow(cls_meta), seq_per_class)  # todas, até o limite definido
      } else {
        nrow(cls_meta)                       # sem limite: todas as sequências
      }
      .select_by_class(cls_meta, n_tr, "train")
    }))
    
    if (nrow(train_metadata) == 0) {
      stop("Could not select any sequences for training.", call. = FALSE)
    }
    
    selected_indices        <- match(train_metadata$sequence_name, metadata$sequence_name)
    train_dataset           <- kmers[selected_indices, , drop = FALSE]
    rownames(train_dataset) <- train_metadata$sequence_name
    
    if (verbose) {
      cat("Total training sequences selected:", nrow(train_metadata), "\n")
      cat("Sequences per class (train):\n")
      print(table(train_metadata$class))
      cat("\n")
    }
    
    # =========================================================================
    # TEST DATASET
    # =========================================================================
    if (mode == "internal_split") {
      
      if (verbose) cat("Mode: Internal Split\n")
      
      train_sequences       <- train_metadata$sequence_name
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
      # mode == "external_test_only"
      if (verbose) cat("Mode: Internal Train + External Test (k =", k, ")\n")
      
      test_tmp_dir <- tempfile(pattern = "camellia_test_")
      dir.create(test_tmp_dir)
      on.exit(unlink(test_tmp_dir, recursive = TRUE), add = TRUE)
      
      if (!file.copy(external_test_fasta, test_tmp_dir)) {
        stop("Failed to copy test FASTA file to temporary directory: ",
             external_test_fasta, call. = FALSE)
      }
      
      external_data_result <- tryCatch({
        create_data(input = test_tmp_dir, k = k)
      }, error = function(e) {
        stop("Error creating external test dataset: ", e$message, call. = FALSE)
      })
      
      test_metadata              <- external_data_result$metadata
      test_dataset               <- external_data_result$kmers
      test_metadata$dataset      <- "test"
      test_metadata$class        <- .class_from_sequence_name(test_metadata$sequence_name)
      test_dataset$CLASS         <- test_metadata$class
      rownames(test_dataset)     <- test_metadata$sequence_name
      
      if (nrow(test_metadata) == 0) {
        stop("No test sequences found in: ", external_test_fasta, call. = FALSE)
      }
      
      if (verbose) {
        cat("Available sequences:\n")
        print(table(test_metadata$class))
        cat("Total:", nrow(test_metadata), "sequences\n\n")
      }
      
      aligned       <- .align_kmer_columns(train_dataset, test_dataset)
      train_dataset <- aligned$train_dataset
      test_dataset  <- aligned$test_dataset
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
  
  if (mode != "internal_split") {
    n_train_per_class <- table(train_metadata$class)
    adjusted          <- FALSE
  }
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
  attr(result, "mode")          <- mode
  attr(result, "min_size")      <- min_size
  attr(result, "seq_per_class") <- seq_per_class
  attr(result, "train_prop")    <- train_prop
  attr(result, "n_train_class") <- n_train_per_class
  attr(result, "adjusted")      <- adjusted
  attr(result, "k")             <- k
  
  class(result) <- "train_test_selection"
  return(result)
}

