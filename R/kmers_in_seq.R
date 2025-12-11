#' Find K-mer Motifs in Training and Validation Sequences
#'
#' @param motifs Character vector or list of motifs to search for
#' @param cluster_result Object from cluster_dendrogram() containing metadata
#' @param sequences Character vector of sequences (default: NULL, reads from input_dir)
#' @param input_dir Path to training FASTA directory (default: NULL)
#' @param sequence_names Character vector of sequence names (default: NULL)
#' @param external_validation_fasta_dir Path to validation FASTA directory (default: NULL)
#' @param validation_predictions Data.frame with 'class' column for validation (default: NULL)
#' @param verbose Logical (default: TRUE)
#'
#' @return Data.frame with columns: motif, sequence_name, class, sequence_length, position_start, position_end, dataset
#' @export
kmers_in_seq <- function(motifs,
                         cluster_result,
                         sequences = NULL,
                         input_dir = NULL,
                         sequence_names = NULL,
                         external_validation_fasta_dir = NULL,
                         validation_predictions = NULL,
                         verbose = TRUE) {

  # ===== VALIDATION =====
  if (is.null(motifs) || length(motifs) == 0) {
    stop("'motifs' must be a non-empty character vector or list", call. = FALSE)
  }

  if (is.null(cluster_result) || !is.list(cluster_result)) {
    stop("'cluster_result' must be provided and must be a list", call. = FALSE)
  }

  if (is.null(cluster_result$data_result) ||
      is.null(cluster_result$data_result$metadata)) {
    stop("'cluster_result' must contain 'data_result$metadata'", call. = FALSE)
  }

  metadata <- cluster_result$data_result$metadata
  if (!is.data.frame(metadata) ||
      !all(c("sequence_name", "class", "length") %in% colnames(metadata))) {
    stop("'cluster_result$data_result$metadata' must be a data.frame with columns: sequence_name, class, length",
         call. = FALSE)
  }

  if (is.null(sequences) && is.null(input_dir)) {
    stop("Either 'sequences' or 'input_dir' must be provided", call. = FALSE)
  }

  # ===== SETUP n_cores (FORCE 1 DURING TESTING) =====
  # Detect if running in testthat environment
  is_testing <- identical(Sys.getenv("TESTTHAT"), "true") ||
    !is.na(Sys.getenv("TESTTHAT_PKG", unset = NA)) ||
    identical(Sys.getenv("NOT_CRAN"), "true")

  n_cores <- if (is_testing) {
    1  # Force sequential processing during tests
  } else if (requireNamespace("parallel", quietly = TRUE)) {
    max(1, parallel::detectCores() - 1)
  } else {
    1
  }

  # Process motifs
  if (is.list(motifs) && !is.data.frame(motifs)) {
    motifs <- unlist(motifs, use.names = FALSE)
  }
  unique_motifs <- unique(motifs)

  # Get training class and length lookup
  train_class_lookup <- setNames(
    as.character(cluster_result$data_result$metadata$class),
    cluster_result$data_result$metadata$sequence_name
  )

  train_length_lookup <- setNames(
    as.integer(cluster_result$data_result$metadata$length),
    cluster_result$data_result$metadata$sequence_name
  )

  # Check if has external validation
  has_validation <- !is.null(external_validation_fasta_dir)

  # Read training sequences
  if (is.null(sequences)) {
    if (verbose) message(sprintf("Reading TRAINING sequences from: %s", input_dir))
    fasta_data <- .read_fasta_sequences(input_dir, verbose)
    sequences <- fasta_data$sequences
    sequence_names <- fasta_data$names
  }

  if (is.null(sequence_names)) {
    sequence_names <- if (!is.null(names(sequences))) names(sequences) else sprintf("seq_%d", seq_along(sequences))
  }

  n_train <- length(sequences)
  n_val <- 0

  # Print summary
  if (verbose) {
    message("\n=== K-mer Search in Sequences ===")
    message(sprintf("TRAINING sequences: %d", n_train))
  }

  # Read validation sequences if provided
  val_sequences <- NULL
  val_sequence_names <- NULL
  val_class_lookup <- NULL
  val_length_lookup <- NULL

  if (has_validation) {
    if (verbose) message(sprintf("Reading VALIDATION sequences from: %s", external_validation_fasta_dir))
    val_fasta_data <- .read_fasta_sequences(external_validation_fasta_dir, verbose)
    val_sequences <- val_fasta_data$sequences
    val_sequence_names <- val_fasta_data$names
    n_val <- length(val_sequences)

    # Get validation class lookup
    if (!is.null(validation_predictions) && "class" %in% colnames(validation_predictions)) {
      val_class_lookup <- setNames(
        as.character(validation_predictions$class),
        val_sequence_names
      )
    }

    # Calculate validation sequence lengths
    val_length_lookup <- setNames(
      sapply(val_sequences, nchar),
      val_sequence_names
    )

    if (verbose) message(sprintf("VALIDATION sequences: %d", n_val))
  }

  if (verbose) {
    message(sprintf("Motifs to search: %d | Cores: %d\n", length(unique_motifs), n_cores))
  }

  start_time <- Sys.time()

  # Search in training
  if (verbose) message("Searching motifs in TRAINING sequences...")
  train_result <- .find_motifs_parallel(
    sequences, unique_motifs, sequence_names, train_class_lookup, train_length_lookup,
    n_cores, verbose, if (has_validation) "training" else NULL
  )

  # Search in validation
  val_result <- NULL
  if (has_validation) {
    if (verbose) message("\nSearching motifs in VALIDATION sequences...")
    val_result <- .find_motifs_parallel(
      val_sequences, unique_motifs, val_sequence_names, val_class_lookup, val_length_lookup,
      n_cores, verbose, "validation"
    )
  }

  # Combine results
  result_df <- if (!is.null(val_result)) rbind(train_result, val_result) else train_result
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Add attributes
  attr(result_df, "n_train_sequences") <- n_train
  attr(result_df, "n_validation_sequences") <- n_val
  attr(result_df, "n_motifs") <- length(unique_motifs)
  attr(result_df, "n_occurrences") <- nrow(result_df)
  attr(result_df, "elapsed_time") <- elapsed
  attr(result_df, "has_validation") <- has_validation
  class(result_df) <- c("kmers_in_seq_result", "data.frame")

  # Print summary
  if (verbose) {
    message(sprintf("\n=== Complete ==="))
    message(sprintf("Total occurrences: %d | Time: %.2f s", nrow(result_df), elapsed))

    if (nrow(result_df) > 0 && has_validation && "dataset" %in% colnames(result_df)) {
      message("\nOccurrences by dataset:")
      print(table(result_df$dataset))
    }
  }

  return(result_df)
}


# ===== INTERNAL FUNCTIONS =====

#' @keywords internal
.read_fasta_sequences <- function(input_dir, verbose) {
  fasta_files <- list.files(input_dir, pattern = "\\.(fasta|fa|fna|faa)$",
                            full.names = TRUE, ignore.case = TRUE)

  if (verbose) message(sprintf("  Found %d FASTA file(s)", length(fasta_files)))

  all_sequences <- list()
  all_names <- character()

  for (fasta_file in fasta_files) {
    if (verbose) message(sprintf("  Reading: %s", basename(fasta_file)))
    seqs <- seqinr::read.fasta(fasta_file, seqtype = "DNA", as.string = TRUE,
                               forceDNAtolower = FALSE)
    for (seq_name in names(seqs)) {
      all_sequences[[seq_name]] <- toupper(seqs[[seq_name]])
      all_names <- c(all_names, seq_name)
    }
  }

  if (verbose) message(sprintf("  Total sequences read: %d", length(all_sequences)))

  return(list(sequences = unlist(all_sequences, use.names = FALSE), names = all_names))
}


#' @keywords internal
.find_motifs_parallel <- function(sequences, motifs, sequence_names,
                                  class_lookup, length_lookup, n_cores, verbose, dataset_label = NULL) {

  process_motif <- function(motif) {
    results_list <- list()
    counter <- 0

    for (j in seq_along(sequences)) {
      positions <- stringi::stri_locate_all_fixed(sequences[j], motif)[[1]]

      if (!is.na(positions[1, 1])) {
        counter <- counter + 1
        seq_name <- sequence_names[j]
        seq_class <- if (!is.null(class_lookup)) {
          cls <- class_lookup[seq_name]
          if (is.na(cls)) "Unknown" else as.character(cls)
        } else {
          "Unknown"
        }

        seq_length <- if (!is.null(length_lookup)) {
          len <- length_lookup[seq_name]
          if (is.na(len)) NA_integer_ else as.integer(len)
        } else {
          NA_integer_
        }

        temp_df <- data.frame(
          motif = rep(motif, nrow(positions)),
          sequence_name = rep(seq_name, nrow(positions)),
          class = rep(seq_class, nrow(positions)),
          position_start = as.integer(positions[, 1]),
          position_end = as.integer(positions[, 2]),
          sequence_length = rep(seq_length, nrow(positions)),
          stringsAsFactors = FALSE
        )

        if (!is.null(dataset_label)) {
          temp_df$dataset <- rep(dataset_label, nrow(positions))
        }

        results_list[[counter]] <- temp_df
      }
    }

    if (counter > 0) do.call(rbind, results_list) else NULL
  }

  # Parallel or sequential
  if (n_cores > 1) {
    if (verbose) message(sprintf("  Using %d cores...", n_cores))
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("sequences", "sequence_names", "class_lookup", "length_lookup", "dataset_label"),
                            envir = environment())
    parallel::clusterEvalQ(cl, library(stringi))
    results_list <- parallel::parLapply(cl, motifs, process_motif)
  } else {
    results_list <- lapply(motifs, process_motif)
  }

  # Combine
  results_list <- results_list[!sapply(results_list, is.null)]

  if (length(results_list) > 0) {
    result_df <- do.call(rbind, results_list)
    rownames(result_df) <- NULL
  } else {
    result_df <- data.frame(
      motif = character(0),
      sequence_name = character(0),
      class = character(0),
      position_start = integer(0),
      position_end = integer(0),
      sequence_length = integer(0),
      stringsAsFactors = FALSE
    )
    if (!is.null(dataset_label)) result_df$dataset <- character(0)
  }

  if (verbose) message(sprintf("  Found %d occurrences", nrow(result_df)))

  return(result_df)
}


# ===== S3 METHODS =====

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


#' @export
summary.kmers_in_seq_result <- function(object, ...) {
  cat("\n=== K-mer Search Summary ===\n\n")
  print(object)

  if (nrow(object) > 0) {
    if ("dataset" %in% names(object)) {
      cat("\n\nBy dataset:\n")
      print(table(object$dataset))

      if ("class" %in% names(object)) {
        cat("\n\nBy class and dataset:\n")
        print(table(object$class, object$dataset))
      }
    }

    cat("\n\nTop 15 motifs:\n")
    print(head(sort(table(object$motif), decreasing = TRUE), 15))

    if ("class" %in% names(object)) {
      cat("\n\nBy class:\n")
      print(sort(table(object$class), decreasing = TRUE))
    }

    cat("\n\nBy sequence (summary):\n")
    by_seq <- table(object$sequence_name)
    cat(sprintf("  Min: %d | Median: %d | Mean: %.1f | Max: %d\n",
                min(by_seq), median(by_seq), mean(by_seq), max(by_seq)))
  }

  invisible(object)
}
