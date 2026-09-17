#' Find K-mer Motifs in Training and Validation Sequences
#'
#' @param motifs Character vector or list of motifs to search for
#' @param data_result Object of class \code{"kmer_data"} from \code{create_data()},
#'   containing \code{$metadata} with class and length for all labeled sequences.
#' @param labeled_sequences Path to directory containing training FASTA files.
#'   The function will read all FASTA files from this directory.
#' @param external_test_fasta_dir Path to directory containing test FASTA files (default: NULL)
#' @param test_predictions Data.frame with 'class' column for test set (default: NULL)
#' @param verbose Logical (default: TRUE)
#'
#' @return Data.frame with columns: motif, sequence_name, class, sequence_length,
#'   position_start, position_end, dataset
#'
#' @details
#' \itemize{
#'  \item Searches for exact matches of k-mer motifs in DNA sequences
#'  \item Uses parallel processing (automatically detects available CPU cores)
#'  \item Supports both training and test datasets
#'  \item Training sequences require \code{data_result} for metadata (class and length)
#'  \item If a sequence is not found in metadata, class and length will be \code{NA}
#'  \item Test dataset is optional and requires external FASTA files
#'  \item Returns positions (start/end) of all motif occurrences
#' }
#'
#' @note
#' \itemize{
#'   \item Requires .find_motifs_parallel and .read_fasta_sequences. See \code{internal-functions.R} for details.
#'   \item S3 methods available. See \code{methods.R} for details.
#' }
#'
#' @examples
#' \dontrun{
#' # Using directory paths for training sequences
#' result_kmers_in_seq <- kmers_in_seq(
#'   motifs            = motifs,
#'   data_result       = result_create_data,
#'   labeled_sequences = "path/to/fasta/directory"
#' )
#'
#' # With external test dataset
#' result_kmers_in_seq <- kmers_in_seq(
#'   motifs                  = motifs,
#'   data_result             = result_create_data,
#'   labeled_sequences       = "E:/TATIANA/CAMELLIA-main/inst/extdata",
#'   external_test_fasta_dir = "E:/path/to/test/fasta",
#'   test_predictions        = result_models$predictions_test_ext_xgb
#' )
#' }
#'
#' @importFrom parallel detectCores makeCluster stopCluster clusterEvalQ parLapply clusterExport
#' @importFrom stringi stri_locate_all_fixed
#' @importFrom seqinr read.fasta
#'
#'
#' @export
kmers_in_seq <- function(motifs,
                         data_result,
                         labeled_sequences,
                         external_test_fasta_dir = NULL,
                         test_predictions        = NULL,
                         verbose                 = TRUE) {

  if (!inherits(data_result, "kmer_data")) {
    stop("'data_result' must be an object of class 'kmer_data' from create_data().")
  }

  metadata <- data_result$metadata
  if (is.null(metadata) || nrow(metadata) == 0L) {
    stop("'data_result$metadata' is NULL or empty.")
  }


  n_cores <- max(1L, parallel::detectCores() - 1L)


  # Process motifs
  motifs        <- unlist(motifs, use.names = FALSE)
  unique_motifs <- unique(motifs)

  # Get training class and length lookup
  train_class_lookup  <- setNames(as.character(metadata$class), metadata$sequence_name)
  train_length_lookup <- setNames(as.integer(metadata$length),  metadata$sequence_name)

  has_test <- !is.null(external_test_fasta_dir)

  # Read training sequences
  if (verbose) message("Reading training sequences from: ", labeled_sequences)
  fasta_data     <- .read_fasta_sequences(labeled_sequences)
  sequences      <- fasta_data$sequences
  sequence_names <- fasta_data$names
  n_train        <- length(sequences)

  if (verbose) message("Loaded ", n_train, " training sequences.")

  # Read validation sequences if provided
  tst_sequences      <- NULL
  tst_sequence_names <- NULL
  tst_class_lookup   <- NULL
  tst_length_lookup  <- NULL
  n_test             <- 0L

  if (has_test) {
    if (verbose) message("Reading test sequences from: ", external_test_fasta_dir)
    tst_fasta_data     <- .read_fasta_sequences(external_test_fasta_dir)
    tst_sequences      <- tst_fasta_data$sequences
    tst_sequence_names <- tst_fasta_data$names
    n_test             <- length(tst_sequences)

    if (verbose) message("Loaded ", n_test, " test sequences.")

    # Class lookup from test predictions
    if (!is.null(test_predictions) &&
        "class" %in% colnames(test_predictions)) {
      tst_class_lookup <- setNames(
        as.character(test_predictions$class),
        tst_sequence_names
      )
    } else {
      if (verbose) {
        message("'test_predictions' not provided: class will be NA for test sequences.")
      }
      tst_class_lookup <- setNames(
        rep(NA_character_, length(tst_sequence_names)),
        tst_sequence_names
      )
    }

    # Length lookup from test sequences
    tst_length_lookup <- setNames(
      sapply(tst_sequences, nchar),
      tst_sequence_names
    )
  }

  start_time <- Sys.time()

  # Search in training
  if (verbose) message("\nSearching motifs in TRAINING sequences...")
  train_result <- .find_motifs_parallel(
    sequences, unique_motifs, sequence_names,
    train_class_lookup, train_length_lookup,
    n_cores,
    dataset = if (has_test) "training" else NULL
  )

  # Search motifs in test sequences
  tst_result <- NULL
  if (has_test) {
    if (verbose) message("Searching motifs in TEST sequences...")
    tst_result <- .find_motifs_parallel(
      tst_sequences, unique_motifs, tst_sequence_names,
      tst_class_lookup, tst_length_lookup,
      n_cores,
      dataset = "test"
    )
  }

  # Combine results
  result_df <- if (!is.null(tst_result)) rbind(train_result, tst_result) else train_result
  elapsed   <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

  # Add attributes
  attr(result_df, "n_train_sequences") <- n_train
  attr(result_df, "n_test_sequences")  <- n_test
  attr(result_df, "n_motifs")          <- length(unique_motifs)
  attr(result_df, "n_occurrences")     <- nrow(result_df)
  attr(result_df, "elapsed_time")      <- elapsed
  attr(result_df, "has_test")          <- has_test
  class(result_df) <- c("kmers_in_seq_result", "data.frame")

  # Print summary
  if (verbose) {
    message("Complete.")
    message(sprintf("Total occurrences: %d | Time: %.2f s", nrow(result_df), elapsed))

    if (nrow(result_df) > 0L && has_test && "dataset" %in% colnames(result_df)) {
      message("Occurrences by dataset:")
      message(paste(capture.output(print(table(result_df$dataset))), collapse = "\n"))
    }
  }

  return(result_df)
}
