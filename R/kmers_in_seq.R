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
#'
#' @details
#' \itemize{
#'  \item Searches for exact matches of k-mer motifs in DNA sequences
#'  \item Uses parallel processing (automatically detects available CPU cores)
#'  \item Supports both training and validation datasets
#'  \item Training sequences require cluster_result object for metadata (class and length)
#'  \item Validation dataset is optional and requires external FASTA files
#'  \item Returns positions (start/end) of all motif occurrences
#'}
#'
#' @note
#' \itemize{
#'   \item Requires .find_motifs_parallel and .read_fasta_sequences. See \code{internal-functions.R} for details.
#'   \item S3 methods available. See \code{methods.R} for details.
#' }
#'
#' @examples
#'\dontrun{
#' # With validation dataset
#' result_kmers_in_seq <- kmers_in_seq(
#'   motifs = motifs,
#'   cluster_result = result_cluster_dendrogram,
#'   input_dir = "path/to/training",
#'   external_validation_fasta_dir = "path/to/validation",
#'   validation_predictions = result_models$predictions_validation_xgb #or predictions_validation_rf
#' )
#'
#' @importFrom parallel detectCores
#' @importFrom stringi stri_locate_all_fixed
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel clusterEvalQ
#' @importFrom parallel parLapply
#' @importFrom parallel clusterExport
#' @importFrom seqinr read.fasta
#'
#'
#' @export
kmers_in_seq <- function(motifs,
                         cluster_result,
                         sequences = NULL,
                         input_dir = NULL,
                         sequence_names = NULL,
                         external_validation_fasta_dir = NULL,
                         validation_predictions = NULL,
                         verbose = TRUE) {

  metadata <- cluster_result$data_result$metadata

  n_cores <- max(1, parallel::detectCores() - 1)

  # Process motifs
  motifs <- unlist(motifs, use.names = FALSE)
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
    fasta_data <- .read_fasta_sequences(input_dir)
    sequences <- fasta_data$sequences
    sequence_names <- fasta_data$names
  }

  if (is.null(sequence_names)) {
    sequence_names <- if (!is.null(names(sequences))) names(sequences) else sprintf("seq_%d", seq_along(sequences))
  }

  n_train <- length(sequences)
  n_val <- 0

  # Read validation sequences if provided
  val_sequences <- NULL
  val_sequence_names <- NULL
  val_class_lookup <- NULL
  val_length_lookup <- NULL

  if (has_validation) {
    val_fasta_data <- .read_fasta_sequences(external_validation_fasta_dir)
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
  }

  start_time <- Sys.time()

  # Search in training
  if (verbose) message("Searching motifs in TRAINING sequences...")
  train_result <- .find_motifs_parallel(
    sequences, unique_motifs, sequence_names, train_class_lookup, train_length_lookup,
    n_cores,  dataset = if (has_validation) "training" else NULL
  )

  # Search in validation
  val_result <- NULL
  if (has_validation) {
    if (verbose) message("\nSearching motifs in VALIDATION sequences...")
    val_result <- .find_motifs_parallel(
      val_sequences, unique_motifs, val_sequence_names, val_class_lookup, val_length_lookup,
      n_cores, dataset = "validation"
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
    message(sprintf("\n Complete."))
    message(sprintf("Total occurrences: %d | Time: %.2f s", nrow(result_df), elapsed))

    if (nrow(result_df) > 0 && has_validation && "dataset" %in% colnames(result_df)) {
      message("\nOccurrences by dataset:")
      print(table(result_df$dataset))
    }
  }

  return(result_df)
}
