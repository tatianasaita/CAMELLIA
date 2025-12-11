#' Create K-mer Count Data from FASTA Files
#'
#' Reads DNA sequences from FASTA files in a directory, counts k-mer occurrences,
#' and returns a list containing k-mer counts and sequence metadata.
#'
#' @param input Character string. Path to directory containing FASTA files
#'   with .fasta extension. The directory must exist and contain at least
#'   one FASTA file.
#'
#' @param k Integer. Length of k-mers to count. Must be a positive integer,
#'   typically between 4 and 8 for DNA sequences. Default: 6.
#'
#' @param alphabet Character vector. Alphabet to use for k-mer generation.
#'   Default: First 4 letters of \code{Biostrings::DNA_ALPHABET}
#'   (A, C, G, T - standard DNA nucleotides).
#'
#' @param verbose Logical. If \code{TRUE} (default), print progress information
#'   to the console. If \code{FALSE}, run silently.
#'
#' @return An object of class \code{"kmer_data"}, which is a list with two elements:
#'   \describe{
#'     \item{kmers}{Data frame with k-mer counts. Rows represent sequences,
#'       columns represent k-mers plus a 'CLASS' column with class labels.
#'       Row names correspond to sequence identifiers from FASTA files.}
#'     \item{metadata}{Data frame with sequence metadata containing:
#'       \describe{
#'         \item{sequence_name}{Sequence identifier from FASTA header}
#'         \item{length}{Sequence length in nucleotides}
#'         \item{class}{Class name extracted from FASTA filename (without extension)}
#'       }}
#'   }
#'
#' @details
#' The function processes all FASTA files in the specified directory in the following steps:
#' \enumerate{
#'   \item Validates input parameters
#'   \item Discovers all .fasta files in the input directory
#'   \item For each FASTA file, reads sequences and counts k-mers
#'   \item Combines results into a single data frame
#'   \item Returns both k-mer counts and metadata
#' }
#'
#' Expected directory structure:
#' \preformatted{
#' input_directory/
#'   ├── class1.fasta
#'   ├── class2.fasta
#'   └── class3.fasta
#' }
#'
#' @examples
#' \dontrun{
#' # Create k-mer data from FASTA files
#' result <- create_data(
#'   input = "/path/to/fasta/directory",
#'   k = 6
#' )
#'
#' # Access k-mer counts
#' kmer_df <- result$kmers
#' head(kmer_df[, 1:10])
#'
#' # Access metadata
#' metadata <- result$metadata
#' head(metadata)
#'
#' # Check class distribution
#' table(metadata$class)
#' }
#'
#' @seealso
#'   \code{\link{count_kmers}} for k-mer counting,
#'   \code{\link[Biostrings]{readDNAStringSet}} for reading FASTA files,
#'   \code{\link[Biostrings]{DNA_ALPHABET}} for available DNA alphabets.
#'
#' @importFrom Biostrings readDNAStringSet mkAllStrings DNA_ALPHABET
#' @importFrom tools file_path_sans_ext
#'
#' @export
create_data <- function(input,
                        k = 6L,
                        alphabet = Biostrings::DNA_ALPHABET[1:4],
                        verbose = TRUE) {

  # ===== PARAMETER VALIDATION =====
  # Validate all parameters BEFORE any file system operations

  # Validate k: check length first
  if (!is.atomic(k) || length(k) != 1L) {
    stop("'k' must be a single positive integer", call. = FALSE)
  }

  # Convert to integer, suppressing warning from non-numeric strings
  k <- suppressWarnings(as.integer(k))

  # Check if k is NA (happens with non-numeric strings or NULL)
  if (is.na(k)) {
    stop("'k' must be a single positive integer", call. = FALSE)
  }

  # Check if k is positive
  if (k <= 0L) {
    stop("'k' must be a single positive integer", call. = FALSE)
  }

  if (k > 10L) {
    warning("Large k value (", k, ") may consume significant memory",
            call. = FALSE)
  }

  # Validate alphabet
  if (!is.character(alphabet) || length(alphabet) == 0L) {
    stop("'alphabet' must be a non-empty character vector", call. = FALSE)
  }

  if (length(alphabet) != length(unique(alphabet))) {
    stop("'alphabet' contains duplicate values", call. = FALSE)
  }

  # Validate verbose
  if (!is.logical(verbose) || length(verbose) != 1L) {
    stop("'verbose' must be a single logical value", call. = FALSE)
  }

  # ===== INPUT VALIDATION =====

  if (!is.character(input) || length(input) != 1L) {
    stop("'input' must be a single character string", call. = FALSE)
  }

  if (!dir.exists(input)) {
    stop("Input directory does not exist: ", input, call. = FALSE)
  }

  # ===== DISCOVER FASTA FILES =====

  fasta_files <- list.files(
    input,
    pattern = "\\.fasta$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(fasta_files) == 0L) {
    stop("No .fasta files found in: ", input, call. = FALSE)
  }

  if (verbose) {
    cat("Found", length(fasta_files), "FASTA file(s)\n")
  }

  # ===== INITIALIZE STORAGE =====

  all_kmers <- Biostrings::mkAllStrings(alphabet, k)
  n_kmers <- length(all_kmers)

  if (verbose) {
    cat("Total possible k-mers (k=", k, "):", n_kmers, "\n", sep = "")
  }

  # Pre-allocate with exact size
  kmer_list <- vector("list", length = length(fasta_files))
  metadata_list <- vector("list", length = length(fasta_files))

  if (verbose) {
    cat("\nProcessing FASTA files:\n")
  }

  # ===== PROCESS EACH FASTA FILE =====

  for (file_idx in seq_along(fasta_files)) {
    fasta_file <- fasta_files[file_idx]
    class_name <- tools::file_path_sans_ext(basename(fasta_file))

    if (verbose) {
      cat("[", file_idx, "/", length(fasta_files), "] ",
          sprintf("%-20s", class_name), " ... ", sep = "")
    }

    tryCatch(
      {
        # Read FASTA file
        sequences <- Biostrings::readDNAStringSet(fasta_file)

        if (length(sequences) == 0L) {
          if (verbose) cat("0 sequences\n")
          next
        }

        n_seq <- length(sequences)

        # Process all sequences at once
        result <- .process_sequences(sequences, all_kmers, k, alphabet, class_name)

        kmer_list[[file_idx]] <- result$kmers
        metadata_list[[file_idx]] <- result$metadata

        if (verbose) {
          cat(n_seq, "sequences\n")
        }
      },
      error = function(e) {
        if (verbose) {
          cat("ERROR: ", conditionMessage(e), "\n", sep = "")
        }
        warning("Error processing file ", basename(fasta_file), ": ",
                conditionMessage(e), call. = FALSE)
      }
    )
  }

  # ===== VALIDATE AND COMBINE RESULTS =====

  # Remove NULL entries
  kmer_list <- kmer_list[!sapply(kmer_list, is.null)]
  metadata_list <- metadata_list[!sapply(metadata_list, is.null)]

  if (length(kmer_list) == 0L) {
    stop("No sequences were successfully processed from any FASTA file",
         call. = FALSE)
  }

  # Combine results
  kmer_df <- do.call(rbind, kmer_list)
  rownames(kmer_df) <- NULL

  metadata_df <- do.call(rbind, metadata_list)
  rownames(metadata_df) <- NULL

  if (verbose) {
    cat("\nTotal sequences processed:", nrow(kmer_df), "\n")
    cat("K-mer matrix dimensions:", nrow(kmer_df), "rows x",
        ncol(kmer_df) - 1L, "columns\n")
    cat("Class distribution:\n")
    print(table(metadata_df$class))
    cat("\n")
  }

  # ===== RETURN RESULTS =====

  result <- list(
    kmers = kmer_df,
    metadata = metadata_df
  )

  class(result) <- c("kmer_data", "list")

  invisible(result)
}


# ==============================================================================
# INTERNAL HELPER FUNCTION
# ==============================================================================

#' Process Sequences from a Single FASTA File
#'
#' Internal helper function that processes all sequences from a FASTA file
#' and returns k-mer counts and metadata.
#'
#' @param sequences A DNAStringSet object containing sequences
#' @param all_kmers Character vector of all possible k-mers
#' @param k Integer k-mer length
#' @param alphabet Character vector of nucleotides
#' @param class_name Character string with the class label
#'
#' @return A list with 'kmers' (data frame) and 'metadata' (data frame)
#'
#' @keywords internal
.process_sequences <- function(sequences, all_kmers, k, alphabet, class_name) {

  n_seq <- length(sequences)
  n_kmers <- length(all_kmers)

  # Pre-allocate matrix
  kmer_matrix <- matrix(0L, nrow = n_seq, ncol = n_kmers + 1L,
                        dimnames = list(NULL, c(all_kmers, "CLASS")))

  seq_names <- character(n_seq)
  seq_lengths <- integer(n_seq)

  # ===== PROCESS EACH SEQUENCE =====

  for (seq_idx in seq_len(n_seq)) {
    seq_name <- names(sequences)[seq_idx]
    seq_obj <- sequences[[seq_idx]]
    seq_length <- length(seq_obj)

    # Count k-mers
    kmer_counts <- count_kmers(seq_obj, k, alphabet)

    # Store in matrix
    kmer_matrix[seq_idx, seq(along = all_kmers)] <- as.integer(kmer_counts[all_kmers])
    kmer_matrix[seq_idx, "CLASS"] <- class_name

    seq_names[seq_idx] <- seq_name
    seq_lengths[seq_idx] <- seq_length
  }

  # ===== CONVERT TO DATA FRAME =====

  kmer_df <- as.data.frame(kmer_matrix, stringsAsFactors = FALSE)
  kmer_df[, seq(along = all_kmers)] <- lapply(kmer_df[, seq(along = all_kmers)],
                                              as.integer)

  # ===== CREATE METADATA =====

  metadata_df <- data.frame(
    sequence_name = seq_names,
    length = seq_lengths,
    class = rep(class_name, n_seq),
    stringsAsFactors = FALSE
  )

  # ===== RETURN RESULTS =====

  list(kmers = kmer_df, metadata = metadata_df)
}


# ==============================================================================
# S3 PRINT METHOD
# ==============================================================================

#' Print Method for kmer_data Objects
#'
#' Displays a summary of the k-mer data created by \code{\link{create_data}}.
#'
#' @param x A kmer_data object created by \code{\link{create_data}}.
#' @param ... Additional arguments (not currently used).
#'
#' @return The input object \code{x}, invisibly.
#'
#' @export
#' @method print kmer_data
print.kmer_data <- function(x, ...) {

  cat("K-mer Data Object\n")
  cat("=================\n\n")

  # ===== K-MER MATRIX SUMMARY =====

  cat("K-mer Matrix:\n")
  cat("  Dimensions:", nrow(x$kmers), "sequences x",
      ncol(x$kmers) - 1L, "k-mers + 1 CLASS column\n")

  kmer_cols <- which(colnames(x$kmers) != "CLASS")
  min_count <- min(x$kmers[, kmer_cols])
  max_count <- max(x$kmers[, kmer_cols])

  cat("  Range of counts: [", min_count, ", ", max_count, "]\n", sep = "")

  # ===== METADATA SUMMARY =====

  cat("\nMetadata:\n")
  cat("  Sequences:", nrow(x$metadata), "\n")
  cat("  Mean sequence length:", round(mean(x$metadata$length), 2), "bp\n")
  cat("  Min sequence length:", min(x$metadata$length), "bp\n")
  cat("  Max sequence length:", max(x$metadata$length), "bp\n")

  # ===== CLASS DISTRIBUTION =====

  cat("\nClass Distribution:\n")
  class_dist <- table(x$metadata$class)

  for (i in seq_along(class_dist)) {
    cat("  ", names(class_dist)[i], ": ", class_dist[i], " sequences\n", sep = "")
  }

  invisible(x)
}
