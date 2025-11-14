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
#'   Default: \code{Biostrings::DNA_ALPHABET[1:4]} (A, C, G, T).
#'
#' @return A list with two elements:
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
#' The function processes all FASTA files in the specified directory:
#' \enumerate{
#'   \item Validates input parameters
#'   \item Discovers all .fasta files in the input directory
#'   \item For each FASTA file:
#'     \itemize{
#'       \item Extracts class name from filename
#'       \item Reads DNA sequences using \code{\link[Biostrings]{readDNAStringSet}}
#'       \item Counts k-mers for each sequence
#'       \item Stores metadata and k-mer counts
#'     }
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
#'
#' # Filter sequences by class
#' hiv_indices <- which(metadata$class == "HIV")
#' hiv_kmers <- result$kmers[hiv_indices, ]
#' }
#'
#' @seealso
#'   \code{\link[Biostrings]{readDNAStringSet}} for reading FASTA files,
#'   \code{\link[Biostrings]{DNA_ALPHABET}} for available alphabets.
#'
#' @importFrom Biostrings readDNAStringSet DNA_ALPHABET mkAllStrings
#' @importFrom tools file_path_sans_ext
#'
#' @export
create_data <- function(input,
                        k = 6L,
                        alphabet = Biostrings::DNA_ALPHABET[1:4]) {

  # ===== INPUT VALIDATION =====

  if (!is.character(input) || length(input) != 1L) {
    stop("'input' must be a single character string", call. = FALSE)
  }

  if (!dir.exists(input)) {
    stop("Input directory does not exist: ", input, call. = FALSE)
  }

  k <- as.integer(k)

  if (length(k) != 1L || k <= 0L) {
    stop("'k' must be a single positive integer", call. = FALSE)
  }

  if (k > 10L) {
    warning("Large k value (", k, ") may consume significant memory",
            call. = FALSE)
  }

  if (!is.character(alphabet) || length(alphabet) == 0L) {
    stop("'alphabet' must be a non-empty character vector", call. = FALSE)
  }

  if (length(alphabet) != length(unique(alphabet))) {
    stop("'alphabet' contains duplicate values", call. = FALSE)
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

  cat("Found", length(fasta_files), "FASTA file(s)\n")

  # ===== INITIALIZE STORAGE =====

  all_kmers <- Biostrings::mkAllStrings(alphabet, k)
  n_kmers <- length(all_kmers)

  cat("Total possible k-mers (k=", k, "):", n_kmers, "\n", sep = "")

  # Pre-allocate lists for efficiency
  all_kmer_counts <- vector("list", length = length(fasta_files))
  all_names <- vector("character")
  all_lengths <- vector("integer")
  all_classes <- vector("character")

  seq_count <- 0L

  # ===== PROCESS EACH FASTA FILE =====

  cat("\nProcessing FASTA files:\n")

  for (file_idx in seq_along(fasta_files)) {
    fasta_file <- fasta_files[file_idx]
    class_name <- tools::file_path_sans_ext(basename(fasta_file))

    cat("[", file_idx, "/", length(fasta_files), "] ", class_name,
        " ... ", sep = "")

    tryCatch(
      {
        # Read FASTA file
        sequences <- Biostrings::readDNAStringSet(fasta_file)

        if (length(sequences) == 0L) {
          warning("No sequences found in: ", fasta_file, call. = FALSE)
          cat("0 sequences\n")
          next
        }

        n_seq <- length(sequences)
        file_kmer_matrix <- matrix(0L, nrow = n_seq, ncol = n_kmers,
                                   dimnames = list(NULL, all_kmers))

        file_names <- character(n_seq)
        file_lengths <- integer(n_seq)
        file_classes <- rep(class_name, n_seq)

        # ===== PROCESS EACH SEQUENCE =====

        for (seq_idx in seq_len(n_seq)) {
          seq_name <- names(sequences)[seq_idx]
          seq_obj <- sequences[[seq_idx]]
          seq_length <- length(seq_obj)

          # Count k-mers for this sequence
          kmer_counts <- count_kmers(seq_obj, k, alphabet)

          # Store results
          file_names[seq_idx] <- seq_name
          file_lengths[seq_idx] <- seq_length
          file_kmer_matrix[seq_idx, ] <- as.integer(kmer_counts[all_kmers])
        }

        # Store file results
        file_kmer_df <- as.data.frame(file_kmer_matrix, stringsAsFactors = FALSE)
        file_kmer_df$CLASS <- file_classes
        rownames(file_kmer_df) <- file_names

        all_kmer_counts[[file_idx]] <- file_kmer_df

        all_names <- c(all_names, file_names)
        all_lengths <- c(all_lengths, file_lengths)
        all_classes <- c(all_classes, file_classes)

        seq_count <- seq_count + n_seq

        cat(n_seq, "sequences\n")
      },
      error = function(e) {
        warning("Error processing file ", basename(fasta_file), ": ",
                conditionMessage(e), call. = FALSE)
        cat("ERROR\n")
      }
    )
  }

  # ===== VALIDATE RESULTS =====

  if (seq_count == 0L) {
    stop("No sequences were successfully processed from any FASTA file",
         call. = FALSE)
  }

  cat("\nTotal sequences processed:", seq_count, "\n")

  # ===== COMBINE RESULTS =====

  # Remove NULL entries from failed files
  all_kmer_counts <- all_kmer_counts[!sapply(all_kmer_counts, is.null)]

  # Combine k-mer matrices
  kmer_df <- do.call(rbind, all_kmer_counts)
  rownames(kmer_df) <- all_names

  # Ensure numeric columns (except CLASS)
  kmer_cols <- setdiff(colnames(kmer_df), "CLASS")
  kmer_df[, kmer_cols] <- lapply(kmer_df[, kmer_cols], as.numeric)

  cat("K-mer matrix dimensions:", nrow(kmer_df), "rows x",
      ncol(kmer_df), "columns\n")

  # ===== CREATE METADATA =====

  metadata <- data.frame(
    sequence_name = all_names,
    length = all_lengths,
    class = all_classes,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  cat("Class distribution:\n")
  print(table(metadata$class))
  cat("\n")

  # ===== RETURN RESULTS =====

  return(list(
    kmers = kmer_df,
    metadata = metadata
  ))
}
