#' Create K-mer Count Data from FASTA Files
#'
#' Reads DNA sequences from FASTA files in a directory, counts k-mer occurrences,
#' and returns a list containing k-mer counts and sequence metadata. Used for training model.
#'
#' @param input Character string. Path to directory containing FASTA files
#'   with .fasta extension. The directory must exist and contain at least
#'   one FASTA file.
#' @param k Integer. Length of k-mers to count. Must be a positive integer,
#'   typically between 4 and 8 for DNA sequences. Default: 6.
#' @param alphabet Character vector. Alphabet to use for k-mer generation.
#'   Default: First 4 letters of \code{Biostrings::DNA_ALPHABET}
#'   (A, C, G, T - standard DNA nucleotides).
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
#' @details The function processes all FASTA files in the specified directory in the following steps:
#' \enumerate{
#'   \item Discovers all .fasta files in the input directory
#'   \item For each FASTA file, reads sequences and counts k-mers
#'   \item Combines results into a single data frame
#'   \item Returns both k-mer counts and metadata
#' }
#'
#' @note
#' \itemize{
#'   \item Requires .process_sequences and .count_kmers. See \code{internal-functions.R} for details.
#'   \item S3 methods available. See \code{methods.R} for details.
#' }
#'
#' @examples
#' \dontrun{
#' # Create k-mer data from FASTA files
#' result <- create_data(
#'   input = "/path/to/fasta/directory",
#'   k = 6
#' )
#'}
#'
#' @importFrom Biostrings DNA_ALPHABET
#' @importFrom Biostrings mkAllStrings
#' @importFrom Biostrings readDNAStringSet
#' @importFrom tools file_path_sans_ext
#' @importFrom methods is
#'
#' @export
create_data <- function(input,
                        k = 6L,
                        alphabet = Biostrings::DNA_ALPHABET[1:4], # "A", "C", "G", "T"
                        verbose = TRUE) {

  # Discover and identify FASTA files in the directory
  fasta_files <- list.files(
    input,
    pattern = "\\.fasta$",
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(fasta_files) == 0L) {
    stop("No fasta files found in: ", input, call. = FALSE)
  }

  if (verbose) {
    cat("Found", length(fasta_files), "FASTA file(s)\n")
  }

  # All k-mers
  all_kmers <- Biostrings::mkAllStrings(alphabet, k)  #Generate all possible k-mers for parameters k and alphabet

  n_kmers <- length(all_kmers) #all_kmers count

  kmer_list <- vector("list", length = length(fasta_files)) #create list: kmer_list
  metadata_list <- vector("list", length = length(fasta_files)) #create list:metadata_list

  # Read FASTA file
  for (file_idx in seq_along(fasta_files)) {
    fasta_file <- fasta_files[file_idx]
    class_name <- tools::file_path_sans_ext(basename(fasta_file))

    if (verbose) {
      cat("[", file_idx, "/", length(fasta_files), "] ",
          sprintf("%-20s", class_name), " ... ", sep = "")
    } # Class of fasta files

    sequences <- Biostrings::readDNAStringSet(fasta_file)

    if (length(sequences) == 0L) {
      if (verbose) cat("0 sequences\n")
      next
    }

    n_seq <- length(sequences)

    # Process all sequences at once
    result <- .process_sequences(sequences, all_kmers, k, alphabet, class_name) #internal help function, see internal-functions.R

    kmer_list[[file_idx]] <- result$kmers
    metadata_list[[file_idx]] <- result$metadata

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

    # Return results
    result <- list(
    kmers = kmer_df,
    metadata = metadata_df
    )

    class(result) <- "kmer_data"
    invisible(result)
}

