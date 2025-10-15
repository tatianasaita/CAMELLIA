#' Create k-mer count data from FASTA files
#'
#' @param input Character string. Path to directory containing FASTA files
#' @param output Character string. Path for output file (currently not used)
#' @param k Integer. Length of k-mers to count (default: 3)
#' @param alphabet Character vector. Alphabet to use for k-mer generation (default: DNA_ALPHABET[1:4])
#'
#' @return A list with two elements: kmers (matrix of k-mer counts) and metadata (data frame)
#'
#' @importFrom Biostrings readDNAStringSet DNA_ALPHABET mkAllStrings
#' @importFrom utils write.csv
#' @importFrom tools file_path_sans_ext
#'
#' @keywords internal
#' @noRd
create_data <- function(input, output, k = 3, alphabet = Biostrings::DNA_ALPHABET[1:4]) {

  # Read all fasta files in the input directory
  fasta_files <- list.files(input, pattern = "\\.fasta$", full.names = TRUE)

  if (length(fasta_files) == 0) {
    stop("No .fasta files found in the input directory")
  }

  # Initialize lists to store results
  all_sequences <- list()
  all_classes <- character()
  all_names <- character()
  all_lengths <- integer()
  all_kmer_counts <- list()

  # Generate all possible kmers once
  all_kmers <- Biostrings::mkAllStrings(alphabet, k)

  cat("Processing", length(fasta_files), "fasta files...\n")

  # Process each fasta file
  for (fasta_file in fasta_files) {
    # Extract class name from filename (remove path and extension)
    class_name <- tools::file_path_sans_ext(basename(fasta_file))

    cat("Processing class:", class_name, "\n")

    # Read sequences from fasta file
    sequences <- Biostrings::readDNAStringSet(fasta_file)

    # Process each sequence in the file
    for (i in seq_along(sequences)) {
      seq_name <- names(sequences)[i]
      seq_obj <- sequences[[i]]
      seq_length <- length(seq_obj)

      # Count kmers for this sequence (internal function)
      kmer_counts <- count_kmers(seq_obj, k, alphabet)

      # Store results
      all_names <- c(all_names, seq_name)
      all_lengths <- c(all_lengths, seq_length)
      all_classes <- c(all_classes, class_name)
      all_kmer_counts[[length(all_kmer_counts) + 1]] <- kmer_counts
    }
  }

  cat("Total sequences processed:", length(all_names), "\n")

  # Create metadata dataframe
  metadata <- data.frame(
    sequence_name = all_names,
    length = all_lengths,
    class = all_classes,
    stringsAsFactors = FALSE
  )

  # Create kmer count matrix
  kmer_matrix <- do.call(rbind, all_kmer_counts)
  rownames(kmer_matrix) <- all_names

  # Add CLASS column
  kmer_matrix <- cbind(kmer_matrix, CLASS = all_classes)

  cat("Kmer matrix dimensions:", nrow(kmer_matrix), "x", ncol(kmer_matrix), "\n")
  cat("Done!\n")

  return(list(kmers = kmer_matrix,
              metadata = metadata))
}
