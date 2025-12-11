#' Count K-mers in a DNA Sequence
#'
#' @param sequence A character string, Biostrings XString, or DNAStringSet
#' @param k Integer. Length of k-mers to count
#' @param alphabet Character vector. Default: c("A", "C", "G", "T")
#'
#' @importFrom Biostrings mkAllStrings
#' @importFrom utils head tail
#'
#' @export
count_kmers <- function(sequence,
                        k,
                        alphabet = c("A", "C", "G", "T")) {

  # === PARAMETER VALIDATION ===

  # Validate k parameter
  if (!is.numeric(k) || length(k) != 1) {
    stop("'k' must be a single numeric value", call. = FALSE)
  }

  if (k <= 0 || k != as.integer(k)) {
    stop("'k' must be a positive integer", call. = FALSE)
  }

  # Validate alphabet parameter
  if (!is.character(alphabet)) {
    stop("'alphabet' must be a character vector", call. = FALSE)
  }

  if (length(alphabet) == 0) {
    stop("'alphabet' cannot be empty", call. = FALSE)
  }

  if (length(alphabet) != length(unique(alphabet))) {
    stop("'alphabet' contains duplicate nucleotides", call. = FALSE)
  }

  # === SEQUENCE CONVERSION ===

  # Handle Biostrings objects
  if (methods::is(sequence, "XString")) {
    # Convert Biostrings XString to character
    sequence <- as.character(sequence)
  } else if (methods::is(sequence, "DNAStringSet")) {
    # If DNAStringSet with multiple sequences, take first one
    if (length(sequence) > 1) {
      warning("Multiple sequences in DNAStringSet. Using first sequence.",
              call. = FALSE)
    }
    sequence <- as.character(sequence[[1]])
  } else if (!is.character(sequence)) {
    stop("'sequence' must be a character string or Biostrings object",
         call. = FALSE)
  }

  # === SEQUENCE VALIDATION ===

  if (length(sequence) != 1) {
    stop("'sequence' must be a single character string", call. = FALSE)
  }

  if (nchar(sequence) == 0) {
    stop("'sequence' cannot be empty", call. = FALSE)
  }

  # === SEQUENCE PROCESSING ===

  # Convert to uppercase for consistency
  sequence <- toupper(sequence)

  # Get sequence length
  seq_length <- nchar(sequence)

  # Validate k against sequence length
  if (k > seq_length) {
    stop(paste("'k' (", k, ") cannot be greater than sequence length (",
               seq_length, ")", sep = ""),
         call. = FALSE)
  }

  # === NUCLEOTIDE VALIDATION ===

  # Check if sequence contains only valid nucleotides
  sequence_chars <- unique(strsplit(sequence, "")[[1]])
  invalid_chars <- setdiff(sequence_chars, alphabet)

  if (length(invalid_chars) > 0) {
    stop(paste("Sequence contains invalid nucleotides not in alphabet:",
               paste(invalid_chars, collapse = ", ")),
         call. = FALSE)
  }

  # === K-MER GENERATION ===

  # Generate all possible k-mers from alphabet
  all_kmers <- Biostrings::mkAllStrings(alphabet, k)

  # === K-MER EXTRACTION ===

  # Extract all k-mers from sequence
  n_kmers_in_seq <- seq_length - k + 1
  extracted_kmers <- character(n_kmers_in_seq)

  for (i in seq_len(n_kmers_in_seq)) {
    extracted_kmers[i] <- substr(sequence, i, i + k - 1)
  }

  # === K-MER COUNTING ===

  # Count occurrences of each extracted k-mer using table
  kmer_table <- table(extracted_kmers)

  # === ENSURE COMPLETE REPRESENTATION ===

  # Initialize vector with all possible k-mers set to 0
  complete_counts <- rep(0L, length(all_kmers))
  names(complete_counts) <- all_kmers

  # Update with observed counts from the sequence
  matching_kmers <- intersect(names(kmer_table), all_kmers)

  if (length(matching_kmers) > 0) {
    complete_counts[matching_kmers] <- as.integer(kmer_table[matching_kmers])
  }

  # === CREATE S3 OBJECT ===

  # Set S3 class
  class(complete_counts) <- c("kmer_counts", "integer")

  # Add attributes for metadata
  attr(complete_counts, "k") <- as.integer(k)
  attr(complete_counts, "alphabet") <- alphabet
  attr(complete_counts, "sequence_length") <- as.integer(seq_length)
  attr(complete_counts, "total_kmers") <- as.integer(sum(complete_counts))

  return(complete_counts)
}

#' Print method for kmer_counts objects
#'
#' @param x A kmer_counts object
#' @param ... Additional arguments (ignored)
#' @export
print.kmer_counts <- function(x, ...) {
  k <- attr(x, "k")
  seq_len <- attr(x, "sequence_length")
  total <- attr(x, "total_kmers")

  cat("K-mer Counts\n")
  cat("============\n\n")
  cat("K-mer length (k):", k, "\n")
  cat("Sequence length:", seq_len, "\n")
  cat("Total k-mers:", total, "\n\n")

  # Get top k-mers
  sorted_counts <- sort(x, decreasing = TRUE)
  top_n <- min(10, length(sorted_counts))
  top_kmers <- head(sorted_counts, top_n)

  cat("Top 10 k-mers by frequency:\n")
  for (i in seq_along(top_kmers)) {
    cat(sprintf("  %s: %d\n", names(top_kmers)[i], top_kmers[i]))
  }

  invisible(x)
}

#' Summary method for kmer_counts objects
#'
#' @param object A kmer_counts object
#' @param ... Additional arguments (ignored)
#' @export
summary.kmer_counts <- function(object, ...) {
  k <- attr(object, "k")
  seq_len <- attr(object, "sequence_length")
  total <- attr(object, "total_kmers")
  alphabet <- attr(object, "alphabet")

  cat("K-mer Counts Summary\n")  # Adicionado 's' em "Counts"
  cat("====================\n\n")  # Ajustado o número de '=' para combinar
  cat("K-mer length (k):", k, "\n")
  cat("Alphabet:", paste(alphabet, collapse = ", "), "\n")
  cat("Sequence length:", seq_len, "\n")
  cat("Total k-mers counted:", total, "\n")
  cat("Possible k-mers:", length(object), "\n")
  cat("Observed k-mers:", sum(object > 0), "\n")
  cat("Zero-count k-mers:", sum(object == 0), "\n\n")

  cat("Frequency distribution:\n")
  cat("  Min:", min(object), "\n")
  cat("  Max:", max(object), "\n")
  cat("  Mean:", round(mean(object), 2), "\n")
  cat("  Median:", median(object), "\n")

  invisible(object)
}

#' Subset method for kmer_counts objects
#'
#' @param x A kmer_counts object
#' @param i Index or name for subsetting
#' @param ... Additional arguments (ignored)
#' @export
`[.kmer_counts` <- function(x, i, ...) {
  result <- NextMethod("[")

  # Preserve class and attributes
  class(result) <- c("kmer_counts", "integer")
  attr(result, "k") <- attr(x, "k")
  attr(result, "alphabet") <- attr(x, "alphabet")
  attr(result, "sequence_length") <- attr(x, "sequence_length")
  attr(result, "total_kmers") <- sum(result)

  return(result)
}

#' Subset method for kmer_counts with pattern matching
#'
#' @param x A kmer_counts object
#' @param pattern Regular expression pattern to match k-mer names
#' @param ... Additional arguments (ignored)
#' @export
subset.kmer_counts <- function(x, pattern = NULL, ...) {
  if (is.null(pattern)) {
    stop("'pattern' argument is required for subset.kmer_counts", call. = FALSE)
  }

  # Find k-mers matching the pattern
  matching_indices <- grep(pattern, names(x))

  if (length(matching_indices) == 0) {
    warning("No k-mers match the pattern", call. = FALSE)
  }

  # Use [ method which preserves class
  result <- x[matching_indices]

  return(result)
}

#' Sort method for kmer_counts objects
#'
#' @param x A kmer_counts object
#' @param decreasing Logical. Sort in decreasing order?
#' @param ... Additional arguments (ignored)
#' @export
sort.kmer_counts <- function(x, decreasing = FALSE, ...) {
  result <- x[order(x, decreasing = decreasing)]
  return(result)
}

#' Head method for kmer_counts objects
#'
#' @param x A kmer_counts object
#' @param n Number of elements to return
#' @param ... Additional arguments (ignored)
#' @export
head.kmer_counts <- function(x, n = 6L, ...) {
  n <- min(n, length(x))
  result <- x[seq_len(n)]
  return(result)
}

#' Tail method for kmer_counts objects
#'
#' @param x A kmer_counts object
#' @param n Number of elements to return
#' @param ... Additional arguments (ignored)
#' @export
tail.kmer_counts <- function(x, n = 6L, ...) {
  len <- length(x)
  n <- min(n, len)
  result <- x[seq(len - n + 1, len)]
  return(result)
}

#' Generate all possible k-mers from an alphabet
#'
#' @param alphabet Character vector of nucleotides
#' @param k Integer. Length of k-mers
#' @return Character vector of all possible k-mers
#' @importFrom Biostrings mkAllStrings
#' @export
generate_all_kmers <- function(alphabet, k) {
  if (!is.character(alphabet) || length(alphabet) == 0) {
    stop("'alphabet' must be a non-empty character vector", call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 1 || k <= 0 || k != as.integer(k)) {
    stop("'k' must be a positive integer", call. = FALSE)
  }

  return(as.character(Biostrings::mkAllStrings(alphabet, k)))
}
