#' Count k-mers in a DNA sequence
#'
#' Internal helper function to count occurrences of all possible k-mers
#' in a DNA sequence. Ensures all possible k-mers are represented in the
#' output, even if their count is zero.
#'
#' @param sequence A DNAString object or character string representing a DNA sequence
#' @param k Integer. Length of k-mers to count
#' @param alphabet Character vector. Alphabet to use for k-mer generation.
#'   Default is DNA_ALPHABET[1:4] (A, C, G, T)
#'
#' @return A named integer vector with counts for all possible k-mers.
#'   Names are the k-mer sequences, values are their counts.
#'
#' @importFrom Biostrings DNAString mkAllStrings oligonucleotideFrequency DNA_ALPHABET
#'
#' @keywords internal
#' @noRd
count_kmers <- function(sequence, k, alphabet = Biostrings::DNA_ALPHABET[1:4]) {
  # Convert to DNAString for efficient processing
  if (!inherits(sequence, "DNAString")) {
    sequence <- Biostrings::DNAString(sequence)
  }

  # Generate all possible kmers
  all_kmers <- Biostrings::mkAllStrings(alphabet, k)

  # Count kmers efficiently using Biostrings
  kmer_counts <- Biostrings::oligonucleotideFrequency(sequence, width = k, as.prob = FALSE)

  # Ensure all kmers are present
  complete_counts <- setNames(rep(0L, length(all_kmers)), all_kmers)
  complete_counts[names(kmer_counts)] <- kmer_counts

  return(complete_counts)
}
