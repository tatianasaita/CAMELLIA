# tests/testthat/test-count_kmers.R
#
# Essential tests for count_kmers function following CRAN standards
# These tests cover:
# - Parameter validation
# - Basic functionality
# - Edge cases
# - S3 methods
# - Consistency

# ===== PARAMETER VALIDATION TESTS =====

test_that("count_kmers validates k parameter", {
  expect_error(
    count_kmers(sequence = "ACGTACGT", k = "2"),
    "'k' must be a single numeric value"
  )

  expect_error(
    count_kmers(sequence = "ACGTACGT", k = c(1, 2)),
    "'k' must be a single numeric value"
  )

  expect_error(
    count_kmers(sequence = "ACGTACGT", k = 0),
    "'k' must be a positive integer"
  )

  expect_error(
    count_kmers(sequence = "ACGTACGT", k = -1),
    "'k' must be a positive integer"
  )

  expect_error(
    count_kmers(sequence = "ACGTACGT", k = 2.5),
    "'k' must be a positive integer"
  )
})

test_that("count_kmers validates alphabet parameter", {
  expect_error(
    count_kmers(sequence = "ACGTACGT", k = 2, alphabet = 123),
    "'alphabet' must be a character vector"
  )

  expect_error(
    count_kmers(sequence = "ACGTACGT", k = 2, alphabet = character(0)),
    "'alphabet' cannot be empty"
  )

  expect_error(
    count_kmers(sequence = "ACGTACGT", k = 2, alphabet = c("A", "C", "G", "A")),
    "'alphabet' contains duplicate nucleotides"
  )
})

test_that("count_kmers validates sequence parameter", {
  expect_error(
    count_kmers(sequence = 123, k = 2),
    "'sequence' must be a character string"
  )

  expect_error(
    count_kmers(sequence = c("ACGT", "TGCA"), k = 2),
    "'sequence' must be a single character string"
  )

  expect_error(
    count_kmers(sequence = "", k = 2),
    "'sequence' cannot be empty"
  )

  expect_error(
    count_kmers(sequence = "ACGT", k = 5),
    "'k' .* cannot be greater than sequence length"
  )
})

test_that("count_kmers detects invalid nucleotides", {
  expect_error(
    count_kmers(sequence = "ACGTN", k = 2),
    "Sequence contains invalid nucleotides not in alphabet"
  )

  expect_error(
    count_kmers(sequence = "ACGT", k = 2, alphabet = c("A", "C", "G")),
    "Sequence contains invalid nucleotides not in alphabet"
  )
})

# ===== BASIC FUNCTIONALITY TESTS =====

test_that("count_kmers returns all possible k-mers", {
  result <- count_kmers(sequence = "ACGTACGT", k = 2)

  # For k=2 with DNA alphabet: 4^2 = 16 possible k-mers
  expect_equal(length(result), 16L)

  result_k3 <- count_kmers(sequence = "ACGTACGT", k = 3)
  # For k=3 with DNA alphabet: 4^3 = 64 possible k-mers
  expect_equal(length(result_k3), 64L)
})

test_that("count_kmers correctly counts k-mers", {
  result <- count_kmers(sequence = "AAAA", k = 2)

  # "AAAA" contains 3 overlapping "AA" k-mers
  expect_equal(result[["AA"]], 3L)
  # Other k-mers should have 0 count
  expect_equal(result[["AC"]], 0L)
  expect_equal(result[["AG"]], 0L)
  expect_equal(result[["AT"]], 0L)
})

test_that("count_kmers handles k=1", {
  result <- count_kmers(sequence = "ACGTACGT", k = 1)

  expect_equal(length(result), 4L)
  expect_equal(result[["A"]], 2L)
  expect_equal(result[["C"]], 2L)
  expect_equal(result[["G"]], 2L)
  expect_equal(result[["T"]], 2L)
})

test_that("count_kmers handles k equal to sequence length", {
  result <- count_kmers(sequence = "ACGT", k = 4)

  # Only one k-mer possible
  expect_equal(sum(result), 1L)
  expect_equal(result[["ACGT"]], 1L)
  # All other k-mers should be 0
  expect_equal(sum(result == 0), length(result) - 1L)
})

# ===== ATTRIBUTE TESTS =====

test_that("count_kmers sets correct attributes", {
  result <- count_kmers(sequence = "ACGTACGT", k = 2)

  expect_equal(attr(result, "k"), 2L)
  expect_equal(attr(result, "alphabet"), c("A", "C", "G", "T"))
  expect_equal(attr(result, "sequence_length"), 8L)
  expect_equal(attr(result, "total_kmers"), 7L)
})

test_that("count_kmers total_kmers attribute is correct", {
  result <- count_kmers(sequence = "ACGTACGT", k = 2)

  # total_kmers should equal sequence_length - k + 1
  expected_total <- nchar("ACGTACGT") - 2 + 1
  expect_equal(attr(result, "total_kmers"), expected_total)
  expect_equal(sum(result), expected_total)
})

# ===== CUSTOM ALPHABET TESTS =====

test_that("count_kmers handles custom alphabet", {
  result <- count_kmers(sequence = "ATAT", k = 2, alphabet = c("A", "T"))

  # For k=2 with custom alphabet [A,T]: 2^2 = 4 possible k-mers
  expect_equal(length(result), 4L)
  expect_equal(result[["AT"]], 2L)
  expect_equal(result[["TA"]], 1L)
})

test_that("count_kmers with single nucleotide alphabet", {
  result <- count_kmers(sequence = "AAAA", k = 2, alphabet = c("A"))

  # Only one k-mer possible (AA)
  expect_equal(length(result), 1L)
  expect_equal(result[["AA"]], 3L)
})

# ===== CASE INSENSITIVITY TESTS =====

test_that("count_kmers handles lowercase and mixed case", {
  result_upper <- count_kmers(sequence = "ACGTACGT", k = 2)
  result_lower <- count_kmers(sequence = "acgtacgt", k = 2)
  result_mixed <- count_kmers(sequence = "AcGtAcGt", k = 2)

  expect_identical(result_upper, result_lower)
  expect_identical(result_upper, result_mixed)
})

# ===== EDGE CASES TESTS =====

test_that("count_kmers handles homopolymer sequences", {
  result <- count_kmers(sequence = "AAAAAAAAAA", k = 2)

  # 10 A's contain 9 overlapping AA k-mers
  expect_equal(result[["AA"]], 9L)
  expect_equal(sum(result), 9L)
  # All other k-mers should have 0 count
  expect_equal(sum(result == 0), length(result) - 1L)
})

test_that("count_kmers handles alternating sequences", {
  result <- count_kmers(sequence = "ACACACAC", k = 2)

  expect_equal(sum(result), 7L)
  expect_equal(result[["AC"]], 4L)
  expect_equal(result[["CA"]], 3L)
})

test_that("count_kmers handles minimum sequence length", {
  result <- count_kmers(sequence = "A", k = 1)

  expect_equal(result[["A"]], 1L)
  expect_equal(sum(result), 1L)
  expect_equal(sum(result == 0), 3L)  # C, G, T have 0 count
})

# ===== S3 METHOD TESTS =====

test_that("print.kmer_counts works correctly", {
  result <- count_kmers(sequence = "ACGTACGT", k = 2)

  expect_output(
    print(result),
    "K-mer Counts"
  )

  expect_output(
    print(result),
    "K-mer length \\(k\\): 2"
  )

  expect_output(
    print(result),
    "Top 10 k-mers by frequency"
  )
})

test_that("summary.kmer_counts works correctly", {
  result <- count_kmers(sequence = "ACGTACGT", k = 2)

  expect_output(
    summary(result),
    "K-mer Counts Summary"
  )

  expect_output(
    summary(result),
    "Mean:"
  )
})

test_that("[.kmer_counts extraction preserves class", {
  result <- count_kmers(sequence = "ACGTACGT", k = 2)

  # Extract by index
  subset_index <- result[1:5]
  expect_s3_class(subset_index, "kmer_counts")
  expect_equal(length(subset_index), 5L)

  # Extract by name
  subset_name <- result[c("AC", "GT")]
  expect_s3_class(subset_name, "kmer_counts")
  expect_equal(length(subset_name), 2L)
})

test_that("subset.kmer_counts works with patterns", {
  result <- count_kmers(sequence = "ACGTACGT", k = 2)

  # Subset k-mers starting with A
  subset_a <- subset(result, pattern = "^A")
  expect_s3_class(subset_a, "kmer_counts")
  expect_true(all(startsWith(names(subset_a), "A")))

  # Subset k-mers ending with T
  subset_t <- subset(result, pattern = "T$")
  expect_s3_class(subset_t, "kmer_counts")
  expect_true(all(endsWith(names(subset_t), "T")))
})

test_that("sort.kmer_counts works correctly", {
  result <- count_kmers(sequence = "AAACCCGGG", k = 1)

  # Sort decreasing (default)
  sorted_dec <- sort(result, decreasing = TRUE)
  expect_s3_class(sorted_dec, "kmer_counts")
  expect_true(sorted_dec[[1]] >= sorted_dec[[2]])

  # Sort increasing
  sorted_inc <- sort(result, decreasing = FALSE)
  expect_s3_class(sorted_inc, "kmer_counts")
  expect_true(sorted_inc[[1]] <= sorted_inc[[2]])
})

test_that("head.kmer_counts returns top k-mers", {
  result <- count_kmers(sequence = "AAACCCGGG", k = 1)

  top3 <- head(result, n = 3)
  expect_s3_class(top3, "kmer_counts")
  expect_equal(length(top3), 3L)
  # Head should return the highest counts
  expect_true(top3[[1]] >= top3[[3]])
})

test_that("tail.kmer_counts returns bottom k-mers", {
  result <- count_kmers(sequence = "AAACCCGGG", k = 1)

  bottom3 <- tail(result, n = 3)
  expect_s3_class(bottom3, "kmer_counts")
  expect_equal(length(bottom3), 3L)
})

# ===== MATHEMATICAL PROPERTIES TESTS =====

test_that("count_kmers sum equals correct total", {
  seq <- "ACGTACGTACGTACGT"

  for (k_val in 1:4) {
    result <- count_kmers(sequence = seq, k = k_val)
    expected_total <- nchar(seq) - k_val + 1
    expect_equal(sum(result), expected_total)
    expect_equal(attr(result, "total_kmers"), expected_total)
  }
})

test_that("count_kmers k-mers are sorted alphabetically", {
  result <- count_kmers(sequence = "ACGTACGT", k = 2)

  # Names should be in alphabetical order
  expect_equal(names(result), sort(names(result)))
})

test_that("count_kmers includes zero-count k-mers", {
  result <- count_kmers(sequence = "ACACAC", k = 1)

  # G and T should have 0 counts
  expect_equal(result[["G"]], 0L)
  expect_equal(result[["T"]], 0L)
  # A and C should have positive counts
  expect_true(result[["A"]] > 0)
  expect_true(result[["C"]] > 0)
})

# ===== CONSISTENCY TESTS =====

test_that("count_kmers is consistent across runs", {
  result1 <- count_kmers(sequence = "ACGTACGTACGT", k = 2)
  result2 <- count_kmers(sequence = "ACGTACGTACGT", k = 2)

  expect_identical(result1, result2)
})

test_that("count_kmers consistency with different k values", {
  seq <- "ACGTACGTACGTACGT"

  for (k_val in 1:3) {
    result1 <- count_kmers(sequence = seq, k = k_val)
    result2 <- count_kmers(sequence = seq, k = k_val)

    expect_equal(length(result1), 4L^k_val)
    expect_identical(result1, result2)
  }
})

# ===== HELPER FUNCTION TESTS =====

test_that("generate_all_kmers produces correct output", {
  # Test k=1
  kmers_1 <- generate_all_kmers(c("A", "C"), 1)
  expect_equal(sort(kmers_1), c("A", "C"))
  expect_equal(length(kmers_1), 2L)

  # Test k=2
  kmers_2 <- generate_all_kmers(c("A", "C"), 2)
  expect_equal(length(kmers_2), 4L)
  expect_true(all(c("AA", "AC", "CA", "CC") %in% kmers_2))
})

test_that("generate_all_kmers output is sorted", {
  kmers <- generate_all_kmers(c("A", "C", "G", "T"), 2)
  expect_equal(kmers, sort(kmers))
})

