# tests/testthat/test-count_kmers.R
#
# Essential tests for count_kmers function following CRAN standards
# These tests cover:
# - Parameter validation
# - Basic functionality
# - Edge cases
# - Consistency
# - Mathematical properties

library(testthat)

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
})

test_that("count_kmers validates alphabet parameter", {
  # Test 1: alphabet não é character
  expect_error(
    count_kmers(sequence = "ACGTACGT", k = 2, alphabet = 123),
    "'alphabet' must be a character vector"
  )

  # Test 2: alphabet vazio
  expect_error(
    count_kmers(sequence = "ACGTACGT", k = 2, alphabet = character(0)),
    "'alphabet' cannot be empty"
  )

  # Test 3: alphabet com duplicatas
  expect_error(
    count_kmers(sequence = "ACGTACGT", k = 2, alphabet = c("A", "C", "G", "T", "A")),
    "'alphabet' contains duplicate nucleotides"
  )

  # Test 4: sequência contém nucleotídeos não no alphabet
  expect_error(
    count_kmers(sequence = "ACGTN", k = 2, alphabet = c("A", "C", "G", "T")),
    "Sequence contains invalid nucleotides not in alphabet"
  )
})

test_that("count_kmers validates sequence parameter", {
  expect_error(
    count_kmers(sequence = 123, k = 2),
    "sequence.*must be.*character.*Biostrings"
  )

  expect_error(
    count_kmers(sequence = c("ACGT", "TGCA"), k = 2),
    "sequence.*must be.*single character string"
  )

  expect_error(
    count_kmers(sequence = "", k = 2),
    "sequence.*cannot be empty"
  )

  expect_error(
    count_kmers(sequence = "ACGT", k = 5),
    "k.*cannot be greater than sequence length"
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

test_that("count_kmers attribute k is integer", {
  result <- count_kmers(sequence = "ACGTACGT", k = 2)

  expect_true(is.integer(attr(result, "k")))
  expect_true(is.integer(attr(result, "sequence_length")))
  expect_true(is.integer(attr(result, "total_kmers")))
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

test_that("count_kmers respects custom alphabet size", {
  test_cases <- list(
    list(seq = "ACGTACGT", alphabet = c("A", "C", "G", "T"), expected_length = 16),
    list(seq = "ACGACGACG", alphabet = c("A", "C", "G"), expected_length = 9),
    list(seq = "ACACACA", alphabet = c("A", "C"), expected_length = 4),
    list(seq = "ATAT", alphabet = c("A", "T"), expected_length = 4),
    list(seq = "AAAA", alphabet = c("A"), expected_length = 1)
  )

  for (test_case in test_cases) {
    result <- count_kmers(
      sequence = test_case$seq,
      k = 2,
      alphabet = test_case$alphabet
    )
    expect_equal(
      length(result),
      test_case$expected_length,
      info = sprintf("Failed for alphabet: %s", paste(test_case$alphabet, collapse = ", "))
    )
  }
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

test_that("count_kmers handles long sequences", {
  # Create a longer sequence
  long_seq <- paste(rep("ACGT", 100), collapse = "")
  result <- count_kmers(sequence = long_seq, k = 2)

  # Should still have all 16 possible 2-mers
  expect_equal(length(result), 16L)
  # Total k-mers should equal length - k + 1
  expected_total <- nchar(long_seq) - 2 + 1
  expect_equal(sum(result), expected_total)
})

# ===== RETURN TYPE AND STRUCTURE TESTS =====

test_that("count_kmers returns named integer vector", {
  result <- count_kmers(sequence = "ACGTACGT", k = 2)

  # Check that it's a numeric/integer type
  expect_true(is.integer(result) || is.numeric(result))

  # Check that it's not a matrix or array
  expect_false(is.matrix(result))
  expect_false(is.array(result))

  # Check that it's a vector-like object
  expect_true(is.atomic(result))

  # Check that all values are named
  expect_true(!is.null(names(result)))
  expect_equal(length(names(result)), length(result))

  # Check that names are character strings
  expect_true(is.character(names(result)))
})

test_that("count_kmers names are valid k-mers", {
  result <- count_kmers(sequence = "ACGTACGT", k = 2)

  # All names should be character strings of length k
  for (name in names(result)) {
    expect_equal(nchar(name), 2L)
    expect_true(all(strsplit(name, "")[[1]] %in% c("A", "C", "G", "T")))
  }
})

test_that("count_kmers all counts are non-negative integers", {
  result <- count_kmers(sequence = "ACGTACGTACGT", k = 2)

  expect_true(all(result >= 0))
  expect_true(all(is.integer(result)))
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

test_that("count_kmers contains all possible k-mers", {
  result <- count_kmers(sequence = "AAAA", k = 1)

  # All 4 nucleotides should be present
  expect_equal(length(result), 4L)
  expect_true(all(c("A", "C", "G", "T") %in% names(result)))

  result_k2 <- count_kmers(sequence = "ACGTACGT", k = 2)
  # All 16 possible 2-mers should be present
  expect_equal(length(result_k2), 16L)
})

# ===== CONSISTENCY AND REPRODUCIBILITY TESTS =====

test_that("count_kmers is reproducible", {
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

test_that("count_kmers order independence", {
  # Same nucleotides in different order should give different k-mer profiles
  seq1 <- "AAACCCGGG"
  seq2 <- "ACGACGACG"

  result1 <- count_kmers(sequence = seq1, k = 2)
  result2 <- count_kmers(sequence = seq2, k = 2)

  # Results should be different
  expect_false(identical(result1, result2))
})

# ===== BOUNDARY VALUE TESTS =====

test_that("count_kmers with k=1 and k=sequence_length", {
  seq <- "ACGTACGT"

  # k=1
  result_k1 <- count_kmers(sequence = seq, k = 1)
  expect_equal(sum(result_k1), nchar(seq))

  # k=sequence_length
  result_k_max <- count_kmers(sequence = seq, k = nchar(seq))
  expect_equal(sum(result_k_max), 1L)
  expect_equal(result_k_max[[seq]], 1L)
})

test_that("count_kmers handles large k values", {
  seq <- "ACGTACGTACGTACGT"

  # k = 10 (mais realista e rápido que k = 15)
  result <- count_kmers(sequence = seq, k = 10)
  expected_count <- nchar(seq) - 10 + 1
  expect_equal(sum(result), expected_count)
  expect_equal(attr(result, "total_kmers"), expected_count)
})

# ===== REPEATED PATTERN TESTS =====

test_that("count_kmers with repeating patterns", {
  # Simple repeat
  result <- count_kmers(sequence = "ATATAT", k = 2)

  expect_equal(result[["AT"]], 3L)
  expect_equal(result[["TA"]], 2L)
  expect_equal(sum(result), 5L)  # 6 - 2 + 1
})

test_that("count_kmers with tandem repeats", {
  # Tandem repeat: ACGT repeated 3 times
  result <- count_kmers(sequence = "ACGTACGTACGT", k = 4)

  # "ACGT" should appear 3 times
  expect_equal(result[["ACGT"]], 3L)
  expect_equal(sum(result), 9L)  # 12 - 4 + 1
})

# ===== DEFAULT PARAMETER TESTS =====

test_that("count_kmers with explicit k parameter", {
  seq <- paste(rep("ACGTACGTACGTACGT", 2), collapse = "")

  # Test with k=6 explicitly
  result <- count_kmers(sequence = seq, k = 6)

  expect_equal(attr(result, "k"), 6L)
  expect_equal(length(result), 4096L)  # 4^6
})

test_that("count_kmers default alphabet parameter", {
  result <- count_kmers(sequence = "ACGTACGT", k = 2)

  # Default alphabet should be A, C, G, T
  expected_kmers <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT",
                      "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
  expect_equal(sort(names(result)), sort(expected_kmers))
})

# ===== XString OBJECT SUPPORT TEST =====

test_that("count_kmers handles Biostrings XString objects", {
  skip_if_not_installed("Biostrings")

  dna_string <- Biostrings::DNAString("ACGTACGT")
  result_xstring <- count_kmers(sequence = dna_string, k = 2)
  result_char <- count_kmers(sequence = "ACGTACGT", k = 2)

  expect_identical(result_xstring, result_char)
})

