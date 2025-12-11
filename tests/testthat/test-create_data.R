# tests/testthat/test-create_data.R
#
# Essential tests for create_data function following CRAN standards
# These tests cover:
# - Parameter validation
# - Basic functionality
# - Error handling
# - Return structure
# - Edge cases
library(testthat)
# ===== SETUP: CREATE TEMPORARY TEST DATA =====

# Helper function to create temporary FASTA files
create_test_fasta <- function(temp_dir, class_name, sequences) {
  fasta_path <- file.path(temp_dir, paste0(class_name, ".fasta"))

  # Create simple FASTA content
  fasta_content <- character()
  for (i in seq_along(sequences)) {
    fasta_content <- c(
      fasta_content,
      paste0(">", class_name, "_seq", i),
      sequences[i]
    )
  }

  writeLines(fasta_content, fasta_path)
  invisible(fasta_path)
}

# ===== PARAMETER VALIDATION TESTS =====

test_that("create_data validates input parameter", {
  # input must be character
  expect_error(
    create_data(input = 123),
    "'input' must be a single character string"
  )

  # input must be single string
  expect_error(
    create_data(input = c("/path1", "/path2")),
    "'input' must be a single character string"
  )

  # input directory must exist
  expect_error(
    create_data(input = "/nonexistent/path/to/directory"),
    "Input directory does not exist"
  )
})

test_that("create_data validates k parameter", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, paste0("test_k_", format(Sys.time(), "%s")))

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir, showWarnings = FALSE)

  writeLines(c(">seq1", "ACGTACGT"), file.path(test_dir, "dummy.fasta"))

  # Test k = 0 (zero is invalid)
  expect_error(
    create_data(input = test_dir, k = 0, verbose = FALSE),
    "'k' must be a single positive integer"
  )

  # Test k = -1 (negative is invalid)
  expect_error(
    create_data(input = test_dir, k = -1, verbose = FALSE),
    "'k' must be a single positive integer"
  )

  # Test k = c(5, 6) (multiple values are invalid)
  expect_error(
    create_data(input = test_dir, k = c(5, 6), verbose = FALSE),
    "'k' must be a single positive integer"
  )

  unlink(test_dir, recursive = TRUE)
})

test_that("create_data validates alphabet parameter", {
  temp_dir <- tempdir()

  expect_error(
    create_data(input = temp_dir, alphabet = 123),
    "'alphabet' must be a non-empty character vector"
  )

  expect_error(
    create_data(input = temp_dir, alphabet = character(0)),
    "'alphabet' must be a non-empty character vector"
  )

  expect_error(
    create_data(input = temp_dir, alphabet = c("A", "C", "G", "A")),
    "'alphabet' contains duplicate values"
  )
})

test_that("create_data validates verbose parameter", {
  temp_dir <- tempdir()

  expect_error(
    create_data(input = temp_dir, verbose = "TRUE"),
    "'verbose' must be a single logical value"
  )

  expect_error(
    create_data(input = temp_dir, verbose = c(TRUE, FALSE)),
    "'verbose' must be a single logical value"
  )
})

test_that("create_data warns about large k values", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, paste0("test_warn_", format(Sys.time(), "%s")))

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir, showWarnings = FALSE)

  # Create FASTA file with reasonable length sequence
  writeLines(
    c(">seq1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"),  # 40 bp
    file.path(test_dir, "dummy.fasta")
  )

  # k > 10 should produce warning
  expect_warning(
    create_data(input = test_dir, k = 11, verbose = FALSE),
    "Large k value"
  )

  unlink(test_dir, recursive = TRUE)
})

# ===== BASIC FUNCTIONALITY TESTS =====

test_that("create_data fails with no FASTA files", {
  temp_dir <- tempdir()
  empty_dir <- file.path(temp_dir, "empty_dir_test")

  if (dir.exists(empty_dir)) {
    unlink(empty_dir, recursive = TRUE)
  }
  dir.create(empty_dir)

  expect_error(
    create_data(input = empty_dir, verbose = FALSE),
    "No .fasta files found"
  )

  unlink(empty_dir, recursive = TRUE)
})

test_that("create_data processes single FASTA file", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_single_fasta")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  # Create test FASTA file
  create_test_fasta(test_dir, "class1", c(
    "ACGTACGTACGTACGT",
    "TGCATGCATGCATGCA"
  ))

  result <- create_data(input = test_dir, k = 2, verbose = FALSE)

  expect_s3_class(result, "kmer_data")
  expect_equal(length(result), 2L)
  expect_true("kmers" %in% names(result))
  expect_true("metadata" %in% names(result))

  unlink(test_dir, recursive = TRUE)
})

test_that("create_data processes multiple FASTA files", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_multiple_fasta")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  # Create multiple FASTA files
  create_test_fasta(test_dir, "class1", c("ACGTACGT", "TGCATGCA"))
  create_test_fasta(test_dir, "class2", c("GGGGGGGG", "CCCCCCCC"))

  result <- create_data(input = test_dir, k = 2, verbose = FALSE)

  # Should have 4 sequences total
  expect_equal(nrow(result$metadata), 4L)
  expect_equal(nrow(result$kmers), 4L)

  # Should have 2 classes
  expect_equal(length(unique(result$metadata$class)), 2L)

  unlink(test_dir, recursive = TRUE)
})

# ===== RETURN STRUCTURE TESTS =====

test_that("create_data returns correct structure", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_structure")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  create_test_fasta(test_dir, "test_class", c("ACGTACGTACGT"))

  result <- create_data(input = test_dir, k = 2, verbose = FALSE)

  # Check kmers data frame
  expect_true(is.data.frame(result$kmers))
  expect_equal(nrow(result$kmers), 1L)
  expect_true("CLASS" %in% colnames(result$kmers))

  # Check metadata data frame
  expect_true(is.data.frame(result$metadata))
  expect_equal(nrow(result$metadata), 1L)
  expect_equal(colnames(result$metadata),
               c("sequence_name", "length", "class"))

  unlink(test_dir, recursive = TRUE)
})

test_that("create_data kmers have correct columns", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_columns")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  create_test_fasta(test_dir, "test", c("ACGTACGT"))

  result <- create_data(input = test_dir, k = 2, verbose = FALSE)

  # For k=2, should have 16 k-mers + CLASS column
  expected_cols <- 16L + 1L  # 4^2 + CLASS
  expect_equal(ncol(result$kmers), expected_cols)

  # All k-mer columns should be integer
  kmer_cols <- which(colnames(result$kmers) != "CLASS")
  for (col in kmer_cols) {
    expect_true(is.integer(result$kmers[[col]]))
  }

  unlink(test_dir, recursive = TRUE)
})

test_that("create_data metadata is complete", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_metadata")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  create_test_fasta(test_dir, "myclass", c(
    "ACGTACGTACGT",
    "TGCATGCATGCATGCA"
  ))

  result <- create_data(input = test_dir, k = 2, verbose = FALSE)

  # Check metadata content
  expect_equal(nrow(result$metadata), 2L)
  expect_true(all(result$metadata$class == "myclass"))
  expect_equal(result$metadata$length, c(12L, 16L))
  expect_true(all(!is.na(result$metadata$sequence_name)))

  unlink(test_dir, recursive = TRUE)
})

# ===== K-MER COUNTING VALIDATION TESTS =====

test_that("create_data counts k-mers correctly", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_counting")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  # Create sequence with known k-mer content
  create_test_fasta(test_dir, "test", c("AAAA"))

  result <- create_data(input = test_dir, k = 2, verbose = FALSE)

  # "AAAA" should have 3 AA k-mers
  aa_col <- which(colnames(result$kmers) == "AA")
  expect_equal(result$kmers[1, aa_col], 3L)

  unlink(test_dir, recursive = TRUE)
})

test_that("create_data respects custom alphabet", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_alphabet")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  create_test_fasta(test_dir, "test", c("ATAT"))

  # Custom alphabet with 2 nucleotides
  result <- create_data(
    input = test_dir,
    k = 2,
    alphabet = c("A", "T"),
    verbose = FALSE
  )

  # Should have 2^2 = 4 k-mers + CLASS
  expect_equal(ncol(result$kmers), 5L)

  unlink(test_dir, recursive = TRUE)
})

# ===== CLASS EXTRACTION TESTS =====

test_that("create_data extracts class from filename", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_class")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  create_test_fasta(test_dir, "bacteria", c("ACGTACGT"))
  create_test_fasta(test_dir, "virus", c("TGCATGCA"))

  result <- create_data(input = test_dir, k = 2, verbose = FALSE)

  # Check that classes are correctly extracted
  expect_equal(sort(unique(result$metadata$class)),
               c("bacteria", "virus"))
  expect_equal(result$kmers$CLASS[1], "bacteria")
  expect_equal(result$kmers$CLASS[2], "virus")

  unlink(test_dir, recursive = TRUE)
})

# ===== ERROR HANDLING TESTS =====

test_that("create_data handles invalid FASTA gracefully", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_invalid")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  # Create valid and invalid FASTA files
  create_test_fasta(test_dir, "valid", c("ACGTACGT"))

  # Create invalid FASTA file
  writeLines("This is not a valid FASTA file",
             file.path(test_dir, "invalid.fasta"))

  # Should handle error gracefully with warning
  expect_warning(
    result <- create_data(input = test_dir, k = 2, verbose = FALSE),
    "Error processing file"
  )

  # Should still return valid data from the valid file
  expect_true(!is.null(result))

  unlink(test_dir, recursive = TRUE)
})

test_that("create_data stops if no sequences process successfully", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_no_sequences")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  # Create FASTA file with only header, no sequences
  writeLines(">empty_sequence",
             file.path(test_dir, "empty.fasta"))

  # ← FIXED: Suppress the warning to avoid it appearing in test output
  expect_error(
    suppressWarnings(
      create_data(input = test_dir, k = 2, verbose = FALSE)
    ),
    "No sequences were successfully processed"
  )

  unlink(test_dir, recursive = TRUE)
})

# ===== S3 METHOD TESTS =====

test_that("print.kmer_data works correctly", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_print")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  create_test_fasta(test_dir, "class1", c("ACGTACGTACGT"))
  result <- create_data(input = test_dir, k = 2, verbose = FALSE)

  expect_output(
    print(result),
    "K-mer Data Object"
  )

  expect_output(
    print(result),
    "K-mer Matrix"
  )

  expect_output(
    print(result),
    "Class Distribution"
  )

  unlink(test_dir, recursive = TRUE)
})

# ===== CONSISTENCY TESTS =====

test_that("create_data is reproducible", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_repro")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  create_test_fasta(test_dir, "test", c(
    "ACGTACGTACGT",
    "TGCATGCATGCA"
  ))

  result1 <- create_data(input = test_dir, k = 2, verbose = FALSE)
  result2 <- create_data(input = test_dir, k = 2, verbose = FALSE)

  expect_identical(result1$kmers, result2$kmers)
  expect_identical(result1$metadata, result2$metadata)

  unlink(test_dir, recursive = TRUE)
})

test_that("create_data produces consistent results for different k", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_k_values")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  create_test_fasta(test_dir, "test", c("ACGTACGTACGT"))

  for (k_val in c(2, 3, 4)) {
    result <- create_data(input = test_dir, k = k_val, verbose = FALSE)

    # Check that number of k-mer columns is 4^k
    expected_cols <- 4L^k_val + 1L  # +1 for CLASS
    expect_equal(ncol(result$kmers), expected_cols)

    # Check metadata is always the same
    expect_equal(nrow(result$metadata), 1L)
  }

  unlink(test_dir, recursive = TRUE)
})

# ===== VERBOSE OUTPUT TEST =====

test_that("create_data verbose option works", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  test_dir <- file.path(temp_dir, "test_verbose")

  if (dir.exists(test_dir)) {
    unlink(test_dir, recursive = TRUE)
  }
  dir.create(test_dir)

  create_test_fasta(test_dir, "test", c("ACGTACGT"))

  # With verbose = TRUE, should produce output
  expect_output(
    create_data(input = test_dir, k = 2, verbose = TRUE),
    "Found"
  )

  # With verbose = FALSE, no output
  expect_silent(
    create_data(input = test_dir, k = 2, verbose = FALSE)
  )

  unlink(test_dir, recursive = TRUE)
})

