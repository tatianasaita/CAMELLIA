# tests/testthat/test-create_data.R

test_that("create_data requires valid input directory", {
  expect_error(
    create_data(input = NULL),
    "'input' must be a single character string"
  )

  expect_error(
    create_data(input = c("dir1", "dir2")),
    "'input' must be a single character string"
  )

  expect_error(
    create_data(input = "/nonexistent/directory"),
    "Input directory does not exist"
  )
})

test_that("create_data requires k to be positive integer", {
  # Create temporary directory for testing
  temp_dir <- tempdir()

  expect_error(
    create_data(input = temp_dir, k = -1L),
    "'k' must be a single positive integer"
  )

  expect_error(
    create_data(input = temp_dir, k = 0L),
    "'k' must be a single positive integer"
  )

  expect_error(
    create_data(input = temp_dir, k = c(5L, 6L)),
    "'k' must be a single positive integer"
  )
})

test_that("create_data requires valid alphabet", {
  temp_dir <- tempdir()

  expect_error(
    create_data(input = temp_dir, alphabet = character(0L)),
    "'alphabet' must be a non-empty character vector"
  )

  expect_error(
    create_data(input = temp_dir, alphabet = c("A", "C", "G", "A")),
    "'alphabet' contains duplicate values"
  )

  expect_error(
    create_data(input = temp_dir, alphabet = 123),
    "'alphabet' must be a non-empty character vector"
  )
})

test_that("create_data detects missing FASTA files", {
  # Create empty temporary directory
  temp_dir <- tempdir()

  expect_error(
    create_data(input = temp_dir),
    "No .fasta files found in"
  )
})

test_that("create_data processes single FASTA file correctly", {
  skip_if_not_installed("Biostrings")

  # Create temporary directory with test FASTA file
  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "test_class.fasta")

  # Create a simple FASTA file
  fasta_content <- ">seq1\nACGTACGT\n>seq2\nTGCATGCA\n"
  writeLines(fasta_content, fasta_file)

  # Suppress console output during test
  suppressMessages({
    result <- create_data(input = temp_dir, k = 2L)
  })

  # Verify structure
  expect_type(result, "list")
  expect_length(result, 2L)
  expect_named(result, c("kmers", "metadata"))

  # Verify kmers component
  expect_s3_class(result$kmers, "data.frame")
  expect_true(nrow(result$kmers) >= 2L)
  expect_true("CLASS" %in% colnames(result$kmers))

  # Verify metadata component
  expect_s3_class(result$metadata, "data.frame")
  expect_equal(nrow(result$metadata), nrow(result$kmers))
  expect_named(result$metadata, c("sequence_name", "length", "class"))

  # Cleanup
  unlink(fasta_file)
})

test_that("create_data assigns correct class labels", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()

  # Create multiple FASTA files with different classes
  fasta_file_1 <- file.path(temp_dir, "HIV.fasta")
  fasta_file_2 <- file.path(temp_dir, "HCV.fasta")

  fasta_content_1 <- ">HIV_seq1\nACGTACGT\n"
  fasta_content_2 <- ">HCV_seq1\nTGCATGCA\n"

  writeLines(fasta_content_1, fasta_file_1)
  writeLines(fasta_content_2, fasta_file_2)

  suppressMessages({
    result <- create_data(input = temp_dir, k = 2L)
  })

  # Check class assignment
  expect_true("HIV" %in% result$metadata$class)
  expect_true("HCV" %in% result$metadata$class)

  # Verify CLASS column matches metadata class
  expect_identical(result$kmers$CLASS, result$metadata$class)

  # Cleanup
  unlink(c(fasta_file_1, fasta_file_2))
})

test_that("create_data correctly counts k-mers", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "kmer_test.fasta")

  # Create FASTA with known k-mer composition
  # Sequence: AAAA -> contains 3x 'AA' k-mers (k=2)
  fasta_content <- ">seq_test\nAAAA\n"
  writeLines(fasta_content, fasta_file)

  suppressMessages({
    result <- create_data(input = temp_dir, k = 2L)
  })

  # Check that k-mer counts are numeric
  kmer_cols <- setdiff(colnames(result$kmers), "CLASS")
  expect_true(all(sapply(result$kmers[, kmer_cols], is.numeric)))

  # Check k-mer counts are non-negative
  expect_true(all(as.matrix(result$kmers[, kmer_cols]) >= 0L))

  # Cleanup
  unlink(fasta_file)
})

test_that("create_data handles sequences of different lengths", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "length_test.fasta")

  # Create FASTA with sequences of different lengths
  fasta_content <- paste(
    ">short_seq",
    "ACGT",
    ">long_seq",
    "ACGTACGTACGTACGTACGT",
    sep = "\n"
  )
  writeLines(fasta_content, fasta_file)

  suppressMessages({
    result <- create_data(input = temp_dir, k = 2L)
  })

  # Verify lengths are recorded correctly
  expect_equal(result$metadata$length[1L], 4L)
  expect_equal(result$metadata$length[2L], 20L)

  # Verify both sequences are processed
  expect_equal(nrow(result$metadata), 2L)

  # Cleanup
  unlink(fasta_file)
})

test_that("create_data preserves sequence names", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "names_test.fasta")

  # Create FASTA with specific sequence names
  fasta_content <- ">sequence_001|description\nACGT\n>sequence_002|other_info\nTGCA\n"
  writeLines(fasta_content, fasta_file)

  suppressMessages({
    result <- create_data(input = temp_dir, k = 2L)
  })

  # Verify sequence names are preserved
  expect_true("sequence_001|description" %in% result$metadata$sequence_name)
  expect_true("sequence_002|other_info" %in% result$metadata$sequence_name)

  # Verify row names match sequence names
  expect_identical(rownames(result$kmers), result$metadata$sequence_name)

  # Cleanup
  unlink(fasta_file)
})

test_that("create_data returns correct k-mer count dimensions", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "dim_test.fasta")

  fasta_content <- ">seq1\nACGT\n>seq2\nTGCA\n"
  writeLines(fasta_content, fasta_file)

  suppressMessages({
    result <- create_data(input = temp_dir, k = 2L)
  })

  # For k=2 with DNA alphabet: 4^2 = 16 possible k-mers + 1 CLASS = 17
  expect_equal(ncol(result$kmers), 17L)

  # Cleanup
  unlink(fasta_file)
})

test_that("create_data handles custom alphabet", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "alphabet_test.fasta")

  fasta_content <- ">seq1\nAC\n"
  writeLines(fasta_content, fasta_file)

  suppressMessages({
    result <- create_data(
      input = temp_dir,
      k = 1L,
      alphabet = c("A", "C")
    )
  })

  # For k=1 with alphabet [A,C]: 2^1 = 2 possible k-mers + 1 CLASS = 3
  expect_equal(ncol(result$kmers), 3L)

  # Cleanup
  unlink(fasta_file)
})

test_that("create_data metadata matches kmers", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "consistency_test.fasta")

  fasta_content <- ">seq1\nACGT\n>seq2\nTGCA\n>seq3\nGGGG\n"
  writeLines(fasta_content, fasta_file)

  suppressMessages({
    result <- create_data(input = temp_dir, k = 2L)
  })

  # Number of rows should match
  expect_equal(nrow(result$kmers), nrow(result$metadata))

  # Sequence names should match
  expect_identical(rownames(result$kmers), result$metadata$sequence_name)

  # CLASS column should match class column in metadata
  expect_identical(result$kmers$CLASS, result$metadata$class)

  # Cleanup
  unlink(fasta_file)
})

test_that("create_data warns about large k values", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "warning_test.fasta")

  # Create a long enough sequence for k=11
  # Generate a sequence of at least 11 bp
  long_sequence <- paste(rep("ACGTACGTAC", 2), collapse = "")  # 20 bp
  fasta_content <- paste(">seq1", long_sequence, sep = "\n")
  writeLines(fasta_content, fasta_file)

  # Should warn for k > 10
  expect_warning(
    suppressMessages({
      create_data(input = temp_dir, k = 11L)
    }),
    "Large k value"
  )

  # Cleanup
  unlink(fasta_file)
})

test_that("create_data converts k to integer", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "int_test.fasta")

  # Create sequence long enough for k=6 (need at least 6 bp)
  fasta_content <- ">seq1\nACGTACGTACGT\n"  # 12 bp
  writeLines(fasta_content, fasta_file)

  # Should accept numeric k and convert to integer
  suppressMessages({
    result <- create_data(input = temp_dir, k = 6.9)
  })

  expect_type(result, "list")
  expect_true(nrow(result$kmers) > 0L)

  # Cleanup
  unlink(fasta_file)
})

test_that("create_data case-insensitive for .fasta extension", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()

  # Test with uppercase extension
  fasta_file <- file.path(temp_dir, "case_test.FASTA")
  fasta_content <- ">seq1\nACGT\n"
  writeLines(fasta_content, fasta_file)

  suppressMessages({
    result <- create_data(input = temp_dir, k = 2L)
  })

  expect_type(result, "list")
  expect_true(nrow(result$kmers) > 0L)

  # Cleanup
  unlink(fasta_file)
})

test_that("create_data returns data frames with correct classes", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "class_test.fasta")

  fasta_content <- ">seq1\nACGT\n"
  writeLines(fasta_content, fasta_file)

  suppressMessages({
    result <- create_data(input = temp_dir, k = 2L)
  })

  expect_s3_class(result$kmers, "data.frame")
  expect_s3_class(result$metadata, "data.frame")

  # Cleanup
  unlink(fasta_file)
})

test_that("create_data handles empty sequences gracefully", {
  skip_if_not_installed("Biostrings")

  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "empty_test.fasta")

  # Create FASTA with one valid and verify behavior
  fasta_content <- ">seq1\nACGT\n"
  writeLines(fasta_content, fasta_file)

  suppressMessages({
    result <- create_data(input = temp_dir, k = 2L)
  })

  # Should process successfully
  expect_equal(nrow(result$kmers), 1L)
  expect_equal(nrow(result$metadata), 1L)

  # Cleanup
  unlink(fasta_file)
})

