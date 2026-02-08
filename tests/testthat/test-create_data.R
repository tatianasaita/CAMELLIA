# tests/testthat/test-create_data.R

library(testthat)

test_that("create_data works with valid input", {
  # Setup: Create temporary directory with test FASTA files
  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "test_class.fasta")
  
  # Write a simple FASTA file
  writeLines(c(
    ">seq1",
    "ACGTACGT",
    ">seq2",
    "TGCATGCA"
  ), fasta_file)
  
  # Test execution
  result <- create_data(
    input = temp_dir,
    k = 3,
    verbose = FALSE
  )
  
  # Assertions
  expect_s3_class(result, "kmer_data")
  expect_type(result, "list")
  expect_named(result, c("kmers", "metadata"))
  
  # Check kmers data frame
  expect_s3_class(result$kmers, "data.frame")
  expect_true("CLASS" %in% colnames(result$kmers))
  expect_equal(nrow(result$kmers), 2)
  expect_equal(ncol(result$kmers), 4^3 + 1) # 64 k-mers + CLASS column
  
  # Check metadata data frame
  expect_s3_class(result$metadata, "data.frame")
  expect_named(result$metadata, c("sequence_name", "length", "class"))
  expect_equal(nrow(result$metadata), 2)
  expect_equal(result$metadata$class, c("test_class", "test_class"))
  expect_equal(result$metadata$length, c(8, 8))
  
  # Cleanup
  unlink(fasta_file)
})

test_that("create_data fails with no FASTA files", {
  temp_dir <- tempdir()
  
  expect_error(
    create_data(input = temp_dir, k = 3, verbose = FALSE),
    "No fasta files found"
  )
})

test_that("create_data handles custom k and alphabet", {
  temp_dir <- tempdir()
  fasta_file <- file.path(temp_dir, "test.fasta")
  
  writeLines(c(">seq1", "ACGT"), fasta_file)
  
  result <- create_data(
    input = temp_dir,
    k = 2,
    alphabet = c("A", "C"),
    verbose = FALSE
  )
  
  expect_equal(ncol(result$kmers), 2^2 + 1) # 4 k-mers + CLASS
  
  unlink(fasta_file)
})
