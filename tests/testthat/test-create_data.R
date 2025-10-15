test_that("create_data processes fasta files correctly", {
  # Create temporary directory with test fasta files
  temp_dir <- tempdir()
  test_input <- file.path(temp_dir, "test_input")
  dir.create(test_input, showWarnings = FALSE)

  # Create test fasta file 1
  fasta1_path <- file.path(test_input, "class1.fasta")
  writeLines(c(
    ">seq1",
    "ACGTACGT",
    ">seq2",
    "AAACCCGGGTTT"
  ), fasta1_path)

  # Create test fasta file 2
  fasta2_path <- file.path(test_input, "class2.fasta")
  writeLines(c(
    ">seq3",
    "TTTTGGGG"
  ), fasta2_path)

  # Run function
  result <- create_data(
    input = test_input,
    output = file.path(temp_dir, "output.csv"),
    k = 2
  )

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("kmers", "metadata"))

  # Check metadata
  expect_s3_class(result$metadata, "data.frame")
  expect_equal(nrow(result$metadata), 3)  # 3 sequences total
  expect_equal(result$metadata$class, c("class1", "class1", "class2"))
  expect_equal(result$metadata$sequence_name, c("seq1", "seq2", "seq3"))

  # Check kmer matrix
  expect_equal(nrow(result$kmers), 3)  # 3 sequences
  expect_equal(ncol(result$kmers), 17)  # 16 2-mers + CLASS column
  expect_true("CLASS" %in% colnames(result$kmers))

  # Clean up
  unlink(test_input, recursive = TRUE)
})

test_that("create_data stops when no fasta files found", {
  # Create empty directory
  temp_dir <- tempdir()
  empty_dir <- file.path(temp_dir, "empty_test")
  dir.create(empty_dir, showWarnings = FALSE)

  # Should throw error
  expect_error(
    create_data(empty_dir, "output.csv"),
    "No .fasta files found"
  )

  # Clean up
  unlink(empty_dir, recursive = TRUE)
})

test_that("create_data works with different k values", {
  # Create temporary directory with test fasta file
  temp_dir <- tempdir()
  test_input <- file.path(temp_dir, "test_k_input")
  dir.create(test_input, showWarnings = FALSE)

  fasta_path <- file.path(test_input, "test.fasta")
  writeLines(c(
    ">seq1",
    "ACGTACGTACGT"
  ), fasta_path)

  # Test with k=3
  result_k3 <- create_data(test_input, "output.csv", k = 3)
  expect_equal(ncol(result_k3$kmers), 65)  # 64 3-mers + CLASS

  # Test with k=1
  result_k1 <- create_data(test_input, "output.csv", k = 1)
  expect_equal(ncol(result_k1$kmers), 5)  # 4 1-mers + CLASS

  # Clean up
  unlink(test_input, recursive = TRUE)
})

test_that("create_data handles multiple sequences per file", {
  temp_dir <- tempdir()
  test_input <- file.path(temp_dir, "test_multi_input")
  dir.create(test_input, showWarnings = FALSE)

  fasta_path <- file.path(test_input, "multi.fasta")
  writeLines(c(
    ">seq1",
    "AAAA",
    ">seq2",
    "CCCC",
    ">seq3",
    "GGGG",
    ">seq4",
    "TTTT"
  ), fasta_path)

  result <- create_data(test_input, "output.csv", k = 2)

  expect_equal(nrow(result$metadata), 4)
  expect_equal(nrow(result$kmers), 4)
  expect_true(all(result$metadata$class == "multi"))

  # Clean up
  unlink(test_input, recursive = TRUE)
})

test_that("create_data assigns correct class names from filenames", {
  temp_dir <- tempdir()
  test_input <- file.path(temp_dir, "test_class_input")
  dir.create(test_input, showWarnings = FALSE)

  # Create files with specific names
  writeLines(c(">s1", "ACGT"), file.path(test_input, "positive.fasta"))
  writeLines(c(">s2", "TGCA"), file.path(test_input, "negative.fasta"))

  result <- create_data(test_input, "output.csv", k = 2)

  classes <- unique(result$metadata$class)
  expect_true("positive" %in% classes)
  expect_true("negative" %in% classes)

  # Clean up
  unlink(test_input, recursive = TRUE)
})
