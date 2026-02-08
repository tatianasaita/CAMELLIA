# tests/testthat/test-kmers_in_seq.R

library(testthat)

test_that("kmers_in_seq returns correct structure", {
  # Mock data
  mock_sequences <- c("ATGCGATAGC", "GCTAGCATGC")
  mock_names <- c("seq1", "seq2")

  mock_cluster <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = mock_names,
        class = c("ClassA", "ClassB"),
        length = c(10, 10),
        stringsAsFactors = FALSE
      )
    )
  )

  motifs <- c("ATG", "GCT")

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = mock_cluster,
    sequences = mock_sequences,
    sequence_names = mock_names,
    verbose = FALSE
  )

  # Test return type
  expect_s3_class(result, "kmers_in_seq_result")
  expect_s3_class(result, "data.frame")

  # Test columns
  expect_true(all(c("motif", "sequence_name", "class", "position_start",
                    "position_end", "sequence_length") %in% colnames(result)))

  # Test attributes
  expect_equal(attr(result, "n_train_sequences"), 2)
  expect_equal(attr(result, "n_validation_sequences"), 0)
  expect_equal(attr(result, "n_motifs"), 2)
  expect_false(attr(result, "has_validation"))
})

test_that("kmers_in_seq finds correct motif positions", {
  mock_sequences <- c("ATGCCCATG")
  mock_names <- c("seq1")

  mock_cluster <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = mock_names,
        class = "ClassA",
        length = 9,
        stringsAsFactors = FALSE
      )
    )
  )

  result <- kmers_in_seq(
    motifs = "ATG",
    cluster_result = mock_cluster,
    sequences = mock_sequences,
    sequence_names = mock_names,
    verbose = FALSE
  )

  # Should find 2 occurrences of ATG
  expect_equal(nrow(result), 2)
  expect_equal(result$position_start, c(1, 7))
  expect_equal(result$position_end, c(3, 9))
  expect_equal(result$motif, c("ATG", "ATG"))
})

test_that("kmers_in_seq handles empty results", {
  mock_sequences <- c("AAAAAAAAAA")
  mock_names <- c("seq1")

  mock_cluster <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = mock_names,
        class = "ClassA",
        length = 10,
        stringsAsFactors = FALSE
      )
    )
  )

  result <- kmers_in_seq(
    motifs = "CGT",
    cluster_result = mock_cluster,
    sequences = mock_sequences,
    sequence_names = mock_names,
    verbose = FALSE
  )

  # Should return empty data.frame with correct structure
  expect_equal(nrow(result), 0)
  expect_true(all(c("motif", "sequence_name", "class") %in% colnames(result)))
})

test_that("kmers_in_seq handles multiple motifs", {
  mock_sequences <- c("ATGCGTTAA")
  mock_names <- c("seq1")

  mock_cluster <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = mock_names,
        class = "ClassA",
        length = 9,
        stringsAsFactors = FALSE
      )
    )
  )

  motifs <- c("ATG", "CGT", "TAA")

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = mock_cluster,
    sequences = mock_sequences,
    sequence_names = mock_names,
    verbose = FALSE
  )

  # Should find all three motifs
  expect_equal(nrow(result), 3)
  expect_setequal(unique(result$motif), motifs)
})

test_that("kmers_in_seq removes duplicate motifs", {
  mock_sequences <- c("ATGCCC")
  mock_names <- c("seq1")

  mock_cluster <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = mock_names,
        class = "ClassA",
        length = 6,
        stringsAsFactors = FALSE
      )
    )
  )

  # Duplicate motifs
  motifs <- c("ATG", "ATG", "CCC")

  result <- kmers_in_seq(
    motifs = motifs,
    cluster_result = mock_cluster,
    sequences = mock_sequences,
    sequence_names = mock_names,
    verbose = FALSE
  )

  # Should only find unique motifs
  expect_equal(attr(result, "n_motifs"), 2)
})

test_that("kmers_in_seq handles validation dataset", {
  skip_if_not_installed("stringi")

  mock_train_names <- c("train1")
  mock_val_names <- c("val1")

  mock_cluster <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = mock_train_names,
        class = "ClassA",
        length = 6,
        stringsAsFactors = FALSE
      )
    )
  )

  # Create SEPARATE temporary directories
  temp_base <- tempdir()
  train_dir <- file.path(temp_base, "train_test")
  val_dir <- file.path(temp_base, "val_test")

  dir.create(train_dir, showWarnings = FALSE)
  dir.create(val_dir, showWarnings = FALSE)

  train_fasta <- file.path(train_dir, "train.fasta")
  val_fasta <- file.path(val_dir, "val.fasta")

  writeLines(c(">train1", "ATGCCC"), train_fasta)
  writeLines(c(">val1", "CCCATG"), val_fasta)

  validation_pred <- data.frame(class = "ClassB", stringsAsFactors = FALSE)

  result <- kmers_in_seq(
    motifs = "ATG",
    cluster_result = mock_cluster,
    input_dir = train_dir,
    external_validation_fasta_dir = val_dir,
    validation_predictions = validation_pred,
    verbose = FALSE
  )

  # Should have dataset column
  expect_true("dataset" %in% colnames(result))
  expect_true(attr(result, "has_validation"))
  expect_equal(attr(result, "n_validation_sequences"), 1)

  # Cleanup
  unlink(train_dir, recursive = TRUE)
  unlink(val_dir, recursive = TRUE)
})

test_that("kmers_in_seq preserves class information", {
  mock_sequences <- c("ATGCCC", "TAAGGC")
  mock_names <- c("seq1", "seq2")

  mock_cluster <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = mock_names,
        class = c("ClassA", "ClassB"),
        length = c(6, 6),
        stringsAsFactors = FALSE
      )
    )
  )

  result <- kmers_in_seq(
    motifs = c("ATG", "TAA"),
    cluster_result = mock_cluster,
    sequences = mock_sequences,
    sequence_names = mock_names,
    verbose = FALSE
  )

  # Check class assignment
  expect_equal(result$class[result$sequence_name == "seq1"], "ClassA")
  expect_equal(result$class[result$sequence_name == "seq2"], "ClassB")
})

