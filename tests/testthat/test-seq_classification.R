# tests/testthat/test-seq_classification.R

library(testthat)

# ===== HELPER FUNCTION =====
create_test_sequences <- function(n_per_class = 30) {
  temp_dir <- file.path(tempdir(), paste0("test_seq_", sample(1:10000, 1)))
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)

  # Create sequences with realistic variation
  set.seed(123)

  for (class_name in c("ClassA", "ClassB")) {
    class_file <- file.path(temp_dir, paste0(class_name, ".fasta"))

    sequences <- character()
    for (i in 1:n_per_class) {
      # Create slightly different sequences for each class
      if (class_name == "ClassA") {
        base_seq <- paste(rep("ATGCATGC", 100), collapse = "")
      } else {
        base_seq <- paste(rep("GCTAGCTA", 100), collapse = "")
      }

      # Add some random variation
      seq_chars <- strsplit(base_seq, "")[[1]]
      n_mutations <- sample(5:15, 1)
      mutation_pos <- sample(1:length(seq_chars), n_mutations)
      seq_chars[mutation_pos] <- sample(c("A", "T", "G", "C"), n_mutations, replace = TRUE)

      sequences <- c(sequences,
                     sprintf(">%s_seq%d", class_name, i),
                     paste(seq_chars, collapse = ""))
    }

    writeLines(sequences, class_file)
  }

  return(temp_dir)
}

# Helper to create valid cluster_result
create_mock_cluster_result <- function() {
  metadata <- data.frame(
    sequence_name = c("seq1", "seq2", "seq3", "seq4", "seq5", "seq6"),
    class = c("A", "A", "A", "B", "B", "B"),
    length = c(500, 900, 1000, 600, 850, 1100),
    stringsAsFactors = FALSE
  )

  kmers <- data.frame(
    kmer1 = c(1, 2, 3, 4, 5, 6),
    kmer2 = c(5, 6, 7, 8, 9, 10),
    kmer3 = c(9, 10, 11, 12, 13, 14)
  )
  rownames(kmers) <- metadata$sequence_name

  list(
    data_result = list(
      metadata = metadata,
      kmers = kmers,
      k = 3
    )
  )
}

# ===== VALIDATION TESTS =====

test_that("select_sample_train_validation validates cluster_result structure", {
  # Missing data_result
  expect_error(
    select_sample_train_validation(
      cluster_result = list(wrong = "structure"),
      k_per_class = 1
    ),
    "'cluster_result' must be a list from cluster_dendrogram\\(\\) with 'data_result'"
  )

  # data_result not a list
  expect_error(
    select_sample_train_validation(
      cluster_result = list(data_result = "not a list"),
      k_per_class = 1
    ),
    "'data_result' must be a list"
  )

  # Missing metadata or kmers
  expect_error(
    select_sample_train_validation(
      cluster_result = list(data_result = list(metadata = data.frame())),
      k_per_class = 1
    ),
    "'data_result' must contain 'metadata' and 'kmers'"
  )
})

test_that("select_sample_train_validation validates metadata structure", {
  expect_error(
    select_sample_train_validation(
      cluster_result = list(
        data_result = list(
          metadata = data.frame(wrong = "columns"),
          kmers = data.frame(k1 = 1)
        )
      ),
      k_per_class = 1
    ),
    "'metadata' missing required columns"
  )
})

test_that("select_sample_train_validation validates k_per_class", {
  cluster_result <- create_mock_cluster_result()

  expect_error(
    select_sample_train_validation(cluster_result, k_per_class = "invalid"),
    "'k_per_class' must be numeric"
  )

  expect_error(
    select_sample_train_validation(cluster_result, k_per_class = -1),
    "'k_per_class' must be positive"
  )

  expect_error(
    select_sample_train_validation(cluster_result, k_per_class = 1.5),
    "'k_per_class' must be an integer"
  )

  expect_error(
    select_sample_train_validation(cluster_result, k_per_class = c(1, 2)),
    "'k_per_class' must be a single value"
  )
})

test_that("select_sample_train_validation validates min_size", {
  cluster_result <- create_mock_cluster_result()

  expect_error(
    select_sample_train_validation(cluster_result, k_per_class = 1, min_size = "invalid"),
    "'min_size' must be numeric"
  )

  expect_error(
    select_sample_train_validation(cluster_result, k_per_class = 1, min_size = -100),
    "'min_size' must be positive"
  )

  expect_error(
    select_sample_train_validation(cluster_result, k_per_class = 1, min_size = 100.5),
    "'min_size' must be an integer"
  )
})

test_that("select_sample_train_validation validates external_validation_fasta_dir", {
  cluster_result <- create_mock_cluster_result()

  expect_error(
    select_sample_train_validation(
      cluster_result,
      k_per_class = 1,
      external_validation_fasta_dir = 123
    ),
    "'external_validation_fasta_dir' must be a character string or NULL"
  )

  expect_error(
    select_sample_train_validation(
      cluster_result,
      k_per_class = 1,
      external_validation_fasta_dir = "/nonexistent/path"
    ),
    "External validation directory does not exist"
  )
})

test_that("select_sample_train_validation validates k_for_external_validation", {
  skip_if_not_installed("seqinr")

  cluster_result <- create_mock_cluster_result()
  temp_dir <- create_test_sequences(5)
  on.exit(unlink(temp_dir, recursive = TRUE))

  expect_error(
    select_sample_train_validation(
      cluster_result,
      k_per_class = 1,
      external_validation_fasta_dir = temp_dir,
      k_for_external_validation = "invalid"
    ),
    "'k_for_external_validation' must be numeric"
  )

  expect_error(
    select_sample_train_validation(
      cluster_result,
      k_per_class = 1,
      external_validation_fasta_dir = temp_dir,
      k_for_external_validation = -1
    ),
    "'k_for_external_validation' must be positive"
  )
})

test_that("select_sample_train_validation validates row count match", {
  metadata <- data.frame(
    sequence_name = c("seq1", "seq2"),
    class = c("A", "B"),
    length = c(900, 1000)
  )

  kmers <- data.frame(kmer1 = c(1, 2, 3))  # Different number of rows

  expect_error(
    select_sample_train_validation(
      cluster_result = list(data_result = list(metadata = metadata, kmers = kmers)),
      k_per_class = 1
    ),
    "'metadata' and 'kmers' must have the same number of rows"
  )
})

test_that("select_sample_train_validation fails when no sequences meet min_size", {
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = c("seq1", "seq2"),
        class = c("A", "B"),
        length = c(500, 600),
        stringsAsFactors = FALSE
      ),
      kmers = data.frame(
        kmer1 = c(1, 2),
        kmer2 = c(3, 4)
      ),
      k = 3
    )
  )
  rownames(cluster_result$data_result$kmers) <- cluster_result$data_result$metadata$sequence_name

  expect_error(
    suppressMessages(
      select_sample_train_validation(cluster_result, k_per_class = 1, min_size = 800)
    ),
    "No sequences available with length >="
  )
})

# ===== FUNCTIONALITY TESTS =====

test_that("select_sample_train_validation returns correct structure", {
  cluster_result <- create_mock_cluster_result()

  result <- suppressMessages(
    select_sample_train_validation(cluster_result, k_per_class = 2, min_size = 800)
  )

  expect_s3_class(result, "sample_selection")
  expect_true(is.list(result))

  expected_names <- c("classification_dataset", "validation_dataset",
                      "classification_metadata", "validation_metadata",
                      "true_labels_classification", "true_labels_validation")
  expect_true(all(expected_names %in% names(result)))

  # Check attributes
  expect_true(!is.null(attr(result, "n_classification")))
  expect_true(!is.null(attr(result, "n_validation")))
  expect_true(!is.null(attr(result, "validation_type")))
  expect_equal(attr(result, "min_size"), 800)
  expect_equal(attr(result, "k_per_class"), 2)
})

test_that("select_sample_train_validation handles min_size correctly", {
  cluster_result <- create_mock_cluster_result()

  result <- suppressMessages(
    select_sample_train_validation(cluster_result, k_per_class = 2, min_size = 800)
  )

  # All classification sequences should have length >= 800
  expect_true(all(result$classification_metadata$length >= 800))

  # Should select up to k_per_class sequences per class
  class_counts <- table(result$classification_metadata$class)
  expect_true(all(class_counts <= 2))
})

test_that("select_sample_train_validation creates non-overlapping datasets", {
  cluster_result <- create_mock_cluster_result()

  result <- suppressMessages(
    select_sample_train_validation(cluster_result, k_per_class = 1, min_size = 800)
  )

  # No sequence should appear in both datasets
  overlap <- intersect(
    result$classification_metadata$sequence_name,
    result$validation_metadata$sequence_name
  )
  expect_equal(length(overlap), 0)
})

test_that("select_sample_train_validation handles insufficient sequences per class", {
  cluster_result <- create_mock_cluster_result()

  # Request more sequences than available (after filtering by min_size)
  result <- suppressMessages(
    select_sample_train_validation(cluster_result, k_per_class = 10, min_size = 800)
  )

  # Should select all available sequences per class
  expect_true(nrow(result$classification_metadata) > 0)
  expect_true(nrow(result$classification_metadata) < 20)  # Less than 10 per class * 2 classes
})

test_that("select_sample_train_validation print method works", {
  cluster_result <- create_mock_cluster_result()

  result <- suppressMessages(
    select_sample_train_validation(cluster_result, k_per_class = 1, min_size = 800)
  )

  expect_output(print(result), "Sample Selection Summary")
  expect_output(print(result), "Classification sequences")
  expect_output(print(result), "Validation sequences")
})

test_that("select_sample_train_validation summary method works", {
  cluster_result <- create_mock_cluster_result()

  result <- suppressMessages(
    select_sample_train_validation(cluster_result, k_per_class = 1, min_size = 800)
  )

  expect_output(summary(result), "Detailed Sample Selection Summary")
  expect_output(summary(result), "Classification Dataset")
  expect_output(summary(result), "Validation Dataset")
})

test_that("select_sample_train_validation saves files correctly", {
  cluster_result <- create_mock_cluster_result()

  # Create temporary directory for test files
  old_wd <- getwd()
  temp_test_dir <- file.path(tempdir(), paste0("test_save_", sample(1:10000, 1)))
  dir.create(temp_test_dir, showWarnings = FALSE)
  setwd(temp_test_dir)
  on.exit({
    setwd(old_wd)
    unlink(temp_test_dir, recursive = TRUE)
  })

  result <- suppressMessages(
    select_sample_train_validation(cluster_result, k_per_class = 1, min_size = 800)
  )

  # Check that files were created
  expect_true(file.exists("true_labels_classification.RData"))
  expect_true(file.exists("classification_dataset.RData"))
  expect_true(file.exists("true_labels_validation.RData"))
  expect_true(file.exists("validation_dataset.RData"))
  expect_true(file.exists("backup_checkpoint1.RData"))

  # Verify file contents
  load("true_labels_classification.RData")
  expect_true(exists("true_labels_classification"))
  expect_s3_class(true_labels_classification, "data.frame")
})

test_that("select_sample_train_validation handles empty validation set", {
  # Create scenario where all sequences are used for classification
  cluster_result <- list(
    data_result = list(
      metadata = data.frame(
        sequence_name = c("seq1", "seq2"),
        class = c("A", "B"),
        length = c(900, 1000),
        stringsAsFactors = FALSE
      ),
      kmers = data.frame(
        kmer1 = c(1, 2),
        kmer2 = c(3, 4)
      ),
      k = 3
    )
  )
  rownames(cluster_result$data_result$kmers) <- cluster_result$data_result$metadata$sequence_name

  expect_warning(
    result <- suppressMessages(
      select_sample_train_validation(cluster_result, k_per_class = 5, min_size = 800)
    ),
    "No sequences available for internal validation"
  )

  expect_equal(nrow(result$validation_metadata), 0)
  expect_equal(nrow(result$validation_dataset), 0)
})

