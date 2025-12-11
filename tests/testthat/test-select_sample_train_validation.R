# tests/testthat/test-select_sample_train_validation.R

# ===== HELPER FUNCTIONS =====

#' Create mock cluster result for testing
#' @keywords internal
create_mock_cluster_result <- function(n_seqs = 20, n_classes = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  metadata <- data.frame(
    sequence_name = paste0("seq", 1:n_seqs),
    length = sample(800:1200, n_seqs, replace = TRUE),
    class = rep(paste0("Class", LETTERS[1:n_classes]), each = n_seqs / n_classes),
    stringsAsFactors = FALSE
  )

  kmers <- as.data.frame(matrix(
    runif(n_seqs * 10),
    nrow = n_seqs,
    ncol = 10
  ))
  colnames(kmers) <- paste0("kmer", 1:10)
  rownames(kmers) <- metadata$sequence_name

  list(
    data_result = list(
      metadata = metadata,
      kmers = kmers
    )
  )
}


# ===== INPUT VALIDATION TESTS =====

test_that("rejects invalid cluster_result", {
  expect_error(
    select_sample_train_validation(
      cluster_result = list(),
      k_per_class = 3,
      min_size = 800
    ),
    "'cluster_result' must be a list from cluster_dendrogram\\(\\) with 'data_result'"
  )

  expect_error(
    select_sample_train_validation(
      cluster_result = list(data_result = list()),
      k_per_class = 3,
      min_size = 800
    ),
    "'data_result' must contain 'metadata' and 'kmers'"
  )

  expect_error(
    select_sample_train_validation(
      cluster_result = list(data_result = "not a list"),
      k_per_class = 3,
      min_size = 800
    ),
    "'data_result' must be a list"
  )
})

test_that("rejects invalid k_per_class", {
  mock_data <- create_mock_cluster_result(seed = 123)

  expect_error(
    select_sample_train_validation(
      mock_data,
      k_per_class = -1,
      min_size = 800
    ),
    "'k_per_class' must be positive"
  )

  expect_error(
    select_sample_train_validation(
      mock_data,
      k_per_class = 2.5,
      min_size = 800
    ),
    "'k_per_class' must be an integer"
  )

  expect_error(
    select_sample_train_validation(
      mock_data,
      k_per_class = 0,
      min_size = 800
    ),
    "'k_per_class' must be positive"
  )

  expect_error(
    select_sample_train_validation(
      mock_data,
      k_per_class = "invalid",
      min_size = 800
    ),
    "'k_per_class' must be numeric"
  )

  expect_error(
    select_sample_train_validation(
      mock_data,
      k_per_class = c(1, 2),
      min_size = 800
    ),
    "'k_per_class' must be a single value"
  )

  # NA numérico - testa a checagem específica de NA
  expect_error(
    select_sample_train_validation(
      mock_data,
      k_per_class = NA_real_,  # NA numérico explícito
      min_size = 800
    ),
    "'k_per_class' cannot be NA"
  )
})

test_that("rejects invalid min_size", {
  mock_data <- create_mock_cluster_result(seed = 123)

  expect_error(
    select_sample_train_validation(
      mock_data,
      k_per_class = 3,
      min_size = -100
    ),
    "'min_size' must be positive"
  )

  expect_error(
    select_sample_train_validation(
      mock_data,
      k_per_class = 3,
      min_size = 1.5
    ),
    "'min_size' must be an integer"
  )

  expect_error(
    select_sample_train_validation(
      mock_data,
      k_per_class = 3,
      min_size = 0
    ),
    "'min_size' must be positive"
  )

  expect_error(
    select_sample_train_validation(
      mock_data,
      k_per_class = 3,
      min_size = "invalid"
    ),
    "'min_size' must be numeric"
  )
})

test_that("rejects invalid external_validation_fasta_dir", {
  mock_data <- create_mock_cluster_result(seed = 123)

  expect_error(
    select_sample_train_validation(
      cluster_result = mock_data,
      k_per_class = 3,
      min_size = 800,
      external_validation_fasta_dir = "/nonexistent/path"
    ),
    "External validation directory does not exist"
  )

  expect_error(
    select_sample_train_validation(
      cluster_result = mock_data,
      k_per_class = 3,
      min_size = 800,
      external_validation_fasta_dir = 123
    ),
    "'external_validation_fasta_dir' must be a character string or NULL"
  )
})

test_that("rejects invalid k_for_external_validation", {
  skip_if_not_installed("seqinr")

  mock_data <- create_mock_cluster_result(seed = 123)

  # Create temp directory with dummy fasta
  temp_dir <- file.path(tempdir(), paste0("test_", as.numeric(Sys.time())))
  dir.create(temp_dir, recursive = TRUE)
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  expect_error(
    select_sample_train_validation(
      cluster_result = mock_data,
      k_per_class = 3,
      min_size = 800,
      external_validation_fasta_dir = temp_dir,
      k_for_external_validation = "invalid"
    ),
    "'k_for_external_validation' must be numeric"
  )

  expect_error(
    select_sample_train_validation(
      cluster_result = mock_data,
      k_per_class = 3,
      min_size = 800,
      external_validation_fasta_dir = temp_dir,
      k_for_external_validation = -1
    ),
    "'k_for_external_validation' must be positive"
  )

  expect_error(
    select_sample_train_validation(
      cluster_result = mock_data,
      k_per_class = 3,
      min_size = 800,
      external_validation_fasta_dir = temp_dir,
      k_for_external_validation = 2.5
    ),
    "'k_for_external_validation' must be an integer"
  )
})


# ===== CORE FUNCTIONALITY TESTS =====

test_that("selects classification dataset correctly with internal validation", {
  # Create fresh mock data
  mock_data <- create_mock_cluster_result(n_seqs = 20, n_classes = 2, seed = 456)

  # Create isolated test directory
  test_dir <- file.path(tempdir(), paste0("test_", as.numeric(Sys.time())))
  dir.create(test_dir, recursive = TRUE)
  old_wd <- getwd()
  setwd(test_dir)

  on.exit({
    setwd(old_wd)
    unlink(test_dir, recursive = TRUE)
  }, add = TRUE)

  result <- suppressMessages(
    select_sample_train_validation(
      mock_data,
      k_per_class = 3,
      min_size = 800
    )
  )

  # Check S3 class
  expect_s3_class(result, "sample_selection")
  expect_true(inherits(result, "list"))

  # Check metadata
  expect_true(!is.null(result$classification_metadata))
  expect_equal(nrow(result$classification_metadata), 6)  # 3 per class * 2 classes
  expect_true(all(c("sequence_name", "length", "class") %in%
                    colnames(result$classification_metadata)))

  expect_true(!is.null(result$validation_metadata))
  expect_gt(nrow(result$validation_metadata), 0)

  # Check datasets
  expect_equal(nrow(result$classification_dataset),
               nrow(result$classification_metadata))
  expect_equal(nrow(result$validation_dataset),
               nrow(result$validation_metadata))

  # Check attributes
  expect_equal(attr(result, "n_classification"), 6)
  expect_equal(attr(result, "validation_type"), "internal")
  expect_equal(attr(result, "min_size"), 800)
  expect_equal(attr(result, "k_per_class"), 3)

  # Check no overlap
  overlap <- intersect(result$classification_metadata$sequence_name,
                       result$validation_metadata$sequence_name)
  expect_equal(length(overlap), 0)
})

test_that("applies min_size filter only to classification dataset", {
  # Create data with varied lengths
  metadata <- data.frame(
    sequence_name = paste0("seq", 1:20),
    length = c(rep(1000, 10), rep(500, 10)),  # Half below 800
    class = rep(c("ClassA", "ClassB"), each = 10),
    stringsAsFactors = FALSE
  )

  kmers <- as.data.frame(matrix(runif(200), nrow = 20, ncol = 10))
  colnames(kmers) <- paste0("kmer", 1:10)
  rownames(kmers) <- metadata$sequence_name

  cluster_result <- list(
    data_result = list(metadata = metadata, kmers = kmers)
  )

  # Create isolated test directory
  test_dir <- file.path(tempdir(), paste0("test_", as.numeric(Sys.time())))
  dir.create(test_dir, recursive = TRUE)
  old_wd <- getwd()
  setwd(test_dir)

  on.exit({
    setwd(old_wd)
    unlink(test_dir, recursive = TRUE)
  }, add = TRUE)

  result <- suppressMessages(
    select_sample_train_validation(
      cluster_result,
      k_per_class = 3,
      min_size = 800
    )
  )

  # Classification should have only length >= 800
  expect_true(all(result$classification_metadata$length >= 800))

  # Validation can have any length
  expect_true(any(result$validation_metadata$length < 800))
})

test_that("handles insufficient sequences gracefully", {
  mock_data <- create_mock_cluster_result(n_seqs = 10, n_classes = 2, seed = 789)

  # Create isolated test directory
  test_dir <- file.path(tempdir(), paste0("test_", as.numeric(Sys.time())))
  dir.create(test_dir, recursive = TRUE)
  old_wd <- getwd()
  setwd(test_dir)

  on.exit({
    setwd(old_wd)
    unlink(test_dir, recursive = TRUE)
  }, add = TRUE)

  # Suppress both messages AND warnings
  result <- suppressWarnings(suppressMessages(
    select_sample_train_validation(
      mock_data,
      k_per_class = 50,  # More than available
      min_size = 800
    )
  ))

  # Should select all available sequences per class
  expect_gt(nrow(result$classification_metadata), 0)
  expect_lte(nrow(result$classification_metadata), 10)  # Max available

  # Each class should have at most 5 sequences
  class_counts <- table(result$classification_metadata$class)
  expect_true(all(class_counts <= 5))
})


test_that("k-mer dataset matches metadata order", {
  mock_data <- create_mock_cluster_result(n_seqs = 20, n_classes = 2, seed = 111)

  # Create isolated test directory
  test_dir <- file.path(tempdir(), paste0("test_", as.numeric(Sys.time())))
  dir.create(test_dir, recursive = TRUE)
  old_wd <- getwd()
  setwd(test_dir)

  on.exit({
    setwd(old_wd)
    unlink(test_dir, recursive = TRUE)
  }, add = TRUE)

  result <- suppressMessages(
    select_sample_train_validation(
      mock_data,
      k_per_class = 3,
      min_size = 800
    )
  )

  # Row names should match sequence names
  expect_equal(
    rownames(result$classification_dataset),
    result$classification_metadata$sequence_name
  )

  expect_equal(
    rownames(result$validation_dataset),
    result$validation_metadata$sequence_name
  )
})


# ===== FILE SAVING TESTS =====

test_that("saves required files", {
  mock_data <- create_mock_cluster_result(n_seqs = 20, n_classes = 2, seed = 222)

  # Create isolated test directory
  test_dir <- file.path(tempdir(), paste0("test_", as.numeric(Sys.time())))
  dir.create(test_dir, recursive = TRUE)
  old_wd <- getwd()
  setwd(test_dir)

  files_to_remove <- c(
    "true_labels_classification.RData",
    "classification_dataset.RData",
    "true_labels_validation.RData",
    "validation_dataset.RData",
    "backup_checkpoint1.RData"
  )

  on.exit({
    file.remove(files_to_remove[file.exists(files_to_remove)])
    setwd(old_wd)
    unlink(test_dir, recursive = TRUE)
  }, add = TRUE)

  result <- suppressMessages(
    select_sample_train_validation(
      mock_data,
      k_per_class = 3,
      min_size = 800
    )
  )

  # Check files exist
  expect_true(file.exists("true_labels_classification.RData"))
  expect_true(file.exists("classification_dataset.RData"))
  expect_true(file.exists("true_labels_validation.RData"))
  expect_true(file.exists("validation_dataset.RData"))
  expect_true(file.exists("backup_checkpoint1.RData"))

  # Verify backup content
  load("backup_checkpoint1.RData")
  expect_true(exists("backup_data"))
  expect_true(all(c("true_labels_classification", "classification_dataset",
                    "true_labels_validation", "validation_dataset") %in%
                    names(backup_data)))

  # Verify loaded data structure
  load("true_labels_classification.RData")
  expect_s3_class(true_labels_classification, "data.frame")
  expect_true(all(c("sequence_name", "length", "class") %in% colnames(true_labels_classification)))
})


# ===== S3 METHODS TESTS =====

test_that("print.sample_selection works correctly", {
  mock_data <- create_mock_cluster_result(n_seqs = 20, n_classes = 2, seed = 333)

  # Create isolated test directory
  test_dir <- file.path(tempdir(), paste0("test_", as.numeric(Sys.time())))
  dir.create(test_dir, recursive = TRUE)
  old_wd <- getwd()
  setwd(test_dir)

  on.exit({
    setwd(old_wd)
    unlink(test_dir, recursive = TRUE)
  }, add = TRUE)

  result <- suppressMessages(
    select_sample_train_validation(
      mock_data,
      k_per_class = 3,
      min_size = 800
    )
  )

  expect_output(print(result), "Sample Selection Summary")
  expect_output(print(result), "Classification sequences")
  expect_output(print(result), "Validation sequences")
  expect_output(print(result), "internal")
  expect_output(print(result), "Class distribution")
})

test_that("summary.sample_selection works correctly", {
  mock_data <- create_mock_cluster_result(n_seqs = 20, n_classes = 2, seed = 444)

  # Create isolated test directory
  test_dir <- file.path(tempdir(), paste0("test_", as.numeric(Sys.time())))
  dir.create(test_dir, recursive = TRUE)
  old_wd <- getwd()
  setwd(test_dir)

  on.exit({
    setwd(old_wd)
    unlink(test_dir, recursive = TRUE)
  }, add = TRUE)

  result <- suppressMessages(
    select_sample_train_validation(
      mock_data,
      k_per_class = 3,
      min_size = 800
    )
  )

  expect_output(summary(result), "Detailed Sample Selection Summary")
  expect_output(summary(result), "Classification Dataset")
  expect_output(summary(result), "Validation Dataset")
  expect_output(summary(result), "Sequence Length Statistics")
})


# ===== EDGE CASES =====

test_that("handles single class data", {
  single_class_data <- create_mock_cluster_result(n_seqs = 10, n_classes = 1, seed = 555)

  # Create isolated test directory
  test_dir <- file.path(tempdir(), paste0("test_", as.numeric(Sys.time())))
  dir.create(test_dir, recursive = TRUE)
  old_wd <- getwd()
  setwd(test_dir)

  on.exit({
    setwd(old_wd)
    unlink(test_dir, recursive = TRUE)
  }, add = TRUE)

  result <- suppressMessages(
    select_sample_train_validation(
      single_class_data,
      k_per_class = 3,
      min_size = 800
    )
  )

  expect_equal(length(unique(result$classification_metadata$class)), 1)
  expect_equal(nrow(result$classification_metadata), 3)
})

test_that("handles case with no validation sequences remaining", {
  small_data <- create_mock_cluster_result(n_seqs = 4, n_classes = 2, seed = 662)

  # Create isolated test directory
  test_dir <- file.path(tempdir(), paste0("test_", as.numeric(Sys.time())))
  dir.create(test_dir, recursive = TRUE)
  old_wd <- getwd()
  setwd(test_dir)

  on.exit({
    setwd(old_wd)
    unlink(test_dir, recursive = TRUE)
  }, add = TRUE)

  result <- suppressWarnings(suppressMessages(
    select_sample_train_validation(
      small_data,
      k_per_class = 2,
      min_size = 800
    )
  ))

  expect_equal(nrow(result$classification_metadata), 4)
  expect_equal(attr(result, "n_validation"), 0)
  expect_equal(nrow(result$validation_metadata), 0)
})

test_that("returns correct aliases for backward compatibility", {
  mock_data <- create_mock_cluster_result(n_seqs = 20, n_classes = 2, seed = 777)

  test_dir <- file.path(tempdir(), paste0("test_", as.numeric(Sys.time())))
  dir.create(test_dir, recursive = TRUE)
  old_wd <- getwd()
  setwd(test_dir)

  on.exit({
    setwd(old_wd)
    unlink(test_dir, recursive = TRUE)
  }, add = TRUE)

  result <- suppressMessages(
    select_sample_train_validation(mock_data, k_per_class = 3, min_size = 800)
  )

  # Check aliases exist and match
  expect_identical(result$true_labels_classification, result$classification_metadata)
  expect_identical(result$true_labels_validation, result$validation_metadata)
})

