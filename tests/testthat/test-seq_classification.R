# tests/testthat/test-seq_classification.R

# ===== HELPER FUNCTION =====

create_test_fasta <- function(dir, n_seqs = 10, seq_length = 1000, n_classes = 2) {
  if (dir.exists(dir)) {
    unlink(dir, recursive = TRUE)
  }
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)

  classes <- paste0("Class", LETTERS[1:n_classes])
  seq_per_class <- n_seqs %/% n_classes

  seq_counter <- 1
  for (class_idx in 1:n_classes) {
    class_name <- classes[class_idx]
    file_path <- file.path(dir, paste0(class_name, ".fasta"))

    for (i in 1:seq_per_class) {
      seq <- paste(sample(c("A", "T", "G", "C"), seq_length, replace = TRUE), collapse = "")
      cat(paste0(">", class_name, "_seq", seq_counter, " ", class_name, "\n"),
          seq, "\n",
          file = file_path,
          append = file.exists(file_path))
      seq_counter <- seq_counter + 1
    }
  }

  invisible(NULL)
}

# ===== TEST 1: PARAMETER VALIDATION =====

test_that("seq_classification validates parameters", {
  skip_if_not_installed("seqinr")

  temp_dir <- file.path(tempdir(), "test_validation")
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  create_test_fasta(temp_dir, n_seqs = 10, n_classes = 2)

  # Invalid k
  expect_error(
    seq_classification(input_dir = temp_dir, k = -1),
    "'k' must be a single positive integer"
  )

  # Invalid hom_thresh
  expect_error(
    seq_classification(input_dir = temp_dir, hom_thresh = 1.5),
    "between 0 and 1"
  )

  # Invalid dist_method
  expect_error(
    seq_classification(input_dir = temp_dir, dist_method = "invalid"),
    "must be one of"
  )

  # Invalid directory
  expect_error(
    seq_classification(input_dir = "/nonexistent/path"),
    "Directory does not exist"
  )
})

# ===== TEST 2: COMPLETE PIPELINE =====

test_that("seq_classification completes full pipeline", {
  skip_if_not_installed("seqinr")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("caret")

  test_id <- format(Sys.time(), "%Y%m%d_%H%M%OS3")
  temp_dir <- file.path(tempdir(), paste0("seq_test_", test_id))
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  # Create 60 sequences (30 per class)
  create_test_fasta(temp_dir, n_seqs = 60, seq_length = 1200, n_classes = 2)

  # Verify 2 files created
  fasta_files <- list.files(temp_dir, pattern = "\\.fasta$")
  expect_equal(length(fasta_files), 2)

  # Run pipeline
  result <- suppressMessages(
    seq_classification(
      input_dir = temp_dir,
      k = 4,
      k_per_class = 15,
      min_size = 800,
      n_motifs = 30,
      prop_train = 0.7,
      cv_folds = 3,
      verbose = FALSE
    )
  )

  # Check result is a list
  expect_type(result, "list")

  # Check essential components exist
  expect_true("classification_results" %in% names(result))
  expect_true("best_model" %in% names(result))
  expect_true("validation_accuracy" %in% names(result))
  expect_true("confusion_matrix" %in% names(result))
  expect_true("parameters" %in% names(result))

  # Validate best_model
  expect_true(result$best_model %in% c("Random Forest", "XGBoost"))

  # Validate accuracy range
  expect_true(result$validation_accuracy >= 0)
  expect_true(result$validation_accuracy <= 1)

  # Check parameters stored correctly
  expect_equal(result$parameters$k, 4)
  expect_equal(result$parameters$k_per_class, 15)
  expect_equal(result$parameters$dist_method, "euclidean")
})

# ===== TEST 3: OUTPUT STRUCTURE =====

test_that("seq_classification returns correct structure", {
  skip_if_not_installed("seqinr")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("caret")

  test_id <- format(Sys.time(), "%Y%m%d_%H%M%OS3")
  temp_dir <- file.path(tempdir(), paste0("seq_struct_", test_id))
  on.exit(unlink(temp_dir, recursive = TRUE), add = TRUE)

  create_test_fasta(temp_dir, n_seqs = 60, seq_length = 1000, n_classes = 2)

  result <- suppressMessages(
    seq_classification(
      input_dir = temp_dir,
      k = 4,
      k_per_class = 15,
      min_size = 500,
      n_motifs = 30,
      cv_folds = 3,
      verbose = FALSE
    )
  )

  # Check all required list elements
  expected_names <- c(
    "classification_results",
    "best_model",
    "validation_predictions",
    "validation_actuals",
    "validation_accuracy",
    "confusion_matrix",
    "model_comparison",
    "kmer_analysis",
    "processing_time",
    "parameters",
    "timestamp"
  )

  for (name in expected_names) {
    expect_true(name %in% names(result), info = paste("Missing:", name))
  }

  # Check data types
  expect_type(result$classification_results, "list")
  expect_type(result$best_model, "character")
  expect_type(result$validation_accuracy, "double")
  expect_s3_class(result$confusion_matrix, "confusionMatrix")
  expect_s3_class(result$processing_time, "difftime")
  expect_s3_class(result$timestamp, "POSIXct")
})

# ===== TEST 4: EXTERNAL VALIDATION =====

test_that("seq_classification handles external validation", {
  skip_if_not_installed("seqinr")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("caret")

  test_id <- format(Sys.time(), "%Y%m%d_%H%M%OS3")
  train_dir <- file.path(tempdir(), paste0("seq_train_", test_id))
  val_dir <- file.path(tempdir(), paste0("seq_val_", test_id))

  on.exit({
    unlink(train_dir, recursive = TRUE)
    unlink(val_dir, recursive = TRUE)
  }, add = TRUE)

  # Create training and validation data
  create_test_fasta(train_dir, n_seqs = 60, seq_length = 1000, n_classes = 2)
  create_test_fasta(val_dir, n_seqs = 30, seq_length = 1000, n_classes = 2)

  # Run with external validation
  result <- suppressMessages(
    seq_classification(
      input_dir = train_dir,
      k = 4,
      k_per_class = 15,
      min_size = 500,
      n_motifs = 30,
      external_validation_fasta_dir = val_dir,
      k_for_external_validation = 4,
      cv_folds = 3,
      verbose = FALSE
    )
  )

  expect_type(result, "list")
  expect_equal(result$parameters$external_validation_fasta_dir, val_dir)
  expect_equal(result$parameters$k_for_external_validation, 4)
})
