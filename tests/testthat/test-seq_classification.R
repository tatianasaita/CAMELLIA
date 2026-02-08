# tests/testthat/test-seq_classification.R

library("testthat")

# Helper function to create minimal test data
create_test_fasta <- function(dir, n_seqs = 40) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)

  # Class A sequences
  classA_seqs <- paste0(
    ">seq", 1:n_seqs, "_classA\n",
    vapply(1:n_seqs, function(i) {
      paste(sample(c("A", "T", "C", "G"), 300, replace = TRUE), collapse = "")
    }, character(1))
  )
  writeLines(classA_seqs, file.path(dir, "classA.fasta"))

  # Class B sequences
  classB_seqs <- paste0(
    ">seq", 1:n_seqs, "_classB\n",
    vapply(1:n_seqs, function(i) {
      paste(sample(c("A", "T", "C", "G"), 300, replace = TRUE), collapse = "")
    }, character(1))
  )
  writeLines(classB_seqs, file.path(dir, "classB.fasta"))
}

test_that("seq_classification requires valid input directory", {
  expect_error(
    seq_classification(input_dir = NULL),
    "path"
  )

  expect_error(
    seq_classification(input_dir = "nonexistent_directory_xyz123"),
    "No fasta files found"
  )
})

test_that("seq_classification returns correct structure with external validation", {
  skip_on_cran()  # Skip on CRAN due to computation time

  temp_dir <- tempfile("test_fasta_")
  temp_val_dir <- tempfile("test_val_fasta_")
  on.exit(unlink(c(temp_dir, temp_val_dir), recursive = TRUE), add = TRUE)

  create_test_fasta(temp_dir, n_seqs = 20)
  create_test_fasta(temp_val_dir, n_seqs = 10)

  result <- suppressMessages(suppressWarnings(
    seq_classification(
      input_dir = temp_dir,
      k = 3,
      seq_per_class = 15,
      min_size = 200,
      n_motifs = 10,
      prop_train = 0.7,
      cv_folds = 3,
      external_validation_fasta_dir = temp_val_dir,
      k_for_external_validation = 3,
      verbose = FALSE
    )
  ))

  expect_type(result, "list")

  # Check required components
  expect_true("classification_results" %in% names(result))
  expect_true("best_model" %in% names(result))
  expect_true("validation_predictions" %in% names(result))
  expect_true("validation_actuals" %in% names(result))
  expect_true("validation_accuracy" %in% names(result))
  expect_true("confusion_matrix" %in% names(result))
  expect_true("kmer_analysis" %in% names(result))
  expect_true("processing_time" %in% names(result))
  expect_true("parameters" %in% names(result))
  expect_true("timestamp" %in% names(result))
})

test_that("seq_classification stores parameters correctly", {
  skip_on_cran()

  temp_dir <- tempfile("test_fasta_")
  temp_val_dir <- tempfile("test_val_fasta_")
  on.exit(unlink(c(temp_dir, temp_val_dir), recursive = TRUE), add = TRUE)

  create_test_fasta(temp_dir, n_seqs = 20)
  create_test_fasta(temp_val_dir, n_seqs = 10)

  result <- suppressMessages(suppressWarnings(
    seq_classification(
      input_dir = temp_dir,
      k = 3,
      hom_thresh = 0.75,
      n_motifs = 10,
      seq_per_class = 15,
      min_size = 200,
      cv_folds = 3,
      external_validation_fasta_dir = temp_val_dir,
      k_for_external_validation = 3,
      verbose = FALSE
    )
  ))

  expect_equal(result$parameters$k, 3)
  expect_equal(result$parameters$hom_thresh, 0.75)
  expect_equal(result$parameters$n_motifs, 10)
  expect_equal(result$parameters$min_size, 200)
  expect_equal(result$parameters$dist_method, "euclidean")
  expect_equal(result$parameters$external_validation_fasta_dir, temp_val_dir)
})

test_that("seq_classification best_model is valid", {
  skip_on_cran()

  temp_dir <- tempfile("test_fasta_")
  temp_val_dir <- tempfile("test_val_fasta_")
  on.exit(unlink(c(temp_dir, temp_val_dir), recursive = TRUE), add = TRUE)

  create_test_fasta(temp_dir, n_seqs = 20)
  create_test_fasta(temp_val_dir, n_seqs = 10)

  result <- suppressMessages(suppressWarnings(
    seq_classification(
      input_dir = temp_dir,
      k = 3,
      seq_per_class = 15,
      min_size = 200,
      cv_folds = 3,
      external_validation_fasta_dir = temp_val_dir,
      k_for_external_validation = 3,
      verbose = FALSE
    )
  ))

  expect_true(result$best_model %in% c("Random Forest", "XGBoost"))
  expect_type(result$best_model, "character")
})

test_that("seq_classification validation accuracy is numeric and bounded", {
  skip_on_cran()

  temp_dir <- tempfile("test_fasta_")
  temp_val_dir <- tempfile("test_val_fasta_")
  on.exit(unlink(c(temp_dir, temp_val_dir), recursive = TRUE), add = TRUE)

  create_test_fasta(temp_dir, n_seqs = 20)
  create_test_fasta(temp_val_dir, n_seqs = 10)

  result <- suppressMessages(suppressWarnings(
    seq_classification(
      input_dir = temp_dir,
      k = 3,
      seq_per_class = 15,
      min_size = 200,
      cv_folds = 3,
      external_validation_fasta_dir = temp_val_dir,
      k_for_external_validation = 3,
      verbose = FALSE
    )
  ))

  expect_type(result$validation_accuracy, "double")
  expect_true(result$validation_accuracy >= 0 && result$validation_accuracy <= 1)
})

test_that("seq_classification confusion matrix has correct structure", {
  skip_on_cran()

  temp_dir <- tempfile("test_fasta_")
  temp_val_dir <- tempfile("test_val_fasta_")
  on.exit(unlink(c(temp_dir, temp_val_dir), recursive = TRUE), add = TRUE)

  create_test_fasta(temp_dir, n_seqs = 20)
  create_test_fasta(temp_val_dir, n_seqs = 10)

  result <- suppressMessages(suppressWarnings(
    seq_classification(
      input_dir = temp_dir,
      k = 3,
      seq_per_class = 15,
      min_size = 200,
      cv_folds = 3,
      external_validation_fasta_dir = temp_val_dir,
      k_for_external_validation = 3,
      verbose = FALSE
    )
  ))

  expect_type(result$confusion_matrix, "list")
  expect_true("table" %in% names(result$confusion_matrix))
  expect_true("overall" %in% names(result$confusion_matrix))
  expect_true("Accuracy" %in% names(result$confusion_matrix$overall))
  expect_true("Kappa" %in% names(result$confusion_matrix$overall))
})

 test_that("seq_classification predictions match actuals length", {
  skip_on_cran()

  temp_dir <- tempfile("test_fasta_")
  temp_val_dir <- tempfile("test_val_fasta_")
  on.exit(unlink(c(temp_dir, temp_val_dir), recursive = TRUE), add = TRUE)

  create_test_fasta(temp_dir, n_seqs = 20)
  create_test_fasta(temp_val_dir, n_seqs = 10)

  result <- suppressMessages(suppressWarnings(
    seq_classification(
      input_dir = temp_dir,
      k = 3,
      seq_per_class = 15,
      min_size = 200,
      cv_folds = 3,
      external_validation_fasta_dir = temp_val_dir,
      k_for_external_validation = 3,
      verbose = FALSE
    )
  ))

  expect_equal(
    length(result$validation_predictions),
    length(result$validation_actuals)
  )
  expect_true(length(result$validation_predictions) > 0)
})

test_that("seq_classification verbose parameter controls intermediate results", {
  skip_on_cran()

  temp_dir <- tempfile("test_fasta_")
  temp_val_dir <- tempfile("test_val_fasta_")
  on.exit(unlink(c(temp_dir, temp_val_dir), recursive = TRUE), add = TRUE)

  create_test_fasta(temp_dir, n_seqs = 20)
  create_test_fasta(temp_val_dir, n_seqs = 10)

  result_quiet <- suppressMessages(suppressWarnings(
    seq_classification(
      input_dir = temp_dir,
      k = 3,
      seq_per_class = 15,
      min_size = 200,
      cv_folds = 3,
      external_validation_fasta_dir = temp_val_dir,
      k_for_external_validation = 3,
      verbose = FALSE
    )
  ))

  expect_null(result_quiet$intermediate_results)
})

test_that("seq_classification timestamp is POSIXct", {
  skip_on_cran()

  temp_dir <- tempfile("test_fasta_")
  temp_val_dir <- tempfile("test_val_fasta_")
  on.exit(unlink(c(temp_dir, temp_val_dir), recursive = TRUE), add = TRUE)

  create_test_fasta(temp_dir, n_seqs = 20)
  create_test_fasta(temp_val_dir, n_seqs = 10)

  result <- suppressMessages(suppressWarnings(
    seq_classification(
      input_dir = temp_dir,
      k = 3,
      seq_per_class = 15,
      min_size = 200,
      cv_folds = 3,
      external_validation_fasta_dir = temp_val_dir,
      k_for_external_validation = 3,
      verbose = FALSE
    )
  ))

  expect_s3_class(result$timestamp, "POSIXct")
})

test_that("seq_classification processing_time is valid", {
  skip_on_cran()

  temp_dir <- tempfile("test_fasta_")
  temp_val_dir <- tempfile("test_val_fasta_")
  on.exit(unlink(c(temp_dir, temp_val_dir), recursive = TRUE), add = TRUE)

  create_test_fasta(temp_dir, n_seqs = 20)
  create_test_fasta(temp_val_dir, n_seqs = 10)

  result <- suppressMessages(suppressWarnings(
    seq_classification(
      input_dir = temp_dir,
      k = 3,
      seq_per_class = 15,
      min_size = 200,
      cv_folds = 3,
      external_validation_fasta_dir = temp_val_dir,
      k_for_external_validation = 3,
      verbose = FALSE
    )
  ))

  expect_s3_class(result$processing_time, "difftime")
  expect_true(as.numeric(result$processing_time) > 0)
})

test_that("seq_classification kmer_analysis has expected structure", {
  skip_on_cran()

  temp_dir <- tempfile("test_fasta_")
  temp_val_dir <- tempfile("test_val_fasta_")
  on.exit(unlink(c(temp_dir, temp_val_dir), recursive = TRUE), add = TRUE)

  create_test_fasta(temp_dir, n_seqs = 20)
  create_test_fasta(temp_val_dir, n_seqs = 10)

  result <- suppressMessages(suppressWarnings(
    seq_classification(
      input_dir = temp_dir,
      k = 3,
      seq_per_class = 15,
      min_size = 200,
      cv_folds = 3,
      external_validation_fasta_dir = temp_val_dir,
      k_for_external_validation = 3,
      verbose = FALSE
    )
  ))

  expect_type(result$kmer_analysis, "list")
  expect_true("unique_cluster_motifs" %in% names(result$kmer_analysis))
  expect_true("unique_class_motifs" %in% names(result$kmer_analysis))
})
