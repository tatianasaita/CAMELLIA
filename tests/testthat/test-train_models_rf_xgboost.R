# tests/testthat/test-train_models_rf_xgboost.R
#
# Essential tests for train_models_rf_xgboost function following CRAN standards
# These tests cover:
# - Input validation (dataset_traintest, selected_motifs, prop_train, cv_folds, seed)
# - Data preparation (column standardization, factor levels consistency)
# - Model training (Random Forest and XGBoost)
# - Predictions and confusion matrices
# - S3 methods (print, summary)
# - File I/O operations
#
# Author: [Your Name]
# Date: 2025-11-26

# ===== Test Setup =====

test_that("setup: create mock data for tests", {
  # Create mock classification dataset
  set.seed(123)
  n_samples <- 60
  n_motifs <- 10

  mock_classification <- data.frame(
    matrix(runif(n_samples * n_motifs), nrow = n_samples, ncol = n_motifs)
  )
  colnames(mock_classification) <- paste0("motif", 1:n_motifs)
  mock_classification$class <- factor(rep(c("ClassA", "ClassB", "ClassC"), each = 20))

  # Create mock validation dataset
  n_val <- 30
  mock_validation <- data.frame(
    matrix(runif(n_val * n_motifs), nrow = n_val, ncol = n_motifs)
  )
  colnames(mock_validation) <- paste0("motif", 1:n_motifs)
  mock_validation$class <- factor(rep(c("ClassA", "ClassB", "ClassC"), each = 10))

  # Create mock dataset_traintest
  mock_dataset_traintest <- list(
    classification_dataset = mock_classification,
    validation_dataset = mock_validation
  )

  # Create mock selected_motifs
  mock_selected_motifs <- list(
    ClassA = c("motif1", "motif2", "motif3"),
    ClassB = c("motif4", "motif5", "motif6"),
    ClassC = c("motif7", "motif8", "motif9")
  )

  assign("mock_dataset_traintest", mock_dataset_traintest, envir = .GlobalEnv)
  assign("mock_selected_motifs", mock_selected_motifs, envir = .GlobalEnv)

  expect_true(exists("mock_dataset_traintest"))
  expect_true(exists("mock_selected_motifs"))
})


# ===== Input Validation Tests =====

test_that("rejects invalid dataset_traintest", {
  expect_error(
    train_models_rf_xgboost(NULL, mock_selected_motifs),
    "must be a list"
  )

  expect_error(
    train_models_rf_xgboost(data.frame(), mock_selected_motifs),
    "must be a list.*not a data.frame"
  )

  expect_error(
    train_models_rf_xgboost(list(), mock_selected_motifs),
    "classification_dataset.*validation_dataset"
  )

  expect_error(
    train_models_rf_xgboost(
      list(classification_dataset = data.frame()),
      mock_selected_motifs
    ),
    "validation_dataset"
  )
})


test_that("rejects invalid selected_motifs", {
  expect_error(
    train_models_rf_xgboost(mock_dataset_traintest, NULL),
    "must be a non-empty list"
  )

  expect_error(
    train_models_rf_xgboost(mock_dataset_traintest, list()),
    "must be a non-empty list"
  )

  expect_error(
    train_models_rf_xgboost(mock_dataset_traintest, "motif1"),
    "must be a non-empty list"
  )
})


test_that("rejects invalid prop_train", {
  expect_error(
    train_models_rf_xgboost(mock_dataset_traintest, mock_selected_motifs, prop_train = -0.1),
    "between 0 and 1"
  )

  expect_error(
    train_models_rf_xgboost(mock_dataset_traintest, mock_selected_motifs, prop_train = 0),
    "between 0 and 1"
  )

  expect_error(
    train_models_rf_xgboost(mock_dataset_traintest, mock_selected_motifs, prop_train = 1),
    "between 0 and 1"
  )

  expect_error(
    train_models_rf_xgboost(mock_dataset_traintest, mock_selected_motifs, prop_train = 1.5),
    "between 0 and 1"
  )
})


test_that("rejects invalid cv_folds", {
  expect_error(
    train_models_rf_xgboost(mock_dataset_traintest, mock_selected_motifs, cv_folds = 0),
    "positive integer"
  )

  expect_error(
    train_models_rf_xgboost(mock_dataset_traintest, mock_selected_motifs, cv_folds = -5),
    "positive integer"
  )

  expect_error(
    train_models_rf_xgboost(mock_dataset_traintest, mock_selected_motifs, cv_folds = 2.5),
    "positive integer"
  )
})


test_that("rejects invalid seed", {
  expect_error(
    train_models_rf_xgboost(mock_dataset_traintest, mock_selected_motifs, seed = 123.5),
    "must be an integer"
  )

  expect_error(
    train_models_rf_xgboost(mock_dataset_traintest, mock_selected_motifs, seed = "abc"),
    "must be an integer"
  )
})


test_that("checks for required packages", {
  skip_if_not_installed("caret")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")

  expect_true(requireNamespace("caret", quietly = TRUE))
  expect_true(requireNamespace("randomForest", quietly = TRUE))
  expect_true(requireNamespace("xgboost", quietly = TRUE))
})


# ===== Data Preparation Tests =====

test_that("validates required columns in datasets", {
  bad_dataset <- mock_dataset_traintest
  bad_dataset$classification_dataset <- bad_dataset$classification_dataset[, 1:5]

  expect_error(
    train_models_rf_xgboost(bad_dataset, mock_selected_motifs),
    "Missing columns"
  )
})


test_that("standardizes CLASS to class column name", {
  skip_on_cran()
  skip_if_not_installed("caret")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")

  dataset_with_CLASS <- mock_dataset_traintest
  colnames(dataset_with_CLASS$classification_dataset)[
    colnames(dataset_with_CLASS$classification_dataset) == "class"
  ] <- "CLASS"
  colnames(dataset_with_CLASS$validation_dataset)[
    colnames(dataset_with_CLASS$validation_dataset) == "class"
  ] <- "CLASS"

  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd))

  expect_no_error(
    result <- train_models_rf_xgboost(
      dataset_with_CLASS,
      mock_selected_motifs,
      cv_folds = 2,
      seed = 123
    )
  )
})


# ===== Model Training Tests =====

test_that("trains models successfully with valid inputs", {
  skip_on_cran()
  skip_if_not_installed("caret")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")

  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd))

  result <- train_models_rf_xgboost(
    mock_dataset_traintest,
    mock_selected_motifs,
    prop_train = 0.7,
    cv_folds = 2,
    seed = 123
  )

  expect_s3_class(result, "train_models_rf_xgboost")
  expect_s3_class(result, "list")

  # Check model objects
  expect_true(!is.null(result$model_rf))
  expect_true(!is.null(result$model_xgb))

  # Check data splits
  expect_true(!is.null(result$train_data))
  expect_true(!is.null(result$test_data))
  expect_true(!is.null(result$validation_data))

  # Check predictions
  expect_true(!is.null(result$predictions_test_rf))
  expect_true(!is.null(result$predictions_validation_rf))
  expect_true(!is.null(result$predictions_test_xgb))
  expect_true(!is.null(result$predictions_validation_xgb))

  # Check confusion matrices
  expect_true(!is.null(result$confusion_matrix_test_rf))
  expect_true(!is.null(result$confusion_matrix_validation_rf))
  expect_true(!is.null(result$confusion_matrix_test_xgb))
  expect_true(!is.null(result$confusion_matrix_validation_xgb))

  # Check comparison
  expect_true(!is.null(result$model_comparison))
  expect_equal(nrow(result$model_comparison), 2)
  expect_true(all(c("Model", "Accuracy_Test", "Accuracy_Validation") %in%
                    colnames(result$model_comparison)))
})


test_that("data partition respects prop_train", {
  skip_on_cran()
  skip_if_not_installed("caret")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")

  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd))

  result <- train_models_rf_xgboost(
    mock_dataset_traintest,
    mock_selected_motifs,
    prop_train = 0.8,
    cv_folds = 2,
    seed = 123
  )

  total_classification <- nrow(result$train_data) + nrow(result$test_data)
  actual_prop <- nrow(result$train_data) / total_classification

  # Allow small tolerance for stratification
  expect_true(abs(actual_prop - 0.8) < 0.1)
})


test_that("factor levels are consistent across datasets", {
  skip_on_cran()
  skip_if_not_installed("caret")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")

  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd))

  result <- train_models_rf_xgboost(
    mock_dataset_traintest,
    mock_selected_motifs,
    cv_folds = 2,
    seed = 123
  )

  train_levels <- levels(result$train_data$class)
  test_levels <- levels(result$test_data$class)
  validation_levels <- levels(result$validation_data$class)

  expect_equal(train_levels, test_levels)
  expect_equal(train_levels, validation_levels)
})


test_that("predictions have correct length", {
  skip_on_cran()
  skip_if_not_installed("caret")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")

  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd))

  result <- train_models_rf_xgboost(
    mock_dataset_traintest,
    mock_selected_motifs,
    cv_folds = 2,
    seed = 123
  )

  expect_equal(length(result$predictions_test_rf), nrow(result$test_data))
  expect_equal(length(result$predictions_test_xgb), nrow(result$test_data))
  expect_equal(length(result$predictions_validation_rf), nrow(result$validation_data))
  expect_equal(length(result$predictions_validation_xgb), nrow(result$validation_data))
})


# ===== File I/O Tests =====

test_that("saves required files", {
  skip_on_cran()
  skip_if_not_installed("caret")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")

  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit({
    setwd(old_wd)
    unlink(file.path(temp_dir, "*.RData"))
  })

  result <- train_models_rf_xgboost(
    mock_dataset_traintest,
    mock_selected_motifs,
    cv_folds = 2,
    seed = 123
  )

  expect_true(file.exists("model_rf.RData"))
  expect_true(file.exists("model_xgboost.RData"))
  expect_true(file.exists("model_comparison.RData"))
})


# ===== S3 Methods Tests =====

test_that("print.train_models_rf_xgboost works correctly", {
  skip_on_cran()
  skip_if_not_installed("caret")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")

  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd))

  result <- train_models_rf_xgboost(
    mock_dataset_traintest,
    mock_selected_motifs,
    cv_folds = 2,
    seed = 123
  )

  expect_output(print(result), "Train Models RF/XGBoost Summary")
  expect_output(print(result), "CV folds:")
  expect_output(print(result), "Train proportion:")
  expect_output(print(result), "Performance comparison:")
  expect_output(print(result), "Best model")
})


test_that("summary.train_models_rf_xgboost works correctly", {
  skip_on_cran()
  skip_if_not_installed("caret")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")

  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd))

  result <- train_models_rf_xgboost(
    mock_dataset_traintest,
    mock_selected_motifs,
    cv_folds = 2,
    seed = 123
  )

  expect_output(summary(result), "Detailed Model Summary")
  expect_output(summary(result), "RANDOM FOREST:")
  expect_output(summary(result), "XGBOOST:")
  expect_output(summary(result), "Test Accuracy:")
  expect_output(summary(result), "Validation Accuracy:")
  expect_output(summary(result), "Confusion Matrix")
})


# ===== Edge Cases =====

test_that("handles binary classification", {
  skip_on_cran()
  skip_if_not_installed("caret")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")

  # Create binary dataset
  binary_classification <- mock_dataset_traintest$classification_dataset[
    mock_dataset_traintest$classification_dataset$class %in% c("ClassA", "ClassB"),
  ]
  binary_classification$class <- droplevels(binary_classification$class)

  binary_validation <- mock_dataset_traintest$validation_dataset[
    mock_dataset_traintest$validation_dataset$class %in% c("ClassA", "ClassB"),
  ]
  binary_validation$class <- droplevels(binary_validation$class)

  binary_dataset <- list(
    classification_dataset = binary_classification,
    validation_dataset = binary_validation
  )

  binary_motifs <- list(
    ClassA = c("motif1", "motif2"),
    ClassB = c("motif3", "motif4")
  )

  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd))

  expect_no_error(
    result <- train_models_rf_xgboost(
      binary_dataset,
      binary_motifs,
      cv_folds = 2,
      seed = 123
    )
  )

  expect_equal(nlevels(result$train_data$class), 2)
})


test_that("seed ensures reproducibility", {
  skip_on_cran()
  skip_if_not_installed("caret")
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")

  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd))

  result1 <- train_models_rf_xgboost(
    mock_dataset_traintest,
    mock_selected_motifs,
    cv_folds = 2,
    seed = 456
  )

  result2 <- train_models_rf_xgboost(
    mock_dataset_traintest,
    mock_selected_motifs,
    cv_folds = 2,
    seed = 456
  )

  # Should have same data splits
  expect_equal(rownames(result1$train_data), rownames(result2$train_data))
  expect_equal(rownames(result1$test_data), rownames(result2$test_data))
})


# ===== Cleanup =====

test_that("cleanup: remove mock objects", {
  if (exists("mock_dataset_traintest", envir = .GlobalEnv)) {
    rm("mock_dataset_traintest", envir = .GlobalEnv)
  }
  if (exists("mock_selected_motifs", envir = .GlobalEnv)) {
    rm("mock_selected_motifs", envir = .GlobalEnv)
  }

  expect_false(exists("mock_dataset_traintest", envir = .GlobalEnv))
  expect_false(exists("mock_selected_motifs", envir = .GlobalEnv))
})

