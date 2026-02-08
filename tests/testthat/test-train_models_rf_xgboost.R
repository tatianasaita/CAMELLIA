# tests/testthat/test-train_models_rf_xgboost.R

library(testthat)

test_that("train_models_rf_xgboost returns correct structure", {
  skip_on_cran()
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("caret")
  
  # Create minimal mock data
  set.seed(123)
  classification_dataset <- data.frame(
    motif1 = rnorm(100),
    motif2 = rnorm(100),
    class = factor(rep(c("A", "B"), each = 50))
  )
  
  validation_dataset <- data.frame(
    motif1 = rnorm(50),
    motif2 = rnorm(50),
    class = factor(rep(c("A", "B"), each = 25))
  )
  
  dataset_traintest <- list(
    classification_dataset = classification_dataset,
    validation_dataset = validation_dataset
  )
  
  selected_motifs <- list(
    A = c("motif1", "motif2"),
    B = c("motif1", "motif2")
  )
  
  # Run function with minimal parameters
  result <- train_models_rf_xgboost(
    dataset_traintest = dataset_traintest,
    selected_motifs = selected_motifs,
    prop_train = 0.7,
    cv_folds = 2,
    verbose = FALSE
  )
  
  # Test class
  expect_s3_class(result, "train_models_rf_xgboost")
  expect_type(result, "list")
  
  # Test essential components
  expect_named(result, c(
    "model_rf", "model_xgb", "train_data", "test_data", "validation_data",
    "predictions_test_rf", "predictions_validation_rf",
    "predictions_test_xgb", "predictions_validation_xgb",
    "actuals_test", "actuals_validation",
    "confusion_matrix_test_rf", "confusion_matrix_validation_rf",
    "confusion_matrix_test_xgb", "confusion_matrix_validation_xgb",
    "model_comparison", "best_model_test", "best_model_validation",
    "motifs_used", "cv_folds", "prop_train"
  ))
  
  # Test model_comparison structure
  expect_s3_class(result$model_comparison, "data.frame")
  expect_equal(nrow(result$model_comparison), 2)
  expect_named(result$model_comparison, c(
    "Model", "Accuracy_Test", "Kappa_Test",
    "Accuracy_Validation", "Kappa_Validation", "Time_seconds"
  ))
})

test_that("train_models_rf_xgboost validates input parameters", {
  skip_on_cran()
  
  # Invalid dataset_traintest structure - missing required columns
  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = list(
        classification_dataset = data.frame(x = 1:10),
        validation_dataset = data.frame(x = 1:10)
      ),
      selected_motifs = list(A = "motif1")
    ),
    "Missing columns"
  )
  
  # Invalid prop_train
  classification_dataset <- data.frame(
    motif1 = rnorm(50),
    class = factor(rep(c("A", "B"), each = 25))
  )
  
  validation_dataset <- data.frame(
    motif1 = rnorm(50),
    class = factor(rep(c("A", "B"), each = 25))
  )
  
  dataset_traintest <- list(
    classification_dataset = classification_dataset,
    validation_dataset = validation_dataset
  )
  
  selected_motifs <- list(A = "motif1", B = "motif1")
  
  # prop_train > 1 should cause error in createDataPartition
  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = dataset_traintest,
      selected_motifs = selected_motifs,
      prop_train = 1.5,
      verbose = FALSE
    )
  )
  
  # prop_train = 0 should cause error
  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = dataset_traintest,
      selected_motifs = selected_motifs,
      prop_train = 0,
      verbose = FALSE
    )
  )
})

test_that("train_models_rf_xgboost handles missing columns", {
  skip_on_cran()
  
  # Dataset without required motif columns
  classification_dataset <- data.frame(
    wrong_motif = rnorm(50),
    class = factor(rep(c("A", "B"), each = 25))
  )
  
  validation_dataset <- data.frame(
    wrong_motif = rnorm(50),
    class = factor(rep(c("A", "B"), each = 25))
  )
  
  dataset_traintest <- list(
    classification_dataset = classification_dataset,
    validation_dataset = validation_dataset
  )
  
  selected_motifs <- list(
    A = "motif1",
    B = "motif1"
  )
  
  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = dataset_traintest,
      selected_motifs = selected_motifs,
      verbose = FALSE
    ),
    "Missing columns"
  )
})

test_that("train_models_rf_xgboost creates output files", {
  skip_on_cran()
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("caret")
  
  # Setup temporary directory
  temp_dir <- tempdir()
  old_wd <- getwd()
  setwd(temp_dir)
  on.exit(setwd(old_wd), add = TRUE)
  
  # Create minimal mock data
  set.seed(456)
  classification_dataset <- data.frame(
    motif1 = rnorm(80),
    motif2 = rnorm(80),
    class = factor(rep(c("A", "B"), each = 40))
  )
  
  validation_dataset <- data.frame(
    motif1 = rnorm(40),
    motif2 = rnorm(40),
    class = factor(rep(c("A", "B"), each = 20))
  )
  
  dataset_traintest <- list(
    classification_dataset = classification_dataset,
    validation_dataset = validation_dataset
  )
  
  selected_motifs <- list(
    A = c("motif1", "motif2"),
    B = c("motif1", "motif2")
  )
  
  # Run function
  result <- train_models_rf_xgboost(
    dataset_traintest = dataset_traintest,
    selected_motifs = selected_motifs,
    prop_train = 0.7,
    cv_folds = 2,
    verbose = FALSE
  )
  
  # Check that files were created
  expect_true(file.exists("model_rf.RData"))
  expect_true(file.exists("model_xgboost.RData"))
  expect_true(file.exists("model_comparison.RData"))
  
  # Clean up
  file.remove("model_rf.RData")
  file.remove("model_xgboost.RData")
  file.remove("model_comparison.RData")
})

test_that("train_models_rf_xgboost handles CLASS column renaming", {
  skip_on_cran()
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("caret")
  
  # Create data with uppercase CLASS column
  set.seed(789)
  classification_dataset <- data.frame(
    motif1 = rnorm(60),
    CLASS = factor(rep(c("A", "B"), each = 30))
  )
  
  validation_dataset <- data.frame(
    motif1 = rnorm(40),
    CLASS = factor(rep(c("A", "B"), each = 20))
  )
  
  dataset_traintest <- list(
    classification_dataset = classification_dataset,
    validation_dataset = validation_dataset
  )
  
  selected_motifs <- list(
    A = "motif1",
    B = "motif1"
  )
  
  # Should not throw error
  expect_silent(
    result <- train_models_rf_xgboost(
      dataset_traintest = dataset_traintest,
      selected_motifs = selected_motifs,
      prop_train = 0.7,
      cv_folds = 2,
      verbose = FALSE
    )
  )
  
  # Check that class column exists in output
  expect_true("class" %in% colnames(result$train_data))
  expect_true("class" %in% colnames(result$test_data))
  expect_true("class" %in% colnames(result$validation_data))
})

test_that("train_models_rf_xgboost verbose parameter works", {
  skip_on_cran()
  skip_if_not_installed("randomForest")
  skip_if_not_installed("xgboost")
  skip_if_not_installed("caret")
  
  set.seed(321)
  classification_dataset <- data.frame(
    motif1 = rnorm(60),
    class = factor(rep(c("A", "B"), each = 30))
  )
  
  validation_dataset <- data.frame(
    motif1 = rnorm(30),
    class = factor(rep(c("A", "B"), each = 15))
  )
  
  dataset_traintest <- list(
    classification_dataset = classification_dataset,
    validation_dataset = validation_dataset
  )
  
  selected_motifs <- list(A = "motif1", B = "motif1")
  
  # Test verbose = FALSE produces no output
  expect_silent(
    train_models_rf_xgboost(
      dataset_traintest = dataset_traintest,
      selected_motifs = selected_motifs,
      cv_folds = 2,
      verbose = FALSE
    )
  )
  
  # Test verbose = TRUE produces output
  expect_output(
    train_models_rf_xgboost(
      dataset_traintest = dataset_traintest,
      selected_motifs = selected_motifs,
      cv_folds = 2,
      verbose = TRUE
    ),
    "Data summary"
  )
})
