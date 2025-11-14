# tests/testthat/test_train_models_rf_xgboost.R

library(testthat)

# ============================================================================
# Setup: Create test data
# ============================================================================

create_test_data <- function() {
  set.seed(123)

  # Create synthetic classification dataset
  n_samples <- 100
  n_motifs <- 20

  # Generate random motif data (0s and 1s)
  motif_data <- matrix(
    sample(0:1, n_samples * n_motifs, replace = TRUE),
    nrow = n_samples,
    ncol = n_motifs
  )

  # Create motif names
  motif_names <- paste0("MOTIF_", seq_len(n_motifs))
  colnames(motif_data) <- motif_names

  # Create classification dataset
  classification_df <- as.data.frame(motif_data)
  classification_df$class <- factor(
    rep(c("ClassA", "ClassB"), length.out = n_samples)
  )

  # Create validation dataset
  validation_df <- as.data.frame(motif_data)
  validation_df$class <- factor(
    rep(c("ClassA", "ClassB"), length.out = n_samples)
  )

  list(
    classification_dataset = classification_df,
    validation_dataset = validation_df
  )
}

# Create selected motifs list
create_selected_motifs <- function() {
  list(
    ClassA = c("MOTIF_1", "MOTIF_2", "MOTIF_3"),
    ClassB = c("MOTIF_4", "MOTIF_5", "MOTIF_6")
  )
}

# ============================================================================
# Test Suite
# ============================================================================

test_that("Input validation: dataset_traintest must be a list", {
  selected_motifs <- create_selected_motifs()

  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = data.frame(x = 1),
      selected_motifs = selected_motifs
    ),
    "must be a list",
    fixed = FALSE
  )
})

test_that("Input validation: dataset_traintest must contain 'classification_dataset'", {
  selected_motifs <- create_selected_motifs()
  invalid_data <- list(validation_dataset = data.frame(x = 1))

  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = invalid_data,
      selected_motifs = selected_motifs
    ),
    "classification_dataset",
    fixed = FALSE
  )
})

test_that("Input validation: dataset_traintest must contain 'validation_dataset'", {
  selected_motifs <- create_selected_motifs()
  datasets <- create_test_data()

  # Remove validation_dataset
  invalid_data <- list(classification_dataset = datasets$classification_dataset)

  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = invalid_data,
      selected_motifs = selected_motifs
    ),
    "validation_dataset",
    fixed = FALSE
  )
})

test_that("Input validation: selected_motifs must be non-empty list", {
  datasets <- create_test_data()

  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = datasets,
      selected_motifs = list()
    ),
    "non-empty list",
    fixed = FALSE
  )
})

test_that("Input validation: prop_train must be between 0 and 1", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = datasets,
      selected_motifs = selected_motifs,
      prop_train = 1.5
    ),
    "prop_train",
    fixed = FALSE
  )
})

test_that("Input validation: prop_train cannot be 0", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = datasets,
      selected_motifs = selected_motifs,
      prop_train = 0
    ),
    "prop_train",
    fixed = FALSE
  )
})

test_that("Input validation: prop_train cannot be 1", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = datasets,
      selected_motifs = selected_motifs,
      prop_train = 1
    ),
    "prop_train",
    fixed = FALSE
  )
})

test_that("Input validation: cv_folds must be positive integer", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = datasets,
      selected_motifs = selected_motifs,
      cv_folds = -1
    ),
    "cv_folds",
    fixed = FALSE
  )
})

test_that("Input validation: cv_folds must be integer (not decimal)", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = datasets,
      selected_motifs = selected_motifs,
      cv_folds = 2.5
    ),
    "cv_folds",
    fixed = FALSE
  )
})

test_that("Input validation: seed must be integer", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  expect_error(
    train_models_rf_xgboost(
      dataset_traintest = datasets,
      selected_motifs = selected_motifs,
      seed = 123.5
    ),
    "seed",
    fixed = FALSE
  )
})

test_that("Function returns object of class 'train_models_rf_xgboost'", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  expect_s3_class(result, "train_models_rf_xgboost")
})

test_that("Result contains all required list elements", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  required_elements <- c(
    "model_rf", "model_xgb",
    "train_data", "test_data", "validation_data",
    "predictions_test_rf", "predictions_validation_rf",
    "predictions_test_xgb", "predictions_validation_xgb",
    "actuals_test", "actuals_validation",
    "confusion_matrix_test_rf", "confusion_matrix_validation_rf",
    "confusion_matrix_test_xgb", "confusion_matrix_validation_xgb",
    "model_comparison",
    "best_model_test",
    "best_model_validation",
    "motifs_used", "motifs_per_class",
    "cv_folds", "prop_train", "seed",
    "time_rf", "time_xgb"
  )

  for (element in required_elements) {
    expect_true(
      element %in% names(result),
      label = paste("Result contains", element)
    )
  }
})

test_that("Models are trained (caret train objects)", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  expect_s3_class(result$model_rf, "train")
  expect_s3_class(result$model_xgb, "train")
})

test_that("Predictions have correct length and class", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  # Test set predictions
  expect_equal(
    length(result$predictions_test_rf),
    nrow(result$test_data)
  )
  expect_equal(
    length(result$predictions_test_xgb),
    nrow(result$test_data)
  )

  # Validation set predictions
  expect_equal(
    length(result$predictions_validation_rf),
    nrow(result$validation_data)
  )
  expect_equal(
    length(result$predictions_validation_xgb),
    nrow(result$validation_data)
  )
})

test_that("Confusion matrices are confusionMatrix objects", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  expect_s3_class(result$confusion_matrix_test_rf, "confusionMatrix")
  expect_s3_class(result$confusion_matrix_test_xgb, "confusionMatrix")
  expect_s3_class(result$confusion_matrix_validation_rf, "confusionMatrix")
  expect_s3_class(result$confusion_matrix_validation_xgb, "confusionMatrix")
})

test_that("Accuracy values are between 0 and 1", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  acc_test_rf <- result$confusion_matrix_test_rf$overall["Accuracy"]
  acc_test_xgb <- result$confusion_matrix_test_xgb$overall["Accuracy"]
  acc_val_rf <- result$confusion_matrix_validation_rf$overall["Accuracy"]
  acc_val_xgb <- result$confusion_matrix_validation_xgb$overall["Accuracy"]

  expect_true(acc_test_rf >= 0 && acc_test_rf <= 1)
  expect_true(acc_test_xgb >= 0 && acc_test_xgb <= 1)
  expect_true(acc_val_rf >= 0 && acc_val_rf <= 1)
  expect_true(acc_val_xgb >= 0 && acc_val_xgb <= 1)
})

test_that("Comparison dataframe has correct structure", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  comparison <- result$model_comparison

  expect_s3_class(comparison, "data.frame")
  expect_equal(nrow(comparison), 2)
  expect_equal(
    colnames(comparison),
    c("Model", "Accuracy_Test", "Kappa_Test",
      "Accuracy_Validation", "Kappa_Validation", "Training_Time_s")
  )
})

test_that("Best models are identified correctly", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  expect_true(
    result$best_model_test %in% c("Random Forest", "XGBoost")
  )
  expect_true(
    result$best_model_validation %in% c("Random Forest", "XGBoost")
  )
})

test_that("Training times are positive numeric values", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  expect_true(is.numeric(result$time_rf))
  expect_true(is.numeric(result$time_xgb))
  expect_gt(result$time_rf, 0)
  expect_gt(result$time_xgb, 0)
})

test_that("print method works correctly", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  expect_output(
    print(result),
    "Train Models RF XGBoost Object"
  )
})

test_that("summary method works correctly", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  expect_output(
    summary(result),
    "Model Summary"
  )
})

test_that("Reproducibility with same seed", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result1 <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 456
      )

      result2 <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 456
      )
    })
  )

  # Check if predictions are identical
  expect_equal(
    result1$predictions_test_rf,
    result2$predictions_test_rf
  )
})

test_that("Data partition respects prop_train parameter", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  n_total <- nrow(datasets$classification_dataset)
  expected_train <- ceiling(n_total * 0.7)

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  n_train <- nrow(result$train_data)
  n_test <- nrow(result$test_data)

  # Allow small deviation due to rounding
  expect_true(abs(n_train - expected_train) <= 1)
  expect_equal(n_train + n_test, n_total)
})

test_that("Selected motifs are correctly used", {
  datasets <- create_test_data()
  selected_motifs <- create_selected_motifs()

  suppressWarnings(
    suppressMessages({
      result <- train_models_rf_xgboost(
        dataset_traintest = datasets,
        selected_motifs = selected_motifs,
        prop_train = 0.7,
        cv_folds = 2,
        seed = 123
      )
    })
  )

  # Check that motifs_used contains only selected motifs
  expected_motifs <- unique(unlist(selected_motifs))
  expect_equal(
    sort(result$motifs_used),
    sort(expected_motifs)
  )
})
