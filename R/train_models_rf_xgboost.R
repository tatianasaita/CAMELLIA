#' Train Random Forest and XGBoost Models
#'
#' @param dataset_traintest List from select_sample_train_validation() containing
#'   classification_dataset and validation_dataset.
#' @param selected_motifs List from select_motifs() with character vectors of motifs per class.
#' @param prop_train Numeric. Proportion for training (0-1). Default: 0.7.
#' @param cv_folds Integer. Number of cross-validation folds. Default: 5.
#' @param verbose Logical. If TRUE, prints detailed training information. Default: TRUE.
#'
#' @return Object of class "train_models_rf_xgboost" with:
#'   \itemize{
#'     \item Model_rf: Trained Random Forest model
#'     \item Model_xgb: Trained XGBoost model
#'     \item Model_comparison: Data frame with performance metrics
#'     \item All predictions and confusion matrices
#'   }
#'
#' @details
#' Trains Random Forest (500 trees) and XGBoost models using k-fold cross-validation.
#' The classification dataset is split into training and test sets according to prop_train.
#' Both models are evaluated on test and validation sets. Performance metrics (accuracy, kappa)
#' and confusion matrices are computed. Models are saved as .RData files in the working directory.
#'
#' @note
#' \itemize{
#'   \item Requires .rename_class_column and .validate_columns. See \code{internal-functions.R} for details.
#'   \item S3 methods available. See \code{methods.R} for details.
#' }
#'
#' @examples
#' \dontrun{
#' result_models <- train_models_rf_xgboost(
#'   dataset_traintest = result_sample_internal,
#'   selected_motifs = result_selected_motifs,
#'   prop_train = 0.8,
#'   cv_folds = 10,
#'   verbose = TRUE
#' )
#'}
#'
#' @importFrom caret createDataPartition
#' @importFrom caret trainControl
#' @importFrom caret train
#' @importFrom caret confusionMatrix
#' @importFrom stats predict
#'
#' @export
train_models_rf_xgboost <- function(dataset_traintest,
                                    selected_motifs,
                                    prop_train = 0.7,
                                    cv_folds = 5,
                                    verbose = TRUE) {

# Check required packages
  if (!requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' required for Random Forest training.\n",
         "Install with: install.packages('randomForest')")
  }

  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' required for XGBoost training.\n",
         "Install with: install.packages('xgboost')")
  }

  set.seed(123)

  # Data preparation
  classification_dataset <- dataset_traintest$classification_dataset
  validation_dataset <- dataset_traintest$validation_dataset

  # Standardize column name to lowercase
  classification_dataset <- .rename_class_column(classification_dataset)
  validation_dataset <- .rename_class_column(validation_dataset)

  # Extract motifs
  motifs_all <- unique(unlist(selected_motifs))
  required_columns <- c(motifs_all, "class")

  # Validate columns
  .validate_columns(classification_dataset, required_columns, "classification_dataset")
  .validate_columns(validation_dataset, required_columns, "validation_dataset")

  # Partition data
  train_partition <- caret::createDataPartition(
    classification_dataset$class,
    p = prop_train,
    list = FALSE
  )

  traindata <- classification_dataset[train_partition, required_columns, drop = FALSE]
  testdata <- classification_dataset[-train_partition, required_columns, drop = FALSE]
  validationdata <- validation_dataset[, required_columns, drop = FALSE]

  # Convert to factors with consistent levels
  traindata$class <- factor(traindata$class)
  testdata$class <- factor(testdata$class, levels = levels(traindata$class))
  validationdata$class <- factor(validationdata$class, levels = levels(traindata$class))

  if (verbose) {
    cat("Data summary\n")
    cat("Training:", nrow(traindata), "| Test:", nrow(testdata),
        "| Validation:", nrow(validationdata), "\n")
    cat("Motifs:", length(motifs_all), "| Classes:", nlevels(traindata$class), "\n\n")
  }

  # CV control
  control <- caret::trainControl(
    method = "cv",
    number = cv_folds,
    savePredictions = "final",
    classProbs = TRUE
  )

  # Train Randon Forest
  start_time_rf <- Sys.time()

  model_rf <- suppressWarnings(
    suppressMessages(
      caret::train(
        class ~ .,
        data = traindata,
        method = "rf",
        trControl = control,
        ntree = 500,
        importance = TRUE,
        verbose = FALSE
      )
    )
  )

  time_rf <- as.numeric(difftime(Sys.time(), start_time_rf, units = "secs"))

  # Predictions
  pred_test_rf <- predict(model_rf, newdata = testdata)
  pred_validation_rf <- predict(model_rf, newdata = validationdata)

  cm_test_rf <- caret::confusionMatrix(pred_test_rf, testdata$class)
  cm_validation_rf <- caret::confusionMatrix(pred_validation_rf, validationdata$class)

  if (verbose) {
    cat("Training Random Forest\n")
    cat("Completed in", round(time_rf, 2), "seconds\n")
    cat("Test Accuracy:", round(cm_test_rf$overall["Accuracy"], 4), "\n")
    cat("Validation Accuracy:", round(cm_validation_rf$overall["Accuracy"], 4), "\n\n")
  }

  # Training XGBoost
  start_time_xgb <- Sys.time()

  invisible(
    capture.output(
      {
        model_xgb <- suppressWarnings(
          suppressMessages(
            caret::train(
              class ~ .,
              data = traindata,
              method = "xgbTree",
              trControl = control,
              verbose = FALSE
            )
          )
        )
      },
      file = NULL
    )
  )

  time_xgb <- as.numeric(difftime(Sys.time(), start_time_xgb, units = "secs"))

  # Predictions
  pred_test_xgb <- predict(model_xgb, newdata = testdata)
  pred_validation_xgb <- predict(model_xgb, newdata = validationdata)

  cm_test_xgb <- caret::confusionMatrix(pred_test_xgb, testdata$class)
  cm_validation_xgb <- caret::confusionMatrix(pred_validation_xgb, validationdata$class)

  if (verbose) {
    cat("Training XGBoost\n")
    cat("Completed in", round(time_xgb, 2), "seconds\n")
    cat("Test Accuracy:", round(cm_test_xgb$overall["Accuracy"], 4), "\n")
    cat("Validation Accuracy:", round(cm_validation_xgb$overall["Accuracy"], 4), "\n\n")
  }

  # Model Comparison
  model_comparison <- data.frame(
    Model = c("Random Forest", "XGBoost"),
    Accuracy_Test = round(c(cm_test_rf$overall["Accuracy"],
                            cm_test_xgb$overall["Accuracy"]), 4),
    Kappa_Test = round(c(cm_test_rf$overall["Kappa"],
                         cm_test_xgb$overall["Kappa"]), 4),
    Accuracy_Validation = round(c(cm_validation_rf$overall["Accuracy"],
                                  cm_validation_xgb$overall["Accuracy"]), 4),
    Kappa_Validation = round(c(cm_validation_rf$overall["Kappa"],
                               cm_validation_xgb$overall["Kappa"]), 4),
    Time_seconds = round(c(time_rf, time_xgb), 2),
    stringsAsFactors = FALSE
  )

  best_test <- if (model_comparison$Accuracy_Test[1] > model_comparison$Accuracy_Test[2]) {
    "Random Forest"
  } else {
    "XGBoost"
  }

  best_validation <- if (model_comparison$Accuracy_Validation[1] > model_comparison$Accuracy_Validation[2]) {
    "Random Forest"
  } else {
    "XGBoost"
  }

  if (verbose) {
    cat("Model comparison\n")
    print(model_comparison)
    cat("Best on test:", best_test, "\n")
    cat("Best on validation:", best_validation, "\n\n")
  }

  # Confusion Matrices
  if (verbose) {
    cat("Confusion Matrix Random Forest (Validation set)\n\n")
    print(cm_validation_rf$table)
    cat("\nOverall Statistics:\n")
    cat("  Accuracy:", round(cm_validation_rf$overall["Accuracy"], 4), "\n")
    cat("  Kappa:", round(cm_validation_rf$overall["Kappa"], 4), "\n")
    cat("  95% CI: (", round(cm_validation_rf$overall["AccuracyLower"], 4), ", ",
        round(cm_validation_rf$overall["AccuracyUpper"], 4), ")\n\n")

    cat("Confusion Matrix XGBoost (Validation set)\n\n")
    print(cm_validation_xgb$table)
    cat("\nOverall Statistics:\n")
    cat("  Accuracy:", round(cm_validation_xgb$overall["Accuracy"], 4), "\n")
    cat("  Kappa:", round(cm_validation_xgb$overall["Kappa"], 4), "\n")
    cat("  95% CI: (", round(cm_validation_xgb$overall["AccuracyLower"], 4), ", ",
        round(cm_validation_xgb$overall["AccuracyUpper"], 4), ")\n\n")
  }

  # Save Results
  save(model_rf, file = "model_rf.RData")
  save(model_xgb, file = "model_xgboost.RData")
  save(model_comparison, file = "model_comparison.RData")
  if (verbose) {
    cat("Files saved: model_rf.RData, model_xgboost.RData, model_comparison.RData\n\n")
  }

  # Return result
  result <- structure(
    list(
      model_rf = model_rf,
      model_xgb = model_xgb,
      train_data = traindata,
      test_data = testdata,
      validation_data = validationdata,
      predictions_test_rf = pred_test_rf,
      predictions_validation_rf = pred_validation_rf,
      predictions_test_xgb = pred_test_xgb,
      predictions_validation_xgb = pred_validation_xgb,
      actuals_test = testdata$class,
      actuals_validation = validationdata$class,
      confusion_matrix_test_rf = cm_test_rf,
      confusion_matrix_validation_rf = cm_validation_rf,
      confusion_matrix_test_xgb = cm_test_xgb,
      confusion_matrix_validation_xgb = cm_validation_xgb,
      model_comparison = model_comparison,
      best_model_test = best_test,
      best_model_validation = best_validation,
      motifs_used = motifs_all,
      cv_folds = cv_folds,
      prop_train = prop_train
    ),
    class = c("train_models_rf_xgboost", "list")
  )

  return(invisible(result))
}
