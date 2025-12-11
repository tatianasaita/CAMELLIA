#' Train Random Forest and XGBoost Models
#'
#' @param dataset_traintest List from select_sample_train_validation() containing
#'   classification_dataset and validation_dataset.
#' @param selected_motifs List from select_motifs() with character vectors of motifs per class.
#' @param prop_train Numeric. Proportion for training (0-1). Default: 0.7.
#' @param cv_folds Integer. Number of cross-validation folds. Default: 5.
#' @param seed Integer. Random seed for reproducibility. Default: 123.
#'
#' @return Object of class "train_models_rf_xgboost" with:
#'   \itemize{
#'     \item model_rf: Trained Random Forest model
#'     \item model_xgb: Trained XGBoost model
#'     \item model_comparison: Data frame with performance metrics
#'     \item All predictions and confusion matrices
#'   }
#'
#' @export
train_models_rf_xgboost <- function(dataset_traintest,
                                    selected_motifs,
                                    prop_train = 0.7,
                                    cv_folds = 5,
                                    seed = 123) {

  # ===== VALIDATION =====

  if (!is.list(dataset_traintest) || is.data.frame(dataset_traintest)) {
    stop("'dataset_traintest' must be a list (not a data.frame)", call. = FALSE)
  }

  required_names <- c("classification_dataset", "validation_dataset")
  if (!all(required_names %in% names(dataset_traintest))) {
    stop("'dataset_traintest' must contain 'classification_dataset' and 'validation_dataset'",
         call. = FALSE)
  }

  if (!is.list(selected_motifs) || length(selected_motifs) == 0) {
    stop("'selected_motifs' must be a non-empty list", call. = FALSE)
  }

  if (!is.numeric(prop_train) || prop_train <= 0 || prop_train >= 1) {
    stop("'prop_train' must be between 0 and 1", call. = FALSE)
  }

  if (!.is_positive_integer(cv_folds)) {
    stop("'cv_folds' must be a positive integer", call. = FALSE)
  }

  if (!is.numeric(seed) || seed != as.integer(seed)) {
    stop("'seed' must be an integer", call. = FALSE)
  }

  # Check required packages
  required_packages <- c("caret", "randomForest", "xgboost")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' required. Install with: install.packages('", pkg, "')",
           call. = FALSE)
    }
  }

  set.seed(seed)

  # ===== DATA PREPARATION =====

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

  cat("\n=== DATA SUMMARY ===\n")
  cat("Training:", nrow(traindata), "| Test:", nrow(testdata),
      "| Validation:", nrow(validationdata), "\n")
  cat("Motifs:", length(motifs_all), "| Classes:", nlevels(traindata$class), "\n\n")

  # CV control
  control <- caret::trainControl(
    method = "cv",
    number = cv_folds,
    savePredictions = "final",
    classProbs = TRUE
  )

  # ===== TRAIN RANDOM FOREST =====

  cat("=== TRAINING RANDOM FOREST ===\n")

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

  cat("Completed in", round(time_rf, 2), "seconds\n")
  cat("Test Accuracy:", round(cm_test_rf$overall["Accuracy"], 4), "\n")
  cat("Validation Accuracy:", round(cm_validation_rf$overall["Accuracy"], 4), "\n\n")

  # ===== TRAIN XGBOOST =====

  cat("=== TRAINING XGBOOST ===\n")

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

  cat("Completed in", round(time_xgb, 2), "seconds\n")
  cat("Test Accuracy:", round(cm_test_xgb$overall["Accuracy"], 4), "\n")
  cat("Validation Accuracy:", round(cm_validation_xgb$overall["Accuracy"], 4), "\n\n")

  # ===== MODEL COMPARISON =====

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

  cat("=== MODEL COMPARISON ===\n")
  print(model_comparison)
  cat("\n")

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

  cat("Best on test:", best_test, "\n")
  cat("Best on validation:", best_validation, "\n\n")

  # ===== CONFUSION MATRICES - VALIDATION SET =====

  cat("=== CONFUSION MATRIX - RANDOM FOREST (Validation Set) ===\n\n")
  print(cm_validation_rf$table)
  cat("\nOverall Statistics:\n")
  cat("  Accuracy:", round(cm_validation_rf$overall["Accuracy"], 4), "\n")
  cat("  Kappa:", round(cm_validation_rf$overall["Kappa"], 4), "\n")
  cat("  95% CI: (", round(cm_validation_rf$overall["AccuracyLower"], 4), ", ",
      round(cm_validation_rf$overall["AccuracyUpper"], 4), ")\n\n")

  cat("=== CONFUSION MATRIX - XGBOOST (Validation Set) ===\n\n")
  print(cm_validation_xgb$table)
  cat("\nOverall Statistics:\n")
  cat("  Accuracy:", round(cm_validation_xgb$overall["Accuracy"], 4), "\n")
  cat("  Kappa:", round(cm_validation_xgb$overall["Kappa"], 4), "\n")
  cat("  95% CI: (", round(cm_validation_xgb$overall["AccuracyLower"], 4), ", ",
      round(cm_validation_xgb$overall["AccuracyUpper"], 4), ")\n\n")

  # ===== SAVE RESULTS =====

  cat("=== SAVING RESULTS ===\n")
  save(model_rf, file = "model_rf.RData")
  save(model_xgb, file = "model_xgboost.RData")
  save(model_comparison, file = "model_comparison.RData")
  cat("Files saved: model_rf.RData, model_xgboost.RData, model_comparison.RData\n\n")

  # ===== RETURN =====

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
      prop_train = prop_train,
      seed = seed
    ),
    class = c("train_models_rf_xgboost", "list")
  )

  return(invisible(result))
}


# ===== INTERNAL FUNCTIONS =====

#' Check if value is a positive integer
#' @keywords internal
.is_positive_integer <- function(x) {
  is.numeric(x) && length(x) == 1L && x > 0 && x == as.integer(x)
}


#' Rename CLASS column to class
#' @keywords internal
.rename_class_column <- function(data) {
  if ("CLASS" %in% colnames(data) && !("class" %in% colnames(data))) {
    colnames(data)[colnames(data) == "CLASS"] <- "class"
  }
  data
}


#' Validate required columns
#' @keywords internal
.validate_columns <- function(data, required_cols, dataset_name) {
  missing <- setdiff(required_cols, colnames(data))
  if (length(missing) > 0) {
    stop("Missing columns in ", dataset_name, ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}


# ===== S3 METHODS =====

#' @export
print.train_models_rf_xgboost <- function(x, ...) {
  cat("\n=== Train Models RF/XGBoost Summary ===\n\n")
  cat("Training configuration:\n")
  cat("  CV folds:", x$cv_folds, "\n")
  cat("  Train proportion:", x$prop_train, "\n")
  cat("  Motifs used:", length(x$motifs_used), "\n")
  cat("  Classes:", nlevels(x$actuals_test), "\n\n")
  cat("Performance comparison:\n")
  print(x$model_comparison)
  cat("\nBest model (test):", x$best_model_test, "\n")
  cat("Best model (validation):", x$best_model_validation, "\n")
  invisible(x)
}


#' @export
summary.train_models_rf_xgboost <- function(object, ...) {
  cat("\n=== Detailed Model Summary ===\n\n")

  cat("RANDOM FOREST:\n")
  cat("  Test Accuracy:",
      round(object$confusion_matrix_test_rf$overall["Accuracy"], 4), "\n")
  cat("  Test Kappa:",
      round(object$confusion_matrix_test_rf$overall["Kappa"], 4), "\n")
  cat("  Validation Accuracy:",
      round(object$confusion_matrix_validation_rf$overall["Accuracy"], 4), "\n")
  cat("  Validation Kappa:",
      round(object$confusion_matrix_validation_rf$overall["Kappa"], 4), "\n\n")

  cat("XGBOOST:\n")
  cat("  Test Accuracy:",
      round(object$confusion_matrix_test_xgb$overall["Accuracy"], 4), "\n")
  cat("  Test Kappa:",
      round(object$confusion_matrix_test_xgb$overall["Kappa"], 4), "\n")
  cat("  Validation Accuracy:",
      round(object$confusion_matrix_validation_xgb$overall["Accuracy"], 4), "\n")
  cat("  Validation Kappa:",
      round(object$confusion_matrix_validation_xgb$overall["Kappa"], 4), "\n\n")

  cat("COMPARISON:\n")
  print(object$model_comparison)

  cat("\n--- Confusion Matrix: Random Forest (Validation) ---\n")
  print(object$confusion_matrix_validation_rf$table)

  cat("\n--- Confusion Matrix: XGBoost (Validation) ---\n")
  print(object$confusion_matrix_validation_xgb$table)

  invisible(object)
}
