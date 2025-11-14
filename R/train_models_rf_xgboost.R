#' Train Random Forest and XGBoost Models
#'
#' @param dataset_traintest List containing classification_dataset and validation_dataset.
#' @param selected_motifs List of character vectors containing selected motifs per class.
#' @param prop_train Numeric. Proportion for training (0-1). Default: 0.7.
#' @param cv_folds Integer. Number of cross-validation folds. Default: 5.
#' @param seed Integer. Random seed for reproducibility. Default: 123.
#'
#' @return Object of class train_models_rf_xgboost with models and predictions.
#'
#' @export
train_models_rf_xgboost <- function(dataset_traintest,
                                    selected_motifs,
                                    prop_train = 0.7,
                                    cv_folds = 5,
                                    seed = 123) {

  # === INPUT VALIDATION ===

  if (!is.list(dataset_traintest) || is.data.frame(dataset_traintest)) {
    stop("'dataset_traintest' must be a list (not a data.frame)")
  }

  if (!("classification_dataset" %in% names(dataset_traintest))) {
    stop("'dataset_traintest' must contain 'classification_dataset'")
  }

  if (!("validation_dataset" %in% names(dataset_traintest))) {
    stop("'dataset_traintest' must contain 'validation_dataset'")
  }

  if (!is.list(selected_motifs) || length(selected_motifs) == 0) {
    stop("'selected_motifs' must be a non-empty list of motifs per class")
  }

  if (!is.numeric(prop_train) || prop_train <= 0 || prop_train >= 1) {
    stop("'prop_train' must be a number between 0 and 1")
  }

  if (!is.numeric(cv_folds) || cv_folds <= 0 || cv_folds != as.integer(cv_folds)) {
    stop("'cv_folds' must be a positive integer")
  }

  if (!is.numeric(seed) || seed != as.integer(seed)) {
    stop("'seed' must be an integer")
  }

  # === PACKAGE REQUIREMENTS ===

  required_packages <- c("caret", "randomForest", "xgboost")

  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Package '%s' is required but not installed. Install with: install.packages('%s')",
                   pkg, pkg))
    }
  }

  set.seed(seed)

  cat("=== DATA PREPARATION ===\n\n")

  classification_dataset <- dataset_traintest$classification_dataset
  validation_dataset <- dataset_traintest$validation_dataset

  if ("CLASS" %in% colnames(classification_dataset) &&
      !("class" %in% colnames(classification_dataset))) {
    cat("Renaming 'CLASS' to 'class' in classification_dataset...\n")
    colnames(classification_dataset)[colnames(classification_dataset) == "CLASS"] <- "class"
  }

  if ("CLASS" %in% colnames(validation_dataset) &&
      !("class" %in% colnames(validation_dataset))) {
    cat("Renaming 'CLASS' to 'class' in validation_dataset...\n")
    colnames(validation_dataset)[colnames(validation_dataset) == "CLASS"] <- "class"
  }

  motifs_all <- unique(unlist(selected_motifs))

  cat("Total unique motifs selected:", length(motifs_all), "\n")
  cat("Classes in motif selector:", paste(names(selected_motifs), collapse = ", "), "\n\n")

  required_columns <- c(motifs_all, "class")
  missing_class <- setdiff(required_columns, colnames(classification_dataset))
  missing_valid <- setdiff(required_columns, colnames(validation_dataset))

  if (length(missing_class) > 0) {
    stop(sprintf("Missing columns in classification_dataset: %s",
                 paste(missing_class, collapse = ", ")))
  }

  if (length(missing_valid) > 0) {
    stop(sprintf("Missing columns in validation_dataset: %s",
                 paste(missing_valid, collapse = ", ")))
  }

  train_partition <- caret::createDataPartition(
    classification_dataset$class,
    p = prop_train,
    list = FALSE
  )

  traindata <- classification_dataset[train_partition, ]
  testdata <- classification_dataset[-train_partition, ]

  train_selected <- traindata[, required_columns, drop = FALSE]
  test_selected <- testdata[, required_columns, drop = FALSE]
  validation_selected <- validation_dataset[, required_columns, drop = FALSE]

  train_selected$class <- factor(train_selected$class)
  test_selected$class <- factor(test_selected$class, levels = levels(train_selected$class))
  validation_selected$class <- factor(validation_selected$class,
                                      levels = levels(train_selected$class))

  cat("Training set observations:", nrow(train_selected), "\n")
  cat("Test set observations:", nrow(test_selected), "\n")
  cat("Validation set observations:", nrow(validation_selected), "\n")
  cat("Number of selected motifs:", length(motifs_all), "\n\n")

  cat("Class distribution in training set:\n")
  print(table(train_selected$class))

  cat("\nClass distribution in test set:\n")
  print(table(test_selected$class))

  cat("\nClass distribution in validation set:\n")
  print(table(validation_selected$class))

  control <- caret::trainControl(
    method = "cv",
    number = cv_folds,
    savePredictions = "final",
    classProbs = TRUE,
    allowParallel = TRUE
  )

  # ==================== RANDOM FOREST ====================
  cat("\n=== RANDOM FOREST MODEL TRAINING ===\n")
  cat("Training Random Forest model with", cv_folds, "cross-validation folds...\n\n")

  start_time_rf <- Sys.time()

  model_rf <- suppressWarnings(
    suppressMessages(
      caret::train(
        class ~ .,
        data = train_selected,
        method = "rf",
        trControl = control,
        ntree = 500,
        importance = TRUE,
        verbose = FALSE
      )
    )
  )

  time_rf <- as.numeric(difftime(Sys.time(), start_time_rf, units = "secs"))

  cat("Random Forest model trained successfully! Time:", round(time_rf, 2), "seconds\n")
  cat("Top 10 important variables:\n")
  var_imp_rf <- caret::varImp(model_rf)
  print(head(var_imp_rf, 10))

  cat("\n--- Random Forest - Test Set ---\n")
  pred_test_rf <- predict(model_rf, newdata = test_selected)
  actuals_test <- test_selected$class
  confusion_matrix_test_rf <- caret::confusionMatrix(pred_test_rf, actuals_test)

  cat("Accuracy on test set:", round(confusion_matrix_test_rf$overall["Accuracy"], 4), "\n")
  cat("Kappa on test set:", round(confusion_matrix_test_rf$overall["Kappa"], 4), "\n")

  cat("\n--- Random Forest - Validation Set ---\n")
  pred_validation_rf <- predict(model_rf, newdata = validation_selected)
  actuals_validation <- validation_selected$class
  confusion_matrix_validation_rf <- caret::confusionMatrix(pred_validation_rf, actuals_validation)

  cat("Accuracy on validation set:", round(confusion_matrix_validation_rf$overall["Accuracy"], 4), "\n")
  cat("Kappa on validation set:", round(confusion_matrix_validation_rf$overall["Kappa"], 4), "\n")

  # ==================== XGBOOST ====================
  cat("\n=== XGBOOST MODEL TRAINING ===\n")
  cat("Training XGBoost model with", cv_folds, "cross-validation folds...\n\n")

  start_time_xgb <- Sys.time()

  invisible(
    capture.output(
      {
        model_xgb <- suppressWarnings(
          suppressMessages(
            caret::train(
              class ~ .,
              data = train_selected,
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

  cat("XGBoost model trained successfully! Time:", round(time_xgb, 2), "seconds\n")
  cat("Top 10 important variables:\n")
  var_imp_xgb <- caret::varImp(model_xgb)
  print(head(var_imp_xgb, 10))

  cat("\n--- XGBoost - Test Set ---\n")
  pred_test_xgb <- predict(model_xgb, newdata = test_selected)
  confusion_matrix_test_xgb <- caret::confusionMatrix(pred_test_xgb, actuals_test)

  cat("Accuracy on test set:", round(confusion_matrix_test_xgb$overall["Accuracy"], 4), "\n")
  cat("Kappa on test set:", round(confusion_matrix_test_xgb$overall["Kappa"], 4), "\n")

  cat("\n--- XGBoost - Validation Set ---\n")
  pred_validation_xgb <- predict(model_xgb, newdata = validation_selected)
  confusion_matrix_validation_xgb <- caret::confusionMatrix(pred_validation_xgb, actuals_validation)

  cat("Accuracy on validation set:", round(confusion_matrix_validation_xgb$overall["Accuracy"], 4), "\n")
  cat("Kappa on validation set:", round(confusion_matrix_validation_xgb$overall["Kappa"], 4), "\n")

  # ==================== MODEL COMPARISON ====================
  cat("\n=== MODEL COMPARISON ===\n\n")

  acc_test_rf <- confusion_matrix_test_rf$overall["Accuracy"]
  acc_test_xgb <- confusion_matrix_test_xgb$overall["Accuracy"]
  kappa_test_rf <- confusion_matrix_test_rf$overall["Kappa"]
  kappa_test_xgb <- confusion_matrix_test_xgb$overall["Kappa"]

  acc_validation_rf <- confusion_matrix_validation_rf$overall["Accuracy"]
  acc_validation_xgb <- confusion_matrix_validation_xgb$overall["Accuracy"]
  kappa_validation_rf <- confusion_matrix_validation_rf$overall["Kappa"]
  kappa_validation_xgb <- confusion_matrix_validation_xgb$overall["Kappa"]

  model_comparison_table <- data.frame(
    Model = c("Random Forest", "XGBoost"),
    Accuracy_Test = c(round(acc_test_rf, 4), round(acc_test_xgb, 4)),
    Kappa_Test = c(round(kappa_test_rf, 4), round(kappa_test_xgb, 4)),
    Accuracy_Validation = c(round(acc_validation_rf, 4), round(acc_validation_xgb, 4)),
    Kappa_Validation = c(round(kappa_validation_rf, 4), round(kappa_validation_xgb, 4)),
    Training_Time_s = c(round(time_rf, 2), round(time_xgb, 2)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  cat("Performance Metrics Comparison:\n")
  print(model_comparison_table)

  best_test <- if (acc_test_rf > acc_test_xgb) "Random Forest" else "XGBoost"
  best_validation <- if (acc_validation_rf > acc_validation_xgb) "Random Forest" else "XGBoost"

  cat("\n Best model on test set:", best_test,
      "(Accuracy:", round(max(acc_test_rf, acc_test_xgb), 4), ")\n")
  cat(" Best model on validation set:", best_validation,
      "(Accuracy:", round(max(acc_validation_rf, acc_validation_xgb), 4), ")\n")

  cat("\n--- Per-Class Metrics (Random Forest - Test) ---\n")
  print(confusion_matrix_test_rf$byClass)

  cat("\n--- Per-Class Metrics (XGBoost - Test) ---\n")
  print(confusion_matrix_test_xgb$byClass)

  cat("\n=== SAVING RESULTS ===\n\n")

  save(model_rf, file = "model_rf.RData")
  save(model_xgb, file = "model_xgboost.RData")
  save(model_comparison_table, file = "model_comparison.RData")

  cat(" Models saved:\n")
  cat("  - model_rf.RData\n")
  cat("  - model_xgboost.RData\n")
  cat("  - model_comparison.RData\n\n")

  result <- list(
    model_rf = model_rf,
    model_xgb = model_xgb,
    train_data = train_selected,
    test_data = test_selected,
    validation_data = validation_selected,
    predictions_test_rf = pred_test_rf,
    predictions_validation_rf = pred_validation_rf,
    predictions_test_xgb = pred_test_xgb,
    predictions_validation_xgb = pred_validation_xgb,
    actuals_test = actuals_test,
    actuals_validation = actuals_validation,
    confusion_matrix_test_rf = confusion_matrix_test_rf,
    confusion_matrix_validation_rf = confusion_matrix_validation_rf,
    confusion_matrix_test_xgb = confusion_matrix_test_xgb,
    confusion_matrix_validation_xgb = confusion_matrix_validation_xgb,
    model_comparison = model_comparison_table,
    best_model_test = best_test,
    best_model_validation = best_validation,
    motifs_used = motifs_all,
    motifs_per_class = selected_motifs,
    cv_folds = cv_folds,
    prop_train = prop_train,
    seed = seed,
    time_rf = time_rf,
    time_xgb = time_xgb
  )

  class(result) <- c("train_models_rf_xgboost", "list")

  return(result)
}


#' Print method for train_models_rf_xgboost objects
#'
#' @param x Object of class train_models_rf_xgboost.
#' @param ... Additional arguments (unused).
#'
#' @export
print.train_models_rf_xgboost <- function(x, ...) {
  cat("=== Train Models RF XGBoost Object ===\n\n")
  cat("Models trained with", x$cv_folds, "cross-validation folds\n")
  cat("Training proportion:", x$prop_train, "\n")
  cat("Number of motifs used:", length(x$motifs_used), "\n")
  cat("Number of classes:", nlevels(x$actuals_test), "\n\n")
  cat("Performance Summary:\n")
  print(x$model_comparison)
  cat("\nBest model on test set:", x$best_model_test, "\n")
  cat("Best model on validation set:", x$best_model_validation, "\n")
}


#' Summary method for train_models_rf_xgboost objects
#'
#' @param object Object of class train_models_rf_xgboost.
#' @param ... Additional arguments (unused).
#'
#' @export
summary.train_models_rf_xgboost <- function(object, ...) {
  cat("=== Model Summary ===\n\n")
  cat("Random Forest - Test Set Accuracy:",
      round(object$confusion_matrix_test_rf$overall["Accuracy"], 4), "\n")
  cat("XGBoost - Test Set Accuracy:",
      round(object$confusion_matrix_test_xgb$overall["Accuracy"], 4), "\n\n")
  cat("Random Forest - Validation Set Accuracy:",
      round(object$confusion_matrix_validation_rf$overall["Accuracy"], 4), "\n")
  cat("XGBoost - Validation Set Accuracy:",
      round(object$confusion_matrix_validation_xgb$overall["Accuracy"], 4), "\n\n")
  cat("Training time - RF:", round(object$time_rf, 2), "seconds\n")
  cat("Training time - XGBoost:", round(object$time_xgb, 2), "seconds\n")
}
