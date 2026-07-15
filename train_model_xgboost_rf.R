#' Train XGBoost or Random Forest Model with Cross-Validation
#'
#' @param result_train_test List from select_train_test1() containing 
#'   train_dataset and test_dataset.
#' @param result_selected_motifs List from select_motifs() with character
#'   vectors of motifs per class.
#' @param ml_method Character. Machine learning method to use: "xgb" for
#'   XGBoost or "rf" for Random Forest. Default: "xgb".
#' @param verbose Logical. If TRUE, prints detailed information. Default: TRUE.
#'
#' @return Object of class "train_model_xgboost" with:
#'   \itemize{
#'     \item method: Method used ("xgb" or "rf")
#'     \item model: Trained model object
#'     \item cv_result: Cross-validation log (xgb only; NULL for rf)
#'     \item best_nrounds: Optimal boosting rounds from CV (xgb only; NULL for rf)
#'     \item cv_accuracy: Mean CV accuracy across folds (xgb only; NULL for rf)
#'     \item cv_accuracy_sd: SD of CV accuracy across folds (xgb only; NULL for rf)
#'     \item train_data: Training data (motifs + CLASS)
#'     \item test_data: Test data (motifs + CLASS)
#'     \item predictions_train: Predicted classes for training set
#'     \item predictions_test: Predicted classes for test set
#'     \item actuals_train: Actual classes for training set
#'     \item actuals_test: Actual classes for test set
#'     \item confusion_matrix_train: Confusion matrix for training set
#'     \item confusion_matrix_test: Confusion matrix for test set
#'     \item model_metrics: Data frame with accuracy and kappa (Train, CV, Test)
#'     \item motifs_used: Character vector of motifs used
#'     \item params: Parameters used
#'   }
#'
#' @details
#' Trains either an XGBoost or Random Forest model with 10-fold cross-validation.
#'
#' XGBoost ("xgb"):
#' \itemize{
#'   \item Uses xgboost package directly (without caret) — avoids
#'     incompatibility between xgboost >= 2.x and caret 7.x
#'   \item CV finds the optimal nrounds (up to 100)
#'   \item Fixed parameters: cv_folds=10, nrounds=100, max_depth=6, eta=0.1,
#'     gamma=0, colsample_bytree=0.8, min_child_weight=1, subsample=0.8
#'   \item CV accuracy (mean and SD across folds) is reported in model_metrics
#'     and accessible via $cv_accuracy and $cv_accuracy_sd
#' }
#'
#' Random Forest ("rf"):
#' \itemize{
#'   \item Uses randomForest package via caret
#'   \item Fixed parameters: cv_folds=10, ntree=100
#' }
#'
#' ML convention used:
#' \itemize{
#'   \item \strong{Training set}: Fits model parameters (train_dataset)
#'   \item \strong{CV (validation)}: Evaluates model during training (10 folds)
#'   \item \strong{Test set}: Final independent evaluation (test_dataset)
#' }
#'
#' @note
#' \itemize{
#'   \item Requires xgboost package for ml_method = "xgb".
#'   \item Requires randomForest package for ml_method = "rf".
#'   \item Requires caret package for confusionMatrix and RF training.
#'   \item S3 methods available. See \code{methods.R} for details.
#' }
#'
#' @examples
#' \dontrun{
#' # XGBoost (default)
#' result_xgb <- train_model_xgboost_rf(
#'   result_train_test      = result_train_test,
#'   result_selected_motifs = result_selected_motifs,
#'   ml_method              = "xgb",
#'   verbose                = TRUE
#' )
#' # Access CV accuracy directly:
#' result_xgb$cv_accuracy
#' result_xgb$cv_accuracy_sd
#'
#' # Full metrics table (Train / CV / Test):
#' result_xgb$model_metrics
#'
#' # Random Forest
#' result_rf <- train_model_xgboost_rf(
#'   result_train_test      = result_train_test,
#'   result_selected_motifs = result_selected_motifs,
#'   ml_method              = "rf",
#'   verbose                = TRUE
#' )
#' }
#'
#' @importFrom xgboost xgb.DMatrix xgb.cv xgb.train
#' @importFrom randomForest randomForest
#' @importFrom caret trainControl train confusionMatrix
#'
#' @export
train_model_xgboost_rf <- function(result_train_test,
                                   result_selected_motifs,
                                   ml_method = c("xgb", "rf"),
                                   verbose   = TRUE) {
  
  ml_method <- match.arg(ml_method)
  
  # ---------------------------------------------------------------------------
  # 1. Check required packages
  # ---------------------------------------------------------------------------
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' required.\n",
         "Install with: install.packages('caret')")
  }
  
  if (ml_method == "xgb" && !requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' required for ml_method = 'xgb'.\n",
         "Install with: install.packages('xgboost')")
  }
  
  if (ml_method == "rf" && !requireNamespace("randomForest", quietly = TRUE)) {
    stop("Package 'randomForest' required for ml_method = 'rf'.\n",
         "Install with: install.packages('randomForest')")
  }
  
  set.seed(123)
  
  # ---------------------------------------------------------------------------
  # 2. Fixed parameters
  # ---------------------------------------------------------------------------
  
  # Shared
  cv_folds <- 10L
  
  # XGBoost
  nrounds          <- 100L
  max_depth        <- 6L
  eta              <- 0.1
  gamma            <- 0
  colsample_bytree <- 0.8
  min_child_weight <- 1
  subsample        <- 0.8
  
  # Random Forest
  ntree <- 100L
  
  # ---------------------------------------------------------------------------
  # 3. Extract datasets
  # ---------------------------------------------------------------------------
  train_dataset <- result_train_test$train_dataset
  test_dataset  <- result_train_test$test_dataset
  
  if (is.null(train_dataset)) stop("result_train_test$train_dataset is NULL.")
  if (is.null(test_dataset))  stop("result_train_test$test_dataset is NULL.")
  
  # ---------------------------------------------------------------------------
  # 4. Extract selected motifs
  # ---------------------------------------------------------------------------
  motifs_all <- unique(unlist(result_selected_motifs))
  motifs_all <- setdiff(motifs_all, "CLASS")
  motifs_ok  <- motifs_all[motifs_all %in% colnames(train_dataset)]
  
  if (length(motifs_ok) == 0) {
    stop("Nenhum motif selecionado encontrado nas colunas de train_dataset.")
  }
  
  if (verbose) {
    cat("Method              :", toupper(ml_method), "\n")
    cat("Motifs selecionados :", length(motifs_all), "\n")
    cat("Motifs presentes    :", length(motifs_ok),  "\n\n")
  }
  
  # ---------------------------------------------------------------------------
  # 5. Handle motifs missing in test_dataset
  # ---------------------------------------------------------------------------
  motifs_missing_test <- setdiff(motifs_ok, colnames(test_dataset))
  if (length(motifs_missing_test) > 0) {
    warning(sprintf(
      "%d motif(s) ausentes em test_dataset — preenchidos com 0.",
      length(motifs_missing_test)
    ))
    for (col in motifs_missing_test) {
      test_dataset[[col]] <- 0L
    }
  }
  
  # ---------------------------------------------------------------------------
  # 6. Subset columns: motifs + CLASS
  # ---------------------------------------------------------------------------
  traindata <- train_dataset[, c(motifs_ok, "CLASS"), drop = FALSE]
  testdata  <- test_dataset[,  c(motifs_ok, "CLASS"), drop = FALSE]
  
  # ---------------------------------------------------------------------------
  # 7. Convert CLASS to factor
  # ---------------------------------------------------------------------------
  traindata$CLASS <- factor(traindata$CLASS)
  testdata$CLASS  <- factor(testdata$CLASS, levels = levels(traindata$CLASS))
  
  # ---------------------------------------------------------------------------
  # 8. Diagnostics
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("Data summary\n")
    cat("Training :", nrow(traindata), "sequences\n")
    cat("Test     :", nrow(testdata),  "sequences\n")
    cat("Motifs   :", length(motifs_ok), "| Classes:", nlevels(traindata$CLASS), "\n\n")
    
    cat("Distribuicao de classes em traindata:\n")
    print(table(traindata$CLASS))
    cat("Sequencias por fold (estimado):",
        floor(min(table(traindata$CLASS)) / cv_folds), "\n\n")
  }
  
  start_time <- Sys.time()
  
  # Placeholders — overwritten per method
  cv_result    <- NULL
  best_nrounds <- NULL
  cv_accuracy  <- NULL
  cv_accuracy_sd <- NULL
  
  # ---------------------------------------------------------------------------
  # 9a. XGBOOST
  # ---------------------------------------------------------------------------
  if (ml_method == "xgb") {
    
    x_train <- as.matrix(traindata[, motifs_ok, drop = FALSE])
    y_train <- as.integer(traindata$CLASS) - 1L
    x_test  <- as.matrix(testdata[, motifs_ok, drop = FALSE])
    y_test  <- as.integer(testdata$CLASS) - 1L
    
    dtrain <- xgboost::xgb.DMatrix(data = x_train, label = y_train)
    dtest  <- xgboost::xgb.DMatrix(data = x_test,  label = y_test)
    
    params <- list(
      booster          = "gbtree",
      objective        = "multi:softmax",
      num_class        = nlevels(traindata$CLASS),
      max_depth        = max_depth,
      eta              = eta,
      gamma            = gamma,
      colsample_bytree = colsample_bytree,
      min_child_weight = min_child_weight,
      subsample        = subsample
    )
    
    if (verbose) cat("Running", cv_folds, "folds cross-validation (XGBoost)...\n")
    
    cv_result <- xgboost::xgb.cv(
      params  = params,
      data    = dtrain,
      nrounds = nrounds,
      nfold   = cv_folds,
      metrics = "merror",
      verbose = FALSE
    )
    
    best_nrounds   <- which.min(cv_result$evaluation_log$test_merror_mean)
    best_error     <- cv_result$evaluation_log$test_merror_mean[best_nrounds]
    best_error_sd  <- cv_result$evaluation_log$test_merror_std[best_nrounds]
    
    # CV accuracy (mean and SD across folds at best_nrounds)
    cv_accuracy    <- round(1 - best_error,    4)
    cv_accuracy_sd <- round(best_error_sd,     4)
    
    if (verbose) {
      cat("Melhor nrounds :", best_nrounds,   "\n")
      cat("Acuracia CV    :", cv_accuracy,     "(+/-", cv_accuracy_sd, "SD)\n\n")
      cat("Training final XGBoost model...\n")
    }
    
    model <- xgboost::xgb.train(
      params  = params,
      data    = dtrain,
      nrounds = best_nrounds,
      verbose = 0
    )
    
    pred_train_raw <- predict(model, dtrain)
    pred_test_raw  <- predict(model, dtest)
    
    pred_train_classes <- factor(
      levels(traindata$CLASS)[pred_train_raw + 1L],
      levels = levels(traindata$CLASS)
    )
    pred_test_classes <- factor(
      levels(traindata$CLASS)[pred_test_raw + 1L],
      levels = levels(traindata$CLASS)
    )
    
    model_params <- params
  }
  
  # ---------------------------------------------------------------------------
  # 9b. RANDOM FOREST
  # ---------------------------------------------------------------------------
  if (ml_method == "rf") {
    
    params <- list(
      ntree    = ntree,
      cv_folds = cv_folds
    )
    
    control <- caret::trainControl(
      method          = "cv",
      number          = cv_folds,
      savePredictions = "final",
      classProbs      = FALSE,
      allowParallel   = FALSE,
      verboseIter     = FALSE,
      returnData      = FALSE
    )
    
    if (verbose) cat("Training Random Forest (", cv_folds, "folds CV)...\n")
    
    model_rf_tmp <- suppressWarnings(
      suppressMessages(
        caret::train(
          CLASS ~ .,
          data       = traindata,
          method     = "rf",
          trControl  = control,
          ntree      = ntree,
          importance = TRUE,
          verbose    = TRUE
        )
      )
    )
    
    model <- model_rf_tmp
    rm(model_rf_tmp)
    
    # Extract CV accuracy from caret resample (one row per fold)
    # model$resample contains Accuracy and Kappa for each fold
    cv_accuracy    <- round(mean(model$resample$Accuracy), 4)
    cv_accuracy_sd <- round(sd(model$resample$Accuracy),   4)
    cv_kappa_mean  <- round(mean(model$resample$Kappa),    4)
    
    if (verbose) {
      cat("CV Accuracy per fold:\n")
      print(round(model$resample$Accuracy, 4))
      cat("Acuracia CV :", cv_accuracy, "(+/-", cv_accuracy_sd, "SD)\n\n")
    }
    
    pred_train_classes <- predict(model, newdata = traindata)
    pred_test_classes  <- predict(model, newdata = testdata)
    
    model_params <- params
  }
  
  time_total <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)
  
  if (verbose) cat("Completed in", time_total, "seconds\n\n")
  
  # ---------------------------------------------------------------------------
  # 10. Metrics
  # ---------------------------------------------------------------------------
  cm_train <- caret::confusionMatrix(pred_train_classes, traindata$CLASS)
  cm_test  <- caret::confusionMatrix(pred_test_classes,  testdata$CLASS)
  
  # Build model_metrics — includes CV row for XGBoost
  if (ml_method == "xgb") {
    model_metrics <- data.frame(
      Set           = c("Train", "CV (validation)", "Test"),
      Accuracy      = c(
        round(cm_train$overall["Accuracy"], 4),
        cv_accuracy,
        round(cm_test$overall["Accuracy"],  4)
      ),
      Accuracy_SD   = c(NA, cv_accuracy_sd, NA),
      Kappa         = c(
        round(cm_train$overall["Kappa"], 4),
        NA,
        round(cm_test$overall["Kappa"],  4)
      ),
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  } else {
    # RF: CV accuracy comes from model$resample extracted above
    model_metrics <- data.frame(
      Set           = c("Train", "CV (validation)", "Test"),
      Accuracy      = c(
        round(cm_train$overall["Accuracy"], 4),
        cv_accuracy,
        round(cm_test$overall["Accuracy"],  4)
      ),
      Accuracy_SD   = c(NA, cv_accuracy_sd, NA),
      Kappa         = c(
        round(cm_train$overall["Kappa"], 4),
        cv_kappa_mean,
        round(cm_test$overall["Kappa"],  4)
      ),
      stringsAsFactors = FALSE,
      row.names        = NULL
    )
  }
  
  if (verbose) {
    cat("Model metrics\n")
    print(model_metrics)
    cat("\n")
    
    if (ml_method == "xgb") {
      cat("CV log at best_nrounds (", best_nrounds, "):\n")
      print(cv_result$evaluation_log[best_nrounds, ])
      cat("\n")
    }
    
    if (ml_method == "rf") {
      cat("CV Accuracy per fold (RF):\n")
      fold_df <- data.frame(
        Fold     = model$resample$Resample,
        Accuracy = round(model$resample$Accuracy, 4),
        Kappa    = round(model$resample$Kappa,    4)
      )
      print(fold_df)
      cat("\n")
    }
    
    cat("Confusion Matrix — Train set\n\n")
    print(cm_train$table)
    cat("\nOverall Statistics:\n")
    cat("  Accuracy:", round(cm_train$overall["Accuracy"], 4), "\n")
    cat("  Kappa   :", round(cm_train$overall["Kappa"],    4), "\n")
    cat("  95% CI: (",
        round(cm_train$overall["AccuracyLower"], 4), ", ",
        round(cm_train$overall["AccuracyUpper"], 4), ")\n\n")
    
    cat("Confusion Matrix — Test set\n\n")
    print(cm_test$table)
    cat("\nOverall Statistics:\n")
    cat("  Accuracy:", round(cm_test$overall["Accuracy"], 4), "\n")
    cat("  Kappa   :", round(cm_test$overall["Kappa"],    4), "\n")
    cat("  95% CI: (",
        round(cm_test$overall["AccuracyLower"], 4), ", ",
        round(cm_test$overall["AccuracyUpper"], 4), ")\n\n")
  }
  
  # ---------------------------------------------------------------------------
  # 11. Return result
  # ---------------------------------------------------------------------------
  structure(
    list(
      method                 = ml_method,
      model                  = model,
      cv_result              = cv_result,
      best_nrounds           = best_nrounds,
      cv_accuracy            = cv_accuracy,        # mean CV accuracy (xgb & rf)
      cv_accuracy_sd         = cv_accuracy_sd,     # SD of CV accuracy (xgb & rf)
      cv_kappa_mean          = if (ml_method == "rf") cv_kappa_mean else NULL,  # rf only
      train_data             = traindata,
      test_data              = testdata,
      predictions_train      = pred_train_classes,
      predictions_test       = pred_test_classes,
      actuals_train          = traindata$CLASS,
      actuals_test           = testdata$CLASS,
      confusion_matrix_train = cm_train,
      confusion_matrix_test  = cm_test,
      model_metrics          = model_metrics,      # Includes CV row for xgb and rf
      motifs_used            = motifs_ok,
      cv_folds               = cv_folds,
      time_seconds           = time_total,
      params                 = model_params
    ),
    class = c("train_model_xgboost", "list")
  )
}
