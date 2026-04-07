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
#'     \item train_data: Training data (motifs + CLASS)
#'     \item test_data: Test data (motifs + CLASS)
#'     \item predictions_train: Predicted classes for training set
#'     \item predictions_test: Predicted classes for test set
#'     \item actuals_train: Actual classes for training set
#'     \item actuals_test: Actual classes for test set
#'     \item confusion_matrix_train: Confusion matrix for training set
#'     \item confusion_matrix_test: Confusion matrix for test set
#'     \item model_metrics: Data frame with accuracy and kappa
#'     \item motifs_used: Character vector of motifs used
#'     \item params: Parameters used
#'   }
#'
#' @details
#' Trains either an XGBoost or Random Forest model with 5-fold cross-validation.
#'
#' XGBoost ("xgb"):
#' \itemize{
#'   \item Uses xgboost package directly (without caret) — avoids
#'     incompatibility between xgboost >= 2.x and caret 7.x
#'   \item CV finds the optimal nrounds (up to 100)
#'   \item Fixed parameters: cv_folds=5, nrounds=100, max_depth=6, eta=0.1,
#'     gamma=0, colsample_bytree=0.8, min_child_weight=1, subsample=0.8
#' }
#'
#' Random Forest ("rf"):
#' \itemize{
#'   \item Uses randomForest package via caret
#'   \item Fixed parameters: cv_folds=5, ntree=100
#' }
#'
#' ML convention used:
#' \itemize{
#'   \item \strong{Training set}: Fits model parameters (train_dataset)
#'   \item \strong{CV}: Evaluates model during training (5 folds)
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
#' result_xgb <- train_model_xgboost(
#'   result_train_test      = result_train_test,
#'   result_selected_motifs = result_selected_motifs,
#'   ml_method              = "xgb",
#'   verbose                = TRUE
#' )
#'
#' # Random Forest
#' result_rf <- train_model_xgboost(
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

  # Check required packages
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

  # Fixed parameters — shared
  cv_folds <- 5L

  # Fixed parameters — XGBoost
  nrounds          <- 100L
  max_depth        <- 6L
  eta              <- 0.1
  gamma            <- 0
  colsample_bytree <- 0.8
  min_child_weight <- 1
  subsample        <- 0.8

  # Fixed parameters — Random Forest
  ntree <- 100L

  #  Extract datasets
  train_dataset <- result_train_test$train_dataset
  test_dataset  <- result_train_test$test_dataset

  if (is.null(train_dataset)) stop("result_train_test$train_dataset is NULL.")
  if (is.null(test_dataset))  stop("result_train_test$test_dataset is NULL.")

  # Extract selected motifs
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

  # Verificar motifs ausentes no test_dataset
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
  # 3. Subsetar colunas dos motifs + CLASS
  # ---------------------------------------------------------------------------
  traindata <- train_dataset[, c(motifs_ok, "CLASS"), drop = FALSE]
  testdata  <- test_dataset[,  c(motifs_ok, "CLASS"), drop = FALSE]

  # ---------------------------------------------------------------------------
  # 4. Converter CLASS para factor
  # ---------------------------------------------------------------------------
  traindata$CLASS <- factor(traindata$CLASS)
  testdata$CLASS  <- factor(testdata$CLASS, levels = levels(traindata$CLASS))

  # ---------------------------------------------------------------------------
  # 5. Diagnóstico
  # ---------------------------------------------------------------------------
  if (verbose) {
    cat("Data summary\n")
    cat("Training :", nrow(traindata), "sequences\n")
    cat("Test     :", nrow(testdata),  "sequences\n")
    cat("Motifs   :", length(motifs_ok), "| Classes:", nlevels(traindata$CLASS), "\n\n")

    cat("Distribuição de classes em traindata:\n")
    print(table(traindata$CLASS))
    cat("Sequências por fold (estimado):",
        floor(min(table(traindata$CLASS)) / cv_folds), "\n\n")
  }

  start_time <- Sys.time()

  # ---------------------------------------------------------------------------
  # 6a. XGBOOST
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

    best_nrounds <- which.min(cv_result$evaluation_log$test_merror_mean)
    best_error   <- min(cv_result$evaluation_log$test_merror_mean)

    if (verbose) {
      cat("Melhor nrounds :", best_nrounds, "\n")
      cat("Acurácia CV    :", round(1 - best_error, 4), "\n\n")
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

    if (ml_method == "rf"){

    cv_result    <- NULL
    best_nrounds <- NULL

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

    pred_train_classes <- predict(model, newdata = traindata)
    pred_test_classes  <- predict(model, newdata = testdata)

    model_params <- params
  }

  time_total <- round(as.numeric(difftime(Sys.time(), start_time, units = "secs")), 2)

  if (verbose) cat("Completed in", time_total, "seconds\n\n")

  # ---------------------------------------------------------------------------
  # 7. Métricas
  # ---------------------------------------------------------------------------
  cm_train <- caret::confusionMatrix(pred_train_classes, traindata$CLASS)
  cm_test  <- caret::confusionMatrix(pred_test_classes,  testdata$CLASS)

  model_metrics <- data.frame(
    Set      = c("Train", "Test"),
    Accuracy = round(c(cm_train$overall["Accuracy"],
                       cm_test$overall["Accuracy"]), 4),
    Kappa    = round(c(cm_train$overall["Kappa"],
                       cm_test$overall["Kappa"]),    4),
    stringsAsFactors = FALSE
  )

  if (verbose) {
    cat("Model metrics\n")
    print(model_metrics)
    cat("\n")

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
  # 8. Salvar resultados
  # ---------------------------------------------------------------------------
#  file_model   <- paste0("model_", ml_method, ".RData")
#  file_metrics <- paste0("model_metrics_", ml_method, ".RData")

#  save(model,         file = file_model)
#  save(model_metrics, file = file_metrics)

#  if (verbose) {
#    cat("Files saved:", file_model, ",", file_metrics, "\n\n")
#  }

  # ---------------------------------------------------------------------------
  # 9. Retornar resultado
  # ---------------------------------------------------------------------------
   structure(
    list(
      method                 = ml_method,
      model                  = model,
      cv_result              = cv_result,
      best_nrounds           = best_nrounds,
      train_data             = traindata,
      test_data              = testdata,
      predictions_train      = pred_train_classes,
      predictions_test       = pred_test_classes,
      actuals_train          = traindata$CLASS,
      actuals_test           = testdata$CLASS,
      confusion_matrix_train = cm_train,
      confusion_matrix_test  = cm_test,
      model_metrics          = model_metrics,
      motifs_used            = motifs_ok,
      cv_folds               = cv_folds,
      time_seconds           = time_total,
      params                 = model_params
    ),
    class = c("train_model_xgboost", "list")
  )
}


