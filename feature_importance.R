feature_importance_shap <- function(classification_result,
                                    top_n      = 10,
                                    fill_color = "steelblue",
                                    verbose    = TRUE) {
  
 # -- Dependências necessárias ------------------------------------------------
 required_packages <- c("xgboost", "ggplot2", "shapviz")
 missing_packages   <- required_packages[
   !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
 ]

 if (length(missing_packages) > 0) {
   stop(
     "Os seguintes pacotes sao necessarios e nao estao instalados: ",
     paste(missing_packages, collapse = ", "),
     "\nInstale com: install.packages(c(",
     paste(sprintf("'%s'", missing_packages), collapse = ", "), "))"
   )
 }

 library(ggplot2)  # necessario para aes(), geom_bar(), theme_minimal(), etc.


  # Validate input
  if (is.null(classification_result$model)) {
    stop("classification_result$model not found. ",
         "This function requires ml_method = 'xgb'.")
  }
  
  if (classification_result$method != "xgb") {
    stop("This function requires ml_method = 'xgb'. ",
         "Current method: ", classification_result$method)
  }
  
  # -- Identify feature columns -----------------------------------------------
  motif_cols <- NULL
  
  for (field in c("motifs_used", "selected_motifs", "feature_names", "motif_names")) {
    if (!is.null(classification_result[[field]]) &&
        length(classification_result[[field]]) > 0) {
      motif_cols <- classification_result[[field]]
      if (verbose) cat(sprintf("  Feature columns from field '%s': %d\n",
                               field, length(motif_cols)))
      break
    }
  }
  
  # Fallback: derive from train_data excluding CLASS
  if (is.null(motif_cols)) {
    motif_cols <- colnames(classification_result$train_data)[
      colnames(classification_result$train_data) != "CLASS"
    ]
    if (verbose) cat(sprintf("  Feature columns derived from train_data: %d\n",
                             length(motif_cols)))
  }
  
  # -- 9a: XGBoost feature importance -----------------------------------------
importance_xgb <- xgboost::xgb.importance(
    feature_names = motif_cols,
    model         = classification_result$model
  )
  
  plot_importance_top <- function(imp, top_n = 10, fill_color = "steelblue") {
    df         <- imp[order(imp$Gain, decreasing = TRUE), ]
    df         <- head(df, top_n)
    df$Feature <- factor(df$Feature, levels = rev(df$Feature))
    
    ggplot(df, aes(x = Gain, y = Feature)) +
      geom_bar(stat = "identity", fill = fill_color) +
      geom_text(aes(label = round(Gain, 4)),
                hjust = -0.1, size = 3) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(x = "Gain (XGB)", y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y        = element_text(size = 10)
      )
  }
  
  plot_importance <- plot_importance_top(importance_xgb,
                                         top_n      = top_n,
                                         fill_color = fill_color)
  
  if (verbose) {
    print(plot_importance)
    cat("  XGBoost feature importance calculated\n")
  }
  
  # -- 9b: SHAP via shapviz ---------------------------------------------------
  missing_cols <- setdiff(motif_cols, colnames(classification_result$train_data))
  if (length(missing_cols) > 0) {
    warning(sprintf(
      "%d feature column(s) not found in train_data and will be removed: %s",
      length(missing_cols),
      paste(head(missing_cols, 5), collapse = ", ")
    ))
    motif_cols <- intersect(motif_cols, colnames(classification_result$train_data))
  }
  
  x_train_shapviz <- as.matrix(
    classification_result$train_data[, motif_cols, drop = FALSE]
  )
  
  shp <- shapviz(
    object = classification_result$model,
    X_pred = x_train_shapviz,
    X      = x_train_shapviz
  )
  
  plot_shap_beeswarm <- sv_importance(shp, kind = "beeswarm", max_display = top_n)
  
  if (verbose) {
    print(plot_shap_beeswarm)
    cat("  SHAP beeswarm plot generated (shapviz)\n\n")
  }
  
  list(
    importance_table   = importance_xgb,
    plot_importance    = plot_importance,
    shp                = shp,
    plot_shap_beeswarm = plot_shap_beeswarm,
    motif_cols         = motif_cols
  )
    }
