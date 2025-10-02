treinar_modelos_rf_xgboost <- function(classification_dataset, 
                                       validation_dataset, 
                                       motifs_selecionados, 
                                       prop_treino = 0.7,
                                       cv_folds = 5,
                                       seed = 123) {
  
  # Carregar bibliotecas necessárias
  library(caret)
  library(randomForest)
  library(xgboost)
  
  # Definir seed para reprodutibilidade
  set.seed(seed)
  
  cat("=== PREPARAÇÃO DOS DADOS ===\n")
  
  # Verificar se as colunas necessárias existem
  colunas_necessarias <- c(motifs_selecionados, "subtype")
  colunas_faltando_class <- setdiff(colunas_necessarias, colnames(classification_dataset))
  colunas_faltando_valid <- setdiff(colunas_necessarias, colnames(validation_dataset))
  
  if (length(colunas_faltando_class) > 0) {
    stop(paste("Colunas não encontradas no classification_dataset:", 
               paste(colunas_faltando_class, collapse = ", ")))
  }
  
  if (length(colunas_faltando_valid) > 0) {
    stop(paste("Colunas não encontradas no validation_dataset:", 
               paste(colunas_faltando_valid, collapse = ", ")))
  }
  
  # Particionar classification_dataset em treino e teste (70/30)
  train_partition <- createDataPartition(classification_dataset$subtype, 
                                         p = prop_treino, 
                                         list = FALSE)
  
  traindata <- classification_dataset[train_partition, ]
  testdata <- classification_dataset[-train_partition, ]
  
  # Selecionar apenas as colunas especificadas
  train_selected <- traindata[, colunas_necessarias]
  test_selected <- testdata[, colunas_necessarias]
  validation_selected <- validation_dataset[, colunas_necessarias]
  
  # Informações sobre os datasets
  cat("Número de observações no conjunto de treino:", nrow(train_selected), "\n")
  cat("Número de observações no conjunto de teste:", nrow(test_selected), "\n")
  cat("Número de observações no conjunto de validação:", nrow(validation_selected), "\n")
  cat("Número de motifs selecionados:", length(motifs_selecionados), "\n")
  
  cat("\nDistribuição de subtipos no treino:\n")
  print(table(train_selected$subtype))
  
  cat("\nDistribuição de subtipos no teste:\n")
  print(table(test_selected$subtype))
  
  cat("\nDistribuição de subtipos na validação:\n")
  print(table(validation_selected$subtype))
  
  # Configurar controle de validação cruzada
  control <- trainControl(method = "cv", 
                          number = cv_folds, 
                          savePredictions = "final")
  
  # ==================== RANDOM FOREST ====================
  cat("\n=== TREINAMENTO DO MODELO RANDOM FOREST ===\n")
  cat("Treinando modelo Random Forest com", cv_folds, "folds de validação cruzada...\n")
  
  # Treinar modelo Random Forest
  model_rf <- train(subtype ~ ., 
                    data = train_selected, 
                    method = "rf", 
                    trControl = control)
  
  cat("Modelo Random Forest treinado com sucesso!\n")
  
  # Predições RF no conjunto de teste
  cat("\n--- Random Forest - Conjunto de Teste ---\n")
  pred_test_rf <- predict(model_rf, newdata = test_selected)
  actuals_test <- factor(test_selected$subtype, levels = levels(factor(train_selected$subtype)))
  confusion_matrix_test_rf <- confusionMatrix(pred_test_rf, actuals_test)
  
  cat("Matriz de Confusão RF - Conjunto de Teste:\n")
  print(confusion_matrix_test_rf)
  
  # Predições RF no conjunto de validação
  cat("\n--- Random Forest - Conjunto de Validação ---\n")
  pred_validation_rf <- predict(model_rf, newdata = validation_selected)
  actuals_validation <- factor(validation_selected$subtype, levels = levels(factor(train_selected$subtype)))
  confusion_matrix_validation_rf <- confusionMatrix(pred_validation_rf, actuals_validation)
  
  cat("Matriz de Confusão RF - Conjunto de Validação:\n")
  print(confusion_matrix_validation_rf)
  
  # ==================== XGBOOST ====================
  cat("\n=== TREINAMENTO DO MODELO XGBOOST ===\n")
  cat("Treinando modelo XGBoost com", cv_folds, "folds de validação cruzada...\n")
  
  # Treinar modelo XGBoost
  model_xgb <- train(subtype ~ ., 
                     data = train_selected, 
                     method = "xgbTree", 
                     trControl = control)
  
  cat("Modelo XGBoost treinado com sucesso!\n")
  
  # Predições XGBoost no conjunto de teste
  cat("\n--- XGBoost - Conjunto de Teste ---\n")
  pred_test_xgb <- predict(model_xgb, newdata = test_selected, iteration_range = NULL)
  confusion_matrix_test_xgb <- confusionMatrix(pred_test_xgb, actuals_test)
  
  cat("Matriz de Confusão XGBoost - Conjunto de Teste:\n")
  print(confusion_matrix_test_xgb)
  
  # Predições XGBoost no conjunto de validação
  cat("\n--- XGBoost - Conjunto de Validação ---\n")
  pred_validation_xgb <- predict(model_xgb, newdata = validation_selected, iteration_range = NULL)
  confusion_matrix_validation_xgb <- confusionMatrix(pred_validation_xgb, actuals_validation)
  
  cat("Matriz de Confusão XGBoost - Conjunto de Validação:\n")
  print(confusion_matrix_validation_xgb) 
  
  # ==================== COMPARAÇÃO DOS MODELOS ====================
  cat("\n=== COMPARAÇÃO DOS MODELOS ===\n")
  
  # Métricas no conjunto de teste
  acc_test_rf <- confusion_matrix_test_rf$overall["Accuracy"]
  acc_test_xgb <- confusion_matrix_test_xgb$overall["Accuracy"]
  
  # Métricas no conjunto de validação
  acc_validation_rf <- confusion_matrix_validation_rf$overall["Accuracy"]
  acc_validation_xgb <- confusion_matrix_validation_xgb$overall["Accuracy"]
  
  # Tabela comparativa
  comparacao <- data.frame(
    Modelo = c("Random Forest", "XGBoost"),
    Acuracia_Teste = c(round(acc_test_rf, 4), round(acc_test_xgb, 4)),
    Acuracia_Validacao = c(round(acc_validation_rf, 4), round(acc_validation_xgb, 4))
  )
  
  cat("Comparação das Acurácias:\n")
  print(comparacao)
  
  # Determinar melhor modelo
  melhor_teste <- if (acc_test_rf > acc_test_xgb) "Random Forest" else "XGBoost"
  melhor_validacao <- if (acc_validation_rf > acc_validation_xgb) "Random Forest" else "XGBoost"
  
  cat("\nMelhor modelo no conjunto de teste:", melhor_teste, "\n")
  cat("Melhor modelo no conjunto de validação:", melhor_validacao, "\n")
  
  # Retornar resultados
  return(list(
    # Modelos treinados
    model_rf = model_rf,
    model_xgb = model_xgb,
    
    # Datasets utilizados
    train_data = train_selected,
    test_data = test_selected,
    validation_data = validation_selected,
    
    # Predições Random Forest
    predictions_test_rf = pred_test_rf,
    predictions_validation_rf = pred_validation_rf,
    
    # Predições XGBoost
    predictions_test_xgb = pred_test_xgb,
    predictions_validation_xgb = pred_validation_xgb,
    
    # Valores reais
    actuals_test = actuals_test,
    actuals_validation = actuals_validation,
    
    # Matrizes de confusão Random Forest
    confusion_matrix_test_rf = confusion_matrix_test_rf,
    confusion_matrix_validation_rf = confusion_matrix_validation_rf,
    
    # Matrizes de confusão XGBoost
    confusion_matrix_test_xgb = confusion_matrix_test_xgb,
    confusion_matrix_validation_xgb = confusion_matrix_validation_xgb,
    
    # Comparação
    comparacao_modelos = comparacao,
    melhor_modelo_teste = melhor_teste,
    melhor_modelo_validacao = melhor_validacao,
    
    # Informações gerais
    motifs_used = motifs_selecionados,
    cv_folds = cv_folds,
    prop_treino = prop_treino
  ))
}
