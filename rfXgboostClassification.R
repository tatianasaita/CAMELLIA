## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac and PC

setwd("D:\\Quali_complete")

## ---------------------------

options(scipen = 6, digits = 4) # For non-scientific notation

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library("Biostrings") #Sequências
library("seqinr") # Sequências
library("dplyr")
library("parallel")
library("pbapply")
library("xgboost")
library("randomForest")


## ---------------------------

## Load up our functions into memory
source("rfXgboostClassification_function.R")

# Passar como parâmetro
resultado <- treinar_modelos_rf_xgboost(
  classification_dataset = datasets$classification_dataset,
  validation_dataset = datasets$validation_dataset,
  motifs_selecionados = features_selecionadas8,
  prop_treino = 0.7,
  cv_folds = 5,
 #train_control = control,  # Novo parâmetro
  seed = 123
)

resultado8 <- resultado

#resultado_valid <- treinar_modelos_rf_xgboost(
#  classification_dataset = datasets$classification_dataset,
#  validation_dataset = datasets$validation_dataset, 
#  motifs_selecionados = motifs_selecionados1024,
#  prop_treino = 0.7,
#  cv_folds = 5,
 # train_control = control,  # Novo parâmetro
#  seed = 123
#)

# 
# # Acessar os resultados
modelo_rf <- resultado$model_rf
modelo_xgb <- resultado$model_xgb
comparacao <- resultado$comparacao_modelos