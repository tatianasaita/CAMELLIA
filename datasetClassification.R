## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac and PC

setwd("/media/tatianamarisaita/LINUX MINT/Git-Quali-master")

## ---------------------------

options(scipen = 6, digits = 4) # For non-scientific notation

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library("Biostrings") #Sequências
library("seqinr") # Sequências
library("dplyr")
library("parallel")
library("pbapply")

## ---------------------------

## Load up our functions into memory

source("datasetClassification_function.R")
source("train_test_validation_dataset_function.R")
## ---------------------------

cat("=== DATAFRAME PREPARATION ===\n")

#Dataset Classificação

classification_dataset_all <- preparar_dataset_classificacao(
  vectorized_matrices = vectorized_matrices, 
  vertices_ordem = vertices,
  true_labels_hiv = true_labels_hiv2
)

# Conjunto treino teste e validação
datasets <- selecionar_amostra_treino_validacao(
  true_labels_hiv2, 
  classification_dataset_all, 
  k_por_subtipo = 180,        # Número de sequências por subtipo no classification_dataset
#  r_por_subtipo = 20,        # Número de sequências por subtipo no validation_dataset
  tamanho_minimo = 800       # Tamanho mínimo para classification_dataset
)