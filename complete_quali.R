## ---------------------------
##
## Script name: complete_quali.R
##
## Purpose of script: Transform sequences into graphs and adjacency matrices
##
## Author: Saita, T. M.
##
## Date Created: 2025-04-28
##
## Copyright (c) Saita, T. M., 2025
## Email: tatisaita@gmail.com
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac and PC

############### ALTERAR CAMINHO!!!!

setwd("/media/tatianamarisaita/LINUX MINT/Git-Quali-master")

## ---------------------------

options(scipen = 6, digits = 4) # For non-scientific notation

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library("Biostrings") #Sequências
library("seqinr") # Sequências
library("igraph") # Grafo
library("dplyr") # Manipulação de dados
#library("pbmcapply") # Paralelização (com progress bar)
library("parallel")
library("pbapply")
library("dendextend") # Dendrogramas
library("stats") # Explícito para dist
library("ggplot2")
library("randomForest")
library("xgboost")
library("caret")
## ---------------------------

# Create dataframe of sequences information (name, length, subtype, subcluster)
  system.time({source("createDataframe.R")})
  
########################################################
# Create Net (if necessary)
# system.time({source("createNet.R")})

######################################################## 
# Create Adjacency Matrix
  system.time({source("createAdjMatrix.R")})


########################################################
# Distance matrix and dendrogram
  system.time({source("plotDendrogram.R")})


########################################################
# Clusterization
  system.time({source("automaticClusterization.R")})


########################################################
# Feature selection
  system.time({source("featureSelection.R")}) 


########################################################  
# Classification dataset (based on initial sequences )
  system.time({source("datasetClassification.R")})
       

######################################################## 
# Validation Dataset (if necessary) - Caso for usar um dataset diferente do inicial
#  system.time({source("validationDataset.R")})

  
########################################################
# Random Forest and XGBoost 
  system.time({source("rfXgboostClassification.R")})