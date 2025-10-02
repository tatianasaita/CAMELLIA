## ---------------------------
##
## Script name: createAdjMAtrix.R
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

setwd("E:/TATIANA/Git-Quali-master")

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

source("generateCombinations.R")
source("createAdjMatrix_function.R")

## ---------------------------
cat("=== GENERATE LIST OF ADJACENCY MATRICES ===\n")
# Parâmetros 
  word <- 3
  step <- 1
  vertices <- generateCombinations(word)
  seq_name <- "HIV_6964.fasta" # arquivo .fasta com as sequências
  seq <- readBStringSet(seq_name) # Leitura de sequências
  n_cores <- detectCores() - 1
  
  #Adjacency Matrix
  adj_matrices <- pblapply(seq_along(seq), function(i) {
    char_seqs <- as.character(seq[i])
    splitted_seq <- strsplit(char_seqs, split = "", fixed = TRUE)[[1]]
    createAdjMatrix(word, step, splitted_seq, vertices) 
  }, cl = n_cores)
  
  name_adj_matrices <- setNames(adj_matrices, true_labels_hiv$name) # Matrizes de adjacência nomeadas
  vectorized_matrices <- lapply(adj_matrices, function(x) as.vector(t(x))) # Transformar matrizes de adjacência em vetores
  
  print("createAdjMatrix... OK")