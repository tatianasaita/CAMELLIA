## ---------------------------
##
## Script name: createNet.R
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

setwd("/media/tatianamarisaita/LINUX MINT/Git-Quali-master")

## ---------------------------

options(scipen = 6, digits = 4) # For non-scientific notation

## ---------------------------

## load up the packages we will need:  (uncomment as required)

library("Biostrings") #Sequências
library("seqinr") # Sequências
library("parallel")
library("pbapply")
library("dplyr")

## ---------------------------

## Load up our functions into memory
source("generateCombinations.R")
source("createNet_function.R")
## ---------------------------


# USAR APENAS QUANDO FOR NECESSÁRIO VISUALIZAR OS GRAFOS!!!


## ---------------------------

cat("=== GENERATE LIST OF GRAPHS ===\n")

# Parâmetros 
word <- 3 #ALTERAR PARÂMETROS!!
step <- 1 #ALTERAR PARÂMETROS!!
vertices <- generateCombinations(word)

# Arquivos  
seq_name <- "HIV_6964.fasta" # ALTERAR ARQUIVO .FASTA!! 
seq <- readBStringSet(seq_name) # Leitura de sequências

# Gerar e armazenar os grafos e matrizes de adjacência de cada sequência
graphs <- list()

for (i in seq_along(seq)) {
  sequence <- strsplit(toString(seq[i]), split = '')[[1]]
  net <- createNet(word, step, sequence)
  graphs[[i]] <- net
} 
cat("=== createNet ok! ===\n")
 