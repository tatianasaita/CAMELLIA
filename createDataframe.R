## ---------------------------
##
## Script name: sequence-graph-adjmatrices1.R
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

## ---------------------------

## Load up our functions into memory
source("createDframe_function.R")

## ---------------------------
cat("=== DATAFRAME PREPARATION ===\n")

true_labels_hiv <- create_labels_dataframe("HIV_6964.fasta")

cat("=== createDataframe ok! ===\n")
