## ---------------------------
##
## Script name: plotDendrogram.R
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
library("dendextend")
library("ggplot2")
## ---------------------------
cat("=== DENDROGRAM PREPARATION ===\n")
## Load up our functions into memory 

source("plotDendrogram_function.R")

plot_dendrogram(vectorized_matrices, true_labels_hiv)
print("plotDendrogram... OK")