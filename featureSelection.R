## ---------------------------
##
## Script name: featureSelection.R
##
## Purpose of script: 
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
library("dplyr")
library("parallel")
library("pbapply")

## --------------------------- 

## Load up our functions into memory

source("datasetFeatureSelection_function.R")
source("featureSelection_function.R")

## ---------------------------
cat("=== FEATURE SELECTION PREPARATION ===\n")
# Dataset of all motifs
all_motifs <- processar_matrizes_subclusters(true_labels_hiv2, name_adj_matrices)

# Features Selection
select_motifs <- select_motifs(all_motifs, clusters_homogeneity, n = 4096)
motifs_selecionados <- unname(unlist(select_motifs))

print("featureSelection... OK")
