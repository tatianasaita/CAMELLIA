## ---------------------------
##
## Script name: automaticClusterization.R
##
## Purpose of script: Tests to automatically clusterize the data
##
## Author: Saita, T. M.
##
## Date Created: 2025-05-06
##
## Copyright (c) Saita, T. M., 2025
## Email: tatisaita@gmail.com
##
## ---------------------------
##
## Notes: 
##   THIS IS A WORK IN PROGRESS VERSION
##   BE CAREFUL WITH THE CODE
##   IT IS NOT FINALIZED YET
## ---------------------------

## set working directory for Mac and PC

setwd("E:/TATIANA/Git-Quali-master")

## ---------------------------

options(scipen = 6, digits = 4) # For non-scientific notation

## ---------------------------

## load up the packages we will need:

library("dendextend")
library("ggplot2")
library("dplyr")
library("cluster")
library("colorspace")

## ---------------------------

## load up our functions into memory

source("assignSubclusters.R")
source("plotDendrogramClusters_function.R")
source("clustersHomogeneity_function.R")
## ---------------------------

cat("=== GENERATING CLUSTERS ===\n")
# Clusterization
labels_hiv <- true_labels_hiv$subtype
clusters_hiv <- assign_subclusters(dend_euc,
                                   labels_hiv,
                                   hom_thresh = 0.7,
                                   min_size = 0)

contagem_clusters <- table(clusters_hiv, labels_hiv) 
sum(contagem_clusters)

true_labels_hiv2 <- true_labels_hiv %>%
  mutate(subcluster = clusters_hiv)

clusters_homogeneity <- analisar_homogeneidade_subcluster(true_labels_hiv2) 

# Dendrogram with clusters
#dend_euc_subclusters <- plot_dendrogram_clusters(dend_euc, true_labels_hiv2$subcluster)


print("automaticClusterization... OK")