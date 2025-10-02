plot_dendrogram <- function(vectorized_matrices, true_labels_hiv, distance_method = "euclidean", cluster_method = "ward.D2") {
  
  ##### Gerando matriz de distância
  data_matrix <- do.call(rbind, vectorized_matrices)
  dist_matriz_euc <- as.matrix(dist(data_matrix, method = distance_method))
  
  ##### Gerando dendrograma
  hc_euc <- hclust(as.dist(dist_matriz_euc), method = cluster_method)
  dend_euc <- as.dendrogram(hc_euc)
  leaf_order_euc <- order.dendrogram(dend_euc)
  species_ordered_euc <- true_labels_hiv$subtype[leaf_order_euc]
  
  # Configuração das cores
  unique_subtypes <- unique(true_labels_hiv$subtype)  
  colors <- rainbow(length(unique_subtypes))
  species_colors_euc <- colors[as.numeric(factor(species_ordered_euc, levels = unique_subtypes))]
  
  # Aplicando cores ao dendrograma
  dend_euc <- dend_euc %>% set("labels_colors", species_colors_euc)
  
  # Plotando
  plot(dend_euc)
  legend("topright", 
         legend = unique_subtypes,
         col = colors,
         pch = 15,
         title = "Subtipos")
  
  # Atribuindo objetos ao ambiente global
  assign("dend_euc", dend_euc, envir = .GlobalEnv)
  assign("dist_matriz_euc", dist_matriz_euc, envir = .GlobalEnv)
  assign("hc_euc", hc_euc, envir = .GlobalEnv)
  assign("colors", colors, envir = .GlobalEnv)
  assign("unique_subtypes", unique_subtypes, envir = .GlobalEnv)
}