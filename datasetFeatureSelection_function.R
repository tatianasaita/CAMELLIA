processar_matrizes_subclusters <- function(true_labels_hiv, name_adj_matrices) {
  # Identificar subclusters únicos
  subclusters <- sort(unique(true_labels_hiv$subcluster))
  
  # Função de normalização com tratamento de casos especiais
  normalizar_matriz <- function(mat) {
    min_val <- min(mat, na.rm = TRUE)
    max_val <- max(mat, na.rm = TRUE)
    
    if (max_val > min_val) {
      return((mat - min_val) / (max_val - min_val))
    } else {
      return(mat)
    }
  }
  
  # Processamento dos subclusters
  subclusters_df <- lapply(subclusters, function(subcluster_num) {
    # Selecionar índices e nomes do subcluster
    indices_subcluster <- which(true_labels_hiv$subcluster == subcluster_num)
    names_subcluster <- true_labels_hiv$name[indices_subcluster]
    
    # Somar e normalizar matrizes do subcluster
    selected_adj_matrices <- name_adj_matrices[names_subcluster]
    summed_matrix <- Reduce(`+`, selected_adj_matrices)
    normalized_matrix <- normalizar_matriz(summed_matrix)
    
    # Converter matriz normalizada para dataframe
    subcluster_df <- as.data.frame(as.table(normalized_matrix))
    subcluster_df$Row_Column <- paste(subcluster_df$Var1, subcluster_df$Var2, sep = " ")
    subcluster_df <- subcluster_df[, c("Row_Column", "Freq")]
    
    colnames(subcluster_df)[colnames(subcluster_df) == "Freq"] <- paste("Cluster_", subcluster_num, sep = "")
    
    return(subcluster_df[, c("Row_Column", paste("Cluster_", subcluster_num, sep = ""))])
  })
  
  # Mesclar dataframes de subclusters
  all_motifs <- Reduce(function(x, y) merge(x, y, by = "Row_Column", all = TRUE), subclusters_df)
  
  # Calcular métricas adicionais
  cluster_columns <- paste0("Cluster_", subclusters)
  all_motifs$total <- rowSums(all_motifs[, cluster_columns], na.rm = TRUE)
  all_motifs$quant_cluster <- rowSums(all_motifs[, cluster_columns] != 0)
  
  # Reordenar colunas
  column_order <- c("Row_Column", 
                    paste0("Cluster_", subclusters), 
                    "total", 
                    "quant_cluster")
  
  all_motifs <- all_motifs[, column_order]
  
  # Retornar dataframe consolidado e matrizes normalizadas
  return(all_motifs)
}
