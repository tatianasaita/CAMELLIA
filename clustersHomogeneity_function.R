analisar_homogeneidade_subcluster <- function(dataframe) {
  
  homogeneidade_subcluster <- dataframe %>%
    group_by(subcluster) %>%
    summarise(
      # Contar total de elementos no subcluster
      total_elementos = n(),
      
      # Encontrar o subtipo mais frequente e sua quantidade
      subtipo_dominante = names(sort(table(subtype), decreasing = TRUE))[1],
      max_freq = max(table(subtype)),
      
      # Calcular homogeneidade
      homogeneidade = max_freq / total_elementos,
      
      # Listar todos os subtipos presentes
      subtipos_presentes = paste(unique(subtype), collapse = ", "),
      
      # Contar quantos subtipos diferentes
      num_subtipos = length(unique(subtype)),
      
      .groups = "drop"
    ) %>%
    arrange(subcluster)
  
  # Imprimir resultado com número total de subclusters
  total_subclusters <- length(unique(dataframe$subcluster))
  print(homogeneidade_subcluster, n = total_subclusters)
  
  # Retornar o dataframe para uso posterior
  return(homogeneidade_subcluster)
}
