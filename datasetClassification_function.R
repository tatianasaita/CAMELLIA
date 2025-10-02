preparar_dataset_classificacao <- function(vectorized_matrices, vertices_ordem, true_labels_hiv) {
  # Verificar inputs
  if (missing(vectorized_matrices) || missing(vertices_ordem) || missing(true_labels_hiv)) {
    stop("Todos os argumentos devem ser fornecidos: vectorized_matrices, vertices_ordem, true_labels_hiv")
  }
  
  # Criar combinações de vértices
  criar_combinacoes_vertices <- function(vertices) {
    formatted_combinations <- character()
    for (first in vertices) {
      for (second in vertices) {
        formatted_combinations <- c(formatted_combinations, paste(first, second))
      }
    }
    return(formatted_combinations)
  }
  
  # Criar combinações de vértices
  formatted_combinations <- criar_combinacoes_vertices(vertices_ordem)
  
  # Converter matrizes vetorizadas em um único dataframe
  vectorized_df <- do.call(rbind, vectorized_matrices)
  
  # Definir nomes das colunas e linhas
  colnames(vectorized_df) <- formatted_combinations
  rownames(vectorized_df) <- true_labels_hiv$name
  
  # Converter para dataframe e adicionar coluna de subtipo
  vectorized_dframe <- as.data.frame(vectorized_df)
  vectorized_dframe$subtype <- true_labels_hiv$subtype
  
  # Adicionar informações de subcluster, se disponível
  if ("subcluster" %in% names(true_labels_hiv)) {
    vectorized_dframe$subcluster <- true_labels_hiv$subcluster
  }
  
  # Verificações adicionais
  if (nrow(vectorized_dframe) != nrow(true_labels_hiv)) {
    warning("Número de linhas não corresponde entre vectorized_dframe e true_labels_hiv")
  }
  
  return(vectorized_dframe)
}

