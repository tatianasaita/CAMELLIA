select_motifs <- function(all_motifs, clusters_homogeneity, n) {
  # 1- Verificar se n é um número inteiro positivo
  if (!is.numeric(n) || n <= 0 || n != as.integer(n)) {
    stop("O valor de n deve ser um número inteiro positivo")
  }
  
  # 2- Calcular a quantidade de subtipos (k)
  subtipos <- unique(clusters_homogeneity$subtipo_dominante)
  k <- length(subtipos)
  print(paste("Quantidade de subtipos (k):", k))
  
  # 3- Calcular a quantidade total de clusters (m)
  clusters <- unique(clusters_homogeneity$subcluster)
  m <- as.numeric(length(clusters))
  print(paste("Quantidade total de clusters (m):", m))
  
  # 4- Calcular a quantidade de clusters em cada subtipo
  subtipo_count <- table(clusters_homogeneity$subtipo_dominante)
  print("Quantidade de clusters por subtipo:")
  print(subtipo_count)
  
  # Ordenar os subtipos em ordem decrescente de quantidade de clusters
  subtipo_order <- names(sort(subtipo_count, decreasing = TRUE))
  print("Subtipos ordenados por quantidade de clusters (decrescente):")
  for (i in 1:length(subtipo_order)) {
    subtipo <- subtipo_order[i]
    count <- subtipo_count[subtipo]
    print(paste(i, ":", subtipo, "-", count, "clusters"))
  }
  
  # Inicializar lista para armazenar os motifs selecionados
  selected_motifs <- list()
  
  # Verificar se n < k (caso em que não é possível fazer feature selection adequada)
  if (n < k) {
    print("CASO 1: n < k - Não é possível selecionar pelo menos um motif para cada subtipo.")
    return(NULL)
  } 
  
  # Caso k <= n < m
  else if (k <= n && n < m) {
    print("CASO 2: k <= n < m - Distribuindo motifs por subtipo")
    
    # Calcular quantos motifs por subtipo (base)
    motifs_per_subtipo_base <- floor(n / k)
    remainder <- n %% k
    
    # Criar vetor para armazenar a quantidade de motifs por subtipo
    motifs_per_subtipo <- setNames(rep(motifs_per_subtipo_base, length(subtipos)), subtipos)
    
    # Distribuir o resto pelos subtipos com mais clusters
    if (remainder > 0) {
      for (i in 1:remainder) {
        subtipo_to_add <- subtipo_order[((i-1) %% k) + 1]
        motifs_per_subtipo[subtipo_to_add] <- motifs_per_subtipo[subtipo_to_add] + 1
      }
    }
    
    # Mostrar a distribuição de motifs por subtipo
    print("Distribuição de motifs por subtipo:")
    for (subtipo in names(motifs_per_subtipo)) {
      print(paste("  ", subtipo, ":", motifs_per_subtipo[subtipo], "motifs"))
    }
    
    # Inicializar vetor para armazenar todos os motifs selecionados
    all_selected_motifs <- c()
    
    # Para cada subtipo, selecionar motifs de acordo com a quantidade definida
    for (subtipo in subtipos) {
      num_motifs_to_select <- motifs_per_subtipo[subtipo]
      
      # Identificar os clusters deste subtipo
      subtipo_clusters <- clusters_homogeneity$subcluster[clusters_homogeneity$subtipo_dominante == subtipo]
      cluster_cols <- paste0("Cluster_", subtipo_clusters)
      
      # Verificar se temos pelo menos um cluster para este subtipo
      if (length(cluster_cols) == 0) {
        warning(paste("Subtipo", subtipo, "não tem clusters associados"))
        selected_motifs[[subtipo]] <- character(0)
        next
      }
      
      # Criar dataframe para análise dos motifs neste subtipo
      subtipo_df <- data.frame(motif = all_motifs$Row_Column)
      
      # Para cada cluster do subtipo, adicionar coluna com os valores
      for (col in cluster_cols) {
        if (col %in% colnames(all_motifs)) {
          subtipo_df[[col]] <- all_motifs[[col]]
        } else {
          warning(paste("Coluna", col, "não encontrada em all_motifs"))
        }
      }
      
      # Avaliar todos os motifs para este subtipo usando a função de importância
      motif_scores <- lapply(1:nrow(subtipo_df), function(i) {
        evaluate_motif_importance(subtipo_df[i,], cluster_cols, all_motifs)
      })
      
      # Extrair os scores para ordenação
      mean_ranks <- sapply(motif_scores, function(x) x$mean_rank)
      top_counts <- sapply(motif_scores, function(x) x$top_count)
      
      # Adicionar scores ao dataframe
      subtipo_df$mean_rank <- mean_ranks
      subtipo_df$top_count <- top_counts
      
      # Ordenar por top_count (decrescente) e depois por mean_rank (decrescente)
      subtipo_df <- subtipo_df[order(-subtipo_df$top_count, -subtipo_df$mean_rank), ]
      
      # Selecionar os motifs que ainda não foram escolhidos
      subtipo_selected <- c()
      for (i in 1:nrow(subtipo_df)) {
        if (length(subtipo_selected) >= num_motifs_to_select) {
          break
        }
        
        current_motif <- subtipo_df$motif[i]
        if (!(current_motif %in% all_selected_motifs)) {
          subtipo_selected <- c(subtipo_selected, current_motif)
          all_selected_motifs <- c(all_selected_motifs, current_motif)
        }
      }
      
      # Verificar se conseguimos selecionar motifs suficientes
      if (length(subtipo_selected) < num_motifs_to_select) {
        warning(paste("Não foi possível selecionar", num_motifs_to_select, 
                      "motifs para o subtipo", subtipo, 
                      ". Apenas", length(subtipo_selected), "foram selecionados."))
      }
      
      # Armazenar os motifs selecionados para este subtipo
      selected_motifs[[subtipo]] <- subtipo_selected
    }
    
    # Mostrar os motifs selecionados por subtipo (formato simplificado)
    print("Resumo dos motifs selecionados por subtipo:")
    for (subtipo in names(selected_motifs)) {
      motifs <- selected_motifs[[subtipo]]
      print(paste("Subtipo", subtipo, ":", length(motifs), "motifs"))
      if (length(motifs) > 0) {
        for (i in 1:length(motifs)) {
          print(paste("   -", motifs[i]))
        }
      } else {
        print("   Nenhum motif selecionado")
      }
    }
  } 
  
  # Caso n >= m
  else {
    print("CASO 3: n >= m - Distribuindo motifs por cluster")
    
    # Calcular motifs por cluster
    motifs_per_cluster_base <- floor(n / m)
    remainder <- n %% m
    
    print(paste("Base de motifs por cluster:", motifs_per_cluster_base))
    print(paste("Resto a distribuir:", remainder))
    
    # Inicializar estrutura para armazenar motifs por subtipo
    subtipo_motifs <- setNames(vector("list", length(subtipos)), subtipos)
    all_selected_motifs <- c()
    
    # Para cada cluster, selecionar motifs_per_cluster_base motifs
    for (cluster_id in clusters) {
      cluster_col <- paste0("Cluster_", cluster_id)
      subtipo <- clusters_homogeneity$subtipo_dominante[clusters_homogeneity$subcluster == cluster_id]
      
      # Criar um dataframe temporário com nomes dos motifs e seus valores no cluster atual
      temp_df <- data.frame(
        motif = all_motifs$Row_Column,
        value = all_motifs[[cluster_col]]
      )
      
      # Ordenar por valor em ordem decrescente
      temp_df <- temp_df[order(temp_df$value, decreasing = TRUE), ]
      
      # Selecionar os top motifs que ainda não foram selecionados
      cluster_selected <- c()
      i <- 1
      while (length(cluster_selected) < motifs_per_cluster_base && i <= nrow(temp_df)) {
        current_motif <- temp_df$motif[i]
        if (!(current_motif %in% all_selected_motifs)) {
          cluster_selected <- c(cluster_selected, current_motif)
          all_selected_motifs <- c(all_selected_motifs, current_motif)
        }
        i <- i + 1
      }
      
      # Verificar se conseguimos motifs suficientes
      if (length(cluster_selected) < motifs_per_cluster_base) {
        warning(paste("Não foi possível selecionar", motifs_per_cluster_base, 
                      "motifs para o cluster", cluster_id, 
                      ". Apenas", length(cluster_selected), "foram selecionados."))
      }
      
      # Adicionar à lista de motifs do subtipo correspondente
      subtipo_motifs[[subtipo]] <- c(subtipo_motifs[[subtipo]], cluster_selected)
    }
    
    # Se houver resto, distribuir motifs extras entre os subtipos com menos clusters
    if (remainder > 0) {
      # Ordenar subtipos por quantidade de clusters (crescente)
      subtipo_order_asc <- names(sort(subtipo_count))
      
      print("Distribuindo motifs extras para os subtipos (começando pelos com menos clusters):")
      
      for (i in 1:remainder) {
        # Selecionar subtipo em ordem cíclica
        subtipo_to_add <- subtipo_order_asc[((i-1) %% k) + 1]
        print(paste("  Processando motif extra", i, "para subtipo:", subtipo_to_add))
        
        # Identificar clusters deste subtipo
        subtipo_clusters <- clusters_homogeneity$subcluster[clusters_homogeneity$subtipo_dominante == subtipo_to_add]
        cluster_cols <- paste0("Cluster_", subtipo_clusters)
        
        # Criar dataframe para análise dos motifs neste subtipo
        subtipo_df <- data.frame(motif = all_motifs$Row_Column)
        
        # Para cada cluster do subtipo, adicionar coluna com os valores
        for (col in cluster_cols) {
          subtipo_df[[col]] <- all_motifs[[col]]
        }
        
        # Avaliar todos os motifs para este subtipo
        motif_scores <- lapply(1:nrow(subtipo_df), function(i) {
          evaluate_motif_importance(subtipo_df[i,], cluster_cols, all_motifs)
        })
        
        # Extrair os scores para ordenação
        mean_ranks <- sapply(motif_scores, function(x) x$mean_rank)
        top_counts <- sapply(motif_scores, function(x) x$top_count)
        
        # Adicionar scores ao dataframe
        subtipo_df$mean_rank <- mean_ranks
        subtipo_df$top_count <- top_counts
        
        # Ordenar por top_count (decrescente) e depois por mean_rank (decrescente)
        subtipo_df <- subtipo_df[order(-subtipo_df$top_count, -subtipo_df$mean_rank), ]
        
        # Selecionar o primeiro motif que ainda não foi escolhido
        selected_extra <- NULL
        for (j in 1:nrow(subtipo_df)) {
          current_motif <- subtipo_df$motif[j]
          if (!(current_motif %in% all_selected_motifs)) {
            selected_extra <- current_motif
            break
          }
        }
        
        if (!is.null(selected_extra)) {
          print(paste("  Motif extra selecionado:", selected_extra))
          all_selected_motifs <- c(all_selected_motifs, selected_extra)
          subtipo_motifs[[subtipo_to_add]] <- c(subtipo_motifs[[subtipo_to_add]], selected_extra)
        } else {
          warning(paste("Não foi possível encontrar um motif extra para o subtipo", subtipo_to_add))
        }
      }
    }
    
    # Atribuir os resultados à lista final
    selected_motifs <- subtipo_motifs
    
    # Mostrar os motifs selecionados por subtipo
    print("Resumo dos motifs selecionados por subtipo:")
    for (subtipo in names(selected_motifs)) {
      motifs <- selected_motifs[[subtipo]]
      print(paste("Subtipo", subtipo, ":", length(motifs), "motifs"))
      if (length(motifs) > 0) {
        for (i in 1:length(motifs)) {
          print(paste("   -", motifs[i]))
        }
      } else {
        print("   Nenhum motif selecionado")
      }
    }
  } 
  
  return(selected_motifs)
}

# Função auxiliar para avaliar a importância de um motif em múltiplos clusters
evaluate_motif_importance <- function(row, cols, all_motifs) {
  values <- as.numeric(row[cols])
  
  # Para cada cluster, calcular o percentil do valor do motif
  cluster_ranks <- numeric(length(cols))
  for (i in seq_along(cols)) {
    col_values <- all_motifs[[cols[i]]]
    col_values_sorted <- sort(col_values, decreasing = TRUE)
    
    # Quanto menor o índice, melhor o ranking (1 = melhor)
    rank_position <- which(col_values_sorted == values[i])[1]
    if (is.na(rank_position)) rank_position <- length(col_values) + 1
    
    # Normalizar o ranking para um valor entre 0 e 1 (1 = melhor)
    cluster_ranks[i] <- 1 - (rank_position / length(col_values))
  }
  
  # Média dos rankings (quanto mais próximo de 1, melhor)
  mean_rank <- mean(cluster_ranks)
  
  # Contar em quantos clusters o motif está entre os top 10%
  top_count <- sum(cluster_ranks >= 0.9)
  
  # Retornar uma lista com ambas as métricas
  return(list(mean_rank = mean_rank, top_count = top_count))
}




