selecionar_amostra_treino_validacao <- function(true_labels_hiv2,   
                                                classification_dataset_all,   
                                                k_por_subtipo,   
                                                #r_por_subtipo = 20,  
                                                tamanho_minimo) {  
  
  # Função para verificar disponibilidade de sequências  
  verificar_disponibilidade <- function(dados, tamanho = NULL, titulo = "") {  
    if (!is.null(tamanho)) {  
      dados_filtrados <- dados[dados$length >= tamanho, ]  
      cat(titulo, "Disponibilidade de sequências por subtipo (tamanho >=", tamanho, "):\n")  
    } else {  
      dados_filtrados <- dados  
      cat(titulo, "Disponibilidade de sequências por subtipo (qualquer tamanho):\n")  
    }  
    
    # Contar por subtipo  
    disponibilidade <- table(dados_filtrados$subtype)  
    print(disponibilidade)  
    
    return(list(  
      dados_filtrados = dados_filtrados,  
      disponibilidade = disponibilidade  
    ))  
  }  
  
  # Verificar disponibilidade inicial para classification_dataset  
  cat("=== SELEÇÃO DO CLASSIFICATION DATASET ===\n")  
  verificacao_classification <- verificar_disponibilidade(true_labels_hiv2, tamanho_minimo)  
  dados_filtrados_classification <- verificacao_classification$dados_filtrados  
  
  # Função para selecionar sequências por subtipo  
  selecionar_por_subtipo <- function(dados, k, dataset_nome = "", sequencias_excluir = NULL) {  
    # Excluir sequências já selecionadas, se houver  
    if (!is.null(sequencias_excluir)) {  
      dados <- dados[!(dados$name %in% sequencias_excluir), ]  
    }  
    
    # Agrupar por subtipo  
    dados_por_subtipo <- split(dados, dados$subtype)  
    
    # Selecionar k sequências de cada subtipo  
    sequencias_selecionadas <- lapply(dados_por_subtipo, function(subtipo_dados) {  
      n_disponivel <- nrow(subtipo_dados)  
      k_ajustado <- min(k, n_disponivel)  
      
      if (k_ajustado > 0) {  
        if (k_ajustado < k) {  
          warning(paste(dataset_nome, "- Subtipo", unique(subtipo_dados$subtype),   
                        "tem apenas", n_disponivel, "sequências disponíveis.",  
                        "Selecionando", k_ajustado, "sequências."))  
        }  
        
        # Selecionar amostra aleatória  
        indices_selecionados <- sample(nrow(subtipo_dados), k_ajustado)  
        return(subtipo_dados[indices_selecionados, ])  
      } else {  
        warning(paste(dataset_nome, "- Subtipo", unique(subtipo_dados$subtype),   
                      "não tem sequências disponíveis."))  
        return(NULL)  
      }  
    })  
    
    # Combinar resultados, removendo NULLs  
    sequencias_selecionadas <- do.call(rbind, sequencias_selecionadas[!sapply(sequencias_selecionadas, is.null)])  
    
    return(sequencias_selecionadas)  
  }  
  
  # Selecionar sequências para classification_dataset  
  true_labels_classification <- selecionar_por_subtipo(dados_filtrados_classification, k_por_subtipo, "Classification")  
  
  # Verificar se conseguiu selecionar sequências  
  if (is.null(true_labels_classification) || nrow(true_labels_classification) == 0) {  
    stop("Não foi possível selecionar nenhuma sequência para classification_dataset com os critérios especificados.")  
  }  
  
  # Extrair nomes das sequências selecionadas para classification  
  nomes_classification <- true_labels_classification$name  
  
  # Selecionar no classification_dataset_all usando rownames  
  classification_dataset <- classification_dataset_all[rownames(classification_dataset_all) %in% nomes_classification, ]  
  
  cat("\n=== SELEÇÃO DO VALIDATION DATASET ===\n")
  
  # Selecionar todas as sequências que não estão no conjunto de treino
  true_labels_validation <- true_labels_hiv2[!(true_labels_hiv2$name %in% nomes_classification), ]
  
  # Verificar se há sequências disponíveis para validação
  if (nrow(true_labels_validation) == 0) {
    warning("Não há sequências disponíveis para o conjunto de validação após remover as do conjunto de treino.")
    # Neste caso, não há como criar um conjunto de validação sem sobreposição
    true_labels_validation <- NULL
  } else {
    # Imprimir informações sobre o conjunto de validação
    cat("Validation - Disponibilidade de sequências por subtipo:\n")
    print(table(true_labels_validation$subtype))
  }
  
  # Extrair nomes das sequências selecionadas para validation
  nomes_validation <- if (!is.null(true_labels_validation)) true_labels_validation$name else character(0)
  
  # Selecionar no classification_dataset_all usando rownames
  validation_dataset <- classification_dataset_all[rownames(classification_dataset_all) %in% nomes_validation, ]
  
  # Verificar quantas sequências foram encontradas no dataset de características
  sequencias_encontradas <- sum(nomes_validation %in% rownames(classification_dataset_all))
  if (length(nomes_validation) > 0 && sequencias_encontradas < length(nomes_validation)) {
    warning(paste("Apenas", sequencias_encontradas, "de", length(nomes_validation), 
                  "sequências de validação foram encontradas no dataset de características."))
  }  
  
  # Verificações finais  
  cat("\n=== RESUMO DA SELEÇÃO ===\n")  
  
  cat("\nCLASSIFICATION DATASET:\n")  
  cat("Total de sequências selecionadas:", nrow(true_labels_classification), "\n")  
  cat("Distribuição por subtipo:\n")  
  print(table(true_labels_classification$subtype))  
  cat("Sequências encontradas no classification_dataset_all:", nrow(classification_dataset), "\n")  
  
  cat("\nVALIDATION DATASET:\n")  
  cat("Total de sequências selecionadas:", nrow(true_labels_validation), "\n")  
  cat("Distribuição por subtipo:\n")  
  print(table(true_labels_validation$subtype))  
  cat("Sequências encontradas no classification_dataset_all:", nrow(validation_dataset), "\n")  
  
  # Verificar sobreposição entre datasets  
  overlap <- intersect(nomes_classification, nomes_validation)  
  if (length(overlap) > 0) {  
    warning(paste("Há", length(overlap), "sequências sobrepostas entre os datasets!"))  
  } else {  
    cat("\nNão há sobreposição entre os datasets.\n")  
  }  
  
  # Verificar sequências não encontradas  
  nomes_classification_nao_encontrados <- setdiff(nomes_classification, rownames(classification_dataset_all))  
  nomes_validation_nao_encontrados <- setdiff(nomes_validation, rownames(classification_dataset_all))  
  
  if (length(nomes_classification_nao_encontrados) > 0) {  
    warning(paste("Sequências de classification não encontradas:", length(nomes_classification_nao_encontrados)))  
  }  
  
  if (length(nomes_validation_nao_encontrados) > 0) {  
    warning(paste("Sequências de validation não encontradas:", length(nomes_validation_nao_encontrados)))  
  }  
  
  # Salvar objetos individualmente  
  save(true_labels_classification, file = "true_labels_classification.RData")  
  save(true_labels_validation, file = "true_labels_validation.RData")  
  save(classification_dataset, file = "classification_dataset.RData")  
  save(validation_dataset, file = "validation_dataset.RData")  
  
  # Salvar checkpoint  
  save.image("backup_checkpoint1.RData")  
  
  # Retornar objetos para uso imediato  
  return(list(  
    true_labels_classification = true_labels_classification,  
    true_labels_validation = true_labels_validation,  
    classification_dataset = classification_dataset,  
    validation_dataset = validation_dataset  
  ))  
}  
