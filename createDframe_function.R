create_labels_dataframe <- function(seq_name, output_csv = NULL) {
  
  # Parâmetros
  n_cores <- detectCores() - 1  # leave one core free
  seq <- readBStringSet(seq_name) # Leitura de sequências
  seq_names <- names(seq) # Nome das sequências
  lengths <- width(seq) # Tamanho das sequências
  
  # Extrair subtype das primeiras letras antes do primeiro ponto
  subtypes <- subtypes <- as.character(sub("\\..*", "", seq_names))
  
  # Criar dataframe
  true_labels_hiv <- data.frame(
    name = seq_names, 
    length = lengths, 
    subtype = subtypes
  )
  
  # Gerar nome do arquivo CSV se não fornecido
  if (is.null(output_csv)) {
    # Extrair nome base do arquivo (sem extensão)
    base_name <- tools::file_path_sans_ext(basename(seq_name))
    output_csv <- paste0(base_name, "dframe.csv")
  }
  
  # Salvar CSV
  write.csv(true_labels_hiv, output_csv, row.names = FALSE)
  
  # Mensagem informativa
  cat("Dataframe created using", nrow(true_labels_hiv), "sequences\n")
  cat("File saved as:", output_csv, "\n")
  
  # Retornar o dataframe
  return(true_labels_hiv)
}




