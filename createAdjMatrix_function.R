createAdjMatrix <- function(word, step, sequence, vertices) {
  cat("===  ===\n")
  n <- length(sequence)
  adj_matrix <- matrix(0L, nrow = length(vertices), ncol = length(vertices),
                       dimnames = list(vertices, vertices))
  
  max_index <- n - word * 2 + 1#  if (max_index <= 0) return(adj_matrix)
    indices <- seq(1, max_index, by = step)
  
  from_seq <- vapply(indices, function(i) paste(sequence[i:(i + word - 1)], collapse = ""), character(1))
  to_seq   <- vapply(indices, function(i) paste(sequence[(i + word):(i + word * 2 - 1)], collapse = ""), character(1))
  
  for (k in seq_along(from_seq)) {
    adj_matrix[from_seq[k], to_seq[k]] <- adj_matrix[from_seq[k], to_seq[k]] + 1L
  }
  
  adj_matrix
}

#############################################################
#createAdjMatrix <- function(word, step, sequence, vertices) {
 # # Se sequence for uma string, converter para vetor de caracteres
 # if(is.character(sequence) && length(sequence) == 1) {
 #   sequence <- strsplit(sequence, split = "", fixed = TRUE)[[1]]
 # }
  
 # n <- length(sequence)
 # adj_matrix <- matrix(0L, nrow = length(vertices), ncol = length(vertices),
 #                      dimnames = list(vertices, vertices))
 # 
 # max_index <- n - word * 2 + 1
 # if (max_index <= 0) return(adj_matrix)
  
 # indices <- seq(1, max_index, by = step)
  
 # from_seq <- vapply(indices, function(i) paste(sequence[i:(i + word - 1)], collapse = ""), character(1))
 # to_seq   <- vapply(indices, function(i) paste(sequence[(i + word):(i + word * 2 - 1)], collapse = ""), character(1))
  
 # for (k in seq_along(from_seq)) {
 #   adj_matrix[from_seq[k], to_seq[k]] <- adj_matrix[from_seq[k], to_seq[k]] + 1L
 # }
  
 # adj_matrix
#}
