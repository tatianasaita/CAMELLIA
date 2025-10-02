createNet <- function(word, step, sequence){
  aux <- ""
  index <- 1
  position <- 0
  cont <- length(sequence)
  x <- 0
  vector <- c()
  
  while((index-1+(word*2)) <= cont){
    
    while(x < word){
      aux <- paste(aux, sequence[index], sep = "")
      x <- x + 1
      index <- index + 1
    }
    vector <- c(vector, aux)
    aux <- ""
    x <- 0
    while(x < word){
      aux <- paste(aux, sequence[index], sep = "")
      x <- x + 1
      index <- index + 1
    }
    vector <- c(vector, aux)
    aux <- ""
    x <- 0
    position <- position + step
    index <- position + 1
  }
  net <- graph <- graph(edges = vector, directed = TRUE)
  
  return(net)
}
