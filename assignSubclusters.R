assign_subclusters <- function(hc, class_labels, hom_thresh = 0.7, min_size = 1) {
  dend <- as.dendrogram(hc)
  n_samples <- length(class_labels)
  cluster_assignments <- rep(NA, n_samples)
  cluster_id <- 1
  
  # Recursive DFS function
  dfs <- function(node) {
    # Get the indices (leaf labels) under this node
    labels <- sort(as.numeric(labels(node)))
    node_labels <- class_labels[labels]
    
    # Homogeneity check
    freq_table <- table(node_labels)
    max_prop <- max(freq_table) / sum(freq_table)
    homogeneous <- max_prop > hom_thresh
    
    # Density check
    dense <- length(labels) > min_size
    
    if (homogeneous && dense) {
      # Assign a new cluster to these indices
      cluster_assignments[labels] <<- cluster_id
      cluster_id <<- cluster_id + 1
    } else if (!is.leaf(node)) {
      dfs(node[[1]])
      dfs(node[[2]])
    }
  }
  
  dfs(dend)
  
  return(cluster_assignments)
}
