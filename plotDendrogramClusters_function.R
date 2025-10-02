plot_dendrogram_clusters <- function(dend, clusters_hiv) {
  # Color the dendrogram by clusters
  dend_colored <- color_branches(dend, clusters = true_labels_hiv2$subcluster)
  
  # Plot dendrogram
  par(mar = c(5, 5, 2, 10))  # Adjust margins
  plot(dend_colored, main = "Dendrogram with Optimized Subclusters", 
       xlab = "", sub = "", horiz = TRUE)
  # Add legend for clusters (if not too many) 
  if (length(unique(true_labels_hiv2$subcluster)) <= 60) {
    cluster_cols <- get_leaves_branches_col(dend_colored)
    cluster_ids <- sort(unique(true_labels_hiv2$subcluster))
    
    # Encontrar os índices originais das cores
    color_indices <- match(cluster_ids, unique(true_labels_hiv2$subcluster))
    
    legend("topright", 
           legend = paste("Cluster", cluster_ids), 
           fill = cluster_cols[color_indices], 
           inset = c(-0.3, 0), 
           xpd = TRUE)
  }
}
