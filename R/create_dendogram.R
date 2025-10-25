#' Create and plot a dendrogram from the given dataset
#' @param data A data frame containing k-mer counts and a 'CLASS' column.
#' @param dist_method A character string specifying the distance method (default: "euclidean
#' for Euclidean distance).
#' @param hclust_method A character string specifying the hierarchical clustering method
#' (default: "ward.D2").
#' @param output A character string specifying the output file path for saving the dendrogram
#' plot as a PNG image. If NULL, the plot will be displayed in the R graphics device.
#' @return None. The function generates and optionally saves a dendrogram plot.
#' @importFrom stats dist hclust order.dendrogram
#' @importFrom grDevices png dev.off
create_dendogram <- function(data, dist_method = "euclidean", hclust_method = "ward.D2", output = NULL) {
  # Clean the data: remove 'CLASS' column
  data_clean <- data[, !colnames(data) %in% c("CLASS")]
  labels <- data$CLASS

  # Create hclust object
  dist_matrix <- stats::dist(data_clean, method = dist_method)
  hclust_obj <- stats::hclust(dist_matrix, method = hclust_method)

  # Create dendrogram object
  dendrogram_obj <- as.dendrogram(hclust_obj)
  remove(hclust_obj, dist_matrix, data_clean)

  # Order dendrogram labels
  leaf_order_euc <- stats::order.dendrogram(dendrogram_obj)
  species_ordered <- labels[leaf_order_euc]

  # Plot colors configuration
  unique_species <- unique(species_ordered)
  n_species <- length(unique_species)
  species_colors <- grDevices::rainbow(n_species)
  species_color_map <- species_colors[as.numeric(factor(species_ordered,
                                                        levels = unique_species))]

  # Setting dendrogram labels and colors
  dendrogram_obj <- dendrogram_obj %>%
    dendextend::set("labels", species_ordered)

  # Plot dendrogram
  if (is.null(output)) {
    dendrogram_obj %>%
      dendextend::set("labels_col", species_color_map) %>%
      plot(main = "Dendrogram", ylab = "Height")
  } else {
    grDevices::png(filename = output, width = 800, height = 600)
    dendrogram_obj %>%
      dendextend::set("labels_col", species_color_map) %>%
      plot(main = "Dendrogram", ylab = "Height")
    grDevices::dev.off()
  }
  # VERIFICAR QUAIS SERAO O RETURN
}
