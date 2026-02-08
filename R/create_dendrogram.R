#' Create and plot a dendrogram from the given dataset
#'
#' Create a dendrogram using optimized distance calculation methods. 
#' The dendrogram leaves are colored according to their class using R's default
#' palette colors for optimal visual distinction.
#'
#' @param data Data frame containing k-mer counts and a 'CLASS' column. See \code{create_data.R} from more details.
#' @param sequence_names Character vector of sequence names (optional).
#' @param dist_method A character string specifying the distance method (default: "euclidean"). See method options in \code{\link[parallelDist]{parDist}}. 
#' @param hclust_method A character string specifying the hierarchical clustering method
#'   (default: "ward.D2"). See method options in \code{\link[stats]{hclust}}. 
#' @param output A character string specifying the output file path for saving the dendrogram
#'   plot as a PNG image. If NULL, the plot will be displayed in the R graphics device.
#'
#' @return A list (invisibly) containing:  
#'   \describe{  
#'     \item{dendrogram}{The dendrogram object (class "dendrogram") with colored leaves}  
#'     \item{hclust}{The hierarchical clustering object (class "hclust")}  
#'     \item{order}{Integer vector with the order of samples in the dendrogram}  
#'     \item{labels}{Character vector with ordered class labels}  
#'     \item{sequence_names}{Character vector with ordered sequence names (if provided)}  
#'     \item{colors}{Character vector with color mapping for each sample}  
#'     \item{base_colors}{Named character vector with color mapping for each class}  
#'   }
#'
#' @details The function create a dendrogram of sequences based on distance calculation, in the following steps:
#' \enumerate{
#'   \item Distance calculation based on k-mers count.
#'   \item Hierarchical clustering.
#'   \item Converts hclust object to dendrogram and orders samples.
#'   \item Applies colors to dendrogram leaves.
#'   \item Returns list with dendrogram object, clustering results, and color mappings.
#' }
#' @note  S3 methods available. See \code{methods.R} for details.
#'
#' @examples
#' \dontrun{
#' result <- create_dendrogram(
#'   data,
#' # sequence_names = paste0("seq_", 1:100),
#'   output = "dendrogram.png"
#' )
#' }
#'
#' @importFrom parallel detectCores
#' @importFrom parallelDist parDist
#' @importFrom fastcluster hclust
#' @importFrom stats as.dendrogram
#' @importFrom stats order.dendrogram
#' @importFrom grDevices palette.colors
#' @importFrom grDevices rainbow
#' @importFrom dendextend set
#' @importFrom grDevices png
#' @importFrom grDevices dev.off
#' 
#' @export
create_dendrogram <- function(data,
                              sequence_names = NULL,
                              dist_method = "euclidean",
                              hclust_method = "ward.D2",
                              output = NULL){

  # Extract data
  class_labels <- data$CLASS
  data_matrix <- as.matrix(data[, !colnames(data) %in% c("CLASS"), drop = FALSE])
  n_samples <- nrow(data_matrix)

  # Distance calculation
  n_threads <- max(1, parallel::detectCores() - 1)
  
  dist_matrix <- parallelDist::parDist(x = data_matrix, 
                                       method = dist_method, 
                                       threads = n_threads)
  
  rm(data_matrix) 
  invisible(gc(verbose = FALSE)) # Remove data_matrix and free memory

  # Hierarchical clustering

  hc_result <- fastcluster::hclust(dist_matrix, method = hclust_method)
  
  rm(dist_matrix)
  invisible(gc(verbose = FALSE)) # Remove data_matrix and free memory

  # Create dendrogram
  dend <- stats::as.dendrogram(hc_result)
  
  dend_order <- stats::order.dendrogram(dend)
  ordered_labels <- as.character(class_labels[dend_order])
  
  # Map sequence_names if provided
  sequence_names_ordered <- NULL
  if (!is.null(sequence_names)) {
    sequence_names_ordered <- sequence_names[dend_order]
  }
  
  # Create consistent color mapping based on original classes
  unique_classes <- unique(class_labels)
  n_classes <- length(unique_classes)

  # Select palette colors
  if (n_classes <= 12) {
    # Use R's "Paired" palette for better color distinction
    base_colors <- grDevices::palette.colors(n = n_classes, palette = "Paired") # Use palette.colors() for up to 12 classes
  } else {
    # Use rainbow for large number of classes
    base_colors <- grDevices::rainbow(n_classes) #Use and rainbow() for more than 12 classes
  }
  
  names(base_colors) <- unique_classes
  label_colors <- base_colors[as.character(class_labels[dend_order])]
  dend_height <- attr(dend, "height")

  # Apply colors to dendrogram leaves based on class
  dend <- dendextend::set(dend, "labels_col", label_colors)

  # Plot dendrogram
  if (is.null(output)) {
    # Plot to screen
    plot(dend, main = "Dendrogram", ylab = "Height")
    
    # Add legend with consistent color order
    legend("topright",
           legend = names(base_colors),
           col = base_colors,
           pch = 15,
           pt.cex = 2,
           cex = 0.8,
           title = "Classes",
           bg = "white",
           box.lty = 1)
    
  } else {
    # Save to file
    grDevices::png(filename = output, width = 800, height = 600)
    plot(dend, main = "Dendrogram", ylab = "Height")
    
    # Add legend with consistent color order
    legend("topright",
           legend = names(base_colors),
           col = base_colors,
           pch = 15,
           pt.cex = 2,
           cex = 0.8,
           title = "Classes",
           bg = "white",
           box.lty = 1)
    
    grDevices::dev.off()
  }

  # Labels and sequences
  numeric_data <- data[, !colnames(data) %in% c("CLASS"), drop = FALSE]
  
  labels_vector <- rownames(numeric_data)
  if (is.null(labels_vector)) {
    labels_vector <- paste0("S", 1:nrow(numeric_data))
  }
  
  sequence_names_result <- if (!is.null(sequence_names)) sequence_names_ordered else labels_vector

  # Return results
  result <- list(
    dendrogram = dend,
    hclust = hc_result,
    order = hc_result$order,
    labels = ordered_labels,
    sequence_names = sequence_names_result,
    colors = label_colors,
    base_colors = base_colors
  )

  # Add S3 class
  class(result) <- c("dendrogram_result", "list")

  invisible(result)
} 
