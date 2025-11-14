# ==============================================================================
# File: R/create_dendrogram.R
# ==============================================================================

# Silence R CMD check NOTEs about ggplot2 NSE variables
utils::globalVariables(c("x", "y", "xend", "yend", "label", "color"))


# ==============================================================================
# MAIN EXPORTED FUNCTION
# ==============================================================================

#' Create and plot a dendrogram from the given dataset
#'
#' This function creates a dendrogram using optimized distance calculation methods.
#' It automatically selects the most efficient algorithm based on data size and
#' available packages. For large datasets (>100 samples), it uses parallelized
#' distance calculation when available.
#'
#' The dendrogram leaves are colored according to their class using R's default
#' palette colors for optimal visual distinction.
#'
#' @param data A data frame containing k-mer counts and a 'CLASS' column.
#' @param sequence_names Character vector of sequence names (optional). If provided,
#'   these names will appear on dendrogram leaves instead of class labels.
#' @param dist_method A character string specifying the distance method (default: "euclidean").
#'   Options: "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski".
#' @param hclust_method A character string specifying the hierarchical clustering method
#'   (default: "ward.D2"). Options: "ward.D", "ward.D2", "single", "complete",
#'   "average", "mcquitty", "median", "centroid".
#' @param output A character string specifying the output file path for saving the dendrogram
#'   plot as a PNG image. If NULL, the plot will be displayed in the R graphics device.
#' @param use_parallel Logical. Use parallel processing for distance calculation when
#'   available (default: TRUE). Recommended for datasets with >100 samples.
#' @param n_threads Integer. Number of threads to use for parallel processing. If NULL,
#'   uses \code{detectCores() - 1}.
#' @param plot_title Character. Title for the dendrogram plot (default: "Dendrogram").
#' @param show_labels Logical. Show sample labels on dendrogram (default: TRUE).
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
#'     \item{height}{Numeric value representing the dendrogram height}
#'     \item{n_samples}{Integer, number of samples}
#'     \item{n_classes}{Integer, number of classes}
#'   }
#'
#' @details
#' The function uses different distance calculation methods based on data size:
#' \itemize{
#'   \item Small datasets (<100 rows): \code{stats::dist()}
#'   \item Medium datasets (100-500 rows): \code{parallelDist::parDist()} (if available)
#'   \item Large datasets (>500 rows): \code{parallelDist::parDist()} with optimized parameters
#' }
#'
#' For clustering, it uses \code{fastcluster::hclust()} when available (faster than
#' \code{stats::hclust()}), otherwise falls back to base R implementation.
#'
#' Color selection uses R's default palette and qualified color schemes for optimal
#' visualization with up to 12 classes. For more classes, rainbow colors are used.
#'
#' All plots are generated using base R graphics, with optional colored leaves when
#' the dendextend package is available.
#'
#' @examples
#' # Create sample data
#' data <- data.frame(
#'   kmer1 = rnorm(100),
#'   kmer2 = rnorm(100),
#'   CLASS = rep(c("A", "B"), each = 50)
#' )
#'
#' # Basic usage
#' result <- create_dendrogram(data)
#'
#' # With sequence names and save to file
#' \dontrun{
#' result <- create_dendrogram(
#'   data,
#'   sequence_names = paste0("seq_", 1:100),
#'   output = "dendrogram.png"
#' )
#' }
#'
#' @importFrom stats dist hclust order.dendrogram as.dendrogram
#' @importFrom grDevices png dev.off rainbow palette.colors
#' @importFrom graphics legend plot
#' @importFrom parallel detectCores
#' @export
create_dendrogram <- function(data,
                              sequence_names = NULL,
                              dist_method = "euclidean",
                              hclust_method = "ward.D2",
                              output = NULL,
                              use_parallel = TRUE,
                              n_threads = NULL,
                              plot_title = "Dendrogram",
                              show_labels = TRUE) {

  # === INPUT VALIDATION ===

  if (!is.data.frame(data)) {
    stop("'data' must be a data.frame")
  }

  if (nrow(data) == 0) {
    stop("'data' must contain at least one row")
  }

  if (!"CLASS" %in% colnames(data)) {
    stop("'data' must contain a 'CLASS' column")
  }

  if (nrow(data) < 2) {
    stop("'data' must contain at least 2 rows for clustering")
  }

  # Validate show_labels
  if (!is.logical(show_labels)) {
    stop("'show_labels' must be a logical value (TRUE or FALSE)")
  }

  # Validate plot_title
  if (!is.character(plot_title)) {
    stop("'plot_title' must be a character string")
  }

  if (length(plot_title) != 1) {
    stop("'plot_title' must be a single character string")
  }

  # Add numeric column validation
  data_numeric <- data[, !colnames(data) %in% c("CLASS"), drop = FALSE]

  # Check if there is at least one numeric column
  if (ncol(data_numeric) == 0) {
    stop("'data' must contain at least one numeric column (besides 'CLASS')")
  }

  # Check if all columns are numeric
  non_numeric_cols <- sapply(data_numeric, function(x) !is.numeric(x))
  if (any(non_numeric_cols)) {
    stop(sprintf("'data' must contain numeric columns. Non-numeric columns: %s",
                 paste(names(non_numeric_cols)[non_numeric_cols], collapse = ", ")))
  }

  # Validate sequence_names
  if (!is.null(sequence_names)) {
    if (length(sequence_names) != nrow(data)) {
      stop(sprintf("'sequence_names' must have length %d (same as nrow(data))", nrow(data)))
    }
    if (!is.character(sequence_names)) {
      stop("'sequence_names' must be a character vector")
    }
  }

  # Check for NA values
  if (any(is.na(data_numeric))) {
    stop("'data' contains NA values. Please remove or impute them before clustering")
  }

  # Check for Inf values
  if (any(is.infinite(as.matrix(data_numeric)))) {
    stop("'data' contains Inf values. Please remove them before clustering")
  }

  # Validate distance method
  valid_dist_methods <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  if (!dist_method %in% valid_dist_methods) {
    stop(sprintf("'dist_method' must be one of: %s", paste(valid_dist_methods, collapse = ", ")))
  }

  # Validate hierarchical clustering method
  valid_hclust_methods <- c("ward.D", "ward.D2", "single", "complete", "average",
                            "mcquitty", "median", "centroid")
  if (!hclust_method %in% valid_hclust_methods) {
    stop(sprintf("'hclust_method' must be one of: %s", paste(valid_hclust_methods, collapse = ", ")))
  }

  # === EXTRACT DATA ===

  class_labels <- data$CLASS
  data_matrix <- as.matrix(data[, !colnames(data) %in% c("CLASS"), drop = FALSE])
  n_samples <- nrow(data_matrix)

  # === DISTANCE CALCULATION ===

  if (use_parallel && n_samples >= 100 && requireNamespace("parallelDist", quietly = TRUE)) {
    if (is.null(n_threads)) {
      n_threads <- max(1, parallel::detectCores() - 1)
    }
    dist_matrix <- parallelDist::parDist(x = data_matrix, method = dist_method, threads = n_threads)
  } else if (n_samples >= 500 && requireNamespace("amap", quietly = TRUE)) {
    n_cores <- if (use_parallel) max(1, parallel::detectCores() - 1) else 1
    dist_matrix <- amap::Dist(x = data_matrix, method = dist_method, nbproc = n_cores)
  } else if (n_samples >= 200 && n_samples < 1000 && requireNamespace("proxy", quietly = TRUE)) {
    dist_matrix <- proxy::dist(x = data_matrix, method = dist_method)
  } else {
    dist_matrix <- stats::dist(data_matrix, method = dist_method)
  }

  rm(data_matrix)
  invisible(gc(verbose = FALSE))

  # === HIERARCHICAL CLUSTERING ===

  if (requireNamespace("fastcluster", quietly = TRUE)) {
    hc_result <- fastcluster::hclust(dist_matrix, method = hclust_method)
  } else {
    hc_result <- stats::hclust(dist_matrix, method = hclust_method)
  }

  rm(dist_matrix)
  invisible(gc(verbose = FALSE))

  # === CREATE DENDROGRAM ===

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

  # === SELECT R DEFAULT PALETTE COLORS ===
  # Use R's palette() for up to 8 classes, palette.colors() for 9-12 classes,
  # and rainbow() for more than 12 classes

  if (n_classes <= 8) {
    # Use R's default qualitative palette
    base_colors <- grDevices::palette()[1:n_classes]
  } else if (n_classes <= 12) {
    # Use R's "Paired" palette for better color distinction
    base_colors <- grDevices::palette.colors(n = n_classes, palette = "Paired")
  } else {
    # Use rainbow for large number of classes
    base_colors <- grDevices::rainbow(n_classes)
  }

  names(base_colors) <- unique_classes
  label_colors <- base_colors[as.character(class_labels[dend_order])]

  dend_height <- attr(dend, "height")

  # === APPLY COLORS TO DENDROGRAM LEAVES ===
  # Check if dendextend is available to modify labels and colors
  if (requireNamespace("dendextend", quietly = TRUE)) {
    # Use sequence_names on dendrogram leaves if provided
    if (!is.null(sequence_names_ordered)) {
      dend <- dendextend::set(dend, "labels", sequence_names_ordered)
    } else {
      dend <- dendextend::set(dend, "labels", ordered_labels)
    }

    # Apply colors to dendrogram leaves based on class
    dend <- dendextend::set(dend, "labels_col", label_colors)

  } else {
    # If dendextend is not available, use base R methods
    # Set labels using attr() - this is the safe way without dendextend
    if (!is.null(sequence_names_ordered)) {
      attr(dend, "labels") <- sequence_names_ordered
    } else {
      attr(dend, "labels") <- ordered_labels
    }

    # Warning: Colors cannot be applied without dendextend
    # The dendrogram will be plotted with default colors
    warning("Package 'dendextend' not available. Dendrogram leaves will not be colored. ",
            "Install dendextend with: install.packages('dendextend')",
            call. = FALSE)
  }

  # === PLOT DENDROGRAM (BASE R WITH COLORED LEAVES) ===

  if (is.null(output)) {
    # Plot to screen
    plot(dend, main = plot_title, ylab = "Height")

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
    plot(dend, main = plot_title, ylab = "Height")

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

  # === LABELS AND SEQUENCES ===

  numeric_data <- data[, !colnames(data) %in% c("CLASS"), drop = FALSE]

  labels_vector <- rownames(numeric_data)
  if (is.null(labels_vector)) {
    labels_vector <- paste0("S", 1:nrow(numeric_data))
  }

  sequence_names_result <- if (!is.null(sequence_names)) sequence_names_ordered else labels_vector

  # === RETURN RESULTS ===

  result <- list(
    dendrogram = dend,
    hclust = hc_result,
    order = hc_result$order,
    labels = ordered_labels,
    sequence_names = sequence_names_result,
    colors = label_colors,
    base_colors = base_colors,
    height = dend_height,
    n_samples = nrow(data),
    n_classes = n_classes
  )

  # Add S3 class
  class(result) <- c("dendrogram_result", "list")

  invisible(result)
}


# ==============================================================================
# S3 PRINT METHOD
# ==============================================================================

#' Print method for dendrogram_result objects
#'
#' Displays a summary of the dendrogram analysis result.
#'
#' @param x A dendrogram_result object created by \code{\link{create_dendrogram}}.
#' @param ... Additional arguments (not currently used).
#'
#' @return The input object \code{x}, invisibly.
#'
#' @export
#' @method print dendrogram_result
print.dendrogram_result <- function(x, ...) {
  cat("Dendrogram Analysis Result\n")
  cat("==========================\n\n")
  cat("Components:\n")
  cat("  - dendrogram: hierarchical clustering dendrogram (with colored leaves)\n")
  cat("  - hclust: hclust object\n")
  cat("  - order: dendrogram leaf order\n")
  cat("  - labels: ordered class labels\n")
  if (!is.null(x$sequence_names)) {
    cat("  - sequence_names: ordered sequence names\n")
  }
  cat("  - colors: color mapping for each sample\n")
  cat("  - base_colors: color mapping for each class\n")
  cat("  - height: maximum tree height\n")
  cat("\n")
  cat("Number of samples:", x$n_samples, "\n")
  cat("Number of classes:", x$n_classes, "\n")
  cat("Tree height:", round(x$height, 4), "\n\n")

  cat("Color mapping by class:\n")
  for (i in seq_along(x$base_colors)) {
    class_name <- names(x$base_colors)[i]
    class_color <- x$base_colors[i]
    cat(sprintf("  %s: %s\n", class_name, class_color))
  }

  invisible(x)
}
