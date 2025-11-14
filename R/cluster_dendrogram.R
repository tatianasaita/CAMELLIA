# ==============================================================================
# File: R/cluster_dendrogram.R
# ==============================================================================

#' Cluster Dendrogram Leaves with Homogeneity and Size Constraints
#'
#' Creates clusters from dendrogram leaves following their hierarchical order,
#' applying homogeneity and minimum size constraints. The function optimizes
#' cluster formation by prioritizing complete class sequences and maximizing
#' cluster sizes while maintaining the specified homogeneity threshold.
#'
#' @param dendrogram A dendrogram object (typically from \code{\link[stats]{as.dendrogram}}
#'   or \code{\link[stats]{hclust}})
#' @param class_labels Character vector of class labels for each dendrogram element,
#'   in dendrogram order (matching the leaf order from \code{\link[stats]{order.dendrogram}})
#' @param hom_thresh Numeric value between 0 and 1. Minimum homogeneity threshold.
#'   Homogeneity is calculated as the proportion of the most frequent class:
#'   \code{max(class_counts) / total_elements}
#' @param min_size Positive integer. Minimum number of elements required per cluster
#' @param allow_flexible_start Logical. If TRUE, allows starting new clusters even when
#'   current element cannot be added to previous cluster. Default is TRUE.
#' @param verbose Logical. If TRUE, prints progress messages. Default is TRUE.
#' @param dendro_order Integer vector. Dendrogram order from create_dendrogram output.
#'   If provided, enables mapping back to original dataset order. Default is NULL.
#' @param sequence_names Character vector. Sequence names in original order.
#'   Used for detailed element assignment. Default is NULL.
#' @param data_result List containing 'kmers' and 'metadata' data.frames.
#'   If provided, adds cluster assignments to metadata. Default is NULL.
#'
#' @return A list with the following components:
#'   \describe{
#'     \item{cluster_summary}{Data.frame with cluster summary information containing:
#'       \itemize{
#'         \item cluster_id: Sequential cluster identification number
#'         \item n_elements: Total number of elements in the cluster
#'         \item n_classes: Number of distinct classes present
#'         \item dominant_class: Class with highest frequency
#'         \item homogeneity: Homogeneity value (0-1)
#'         \item class_distribution: Detailed class composition
#'       }
#'     }
#'     \item{element_assignment}{Data.frame mapping each element to its cluster (if dendro_order provided):
#'       \itemize{
#'         \item original_index: Position in original dataset
#'         \item sequence_name: Name of the sequence (if provided)
#'         \item dendro_index: Position in dendrogram order
#'         \item class: Original class label
#'         \item cluster: Assigned cluster ID
#'         \item dominant_class: Dominant class in assigned cluster
#'         \item homogeneity: Homogeneity of assigned cluster
#'       }
#'     }
#'     \item{data_result}{Updated data_result with cluster column added to metadata (if data_result provided)}
#'     \item{cluster_assignment_dendro_order}{Integer vector with cluster IDs in dendrogram order}
#'     \item{cluster_assignment_original_order}{Integer vector with cluster IDs in original order (if dendro_order provided)}
#'   }
#'
#' @details
#' The algorithm proceeds in the following steps:
#' \enumerate{
#'   \item \strong{Priority for complete classes:} If a consecutive sequence contains
#'     ALL elements of a particular class, creates a pure cluster (homogeneity = 1.0)
#'   \item \strong{Greedy expansion:} For mixed sequences, builds the largest possible
#'     clusters while maintaining homogeneity >= \code{hom_thresh}
#'   \item \strong{Smart element assignment:} When an element cannot be added to the
#'     previous cluster, starts a new cluster if possible
#'   \item \strong{Constraint enforcement:} Ensures all clusters satisfy both
#'     \code{min_size} and \code{hom_thresh} through iterative merging
#'   \item \strong{Complete assignment:} Guarantees every element is assigned to
#'     exactly one cluster
#'   \item \strong{Order mapping:} If dendro_order is provided, maps cluster assignments
#'     back to original dataset order
#'   \item \strong{Data integration:} If data_result is provided, adds cluster assignments
#'     to the metadata data.frame
#' }
#'
#' If it's impossible to satisfy both constraints simultaneously, the function will
#' attempt to create the best possible clustering and issue a warning.
#'
#' @examples
#' \dontrun{
#' # Create sample data
#' set.seed(123)
#' data_matrix <- matrix(rnorm(100 * 10), nrow = 100)
#' classes <- rep(c("A", "B", "C", "D"), each = 25)
#'
#' # Create dendrogram
#' dist_mat <- dist(data_matrix)
#' hc <- hclust(dist_mat, method = "ward.D2")
#' dend <- as.dendrogram(hc)
#'
#' # Prepare class labels in dendrogram order
#' dend_order <- order.dendrogram(dend)
#' class_labels <- classes[dend_order]
#'
#' # Create clusters
#' cluster_result <- cluster_dendrogram(
#'   dendrogram = dend,
#'   class_labels = class_labels,
#'   hom_thresh = 0.7,
#'   min_size = 5,
#'   dendro_order = dend_order
#' )
#'
#' # View results
#' print(cluster_result)
#' summary(cluster_result)
#' }
#'
#' @importFrom stats is.leaf
#' @export
cluster_dendrogram <- function(dendrogram, class_labels, hom_thresh, min_size,
                               allow_flexible_start = TRUE, verbose = TRUE,
                               dendro_order = NULL, sequence_names = NULL,
                               data_result = NULL) {

  # === Helper: Count dendrogram leaves correctly ===

  count_dendro_leaves <- function(dend) {
    if (stats::is.leaf(dend)) return(1L)
    sum(vapply(dend, count_dendro_leaves, FUN.VALUE = integer(1)))
  }

  # === INPUT VALIDATION ===

  if (!inherits(dendrogram, "dendrogram")) {
    stop("'dendrogram' must be a dendrogram object")
  }

  if (!is.character(class_labels)) {
    stop("'class_labels' must be a character vector")
  }

  # Count leaves correctly
  n_leaves <- count_dendro_leaves(dendrogram)
  if (length(class_labels) != n_leaves) {
    stop(sprintf(
      "'class_labels' must have length %d (same as dendrogram leaves)",
      n_leaves
    ))
  }

  if (any(is.na(class_labels))) {
    stop("'class_labels' contains NA values")
  }

  if (any(class_labels == "")) {
    stop("'class_labels' contains empty strings")
  }

  if (!is.numeric(hom_thresh) || length(hom_thresh) != 1 || is.na(hom_thresh) ||
      hom_thresh < 0 || hom_thresh > 1) {
    stop("'hom_thresh' must be a numeric value between 0 and 1")
  }

  if (!is.numeric(min_size) || length(min_size) != 1 || is.na(min_size) ||
      min_size < 1 || min_size != as.integer(min_size)) {
    stop("'min_size' must be an integer >= 1")
  }

  # Validate optional parameters
  if (!is.null(dendro_order)) {
    if (!is.integer(dendro_order) && !is.numeric(dendro_order)) {
      stop("'dendro_order' must be an integer vector")
    }
    if (length(dendro_order) != n_leaves) {
      stop(sprintf("'dendro_order' must have length %d", n_leaves))
    }
  }

  if (!is.null(sequence_names)) {
    if (!is.character(sequence_names)) {
      stop("'sequence_names' must be a character vector")
    }
    if (length(sequence_names) != n_leaves) {
      stop(sprintf("'sequence_names' must have length %d", n_leaves))
    }
  }

  # Validate data_result if provided
  if (!is.null(data_result)) {
    if (!is.list(data_result)) {
      stop("'data_result' must be a list")
    }

    if (!all(c("kmers", "metadata") %in% names(data_result))) {
      stop("'data_result' must contain 'kmers' and 'metadata' elements")
    }

    if (!is.data.frame(data_result$metadata)) {
      stop("'data_result$metadata' must be a data.frame")
    }

    if (nrow(data_result$metadata) != n_leaves) {
      stop(sprintf(
        "data_result has %d sequences but dendrogram has %d elements",
        nrow(data_result$metadata), n_leaves
      ))
    }

    if (is.null(dendro_order)) {
      stop("'dendro_order' must be provided when 'data_result' is provided")
    }
  }

  # === Pre-compute: Total count of each class ===

  class_total_counts <- table(class_labels)

  if (verbose) {
    message("=== Class Distribution ===")
    for (cls in names(class_total_counts)) {
      message(sprintf("%s: %d elements", cls, class_total_counts[cls]))
    }
    message("")
  }

  # === Helper functions ===

  # Optimized homogeneity calculation
  calculate_homogeneity <- function(classes) {
    if (length(classes) == 0) return(NA_real_)
    class_counts <- table(classes)
    max_count <- max(class_counts)
    return(max_count / length(classes))
  }

  # Optimized dominant class calculation
  get_dominant_class <- function(classes) {
    if (length(classes) == 0) return(NA_character_)
    class_counts <- table(classes)
    return(names(class_counts)[which.max(class_counts)])
  }

  # OPTIMIZATION B: Use rle() for faster sequence detection
  find_complete_class_sequence <- function(labels, start_pos, total_counts) {
    if (start_pos > length(labels)) {
      return(list(found = FALSE, end_pos = start_pos - 1, class_name = NA))
    }

    class_names <- names(total_counts)
    class_needed <- as.integer(total_counts[class_names])
    names(class_needed) <- class_names

    for (i in seq_along(class_names)) {
      cls <- class_names[i]
      total_needed <- class_needed[i]

      if (is.na(total_needed) || total_needed <= 0) next
      if (start_pos + total_needed - 1 > length(labels)) next

      end_check <- start_pos + total_needed - 1
      candidate_seq <- labels[start_pos:end_check]

      # Use rle for faster detection of uniform sequences
      run_lengths <- rle(candidate_seq)
      if (length(run_lengths$lengths) == 1 && run_lengths$values[1] == cls) {
        return(list(
          found = TRUE,
          end_pos = end_check,
          class_name = cls
        ))
      }
    }

    return(list(found = FALSE, end_pos = start_pos - 1, class_name = NA))
  }

  # OPTIMIZATION A: Use binary search for cluster validity
  can_create_valid_cluster <- function(labels, start_pos, end_pos_limit) {
    if (start_pos + min_size - 1 > end_pos_limit) {
      return(list(valid = FALSE, end_pos = NA))
    }

    # Binary search for the maximum valid end position
    left <- start_pos + min_size - 1
    right <- end_pos_limit
    best_end <- NA

    while (left <= right) {
      mid <- as.integer((left + right) / 2)
      candidate_classes <- labels[start_pos:mid]
      candidate_hom <- calculate_homogeneity(candidate_classes)

      if (candidate_hom >= hom_thresh) {
        best_end <- mid
        left <- mid + 1
      } else {
        right <- mid - 1
      }
    }

    if (!is.na(best_end)) {
      return(list(valid = TRUE, end_pos = best_end))
    } else {
      return(list(valid = FALSE, end_pos = NA))
    }
  }

  # === Main clustering algorithm ===

  n_elements <- length(class_labels)
  clusters <- list()
  cluster_id <- 0
  i <- 1

  while (i <= n_elements) {

    # Step 1: Check if there's a COMPLETE class sequence starting at i
    complete_check <- find_complete_class_sequence(class_labels, i, class_total_counts)

    if (complete_check$found) {
      cluster_id <- cluster_id + 1
      clusters[[cluster_id]] <- list(
        indices = i:complete_check$end_pos,
        classes = class_labels[i:complete_check$end_pos]
      )

      if (verbose) {
        message(sprintf(">> Created PURE cluster for class '%s' with ALL %d elements (homog=1.0)",
                        complete_check$class_name,
                        complete_check$end_pos - i + 1))
      }

      i <- complete_check$end_pos + 1
      next
    }

    # Step 2: Check if we have enough elements remaining
    remaining_elements <- n_elements - i + 1

    if (remaining_elements < min_size) {
      if (cluster_id > 0) {
        remaining <- i:n_elements
        merged_classes <- c(clusters[[cluster_id]]$classes, class_labels[remaining])
        merged_hom <- calculate_homogeneity(merged_classes)

        if (merged_hom >= hom_thresh) {
          clusters[[cluster_id]]$indices <- c(clusters[[cluster_id]]$indices, remaining)
          clusters[[cluster_id]]$classes <- merged_classes
          if (verbose) {
            message(sprintf(">> Appended %d remaining elements to cluster %d (new homog=%.3f)",
                            length(remaining), cluster_id, merged_hom))
          }
        } else {
          warning(sprintf(
            "Last %d elements cannot satisfy homogeneity constraint. Appending anyway with homogeneity=%.3f (threshold=%.3f)",
            length(remaining), merged_hom, hom_thresh
          ))
          clusters[[cluster_id]]$indices <- c(clusters[[cluster_id]]$indices, remaining)
          clusters[[cluster_id]]$classes <- merged_classes
        }
      } else {
        cluster_id <- 1
        clusters[[cluster_id]] <- list(
          indices = i:n_elements,
          classes = class_labels[i:n_elements]
        )
        warning(sprintf(
          "Creating single cluster with %d elements (less than min_size=%d)",
          remaining_elements, min_size
        ))
      }
      break
    }

    # Step 3: Try to create a new valid cluster starting at position i
    cluster_check <- can_create_valid_cluster(class_labels, i, n_elements)

    if (cluster_check$valid) {
      cluster_id <- cluster_id + 1
      clusters[[cluster_id]] <- list(
        indices = i:cluster_check$end_pos,
        classes = class_labels[i:cluster_check$end_pos]
      )

      hom_val <- calculate_homogeneity(class_labels[i:cluster_check$end_pos])
      if (verbose) {
        message(sprintf(">> Created cluster %d: elements %d-%d (size=%d, homog=%.3f)",
                        cluster_id, i, cluster_check$end_pos,
                        cluster_check$end_pos - i + 1, hom_val))
      }

      i <- cluster_check$end_pos + 1
    } else {
      if (cluster_id > 0) {
        merged_classes <- c(clusters[[cluster_id]]$classes, class_labels[i])
        merged_hom <- calculate_homogeneity(merged_classes)

        if (merged_hom >= hom_thresh) {
          clusters[[cluster_id]]$indices <- c(clusters[[cluster_id]]$indices, i)
          clusters[[cluster_id]]$classes <- merged_classes
          i <- i + 1
        } else {
          if (allow_flexible_start) {
            found_valid <- FALSE
            for (look_ahead in (i + 1):min(i + 20, n_elements)) {
              ahead_check <- can_create_valid_cluster(class_labels, look_ahead, n_elements)
              if (ahead_check$valid) {
                problematic <- i:(look_ahead - 1)
                merged_classes <- c(clusters[[cluster_id]]$classes, class_labels[problematic])
                merged_hom <- calculate_homogeneity(merged_classes)

                warning(sprintf(
                  "Elements %d-%d added to cluster %d with reduced homogeneity %.3f (threshold=%.3f)",
                  i, look_ahead - 1, cluster_id, merged_hom, hom_thresh
                ))

                clusters[[cluster_id]]$indices <- c(clusters[[cluster_id]]$indices, problematic)
                clusters[[cluster_id]]$classes <- merged_classes
                i <- look_ahead
                found_valid <- TRUE
                break
              }
            }

            if (!found_valid) {
              warning(sprintf(
                "Element %d (class '%s') added to cluster %d, resulting homogeneity=%.3f (threshold=%.3f)",
                i, class_labels[i], cluster_id, merged_hom, hom_thresh
              ))
              clusters[[cluster_id]]$indices <- c(clusters[[cluster_id]]$indices, i)
              clusters[[cluster_id]]$classes <- merged_classes
              i <- i + 1
            }
          } else {
            stop(sprintf(
              "Cannot assign element at position %d without violating homogeneity constraint. Current element class: '%s'. Previous cluster homogeneity would become %.3f (threshold: %.3f). Consider reducing hom_thresh or min_size, or set allow_flexible_start=TRUE.",
              i, class_labels[i], merged_hom, hom_thresh
            ))
          }
        }
      } else {
        cluster_id <- 1
        end_pos <- min(i + min_size - 1, n_elements)
        clusters[[cluster_id]] <- list(
          indices = i:end_pos,
          classes = class_labels[i:end_pos]
        )

        hom_val <- calculate_homogeneity(class_labels[i:end_pos])
        if (hom_val < hom_thresh) {
          warning(sprintf(
            "First cluster has homogeneity %.3f (threshold=%.3f)",
            hom_val, hom_thresh
          ))
        }

        i <- end_pos + 1
      }
    }
  }

  # === Post-processing: Merge small clusters ===

  k <- 1
  while (k <= length(clusters)) {
    current_size <- length(clusters[[k]]$indices)
    current_hom <- calculate_homogeneity(clusters[[k]]$classes)

    if (current_size < min_size) {
      if (k > 1) {
        merged_classes <- c(clusters[[k - 1]]$classes, clusters[[k]]$classes)
        merged_hom <- calculate_homogeneity(merged_classes)

        clusters[[k - 1]]$indices <- c(clusters[[k - 1]]$indices, clusters[[k]]$indices)
        clusters[[k - 1]]$classes <- merged_classes
        clusters[[k]] <- NULL
        clusters <- clusters[!vapply(clusters, is.null, FUN.VALUE = logical(1))]

        if (merged_hom < hom_thresh) {
          warning(sprintf(
            "Merged cluster %d has homogeneity %.3f (threshold=%.3f)",
            k - 1, merged_hom, hom_thresh
          ))
        }
      } else if (k < length(clusters)) {
        merged_classes <- c(clusters[[k]]$classes, clusters[[k + 1]]$classes)
        merged_hom <- calculate_homogeneity(merged_classes)

        clusters[[k + 1]]$indices <- c(clusters[[k]]$indices, clusters[[k + 1]]$indices)
        clusters[[k + 1]]$classes <- merged_classes
        clusters[[k]] <- NULL
        clusters <- clusters[!vapply(clusters, is.null, FUN.VALUE = logical(1))]

        if (merged_hom < hom_thresh) {
          warning(sprintf(
            "Merged cluster %d has homogeneity %.3f (threshold=%.3f)",
            k, merged_hom, hom_thresh
          ))
        }
      } else {
        warning(sprintf(
          "Single cluster with size %d (min_size=%d)",
          current_size, min_size
        ))
        k <- k + 1
      }
    } else {
      k <- k + 1
    }
  }

  # === OPTIMIZATION C: Pre-allocate data.frame using vectorized operations ===

  n_clusters <- length(clusters)

  cluster_ids <- integer(n_clusters)
  n_elements_vec <- integer(n_clusters)
  n_classes_vec <- integer(n_clusters)
  dominant_classes <- character(n_clusters)
  homogeneities <- numeric(n_clusters)
  class_distributions <- character(n_clusters)

  for (k in seq_len(n_clusters)) {
    cluster_info <- clusters[[k]]
    classes <- cluster_info$classes
    n_elem <- length(classes)
    dom_class <- get_dominant_class(classes)
    homogeneity <- calculate_homogeneity(classes)

    class_counts <- table(classes)
    n_classes <- length(class_counts)

    class_distribution <- paste(
      vapply(names(class_counts), function(cls) {
        sprintf("%s(%d)", cls, class_counts[cls])
      }, FUN.VALUE = character(1)),
      collapse = ", "
    )

    cluster_ids[k] <- k
    n_elements_vec[k] <- n_elem
    n_classes_vec[k] <- n_classes
    dominant_classes[k] <- dom_class
    homogeneities[k] <- homogeneity
    class_distributions[k] <- class_distribution
  }

  cluster_summary <- data.frame(
    cluster_id = cluster_ids,
    n_elements = n_elements_vec,
    n_classes = n_classes_vec,
    dominant_class = dominant_classes,
    homogeneity = homogeneities,
    class_distribution = class_distributions,
    stringsAsFactors = FALSE
  )

  # === Final Validation ===

  total_elements <- sum(cluster_summary$n_elements)
  if (total_elements != n_elements) {
    stop(sprintf(
      "INTERNAL ERROR: Assignment mismatch - expected %d elements, got %d",
      n_elements, total_elements
    ))
  }

  # === Build cluster assignment vectors ===

  cluster_assignment_dendro_order <- rep(NA_integer_, n_elements)
  current_pos <- 1

  for (row_idx in seq_len(nrow(cluster_summary))) {
    cluster_id <- cluster_summary$cluster_id[row_idx]
    n_elem <- cluster_summary$n_elements[row_idx]
    end_pos <- current_pos + n_elem - 1
    cluster_assignment_dendro_order[current_pos:end_pos] <- cluster_id
    current_pos <- end_pos + 1
  }

  # Map back to original order if dendro_order provided
  cluster_assignment_original_order <- NULL
  if (!is.null(dendro_order)) {
    cluster_assignment_original_order <- rep(NA_integer_, n_elements)
    cluster_assignment_original_order[dendro_order] <- cluster_assignment_dendro_order
  }

  # === Create detailed element assignment ===

  element_assignment <- NULL
  if (!is.null(dendro_order)) {
    element_assignment <- data.frame(
      original_index = 1:n_elements,
      sequence_name = if (!is.null(sequence_names)) sequence_names else paste0("seq_", 1:n_elements),
      dendro_index = dendro_order,
      class = class_labels,
      cluster = cluster_assignment_original_order,
      stringsAsFactors = FALSE
    )

    # OPTIMIZATION D: Use match() instead of merge() for better performance
    idx <- match(element_assignment$cluster, cluster_summary$cluster_id)
    element_assignment$dominant_class <- cluster_summary$dominant_class[idx]
    element_assignment$homogeneity <- cluster_summary$homogeneity[idx]

    element_assignment <- element_assignment[order(element_assignment$original_index), ]
    rownames(element_assignment) <- NULL
  }

  # === Add cluster to data_result metadata ===

  if (!is.null(data_result)) {
    data_result$metadata$cluster <- cluster_assignment_original_order
  }

  # === Summary messages ===

  n_below_minsize <- sum(cluster_summary$n_elements < min_size)
  n_below_hom <- sum(cluster_summary$homogeneity < hom_thresh)

  if (verbose) {
    message("\n=== Dendrogram Clustering Summary ===")
    message(sprintf("Total elements: %d", n_elements))
    message(sprintf("Total clusters: %d", nrow(cluster_summary)))
    message(sprintf("Pure clusters (homog=1.0): %d", sum(cluster_summary$homogeneity == 1.0)))
    message(sprintf("Homogeneity threshold: %.2f", hom_thresh))
    message(sprintf("Minimum cluster size: %d", min_size))
    message(sprintf("Average homogeneity: %.4f", round(mean(cluster_summary$homogeneity), 4)))
    message(sprintf("Min/Max cluster sizes: %d / %d",
                    min(cluster_summary$n_elements), max(cluster_summary$n_elements)))
    message(sprintf("Min/Max homogeneity: %.4f / %.4f",
                    min(cluster_summary$homogeneity), max(cluster_summary$homogeneity)))

    if (n_below_hom > 0) {
      message(sprintf("WARNING: %d cluster(s) below homogeneity threshold", n_below_hom))
    }
    if (n_below_minsize > 0) {
      message(sprintf("WARNING: %d cluster(s) below minimum size", n_below_minsize))
    }

    if (!is.null(dendro_order)) {
      message("\n=== Cluster Assignment Summary ===")
      message(sprintf("Total sequences assigned: %d", n_elements))
      message("Cluster distribution:")
      cluster_dist <- table(cluster_assignment_original_order)
      for (cid in names(cluster_dist)) {
        message(sprintf("  Cluster %s: %d sequences", cid, cluster_dist[cid]))
      }
    }

    if (!is.null(data_result)) {
      message("\nCluster column successfully added to data_result$metadata")
    }

    message("")
  }

  # === Return results ===

  result <- list(
    cluster_summary = cluster_summary,
    element_assignment = element_assignment,
    data_result = data_result,
    cluster_assignment_dendro_order = cluster_assignment_dendro_order,
    cluster_assignment_original_order = cluster_assignment_original_order
  )

  class(result) <- c("cluster_dendrogram_result", "list")

  invisible(result)
}


# ==============================================================================
# S3 PRINT METHOD
# ==============================================================================

#' Print method for cluster_dendrogram_result objects
#'
#' Displays a summary of the dendrogram clustering analysis result.
#'
#' @param x A cluster_dendrogram_result object created by \code{\link{cluster_dendrogram}}.
#' @param ... Additional arguments (not currently used).
#'
#' @return The input object \code{x}, invisibly.
#'
#' @method print cluster_dendrogram_result
#' @export
print.cluster_dendrogram_result <- function(x, ...) {
  cat("Dendrogram Clustering Result\n")
  cat("============================\n\n")

  cat("Number of clusters: ", nrow(x$cluster_summary), "\n")
  cat("Total elements:     ", sum(x$cluster_summary$n_elements), "\n")

  if (!is.null(x$cluster_summary)) {
    mean_hom <- mean(x$cluster_summary$homogeneity, na.rm = TRUE)
    cat("Mean homogeneity:   ", round(mean_hom, 3), "\n")
  }

  cat("\nCluster Summary:\n")
  cat(strrep("-", 80), "\n")

  # Imprimir sem parâmetro 'n'
  print(x$cluster_summary, na.print = "")

  cat(strrep("-", 80), "\n")
  invisible(x)
}


# ==============================================================================
# S3 SUMMARY METHOD
# ==============================================================================

#' Summary method for cluster_dendrogram_result objects
#'
#' Provides detailed statistical summary of the dendrogram clustering result.
#'
#' @param object A cluster_dendrogram_result object created by \code{\link{cluster_dendrogram}}.
#' @param ... Additional arguments (not currently used).
#'
#' @return The object invisibly with additional statistics printed.
#'
#' @export
#' @method summary cluster_dendrogram_result
summary.cluster_dendrogram_result <- function(object, ...) {
  cat("Dendrogram Clustering Analysis Summary\n")
  cat("======================================\n\n")

  cat("Cluster Statistics:\n")
  cat("  Number of clusters: ", nrow(object$cluster_summary), "\n")
  cat("  Total elements: ", sum(object$cluster_summary$n_elements), "\n")
  cat("  Pure clusters (homogeneity = 1.0): ",
      sum(object$cluster_summary$homogeneity == 1.0), "\n\n")

  cat("Cluster Size Distribution:\n")
  size_summary <- summary(object$cluster_summary$n_elements)
  print(size_summary)

  cat("\nHomogeneity Distribution:\n")
  hom_summary <- summary(object$cluster_summary$homogeneity)
  print(hom_summary)

  cat("\nClass Distribution Across Clusters:\n")
  if (!is.null(object$element_assignment)) {
    class_table <- table(object$element_assignment$class, object$element_assignment$cluster)
    print(class_table)
  } else {
    cat("  (Not available - dendro_order was not provided)\n")
  }

  cat("\nCluster Details:\n")
  for (i in seq_len(nrow(object$cluster_summary))) {
    row <- object$cluster_summary[i, ]
    cat(sprintf("  Cluster %d: %d elements, %d class(es), %s (homog=%.3f)\n",
                row$cluster_id, row$n_elements, row$n_classes,
                row$dominant_class, row$homogeneity))
  }

  invisible(object)
}
