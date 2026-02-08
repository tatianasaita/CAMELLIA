#' Cluster Dendrogram Leaves with Homogeneity and Size Constraints
#'
#' Creates clusters from dendrogram leaves following their hierarchical order,
#' applying homogeneity and size constraints while respecting dendrogram contiguity.
#'
#'
#' @param dendrogram A dendrogram object
#' @param class_labels Character vector of class labels in dendrogram order
#' @param hom_thresh Minimum homogeneity threshold (0-1)
#' @param min_size_cluster Minimum cluster size (default: 3)
#' @param verbose Print progress messages (default: TRUE)
#' @param sequence_names Character vector of sequence names
#' @param data_result List with 'kmers' and 'metadata' for adding assignments
#'
#' @return List with clustering results
#'
#' @details
#' The clustering algorithm sequentially applies the following operations:
#' \itemize{
#'   \item Creates pure clusters for rare classes (total count < \code{min_size_cluster})
#'   \item Identifies contiguous pure-class segments of size >= \code{min_size_cluster}
#'   \item Extends clusters or forms new ones maintaining homogeneity >= \code{hom_thresh}
#'   \item Processes remaining fragments by merging with adjacent clusters or forming new ones
#'   \item Forces assignment of any unassigned elements into pure-class clusters
#'   \item Merges adjacent clusters iteratively with protection: clusters containing
#'     100\% of a class's total elements cannot be contaminated with other classes
#' }
#'
#' @note
#' \itemize{
#'   \item Requires .process_sequences and .count_kmers. See \code{internal-functions.R} for details.
#'   \item S3 methods available. See \code{methods.R} for details.
#' }
#'
#' @examples
#' \dontrun{
#' result_cluster_dendrogram <- cluster_dendrogram(
#'   dendrogram = dend,
#'   class_labels = classes,
#'   hom_thresh = 0.7,
#'   min_size_cluster = 3
#' )
#'}
#'
#' @importFrom stats is.leaf
#'
#' @export
cluster_dendrogram <- function(dendrogram,
                               class_labels,
                               hom_thresh,
                               min_size_cluster = 3L,
                               verbose = TRUE,
                               sequence_names = NULL,
                               data_result = NULL) {
# Intialization
  n_elements <- length(class_labels)
  clusters <- list()
  cluster_counter <- 0L
  assigned <- rep(FALSE, n_elements)
  class_total_counts <- table(class_labels)

  if (is.null(sequence_names)) {
    sequence_names <- paste0("seq_", seq_len(n_elements))
  }

# Handle rare classes (total count < min_size_cluster): Homogeneous cluster with classes < min_size_cluster
  rare_classes <- names(class_total_counts[class_total_counts < min_size_cluster])

    for (rare_class in rare_classes) {
      rare_indices <- which(class_labels == rare_class & !assigned)
      if (length(rare_indices) == 0) next

      segments <- .find_segments(rare_indices)

      for (seg in segments) {
        cluster_counter <- cluster_counter + 1L
        clusters[[cluster_counter]] <- .create_cluster(seg, rep(rare_class, length(seg)), cluster_counter)
        assigned[seg] <- TRUE
      }
    }

# Create pure class segments (homogeneity = 1.0): Homogeneous clusters with classes > min_size_cluster
  i <- 1L
  while (i <= n_elements) {
    if (assigned[i]) {
      i <- i + 1L
      next
    }

    current_class <- class_labels[i]
    j <- i + 1L
    while (j <= n_elements && !assigned[j] && class_labels[j] == current_class) {
      j <- j + 1L
    }

    segment_length <- j - i
    if (segment_length >= min_size_cluster) {
      cluster_counter <- cluster_counter + 1L
      idx <- i:(j - 1L)
      clusters[[cluster_counter]] <- .create_cluster(idx, class_labels[idx], cluster_counter)
      assigned[idx] <- TRUE
    }
    i <- j
  }

# Greedy expansion with homogeneity constraint: Extending last cluster or forming new cluster
  i <- 1L
  while (i <= n_elements) {
    if (assigned[i]) {
      i <- i + 1L
      next
    }

    # Try extending last cluster
    if (cluster_counter > 0L) {
      last_cluster <- clusters[[cluster_counter]]
      if (max(last_cluster$indices) == i - 1L) {
        temp_classes <- c(last_cluster$classes, class_labels[i])
        temp_hom <- .calc_homogeneity(temp_classes)

        if (temp_hom >= hom_thresh) {
          clusters[[cluster_counter]]$indices <- c(last_cluster$indices, i)
          clusters[[cluster_counter]]$classes <- temp_classes
          assigned[i] <- TRUE

          i <- i + 1L
          next
        }
      }
    }

    # Try forming new cluster
    best_end <- i
    best_hom <- 1.0

    if (i < n_elements) {
      for (j in (i + 1L):n_elements) {
        if (assigned[j]) break
        hom <- .calc_homogeneity(class_labels[i:j])
        if (hom >= hom_thresh) {
          best_end <- j
          best_hom <- hom
        } else {
          break
        }
      }
    }

    segment_size <- best_end - i + 1L
    if (segment_size >= min_size_cluster && best_hom >= hom_thresh) {
      cluster_counter <- cluster_counter + 1L
      idx <- i:best_end
      clusters[[cluster_counter]] <- .create_cluster(idx, class_labels[idx], cluster_counter)
      assigned[idx] <- TRUE
    }
    i <- best_end + 1L
  }

# Process remaining elements
  remaining <- which(!assigned)
  if (length(remaining) > 0) {
    segments <- .find_segments(remaining)

    for (seg_idx in seq_along(segments)) {
      seg <- segments[[seg_idx]]
      seg_classes <- class_labels[seg]
      n_seg <- length(seg)

      # Try merging with adjacent cluster
      merged <- FALSE
      if (cluster_counter > 0L) {
        for (cid in seq_along(clusters)) {
          if (is.null(clusters[[cid]])) next
          if (max(clusters[[cid]]$indices) == seg[1] - 1L) {
            temp_classes <- c(clusters[[cid]]$classes, seg_classes)
            temp_hom <- .calc_homogeneity(temp_classes)

            if (temp_hom >= hom_thresh) {
              clusters[[cid]]$indices <- c(clusters[[cid]]$indices, seg)
              clusters[[cid]]$classes <- temp_classes
              assigned[seg] <- TRUE
              merged <- TRUE
              break
            }
          }
        }
      }

      if (merged) next

      # Try forming clusters within segment
      if (n_seg >= min_size_cluster) {
        i <- 1L
        while (i <= n_seg) {
          if (assigned[seg[i]]) {
            i <- i + 1L
            next
          }

          best_end <- i
          best_hom <- .calc_homogeneity(seg_classes[i])

          if (i < n_seg) {
            for (j in (i + 1L):n_seg) {
              candidate_hom <- .calc_homogeneity(seg_classes[i:j])
              if (candidate_hom >= hom_thresh) {
                best_end <- j
                best_hom <- candidate_hom
              } else {
                break
              }
            }
          }

          cluster_size <- best_end - i + 1L
          if (best_hom >= hom_thresh && cluster_size >= min_size_cluster) {
            cluster_counter <- cluster_counter + 1L
            cluster_indices <- seg[i:best_end]
            clusters[[cluster_counter]] <- .create_cluster(cluster_indices, seg_classes[i:best_end], cluster_counter)
            assigned[cluster_indices] <- TRUE
          }
          i <- best_end + 1L
        }
      }
    }
  }

# Final merge attempt
  remaining <- which(!assigned)
  if (length(remaining) > 0) {
    segments <- .find_segments(remaining)

    for (seg in segments) {
      if (all(assigned[seg])) next
      seg_classes <- class_labels[seg]

      best_cid <- NULL
      best_hom <- -1

      for (cid in seq_along(clusters)) {
        if (is.null(clusters[[cid]])) next
        if (max(clusters[[cid]]$indices) == seg[1] - 1L) {
          temp_hom <- .calc_homogeneity(c(clusters[[cid]]$classes, seg_classes))
          if (temp_hom > best_hom) {
            best_cid <- cid
            best_hom <- temp_hom
          }
        }
      }

      if (!is.null(best_cid) && best_hom >= hom_thresh) {
        clusters[[best_cid]]$indices <- c(clusters[[best_cid]]$indices, seg)
        clusters[[best_cid]]$classes <- c(clusters[[best_cid]]$classes, seg_classes)
        assigned[seg] <- TRUE
      }
    }
  }

# Force cluster creation for unassigned elements
  remaining <- which(!assigned)
  if (length(remaining) > 0) {

    for (cls in unique(class_labels[remaining])) {
      cls_indices <- which(class_labels == cls & !assigned)
      if (length(cls_indices) == 0) next

      segments <- .find_segments(cls_indices)

      for (seg in segments) {
        cluster_counter <- cluster_counter + 1L
        clusters[[cluster_counter]] <- .create_cluster(seg, class_labels[seg], cluster_counter)
        assigned[seg] <- TRUE
      }
    }
  }

# Cluster merging: Clusters with 1.0 homogeneity are not merged
  if (length(clusters) > 1) {
    merged_any <- TRUE
    iteration <- 0
    max_iterations <- 10

    while (merged_any && iteration < max_iterations) {
      merged_any <- FALSE
      iteration <- iteration + 1

      i <- 1
      while (i < length(clusters)) {
        if (is.null(clusters[[i]]) || is.null(clusters[[i+1]])) {
          i <- i + 1
          next
        }

        cluster_i <- clusters[[i]]
        cluster_j <- clusters[[i+1]]

        if (max(cluster_i$indices) + 1 != min(cluster_j$indices)) {
          i <- i + 1
          next
        }

        is_i_complete <- .is_complete_class_cluster(cluster_i$classes, class_total_counts)
        is_j_complete <- .is_complete_class_cluster(cluster_j$classes, class_total_counts)

        merged_classes <- c(cluster_i$classes, cluster_j$classes)
        merged_hom <- .calc_homogeneity(merged_classes)

        can_merge <- FALSE

        if (is_i_complete && is_j_complete) {
          dom_i <- .get_dominant(cluster_i$classes)
          dom_j <- .get_dominant(cluster_j$classes)
          can_merge <- (dom_i == dom_j)

        } else if (is_i_complete) {
          dom_i <- .get_dominant(cluster_i$classes)
          can_merge <- all(cluster_j$classes == dom_i)

        } else if (is_j_complete) {
          dom_j <- .get_dominant(cluster_j$classes)
          can_merge <- all(cluster_i$classes == dom_j)

        } else {
          can_merge <- (merged_hom >= hom_thresh)
        }

        if (can_merge) {
          old_id_i <- cluster_i$id
          old_id_j <- cluster_j$id
          clusters[[i]]$indices <- c(cluster_i$indices, cluster_j$indices)
          clusters[[i]]$classes <- merged_classes
          clusters[[i+1]] <- NULL
          merged_any <- TRUE
          next
        }

        i <- i + 1
      }

      clusters <- clusters[!sapply(clusters, is.null)]
    }
  }

# Finalization
  clusters <- clusters[!sapply(clusters, is.null)]

  if (length(clusters) > 0) {
    order_idx <- order(sapply(clusters, function(cl) min(cl$indices)))
    clusters <- clusters[order_idx]
    for (k in seq_along(clusters)) clusters[[k]]$id <- k
  }

  n_clusters <- length(clusters)
  n_unassigned <- sum(!assigned)

  if (n_unassigned > 0) {
    warning("Some elements remain unassigned - this should not happen!", call. = FALSE)
  }

# Outputs
  cluster_assignment <- rep(NA_integer_, n_elements)
  for (cid in seq_along(clusters)) {
    cluster_assignment[clusters[[cid]]$indices] <- clusters[[cid]]$id
  }

  cluster_orig <- cluster_assignment

  cluster_summary <- data.frame(
    cluster_id = integer(n_clusters),
    n_elements = integer(n_clusters),
    dominant_class = character(n_clusters),
    homogeneity = numeric(n_clusters),
    n_classes = integer(n_clusters),
    class_composition = character(n_clusters),
    is_complete_class = logical(n_clusters),
    stringsAsFactors = FALSE
  )

  if (n_clusters > 0) {
    for (k in seq_len(n_clusters)) {
      cl <- clusters[[k]]
      class_table <- table(cl$classes)
      dominant_class <- names(which.max(class_table))
      is_complete <- class_table[dominant_class] == class_total_counts[dominant_class]

      cluster_summary$cluster_id[k] <- k
      cluster_summary$n_elements[k] <- length(cl$classes)
      cluster_summary$dominant_class[k] <- dominant_class
      cluster_summary$homogeneity[k] <- .calc_homogeneity(cl$classes)
      cluster_summary$n_classes[k] <- length(class_table)
      cluster_summary$class_composition[k] <- paste(names(class_table), class_table, sep = ":", collapse = "; ")
      cluster_summary$is_complete_class[k] <- is_complete
    }
  }

  element_assignment <- data.frame(
    sequence_name = sequence_names,
    dendro_index = seq_len(n_elements),
    class = class_labels,
    cluster = cluster_assignment,
    stringsAsFactors = FALSE
  )

  if (!is.null(data_result)) {
    data_result$metadata$cluster <- cluster_orig
  }

  structure(
    list(
      dendrogram = dendrogram,
      clusters = clusters,
      cluster_summary = cluster_summary,
      element_assignment = element_assignment,
      data_result = data_result,
      cluster_assignment_dendro_order = cluster_assignment,
      unassigned_elements = which(!assigned),
      n_unassigned = n_unassigned,
      hom_thresh = hom_thresh,
      min_size_cluster = min_size_cluster,
      n_elements = n_elements,
      valid_clusters = clusters
    ),
    class = "cluster_dendrogram_result"
  )
}

