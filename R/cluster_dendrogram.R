#' Cluster Dendrogram Leaves with Homogeneity and Size Constraints
#'
#' Creates clusters from dendrogram leaves following their hierarchical order,
#' applying homogeneity and size constraints while respecting dendrogram contiguity.
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
#' @export
cluster_dendrogram <- function(dendrogram = NULL,
                               class_labels,
                               hom_thresh = NULL,
                               min_size_cluster = 3L,
                               verbose = TRUE,
                               sequence_names = NULL,
                               data_result = NULL) {
  
  min_size_input <- as.integer(min_size_cluster)
  hom_thresh_input <- hom_thresh
  
  # ===== DENDROGRAM METHOD =====
  n_elements <- length(class_labels)
  clusters_dend <- list()
  cluster_counter_dend <- 0L
  assigned_dend <- rep(FALSE, n_elements)
  class_total_counts <- table(class_labels)
  
  if (is.null(sequence_names)) {
    sequence_names <- paste0("seq_", seq_len(n_elements))
  }
  
  # Handle rare classes
  rare_classes_dend <- names(class_total_counts[class_total_counts < min_size_input])
  
  for (rare_class in rare_classes_dend) {
    rare_indices_dend <- which(class_labels == rare_class & !assigned_dend)
    if (length(rare_indices_dend) == 0) next
    
    segments_dend <- .find_segments(rare_indices_dend)
    
    for (seg in segments_dend) {
      cluster_counter_dend <- cluster_counter_dend + 1L
      clusters_dend[[cluster_counter_dend]] <- .create_cluster(seg, rep(rare_class, length(seg)), cluster_counter_dend)
      assigned_dend[seg] <- TRUE
    }
  }
  
  # Create pure class segments
  i <- 1L
  while (i <= n_elements) {
    if (assigned_dend[i]) {
      i <- i + 1L
      next
    }
    
    current_class_dend <- class_labels[i]
    j <- i + 1L
    while (j <= n_elements && !assigned_dend[j] && class_labels[j] == current_class_dend) {
      j <- j + 1L
    }
    
    segment_length_dend <- j - i
    if (segment_length_dend >= min_size_input) {
      cluster_counter_dend <- cluster_counter_dend + 1L
      idx <- i:(j - 1L)
      clusters_dend[[cluster_counter_dend]] <- .create_cluster(idx, class_labels[idx], cluster_counter_dend)
      assigned_dend[idx] <- TRUE
    }
    i <- j
  }
  
  # Greedy expansion with homogeneity constraint
  i <- 1L
  while (i <= n_elements) {
    if (assigned_dend[i]) {
      i <- i + 1L
      next
    }
    
    if (cluster_counter_dend > 0L) {
      last_cluster_dend <- clusters_dend[[cluster_counter_dend]]
      if (max(last_cluster_dend$indices) == i - 1L) {
        temp_classes_dend <- c(last_cluster_dend$classes, class_labels[i])
        temp_hom_dend <- .calc_homogeneity(temp_classes_dend)
        
        if (temp_hom_dend >= hom_thresh_input) {
          clusters_dend[[cluster_counter_dend]]$indices <- c(last_cluster_dend$indices, i)
          clusters_dend[[cluster_counter_dend]]$classes <- temp_classes_dend
          assigned_dend[i] <- TRUE
          
          i <- i + 1L
          next
        }
      }
    }
    
    best_end_dend <- i
    best_hom_dend <- 1.0
    
    if (i < n_elements) {
      for (j in (i + 1L):n_elements) {
        if (assigned_dend[j]) break
        hom_dend <- .calc_homogeneity(class_labels[i:j])
        if (hom_dend >= hom_thresh_input) {
          best_end_dend <- j
          best_hom_dend <- hom_dend
        } else {
          break
        }
      }
    }
    
    segment_size_dend <- best_end_dend - i + 1L
    if (segment_size_dend >= min_size_input && best_hom_dend >= hom_thresh_input) {
      cluster_counter_dend <- cluster_counter_dend + 1L
      idx <- i:best_end_dend
      clusters_dend[[cluster_counter_dend]] <- .create_cluster(idx, class_labels[idx], cluster_counter_dend)
      assigned_dend[idx] <- TRUE
    }
    i <- best_end_dend + 1L
  }
  
  # Process remaining elements
  remaining_dend <- which(!assigned_dend)
  if (length(remaining_dend) > 0) {
    segments_dend <- .find_segments(remaining_dend)
    
    for (seg_idx in seq_along(segments_dend)) {
      seg <- segments_dend[[seg_idx]]
      seg_classes_dend <- class_labels[seg]
      n_seg_dend <- length(seg)
      
      merged_dend <- FALSE
      if (cluster_counter_dend > 0L) {
        for (cid in seq_along(clusters_dend)) {
          if (is.null(clusters_dend[[cid]])) next
          if (max(clusters_dend[[cid]]$indices) == seg[1] - 1L) {
            temp_classes_dend <- c(clusters_dend[[cid]]$classes, seg_classes_dend)
            temp_hom_dend <- .calc_homogeneity(temp_classes_dend)
            
            if (temp_hom_dend >= hom_thresh_input) {
              clusters_dend[[cid]]$indices <- c(clusters_dend[[cid]]$indices, seg)
              clusters_dend[[cid]]$classes <- temp_classes_dend
              assigned_dend[seg] <- TRUE
              merged_dend <- TRUE
              break
            }
          }
        }
      }
      
      if (merged_dend) next
      
      if (n_seg_dend >= min_size_input) {
        i <- 1L
        while (i <= n_seg_dend) {
          if (assigned_dend[seg[i]]) {
            i <- i + 1L
            next
          }
          
          best_end_dend <- i
          best_hom_dend <- .calc_homogeneity(seg_classes_dend[i])
          
          if (i < n_seg_dend) {
            for (j in (i + 1L):n_seg_dend) {
              candidate_hom_dend <- .calc_homogeneity(seg_classes_dend[i:j])
              if (candidate_hom_dend >= hom_thresh_input) {
                best_end_dend <- j
                best_hom_dend <- candidate_hom_dend
              } else {
                break
              }
            }
          }
          
          cluster_size_dend <- best_end_dend - i + 1L
          if (best_hom_dend >= hom_thresh_input && cluster_size_dend >= min_size_input) {
            cluster_counter_dend <- cluster_counter_dend + 1L
            cluster_indices_dend <- seg[i:best_end_dend]
            clusters_dend[[cluster_counter_dend]] <- .create_cluster(cluster_indices_dend, seg_classes_dend[i:best_end_dend], cluster_counter_dend)
            assigned_dend[cluster_indices_dend] <- TRUE
          }
          i <- best_end_dend + 1L
        }
      }
    }
  }
  
  # Final merge attempt
  remaining_dend <- which(!assigned_dend)
  if (length(remaining_dend) > 0) {
    segments_dend <- .find_segments(remaining_dend)
    
    for (seg in segments_dend) {
      if (all(assigned_dend[seg])) next
      seg_classes_dend <- class_labels[seg]
      
      best_cid_dend <- NULL
      best_hom_dend <- -1
      
      for (cid in seq_along(clusters_dend)) {
        if (is.null(clusters_dend[[cid]])) next
        if (max(clusters_dend[[cid]]$indices) == seg[1] - 1L) {
          temp_hom_dend <- .calc_homogeneity(c(clusters_dend[[cid]]$classes, seg_classes_dend))
          if (temp_hom_dend > best_hom_dend) {
            best_cid_dend <- cid
            best_hom_dend <- temp_hom_dend
          }
        }
      }
      
      if (!is.null(best_cid_dend) && best_hom_dend >= hom_thresh_input) {
        clusters_dend[[best_cid_dend]]$indices <- c(clusters_dend[[best_cid_dend]]$indices, seg)
        clusters_dend[[best_cid_dend]]$classes <- c(clusters_dend[[best_cid_dend]]$classes, seg_classes_dend)
        assigned_dend[seg] <- TRUE
      }
    }
  }
  
  # Force cluster creation for unassigned elements
  remaining_dend <- which(!assigned_dend)
  if (length(remaining_dend) > 0) {
    
    for (cls in unique(class_labels[remaining_dend])) {
      cls_indices_dend <- which(class_labels == cls & !assigned_dend)
      if (length(cls_indices_dend) == 0) next
      
      segments_dend <- .find_segments(cls_indices_dend)
      
      for (seg in segments_dend) {
        cluster_counter_dend <- cluster_counter_dend + 1L
        clusters_dend[[cluster_counter_dend]] <- .create_cluster(seg, class_labels[seg], cluster_counter_dend)
        assigned_dend[seg] <- TRUE
      }
    }
  }
  
  # Cluster merging
  if (length(clusters_dend) > 1) {
    merged_any_dend <- TRUE
    iteration_dend <- 0L
    max_iterations_dend <- 10L
    
    while (merged_any_dend && iteration_dend < max_iterations_dend) {
      merged_any_dend <- FALSE
      iteration_dend <- iteration_dend + 1L
      
      i <- 1L
      while (i < length(clusters_dend)) {
        if (is.null(clusters_dend[[i]]) || is.null(clusters_dend[[i+1L]])) {
          i <- i + 1L
          next
        }
        
        cluster_i_dend <- clusters_dend[[i]]
        cluster_j_dend <- clusters_dend[[i+1L]]
        
        if (max(cluster_i_dend$indices) + 1L != min(cluster_j_dend$indices)) {
          i <- i + 1L
          next
        }
        
        is_i_complete_dend <- .is_complete_class_cluster(cluster_i_dend$classes, class_total_counts)
        is_j_complete_dend <- .is_complete_class_cluster(cluster_j_dend$classes, class_total_counts)
        
        merged_classes_dend <- c(cluster_i_dend$classes, cluster_j_dend$classes)
        merged_hom_dend <- .calc_homogeneity(merged_classes_dend)
        
        can_merge_dend <- FALSE
        
        if (is_i_complete_dend && is_j_complete_dend) {
          dom_i_dend <- .get_dominant(cluster_i_dend$classes)
          dom_j_dend <- .get_dominant(cluster_j_dend$classes)
          can_merge_dend <- (dom_i_dend == dom_j_dend)
          
        } else if (is_i_complete_dend) {
          dom_i_dend <- .get_dominant(cluster_i_dend$classes)
          can_merge_dend <- all(cluster_j_dend$classes == dom_i_dend)
          
        } else if (is_j_complete_dend) {
          dom_j_dend <- .get_dominant(cluster_j_dend$classes)
          can_merge_dend <- all(cluster_i_dend$classes == dom_j_dend)
          
        } else {
          can_merge_dend <- (merged_hom_dend >= hom_thresh_input)
        }
        
        if (can_merge_dend) {
          clusters_dend[[i]]$indices <- c(cluster_i_dend$indices, cluster_j_dend$indices)
          clusters_dend[[i]]$classes <- merged_classes_dend
          clusters_dend[[i+1]] <- NULL
          merged_any_dend <- TRUE
          next
        }
        
        i <- i + 1L
      }
      
      clusters_dend <- clusters_dend[!sapply(clusters_dend, is.null)]
    }
  }
  
  # Finalization
  clusters_dend <- clusters_dend[!sapply(clusters_dend, is.null)]
  
  if (length(clusters_dend) > 0) {
    order_idx_dend <- order(sapply(clusters_dend, function(cl) min(cl$indices)))
    clusters_dend <- clusters_dend[order_idx_dend]
    for (k in seq_along(clusters_dend)) clusters_dend[[k]]$id <- k
  }
  
  n_clusters_dend <- length(clusters_dend)
  n_unassigned_dend <- sum(!assigned_dend)
  
  if (n_unassigned_dend > 0) {
    warning("Some elements remain unassigned - this should not happen!", call. = FALSE)
  }
  
  cluster_assignment_dend <- rep(NA_integer_, n_elements)
  for (cid in seq_along(clusters_dend)) {
    cluster_assignment_dend[clusters_dend[[cid]]$indices] <- clusters_dend[[cid]]$id
  }
  
  cluster_summary_dend <- data.frame(
    cluster_id = integer(n_clusters_dend),
    n_elements = integer(n_clusters_dend),
    dominant_class = character(n_clusters_dend),
    homogeneity = numeric(n_clusters_dend),
    n_classes = integer(n_clusters_dend),
    class_composition = character(n_clusters_dend),
    is_complete_class = logical(n_clusters_dend),
    stringsAsFactors = FALSE
  )
  
  if (n_clusters_dend > 0) {
    for (k in seq_len(n_clusters_dend)) {
      cl <- clusters_dend[[k]]
      class_table_dend <- table(cl$classes)
      dominant_class_dend <- names(which.max(class_table_dend))
      is_complete_dend <- class_table_dend[dominant_class_dend] == class_total_counts[dominant_class_dend]
      
      cluster_summary_dend$cluster_id[k] <- k
      cluster_summary_dend$n_elements[k] <- length(cl$classes)
      cluster_summary_dend$dominant_class[k] <- dominant_class_dend
      cluster_summary_dend$homogeneity[k] <- .calc_homogeneity(cl$classes)
      cluster_summary_dend$n_classes[k] <- length(class_table_dend)
      cluster_summary_dend$class_composition[k] <- paste(names(class_table_dend), class_table_dend, sep = ":", collapse = "; ")
      cluster_summary_dend$is_complete_class[k] <- is_complete_dend
    }
  }
  
  hom_thresh_observed_dend <- if (n_clusters_dend > 0) min(cluster_summary_dend$homogeneity) else NA_real_
  min_size_observed_dend <- if (n_clusters_dend > 0) min(cluster_summary_dend$n_elements) else NA_integer_
  
  element_assignment_dend <- data.frame(
    sequence_name = sequence_names,
    dendro_index = seq_len(n_elements),
    class = class_labels,
    cluster = cluster_assignment_dend,
    stringsAsFactors = FALSE
  )
  
  if (!is.null(data_result)) {
    data_result$metadata$cluster <- cluster_assignment_dend
  }
  
  structure(
    list(
      dendrogram = dendrogram,
      clusters = clusters_dend,
      cluster_summary = cluster_summary_dend,
      element_assignment = element_assignment_dend,
      data_result = data_result,
      cluster_assignment_dendro_order = cluster_assignment_dend,
      cluster_assignment = cluster_assignment_dend,
      unassigned_elements = which(!assigned_dend),
      n_unassigned = n_unassigned_dend,
      hom_thresh_input = hom_thresh_input,
      hom_thresh_observed = hom_thresh_observed_dend,
      min_size_input = min_size_input,
      min_size_observed = min_size_observed_dend,
      n_elements = n_elements,
      valid_clusters = clusters_dend
    ),
    class = "cluster_dendrogram_result"
  )
}