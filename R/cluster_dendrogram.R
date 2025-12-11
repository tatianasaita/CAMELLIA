#' Cluster Dendrogram Leaves with Homogeneity and Size Constraints
#'
#' Creates clusters from dendrogram leaves following their hierarchical order,
#' applying homogeneity and size constraints while respecting dendrogram contiguity.
#'
#' @param dendrogram A dendrogram object
#' @param class_labels Character vector of class labels in dendrogram order
#' @param hom_thresh Minimum homogeneity threshold (0-1)
#' @param min_size Minimum cluster size (default: 3)
#' @param verbose Print progress messages (default: TRUE)
#' @param dendro_order Integer vector for mapping back to original order
#' @param sequence_names Character vector of sequence names
#' @param data_result List with 'kmers' and 'metadata' for adding assignments
#'
#' @return List with clustering results
#' @export
cluster_dendrogram <- function(dendrogram,
                               class_labels,
                               hom_thresh,
                               min_size = 3L,
                               verbose = TRUE,
                               dendro_order = NULL,
                               sequence_names = NULL,
                               data_result = NULL) {

  # === HELPER FUNCTIONS ===

  count_leaves <- function(dend) {
    if (stats::is.leaf(dend)) return(1L)
    sum(vapply(dend, count_leaves, FUN.VALUE = integer(1)))
  }

  calc_homogeneity <- function(classes) {
    if (length(classes) == 0) return(NA_real_)
    max(table(classes)) / length(classes)
  }

  get_dominant <- function(classes) {
    if (length(classes) == 0) return(NA_character_)
    names(which.max(table(classes)))
  }

  find_segments <- function(indices) {
    if (length(indices) <= 1) return(list(indices))
    segs <- list()
    current <- indices[1]
    for (i in 2:length(indices)) {
      if (indices[i] == indices[i-1] + 1) {
        current <- c(current, indices[i])
      } else {
        segs[[length(segs) + 1]] <- current
        current <- indices[i]
      }
    }
    segs[[length(segs) + 1]] <- current
    segs
  }

  is_complete_class_cluster <- function(cluster_classes, total_counts) {
    class_table <- table(cluster_classes)
    dominant_class <- names(which.max(class_table))
    class_table[dominant_class] == total_counts[dominant_class]
  }

  create_cluster <- function(indices, classes, counter) {
    list(indices = indices, classes = classes, id = counter)
  }

  log_msg <- function(...) if (verbose) message(...)

  # === VALIDATION ===

  if (!inherits(dendrogram, "dendrogram")) {
    stop("'dendrogram' must be a dendrogram object", call. = FALSE)
  }

  n_elements <- count_leaves(dendrogram)

  if (length(class_labels) != n_elements) {
    stop(sprintf("'class_labels' length must be %d", n_elements), call. = FALSE)
  }

  if (anyNA(class_labels) || any(class_labels == "")) {
    stop("'class_labels' contains NA or empty values", call. = FALSE)
  }

  min_size <- as.integer(min_size)

  log_msg(sprintf("Starting clustering (n=%d, min_size=%d, hom_thresh=%.2f)",
                  n_elements, min_size, hom_thresh))

  # === INITIALIZATION ===

  clusters <- list()
  cluster_counter <- 0L
  assigned <- rep(FALSE, n_elements)
  class_total_counts <- table(class_labels)

  # === PHASE 0: Handle rare classes (total count < min_size) ===

  log_msg("Phase 0: Identifying rare classes (total count < min_size)")

  rare_classes <- names(class_total_counts[class_total_counts < min_size])

  if (length(rare_classes) > 0) {
    log_msg(sprintf("  Found %d rare classes: %s",
                    length(rare_classes), paste(rare_classes, collapse = ", ")))

    for (rare_class in rare_classes) {
      rare_indices <- which(class_labels == rare_class & !assigned)
      if (length(rare_indices) == 0) next

      segments <- find_segments(rare_indices)
      log_msg(sprintf("  Class '%s' (%d total elements, %d segments)",
                      rare_class, length(rare_indices), length(segments)))

      for (seg in segments) {
        cluster_counter <- cluster_counter + 1L
        clusters[[cluster_counter]] <- create_cluster(seg, rep(rare_class, length(seg)), cluster_counter)
        assigned[seg] <- TRUE

        log_msg(sprintf("    Created %s cluster %d at position%s %s",
                        if (length(seg) == 1) "singleton" else "",
                        cluster_counter,
                        if (length(seg) > 1) "s" else "",
                        if (length(seg) == 1) seg else sprintf("%d-%d", min(seg), max(seg))))
      }
    }
    log_msg(sprintf("  Phase 0 complete: %d clusters created for rare classes\n", cluster_counter))
  } else {
    log_msg("  No rare classes found\n")
  }

  # === PHASE 1: Create pure class segments (homogeneity = 1.0) ===

  log_msg("Phase 1: Creating pure class segments")

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
    if (segment_length >= min_size) {
      cluster_counter <- cluster_counter + 1L
      idx <- i:(j - 1L)
      clusters[[cluster_counter]] <- create_cluster(idx, class_labels[idx], cluster_counter)
      assigned[idx] <- TRUE
      log_msg(sprintf("  Created cluster %d: %d elements of class '%s' (positions %d-%d)",
                      cluster_counter, segment_length, current_class, i, j-1))
    }
    i <- j
  }

  # === PHASE 2: Greedy expansion with homogeneity constraint ===

  log_msg("\nPhase 2: Greedy expansion")

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
        temp_hom <- calc_homogeneity(temp_classes)

        if (temp_hom >= hom_thresh) {
          clusters[[cluster_counter]]$indices <- c(last_cluster$indices, i)
          clusters[[cluster_counter]]$classes <- temp_classes
          assigned[i] <- TRUE
          log_msg(sprintf("  Extended cluster %d (size=%d, hom=%.3f)",
                          cluster_counter, length(temp_classes), temp_hom))
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
        hom <- calc_homogeneity(class_labels[i:j])
        if (hom >= hom_thresh) {
          best_end <- j
          best_hom <- hom
        } else {
          break
        }
      }
    }

    segment_size <- best_end - i + 1L
    if (segment_size >= min_size && best_hom >= hom_thresh) {
      cluster_counter <- cluster_counter + 1L
      idx <- i:best_end
      clusters[[cluster_counter]] <- create_cluster(idx, class_labels[idx], cluster_counter)
      assigned[idx] <- TRUE
      log_msg(sprintf("  Created cluster %d: %d elements, hom=%.3f (positions %d-%d)",
                      cluster_counter, segment_size, best_hom, i, best_end))
    }
    i <- best_end + 1L
  }

  # === PHASE 3: Process remaining elements ===

  remaining <- which(!assigned)
  if (length(remaining) > 0) {
    log_msg(sprintf("\nPhase 3: Processing %d remaining elements", length(remaining)))
    segments <- find_segments(remaining)

    for (seg_idx in seq_along(segments)) {
      seg <- segments[[seg_idx]]
      seg_classes <- class_labels[seg]
      n_seg <- length(seg)

      log_msg(sprintf("  Segment %d: %d elements at positions %d-%d",
                      seg_idx, n_seg, min(seg), max(seg)))

      # Try merging with adjacent cluster
      merged <- FALSE
      if (cluster_counter > 0L) {
        for (cid in seq_along(clusters)) {
          if (is.null(clusters[[cid]])) next
          if (max(clusters[[cid]]$indices) == seg[1] - 1L) {
            temp_classes <- c(clusters[[cid]]$classes, seg_classes)
            temp_hom <- calc_homogeneity(temp_classes)

            if (temp_hom >= hom_thresh) {
              clusters[[cid]]$indices <- c(clusters[[cid]]$indices, seg)
              clusters[[cid]]$classes <- temp_classes
              assigned[seg] <- TRUE
              merged <- TRUE
              log_msg(sprintf("    [OK] Merged into cluster %d (size=%d, hom=%.3f)",
                              clusters[[cid]]$id, length(temp_classes), temp_hom))
              break
            }
          }
        }
      }

      if (merged) next

      # Try forming clusters within segment
      if (n_seg >= min_size) {
        i <- 1L
        while (i <= n_seg) {
          if (assigned[seg[i]]) {
            i <- i + 1L
            next
          }

          best_end <- i
          best_hom <- calc_homogeneity(seg_classes[i])

          if (i < n_seg) {
            for (j in (i + 1L):n_seg) {
              candidate_hom <- calc_homogeneity(seg_classes[i:j])
              if (candidate_hom >= hom_thresh) {
                best_end <- j
                best_hom <- candidate_hom
              } else {
                break
              }
            }
          }

          cluster_size <- best_end - i + 1L
          if (best_hom >= hom_thresh && cluster_size >= min_size) {
            cluster_counter <- cluster_counter + 1L
            cluster_indices <- seg[i:best_end]
            clusters[[cluster_counter]] <- create_cluster(cluster_indices, seg_classes[i:best_end], cluster_counter)
            assigned[cluster_indices] <- TRUE
            log_msg(sprintf("    [OK] Created cluster %d: %d elements, hom=%.3f",
                            cluster_counter, cluster_size, best_hom))
          }
          i <- best_end + 1L
        }
      } else {
        log_msg(sprintf("    -> Segment too small (%d < min_size=%d), will try Phase 4", n_seg, min_size))
      }
    }
  }

  # === PHASE 4: Final merge attempt ===

  remaining <- which(!assigned)
  if (length(remaining) > 0) {
    log_msg(sprintf("\nPhase 4: Final merge attempt for %d elements", length(remaining)))
    segments <- find_segments(remaining)

    for (seg in segments) {
      if (all(assigned[seg])) next
      seg_classes <- class_labels[seg]

      log_msg(sprintf("  Attempting to merge %d elements at positions %d-%d",
                      length(seg), min(seg), max(seg)))

      best_cid <- NULL
      best_hom <- -1

      for (cid in seq_along(clusters)) {
        if (is.null(clusters[[cid]])) next
        if (max(clusters[[cid]]$indices) == seg[1] - 1L) {
          temp_hom <- calc_homogeneity(c(clusters[[cid]]$classes, seg_classes))
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
        log_msg(sprintf("    [OK] Merged into cluster %d (hom=%.3f)", clusters[[best_cid]]$id, best_hom))
      } else {
        log_msg(if (is.null(best_cid)) "    [FAIL] No adjacent cluster found" else
          sprintf("    [FAIL] Best merge would have hom=%.3f (< %.3f)", best_hom, hom_thresh))
      }
    }
  }

  # === PHASE 5: Force cluster creation for unassigned elements ===

  remaining <- which(!assigned)
  if (length(remaining) > 0) {
    log_msg(sprintf("\nPhase 5: Creating clusters for %d unassigned elements by class", length(remaining)))

    for (cls in unique(class_labels[remaining])) {
      cls_indices <- which(class_labels == cls & !assigned)
      if (length(cls_indices) == 0) next

      segments <- find_segments(cls_indices)
      log_msg(sprintf("  Class '%s': %d unassigned elements in %d segment(s)",
                      cls, length(cls_indices), length(segments)))

      for (seg in segments) {
        cluster_counter <- cluster_counter + 1L
        clusters[[cluster_counter]] <- create_cluster(seg, class_labels[seg], cluster_counter)
        assigned[seg] <- TRUE

        log_msg(sprintf("    Created %s cluster %d at position%s %s",
                        if (length(seg) == 1) "singleton" else "",
                        cluster_counter,
                        if (length(seg) > 1) "s" else "",
                        if (length(seg) == 1) seg else sprintf("%d-%d", min(seg), max(seg))))
      }
    }
    log_msg("  Phase 5 complete: All elements now assigned\n")
  }

  # === PHASE 6: Intelligent cluster merging ===

  if (length(clusters) > 1) {
    log_msg("\nPhase 6: Intelligent cluster merging")
    log_msg("  Rule: Protect complete-class clusters from contamination")

    merged_any <- TRUE
    iteration <- 0
    max_iterations <- 10

    while (merged_any && iteration < max_iterations) {
      merged_any <- FALSE
      iteration <- iteration + 1
      log_msg(sprintf("\n  Iteration %d:", iteration))

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

        is_i_complete <- is_complete_class_cluster(cluster_i$classes, class_total_counts)
        is_j_complete <- is_complete_class_cluster(cluster_j$classes, class_total_counts)

        merged_classes <- c(cluster_i$classes, cluster_j$classes)
        merged_hom <- calc_homogeneity(merged_classes)

        can_merge <- FALSE
        reason <- ""

        if (is_i_complete && is_j_complete) {
          dom_i <- get_dominant(cluster_i$classes)
          dom_j <- get_dominant(cluster_j$classes)
          can_merge <- (dom_i == dom_j)
          reason <- if (can_merge) sprintf("Both complete for same class '%s'", dom_i) else
            sprintf("Both complete: '%s' (%d) and '%s' (%d) - cannot contaminate",
                    dom_i, length(cluster_i$classes), dom_j, length(cluster_j$classes))

        } else if (is_i_complete) {
          dom_i <- get_dominant(cluster_i$classes)
          can_merge <- all(cluster_j$classes == dom_i)
          reason <- if (can_merge) sprintf("Cluster %d complete for '%s', cluster %d also only '%s'",
                                           cluster_i$id, dom_i, cluster_j$id, dom_i) else
                                             sprintf("Cluster %d is COMPLETE for '%s' (%d elements) - cannot add other classes",
                                                     cluster_i$id, dom_i, length(cluster_i$classes))

        } else if (is_j_complete) {
          dom_j <- get_dominant(cluster_j$classes)
          can_merge <- all(cluster_i$classes == dom_j)
          reason <- if (can_merge) sprintf("Cluster %d complete for '%s', cluster %d also only '%s'",
                                           cluster_j$id, dom_j, cluster_i$id, dom_j) else
                                             sprintf("Cluster %d is COMPLETE for '%s' (%d elements) - cannot add other classes",
                                                     cluster_j$id, dom_j, length(cluster_j$classes))

        } else {
          can_merge <- (merged_hom >= hom_thresh)
          reason <- if (can_merge) {
            sprintf("Non-complete clusters: hom_i=%.3f, hom_j=%.3f -> merged=%.3f >= %.3f",
                    calc_homogeneity(cluster_i$classes), calc_homogeneity(cluster_j$classes),
                    merged_hom, hom_thresh)
          } else {
            sprintf("Merged homogeneity %.3f < threshold %.3f", merged_hom, hom_thresh)
          }
        }

        if (can_merge) {
          old_id_i <- cluster_i$id
          old_id_j <- cluster_j$id
          clusters[[i]]$indices <- c(cluster_i$indices, cluster_j$indices)
          clusters[[i]]$classes <- merged_classes
          clusters[[i+1]] <- NULL
          merged_any <- TRUE

          log_msg(sprintf("    [OK] Merged clusters %d and %d (n=%d + %d = %d, hom=%.3f)",
                          old_id_i, old_id_j, length(cluster_i$classes),
                          length(cluster_j$classes), length(merged_classes), merged_hom))
          if (reason != "") log_msg(sprintf("      Reason: %s", reason))
          next
        } else {
          log_msg(sprintf("    [FAIL] Cannot merge clusters %d and %d", cluster_i$id, cluster_j$id))
          if (reason != "") log_msg(sprintf("      Reason: %s", reason))
        }
        i <- i + 1
      }

      clusters <- clusters[!sapply(clusters, is.null)]
      if (!merged_any) log_msg("    No more merges possible")
    }

    log_msg(sprintf("\n  Phase 6 complete after %d iteration(s)", iteration))
    log_msg(sprintf("  Final cluster count: %d", length(clusters)))
  }

  # === FINALIZATION ===

  clusters <- clusters[!sapply(clusters, is.null)]

  if (length(clusters) > 0) {
    order_idx <- order(sapply(clusters, function(cl) min(cl$indices)))
    clusters <- clusters[order_idx]
    for (k in seq_along(clusters)) clusters[[k]]$id <- k
    log_msg("\nClusters reordered by dendrogram position")
  }

  n_clusters <- length(clusters)
  n_unassigned <- sum(!assigned)

  log_msg(sprintf("\n=== CLUSTERING COMPLETE ==="))
  log_msg(sprintf("Total clusters: %d", n_clusters))
  log_msg(sprintf("Unassigned elements: %d", n_unassigned))

  if (n_unassigned > 0) {
    warning("Some elements remain unassigned - this should not happen!", call. = FALSE)
  }

  # === BUILD OUTPUTS ===

  cluster_assignment <- rep(NA_integer_, n_elements)
  for (cid in seq_along(clusters)) {
    cluster_assignment[clusters[[cid]]$indices] <- clusters[[cid]]$id
  }

  cluster_orig <- if (!is.null(dendro_order)) {
    cluster_orig <- rep(NA_integer_, n_elements)
    for (i in seq_len(n_elements)) cluster_orig[dendro_order[i]] <- cluster_assignment[i]
    cluster_orig
  } else {
    cluster_assignment
  }

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
      cluster_summary$homogeneity[k] <- calc_homogeneity(cl$classes)
      cluster_summary$n_classes[k] <- length(class_table)
      cluster_summary$class_composition[k] <- paste(names(class_table), class_table, sep = ":", collapse = "; ")
      cluster_summary$is_complete_class[k] <- is_complete
    }
  }

  if (is.null(sequence_names)) {
    sequence_names <- sprintf("seq_%d", seq_len(n_elements))
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
      cluster_assignment_original_order = cluster_orig,
      unassigned_elements = which(!assigned),
      n_unassigned = n_unassigned,
      hom_thresh = hom_thresh,
      min_size = min_size,
      n_elements = n_elements,
      valid_clusters = clusters
    ),
    class = "cluster_dendrogram_result"
  )
}

# === S3 METHODS ===

#' @export
print.cluster_dendrogram_result <- function(x, ...) {
  cat("\nDendrogram Clustering Result\n")
  cat("=============================\n\n")
  cat(sprintf("Total elements:     %d\n", x$n_elements))
  cat(sprintf("Total clusters:     %d\n", nrow(x$cluster_summary)))
  cat(sprintf("Unassigned:         %d\n", x$n_unassigned))
  cat(sprintf("Min size:           %d\n", x$min_size))
  cat(sprintf("Hom threshold:      %.3f\n", x$hom_thresh))

  if (nrow(x$cluster_summary) > 0) {
    cat("\nCluster Summary:\n")
    print(x$cluster_summary, row.names = FALSE)

    complete_clusters <- x$cluster_summary[x$cluster_summary$is_complete_class, ]
    if (nrow(complete_clusters) > 0) {
      cat("\n[OK] Complete-class clusters (protected from contamination):\n")
      for (i in 1:nrow(complete_clusters)) {
        cat(sprintf("  Cluster %d: ALL %d elements of class '%s'\n",
                    complete_clusters$cluster_id[i],
                    complete_clusters$n_elements[i],
                    complete_clusters$dominant_class[i]))
      }
    }
  }

  if (x$n_unassigned > 0) {
    cat("\n[WARNING] Some elements could not be assigned to clusters\n")
  }

  invisible(x)
}

#' @export
summary.cluster_dendrogram_result <- function(object, ...) {
  cat("\nClustering Summary\n")
  cat("------------------\n")
  cat(sprintf("Elements:  %d\n", object$n_elements))
  cat(sprintf("Clusters:  %d\n", nrow(object$cluster_summary)))
  cat(sprintf("Unassigned: %d\n", object$n_unassigned))

  if (nrow(object$cluster_summary) > 0) {
    cat(sprintf("Mean homogeneity: %.3f\n", mean(object$cluster_summary$homogeneity)))
    cat(sprintf("Mean size: %.1f\n", mean(object$cluster_summary$n_elements)))
    cat(sprintf("Size range: %d - %d\n",
                min(object$cluster_summary$n_elements),
                max(object$cluster_summary$n_elements)))
    cat(sprintf("Complete-class clusters: %d\n", sum(object$cluster_summary$is_complete_class)))
  }

  invisible(object)
}
