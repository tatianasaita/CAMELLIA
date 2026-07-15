#' Cluster Dendrogram Leaves with Homogeneity and Size Constraints
#'
#' Creates clusters from dendrogram leaves following their hierarchical order,
#' applying homogeneity and size constraints while respecting dendrogram contiguity.
#'
library(apcluster)
#' @param dendrogram A dendrogram object
#' @param class_labels Character vector of class labels in dendrogram order
#' @param hom_thresh Minimum homogeneity threshold (0-1)
#' @param min_size_cluster Minimum cluster size (default: 3)
#' @param verbose Print progress messages (default: TRUE)
#' @param sequence_names Character vector of sequence names
#' @param data_result List with 'kmers' and 'metadata' for adding assignments
#' @param method Clustering method: "dendrogram" or "apcluster" (default: "dendrogram")
#' @param feature_matrix Numeric matrix for apcluster (required if method = "apcluster")
#' @param ap_r Distance metric parameter for apcluster (default: 2)
#'
#' @return List with clustering results:
#' \itemize{
#'   \item \code{method}: Clustering method used
#'   \item \code{clusters}: List of cluster objects
#'   \item \code{cluster_summary}: Data frame with cluster statistics
#'   \item \code{element_assignment}: Data frame with element-to-cluster assignments
#'   \item \code{min_size_input}: Input parameter value for minimum cluster size
#'   \item \code{min_size_observed}: Actual minimum cluster size in results
#'   \item \code{hom_thresh_input}: Input parameter value (dendrogram only)
#'   \item \code{hom_thresh_observed}: Actual minimum homogeneity observed
#' }
#'
#' @details
#' Dendrogram method (default):
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
#' Affinity Propagation method:
#' \itemize{
#'   \item Computes negative distance similarity matrix
#'   \item Identifies exemplar sequences through message passing
#'   \item Assigns all sequences to their nearest exemplar
#'   \item Computes homogeneity and class composition metrics
#'   \item Uses default parameters: damping=0.9, maxits=1000, convits=100
#' @return List with clustering results:
#' \itemize{
#'   \item \code{method}: Clustering method used
#'   \item \code{clusters}: List of cluster objects
#'   \item \code{cluster_summary}: Data frame with cluster statistics
#'   \item \code{element_assignment}: Data frame with element-to-cluster assignments
#'   \item \code{min_size_input}: Input parameter value for minimum cluster size
#'   \item \code{min_size_observed}: Actual minimum cluster size in results
#'   \item \code{hom_thresh_input}: Input parameter value (dendrogram only)
#'   \item \code{hom_thresh_observed}: Actual minimum homogeneity observed
#' }
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
#' # Dendrogram method
#' result_cluster_dendrogram <- cluster_dendrogram(
#'   dendrogram = dend,
#'   class_labels = classes,
#'   hom_thresh = 0.7,
#'   min_size_cluster = 3,
#'   method = "dendrogram"
#' )
#' # Affinity Propagation method
#' result_cluster_ap <- cluster_dendrogram(
#'   class_labels = classes,
#'   method = "apcluster",
#'   feature_matrix = kmers_matrix,
#'   ap_r = 2
#' )
#'}
#'
#' @importFrom stats is.leaf
#'
#' @export
cluster_dendrogram <- function(dendrogram = NULL,
class_labels,
hom_thresh = NULL,
min_size_cluster = 3L,
verbose = TRUE,
sequence_names = NULL,
data_result = NULL,
method = c("dendrogram", "apcluster"),
feature_matrix = NULL,
ap_r = 2,
seq_lengths    = NULL,
min_seq_length = 800L) {

  method <- match.arg(method)
  min_size_input   <- as.integer(min_size_cluster)
  hom_thresh_input <- hom_thresh

  # Affinity propagation method
  if (method == "apcluster") {

    if (!requireNamespace("apcluster", quietly = TRUE)) {
      stop("Package 'apcluster' is required. Install with: install.packages('apcluster')")
    }

    if (is.null(feature_matrix)) {
      stop("'feature_matrix' must be provided when method = 'apcluster'.")
    }

    n_elements_full <- nrow(feature_matrix)
    n_elements_original <- n_elements_full

    class_labels_orig <- class_labels

    if (is.null(sequence_names)) {
      sequence_names <- paste0("seq_", seq_len(n_elements_full))
    }
    sequence_names_orig <- sequence_names

    # Aplicar filtro de comprimento se fornecido

    pass_length <- seq_len(n_elements_full)

    if (!is.null(seq_lengths)) {
      if (length(seq_lengths) != n_elements_full) {
        stop("'seq_lengths' must have the same length as nrow(feature_matrix).")
      }

      pass_length <- which(seq_lengths >= min_seq_length)
      if (length(pass_length) == 0L) {
        stop(sprintf("No sequences passed the minimum length filter (min_seq_length = %d).", min_seq_length))
      }

      if (verbose) {
        message(sprintf("[Sampling] Step 1 — Length filter (>= %d bp): %d / %d sequences retained.",
                        min_seq_length, length(pass_length), n_elements_full))
      }

      # Aplicar filtro
      feature_matrix <- feature_matrix[pass_length, , drop = FALSE]
      class_labels <- class_labels[pass_length]
      sequence_names <- sequence_names[pass_length]
 #     n_elements_full <- length(pass_length)
    }

    # ETAPA 1: Aplicar apcluster por classe para obter exemplares
    unique_classes <- unique(class_labels)
    n_classes <- length(unique_classes)

    if (verbose) {
      message(sprintf("[Sampling] Step 2 — Applying apcluster within each class (%d classes)", n_classes))
    }

    # Coletar exemplares de cada classe
    selected_indices_all <- integer(0)

    for (cls in unique_classes) {
      cls_indices <- which(class_labels == cls)
      cls_matrix <- feature_matrix[cls_indices, , drop = FALSE]
      n_cls <- nrow(cls_matrix)

      if (verbose) {
        message(sprintf("[Sampling]   Processing class '%s': %d sequences", cls, n_cls))
      }

      # Aplicar apcluster na classe
      sim_matrix_cls <- apcluster::negDistMat(cls_matrix, r = ap_r)

      ap_result_cls <- apcluster::apcluster(
        s = sim_matrix_cls,
        lam = 0.9,
        maxits = 500,
        convits = 50,
        details = FALSE
      )

      n_exemplars <- length(ap_result_cls@exemplars)

      if (verbose) {
        message(sprintf("[Sampling]     Class '%s': selected %d exemplars", cls, n_exemplars))
      }

      # Converter índices locais para globais
      exemplar_local_indices <- cls_indices[ap_result_cls@exemplars]
      selected_indices_all   <- c(selected_indices_all, exemplar_local_indices)
    }

    selected_indices_all <- sort(unique(selected_indices_all))

    if (verbose) {
      message(sprintf("[Sampling] Total exemplars selected: %d sequences", length(selected_indices_all)))
    }

    selected_indices_all_global <- pass_length[selected_indices_all]

    # ETAPA 2: Aplicar apcluster final em todos os exemplares juntos
    if (verbose) {
      message("[Clustering] Step 3 — Final apcluster on all selected exemplars")
    }

    # Subset final para apcluster
    feature_matrix_final <- feature_matrix[selected_indices_all, , drop = FALSE]
    class_labels_final <- class_labels[selected_indices_all]
    seq_names_final <- sequence_names[selected_indices_all]

    if (is.null(rownames(feature_matrix_final))) {
      rownames(feature_matrix_final) <- seq_names_final
    }

    class_total_counts <- table(class_labels)

    # Apcluster final
    sim_matrix_final <- apcluster::negDistMat(feature_matrix_final, r = ap_r)

    if (verbose) message("[Clustering] Computing final similarity matrix and clustering...")

    ap_result_final <- apcluster::apcluster(
      s = sim_matrix_final,
      lam = 0.9,
      maxits = 1000,
      convits = 100,
      details = verbose
    )

    if (verbose) {
      message(sprintf("[Clustering] Final clustering: %d clusters formed", length(ap_result_final@exemplars)))
    }

    # Extrair atribuições de cluster
    n_clusters_ap <- length(ap_result_final@clusters)
    clusters_ap <- vector("list", n_clusters_ap)

    for (k in seq_len(n_clusters_ap)) {
      local_idx <- ap_result_final@clusters[[k]]
      global_idx <- selected_indices_all_global[local_idx]

      clusters_ap[[k]] <- list(
        id = k,
        indices = local_idx,
        global_indices = global_idx,
        classes = class_labels_final[local_idx],
        exemplar_id = ap_result_final@exemplars[k]
      )
    }

    # Atribuições
    n_elements_clustered <- length(selected_indices_all)
    n_filtered           <- nrow(feature_matrix)  # tamanho filtrado (6133)

    # Atribuição no espaço dos exemplares (1..611)
    cluster_assignment_local <- rep(NA_integer_, n_elements_clustered)
    for (cid in seq_along(clusters_ap)) {
      cluster_assignment_local[clusters_ap[[cid]]$indices] <- clusters_ap[[cid]]$id
    }

    # Expandir: exemplares (1..611) → subconjunto filtrado (1..6133)
    cluster_assignment_filtered <- rep(NA_integer_, n_filtered)
    cluster_assignment_filtered[selected_indices_all] <- cluster_assignment_local

    # Expandir: subconjunto filtrado (1..6133) → dataset original (1..6964)
    cluster_assignment_original <- rep(NA_integer_, n_elements_original)
    cluster_assignment_original[pass_length] <- cluster_assignment_filtered


    # Estatísticas dos clusters
    cluster_sizes_ap <- sapply(clusters_ap, function(cl) length(cl$indices))
    min_size_observed_ap <- min(cluster_sizes_ap)

    cluster_summary_ap <- .build_cluster_summary_ap(
      clusters = clusters_ap,
      class_total_counts = class_total_counts,
      n_clusters = n_clusters_ap
    )

    hom_thresh_observed_ap <- min(cluster_summary_ap$homogeneity)

    # Construir vetor de índices globais dos exemplares
    exemplar_global <- selected_indices_all_global[ap_result_final@exemplars]
    
    element_assignment_ap <- data.frame(
      sequence_name = sequence_names_orig,
      index         = seq_len(n_elements_original),
      class         = class_labels_orig,
      cluster       = cluster_assignment_original,
      is_selected   = seq_len(n_elements_original) %in% selected_indices_all_global,
      is_exemplar   = seq_len(n_elements_original) %in% exemplar_global,
      stringsAsFactors = FALSE
    )

    # data_result
    if (!is.null(data_result)) {
      data_result$metadata$cluster <- cluster_assignment_original  # FIX 4
    }

    return(structure(
      list(
        method          = "apcluster",
        ap_result       = ap_result_final,
        dendrogram      = NULL,
        clusters        = clusters_ap,
        cluster_summary = cluster_summary_ap,
        element_assignment = element_assignment_ap,
        data_result     = data_result,
        cluster_assignment_dendro_order = NULL,
        cluster_assignment = cluster_assignment_original,          # FIX 4: tamanho original (6964)
        unassigned_elements = integer(0),
        n_unassigned    = 0L,
        hom_thresh_input    = hom_thresh_input,
        hom_thresh_observed = hom_thresh_observed_ap,
        min_size_input      = min_size_input,
        min_size_observed   = min_size_observed_ap,
        n_elements          = n_elements_original,                 # FIX 1: tamanho original (6964)
        n_elements_clustered = n_elements_clustered,
        valid_clusters      = clusters_ap,
        selected_indices    = selected_indices_all_global,         # FIX 3: índices globais
        ap_parameters = list(
          r          = ap_r,
          damping    = 0.9,
          maxiter    = 1000,
          conviter   = 100,
          preference = NULL
        )
      ),
      class = "cluster_dendrogram_result"
    ))
  }

  # Dendrogram Method
  n_elements <- length(class_labels)
  clusters_dend <- list()
  cluster_counter_dend <- 0L
  assigned_dend <- rep(FALSE, n_elements)
  class_total_counts <- table(class_labels)

  if (is.null(sequence_names)) {
    sequence_names <- paste0("seq_", seq_len(n_elements))
  }

  # Handle rare classes (total count < min_size_cluster): Homogeneous cluster with classes < min_size_cluster
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

  # Create pure class segments (homogeneity = 1.0): Homogeneous clusters with classes > min_size_cluster
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

  # Greedy expansion with homogeneity constraint: Extending last cluster or forming new cluster
  i <- 1L
  while (i <= n_elements) {
    if (assigned_dend[i]) {
      i <- i + 1L
      next
    }

    # Try extending last cluster
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

    # Try forming new cluster
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

      # Try merging with adjacent cluster
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

      # Try forming clusters within segment
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

  # Cluster merging: Clusters with 1.0 homogeneity are not merged
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
          # old_id_i_dend <- cluster_i_dend$id
          # old_id_j_dend <- cluster_j_dend$id
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

  # Outputs
  cluster_assignment_dend <- rep(NA_integer_, n_elements)
  for (cid in seq_along(clusters_dend)) {
    cluster_assignment_dend[clusters_dend[[cid]]$indices] <- clusters_dend[[cid]]$id
  }

  cluster_orig_dend <- cluster_assignment_dend

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

  # Observed minimum homogeneity
  hom_thresh_observed_dend <- if (n_clusters_dend > 0) min(cluster_summary_dend$homogeneity) else NA_real_

  # Observed minimum size
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
      method = "dendrogram",
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
