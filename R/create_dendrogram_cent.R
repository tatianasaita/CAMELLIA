#' Create and plot a dendrogram from the given dataset
#'
#' @param data Data frame containing k-mer counts and a 'CLASS' column.
#' @param sequence_names Character vector of sequence names (optional).
#' @param dist_method Distance method (default: "euclidean"). See \code{\link[parallelDist]{parDist}}.
#' @param hclust_method Hierarchical clustering method (default: "ward.D2"). See \code{\link[stats]{hclust}}.
#' @param output Output file path for PNG. If NULL, plots to screen.
#' @param seq_lengths Integer vector with sequence lengths aligned with rows of \code{data}.
#'   When provided, activates length filter and centroid-distance sampling. Default: NULL.
#' @param min_seq_length Minimum sequence length to retain (default: 800). Only used
#'   when \code{seq_lengths} is provided.
#'
#' @return A list (invisibly) containing: dendrogram, hclust, order, labels,
#'   sequence_names, colors, base_colors, selected_indices, n_per_class_used,
#'   limiting_class.
#'
#' @details
#' \enumerate{
#'   \item (Optional) Filter sequences by \code{min_seq_length}.
#'   \item (Optional) Determine \code{n_per_class} from the limiting class and
#'     sample 70\% from the far tail and 30\% from the near tail of each class's
#'     centroid-distance distribution.
#'   \item Distance calculation based on k-mer counts.
#'   \item Hierarchical clustering.
#'   \item Converts hclust object to dendrogram and orders samples.
#'   \item Applies colors to dendrogram leaves.
#' }
#'
#' @importFrom parallel detectCores
#' @importFrom parallelDist parDist
#' @importFrom fastcluster hclust
#' @importFrom stats as.dendrogram order.dendrogram
#' @importFrom grDevices palette.colors rainbow png dev.off
#' @importFrom dendextend set
#' @export
create_dendrogram_cent <- function(data,
                              sequence_names = NULL,
                              dist_method    = "euclidean",
                              hclust_method  = "ward.D2",
                              output         = NULL,
                              seq_lengths    = NULL,
                              min_seq_length = 800L) {

  # Proportions for centroid-distance sampling (fixed)
  prop_far  <- 0.70   # 70% from the far tail
  prop_near <- 0.30   # 30% from the near tail
  max_total_seqs <- 5000L
  # ---------------------------------------------------------------------------
  # Extract data
  # ---------------------------------------------------------------------------
  class_labels <- data$CLASS
  data_matrix  <- as.matrix(data[, !colnames(data) %in% "CLASS", drop = FALSE])
  n_samples    <- nrow(data_matrix)
  n_classes_total        <- length(unique(class_labels))

  selected_indices <- NULL
  n_per_class_used <- NULL
  limiting_class   <- NULL

  # Centroid-distance sampling (only when seq_lengths is provided)
  if (!is.null(seq_lengths)) {

    if (length(seq_lengths) != n_samples) {
      stop("'seq_lengths' must have the same length as the number of rows in 'data'.")
    }

    # Step 1: filter by minimum sequence length
    pass_length <- which(seq_lengths >= min_seq_length)

    if (length(pass_length) == 0L) {
      stop(sprintf(
        "No sequences passed the minimum length filter (min_seq_length = %d).",
        min_seq_length
      ))
    }

    message(sprintf(
      "[Sampling] Step 1 — Length filter (>= %d bp): %d / %d sequences retained.",
      min_seq_length, length(pass_length), n_samples
    ))

    # Step 2: find the limiting class
    filtered_labels     <- class_labels[pass_length]
    counts_after_filter <- table(filtered_labels)
    n_classes_filtered <- length(counts_after_filter)

    message("[Sampling] Step 2 — Sequences per class after length filter:")
    for (cls_name in names(counts_after_filter)) {
      message(sprintf("             %s : %d", cls_name, counts_after_filter[[cls_name]]))
    }

    n_per_class_used <- as.integer(min(counts_after_filter))
    limiting_class   <- names(which.min(counts_after_filter))

    message(sprintf(
      "[Sampling] Limiting class: '%s' with %d sequences → n_per_class = %d for all classes.",
      limiting_class, n_per_class_used, n_per_class_used
    ))

    # Cap total sequences at 5000

    if (n_per_class_used * n_classes_filtered > max_total_seqs) {
      n_per_class_capped <- floor(max_total_seqs / n_classes_filtered)

      message(sprintf(
        "[Sampling] Cap applied: %d classes x %d = %d > 5000 (max).",
        n_classes_filtered, n_per_class_used, n_per_class_used * n_classes_filtered
      ))
      message(sprintf(
        "[Sampling] n_per_class reduced from %d to %d (= floor(5000 / %d)).",
        n_per_class_used, n_per_class_capped, n_classes_filtered
      ))

      n_per_class_used <- as.integer(n_per_class_capped)
      limiting_class   <- paste0(limiting_class, " [capped at 5000]")
    }

    message(sprintf(
      "[Sampling] Final n_per_class = %d | Expected total <= %d sequences.",
      n_per_class_used, n_per_class_used * n_classes_filtered
    ))

    # Step 3: centroid-distance sampling (70% far / 30% near)
    n_far_per_class  <- round(n_per_class_used * prop_far)
    n_near_per_class <- n_per_class_used - n_far_per_class

    message(sprintf(
      "[Sampling] Step 3 — Centroid-distance sampling: %d far (70%%) + %d near (30%%) per class.",
      n_far_per_class, n_near_per_class
    ))

    selected_indices <- integer(0)

    for (cls in names(counts_after_filter)) {

      cls_pass_idx <- pass_length[filtered_labels == cls]
      n_cls        <- length(cls_pass_idx)

      if (n_cls <= n_per_class_used) {
        selected_indices <- c(selected_indices, cls_pass_idx)
        message(sprintf(
          "[Sampling]   Class '%s': only %d sequences available — keeping all.",
          cls, n_cls
        ))
        next
      }

      cls_mat  <- data_matrix[cls_pass_idx, , drop = FALSE]
      centroid <- colMeans(cls_mat)
      dists    <- sqrt(rowSums(sweep(cls_mat, 2L, centroid)^2))
      ord      <- order(dists)  # ascending: near -> far

      near_local <- ord[seq_len(n_near_per_class)]
      far_local  <- ord[(n_cls - n_far_per_class + 1L):n_cls]

      chosen_global    <- cls_pass_idx[unique(c(near_local, far_local))]
      selected_indices <- c(selected_indices, chosen_global)

      message(sprintf(
        "[Sampling]   Class '%s': selected %d sequences (%d near + %d far).",
        cls, length(chosen_global), n_near_per_class, n_far_per_class
      ))
    }

    selected_indices <- sort(unique(selected_indices))

    message(sprintf(
      "[Sampling] Total sequences selected for dendrogram: %d (hard cap: 5000)",
      length(selected_indices)
    ))

    # Subset to selected sequences
    data_matrix  <- data_matrix[selected_indices, , drop = FALSE]
    class_labels <- class_labels[selected_indices]

    if (!is.null(sequence_names)) {
      sequence_names <- sequence_names[selected_indices]
    }
  }

  # Distance calculation
  n_threads <- max(1L, parallel::detectCores() - 1L)

  dist_matrix <- parallelDist::parDist(
    x       = data_matrix,
    method  = dist_method,
    threads = n_threads
  )

  rm(data_matrix)
  invisible(gc(verbose = FALSE))

  # Hierarchical clustering
  hc_result <- fastcluster::hclust(dist_matrix, method = hclust_method)

  rm(dist_matrix)
  invisible(gc(verbose = FALSE))


  # Build dendrogram
  dend       <- stats::as.dendrogram(hc_result)
  dend_order <- stats::order.dendrogram(dend)

  ordered_labels <- as.character(class_labels[dend_order])

  sequence_names_ordered <- NULL
  if (!is.null(sequence_names)) {
    sequence_names_ordered <- sequence_names[dend_order]
  }

  # Color mapping
  unique_classes <- unique(class_labels)
  n_classes_plot      <- length(unique_classes)

  if (n_classes_plot <= 12L) {
    base_colors <- grDevices::palette.colors(n = n_classes_plot, palette = "Paired")
  } else {
    base_colors <- grDevices::rainbow(n_classes_plot)
  }

  names(base_colors) <- unique_classes
  label_colors       <- base_colors[as.character(class_labels[dend_order])]

  dend <- dendextend::set(dend, "labels_col", label_colors)

  # Plot
  .plot_dendrogram <- function() {
    plot(dend, main = "Dendrogram", ylab = "Height")
    legend("topright",
           legend  = names(base_colors),
           col     = base_colors,
           pch     = 15,
           pt.cex  = 2,
           cex     = 0.8,
           title   = "Classes",
           bg      = "white",
           box.lty = 1)
  }

  if (is.null(output)) {
    .plot_dendrogram()
  } else {
    n_leaves  <- length(class_labels)
    img_width <- max(800L, n_leaves * 2L)   # minimum 800px, ~2px per leaf

    grDevices::png(filename = output, width = img_width, height = 600, res = 150)
    .plot_dendrogram()
    grDevices::dev.off()
  }

  # Build return labels / sequence_names vectors

  labels_vector <- rownames(data)
  if (is.null(labels_vector)) {
    labels_vector <- paste0("S", seq_len(nrow(data)))
  }

  if (!is.null(selected_indices)) {
    labels_vector <- labels_vector[selected_indices]
  }

  labels_vector_ordered <- labels_vector[dend_order]

  sequence_names_result <- if (!is.null(sequence_names_ordered)) {
    sequence_names_ordered
  } else {
    labels_vector_ordered
  }

  # Return
  result <- list(
    dendrogram       = dend,
    hclust           = hc_result,
    order            = hc_result$order,
    labels           = ordered_labels,          # class labels in dendrogram order
    sequence_names   = sequence_names_result,   # sequence names in dendrogram order
    colors           = label_colors,            # per-leaf color in dendrogram order
    base_colors      = base_colors,             # palette: class name -> color
    selected_indices = selected_indices,        # global indices selected (NULL if no sampling)
    n_per_class_used = n_per_class_used,        # samples per class (NULL if no sampling)
    limiting_class   = limiting_class           # class that defined n (NULL if no sampling)
  )

  class(result) <- c("dendrogram_result", "list")

  invisible(result)
}
