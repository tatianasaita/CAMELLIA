#' Select Motifs from Cluster Analysis
#'
#' Selects n motifs based on cluster information, distributing them intelligently
#' across classes according to their distribution across all clusters.
#'
#' @param motif_cluster Data.frame with k-mer motifs by cluster (normalized).
#' @param cluster_result List containing cluster_summary from cluster_dendrogram().
#' @param n Integer. Number of motifs to select.
#' @param verbose Logical. If TRUE, prints detailed information. Default is TRUE.
#'
#' @return An object of class "select_motifs" with selected motifs per class.
#'
#' @export
select_motifs <- function(motif_cluster, cluster_result, n, verbose = TRUE) {

  # ===== INPUT VALIDATION =====

  if (!is.numeric(n) || n <= 0 || n != as.integer(n)) {
    stop("Parameter 'n' must be a positive integer")
  }

  if (!is.data.frame(motif_cluster)) {
    stop("'motif_cluster' must be a data.frame")
  }

  # Validate that rownames exist and are meaningful (not default numeric)
  if (is.null(rownames(motif_cluster)) || length(rownames(motif_cluster)) == 0) {
    stop("'motif_cluster' must have rownames containing k-mer sequences")
  }

  # Check if rownames are only numeric (R default)
  rn <- rownames(motif_cluster)
  if (all(grepl("^[0-9]+$", rn))) {
    stop("'motif_cluster' must have meaningful rownames containing k-mer sequences (not default numeric rownames)")
  }

  if (!all(grepl("^Cluster_", colnames(motif_cluster)))) {
    stop("'motif_cluster' columns must be named 'Cluster_<id>'")
  }

  if (!is.list(cluster_result)) {
    stop("'cluster_result' must be a list from cluster_dendrogram()")
  }

  if (!("cluster_summary" %in% names(cluster_result))) {
    stop("'cluster_result' must contain 'cluster_summary' element")
  }

  cluster_summary <- cluster_result$cluster_summary

  if (!is.data.frame(cluster_summary)) {
    stop("'cluster_result$cluster_summary' must be a data.frame")
  }

  required_cols <- c("cluster_id", "dominant_class", "class_distribution")
  missing_cols <- setdiff(required_cols, colnames(cluster_summary))
  if (length(missing_cols) > 0) {
    stop(paste("'cluster_summary' missing columns:", paste(missing_cols, collapse = ", ")))
  }

  if (!is.logical(verbose)) {
    stop("'verbose' must be logical (TRUE or FALSE)")
  }

  # ===== EXTRACT BASIC INFORMATION =====

  classes <- unique(cluster_summary$dominant_class)
  k <- length(classes)

  clusters <- unique(cluster_summary$cluster_id)
  m <- as.numeric(length(clusters))

  classe_dominant_count <- table(cluster_summary$dominant_class)

  if (verbose) {
    message(sprintf("Number of classes (k): %d", k))
    message(sprintf("Total number of clusters (m): %d", m))
  }

  # ===== COUNT CLASS DISTRIBUTION ACROSS ALL CLUSTERS =====

  classes_separated <- .extract_all_classes(cluster_summary$class_distribution)
  classe_distribution_count <- table(classes_separated)
  classe_distribution_count <- sort(classe_distribution_count, decreasing = TRUE)
  classe_order <- names(classe_distribution_count)

  if (verbose) {
    message("Classes ordered by distribution across all clusters (decreasing):")
    for (i in seq_along(classe_order)) {
      classe <- classe_order[i]
      count <- classe_distribution_count[classe]
      message(sprintf("  %d: %s - %d clusters", i, classe, count))
    }
  }

  # ===== INITIALIZE RESULT STORAGE =====

  selected_motifs <- list()
  case_used <- NULL

  # ===== CASE 1: n < k =====

  if (n < k) {
    if (verbose) {
      message(sprintf(
        "\nCASE 1: n < k\nImpossible to select at least one motif per class (%d < %d)",
        n, k
      ))
    }

    result <- structure(
      list(),
      class = c("select_motifs", "empty_result"),
      n = n,
      k = k,
      m = m,
      case = "CASE_1",
      reason = "n < k: insufficient motifs for all classes"
    )

    return(invisible(result))
  }

  # ===== CASE 2: k <= n < m =====

  else if (k <= n && n < m) {
    if (verbose) {
      message(sprintf(
        "\nCASE 2: k <= n < m\nDistributing motifs by class (%d <= %d < %d)",
        k, n, m
      ))
    }

    case_used <- "CASE_2"

    motifs_per_classe_base <- floor(n / k)
    remainder <- n %% k

    motifs_per_classe <- setNames(
      rep(motifs_per_classe_base, length(classes)),
      classes
    )

    if (remainder > 0) {
      for (i in seq_len(remainder)) {
        classe_to_add <- classe_order[((i - 1) %% k) + 1]
        motifs_per_classe[classe_to_add] <- motifs_per_classe[classe_to_add] + 1
      }
    }

    if (verbose) {
      message("\nMotif distribution per class:")
      for (classe in names(motifs_per_classe)) {
        message(sprintf("  %s: %d motifs", classe, motifs_per_classe[classe]))
      }
    }

    all_selected_motifs <- new.env(hash = TRUE)

    for (classe in classes) {
      num_motifs_to_select <- motifs_per_classe[classe]

      classe_clusters <- cluster_summary$cluster_id[
        cluster_summary$dominant_class == classe
      ]
      cluster_cols <- paste0("Cluster_", classe_clusters)

      if (length(cluster_cols) == 0) {
        warning(sprintf("Class '%s' has no associated clusters", classe))
        selected_motifs[[classe]] <- character(0)
        next
      }

      tryCatch({
        classe_df <- .build_class_dataframe(
          motif_cluster,
          rownames(motif_cluster),
          cluster_cols
        )

        motif_scores <- .calculate_motif_scores(classe_df, cluster_cols, motif_cluster)

        classe_df$mean_rank <- sapply(motif_scores, function(x) x$mean_rank)
        classe_df$top_count <- sapply(motif_scores, function(x) x$top_count)

        classe_df <- classe_df[order(-classe_df$top_count, -classe_df$mean_rank), ]

        classe_selected <- c()
        for (i in seq_len(nrow(classe_df))) {
          if (length(classe_selected) >= num_motifs_to_select) {
            break
          }

          current_motif <- classe_df$motif[i]
          if (!exists(current_motif, envir = all_selected_motifs)) {
            classe_selected <- c(classe_selected, current_motif)
            assign(current_motif, TRUE, envir = all_selected_motifs)
          }
        }

        if (length(classe_selected) < num_motifs_to_select) {
          warning(sprintf(
            "Could not select %d motifs for class '%s'. Only %d selected.",
            num_motifs_to_select, classe, length(classe_selected)
          ))
        }

        selected_motifs[[classe]] <- classe_selected

      }, error = function(e) {
        warning(sprintf(
          "Error processing class '%s': %s",
          classe, e$message
        ))
        selected_motifs[[classe]] <<- character(0)
      })
    }
  }

  # ===== CASE 3: n >= m =====

  else {
    if (verbose) {
      message(sprintf(
        "\nCASE 3: n >= m\nDistributing motifs by cluster (%d >= %d)",
        n, m
      ))
    }

    case_used <- "CASE_3"

    motifs_per_cluster_base <- floor(n / m)
    remainder <- n %% m

    if (verbose) {
      message(sprintf("Base motifs per cluster: %d", motifs_per_cluster_base))
      message(sprintf("Remainder to distribute: %d", remainder))
    }

    classe_motifs <- setNames(vector("list", length(classes)), classes)
    all_selected_motifs <- new.env(hash = TRUE)

    # Select motifs for each cluster
    for (cluster_id in clusters) {
      cluster_col <- paste0("Cluster_", cluster_id)
      classe <- cluster_summary$dominant_class[
        cluster_summary$cluster_id == cluster_id
      ]

      tryCatch({
        temp_df <- data.frame(
          motif = rownames(motif_cluster),
          value = motif_cluster[[cluster_col]],
          stringsAsFactors = FALSE
        )

        temp_df <- temp_df[order(temp_df$value, decreasing = TRUE), ]

        cluster_selected <- c()
        i <- 1
        while (length(cluster_selected) < motifs_per_cluster_base &&
               i <= nrow(temp_df)) {
          current_motif <- temp_df$motif[i]
          if (!exists(current_motif, envir = all_selected_motifs)) {
            cluster_selected <- c(cluster_selected, current_motif)
            assign(current_motif, TRUE, envir = all_selected_motifs)
          }
          i <- i + 1
        }

        if (length(cluster_selected) < motifs_per_cluster_base) {
          warning(sprintf(
            "Could not select %d motifs for cluster %s. Only %d selected.",
            motifs_per_cluster_base, cluster_id, length(cluster_selected)
          ))
        }

        classe_motifs[[classe]] <- c(classe_motifs[[classe]], cluster_selected)

      }, error = function(e) {
        warning(sprintf(
          "Error processing cluster %s: %s",
          cluster_id, e$message
        ))
      })
    }

    # Distribute remainder motifs
    if (remainder > 0) {
      if (verbose) {
        message("\nDistributing extra motifs to classes (by distribution order):")
      }

      for (i in seq_len(remainder)) {
        classe_to_add <- classe_order[((i - 1) %% k) + 1]

        if (verbose) {
          message(sprintf("  Processing extra motif %d for class: %s", i, classe_to_add))
        }

        classe_clusters <- cluster_summary$cluster_id[
          cluster_summary$dominant_class == classe_to_add
        ]
        cluster_cols <- paste0("Cluster_", classe_clusters)

        tryCatch({
          classe_df <- .build_class_dataframe(
            motif_cluster,
            rownames(motif_cluster),
            cluster_cols
          )

          motif_scores <- .calculate_motif_scores(
            classe_df,
            cluster_cols,
            motif_cluster
          )

          classe_df$mean_rank <- sapply(motif_scores, function(x) x$mean_rank)
          classe_df$top_count <- sapply(motif_scores, function(x) x$top_count)

          classe_df <- classe_df[order(-classe_df$top_count, -classe_df$mean_rank), ]

          selected_extra <- NULL
          for (j in seq_len(nrow(classe_df))) {
            current_motif <- classe_df$motif[j]
            if (!exists(current_motif, envir = all_selected_motifs)) {
              selected_extra <- current_motif
              break
            }
          }

          if (!is.null(selected_extra)) {
            if (verbose) {
              message(sprintf("  Extra motif selected: %s", selected_extra))
            }
            assign(selected_extra, TRUE, envir = all_selected_motifs)
            classe_motifs[[classe_to_add]] <- c(
              classe_motifs[[classe_to_add]],
              selected_extra
            )
          } else {
            warning(sprintf(
              "Could not find an extra motif for class '%s'",
              classe_to_add
            ))
          }

        }, error = function(e) {
          warning(sprintf(
            "Error selecting extra motif for class '%s': %s",
            classe_to_add, e$message
          ))
        })
      }
    }

    selected_motifs <- classe_motifs
  }

  # ===== FORMAT OUTPUT =====

  if (verbose) {
    message("\nSummary of selected motifs by class:")
    for (classe in names(selected_motifs)) {
      motifs <- selected_motifs[[classe]]
      message(sprintf("  Class %s: %d motifs", classe, length(motifs)))
    }
    message("")
  }

  class(selected_motifs) <- c("select_motifs", "list")
  attr(selected_motifs, "n") <- n
  attr(selected_motifs, "k") <- k
  attr(selected_motifs, "m") <- m
  attr(selected_motifs, "case") <- case_used
  attr(selected_motifs, "classes_order") <- classe_order
  attr(selected_motifs, "motifs_per_class") <- sapply(selected_motifs, length)

  return(invisible(selected_motifs))
}


# ===== INTERNAL HELPER FUNCTIONS =====

#' Extract All Classes from Class Distribution String
#'
#' @keywords internal
.extract_all_classes <- function(class_dist_vector) {
  class_dist_vector <- class_dist_vector[!is.na(class_dist_vector) &
                                           class_dist_vector != ""]

  if (length(class_dist_vector) == 0) {
    return(character(0))
  }

  all_text <- paste(class_dist_vector, collapse = ",")
  classes <- unique(trimws(gsub("\\([^)]*\\)", "",
                                strsplit(all_text, ",")[[1]])))

  return(classes[classes != ""])
}


#' Build Class Data.frame for Motif Analysis
#'
#' @keywords internal
.build_class_dataframe <- function(motif_cluster, motif_names, cluster_cols) {

  missing_cols <- setdiff(cluster_cols, colnames(motif_cluster))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "Columns not found in motif_cluster: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  classe_df <- data.frame(
    motif = motif_names,
    stringsAsFactors = FALSE
  )

  for (col in cluster_cols) {
    classe_df[[col]] <- motif_cluster[[col]]
  }

  return(classe_df)
}


#' Calculate Motif Importance Scores Across Clusters
#'
#' @keywords internal
.calculate_motif_scores <- function(classe_df, cluster_cols, motif_cluster) {

  motif_scores <- lapply(seq_len(nrow(classe_df)), function(i) {
    .evaluate_motif_importance(classe_df[i, ], cluster_cols, motif_cluster)
  })

  return(motif_scores)
}


#' Evaluate Motif Importance Across Multiple Clusters
#'
#' @keywords internal
.evaluate_motif_importance <- function(row, cluster_cols, motif_cluster) {

  values <- as.numeric(row[cluster_cols])

  cluster_ranks <- numeric(length(cluster_cols))

  for (i in seq_along(cluster_cols)) {
    col_values <- motif_cluster[[cluster_cols[i]]]
    col_values_sorted <- sort(col_values, decreasing = TRUE)

    rank_position <- which(col_values_sorted == values[i])[1]
    if (is.na(rank_position)) {
      rank_position <- length(col_values) + 1
    }

    cluster_ranks[i] <- 1 - (rank_position / length(col_values))
  }

  mean_rank <- mean(cluster_ranks)
  top_count <- sum(cluster_ranks >= 0.9)

  return(list(mean_rank = mean_rank, top_count = top_count))
}


#' Print Method for select_motifs Objects
#'
#' @param x An object of class "select_motifs"
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns x
#'
#' @export
print.select_motifs <- function(x, ...) {

  cat("\n=== Selected Motifs Summary ===\n\n")

  # Check if empty result (CASE 1 failure)
  if (inherits(x, "empty_result")) {
    cat("Result: Empty (insufficient motifs)\n")
    cat("Case: ", attr(x, "case"), "\n")
    cat("Reason: ", attr(x, "reason"), "\n")
    return(invisible(x))
  }

  cat(sprintf("Total motifs selected (n):  %d\n", attr(x, "n")))
  cat(sprintf("Number of classes (k):     %d\n", attr(x, "k")))
  cat(sprintf("Number of clusters (m):    %d\n", attr(x, "m")))
  cat(sprintf("Selection case:            %s\n", attr(x, "case")))
  cat("\n")

  cat("Motifs per class:\n")
  motifs_per_class <- attr(x, "motifs_per_class")
  for (i in seq_along(x)) {
    classe <- names(x)[i]
    count <- motifs_per_class[i]
    cat(sprintf("  %s: %d motif%s\n", classe, count,
                ifelse(count == 1, "", "s")))
  }

  cat("\n")

  return(invisible(x))
}

#' Summary Method for select_motifs Objects
#'
#' @param object An object of class "select_motifs"
#' @param ... Additional arguments (ignored)
#'
#' @return Invisibly returns object
#'
#' @export
summary.select_motifs <- function(object, ...) {

  cat("\n=== Selected Motifs Detailed Summary ===\n\n")

  # Check if empty result (CASE 1 failure)
  if (inherits(object, "empty_result")) {
    cat("Result: Empty (insufficient motifs)\n")
    cat("Reason: ", attr(object, "reason"), "\n")
    return(invisible(object))
  }

  cat(sprintf("Total motifs selected: %d\n", attr(object, "n")))
  cat(sprintf("Number of classes: %d\n", attr(object, "k")))
  cat(sprintf("Number of clusters: %d\n", attr(object, "m")))
  cat(sprintf("Selection case: %s\n\n", attr(object, "case")))

  cat("Selected motifs by class:\n")
  cat("---\n")
  for (classe in names(object)) {
    motifs <- object[[classe]]
    cat(sprintf("\nClass: %s\n", classe))
    cat(sprintf("Count: %d motif%s\n", length(motifs),
                ifelse(length(motifs) == 1, "", "s")))

    if (length(motifs) > 0) {
      cat("Motifs:\n")
      for (j in seq_along(motifs)) {
        cat(sprintf("  [%d] %s\n", j, motifs[j]))
      }
    } else {
      cat("  (none selected)\n")
    }
  }

  cat("\n")

  return(invisible(object))
}
