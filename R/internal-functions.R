################################################################################
#' Internal helper functions of create_data.R
################################################################################
#'
#' .count_kmers: Count K-mers in single sequence
#'
#' @param sequences Sequences of DNAStringSet
#' @param k Integer. Length of k-mers to count
#' @param alphabet Character vector of valid nucleotides
#'
#' @return Named integer vector with counts for all possible k-mers
#'
#' @keywords internal
#'
.count_kmers <- function(sequence, k, alphabet) {

# Handle Biostrings objects
if (methods::is(sequence, "XString")) {
  # Convert Biostrings XString to character
  sequences <- as.character(sequence)
}
# Sequence processing
sequence <- toupper(sequence) # Convert to uppercase
seq_length <- nchar(sequence) # Sequence length

# Check if sequence contains only valid nucleotides
sequence_chars <- unique(strsplit(sequence, "")[[1]])
invalid_chars <- setdiff(sequence_chars, alphabet)

if (length(invalid_chars) > 0) {
  # Remove invalid characters from sequence
  pattern <- paste0("[^", paste(alphabet, collapse = ""), "]")
  sequence <- gsub(pattern, "", sequence)

  # Update sequence length after cleaning
  seq_length <- nchar(sequence)
}

# Generate all possible k-mers from alphabet
all_kmers <- as.character(Biostrings::mkAllStrings(alphabet, k))

# Extract all k-mers from sequence
n_kmers_in_seq <- seq_length - k + 1
extracted_kmers <- character(n_kmers_in_seq)

for (i in seq_len(n_kmers_in_seq)) {
  extracted_kmers[i] <- substr(sequence, i, i + k - 1)
}

# Count occurrences of each extracted k-mer using table
kmer_table <- table(extracted_kmers)


# Initialize vector with all possible k-mers set to 0
complete_counts <- setNames(rep(0L, length(all_kmers)), all_kmers)
matching_kmers <- intersect(names(kmer_table), all_kmers)

if (length(matching_kmers) > 0) {
  complete_counts[matching_kmers] <- as.integer(kmer_table[matching_kmers])
}

# Return results
return(complete_counts)
}

################################################################################
#'
#' .process_sequences: Process Sequences from a Single FASTA File
#'
#' @param sequences A DNAStringSet object containing sequences
#' @param all_kmers Character vector of all possible k-mers
#' @param k Integer k-mer length
#' @param alphabet Character vector of nucleotides
#' @param class_name Character string with the class label
#'
#' @note Require .count_kmers aux. function.
#' @return kmer_df (data frame) and metadata_df (data frame)
#'
#' @keywords internal
#'
.process_sequences <- function(sequences, all_kmers, k, alphabet, class_name) {

  n_seq <- length(sequences) #length of sequences
  n_kmers <- length(all_kmers) #length of k-mers

  # Pre-allocate matrix
  kmer_matrix <- matrix(0L, nrow = n_seq, ncol = n_kmers + 1L,
                        dimnames = list(NULL, c(all_kmers, "CLASS"))) #create matrix w/ n_seq rows, n_kmers cols and class column

  seq_names <- character(n_seq)
  seq_lengths <- integer(n_seq)

  # For each sequence
  for (seq_idx in seq_len(n_seq)) {
    seq_name <- names(sequences)[seq_idx]
    seq_obj <- sequences[[seq_idx]]
    seq_length <- length(seq_obj)

    kmer_counts <- .count_kmers(seq_obj, k, alphabet)

    # Store in matrix
    kmer_matrix[seq_idx, seq(along = all_kmers)] <- as.integer(kmer_counts[all_kmers])
    kmer_matrix[seq_idx, "CLASS"] <- class_name

    seq_names[seq_idx] <- seq_name
    seq_lengths[seq_idx] <- seq_length
  }

  # Convert to dataframe
  kmer_df <- as.data.frame(kmer_matrix, stringsAsFactors = FALSE)
  kmer_df[, seq(along = all_kmers)] <- lapply(kmer_df[, seq(along = all_kmers)],
                                              as.integer)

  # Create metadata
  metadata_df <- data.frame(
    sequence_name = seq_names,
    length = seq_lengths,
    class = rep(class_name, n_seq),
    stringsAsFactors = FALSE
  )

  # Return results
  list(kmers = kmer_df, metadata = metadata_df)
}

################################################################################
#' Internal helper functions of cluster_dendrogram.R
################################################################################
#'
#' .count_leaves: Count Leaves in a Dendrogram
#'
#' @param dend A dendrogram object or sub-dendrogram node
#' @return Integer. The total number of leaves in the dendrogram.
#' @keywords internal
#'
.count_leaves <- function(dend) {
  if (stats::is.leaf(dend)) return(1L)
  sum(vapply(dend, .count_leaves, FUN.VALUE = integer(1)))
}
################################################################################
#'
#'.calc_homogeneity: Calculate Homogeneity of a Class Vector
#'
#' @param classes Character vector of class labels
#' @return Numeric. Homogeneity value between 0 and 1, or NA if the input
#'   vector is empty.
#' @keywords internal
#'
.calc_homogeneity <- function(classes) {
  if (length(classes) == 0) return(NA_real_)
  max(table(classes)) / length(classes)
}
################################################################################
#'
#' .get_dominant: Get Dominant Class from a Class Vector
#'
#' @param classes Character vector of class labels
#' @return Character. The name of the most frequent class, or NA if the input
#'   vector is empty.
#' @keywords internal
#'
.get_dominant <- function(classes) {
  if (length(classes) == 0) return(NA_character_)
  names(which.max(table(classes)))
}
################################################################################
#'
#' .find_segments: Find Contiguous Segments in an Index Vector
#'
#' @param indices Integer vector of indices (not necessarily sorted or contiguous)
#' @return List of integer vectors, where each element is a contiguous segment.
#'   Example: \code{c(1,2,3,7,8)} returns \code{list(c(1,2,3), c(7,8))}.
#' @keywords internal
#'
.find_segments <- function(indices) {
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
################################################################################
#'
#' .is_complete_class_cluster: Check if Cluster Contains All Members of Dominant Class
#'
#' @param cluster_classes Character vector of class labels within the cluster
#' @param total_counts Named integer vector with total counts for each class
#'   in the complete dataset
#' @return Logical. TRUE if the cluster contains all members of its dominant
#'   class, FALSE otherwise.
#' @keywords internal
#'
.is_complete_class_cluster <- function(cluster_classes, total_counts) {
  class_table <- table(cluster_classes)
  dominant_class <- names(which.max(class_table))
  class_table[dominant_class] == total_counts[dominant_class]
}
################################################################################
#'
#' .create_cluster: Create a Cluster Object
#' @param indices Integer vector of indices belonging to the cluster
#' @param classes Character vector of class labels corresponding to the indices
#' @param counter Integer. Unique cluster identifier
#' @return List with three elements:
#'   \itemize{
#'     \item \code{indices}: Integer vector of cluster member indices
#'     \item \code{classes}: Character vector of class labels
#'     \item \code{id}: Integer cluster ID
#'   }
#' @keywords internal
#'
.create_cluster <- function(indices, classes, counter) {
  list(indices = indices, classes = classes, id = counter)
}

################################################################################
#' Internal helper functions of calculate_cluster_motifs.R
################################################################################
#'
#'.normalize_motif_matrix: Normalize Motif Matrix Using Min-Max Normalization
#'
#' @param motif_raw Numeric matrix to normalize (by column)
#' @return Normalized matrix with values between 0 and 1
#' @keywords internal
.normalize_motif_matrix <- function(motif_raw) {

  col_mins <- apply(motif_raw, 2, min)
  col_maxs <- apply(motif_raw, 2, max)
  col_ranges <- col_maxs - col_mins

  # Identify columns with zero range (all values identical)
  zero_range <- col_ranges == 0
  col_ranges[zero_range] <- 1  # Avoid division by zero

  # Apply min-max normalization using sweep (vectorized)
  motif_normalized <- sweep(
    sweep(motif_raw, 2, col_mins, "-"),
    2,
    col_ranges,
    "/"
  )

  # Set zero-range columns to 0.5
  motif_normalized[, zero_range] <- 0.5

  return(motif_normalized)
}

################################################################################
#' Internal helper functions of select_motifs.R
################################################################################
#'
#'.select_by_clas_fast: Select Top Motifs by Class (Fast Method)
#'
#' @param motif_cluster A data.frame or matrix with motifs as rows and clusters
#'   as columns. Row names should be motif identifiers.
#' @param cluster_to_class A named vector mapping cluster names to class labels.
#'   Names should match column names in \code{motif_cluster}.
#' @param classe_order A character vector specifying the order of classes for
#'   motif distribution.
#' @param n An integer specifying the total number of motifs to select.
#'
#' @return A named list where each element corresponds to a class and contains
#'   a character vector of selected motif names. List names match \code{classe_order}.
#'
#' @keywords internal
.select_by_class_fast <- function(motif_cluster, cluster_to_class, classe_order, n, verbose) {

  # Calculate distribution
  n_classes <- length(classe_order)
  motifs_per_classe <- rep(n %/% n_classes, n_classes)
  remainder <- n %% n_classes
  if (remainder > 0) {
    motifs_per_classe[seq_len(remainder)] <- motifs_per_classe[seq_len(remainder)] + 1L
  }
  names(motifs_per_classe) <- classe_order

  # Group columns by class
  class_cols <- split(names(cluster_to_class), cluster_to_class)
  selected_motifs <- setNames(vector("list", n_classes), classe_order)
  all_motifs <- rownames(motif_cluster)
  used_motifs <- setNames(rep(FALSE, length(all_motifs)), all_motifs)

  for (classe in classe_order) {
    cols <- class_cols[[classe]]
    if (length(cols) == 0) {
      selected_motifs[[classe]] <- character(0)
      next
    }

    scores <- rowSums(motif_cluster[, cols, drop = FALSE])
    top_idx <- order(scores, decreasing = TRUE)
    selected <- character(0)

    for (idx in top_idx) {
      if (length(selected) >= motifs_per_classe[classe]) break
      motif <- all_motifs[idx]
      if (!used_motifs[motif]) {
        selected <- c(selected, motif)
        used_motifs[motif] <- TRUE
      }
    }
    selected_motifs[[classe]] <- selected
  }

  selected_motifs
}
################################################################################
#'
#'.select_by_cluster_fast: Select Top Motifs by Cluster (Fast Method)
#' @param motif_cluster A data.frame or matrix with motifs as rows and clusters
#'   as columns. Row names should be motif identifiers.
#' @param cluster_to_class A named vector mapping cluster names to class labels.
#'   Names should match column names in \code{motif_cluster}.
#' @param classe_order A character vector specifying the order of classes for
#'   remainder distribution priority.
#' @param n An integer specifying the total number of motifs to select.
#' @param m An integer specifying the total number of clusters.
#'
#' @return A named list where each element corresponds to a class and contains
#'   a character vector of selected motif names.
#'
#' @keywords internal
.select_by_cluster_fast <- function(motif_cluster, cluster_to_class, classe_order, n, m) {

  motifs_per_cluster <- n %/% m
  remainder <- n %% m

  classe_motifs <- setNames(vector("list", length(classe_order)), classe_order)
  all_motifs <- rownames(motif_cluster)
  used_motifs <- setNames(rep(FALSE, length(all_motifs)), all_motifs)

  # Process each cluster
  for (col in names(cluster_to_class)) {
    classe <- cluster_to_class[col]
    values <- motif_cluster[[col]]
    top_idx <- order(values, decreasing = TRUE)
    selected <- character(0)

    for (idx in top_idx) {
      if (length(selected) >= motifs_per_cluster) break
      motif <- all_motifs[idx]
      if (!used_motifs[motif]) {
        selected <- c(selected, motif)
        used_motifs[motif] <- TRUE
      }
    }
    classe_motifs[[classe]] <- c(classe_motifs[[classe]], selected)
  }

  # Distribute remainder
  if (remainder > 0) {
    for (i in seq_len(remainder)) {
      classe_to_add <- classe_order[((i - 1L) %% length(classe_order)) + 1L]
      cols <- names(cluster_to_class)[cluster_to_class == classe_to_add]
      if (length(cols) == 0) next

      scores <- rowSums(motif_cluster[, cols, drop = FALSE])
      top_idx <- order(scores, decreasing = TRUE)

      for (idx in top_idx) {
        motif <- all_motifs[idx]
        if (!used_motifs[motif]) {
          classe_motifs[[classe_to_add]] <- c(classe_motifs[[classe_to_add]], motif)
          used_motifs[motif] <- TRUE
          break
        }
      }
    }
  }

  classe_motifs
}

################################################################################
#' Internal helper functions of select_sample_train_validation.R
################################################################################
#' Select k sequences per class
#'
#' .select_by_class:  Select k Sequences per Class
#'
#' @param data A data.frame containing at least a 'class' column and sequence
#'   information.
#' @param k Integer. Number of sequences to select per class.
#' @param dataset_name Character. Name identifier for the dataset (used in messages).
#'   Default is empty string.
#'
#' @return A data.frame with selected sequences from all classes combined.
#'   Returns empty data.frame if no sequences can be selected. Row names are reset.
#'
#' @keywords internal
.select_by_class <- function(data, k, dataset_name = "") {
  classes <- unique(data$class)

  sequences_selected <- lapply(classes, function(cls) {
    class_data <- data[data$class == cls, ]
    n_available <- nrow(class_data)
    k_adjusted <- min(k, n_available)

    if (k_adjusted == 0) {
      return(NULL)
    }

    set.seed(123)  # For reproducibility in tests
    selected_indices <- sample(nrow(class_data), k_adjusted)
    class_data[selected_indices, ]
  })

  sequences_selected <- do.call(rbind, sequences_selected[!sapply(sequences_selected, is.null)])
  rownames(sequences_selected) <- NULL

  sequences_selected
}
################################################################################
#' Internal helper functions of train_models_rf_xgboost.R
################################################################################
#'
#' .rename_class_column: Rename CLASS column to class
#' @param data A data frame that may contain a column named "CLASS".
#' @return A data frame with the "CLASS" column renamed to "class" if it exists
#'   and there is no existing "class" column. If "CLASS" does not exist or
#'   "class" already exists, returns the data frame unchanged.
#' @keywords internal
#'
.rename_class_column <- function(data) {
  if ("CLASS" %in% colnames(data) && !("class" %in% colnames(data))) {
    colnames(data)[colnames(data) == "CLASS"] <- "class"
  }
  data
}
################################################################################
#'
#' .validate_columns: Validate required columns
#'
#' @param data A data frame to be validated.
#' @param required_cols A character vector containing the names of columns that
#'   must be present in the data frame.
#' @param dataset_name A character string with the name of the dataset, used in
#'   the error message for identification purposes.
#' @return NULL (invisibly). The function is called for its side effect of
#'   stopping execution with an error message if required columns are missing.
#'   If all required columns are present, the function completes silently.
#' @keywords internal
#'
.validate_columns <- function(data, required_cols, dataset_name) {
  missing <- setdiff(required_cols, colnames(data))
  if (length(missing) > 0) {
    stop("Missing columns in ", dataset_name, ": ",
         paste(missing, collapse = ", "), call. = FALSE)
  }
}

################################################################################
#' Internal helper functions of kmer_analysis.R
################################################################################
#'
#' .create_motifs_rank: Create motifs rank table with sum and count statistics
#'
#' @param data A data frame containing motif columns and a class column.
#' @param class_col A character string specifying the name of the class column.
#'   Default: "CLASS".
#'
#' @return A data frame with motif statistics where:
#'   \itemize{
#'     \item \code{motif}: Character vector of motif names
#'     \item Class columns: Formatted strings "sum|count" for each class
#'     \item \code{mean}: Numeric vector of mean sums across all classes
#'   }
#'   Each class column contains values formatted as "sum|count", where sum is
#'   the total sum of motif values and count is the number of non-zero occurrences
#'   for that motif in that class.
#'
#' @keywords internal
.create_motifs_rank <- function(data, class_col = "CLASS") {
  class_col_idx <- which(colnames(data) == class_col)
  motif_cols <- setdiff(seq_len(ncol(data)), class_col_idx)

  classes <- data[[class_col]]
  motifs <- data[, motif_cols, drop = FALSE]
  motif_names <- colnames(motifs)
  unique_classes <- unique(classes)

  sums_mat <- matrix(0, nrow = length(motif_names), ncol = length(unique_classes),
                     dimnames = list(motif_names, unique_classes))
  counts_mat <- matrix(0L, nrow = length(motif_names), ncol = length(unique_classes),
                       dimnames = list(motif_names, unique_classes))

  for (class_name in unique_classes) {
    class_rows <- classes == class_name
    class_data <- motifs[class_rows, , drop = FALSE]
    sums_mat[, class_name] <- colSums(class_data, na.rm = TRUE)
    counts_mat[, class_name] <- colSums(class_data > 0, na.rm = TRUE)
  }

  value_width <- nchar(as.character(ceiling(max(sums_mat))))
  formatted_mat <- matrix("", nrow = nrow(sums_mat), ncol = ncol(sums_mat),
                          dimnames = dimnames(sums_mat))

  for (i in seq_len(nrow(sums_mat))) {
    for (j in seq_len(ncol(sums_mat))) {
      formatted_mat[i, j] <- sprintf("%0*d|%d", value_width,
                                     round(sums_mat[i, j]), counts_mat[i, j])
    }
  }

  result <- as.data.frame(formatted_mat, stringsAsFactors = FALSE)
  result$motif <- motif_names
  result$mean <- round(rowMeans(sums_mat, na.rm = TRUE), 2)
  result <- result[, c("motif", unique_classes, "mean")]
  rownames(result) <- NULL
  return(result)
}

################################################################################
#' Internal helper functions of kmers_in_seq.R
################################################################################
#'
#' .read_fasta_sequences: Read and consolidate FASTA sequences from directory
#'
#' @param input_dir Character string specifying the directory containing FASTA files
#' @return A list with two elements: 'sequences' (character vector of uppercase DNA sequences) and 'names' (character vector of sequence identifiers)
#' @keywords internal
.read_fasta_sequences <- function(input_dir) {
  fasta_files <- list.files(input_dir, pattern = "\\.(fasta)$",
                            full.names = TRUE, ignore.case = TRUE)

  all_sequences <- list()
  all_names <- character()

  for (fasta_file in fasta_files) {
    seqs <- seqinr::read.fasta(fasta_file, seqtype = "DNA", as.string = TRUE,
                               forceDNAtolower = FALSE)
    for (seq_name in names(seqs)) {
      all_sequences[[seq_name]] <- toupper(seqs[[seq_name]])
      all_names <- c(all_names, seq_name)
    }
  }

  return(list(sequences = unlist(all_sequences, use.names = FALSE), names = all_names))
}

################################################################################
#'
#' .find_motifs_parallel: Find motif occurrences in sequences using parallel processing
#'
#' @param sequences Character vector of DNA sequences
#' @param motifs Character vector of motif patterns to search
#' @param sequence_names Character vector of sequence identifiers
#' @param class_lookup Named vector for mapping sequence names to classes
#' @param length_lookup Named vector for mapping sequence names to lengths
#' @param n_cores Integer specifying number of CPU cores to use
#' @param dataset Character or NULL. Dataset label ("training" or "validation").
#'   Default: NULL.
#'
#' @return A data.frame with columns: motif, sequence_name, class, position_start, position_end, sequence_length, and optionally dataset
#'
#' @keywords internal
.find_motifs_parallel <- function(sequences, motifs, sequence_names,
                                  class_lookup, length_lookup, n_cores,dataset = NULL) {

  process_motif <- function(motif) {
    results_list <- list()
    counter <- 0

    for (j in seq_along(sequences)) {
      positions <- stringi::stri_locate_all_fixed(sequences[j], motif)[[1]]

      if (!is.na(positions[1, 1])) {
        counter <- counter + 1
        seq_name <- sequence_names[j]
        seq_class <- if (!is.null(class_lookup)) {
          cls <- class_lookup[seq_name]
          if (is.na(cls)) "Unknown" else as.character(cls)
        } else {
          "Unknown"
        }

        seq_length <- if (!is.null(length_lookup)) {
          len <- length_lookup[seq_name]
          if (is.na(len)) NA_integer_ else as.integer(len)
        } else {
          NA_integer_
        }

        temp_df <- data.frame(
          motif = rep(motif, nrow(positions)),
          sequence_name = rep(seq_name, nrow(positions)),
          class = rep(seq_class, nrow(positions)),
          position_start = as.integer(positions[, 1]),
          position_end = as.integer(positions[, 2]),
          sequence_length = rep(seq_length, nrow(positions)),
          stringsAsFactors = FALSE
        )
        # Add dataset column if provided
        if (!is.null(dataset)) {
          temp_df$dataset <- rep(dataset, nrow(positions))
        }

        results_list[[counter]] <- temp_df
      }
    }

    if (counter > 0) do.call(rbind, results_list) else NULL
  }

  # Parallel or sequential
  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("sequences", "sequence_names", "class_lookup", "length_lookup", "dataset"),
                            envir = environment())
    parallel::clusterEvalQ(cl, library(stringi))
    results_list <- parallel::parLapply(cl, motifs, process_motif)
  } else {
    results_list <- lapply(motifs, process_motif)
  }

  # Combine
  results_list <- results_list[!sapply(results_list, is.null)]

  if (length(results_list) > 0) {
    result_df <- do.call(rbind, results_list)
    rownames(result_df) <- NULL
  } else {
    result_df <- data.frame(
      motif = character(0),
      sequence_name = character(0),
      class = character(0),
      position_start = integer(0),
      position_end = integer(0),
      sequence_length = integer(0),
      stringsAsFactors = FALSE
    )
    if (!is.null(dataset)) {
      result_df$dataset <- character(0)
    }
  }

  return(result_df)
}
