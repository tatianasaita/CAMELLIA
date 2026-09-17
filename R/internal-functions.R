################################################################################
# Internal helper functions of create_data.R
################################################################################
#' .generateCombinations: Combinations of Nucleotides for K-mers
#'
#' Generates all possible k-mer combinations for a given word length using
#' the DNA alphabet (A, T, C, G).
#'
#' @param word Integer. The k-mer length (word size).
#'
#' @return Character vector of all possible k-mer combinations, sorted
#'   alphabetically.
#'
#' @keywords internal
.generateCombinations <- function(word) {
 nucleotides <- c("A", "T", "C", "G")
 combinations <- do.call(expand.grid,
                       replicate(word, nucleotides, simplify = FALSE))
 sequences <- apply(combinations, 1, paste0, collapse = "")
 sort(sequences)
 }
#'
################################################################################
#'
#' .createNet_fast: Create Nets
#'
#' Builds a directed graph (network) of k-mer transitions from a DNA sequence,
#' used for network-based sequence representation.
#'
#' @param word Integer. K-mer length used to build each node.
#' @param step Integer. Step size between consecutive k-mer windows.
#' @param sequence Character vector. The DNA sequence, split into individual
#'   characters.
#' @param all_vertices Character vector or NULL. Optional set of vertex names
#'   that must be present in the resulting graph, even if not observed in
#'   \code{sequence}. Default: NULL.
#'
#' @return An \code{igraph} directed graph object representing k-mer
#'   transitions.
#'
#' @keywords internal
.createNet_fast <- function(word, step, sequence, all_vertices = NULL) {
   cont <- length(sequence)
   
   starts <- seq(1, cont - (word * 2) + 1, by = step)
   
   make_kmers <- function(idx_starts) {
     sapply(idx_starts, function(i) {
       paste0(sequence[i:(i + word - 1)], collapse = "")
     })
   }
   
   from_kmers <- make_kmers(starts)
   to_kmers   <- make_kmers(starts + word)
   
   edges_vec <- as.vector(rbind(from_kmers, to_kmers))
   
   net <- igraph::graph(edges = edges_vec, directed = TRUE)
   
   if (!is.null(all_vertices)) {
     missing_v <- setdiff(all_vertices, igraph::V(net)$name)
     if (length(missing_v) > 0) {
       net <- igraph::add_vertices(net, length(missing_v), name = missing_v)
     }
   }
   
   return(net)
 }


################################################################################
#'
#' .count_kmers: Count K-mers in a Single Sequence
#'
#' @param sequence A single DNA sequence, either a \code{DNAString} object or
#'   a character string/object coercible to one.
#' @param k Integer. Length of k-mers to count.
#' @param alphabet Character vector of valid nucleotides.
#'
#' @return Named integer vector with counts for all possible k-mers.
#'
#' @keywords internal
.count_kmers <- function(sequence, k, alphabet) {

# Handle Biostrings objects
  if (!methods::is(sequence, "DNAString")) {
    sequence <- Biostrings::DNAString(toupper(as.character(sequence)))
  }

  # Generate complete k-mer vocabulary and initialise all counts to zero
  all_kmers      <- as.character(Biostrings::mkAllStrings(alphabet, k))
  complete_counts <- setNames(rep(0L, length(all_kmers)), all_kmers)

  # Count k-mers using Biostrings C-level routine
  kmer_counts <- Biostrings::oligonucleotideFrequency(sequence, width = k,
                                                      as.prob = FALSE)

  # Fill only the k-mers present in our alphabet vocabulary
  common_kmers <- intersect(names(kmer_counts), all_kmers)
  complete_counts[common_kmers] <- as.integer(kmer_counts[common_kmers])

  return(complete_counts)
}

################################################################################
#'
#' .process_sequences: Process Sequences from a Single FASTA File
#'
#' @param sequences A DNAStringSet object containing sequences.
#' @param all_kmers Character vector of all possible k-mers.
#' @param k Integer. K-mer length.
#' @param alphabet Character vector of nucleotides.
#' @param class_name Character string with the class label.
#'
#' @return A list with two elements: \code{kmers} (data frame of k-mer
#'   counts) and \code{metadata} (data frame of sequence metadata).
#'
#' @keywords internal
.process_sequences <- function(sequences, all_kmers, k, alphabet, class_name) {

    n_seq   <- length(sequences)
    n_kmers <- length(all_kmers)

    # Count k-mers for all sequences at once (returns a matrix: seqs x k-mers)
    # oligonucleotideFrequency accepts DNAStringSet directly — no loop needed
    kmer_matrix_raw <- Biostrings::oligonucleotideFrequency(sequences, width = k,
                                                            as.prob = FALSE)

    # Ensure every column of the vocabulary is present (zero-fill absent k-mers)
    complete_matrix <- matrix(
      0L,
      nrow     = n_seq,
      ncol     = n_kmers,
      dimnames = list(names(sequences), all_kmers)
    )

    common_kmers <- intersect(colnames(kmer_matrix_raw), all_kmers)
    complete_matrix[, common_kmers] <- kmer_matrix_raw[, common_kmers]

    # Convert to data.frame and append CLASS column
    kmer_df       <- as.data.frame(complete_matrix, check.names = FALSE)
    kmer_df$CLASS <- class_name

    # Build metadata using Biostrings::width() — avoids converting to character
    metadata_df <- data.frame(
      sequence_name    = names(sequences),
      length           = Biostrings::width(sequences),
      class            = class_name,
      stringsAsFactors = FALSE,
      row.names        = NULL
    )

    list(kmers = kmer_df, metadata = metadata_df)
  }

################################################################################
# Internal helper functions of create_dendrogram_cent.R
################################################################################
#'
#' .plot_dendrogram: Plot Dendrogram with Class Color Legend
#'
#' Plots a dendrogram with colored leaves and a legend mapping colors to
#' class names.
#'
#' @param dend A dendrogram object, typically produced by
#'   \code{\link[stats]{as.dendrogram}} with leaf colors already applied
#'   (e.g., via \code{dendextend::set(dend, "labels_col", ...)}).
#' @param base_colors Named character vector mapping class names to colors,
#'   used to build the legend.
#'
#' @return NULL (invisibly). Called for its side effect of producing a plot
#'   on the current graphics device.
#'
#' @keywords internal
.plot_dendrogram <- function(dend, base_colors) {
  graphics::plot(dend, main = "Dendrogram", ylab = "Height")
  graphics::legend("topright",
                   legend  = names(base_colors),
                   col     = base_colors,
                   pch     = 15,
                   pt.cex  = 2,
                   cex     = 0.8,
                   title   = "Classes",
                   bg      = "white",
                   box.lty = 1)
}
 
################################################################################
# Internal helper functions of cluster_dendrogram.R
################################################################################
#'
#' .find_segments: Find Contiguous Segments in Index Vector
#'
#' @param indices Integer vector of indices to segment. May contain duplicates
#'   and be unsorted; the function will sort them internally.
#'
#' @return A list of integer vectors, where each element is a contiguous segment
#'   of indices. Returns an empty list if input is empty.
#'
#' @keywords internal
.find_segments <- function(indices) {
  if (length(indices) == 0) return(list())

  indices <- sort(indices[!is.na(indices)])

  if (length(indices) == 0) return(list())
  segments <- list()
  segment_start <- indices[1]
  segment_end <- indices[1]

  for (i in 2:length(indices)) {
    current <- indices[i]
    if (is.na(current)) next

    if (current == segment_end + 1) {
      segment_end <- current
    } else {
      segments[[length(segments) + 1]] <- segment_start:segment_end
      segment_start <- current
      segment_end <- current
    }
  }

  segments[[length(segments) + 1]] <- segment_start:segment_end

  return(segments)
}

################################################################################
#'
#' .create_cluster: Create Cluster Object
#' 
#' @param indices Integer vector of element indices assigned to this cluster.
#'   These are positions in the original dendrogram or data ordering.
#' @param classes Character vector of class labels corresponding to each index.
#'   Must have the same length as \code{indices}.
#' @param id Integer. Unique cluster identifier used for reference and result
#'   reporting.
#'
#' @return A list with three components:
#'   \itemize{
#'     \item \code{id}: Numeric cluster identifier
#'     \item \code{indices}: Integer vector of element positions in the cluster
#'     \item \code{classes}: Character vector of class labels for each element
#'   }
#'
#' @keywords internal
.create_cluster <- function(indices, classes, id) {
  list(
    id = id,
    indices = indices,
    classes = classes
  )
}
################################################################################
#'
#' .calc_homogeneity: Calculate Cluster Homogeneity
#'
#' @param classes Character vector of class labels for cluster elements.
#'
#' @return Numeric value between 0 and 1 representing cluster homogeneity.
#'   Returns 0 for empty input.
#'
#' @keywords internal
.calc_homogeneity <- function(classes) {
  if (length(classes) == 0) return(0)
  class_counts <- table(classes)
  max_count <- max(class_counts)
  return(max_count / length(classes))
}
################################################################################
#'
#' .is_complete_class_cluster: Check if Cluster Contains Complete Class
#' 
#' @param cluster_classes Character vector of class labels in the cluster being
#'   evaluated.
#' @param class_total_counts Named integer vector (typically from \code{table()})
#'   containing the total count of each class across all elements. Names should
#'   match class labels in \code{cluster_classes}.
#'
#' @return Logical. \code{TRUE} if the cluster contains all elements of at least
#'   one class, \code{FALSE} otherwise.
#'
#' @keywords internal
.is_complete_class_cluster <- function(cluster_classes, class_total_counts) {
  class_counts <- table(cluster_classes)
  for (cls in names(class_counts)) {
    if (class_counts[cls] == class_total_counts[cls]) {
      return(TRUE)
    }
  }
  return(FALSE)
}
################################################################################
#'
#' .get_dominant: Get Dominant Class
#' 
#' @param classes Character vector of class labels.
#'
#' @return Character. The name of the dominant (most frequent) class. If multiple
#'   classes have equal maximum frequency, returns the first in lexicographic
#'   order.
#' @keywords internal
.get_dominant <- function(classes) {
  class_counts <- table(classes)
  dominant <- names(which.max(class_counts))
  return(as.character(dominant))
}

################################################################################
# Internal helper functions of calculate_cluster_motifs.R
################################################################################
#'
#' .normalize_motif_matrix: Normalize Motif Matrix Using Min-Max Normalization
#'
#' @param motif_raw Numeric matrix to normalize (by column).
#'
#' @return Normalized matrix with values between 0 and 1.
#'
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
# Internal helper functions of select_motifs.R
################################################################################
#'
#' .select_by_class_fast: Select Top Motifs by Class (Fast Method)
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
.select_by_class_fast <- function(motif_cluster, cluster_to_class, classe_order, n) {

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
#' .select_by_cluster_fast: Select Top Motifs by Cluster (Fast Method)
#'
#' @param motif_cluster A data.frame or matrix with motifs as rows and clusters
#'   as columns. Row names should be motif identifiers.
#' @param cluster_to_class A named vector mapping cluster names to class labels.
#'   Names should match column names in \code{motif_cluster}.
#' @param classe_order A character vector specifying the order of classes for
#'   remainder distribution priority.
#' @param n An integer specifying the total number of motifs to select.
#' @param m An integer specifying the total number of clusters.
#' @param cluster_homogeneity Named numeric vector or NULL. Homogeneity score
#'   per cluster, with names matching the columns of \code{motif_cluster}. Used
#'   to prioritize less homogeneous clusters when distributing the remainder.
#'   Default: NULL.
#'
#' @return A named list where each element corresponds to a class and contains
#'   a character vector of selected motif names.
#'
#' @keywords internal
.select_by_cluster_fast <- function(motif_cluster, cluster_to_class, classe_order, n, m,
                                    cluster_homogeneity = NULL) {
  
  motifs_per_cluster <- n %/% m
  remainder          <- n %% m
  
  classe_motifs <- setNames(vector("list", length(classe_order)), classe_order)
  all_motifs    <- rownames(motif_cluster)
  used_motifs   <- setNames(rep(FALSE, length(all_motifs)), all_motifs)
  
  # Select base motifs for each cluster
  for (col in names(cluster_to_class)) {
    classe <- cluster_to_class[col]
    values <- motif_cluster[, col, drop = TRUE]
    
    if (is.null(values) || length(values) == 0) {
      warning(sprintf(
        ".select_by_cluster_fast: column '%s' not found or empty in motif_cluster - skipped.", col
      ))
      next
    }
    
    if (!is.numeric(values)) values <- suppressWarnings(as.numeric(values))
    
    if (all(is.na(values))) {
      warning(sprintf(
        ".select_by_cluster_fast: column '%s' is entirely NA - skipped.", col
      ))
      next
    }
    
    top_idx  <- order(values, decreasing = TRUE, na.last = TRUE)
    selected <- character(0)
    
    for (idx in top_idx) {
      if (length(selected) >= motifs_per_cluster) break
      motif <- all_motifs[idx]
      if (!used_motifs[motif]) {
        selected        <- c(selected, motif)
        used_motifs[motif] <- TRUE
      }
    }
    classe_motifs[[classe]] <- c(classe_motifs[[classe]], selected)
  }
  
  # Distribui o resto por homogeneidade crescente
  if (remainder > 0 && !is.null(cluster_homogeneity)) {
    
    # cluster_homogeneity: named numeric vector - names match motif_cluster columns
    valid_hom <- cluster_homogeneity[names(cluster_homogeneity) %in% names(cluster_to_class)]
    
    if (length(valid_hom) == 0) {
      warning(".select_by_cluster_fast: 'cluster_homogeneity' has no names matching the clusters - remainder ignored.")
    } else {
      # Sort clusters from least to most homogeneous
      clusters_by_hom <- names(sort(valid_hom, decreasing = FALSE))
      
      for (i in seq_len(remainder)) {
        col    <- clusters_by_hom[[i]]           # i-th cluster in ascending order
        classe <- cluster_to_class[col]
        values <- motif_cluster[, col, drop = TRUE]
        
        if (!is.numeric(values)) values <- suppressWarnings(as.numeric(values))
        
        top_idx <- order(values, decreasing = TRUE, na.last = TRUE)
        
        for (idx in top_idx) {
          motif <- all_motifs[idx]
          if (!used_motifs[motif]) {
            classe_motifs[[classe]] <- c(classe_motifs[[classe]], motif)
            used_motifs[motif]      <- TRUE
            break
          }
        }
      }
    }
    
  } else if (remainder > 0 && is.null(cluster_homogeneity)) {
    warning(".select_by_cluster_fast: 'cluster_homogeneity' not provided - remainder will not be distributed.")
  }
  
  classe_motifs
}
################################################################################
# Internal helper functions of select_train_test.R
################################################################################
#' .select_by_class: Select Sequences for Training (Internal)
#'
#' Randomly samples \code{n_train} rows from \code{data}, which must already
#' be filtered to a single class. The seed is intentionally not set here —
#' reproducibility is the caller's responsibility (set seed before calling
#' \code{select_train_test}).
#'
#' @param data Data frame of metadata for a single class.
#' @param n_train Integer. Number of sequences to select.
#' @param dataset_name Character label written to the \code{dataset} column
#'   (default: \code{"train"}).
#'
#' @return A data frame with \code{n_train} rows (or fewer if \code{data} has
#'   fewer rows) and a \code{dataset} column set to \code{dataset_name}.
#'
#' @keywords internal
.select_by_class <- function(data, n_train, dataset_name = "train") {
  n_available        <- nrow(data)
  n_select           <- min(as.integer(n_train), n_available)
  selected           <- data[sample(n_available, n_select), ]
  selected$dataset   <- dataset_name
  rownames(selected) <- NULL
  selected
}
################################################################################
#'
#' .align_kmer_columns: Align K-mer Columns Between Train and Test Matrices
#'
#' Checks whether the two matrices share the same k-mer columns and, if not,
#' subsets both to their intersection, emitting warnings as appropriate.
#'
#' @param train_dataset Matrix. K-mer counts for the training set.
#' @param test_dataset  Matrix. K-mer counts for the test set.
#'
#' @return A named list with elements \code{train_dataset} and
#'   \code{test_dataset}, both restricted to common columns.
#'
#' @keywords internal
.align_kmer_columns <- function(train_dataset, test_dataset) {
  n_kmers_train <- ncol(train_dataset)
  n_kmers_test  <- ncol(test_dataset)
  
  if (n_kmers_train != n_kmers_test) {
    warning("K-mer count mismatch: train=", n_kmers_train,
            ", test=", n_kmers_test, call. = FALSE)
    
    common_kmers <- intersect(colnames(train_dataset), colnames(test_dataset))
    
    if (length(common_kmers) < n_kmers_train * 0.7) {
      warning("Only ", length(common_kmers), " common k-mers found (",
              round(length(common_kmers) / n_kmers_train * 100, 1), "%)",
              call. = FALSE)
    }
    
    train_dataset <- train_dataset[, common_kmers, drop = FALSE]
    test_dataset  <- test_dataset[,  common_kmers, drop = FALSE]
  }
  
  list(train_dataset = train_dataset, test_dataset = test_dataset)
}

################################################################################
#' Extract Class Label from Sequence Name
#'
#' Extracts the class label from a sequence name by removing everything from
#' the first "." onward (e.g., \code{"type1.seq_042"} becomes \code{"type1"}).
#'
#' @param sequence_names Character vector of sequence names.
#'
#' @return Character vector with the class label extracted from each sequence
#'   name.
#'
#' @keywords internal
#' 
.class_from_sequence_name <- function(sequence_names) {
  sub("\\..*", "", sequence_names)
}

################################################################################
# Internal helper functions of train_models_rf_xgboost.R
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
# Internal helper functions of kmer_analysis.R
################################################################################
#' .create_motifs_rank: Rank Motifs by Class Frequency 
#'
#' Builds a formatted summary table of motif occurrence sums and counts per
#' class.
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
  motif_cols    <- setdiff(seq_len(ncol(data)), class_col_idx)

  classes        <- data[[class_col]]
  motifs         <- data[, motif_cols, drop = FALSE]
  motif_names    <- colnames(motifs)
  unique_classes <- unique(classes)

  sums_mat <- matrix(
    0,
    nrow     = length(motif_names),
    ncol     = length(unique_classes),
    dimnames = list(motif_names, unique_classes)
  )
  counts_mat <- matrix(
    0L,
    nrow     = length(motif_names),
    ncol     = length(unique_classes),
    dimnames = list(motif_names, unique_classes)
  )

  for (class_name in unique_classes) {
    class_rows               <- classes == class_name
    class_data               <- motifs[class_rows, , drop = FALSE]
    sums_mat[, class_name]   <- colSums(class_data, na.rm = TRUE)
    counts_mat[, class_name] <- colSums(class_data > 0, na.rm = TRUE)
  }

  value_width   <- nchar(as.character(ceiling(max(sums_mat))))
  formatted_mat <- matrix(
    "",
    nrow     = nrow(sums_mat),
    ncol     = ncol(sums_mat),
    dimnames = dimnames(sums_mat)
  )

  for (i in seq_len(nrow(sums_mat))) {
    for (j in seq_len(ncol(sums_mat))) {
      formatted_mat[i, j] <- sprintf(
        "%0*d|%d",
        value_width,
        round(sums_mat[i, j]),
        counts_mat[i, j]
      )
    }
  }

  result       <- as.data.frame(formatted_mat, stringsAsFactors = FALSE)
  result$motif <- motif_names
  result$mean  <- round(rowMeans(sums_mat, na.rm = TRUE), 2)
  result       <- result[, c("motif", unique_classes, "mean")]
  rownames(result) <- NULL
  return(result)
}

################################################################################
# Internal helper functions of kmers_in_seq.R
################################################################################
#'
#' .read_fasta_sequences: Read and consolidate FASTA sequences from directory
#'
#' @param input_dir Character string specifying the directory containing FASTA files
#' @return A list with two elements: 'sequences' (character vector of uppercase DNA
#'   sequences) and 'names' (character vector of sequence identifiers)
#' @keywords internal
.read_fasta_sequences <- function(input_dir) {
  fasta_files <- list.files(input_dir, pattern = "\\.(fasta|fa|fna)$",
                            full.names = TRUE, ignore.case = TRUE)

  all_sequences <- list()
  all_names     <- character()

  for (fasta_file in fasta_files) {
    seqs <- seqinr::read.fasta(fasta_file, seqtype = "DNA", as.string = TRUE,
                               forceDNAtolower = FALSE)
    for (seq_name in names(seqs)) {
      all_sequences[[seq_name]] <- as.character(toupper(seqs[[seq_name]]))
      all_names <- c(all_names, seq_name)
    }
  }

  return(list(sequences = unlist(all_sequences, use.names = FALSE), names = all_names))
}

################################################################################

#' .find_motifs_parallel: Find motif occurrences in sequences using parallel processing
#'
#' @param sequences Character vector of DNA sequences
#' @param motifs Character vector of motif patterns to search
#' @param sequence_names Character vector of sequence identifiers
#' @param class_lookup Named vector for mapping sequence names to classes
#' @param length_lookup Named vector for mapping sequence names to lengths
#' @param n_cores Integer specifying number of CPU cores to use
#' @param dataset Character or NULL. Dataset label ("training" or "test").
#'   Default: NULL.
#'
#' @return A data.frame with columns: motif, sequence_name, class, position_start,
#'   position_end, sequence_length, and optionally dataset
#'
#' @keywords internal
.find_motifs_parallel <- function(sequences, motifs, sequence_names,
                                  class_lookup, length_lookup, n_cores, dataset = NULL) {

  process_motif <- function(motif) {
    results_list <- list()
    counter <- 0L

    for (j in seq_along(sequences)) {
      positions <- stringi::stri_locate_all_fixed(sequences[j], motif)[[1]]

      if (!is.na(positions[1L, 1L])) {
        counter  <- counter + 1L
        seq_name <- sequence_names[j]

        seq_class <- if (!is.null(class_lookup)) {
          cls <- class_lookup[seq_name]
          if (is.na(cls)) NA_character_ else as.character(cls)
        } else {
          NA_character_
        }

        seq_length <- if (!is.null(length_lookup)) {
          len <- length_lookup[seq_name]
          if (is.na(len)) NA_integer_ else as.integer(len)
        } else {
          NA_integer_
        }

        temp_df <- data.frame(
          motif           = rep(motif,      nrow(positions)),
          sequence_name   = rep(seq_name,   nrow(positions)),
          class           = rep(seq_class,  nrow(positions)),
          position_start  = as.integer(positions[, 1L]),
          position_end    = as.integer(positions[, 2L]),
          sequence_length = rep(seq_length, nrow(positions)),
          stringsAsFactors = FALSE
        )

        if (!is.null(dataset)) {
          temp_df$dataset <- dataset
        }

        results_list[[counter]] <- temp_df
      }
    }

    if (counter > 0L) do.call(rbind, results_list) else NULL
  }

  # Parallel or sequential
  if (n_cores > 1L) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      c("sequences", "sequence_names", "class_lookup", "length_lookup", "dataset"),
      envir = environment()
    )
    
    results_list <- parallel::parLapply(cl, motifs, process_motif)
  } else {
    results_list <- lapply(motifs, process_motif)
  }

  # Combine
  results_list <- results_list[!sapply(results_list, is.null)]

  if (length(results_list) > 0L) {
    result_df <- do.call(rbind, results_list)
    rownames(result_df) <- NULL
  } else {
    result_df <- data.frame(
      motif           = character(0),
      sequence_name   = character(0),
      class           = character(0),
      position_start  = integer(0),
      position_end    = integer(0),
      sequence_length = integer(0),
      stringsAsFactors = FALSE
    )
    if (!is.null(dataset)) result_df$dataset <- character(0)
  }

  return(result_df)
}

################################################################################
# Internal helper functions of seq_classification_cent.R
################################################################################
#' .print_step: Print Pipeline Step Header
#'
#' Prints a formatted step header (step number, total steps, and title) to
#' the console, followed by an optional details message. Used to report
#' progress through the sequence classification pipeline.
#'
#' @param step Integer. Current step number.
#' @param total_steps Integer. Total number of steps in the pipeline.
#' @param title Character. Short description of the current step.
#' @param verbose Logical. If FALSE, the function returns immediately without
#'   printing anything. Default: TRUE.
#' @param details Character or NULL. Optional additional message printed
#'   after the step header. Default: NULL.
#'
#' @return NULL (invisibly). Called for its side effect of printing progress
#'   messages.
#'
#' @keywords internal
.print_step <- function(step, total_steps, title, verbose = TRUE, details = NULL) {
  if (!verbose) return()
  message(sprintf("STEP %d/%d: %s", step, total_steps, title))
  message(paste(rep("-", 80), collapse = ""))
  if (!is.null(details)) message(details)
} 