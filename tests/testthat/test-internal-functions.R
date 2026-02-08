# tests/testthat/test-internal-functions.R

library(testthat)

test_that(".count_kmers works correctly with valid DNA sequences", {
  # Test basic k-mer counting
  sequence <- "ATCGATCG"
  result <- .count_kmers(sequence, k = 2, alphabet = c("A", "T", "C", "G"))

  expect_type(result, "integer")
  expect_named(result)
  expect_equal(sum(result), 7) # 8 - 2 + 1 = 7 k-mers
  expect_true(result["AT"] == 2L)
  expect_true(result["TC"] == 2L)
  expect_true(result["CG"] == 2L)
})

test_that(".count_kmers handles invalid characters", {
  sequence <- "ATCNGATCG"
  result <- .count_kmers(sequence, k = 2, alphabet = c("A", "T", "C", "G"))

  expect_type(result, "integer")
  # Should remove N and continue
  expect_true(all(result >= 0))
})

test_that(".count_kmers handles Biostrings objects", {
  skip_if_not_installed("Biostrings")

  dna_string <- Biostrings::DNAString("ATCG")
  result <- .count_kmers(dna_string, k = 2, alphabet = c("A", "T", "C", "G"))

  expect_type(result, "integer")
  expect_true(result["AT"] == 1L)
  expect_true(result["TC"] == 1L)
  expect_true(result["CG"] == 1L)
})

test_that(".process_sequences returns correct structure", {
  skip_if_not_installed("Biostrings")

  sequences <- Biostrings::DNAStringSet(c("ATCG", "GCTA"))
  names(sequences) <- c("seq1", "seq2")
  all_kmers <- c("AT", "TC", "CG", "GC", "CT", "TA", "AA", "TT", "CC", "GG",
                 "AC", "AG", "CA", "GA", "TG", "GT")

  result <- .process_sequences(sequences, all_kmers, k = 2,
                               alphabet = c("A", "T", "C", "G"),
                               class_name = "class1")

  expect_type(result, "list")
  expect_named(result, c("kmers", "metadata"))
  expect_s3_class(result$kmers, "data.frame")
  expect_s3_class(result$metadata, "data.frame")
  expect_equal(nrow(result$kmers), 2)
  expect_equal(nrow(result$metadata), 2)
  expect_true("CLASS" %in% colnames(result$kmers))
  expect_equal(result$metadata$class, c("class1", "class1"))
})

test_that(".count_leaves counts dendrogram leaves correctly", {
  # Create simple dendrogram
  hc <- hclust(dist(matrix(1:6, ncol = 2)))
  dend <- as.dendrogram(hc)

  result <- .count_leaves(dend)
  expect_type(result, "integer")
  expect_equal(result, 3L)
})

test_that(".calc_homogeneity calculates correctly", {
  # Perfect homogeneity
  expect_equal(.calc_homogeneity(c("A", "A", "A")), 1.0)

  # Mixed classes
  expect_equal(.calc_homogeneity(c("A", "A", "B")), 2/3)

  # Empty vector
  expect_true(is.na(.calc_homogeneity(character(0))))
})

test_that(".get_dominant returns most frequent class", {
  expect_equal(.get_dominant(c("A", "A", "B")), "A")
  expect_equal(.get_dominant(c("B", "A", "B", "C")), "B")
  expect_true(is.na(.get_dominant(character(0))))
})

test_that(".find_segments identifies contiguous segments", {
  expect_equal(.find_segments(c(1, 2, 3, 7, 8)), list(c(1, 2, 3), c(7, 8)))
  expect_equal(.find_segments(c(1)), list(c(1)))
  expect_equal(.find_segments(c(1, 3, 5)), list(c(1), c(3), c(5)))
})

test_that(".is_complete_class_cluster works correctly", {
  total_counts <- c(A = 3, B = 2)

  # Complete cluster for class A
  expect_true(.is_complete_class_cluster(c("A", "A", "A"), total_counts))

  # Incomplete cluster
  expect_false(.is_complete_class_cluster(c("A", "A", "B"), total_counts))
  expect_false(.is_complete_class_cluster(c("A", "A"), total_counts))
})

test_that(".create_cluster creates correct structure", {
  result <- .create_cluster(c(1, 2, 3), c("A", "A", "B"), counter = 5)

  expect_type(result, "list")
  expect_named(result, c("indices", "classes", "id"))
  expect_equal(result$indices, c(1, 2, 3))
  expect_equal(result$classes, c("A", "A", "B"))
  expect_equal(result$id, 5)
})

test_that(".normalize_motif_matrix normalizes correctly", {
  # Test matrix with varying ranges
  mat <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3, ncol = 2)
  result <- .normalize_motif_matrix(mat)

  expect_true(all(result >= 0 & result <= 1))
  expect_equal(dim(result), dim(mat))
  expect_equal(result[1, 1], 0) # minimum value
  expect_equal(result[3, 1], 1) # maximum value
})

test_that(".normalize_motif_matrix handles zero range columns", {
  # Column with all same values
  mat <- matrix(c(1, 2, 3, 5, 5, 5), nrow = 3, ncol = 2)
  result <- .normalize_motif_matrix(mat)

  expect_equal(result[, 2], c(0.5, 0.5, 0.5))
})

test_that(".select_by_class_fast distributes motifs correctly", {
  motif_cluster <- data.frame(
    cluster1 = c(10, 5, 3, 1),
    cluster2 = c(2, 8, 6, 4),
    row.names = c("motif1", "motif2", "motif3", "motif4")
  )
  cluster_to_class <- c(cluster1 = "A", cluster2 = "B")
  classe_order <- c("A", "B")

  result <- .select_by_class_fast(motif_cluster, cluster_to_class,
                                  classe_order, n = 4, verbose = FALSE)

  expect_type(result, "list")
  expect_named(result, c("A", "B"))
  expect_equal(length(unlist(result)), 4)
})

test_that(".select_by_cluster_fast distributes motifs per cluster", {
  motif_cluster <- data.frame(
    cluster1 = c(10, 5, 3, 1),
    cluster2 = c(2, 8, 6, 4),
    row.names = c("motif1", "motif2", "motif3", "motif4")
  )
  cluster_to_class <- c(cluster1 = "A", cluster2 = "B")
  classe_order <- c("A", "B")

  result <- .select_by_cluster_fast(motif_cluster, cluster_to_class,
                                    classe_order, n = 4, m = 2)

  expect_type(result, "list")
  expect_equal(length(unlist(result)), 4)
})

test_that(".select_by_class selects correct number of sequences", {
  data <- data.frame(
    seq = paste0("seq", 1:10),
    class = rep(c("A", "B"), each = 5),
    stringsAsFactors = FALSE
  )

  result <- .select_by_class(data, k = 3)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 6) # 3 per class
})

test_that(".select_by_class handles insufficient sequences", {
  data <- data.frame(
    seq = paste0("seq", 1:3),
    class = c("A", "A", "B"),
    stringsAsFactors = FALSE
  )

  result <- .select_by_class(data, k = 5)

  expect_equal(nrow(result), 3) # Returns all available
})

test_that(".rename_class_column renames CLASS to class", {
  data <- data.frame(CLASS = c("A", "B"), value = c(1, 2))
  result <- .rename_class_column(data)

  expect_true("class" %in% colnames(result))
  expect_false("CLASS" %in% colnames(result))
})

test_that(".rename_class_column leaves data unchanged if no CLASS column", {
  data <- data.frame(type = c("A", "B"), value = c(1, 2))
  result <- .rename_class_column(data)

  expect_equal(result, data)
})

test_that(".validate_columns stops on missing columns", {
  data <- data.frame(a = 1, b = 2)

  expect_error(
    .validate_columns(data, c("a", "c"), "test_data"),
    "Missing columns in test_data: c"
  )
})

test_that(".validate_columns passes with all required columns", {
  data <- data.frame(a = 1, b = 2, c = 3)

  expect_silent(.validate_columns(data, c("a", "b"), "test_data"))
})

test_that(".create_motifs_rank creates correct structure", {
  data <- data.frame(
    motif1 = c(5, 0, 3),
    motif2 = c(2, 4, 0),
    CLASS = c("A", "B", "A")
  )

  result <- .create_motifs_rank(data, class_col = "CLASS")

  expect_s3_class(result, "data.frame")
  expect_true("motif" %in% colnames(result))
  expect_true("mean" %in% colnames(result))
  expect_true(all(c("A", "B") %in% colnames(result)))
  expect_equal(nrow(result), 2) # 2 motifs
})

test_that(".create_motifs_rank formats values correctly", {
  data <- data.frame(
    motif1 = c(5, 0, 3),
    CLASS = c("A", "A", "B")
  )

  result <- .create_motifs_rank(data, class_col = "CLASS")

  # Check format: sum|count
  expect_match(result$A[1], "\\d+\\|\\d+")
  expect_match(result$B[1], "\\d+\\|\\d+")
})

test_that(".read_fasta_sequences returns correct structure", {
  skip_if_not_installed("seqinr")

  # Create temporary FASTA file
  temp_dir <- tempdir()
  temp_fasta <- file.path(temp_dir, "test.fasta")
  writeLines(c(">seq1", "ATCG", ">seq2", "GCTA"), temp_fasta)

  result <- .read_fasta_sequences(temp_dir)

  expect_type(result, "list")
  expect_named(result, c("sequences", "names"))
  expect_type(result$sequences, "character")
  expect_type(result$names, "character")
  expect_true(length(result$sequences) >= 2)

  # Cleanup
  unlink(temp_fasta)
})

test_that(".find_motifs_parallel finds motifs correctly", {
  skip_if_not_installed("stringi")

  sequences <- c("ATCGATCG", "GCTAGCTA")
  motifs <- c("ATC", "GCT")
  sequence_names <- c("seq1", "seq2")
  class_lookup <- c(seq1 = "A", seq2 = "B")
  length_lookup <- c(seq1 = 8L, seq2 = 8L)

  result <- .find_motifs_parallel(sequences, motifs, sequence_names,
                                  class_lookup, length_lookup, n_cores = 1)

  expect_s3_class(result, "data.frame")
  expect_true(all(c("motif", "sequence_name", "class", "position_start",
                    "position_end", "sequence_length") %in% colnames(result)))
  expect_true(nrow(result) >= 2)
})

test_that(".find_motifs_parallel handles dataset parameter", {
  skip_if_not_installed("stringi")

  sequences <- c("ATCG")
  motifs <- c("ATC")
  sequence_names <- c("seq1")

  result <- .find_motifs_parallel(sequences, motifs, sequence_names,
                                  class_lookup = NULL, length_lookup = NULL,
                                  n_cores = 1, dataset = "training")

  expect_true("dataset" %in% colnames(result))
  expect_equal(result$dataset[1], "training")
})

test_that(".find_motifs_parallel returns empty dataframe when no matches", {
  skip_if_not_installed("stringi")

  sequences <- c("AAAA")
  motifs <- c("GGG")
  sequence_names <- c("seq1")

  result <- .find_motifs_parallel(sequences, motifs, sequence_names,
                                  class_lookup = NULL, length_lookup = NULL,
                                  n_cores = 1)

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 0)
})

